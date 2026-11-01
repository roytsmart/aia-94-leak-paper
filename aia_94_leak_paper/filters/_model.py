"""
How much of the leak, and of the EUV, the thick focal-plane filter passes
relative to the thin one.

The frame through the thick filter, :func:`aia_94_leak_paper.filters.frame`,
is compared with the run through the thin filter after it,
:func:`aia_94_leak_paper.filters.run`, over the disk. Across its regions, the
signal through the thick filter is fit as a multiple :math:`\\tau` of the
signal through the thin one plus a multiple of the leak through the thin
one, :math:`L`,

.. math::

    S_\\mathrm{thick} = \\tau S_\\mathrm{thin} + c L,

where :math:`c L` takes up whatever the thick filter does to the leak
beyond what it does to the EUV. The structure of the corona fixes
:math:`\\tau`, the attenuation of the EUV, and the attenuation of the leak
is :math:`\\tau + c`. A single frame has no noise of its own to measure,
so this is the only measure of the leak through the thick filter. The leak
through the thin one, :func:`leak_thin`, is measured two ways,
:data:`measures`: by the noise of the run, and by the eclipses that bracket
the day.

Light that crossed the metal of the filters would be attenuated by the
extra zirconium of the thick one, :func:`metal`, which the optical
constants of zirconium give at any wavelength: in the EUV as the
calibration expects, and in the visible almost entirely, by every
measurement of the constants, :data:`tables`.
"""

import dataclasses
import itertools

import astropy.units as u
import named_arrays as na
import numpy as np
import optika

import aia_94_leak_paper

__all__ = [
    "Attenuation",
    "attenuation",
    "euv_table",
    "leak_thin",
    "measures",
    "metal",
    "metal_visible",
    "radius",
    "tables",
    "visible",
    "wavelength_euv",
    "zirconium",
]

radius: float = 0.9
"""Regions within this distance of disk center, in solar radii, are compared; nearer the limb they straddle it."""

measures: tuple[str, ...] = ("noise", "eclipses")
"""
The measures of the leak through the thin filter on the day, :func:`leak_thin`:
the noise of the run through it, and the eclipses that bracket the day.
"""

zirconium: dict[str, u.Quantity] = {"thin": 2000 * u.AA, "thick": 3000 * u.AA}
"""
The zirconium of each focal-plane filter of the channel, as built
(Lemen et al. 2012, Table 3); the thick filter is a spare, in case the thin
one develops pinholes.
"""

euv_table: float = 0.193237 / 0.26922
"""
The attenuation of the EUV by the thick filter that the calibration of the
channel expects: the ratio of its effective areas through the thick and
thin focal-plane filters, ``94_THICK`` and ``94_THIN`` in version 10 of the
correction table at JSOC. Every version of the table has the same ratio,
set by the transmissions of the two filters measured before launch.
"""

wavelength_euv: u.Quantity = 94 * u.AA
"""The wavelength of the channel, which its effective areas are given at."""

visible: tuple[u.Quantity, u.Quantity] = (4000 * u.AA, 11000 * u.AA)
"""
The wavelengths of visible light the leak could be made of: from the
violet end of the visible to 1.1 µm, the wavelength of the band gap of
silicon, beyond which the detector is transparent.
"""

tables: tuple[str, ...] = ("querry", "lynch", "palm")
"""
The measurements of the optical constants of zirconium in the visible, as
the names of their tables in optika (:attr:`optika.chemicals.Chemical.table`):
Querry (1987), optika's default; Lynch, Olson & Weaver (1975); and the film
of Palm et al. (2018), the least absorbing of them.
"""


@dataclasses.dataclass(frozen=True)
class Attenuation:
    """
    How much the thick focal-plane filter passes relative to the thin one,
    over the disk.
    """

    leak: float
    """The attenuation of the leak, :math:`\\tau + c`."""

    leak_quadrants: np.ndarray
    """The same, fit in each quadrant of the detector."""

    euv: float
    """The attenuation of the EUV, :math:`\\tau`, from the structure of the corona."""

    center: u.Quantity
    """The leak through the thin filter at disk center, the median within a tenth of a radius."""

    @property
    def leak_spread(self) -> float:
        """The standard deviation of :attr:`leak_quadrants`."""
        return float(np.std(self.leak_quadrants))


def _fit(
    thin: np.ndarray,
    thick: np.ndarray,
    leak: np.ndarray,
) -> tuple[float, float]:
    """
    The attenuation of the EUV, :math:`\\tau`, and the multiple :math:`c` of
    the leak through the thin filter, fit to the signals of some regions.

    Parameters
    ----------
    thin
        The signal of each region through the thin filter.
    thick
        The signal of each region through the thick filter.
    leak
        The leak through the thin filter in each region.
    """
    design = np.stack([thin, leak], axis=-1)
    (tau, c), *_ = np.linalg.lstsq(design, thick, rcond=None)
    return float(tau), float(c)


def _leak_eclipses(run: "aia_94_leak_paper.noise.Run") -> np.ndarray:
    """
    The leak through the thin filter on the day, in the units of the detector
    per frame, in each region of the run, from the two eclipses that bracket
    the day.

    The leak of each eclipse, in level 1, is grown to the day along the
    exponential of the eclipses and taken as a profile in rings of a
    hundredth of a radius, out to 1.3 radii, since the Sun is not the same
    size in each. It is then multiplied by the flat field in each region,
    which is fixed to the detector, the mean of the flat fields of the two
    eclipses, to put it in the units of the detector, as the run is.

    Parameters
    ----------
    run
        The run through the thin filter, :func:`aia_94_leak_paper.filters.run`.
    """
    eclipses = aia_94_leak_paper.eclipses
    day = aia_94_leak_paper.filters.switch.date()
    doubling = eclipses.growth(eclipses.history()).doubling
    edges = np.linspace(0, 1.3, 131)
    middle = (edges[:-1] + edges[1:]) / 2
    axes = ("y", "x")
    shape = run.r.shape
    profiles, flats = [], []
    for date in (
        max(d for d in eclipses.dates if d <= day),
        min(d for d in eclipses.dates if d >= day),
    ):
        eclipse = eclipses.load(date)
        growth = 2 ** ((day - date).days * u.day / doubling).to_value(u.one)
        # the leak with the whole disk up, in level 1
        seen = eclipses.seen(date) / eclipse.flat_center * run.exposure * growth
        seen = seen.to_value(u.DN).ndarray_aligned(axes)
        position = eclipse.position
        distance = np.sqrt(np.square(position.x) + np.square(position.y))
        distance = distance.ndarray_aligned(axes)
        profiles.append(
            [
                np.nanmedian(seen[(distance >= lo) & (distance < hi)])
                for lo, hi in itertools.pairwise(edges)
            ]
        )
        # the bins of the eclipses tile the regions of the run
        flat = eclipse.flat.ndarray_aligned(axes)
        ny, nx = (flat.shape[i] // shape[a] for i, a in enumerate(axes))
        flats.append(
            flat.reshape(shape["y"], ny, shape["x"], nx).mean(axis=(1, 3)),
        )
    r = run.r.ndarray_aligned(axes)
    level_1 = np.interp(r, middle, np.mean(profiles, axis=0), right=np.nan)
    return level_1 * np.mean(flats, axis=0)


def leak_thin(measure: str = "noise") -> na.ScalarArray:
    """
    The leak through the thin filter in each region of the run through it,
    per frame.

    Parameters
    ----------
    measure
        One of :data:`measures`: ``"noise"``, the limb-darkened disk of the
        noiseless signal of the run, :func:`aia_94_leak_paper.noise.disk`,
        which is nothing beyond the limb; or ``"eclipses"``, the leak in the
        eclipses that bracket the day, through the point spread function of
        the eclipses, which spreads some of it beyond the limb.
    """
    thin = aia_94_leak_paper.filters.run()
    axes = ("y", "x")
    r = thin.r.ndarray_aligned(axes)
    if measure == "noise":
        disk = aia_94_leak_paper.noise.disk(thin)
        inside = disk(np.clip(r, 0, 1)).to_value(u.DN)
        leak = np.where(r < 1, inside, 0)
    elif measure == "eclipses":
        leak = _leak_eclipses(thin)
    else:  # pragma: no cover
        raise ValueError(f"measure must be one of {measures}, got {measure!r}")
    return na.ScalarArray(leak, axes=axes) << u.DN


def attenuation(measure: str = "noise") -> Attenuation:
    """
    How much the thick focal-plane filter passes of the leak and of the EUV, relative to the thin one.

    Parameters
    ----------
    measure
        The measure of the leak through the thin filter, one of :data:`measures`.
    """
    filters = aia_94_leak_paper.filters
    noise = aia_94_leak_paper.noise
    thin = filters.run()
    axes = ("y", "x")
    s_thin = thin.image.to_value(u.DN).ndarray_aligned(axes)
    s_thick = filters.frame("clean").image.to_value(u.DN).ndarray_aligned(axes)
    leak = leak_thin(measure).to_value(u.DN).ndarray_aligned(axes)
    r = thin.r.ndarray_aligned(axes)
    use = (r < radius) & np.isfinite(s_thin) & np.isfinite(s_thick)

    tau, c = _fit(s_thin[use], s_thick[use], leak[use])
    quadrants = []
    for quadrant in noise.quadrants(thin):
        w = use & quadrant.ndarray_aligned(axes)
        tau_q, c_q = _fit(s_thin[w], s_thick[w], leak[w])
        quadrants.append(tau_q + c_q)

    return Attenuation(
        leak=tau + c,
        leak_quadrants=np.array(quadrants),
        euv=tau,
        center=np.median(leak[r < 0.1]) * u.DN,
    )


def metal(
    wavelength: u.Quantity | na.AbstractScalar,
    table: None | str = None,
) -> na.AbstractScalar:
    """
    How much the thick focal-plane filter would pass, relative to the thin
    one, of light that crossed the metal: what its extra zirconium transmits,

    .. math::

        e^{-4 \\pi k \\, \\Delta d / \\lambda},

    where :math:`k` is the extinction coefficient of zirconium at the
    wavelength :math:`\\lambda` and :math:`\\Delta d` is how much more
    zirconium the thick filter has. The reflections at the surfaces of the
    two filters, and the oxide on them, are the same, and so are left out
    of the ratio.

    Parameters
    ----------
    wavelength
        The wavelength of the light in vacuum.
    table
        The table of the optical constants of zirconium in optika, one of
        :data:`tables`, or :obj:`None` for optika's default.
    """
    chemical = optika.chemicals.Chemical("Zr", table=table)
    extra = zirconium["thick"] - zirconium["thin"]
    depth = chemical.absorption(na.as_named_array(wavelength)) * extra
    return np.exp(-depth.to(u.one).value)


def metal_visible(table: None | str = None) -> float:
    """
    The most that the thick focal-plane filter would pass, relative to the
    thin one, of light in :data:`visible` that crossed the metal,
    :func:`metal`.

    Parameters
    ----------
    table
        The table of the optical constants of zirconium in optika, one of
        :data:`tables`, or :obj:`None` for optika's default.
    """
    wavelength = na.geomspace(*visible, axis="wavelength", num=1001)
    return float(metal(wavelength, table).max().ndarray)

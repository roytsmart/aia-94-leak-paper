"""
How much of the leak, and of the EUV, the thick focal-plane filter passes
relative to the thin one.

The frame through the thick filter, :func:`aia_94_leak_paper.filters.frame`,
is compared with the run through the thin filter after it,
:func:`aia_94_leak_paper.filters.run`, over the disk. Across its regions, the
signal through the thick filter is fit as a multiple :math:`\\tau` of the
signal through the thin one plus a limb-darkened disk,

.. math::

    S_\\mathrm{thick} = \\tau S_\\mathrm{thin} + c \\left[ 1 - u (1 - \\mu) \\right],

where the disk, with the limb darkening :math:`u` of the leak through the
thin filter, takes up whatever the thick filter does to the leak beyond
what it does to the EUV. The structure of the corona fixes :math:`\\tau`,
the attenuation of the EUV. Since the disk is :math:`c` at its center, the
noiseless signal through the thin filter there, :math:`A` of
:func:`aia_94_leak_paper.noise.disk`, gives the attenuation of the leak,
:math:`\\tau + c / A`. A single frame has no noise of its own to measure,
so this is the only measure of the leak through the thick filter.

Light that crossed the metal of the filters would be attenuated by the
extra zirconium of the thick one, :func:`metal`, which the optical
constants of zirconium give at any wavelength: in the EUV as the
calibration expects, and in the visible almost entirely, by every
measurement of the constants, :data:`tables`.
"""

import dataclasses

import astropy.units as u
import named_arrays as na
import numpy as np
import optika

import aia_94_leak_paper

__all__ = [
    "Attenuation",
    "attenuation",
    "band",
    "band_peak",
    "euv_table",
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
    """The attenuation of the leak, :math:`\\tau + c / A`."""

    leak_quadrants: np.ndarray
    """The same, fit in each quadrant of the detector."""

    euv: float
    """The attenuation of the EUV, :math:`\\tau`, from the structure of the corona."""

    noiseless_center: u.Quantity
    """The noiseless signal through the thin filter at disk center, :math:`A`."""

    @property
    def leak_spread(self) -> float:
        """The standard deviation of :attr:`leak_quadrants`."""
        return float(np.std(self.leak_quadrants))


def _fit(
    thin: np.ndarray,
    thick: np.ndarray,
    r: np.ndarray,
    limb_darkening: float,
) -> tuple[float, float]:
    """
    The attenuation of the EUV, :math:`\\tau`, and the disk, :math:`c`, fit
    to the signals of some regions.

    Parameters
    ----------
    thin
        The signal of each region through the thin filter.
    thick
        The signal of each region through the thick filter.
    r
        The distance of each region from the center of the disk, in solar radii.
    limb_darkening
        The limb darkening :math:`u` of the leak through the thin filter.
    """
    mu = np.sqrt(1 - np.square(r))
    design = np.stack([thin, 1 - limb_darkening * (1 - mu)], axis=-1)
    (tau, c), *_ = np.linalg.lstsq(design, thick, rcond=None)
    return float(tau), float(c)


def attenuation() -> Attenuation:
    """How much the thick focal-plane filter passes of the leak and of the EUV, relative to the thin one."""
    filters = aia_94_leak_paper.filters
    noise = aia_94_leak_paper.noise
    thin = filters.run()
    axes = ("y", "x")
    s_thin = thin.image.to_value(u.DN).ndarray_aligned(axes)
    s_thick = filters.frame("clean").image.to_value(u.DN).ndarray_aligned(axes)
    r = thin.r.ndarray_aligned(axes)
    disk = noise.disk(thin)
    center = disk.amplitude.to_value(u.DN)
    use = (r < radius) & np.isfinite(s_thin) & np.isfinite(s_thick)

    tau, c = _fit(s_thin[use], s_thick[use], r[use], disk.limb_darkening)
    quadrants = []
    for quadrant in noise.quadrants(thin):
        w = use & quadrant.ndarray_aligned(axes)
        tau_q, c_q = _fit(s_thin[w], s_thick[w], r[w], disk.limb_darkening)
        quadrants.append(tau_q + c_q / center)

    return Attenuation(
        leak=tau + c / center,
        leak_quadrants=np.array(quadrants),
        euv=tau,
        noiseless_center=disk.amplitude,
    )


def band() -> na.ScalarArray:
    """
    The band in each region: the frame through the thick filter that was
    read out as the filter wheel moved between the two thick filters, less
    the one read out as it moved to the thin filter of 304 Å.
    """
    filters = aia_94_leak_paper.filters
    return filters.frame("banded").image - filters.frame("clean").image


def band_peak() -> u.Quantity:
    """The most the band adds to a frame: the largest median of :func:`band` along a row of the regions on the disk."""
    axes = ("y", "x")
    b = band().to_value(u.DN).ndarray_aligned(axes)
    disk = aia_94_leak_paper.filters.run().r.ndarray_aligned(axes) < 1
    rows = [np.median(row[d]) for row, d in zip(b, disk) if d.any()]
    return max(rows) * u.DN


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

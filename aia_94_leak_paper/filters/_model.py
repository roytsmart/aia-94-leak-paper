"""
How much of the leak, and of the EUV, the thick focal-plane filter passes
relative to the thin one.

The thick filter does not transmit evenly over the field: relative to the
thin one, the signal through it falls toward the top and bottom of the
detector. Its flat field is that of the thin filter, so the level-1 images
do not correct for this, and the two filters are compared only within a
band across the middle of the disk, :data:`band`, where the thick one
transmits evenly.

Within the band, the noiseless signal of each run,
:func:`aia_94_leak_paper.noise.noiseless`, measures the leak through each
filter. The EUV through each is measured from the structure of the corona:
across the regions of the band, the signal through the thick filter is fit
as a multiple :math:`\\tau` of the signal through the thin one plus a
limb-darkened disk,

.. math::

    S_\\mathrm{thick} = \\tau S_\\mathrm{thin} + c \\left[ 1 - u (1 - \\mu) \\right],

where the disk, with the limb darkening :math:`u` of the leak through the
thin filter, takes up whatever the thick filter does to the leak beyond
what it does to the EUV. The structure of the corona fixes :math:`\\tau`,
the attenuation of the EUV, and :math:`c` gives the attenuation of the leak
from the signals alone, :math:`\\tau + c / p_\\mathrm{thin}`.

Light that crossed the metal of the filters would be attenuated by the
extra zirconium of the thick one, :func:`metal`, which the optical
constants of zirconium give at any wavelength. In the EUV this is what the
calibration expects, and in the visible it is far less than the thick
filter passes of the leak, by every measurement of the constants,
:data:`tables`.
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
    "euv_table",
    "metal",
    "metal_visible",
    "tables",
    "visible",
    "wavelength_euv",
    "within",
    "zirconium",
]

band: float = 0.45
"""
Half the height of the band across the disk, in solar radii, within which
the thick filter transmits evenly.
"""

_radius = 0.9
"""Regions of the band within this distance of disk center are compared; nearer the limb they straddle it."""

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
    within :data:`band`.
    """

    leak: float
    """The noiseless signal through the thick filter over that through the thin one."""

    leak_quadrants: np.ndarray
    """The same, in each quadrant of the detector."""

    leak_signals: float
    """The same, from the signals alone, :math:`\\tau + c / p_\\mathrm{thin}`."""

    euv: float
    """The attenuation of the EUV, :math:`\\tau`, from the structure of the corona."""

    noiseless_thin: u.Quantity
    """The median noiseless signal through the thin filter, :math:`p_\\mathrm{thin}`."""

    noiseless_thick: u.Quantity
    """The median noiseless signal through the thick filter."""

    @property
    def leak_spread(self) -> float:
        """The standard deviation of :attr:`leak_quadrants`."""
        return float(np.std(self.leak_quadrants))


def within(run: "aia_94_leak_paper.noise.Run") -> np.ndarray:
    """
    Which regions of a run are compared: those of the band within
    :data:`_radius` of disk center, as a boolean array with the axes ``y`` and ``x``.

    Parameters
    ----------
    run
        A run, as :func:`aia_94_leak_paper.filters.load` gives it.
    """
    where = (np.abs(run.position.y) < band) & (run.r < _radius)
    return where.ndarray_aligned(("y", "x"))


def attenuation() -> Attenuation:
    """How much the thick focal-plane filter passes of the leak and of the EUV, relative to the thin one."""
    filters = aia_94_leak_paper.filters
    noise = aia_94_leak_paper.noise
    thin, thick = filters.load("thin"), filters.load("thick")
    axes = ("y", "x")
    use = within(thin)
    p = {
        name: noise.noiseless(run).to_value(u.DN).ndarray_aligned(axes)
        for name, run in (("thin", thin), ("thick", thick))
    }
    s = {
        name: run.image.to_value(u.DN).ndarray_aligned(axes)
        for name, run in (("thin", thin), ("thick", thick))
    }
    p_thin = np.nanmedian(p["thin"][use])
    p_thick = np.nanmedian(p["thick"][use])
    quadrants = []
    for quadrant in noise.quadrants(thin):
        w = use & quadrant.ndarray_aligned(axes)
        quadrants.append(np.nanmedian(p["thick"][w]) / np.nanmedian(p["thin"][w]))

    # the EUV from the structure of the corona, and the leak from the disk
    r = thin.r.ndarray_aligned(axes)
    ok = use & np.isfinite(s["thin"]) & np.isfinite(s["thick"])
    mu = np.sqrt(1 - np.square(r[ok]))
    limb_darkening = noise.disk(thin).limb_darkening
    design = np.stack([s["thin"][ok], 1 - limb_darkening * (1 - mu)], axis=-1)
    (tau, c), *_ = np.linalg.lstsq(design, s["thick"][ok], rcond=None)

    return Attenuation(
        leak=float(p_thick / p_thin),
        leak_quadrants=np.array(quadrants),
        leak_signals=float(tau + c / p_thin),
        euv=float(tau),
        noiseless_thin=p_thin * u.DN,
        noiseless_thick=p_thick * u.DN,
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

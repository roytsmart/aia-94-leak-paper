"""
The differential emission measure (DEM) of the corona, with the leak left
in 94 Å and taken out.

The DEM of each bin is inverted from the six coronal channels by the
regularized method of Plowman & Caspi (2020), with the temperature response
of SolarSoft, corrected for the degradation of each channel and normalized
to EVE, including the empirical correction of 94 Å for the lines CHIANTI
lacks, and with the uncertainty of SolarSoft. It is inverted twice: from the
images as they are, and with the leak measured in the eclipse of the same
day subtracted from 94 Å. The leak is the limb-darkened disk of that eclipse
seen through its point spread function. It reaches the detector as the EUV
does, so in level 1, flat fielded with the rest of the image, it keeps the
shape the eclipse measures; the eclipse measures it in the units of the
detector at disk center, so it is divided by the flat field there.
"""

import dataclasses
import datetime
import hashlib
import pathlib

import astropy.time
import astropy.units as u
import joblib
import named_arrays as na
import numpy as np
import scipy.ndimage
import sdo
import utu

import aia_94_leak_paper

from ..eclipses._archive import directory_default
from ._build import wavelengths

__all__ = [
    "DEM",
    "Regions",
    "Summary",
    "date_dem",
    "fewest",
    "history",
    "hot",
    "invert",
    "leak",
    "regions",
    "response",
    "summary",
]

date_dem: datetime.date = datetime.date(2026, 1, 19)
"""
The day whose DEMs are shown: the leak is near its brightest, and a large
coronal hole is on the disk.
"""

hot: float = 6.5
"""The log T above which plasma counts as hot, 3.2 MK."""

_temperature = slice(20, 71)
"""Log T 5.0 to 7.5 on the grid of the response, which starts at 4.0 in steps of 0.05."""

_saturation = 15000 * u.DN
"""A bin this bright is at or near the 16383 DN at which the camera saturates."""

_inner = 0.9
"""Regions are chosen within this distance of the center of the disk, in solar radii."""

_smooth = 5
"""The side, in bins, of the square 193 Å is averaged over to find coronal holes."""

_hole = 0.3
"""A coronal hole is darker in 193 Å, so averaged, than this fraction of the median of the disk."""

_quiet = 0.2
"""The quiet Sun is within this fraction of the median of the disk in 193 Å."""

fewest: int = 500
"""The fewest bins of coronal hole whose DEMs are summarized."""

_memory = joblib.Memory(directory_default / "joblib", verbose=0)

_code = hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest()
"""
A digest of this module, passed to every cached function, since `joblib`
keys only on the source of the function it caches and not on the functions
that one calls.
"""


def leak(date: datetime.date) -> na.ScalarArray:
    """
    The leak in the 94 Å image of one day, in the units of level 1.5 per
    second, with the axes ``y`` and ``x``.

    Parameters
    ----------
    date
        A day of :data:`aia_94_leak_paper.eclipses.dates`.
    """
    eclipses = aia_94_leak_paper.eclipses
    synoptic = aia_94_leak_paper.dems.load(date)
    fit = eclipses.fit(date)
    r = synoptic.r
    mu = np.sqrt(np.clip(1 - np.square(r), 0, 1))
    disk = np.where(r < 1, 1 - fit.limb_darkening * (1 - mu), 0)
    axes = ("y", "x")
    image = eclipses.psf()(disk.ndarray_aligned(axes), synoptic.radius)
    inner = (r < 0.1).ndarray_aligned(axes)
    flat_center = np.median(synoptic.flat.ndarray_aligned(axes)[inner])
    return na.ScalarArray(image, axes=axes) * fit.amplitude / flat_center


def response(time: astropy.time.Time) -> na.FunctionArray:
    """
    The temperature response of the six channels at a time, along the axes
    ``wavelength`` and ``temperature``.

    Parameters
    ----------
    time
        When the images were observed.
    """
    result = sdo.aia.temperature_response(
        wavelength=na.ScalarArray(np.array(wavelengths) * u.AA, axes="wavelength"),
        time=time,
        eve=True,
        chiantifix=True,
    )
    return result[{"temperature": _temperature}]


@dataclasses.dataclass(frozen=True)
class DEM:
    """The DEM of every bin of one day."""

    date: datetime.date
    """The day."""

    corrected: bool
    """Whether the leak was subtracted from 94 Å."""

    temperature: na.ScalarArray
    """The temperatures of the DEM, along ``temperature``."""

    values: na.ScalarArray
    """The DEM per unit log T, with the axes ``y``, ``x``, and ``temperature``."""

    chi2: na.ScalarArray
    """The reduced :math:`\\chi^2` of the fit of each bin, negative where it failed."""

    intensity: na.ScalarArray
    """The intensity each channel was fit to, per second."""

    predicted: na.ScalarArray
    """The intensity of each channel the DEM predicts, per second."""

    @property
    def fraction_hot(self) -> na.ScalarArray:
        """The fraction of the emission measure of each bin above :data:`hot`."""
        logt = np.log10(self.temperature.ndarray.to_value(u.K))
        above = np.flatnonzero(logt >= hot)
        values = self.values.to_value(u.cm**-5)
        with np.errstate(divide="ignore", invalid="ignore"):
            return values[{"temperature": slice(above[0], None)}].sum(
                "temperature"
            ) / values.sum("temperature")


def invert(date: datetime.date, corrected: bool) -> DEM:
    """
    Invert the six channels of one day for the DEM of every bin.

    Parameters
    ----------
    date
        A day of :data:`aia_94_leak_paper.eclipses.dates`.
    corrected
        Whether to subtract the leak from 94 Å first.
    """
    synoptic = aia_94_leak_paper.dems.load(date)
    sigma = sdo.aia.uncertainty(
        counts=synoptic.images,
        wavelength=synoptic.wavelength,
        n_sample=synoptic.samples,
    )
    rate = synoptic.rate
    if corrected:
        is_94 = synoptic.wavelength == 94 * u.AA
        rate = rate - np.where(is_94, leak(date), 0 * u.DN / u.s)
    # a saturated bin is left unfit
    rate = np.where(synoptic.images > _saturation, np.nan, rate)
    per_bin = u.DN / u.s / u.pix
    r = response(synoptic.time)
    dem, chi2 = utu.dem.plowman(
        intensity=(rate / u.pix).to(per_bin),
        uncertainty=(sigma / synoptic.exposure / u.pix).to(per_bin),
        response=r,
        axis_channel="wavelength",
        axis_temperature="temperature",
    )
    # the DEM and the response are both piecewise linear in log T, and the
    # trapezoid rule is the integral of their product to within a few percent
    logt = np.log10(dem.inputs.to_value(u.K))
    predicted = np.trapezoid(
        (dem.outputs * r.outputs).ndarray_aligned(rate.axes + ("temperature",)),
        x=logt.ndarray,
        axis=-1,
    )
    return DEM(
        date=date,
        corrected=corrected,
        temperature=dem.inputs,
        values=dem.outputs,
        chi2=chi2,
        intensity=rate,
        predicted=(na.ScalarArray(predicted, axes=rate.axes) * u.pix).to(u.DN / u.s),
    )


@dataclasses.dataclass(frozen=True)
class Regions:
    """The coronal holes and the quiet Sun of one day, chosen in 193 Å."""

    holes: na.ScalarArray
    """Which bins are coronal hole."""

    quiet: na.ScalarArray
    """Which bins are quiet Sun."""


def regions(date: datetime.date) -> Regions:
    """
    The coronal holes and the quiet Sun of one day, within :data:`_inner` of
    the center of the disk.

    Parameters
    ----------
    date
        A day of :data:`aia_94_leak_paper.eclipses.dates`.
    """
    synoptic = aia_94_leak_paper.dems.load(date)
    axes = ("y", "x")
    rate = synoptic.rate[{"wavelength": wavelengths.index(193)}]
    rate = rate.to_value(u.DN / u.s).ndarray_aligned(axes)
    inner = synoptic.r.ndarray_aligned(axes) < _inner
    median = np.median(rate[inner])
    smooth = scipy.ndimage.uniform_filter(rate, size=_smooth)
    return Regions(
        holes=na.ScalarArray(inner & (smooth < _hole * median), axes=axes),
        quiet=na.ScalarArray(inner & (np.abs(rate / median - 1) < _quiet), axes=axes),
    )


@dataclasses.dataclass(frozen=True)
class Summary:
    """
    The DEMs of the coronal holes and the quiet Sun of one day, as observed
    and with the leak subtracted.
    """

    date: datetime.date
    """The day."""

    holes: int
    """The bins of coronal hole."""

    quiet_observed: float
    """The median fraction of the emission measure of the quiet Sun above :data:`hot`, as observed."""

    quiet_corrected: float
    """The same, with the leak subtracted."""

    holes_observed: float
    """The median fraction of the emission measure of the coronal holes above :data:`hot`, as observed, if there are :data:`fewest` bins of them."""

    holes_corrected: float
    """The same, with the leak subtracted."""

    chi2_observed: float
    """The median reduced :math:`\\chi^2` of the quiet Sun, as observed."""

    chi2_corrected: float
    """The same, with the leak subtracted."""

    ratio_94: float
    """The median, over the quiet Sun, of 94 Å with the leak subtracted over the 94 Å its DEM predicts."""


@_memory.cache
def _summary(
    date: datetime.date,
    leak: tuple[float, float, float, float],
    version: str,
    code: str,
) -> dict:
    """The fields of :class:`Summary`, cached on the leak and the archive."""
    chosen = regions(date)
    result = {"date": date, "holes": int(chosen.holes.ndarray.sum())}
    for corrected, suffix in ((False, "observed"), (True, "corrected")):
        dem = invert(date, corrected)
        good = dem.chi2 >= 0
        fraction = dem.fraction_hot
        for name, where, least in (
            ("quiet", chosen.quiet & good, 1),
            ("holes", chosen.holes & good, fewest),
        ):
            values = fraction.ndarray[where.ndarray]
            result[f"{name}_{suffix}"] = (
                float(np.median(values)) if values.size >= least else np.nan
            )
        quiet = (chosen.quiet & good).ndarray
        result[f"chi2_{suffix}"] = float(np.median(dem.chi2.ndarray[quiet]))
        if corrected:
            index = {"wavelength": wavelengths.index(94)}
            ratio = dem.intensity[index] / dem.predicted[index]
            result["ratio_94"] = float(np.median(ratio.to_value(u.one).ndarray[quiet]))
    return result


def summary(date: datetime.date) -> Summary:
    """
    The DEMs of the coronal holes and the quiet Sun of one day, as observed
    and with the leak subtracted.

    Parameters
    ----------
    date
        A day of :data:`aia_94_leak_paper.eclipses.dates`.
    """
    fit = aia_94_leak_paper.eclipses.fit(date)
    psf = aia_94_leak_paper.eclipses.psf()
    leak = (
        float(fit.amplitude.to_value(u.DN / u.s)),
        float(fit.limb_darkening),
        float(psf.eta),
        float(psf.sigma),
    )
    return Summary(
        **_summary(date, leak, aia_94_leak_paper.dems.version, _code),
    )


def history() -> list[Summary]:
    """The DEMs of every day of :data:`aia_94_leak_paper.eclipses.dates`."""
    return [summary(date) for date in aia_94_leak_paper.eclipses.dates]

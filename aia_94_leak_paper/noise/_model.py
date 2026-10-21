"""
The noiseless part of the 94 Å signal, from the variance of a run.

Photon noise makes the variance of a frame grow with its signal :math:`S`,
on top of the read noise :math:`R^2`. If part of the signal, :math:`p`,
carries no photon noise, the variance grows only with the rest:

.. math::

    \\sigma^2 = R^2 + b (S - p),

so a line fit to the variance against the signal crosses the read noise
not at zero signal but at :math:`S = p`. The read noise is measured in the
corners of the detector beyond the field stop, which see no light.
"""

import dataclasses
import datetime

import astropy.time
import astropy.units as u
import named_arrays as na
import numpy as np

import aia_94_leak_paper

__all__ = [
    "Noise",
    "fit",
    "history",
    "noiseless",
    "pooled",
]

_corners = 1.75
"""Regions this far from the center of the disk, in solar radii, are beyond the field stop."""

_center = 0.5
"""Regions within this distance of the center of the disk are fit for the noiseless signal."""

_signal = (2, 20) * u.DN
"""The range of signal the line is fit over, where the leak is a large part of it."""

_count = 3000
"""The fewest residuals a bin of signal needs to be fit."""


def pooled(run: "aia_94_leak_paper.noise.Run", where: na.ScalarArray) -> dict:
    """
    The count, mean signal, and variance of each bin of signal, pooled over
    the regions selected.

    Parameters
    ----------
    run
        A run of the archive.
    where
        Which regions to pool, a boolean array with the axes ``y`` and ``x``.
    """
    axes = ("y", "x")
    count = np.where(where & (run.count > 0), run.count, 0)
    total = count.sum(axes)
    with np.errstate(divide="ignore", invalid="ignore"):
        signal = np.nan_to_num(run.signal * count).sum(axes) / total
        mean = np.nan_to_num(run.mean * count).sum(axes) / total
        square = (
            np.nan_to_num((1.5 * run.variance + np.square(run.mean)) * count).sum(axes)
            / total
        )
    return {
        "count": total,
        "signal": signal,
        "variance": (square - np.square(mean)) / 1.5,
    }


@dataclasses.dataclass(frozen=True)
class Noise:
    """The noise of the 94 Å channel in one run, in the units of the detector."""

    date: datetime.date
    """The day of the run."""

    time: astropy.time.Time
    """When the first frame was observed."""

    exposure: u.Quantity
    """The exposure time of every frame."""

    read: u.Quantity
    """The variance of the read noise of one frame, from the corners."""

    slope: float
    """How much the variance grows per DN of signal, :math:`b`."""

    intercept: u.Quantity
    """The variance the line reaches at zero signal."""

    noiseless: u.Quantity
    """The noiseless signal near disk center in one frame, :math:`p`."""

    @property
    def rate(self) -> u.Quantity:
        """The noiseless signal near disk center per second."""
        return self.noiseless / self.exposure


def _line(table: dict) -> tuple[float, u.Quantity]:
    """Fit variance against signal over :data:`_signal`, weighted by the uncertainty of each variance."""
    signal = table["signal"].ndarray.to_value(u.DN)
    variance = table["variance"].ndarray.to_value(u.DN**2)
    count = table["count"].ndarray
    use = (
        (signal >= _signal[0].value) & (signal <= _signal[1].value) & (count >= _count)
    )
    weight = np.sqrt(count[use]) / variance[use]
    design = np.stack([np.ones(use.sum()), signal[use]], axis=-1)
    (intercept, slope), *_ = np.linalg.lstsq(
        design * weight[:, None], variance[use] * weight, rcond=None
    )
    return float(slope), intercept * u.DN**2


def fit(date: datetime.date) -> Noise:
    """
    The read noise, and the noiseless signal near disk center, in one run.

    Parameters
    ----------
    date
        A day of :data:`aia_94_leak_paper.eclipses.dates`.
    """
    run = aia_94_leak_paper.noise.load(date)
    r = run.r
    corners = pooled(run, r > _corners)
    count = corners["count"].ndarray
    filled = count > 0
    read = np.average(corners["variance"].ndarray[filled], weights=count[filled])
    slope, intercept = _line(pooled(run, r < _center))
    return Noise(
        date=date,
        time=run.time,
        exposure=run.exposure,
        read=read,
        slope=slope,
        intercept=intercept,
        noiseless=(read - intercept) / slope / u.DN,
    )


def noiseless(date: datetime.date) -> na.ScalarArray:
    """
    The noiseless signal of each region of one run, with the slope of the
    line held at its value near disk center.

    Parameters
    ----------
    date
        A day of :data:`aia_94_leak_paper.eclipses.dates`.
    """
    run = aia_94_leak_paper.noise.load(date)
    result = fit(date)
    use = (
        (run.signal >= _signal[0])
        & (run.signal <= 1.5 * _signal[1])
        & (run.count > 100)
    )
    # the intercept of each region, given the slope, weighted by the residuals in each bin
    weight = np.where(use, run.count, 0)
    with np.errstate(divide="ignore", invalid="ignore"):
        intercept = np.nan_to_num(
            (run.variance - result.slope * run.signal * u.DN) * weight
        ).sum("signal") / weight.sum("signal")
    return (result.read - intercept) / result.slope / u.DN


def history() -> list[Noise]:
    """The noise of the run on every day of :data:`aia_94_leak_paper.eclipses.dates`."""
    return [fit(date) for date in aia_94_leak_paper.eclipses.dates]

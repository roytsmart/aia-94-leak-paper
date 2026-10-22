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

The detector is read out through four amplifiers, one per quadrant, whose
read noise and gain differ, so each quadrant is fit with its own read noise
(from its own corner) and its own slope.
"""

import dataclasses
import datetime

import astropy.time
import astropy.units as u
import named_arrays as na
import numpy as np

import aia_94_leak_paper

__all__ = [
    "Disk",
    "Noise",
    "disk",
    "fit",
    "history",
    "measured",
    "noiseless",
    "pooled",
    "quadrants",
    "spread",
]

_corners = 2600
"""
Regions whose centers are this far from the center of the detector, in
pixels, are beyond its field stop, whose light ends near 2500 pixels.
"""

_center = 0.5
"""Regions within this distance of the center of the disk are fit for the noiseless signal."""

_signal = (2, 20) * u.DN
"""The range of signal the line is fit over, where the leak is a large part of it."""

_count = 3000
"""The fewest residuals a bin of signal needs to be fit."""

_disk = 0.95
"""
Regions within this distance of the center of the disk are fit by
:func:`disk`; nearer the limb the regions straddle it.
"""


def quadrants(run: "aia_94_leak_paper.noise.Run") -> list[na.ScalarArray]:
    """
    Which regions lie in each quadrant of the detector, each read out by its
    own amplifier, as boolean arrays with the axes ``y`` and ``x``.

    Parameters
    ----------
    run
        A run of the archive.
    """
    n = run.count.shape["x"]
    y = na.arange(0, n, axis="y") < n // 2
    x = na.arange(0, n, axis="x") < n // 2
    return [y & x, y & ~x, ~y & x, ~y & ~x]


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
    """
    The noise of the 94 Å channel in one run, in the units of the detector,
    for each quadrant of the detector (lower left, lower right, upper left,
    upper right).
    """

    date: datetime.date
    """The day of the run."""

    time: astropy.time.Time
    """When the first frame was observed."""

    exposure: u.Quantity
    """The exposure time of every frame."""

    read: u.Quantity
    """The variance of the read noise of one frame, from the corner of each quadrant."""

    slope: np.ndarray
    """How much the variance grows per DN of signal near disk center, :math:`b`."""

    intercept: u.Quantity
    """The variance the line of each quadrant reaches at zero signal."""

    @property
    def noiseless_quadrants(self) -> u.Quantity:
        """The noiseless signal near disk center in one frame, :math:`p`, in each quadrant."""
        return (self.read - self.intercept) / self.slope / u.DN

    @property
    def noiseless(self) -> u.Quantity:
        """
        The noiseless signal near disk center in one frame, the mean of the
        quadrants that could be fit.
        """
        return np.nanmean(self.noiseless_quadrants)

    @property
    def spread(self) -> u.Quantity:
        """The standard deviation of the noiseless signal among the quadrants that could be fit."""
        return np.nanstd(self.noiseless_quadrants)

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
    if use.sum() < 3:
        # near solar minimum the quiet Sun of a quadrant can be too faint
        return np.nan, np.nan * u.DN**2
    weight = np.sqrt(count[use]) / variance[use]
    design = np.stack([np.ones(use.sum()), signal[use]], axis=-1)
    (intercept, slope), *_ = np.linalg.lstsq(
        design * weight[:, None], variance[use] * weight, rcond=None
    )
    return float(slope), intercept * u.DN**2


def _corner(run: "aia_94_leak_paper.noise.Run") -> na.ScalarArray:
    """The regions of the corners of the detector, beyond the field stop."""
    n = run.count.shape["x"]
    x, y = ((na.arange(0, n, axis=a) + 0.5 - n / 2) * run.region for a in ("x", "y"))
    return np.sqrt(np.square(x) + np.square(y)) > _corners


def _read(run: "aia_94_leak_paper.noise.Run", where: na.ScalarArray) -> u.Quantity:
    """The variance of the read noise in the corner of the regions selected."""
    corner = pooled(run, _corner(run) & where)
    count = corner["count"].ndarray
    filled = count > 0
    return np.average(corner["variance"].ndarray[filled], weights=count[filled])


def fit(run: "aia_94_leak_paper.noise.Run") -> Noise:
    """
    The read noise, and the noiseless signal near disk center, in each
    quadrant of one run.

    Parameters
    ----------
    run
        A run, as :func:`aia_94_leak_paper.noise.load` gives it.
    """
    reads, slopes, intercepts = [], [], []
    for quadrant in quadrants(run):
        reads.append(_read(run, quadrant))
        slope, intercept = _line(pooled(run, (run.r < _center) & quadrant))
        slopes.append(slope)
        intercepts.append(intercept)
    return Noise(
        date=run.date,
        time=run.time,
        exposure=run.exposure,
        read=u.Quantity(reads),
        slope=np.array(slopes),
        intercept=u.Quantity(intercepts),
    )


def noiseless(run: "aia_94_leak_paper.noise.Run") -> na.ScalarArray:
    """
    The noiseless signal of each region of one run, with the read noise and
    the slope of the line held at their values near disk center in its
    quadrant.

    Parameters
    ----------
    run
        A run, as :func:`aia_94_leak_paper.noise.load` gives it.
    """
    result = fit(run)
    use = (
        (run.signal >= _signal[0])
        & (run.signal <= 1.5 * _signal[1])
        & (run.count > 100)
    )
    # the read noise and slope of the quadrant of each region
    read = 0 * u.DN**2
    slope = 0
    for quadrant, r, b in zip(quadrants(run), result.read, result.slope):
        read = np.where(quadrant, r, read)
        slope = np.where(quadrant, b, slope)
    # the intercept of each region, given the slope, weighted by the residuals in each bin
    weight = np.where(use, run.count, 0)
    with np.errstate(divide="ignore", invalid="ignore"):
        intercept = np.nan_to_num(
            (run.variance - slope * run.signal * u.DN) * weight
        ).sum("signal") / weight.sum("signal")
    return (read - intercept) / slope / u.DN


@dataclasses.dataclass(frozen=True)
class Disk:
    """
    The map of the noiseless signal of one run, fit as a limb-darkened disk
    plus a part that grows with the 94 Å signal of each region,

    .. math::

        p = A \\left[ 1 - u (1 - \\mu) \\right] + k (S - \\bar{S}),

    where :math:`\\mu` is the cosine of the angle from disk center and
    :math:`\\bar{S}` is the median signal of the regions fit. If the
    noiseless signal were a part of the 94 Å signal, :math:`k` would be
    that part; if it is light from the photosphere, :math:`k` is zero and
    :math:`u` is the limb darkening of the photosphere.
    """

    date: datetime.date
    """The day of the run."""

    amplitude: u.Quantity
    """The noiseless signal at disk center at the median signal, :math:`A`."""

    limb_darkening: float
    """The coefficient :math:`u`."""

    per_signal: float
    """How much the noiseless signal grows per DN of the signal of a region, :math:`k`."""

    signal: u.Quantity
    """The median signal of the regions fit, :math:`\\bar{S}`."""

    def __call__(self, r: np.ndarray) -> u.Quantity:
        """
        The limb-darkened disk at the median signal.

        Parameters
        ----------
        r
            The distance from the center of the disk, in solar radii.
        """
        mu = np.sqrt(np.clip(1 - np.square(r), 0, 1))
        return self.amplitude * (1 - self.limb_darkening * (1 - mu))


def disk(run: "aia_94_leak_paper.noise.Run") -> Disk:
    """
    Fit the map of the noiseless signal of one run, short of the limb, as a
    limb-darkened disk plus a part that grows with the 94 Å signal.

    Parameters
    ----------
    run
        A run, as :func:`aia_94_leak_paper.noise.load` gives it.
    """
    axes = ("y", "x")
    r = run.r.ndarray_aligned(axes)
    p = noiseless(run).to_value(u.DN).ndarray_aligned(axes)
    s = run.image.to_value(u.DN).ndarray_aligned(axes)
    use = (r < _disk) & np.isfinite(p)
    median = np.median(s[use])
    mu = np.sqrt(1 - np.square(r[use]))
    design = np.stack([np.ones(use.sum()), mu, s[use] - median], axis=-1)
    (a, b, k), *_ = np.linalg.lstsq(design, p[use], rcond=None)
    return Disk(
        date=run.date,
        amplitude=(a + b) * u.DN,
        limb_darkening=float(b / (a + b)),
        per_signal=float(k),
        signal=median * u.DN,
    )


def history() -> list[Noise]:
    """The noise of the run on every day of :data:`aia_94_leak_paper.eclipses.dates`."""
    load = aia_94_leak_paper.noise.load
    return [fit(load(date)) for date in aia_94_leak_paper.eclipses.dates]


def spread() -> u.Quantity:
    """The typical spread of the noiseless signal among the quadrants, the median over the days."""
    return np.nanmedian(u.Quantity([n.spread for n in history()]))


def measured(
    noises: list[Noise],
    leaks: "list[aia_94_leak_paper.eclipses.Leak]",
) -> list[bool]:
    """
    Whether the noise of each day measures the leak: whether the leak in the
    eclipse of that day is three times the typical :func:`spread` of the method.

    Parameters
    ----------
    noises
        The noise of each day, as :func:`history` gives it.
    leaks
        The leak in the eclipse of each day, as
        :func:`aia_94_leak_paper.eclipses.history` gives it.
    """
    s = spread()
    return [bool(k.center * n.exposure > 3 * s) for n, k in zip(noises, leaks)]

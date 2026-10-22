"""
The visible light that passes the entrance filters, from the images without
a focal-plane filter, against the leak.

If the leak is visible light that passes the entrance filters, of which the
focal-plane filter passes a fixed part, the leak grows as the images
without a focal-plane filter brighten: the ratio of the leak in each eclipse
to the brightness of the images at the time is constant.

Both are compared in a ring near the limb, :data:`ring`, rather than at
disk center, which in the images has passed the range over which the
detector is linear: the short images reach 11,000 DN at disk center in
2024 and 14,000 DN in 2026, and the detector is linear to 1% up to 11,000 DN
(Boerner et al. 2012), while the ring stays below 9,000 DN.
"""

import dataclasses
import datetime
import functools

import astropy.time
import astropy.units as u
import numpy as np

import aia_94_leak_paper

__all__ = [
    "Comparison",
    "brightness",
    "compare",
    "dominance",
    "growth",
    "leak",
    "ring",
]

ring: tuple[float, float] = (0.85, 0.9)
"""The ring the images and the leak are compared in, in solar radii."""

dominance: int = 3
"""
How many times its mean before the onset the leak in the ring must be for
the growing part of it to dominate the glow that came before.
"""


def _median(image: np.ndarray, r: np.ndarray, lo: float, hi: float) -> float:
    """The median of an image between two distances from disk center."""
    return float(np.nanmedian(image[(r >= lo) & (r < hi)]))


def brightness(
    monitor: "aia_94_leak_paper.entrance.Monitor",
) -> u.Quantity:
    """
    The brightness of each image in :data:`ring`, the median of its bins.

    Parameters
    ----------
    monitor
        The images, as :func:`aia_94_leak_paper.entrance.load` gives them.
    """
    axes = ("time", "y", "x")
    unit = monitor.images.unit
    images = monitor.images.to_value(unit).ndarray_aligned(axes)
    r = monitor.r.ndarray_aligned(axes)
    return u.Quantity([_median(i, d, *ring) for i, d in zip(images, r)]) * unit


@functools.cache
def leak(date: datetime.date) -> u.Quantity:
    """
    The leak of one eclipse in :data:`ring`, with the whole disk up, through
    the point spread function, per second.

    Parameters
    ----------
    date
        A day of :data:`aia_94_leak_paper.eclipses.dates`.
    """
    eclipses = aia_94_leak_paper.eclipses
    axes = ("y", "x")
    seen = eclipses.seen(date)
    unit = seen.unit
    position = eclipses.load(date).position
    r = np.sqrt(np.square(position.x) + np.square(position.y))
    return (
        _median(
            seen.to_value(unit).ndarray_aligned(axes), r.ndarray_aligned(axes), *ring
        )
        * unit
    )


def growth(
    monitor: "aia_94_leak_paper.entrance.Monitor",
    since: None | int = None,
) -> "aia_94_leak_paper.eclipses.Growth":
    """
    Fit the brightness of the images with an exponential in time, as
    :func:`aia_94_leak_paper.eclipses.growth` fits the leak.

    Parameters
    ----------
    monitor
        The images, as :func:`aia_94_leak_paper.entrance.load` gives them.
    since
        The first year to fit, :data:`aia_94_leak_paper.eclipses.onset` if
        :obj:`None`.
    """
    if since is None:
        since = aia_94_leak_paper.eclipses.onset
    year = monitor.time.decimalyear
    value = brightness(monitor)
    where = (year >= since) & np.isfinite(value)
    slope, intercept = np.polyfit(year[where] - since, np.log(value[where].value), 1)
    residual = np.log(value[where].value) - (intercept + slope * (year[where] - since))
    return aia_94_leak_paper.eclipses.Growth(
        since=since,
        start=np.exp(intercept) * value.unit,
        doubling=np.log(2) / slope * u.yr,
        scatter=float(np.sqrt(np.mean(np.square(residual)))),
    )


@dataclasses.dataclass(frozen=True)
class Comparison:
    """The leak in each eclipse in :data:`ring` against the brightness of the images there at the time."""

    dates: list[datetime.date]
    """The days of the eclipses since the first image."""

    time: astropy.time.Time
    """When each eclipse began."""

    leak: u.Quantity
    """The leak in each eclipse."""

    brightness: u.Quantity
    """The brightness of the images, interpolated in its logarithm to each eclipse."""

    dominant: np.ndarray
    """
    Whether the leak in each eclipse is :data:`dominance` times its mean
    before the onset, so that its growing part dominates the glow that came
    before.
    """

    @property
    def ratio(self) -> np.ndarray:
        """The leak over the brightness of the images, the part the focal-plane filter passes."""
        return (self.leak / self.brightness).to_value(u.one)

    @property
    def ratio_mean(self) -> float:
        """The geometric mean of :attr:`ratio` where the leak is :attr:`dominant`."""
        return float(np.exp(np.mean(np.log(self.ratio[self.dominant]))))

    @property
    def scatter(self) -> float:
        """The root-mean-square scatter of the natural logarithm of :attr:`ratio` where the leak is :attr:`dominant`."""
        log = np.log(self.ratio[self.dominant])
        return float(np.sqrt(np.mean(np.square(log - log.mean()))))


def compare() -> Comparison:
    """The leak in each eclipse since the first image against the brightness of the images."""
    eclipses = aia_94_leak_paper.eclipses
    monitor = aia_94_leak_paper.entrance.load()
    value = brightness(monitor)
    good = np.isfinite(value)
    year = monitor.time.decimalyear[good]
    log = np.log(value[good].value)
    leaks = eclipses.history()
    before = u.Quantity(
        [leak(k.date) for k in leaks if k.date.year < eclipses.onset]
    ).mean()
    leaks = [k for k in leaks if k.time.decimalyear >= year[0]]
    at = np.array([k.time.decimalyear for k in leaks])
    inside = u.Quantity([leak(k.date) for k in leaks])
    return Comparison(
        dates=[k.date for k in leaks],
        time=astropy.time.Time([k.time for k in leaks]),
        leak=inside,
        brightness=np.exp(np.interp(at, year, log)) * value.unit,
        dominant=(inside > dominance * before),
    )

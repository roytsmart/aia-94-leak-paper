"""
A forward model of the leak through the ingress of an eclipse.

During ingress the atmosphere absorbs the EUV well before the visible disk
sets, so for some tens of seconds each part of the disk is seen in visible
light alone. Each bin is put on a clock started by its own EUV edge (when
171 Å falls through half of its value before the eclipse), and the 94 Å
frames in which the EUV is gone are fit with

.. math::

    I_{94} = \\mathrm{PSF} * \\left( A \\left[1 - u (1 - \\mu)\\right]
        \\tfrac{1}{2} \\mathrm{erfc}\\frac{s - s_\\mathrm{set}}{w} \\right)
        + a P_{94},

a limb-darkened disk of amplitude :math:`A` setting behind the atmosphere
at a time :math:`s_\\mathrm{set}` after the EUV edge, seen through the point
spread function of the leak, plus a fraction :math:`a` of the 94 Å image
before the edge, whatever EUV is left. The point spread function is a sharp
core and a broad Gaussian wing which carries a fraction :math:`\\eta` of the
light; it is fit once to the eclipses in which the leak is brightest and
then held fixed.
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
import scipy.optimize
import scipy.special

import aia_94_leak_paper

from ._archive import directory_default

__all__ = [
    "PSF",
    "Growth",
    "Leak",
    "dates_psf",
    "fit",
    "growth",
    "history",
    "onset",
    "psf",
    "seen",
]

dates_psf: tuple[datetime.date, ...] = (
    datetime.date(2025, 1, 23),
    datetime.date(2025, 7, 23),
    datetime.date(2026, 1, 19),
    datetime.date(2026, 7, 15),
)
"""The eclipses the point spread function of the leak is fit to."""

_window = (-40, -15)
"""The frames before the edge of a bin which give its image before the edge, in seconds."""

_gone = 39
"""How long after the EUV edge of a bin on the disk its EUV is gone, in seconds."""

_gone_corona = 70
"""The same for a bin off the limb, where the corona is fainter but higher."""

_extent = 1.4
"""How far from the center of the disk bins are fit, in solar radii."""

_downsample = 4
"""The factor the images are averaged down by before the wing is computed."""

_memory = joblib.Memory(
    directory_default / "joblib",
    verbose=0,
)

_code = hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest()
"""
A digest of this module, passed to every cached function, since `joblib`
keys only on the source of the function it caches and not on the functions
that one calls.
"""


@dataclasses.dataclass(frozen=True)
class PSF:
    """The point spread function of the leak: a sharp core and a Gaussian wing."""

    eta: float
    """The fraction of the light in the wing."""

    sigma: float
    """The standard deviation of the wing, in solar radii."""

    def __call__(self, images: np.ndarray, radius: float) -> np.ndarray:
        """
        Images seen through this point spread function.

        The wing is broad, so it is computed on images averaged down by
        :data:`_downsample` and repeated back up.

        Parameters
        ----------
        images
            A stack of images, with the last two axes ``y`` and ``x``.
        radius
            The radius of the disk, in bins.
        """
        n = _downsample
        *rest, ny, nx = images.shape
        small = images.reshape(*rest, ny // n, n, nx // n, n).mean(axis=(-3, -1))
        sigma = self.sigma * radius / n
        wing = scipy.ndimage.gaussian_filter(
            small,
            sigma=(0,) * len(rest) + (sigma, sigma),
            mode="constant",
            truncate=4,
        )
        wing = np.repeat(np.repeat(wing, n, axis=-2), n, axis=-1)
        return (1 - self.eta) * images + self.eta * wing


@dataclasses.dataclass
class _Ingress:
    """The arrays of one ingress the model is fit to, aligned as ``(time, y, x)``."""

    date: datetime.date
    s: np.ndarray
    """The time of each frame after the EUV edge of each bin, in seconds."""
    data: np.ndarray
    """94 Å, less the dark, scaled to the median exposure, in DN."""
    before: np.ndarray
    """94 Å before the EUV edge of each bin, less the dark, in DN."""
    r: np.ndarray
    """The distance of each bin from the center of the disk, in solar radii."""
    use: np.ndarray
    """The samples the model is fit to."""
    radius: float
    """The radius of the disk, in bins."""
    exposure: u.Quantity
    """The median exposure time of the 94 Å frames."""
    time: astropy.time.Time
    """When the 171 Å disk mean began to fall."""


def _seconds(times: na.ScalarArray, start: astropy.time.Time) -> na.ScalarArray:
    """Times after a start, in seconds, as a named array of floats."""
    return na.ScalarArray((times.ndarray - start).to_value(u.s), axes=times.axes)


def _ingress(eclipse: "aia_94_leak_paper.eclipses.Eclipse") -> _Ingress:
    """Put the 94 Å frames of an ingress on the clock of the EUV edge of each bin."""
    images, exposures = eclipse.images, eclipse.exposures

    def rate(key: str) -> na.ScalarArray:
        return images[key] / exposures[key]

    # when 171 falls through half its value before the eclipse, in each bin;
    # the corners beyond the field stop see nothing, and divide by zero
    dark171 = rate("UMBRA171").mean("time")
    t171 = _seconds(eclipse.times["INGRESS171"], eclipse.start)
    with np.errstate(divide="ignore", invalid="ignore"):
        transmission = (rate("INGRESS171") - dark171) / (
            rate("PRE171").mean("time") - dark171
        )
        below = transmission < 0.5
        after = np.argmax(below, axis="time")
        before = after | {"time": np.maximum(after["time"] - 1, 0)}
        va, vb = transmission[before], transmission[after]
        ta, tb = t171[{"time": before["time"]}], t171[{"time": after["time"]}]
        edge = ta + (va - 0.5) / (va - vb) * (tb - ta)
    crossed = np.any(below, axis="time") & (after["time"] > 0)

    # the edge is straight across the disk, so fit it with a plane
    position = eclipse.position
    r = np.sqrt(np.square(position.x) + np.square(position.y))
    axes = ("y", "x")
    x, y = position.x.broadcast_to(r.shape), position.y.broadcast_to(r.shape)
    on_disk = (crossed & (r < 0.95) & np.isfinite(edge)).ndarray_aligned(axes)
    design = np.stack(
        [
            np.ones(on_disk.sum()),
            x.ndarray_aligned(axes)[on_disk],
            y.ndarray_aligned(axes)[on_disk],
        ],
        axis=-1,
    )
    (t0, gx, gy), *_ = np.linalg.lstsq(
        design, edge.ndarray_aligned(axes)[on_disk], rcond=None
    )
    plane = t0 + gx * x + gy * y

    # 94 on that clock, less the dark and scaled to the median exposure
    exposure = np.median(exposures["INGRESS94"].ndarray)
    t94 = _seconds(eclipse.times["INGRESS94"], eclipse.start)
    s = t94 - plane
    dark = (images["UMBRA94"] * exposure / exposures["UMBRA94"]).mean("time")
    data = images["INGRESS94"] * exposure / exposures["INGRESS94"] - dark

    # each bin before its edge, from its own frames where there are any
    early = (s >= _window[0]) & (s < _window[1])
    count = np.sum(early, axis="time")
    pre = (images["PRE94"] * exposure / exposures["PRE94"]).mean("time") - dark
    before94 = np.where(
        count > 0,
        np.sum(np.where(early, data, 0 * u.DN), axis="time") / np.maximum(count, 1),
        pre,
    )

    # the bright cores of active regions are left out, where the EUV left matters most
    disk = r < 1
    core = disk & (
        before94
        > 4 * np.median(before94.ndarray_aligned(axes)[(r < 0.9).ndarray_aligned(axes)])
    )
    use = ((s >= _gone) & disk) | ((s >= _gone_corona) & ~disk & (r < _extent))
    use = use & ~core

    axes3 = ("time", "y", "x")
    return _Ingress(
        date=eclipse.date,
        s=s.ndarray_aligned(axes3),
        data=data.to(u.DN).ndarray_aligned(axes3).value,
        before=before94.to(u.DN).ndarray_aligned(axes).value,
        r=r.ndarray_aligned(axes),
        use=use.broadcast_to(s.shape).ndarray_aligned(axes3),
        radius=eclipse.radius,
        exposure=exposure,
        time=eclipse.start,
    )


def _model(q: np.ndarray, ingress: _Ingress, psf: PSF) -> np.ndarray:
    """The model of every frame of an ingress, for parameters ``q = (A, u, s_set, w, a)``."""
    amplitude, u_limb, s_set, w, a = q
    mu = np.sqrt(np.clip(1 - np.square(ingress.r), 0, 1))
    disk = np.where(ingress.r < 1, amplitude * (1 - u_limb * (1 - mu)), 0)
    images = disk * 0.5 * scipy.special.erfc((ingress.s - s_set) / w)
    return psf(images, ingress.radius) + a * ingress.before


_q0 = np.array([1.0, 0.55, 75.0, 5.0, 0.01])
_lower = np.array([0, -0.5, 40, 0.5, -0.2])
_upper = np.array([50, 1.5, 130, 30, 0.5])
_scale = np.array([0.5, 0.1, 5, 2, 0.01])
_loss = {"loss": "soft_l1", "f_scale": 0.3}
"""Robust against the particle hits the despiker missed, which reach hundreds of DN."""


@_memory.cache
def _psf(
    dates: tuple[datetime.date, ...], version: str, code: str
) -> tuple[float, float]:
    ingresses = [_ingress(aia_94_leak_paper.eclipses.load(d)) for d in dates]

    def residual(p: np.ndarray) -> np.ndarray:
        psf = PSF(eta=p[0], sigma=p[1])
        out = []
        for i, ingress in enumerate(ingresses):
            q = p[2 + 5 * i : 7 + 5 * i]
            out.append(_model(q, ingress, psf)[ingress.use] - ingress.data[ingress.use])
        return np.concatenate(out)

    p0 = np.concatenate([[0.5, 1.0]] + [_q0 * [2, 1, 1, 1, 1]] * len(ingresses))
    lower = np.concatenate([[0, 0.2]] + [_lower] * len(ingresses))
    upper = np.concatenate([[0.95, 5]] + [_upper] * len(ingresses))
    scale = np.concatenate([[0.05, 0.1]] + [_scale] * len(ingresses))
    result = scipy.optimize.least_squares(
        residual, p0, bounds=(lower, upper), x_scale=scale, **_loss
    )
    return float(result.x[0]), float(result.x[1])


def psf() -> PSF:
    """The point spread function of the leak, fit jointly to the eclipses of :data:`dates_psf`."""
    eta, sigma = _psf(dates_psf, aia_94_leak_paper.eclipses.version, _code)
    return PSF(eta=eta, sigma=sigma)


@dataclasses.dataclass(frozen=True)
class Leak:
    """The leak measured in one eclipse, in the units of the detector."""

    date: datetime.date
    """The day of the eclipse."""

    time: astropy.time.Time
    """When the 171 Å disk mean began to fall."""

    amplitude: u.Quantity
    """The leak at disk center before the point spread function, per second."""

    limb_darkening: float
    """The coefficient :math:`u` of :math:`1 - u (1 - \\mu)`."""

    sunset: u.Quantity
    """How long after the EUV edge of a bin the visible disk is half set there."""

    sunset_width: u.Quantity
    """The width of the setting, :math:`w`."""

    euv: float
    """The fraction of the 94 Å image before the edge that is left, :math:`a`."""

    center: u.Quantity
    """The leak at disk center with the whole disk up, wing included, per second."""

    glow: u.Quantity
    """The leak just off the limb (1.05 to 1.15 solar radii) with the whole disk up, per second."""

    quiet: u.Quantity
    """The median 94 Å signal within a third of a radius of disk center before the eclipse, per second."""

    rms: u.Quantity
    """The root-mean-square residual of the fit, per bin and frame."""


@_memory.cache
def _fit(
    date: datetime.date,
    eta: float,
    sigma: float,
    setting: None | tuple[float, float],
    version: str,
    code: str,
) -> dict:
    ingress = _ingress(aia_94_leak_paper.eclipses.load(date))
    psf = PSF(eta=eta, sigma=sigma)
    use = ingress.use

    # the time and width of the setting are fit, or held at `setting`
    free = np.array([True, True, setting is None, setting is None, True])
    fixed = _q0.copy()
    if setting is not None:
        fixed[2:4] = setting

    def full(q: np.ndarray) -> np.ndarray:
        result = fixed.copy()
        result[free] = q
        return result

    def residual(q: np.ndarray) -> np.ndarray:
        return _model(full(q), ingress, psf)[use] - ingress.data[use]

    result = scipy.optimize.least_squares(
        residual,
        _q0[free],
        bounds=(_lower[free], _upper[free]),
        x_scale=_scale[free],
        **_loss,
    )
    amplitude, u_limb, s_set, w, a = full(result.x)
    r = ingress.r
    mu = np.sqrt(np.clip(1 - np.square(r), 0, 1))
    seen = psf(np.where(r < 1, amplitude * (1 - u_limb * (1 - mu)), 0), ingress.radius)
    per_second = 1 / ingress.exposure.to_value(u.s)
    return {
        "time": ingress.time.isot,
        "amplitude": amplitude * per_second,
        "limb_darkening": u_limb,
        "sunset": s_set,
        "sunset_width": w,
        "euv": a,
        "center": np.median(seen[r < 0.1]) * per_second,
        "glow": np.median(seen[(r > 1.05) & (r < 1.15)]) * per_second,
        "quiet": np.median(ingress.before[r < 1 / 3]) * per_second,
        "rms": float(np.sqrt(np.mean(np.square(result.fun)))),
    }


onset: int = 2016
"""
The year from which the leak is fit with an exponential.

Before it the leak is a faint diffuse glow, which lies above the
exponential the later eclipses follow, and too faint to time its setting.
It is also the first year with an eclipse in each season.
"""


def _setting(p: PSF) -> tuple[float, float]:
    """The median time and width of the setting in the eclipses since :data:`onset`."""
    fits = [
        _fit(d, p.eta, p.sigma, None, aia_94_leak_paper.eclipses.version, _code)
        for d in aia_94_leak_paper.eclipses.dates
        if d.year >= onset
    ]
    return (
        float(np.median([f["sunset"] for f in fits])),
        float(np.median([f["sunset_width"] for f in fits])),
    )


def fit(date: datetime.date) -> Leak:
    """
    The leak in the ingress of one eclipse, with the point spread function of :func:`psf`.

    Before :data:`onset` the leak is too faint to time its setting, so the
    time and width of the setting are held at their medians since.

    Parameters
    ----------
    date
        A day of :data:`aia_94_leak_paper.eclipses.dates`.
    """
    p = psf()
    setting = None if date.year >= onset else _setting(p)
    result = _fit(
        date, p.eta, p.sigma, setting, aia_94_leak_paper.eclipses.version, _code
    )
    rate = u.DN / u.s
    return Leak(
        date=date,
        time=astropy.time.Time(result["time"], scale="utc"),
        amplitude=result["amplitude"] * rate,
        limb_darkening=result["limb_darkening"],
        sunset=result["sunset"] * u.s,
        sunset_width=result["sunset_width"] * u.s,
        euv=result["euv"],
        center=result["center"] * rate,
        glow=result["glow"] * rate,
        quiet=result["quiet"] * rate,
        rms=result["rms"] * u.DN,
    )


def history() -> list[Leak]:
    """The leak in every eclipse of :data:`aia_94_leak_paper.eclipses.dates`."""
    return [fit(date) for date in aia_94_leak_paper.eclipses.dates]


def seen(date: datetime.date) -> na.ScalarArray:
    """
    The leak of one eclipse with the whole disk up, through the point spread
    function, per second, on the bins of the archive.

    Parameters
    ----------
    date
        A day of :data:`aia_94_leak_paper.eclipses.dates`.
    """
    eclipse = aia_94_leak_paper.eclipses.load(date)
    leak = fit(date)
    position = eclipse.position
    r = np.sqrt(np.square(position.x) + np.square(position.y))
    mu = np.sqrt(np.clip(1 - np.square(r), 0, 1))
    disk = np.where(r < 1, 1 - leak.limb_darkening * (1 - mu), 0)
    axes = ("y", "x")
    image = psf()(disk.ndarray_aligned(axes), eclipse.radius)
    return na.ScalarArray(image, axes=axes) * leak.amplitude


@dataclasses.dataclass(frozen=True)
class Growth:
    """An exponential fit to the leak at disk center against time."""

    since: int
    """The first year of the eclipses fit."""

    start: u.Quantity
    """The leak at disk center at the start of :attr:`since`."""

    doubling: u.Quantity
    """The time the leak takes to double."""

    scatter: float
    """The root-mean-square scatter of the natural logarithm of the leak about the fit."""

    def __call__(self, year: np.ndarray) -> u.Quantity:
        """The fit at some decimal years."""
        return self.start * 2 ** ((year - self.since) / self.doubling.to_value(u.yr))


def growth(leaks: list[Leak], since: int = onset) -> Growth:
    """
    Fit the leak at disk center with an exponential in time.

    Parameters
    ----------
    leaks
        The leak in each eclipse, as :func:`history` gives it.
    since
        The first year to fit.
    """
    year = np.array([leak.time.decimalyear for leak in leaks])
    center = u.Quantity([leak.center for leak in leaks])
    where = year >= since
    slope, intercept = np.polyfit(year[where] - since, np.log(center[where].value), 1)
    residual = np.log(center[where].value) - (intercept + slope * (year[where] - since))
    return Growth(
        since=since,
        start=np.exp(intercept) * center.unit,
        doubling=np.log(2) / slope * u.yr,
        scatter=float(np.sqrt(np.mean(np.square(residual)))),
    )

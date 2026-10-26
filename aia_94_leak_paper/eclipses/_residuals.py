"""
The fit of :mod:`._model` against the frames it is made to.

It is kept apart from the model, whose cached fits are keyed on a digest of
its module.
"""

import dataclasses
import datetime

import astropy.units as u
import named_arrays as na
import numpy as np

import aia_94_leak_paper

from ._model import _ingress, _model, fit, psf

__all__ = [
    "Residuals",
    "before_setting",
    "residuals",
]

before_setting: u.Quantity = 10 * u.s
"""How long before the visible disk sets in a bin the plateau of that bin ends."""

_gradient = 0.9
"""How far from the center of the disk the plateau is fit with a plane, in solar radii."""


@dataclasses.dataclass(frozen=True)
class Residuals:
    """
    The 94 Å frames of one ingress against the fit to them, on the clock of
    the EUV edge of each bin, per second, in the units of
    :class:`aia_94_leak_paper.eclipses.Leak`.
    """

    date: datetime.date
    """The day of the eclipse."""

    s: na.ScalarArray
    """The time of each frame after the EUV edge of each bin."""

    data: na.ScalarArray
    """The frames, less the dark."""

    model: na.ScalarArray
    """The fit to the frames."""

    use: na.ScalarArray
    """The samples the fit is made to."""

    position: na.Cartesian2dVectorArray
    """The center of each bin, in solar radii from the center of the disk."""

    setting: u.Quantity
    """How long after the EUV edge of a bin the visible disk is half set there."""

    @property
    def r(self) -> na.ScalarArray:
        """The distance of each bin from the center of the disk, in solar radii."""
        return np.sqrt(np.square(self.position.x) + np.square(self.position.y))

    @property
    def residual(self) -> na.ScalarArray:
        """The frames less the fit."""
        return self.data - self.model

    @property
    def plateau(self) -> na.ScalarArray:
        """
        The mean residual of each bin on the disk over its plateau, from when
        its EUV is gone until :data:`before_setting` before the visible disk
        sets there, and NaN in the bins with no frame there.
        """
        where = self.use & (self.r < 1) & (self.s < self.setting - before_setting)
        count = np.sum(where, axis="time")
        total = np.sum(np.where(where, self.residual, 0), axis="time")
        return np.where(count > 0, total / np.maximum(count, 1), np.nan)

    @property
    def gradient(self) -> na.Cartesian2dVectorArray:
        """
        The gradient of a plane fit to :attr:`plateau` within
        :data:`_gradient` of the center of the disk, per solar radius.
        """
        axes = ("y", "x")
        plateau = self.plateau
        unit = plateau.unit
        values = plateau.to_value(unit).ndarray_aligned(axes)
        x, y = (
            c.broadcast_to(self.r.shape).ndarray_aligned(axes)
            for c in (self.position.x, self.position.y)
        )
        where = np.isfinite(values) & (self.r < _gradient).ndarray_aligned(axes)
        design = np.stack([np.ones(where.sum()), x[where], y[where]], axis=-1)
        (_, gx, gy), *_ = np.linalg.lstsq(design, values[where], rcond=None)
        return na.Cartesian2dVectorArray(x=gx * unit, y=gy * unit)


def residuals(date: datetime.date) -> Residuals:
    """
    The frames of the ingress of one eclipse against the fit of
    :func:`aia_94_leak_paper.eclipses.fit`.

    Parameters
    ----------
    date
        A day of :data:`aia_94_leak_paper.eclipses.dates`.
    """
    eclipse = aia_94_leak_paper.eclipses.load(date)
    ingress = _ingress(eclipse)
    leak = fit(date)
    exposure = ingress.exposure.to(u.s)
    q = np.array(
        [
            (leak.amplitude * exposure).to_value(u.DN),
            leak.limb_darkening,
            leak.sunset.to_value(u.s),
            leak.sunset_width.to_value(u.s),
            leak.euv,
        ]
    )
    model = _model(q, ingress, psf())
    axes = ("time", "y", "x")
    return Residuals(
        date=date,
        s=na.ScalarArray(ingress.s * u.s, axes=axes),
        data=na.ScalarArray(ingress.data * u.DN / exposure, axes=axes),
        model=na.ScalarArray(model * u.DN / exposure, axes=axes),
        use=na.ScalarArray(ingress.use, axes=axes),
        position=eclipse.position,
        setting=leak.sunset,
    )

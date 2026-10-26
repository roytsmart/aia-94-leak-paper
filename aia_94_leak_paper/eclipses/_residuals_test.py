import astropy.units as u
import numpy as np
import pytest

import aia_94_leak_paper

eclipses = aia_94_leak_paper.eclipses


@pytest.mark.parametrize(
    "date",
    [d for d in eclipses.dates if d.year >= eclipses.onset],
    ids=str,
)
def test_residuals(date):
    result = eclipses.residuals(date)
    assert result.date == date
    center = eclipses.fit(date).center
    residual = result.residual.ndarray[result.use.ndarray]
    assert abs(np.median(residual)) < 0.02 * center
    # the plateau is brighter than the fit toward the right of the detector,
    # whichever way the edge sweeps the disk
    gradient = result.gradient
    assert gradient.x > 0 * u.DN / u.s
    assert abs(gradient.y) < gradient.x

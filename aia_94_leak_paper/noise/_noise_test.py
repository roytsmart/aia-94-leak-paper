import astropy.units as u
import numpy as np
import pytest

import aia_94_leak_paper

noise = aia_94_leak_paper.noise
dates = aia_94_leak_paper.eclipses.dates


@pytest.mark.parametrize("date", dates, ids=str)
def test_load(date):
    run = noise.load(date)
    assert run.date == date
    assert run.count.axes == ("signal", "y", "x")
    assert run.count.sum().ndarray > 1e8
    assert run.exposure > 0 * u.s


@pytest.mark.parametrize("date", dates, ids=str)
def test_fit(date):
    result = noise.fit(date)
    assert result.read.shape == result.slope.shape == (4,)
    assert np.all((1 * u.DN**2 < result.read) & (result.read < 2 * u.DN**2))
    assert 0.3 < np.nanmean(result.slope) < 1
    assert np.isfinite(result.noiseless)
    if date.year >= aia_94_leak_paper.eclipses.onset:
        assert result.noiseless > 0 * u.DN

import astropy.units as u
import numpy as np
import pytest

import aia_94_leak_paper

eclipses = aia_94_leak_paper.eclipses


def test_psf_keeps_the_light():
    """Light near the middle of the frame is spread out, not lost or made."""
    images = np.zeros((2, 256, 256))
    images[:, 120:136, 120:136] = 1
    psf = eclipses.PSF(eta=0.5, sigma=0.2)
    result = psf(images, radius=100)
    assert result.shape == images.shape
    assert np.allclose(result.sum(axis=(-2, -1)), images.sum(axis=(-2, -1)), rtol=1e-3)


@pytest.mark.parametrize("date", eclipses.dates, ids=str)
def test_fit(date):
    result = eclipses.fit(date)
    assert result.date == date
    assert 0 * u.DN / u.s < result.glow < result.center < result.quiet
    assert 40 * u.s < result.sunset < 130 * u.s


def test_growth():
    result = eclipses.growth(eclipses.history())
    assert 2 * u.yr < result.doubling < 5 * u.yr
    assert result.scatter < 0.1

import astropy.units as u
import numpy as np
import pytest

import aia_94_leak_paper

dems = aia_94_leak_paper.dems
dates = aia_94_leak_paper.eclipses.dates


@pytest.mark.parametrize("date", dates, ids=str)
def test_load(date):
    synoptic = dems.load(date)
    assert synoptic.date == date
    assert synoptic.images.axes == ("wavelength", "y", "x")
    assert np.all(synoptic.exposure > 0 * u.s)
    # the flat field is near one on the disk
    flat = synoptic.flat.ndarray[synoptic.r.ndarray_aligned(("y", "x")) < 1]
    assert 0.8 < np.median(flat) < 1.2


def test_leak():
    date = dems.date_dem
    leak = dems.leak(date).to_value(u.DN / u.s).ndarray
    n = leak.shape[0]
    center = leak[n // 2 - 2 : n // 2 + 2, n // 2 - 2 : n // 2 + 2].mean()
    expected = aia_94_leak_paper.eclipses.fit(date).center.to_value(u.DN / u.s)
    assert center == pytest.approx(expected, rel=0.2)


@pytest.mark.parametrize("date", dates, ids=str)
def test_summary(date):
    result = dems.summary(date)
    # the leak only ever adds hot plasma, and the quiet Sun without it is cool
    assert result.quiet_corrected <= result.quiet_observed
    assert result.quiet_corrected < 0.1
    if np.isfinite(result.holes_observed):
        assert result.holes_corrected <= result.holes_observed
        assert result.holes_corrected < 0.2

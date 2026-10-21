import astropy.units as u
import numpy as np
import pytest

import aia_94_leak_paper

eclipses = aia_94_leak_paper.eclipses


@pytest.mark.parametrize("date", eclipses.dates, ids=str)
def test_load(date):
    result = eclipses.load(date)
    assert isinstance(result, eclipses.Eclipse)
    assert result.date == date

    groups = [f"{w}{c}" for c in eclipses.channels for w in ("PRE", "INGRESS", "UMBRA")]
    assert set(result.images) == set(groups)
    for group in groups:
        image = result.images[group]
        assert image.axes == ("time", "y", "x")
        assert image.unit == u.DN
        assert result.times[group].shape == {"time": image.shape["time"]}
        assert np.all(result.exposures[group] > 0 * u.s)
        t = result.times[group].ndarray
        assert np.all(np.diff(t.unix) > 0)

    # the windows come in order, and the ingress spans the fall of 171
    t = {k: v.ndarray for k, v in result.times.items()}
    assert t["PRE94"].max() < result.start < t["INGRESS94"].max()
    assert t["INGRESS94"].min() < result.gone < t["UMBRA94"].min()

    # the disk is near the middle of the frame, and about 100 bins across
    position = result.position
    r = np.sqrt(np.square(position.x) + np.square(position.y))
    assert r.min().ndarray < 0.02
    assert 80 < result.radius < 120

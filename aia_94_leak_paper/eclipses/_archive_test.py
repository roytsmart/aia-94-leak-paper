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

    # the windows come in order, the ingress spans the fall of 171, and the
    # dark follows it within minutes
    t = {k: v.ndarray for k, v in result.times.items()}
    assert t["PRE94"].max() < result.start < t["INGRESS94"].max()
    assert t["INGRESS94"].min() < result.gone < t["INGRESS94"].max()
    assert t["INGRESS94"].max() < t["UMBRA94"].min()
    assert (t["UMBRA94"].max() - result.gone).to_value(u.min) < 6

    # the disk is near the middle of the frame, and about 100 bins across
    position = result.position
    r = np.sqrt(np.square(position.x) + np.square(position.y))
    assert r.min().ndarray < 0.02
    assert 80 < result.radius < 120

    # the flat field is near one on the disk
    axes = ("y", "x")
    assert result.flat.axes == axes
    flat = result.flat.ndarray_aligned(axes)[(r < 1).ndarray_aligned(axes)]
    assert 0.8 < np.median(flat) < 1.2
    assert 0.8 < result.flat_center < 1.2

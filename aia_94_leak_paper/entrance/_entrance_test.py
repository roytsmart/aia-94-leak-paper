import astropy.units as u
import numpy as np

import aia_94_leak_paper

entrance = aia_94_leak_paper.entrance


def test_load():
    monitor = entrance.load()
    assert monitor.images.axes == ("time", "y", "x")
    # about one image a month since the first set, at the end of 2011
    assert monitor.images.shape["time"] > 150
    assert np.all(np.diff(monitor.time.mjd) > 0)
    assert np.all(monitor.exposure < 0.006 * u.s)


def test_brightness():
    monitor = entrance.load()
    brightness = entrance.brightness(monitor)
    assert np.all(np.isfinite(brightness))
    # the entrance filters pass more and more visible light
    assert brightness[-1] > 10 * brightness[0]
    growth = entrance.growth(monitor)
    assert 2 * u.yr < growth.doubling < 5 * u.yr


def test_compare():
    result = entrance.compare()
    assert result.dominant[-1]
    assert not result.dominant[0]
    # once the leak dominates, it is a constant part of the light
    assert result.scatter < 0.1
    assert 1e-7 < result.ratio_mean < 1e-6
    # before, it was more
    assert result.ratio[0] > result.ratio_mean

import astropy.units as u
import pytest

import aia_94_leak_paper

filters = aia_94_leak_paper.filters


@pytest.mark.parametrize("name", filters.filters)
def test_load(name):
    run = filters.load(name)
    assert run.date == filters.switch.date()
    assert run.count.axes == ("signal", "y", "x")
    assert run.count.sum().ndarray > 1e8
    assert run.exposure > 0 * u.s


def test_attenuation():
    result = filters.attenuation()
    # the EUV, measured from the corona, is attenuated as the calibration expects
    assert result.euv == pytest.approx(filters.euv_table, rel=0.1)
    # the leak about as much: visible light through the extra zirconium
    # would be all but gone
    assert 0.5 < result.leak < 1
    assert result.leak == pytest.approx(result.euv, rel=0.2)
    # the noise and the signals alone agree
    assert result.leak_signals == pytest.approx(result.leak, abs=0.1)
    assert result.noiseless_thick < result.noiseless_thin

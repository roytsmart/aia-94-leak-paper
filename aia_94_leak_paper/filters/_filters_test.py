import astropy.units as u
import named_arrays as na
import numpy as np
import optika
import pytest

import aia_94_leak_paper

filters = aia_94_leak_paper.filters


def test_run():
    run = filters.run()
    assert run.date == filters.switch.date()
    assert run.count.axes == ("signal", "y", "x")
    assert run.count.sum().ndarray > 1e8
    assert run.exposure > 0 * u.s


@pytest.mark.parametrize("name", filters.frames)
def test_frame(name):
    frame = filters.frame(name)
    run = filters.run()
    assert frame.time < run.time
    # equal to the millisecond, as the frames of a run are
    assert frame.exposure.to_value(u.s) == pytest.approx(
        run.exposure.to_value(u.s), abs=1e-3
    )
    assert frame.image.axes == ("y", "x")
    assert frame.image.shape == {a: run.count.shape[a] for a in ("y", "x")}
    expected = {"clean": "304_THIN", "banded": "304_THICK"}
    assert frame.followed == expected[name]


def test_band():
    # the band is about as large as the leak through the thin filter, and
    # leaves the top and bottom rows of the detector alone
    peak = filters.band_peak()
    assert peak > 0.5 * filters.attenuation().noiseless_center
    band = filters.band().to_value(u.DN).ndarray_aligned(("y", "x"))
    rows = np.median(band, axis=1)
    assert np.abs(rows[:4]).max() < peak.value / 10
    assert np.abs(rows[-4:]).max() < peak.value / 10


def test_attenuation():
    result = filters.attenuation()
    # the EUV, measured from the corona, is attenuated as the calibration expects
    assert result.euv == pytest.approx(filters.euv_table, rel=0.1)
    # the leak far more: the thick filter all but stops it
    assert abs(result.leak) < 3 * result.leak_spread
    assert result.leak < result.euv / 4
    assert result.noiseless_center > 0 * u.DN


def test_metal():
    # the extra zirconium attenuates the EUV as the calibration expects
    euv = float(filters.metal(filters.wavelength_euv).ndarray)
    assert euv == pytest.approx(filters.euv_table, rel=0.02)
    # and would all but stop visible light that crossed it, by every
    # measurement of the constants, passing no less than the thick filter
    # passes of the leak
    result = filters.attenuation()
    for table in filters.tables:
        most = filters.metal_visible(table)
        assert 0 < most < result.euv / 4
        assert result.leak < most + 3 * result.leak_spread
        for end in filters.visible:
            assert most >= float(filters.metal(end, table).ndarray)
    # the default of optika is the first of them
    assert filters.metal_visible() == filters.metal_visible(filters.tables[0])


@pytest.mark.parametrize("table", filters.tables)
def test_metal_thin_film(table):
    # the reflections at the surfaces of the filters, their oxide, and the
    # mesh they are supported by are left out of the ratio, since they are
    # the same for both
    wavelength = na.geomspace(*filters.visible, axis="wavelength", num=11)
    rays = optika.rays.RayVectorArray(
        wavelength=wavelength,
        direction=na.Cartesian3dVectorArray(0, 0, 1),
    )
    efficiency = {
        name: optika.materials.ThinFilmFilter(
            layer=optika.materials.Layer(
                chemical=optika.chemicals.Chemical("Zr", table=table),
                thickness=thickness,
            ),
            layer_oxide=optika.materials.Layer(chemical="ZrO2", thickness=5 * u.nm),
            mesh=optika.materials.meshes.Mesh(
                chemical="Ni",
                efficiency=1,
                pitch=1 / u.mm,
            ),
        ).efficiency(rays=rays, normal=na.Cartesian3dVectorArray(0, 0, -1))
        for name, thickness in filters.zirconium.items()
    }
    ratio = (efficiency["thick"] / efficiency["thin"]).ndarray
    expected = filters.metal(wavelength, table).ndarray
    assert ratio == pytest.approx(expected, rel=0.01)

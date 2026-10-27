import astropy.units as u
import named_arrays as na
import optika
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


def test_metal():
    # the extra zirconium attenuates the EUV as the calibration expects
    euv = float(filters.metal(filters.wavelength_euv).ndarray)
    assert euv == pytest.approx(filters.euv_table, rel=0.02)
    # and would attenuate visible light that crossed it far more than the
    # thick filter attenuates the leak, by every measurement of the constants
    leak = filters.attenuation().leak
    for table in filters.tables:
        most = filters.metal_visible(table)
        assert 0 < most < leak / 5
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

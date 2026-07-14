#!/usr/bin/env python
#
# test_galaxy_shape_factor.py
#
# Copyright (C) 2026 Sandro Dias Pinto Vitenti <vitenti@uel.br>
#
# numcosmo is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# numcosmo is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program. If not, see <http://www.gnu.org/licenses/>.

"""Unit tests for ``NcGalaxyShapeFactor``'s own (scheme-independent) API.

The per-scheme test files (``test_galaxy_shape_factor_var_add.py`` etc.)
exercise ``eval_marginal``/``eval_ln_marginal`` end-to-end; this file covers
what the base class itself owns and every scheme shares: the ``ellip-conv``
property, ``check_obs``, the ellipticity-transform wrappers
(``apply_shear``/``apply_shear_inv``/``lndet_jac``), and the plain
``NcGalaxyShapeFactorData`` get/set/radius accessors.
"""

import pytest
from numpy.testing import assert_allclose
from gi.repository import GLib

from numcosmo_py import Ncm, Nc
from numcosmo_py.helper import duplicate_via_serialization

Ncm.cfg_init()

_CONVS = [Nc.GalaxyWLObsEllipConv.TRACE, Nc.GalaxyWLObsEllipConv.TRACE_DET]


def _build_mset():
    cosmo = Nc.HICosmoDEXcdm.new()
    dist = Nc.Distance.new(100.0)
    hms = Nc.HaloCMParam.new(Nc.HaloMassSummaryMassDef.MEAN, 200.0)
    dp = Nc.HaloDensityProfileNFW.new(hms)
    hp = Nc.HaloPosition.new(dist)
    smd = Nc.WLSurfaceMassDensity.new(dist)
    pop = Nc.GalaxyShapePopGauss.new()

    hms.param_set_by_name("log10MDelta", 14.0)
    hp.param_set_by_name("z", 0.2)
    pop.param_set_by_name("sigma", 0.3)
    hp.prepare(cosmo)

    mset = Ncm.MSet.empty_new()
    for model in (cosmo, dp, hp, smd, pop):
        mset.set(model)
    mset.set(Nc.GalaxyRedshiftPopLSSTSRD.new_y1_source())
    mset.set(Nc.GalaxyRedshiftObsGauss.new())

    return mset


def _build_data(mset, ellip_conv):
    posf = Nc.GalaxyPositionFactorFlat.new(-0.2, 0.2, -0.2, 0.2)
    pos_data = Nc.GalaxyPositionFactorData.new(posf, mset)
    pos_data.ra = 0.03
    pos_data.dec = 0.02

    zf = Nc.GalaxyRedshiftFactorComposed.new(0.0, 20.0)
    z_data = Nc.GalaxyRedshiftFactorData.new(zf, mset)
    z_data.z = 0.9

    gsf = Nc.GalaxyShapeFactorVarAdd.new(ellip_conv)
    data = Nc.GalaxyShapeFactorData.new(gsf, mset, pos_data, z_data)

    return gsf, data


def test_serialize_deserialize():
    """A round trip through NcmSerialize preserves the "ellip-conv"
    property (checked on a concrete scheme, since the base class itself is
    abstract)."""
    gsf = Nc.GalaxyShapeFactorVarAdd.new(Nc.GalaxyWLObsEllipConv.TRACE_DET)

    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
    gsf2 = duplicate_via_serialization(gsf, ser)

    assert isinstance(gsf2, Nc.GalaxyShapeFactorVarAdd)
    assert gsf2 is not gsf
    assert gsf2.get_ellip_conv() == gsf.get_ellip_conv()


def test_integrand_copy_matches_original():
    """A copied integrand evaluates identically to the original."""
    mset = _build_mset()
    gsf, data = _build_data(mset, Nc.GalaxyWLObsEllipConv.TRACE)
    gsf.data_set(data, 0.05, -0.02, 0.03, 0.005, -0.003, 0.05, Nc.WLEllipticityFrame.CELESTIAL)
    gsf.prepare(mset)
    gsf.prepare_data_array(mset, [data], True, True)

    integ = gsf.integ(mset, False)
    integ_copy = integ.copy()

    assert integ_copy.eval(0.5, data) == integ.eval(0.5, data)


@pytest.mark.parametrize("ellip_conv", _CONVS)
def test_ellip_conv_property_round_trip(ellip_conv):
    """The "ellip-conv" property and get_ellip_conv() agree with the value
    passed at construction."""
    gsf = Nc.GalaxyShapeFactorVarAdd.new(ellip_conv)

    assert gsf.get_ellip_conv() == ellip_conv
    assert gsf.props.ellip_conv == ellip_conv


def test_check_obs_mismatch_raises():
    """check_obs() raises when the observation's convention differs from
    the shape factor's own."""
    mset = _build_mset()
    gsf, data = _build_data(mset, Nc.GalaxyWLObsEllipConv.TRACE)

    cols = Nc.GalaxyShapeFactorData.required_columns(data)
    obs = Nc.GalaxyWLObs.new(
        Nc.GalaxyWLObsEllipConv.TRACE_DET, Nc.WLEllipticityFrame.CELESTIAL, 1, cols
    )

    with pytest.raises(GLib.Error, match="ellip_conv_mismatch|convention"):
        gsf.check_obs(obs)


def test_check_obs_match():
    """check_obs() returns True when the conventions agree."""
    mset = _build_mset()
    gsf, data = _build_data(mset, Nc.GalaxyWLObsEllipConv.TRACE)

    cols = Nc.GalaxyShapeFactorData.required_columns(data)
    obs = Nc.GalaxyWLObs.new(
        Nc.GalaxyWLObsEllipConv.TRACE, Nc.WLEllipticityFrame.CELESTIAL, 1, cols
    )

    assert gsf.check_obs(obs) is True


@pytest.mark.parametrize("ellip_conv", _CONVS)
def test_apply_shear_round_trip(ellip_conv):
    """apply_shear_inv() undoes apply_shear() for both ellipticity
    conventions (exercises both switch branches in each wrapper)."""
    gsf = Nc.GalaxyShapeFactorVarAdd.new(ellip_conv)

    g = Ncm.Complex.new()
    g.set(0.05, -0.02)
    e_intrinsic = Ncm.Complex.new()
    e_intrinsic.set(0.1, 0.03)

    e_obs = Ncm.Complex.new()
    gsf.apply_shear(g, e_intrinsic, e_obs)

    e_recovered = Ncm.Complex.new()
    gsf.apply_shear_inv(g, e_obs, e_recovered)

    assert_allclose(e_recovered.Re(), e_intrinsic.Re(), atol=1.0e-10)
    assert_allclose(e_recovered.Im(), e_intrinsic.Im(), atol=1.0e-10)


@pytest.mark.parametrize("ellip_conv", _CONVS)
def test_lndet_jac_finite(ellip_conv):
    """lndet_jac() returns a finite value for both ellipticity conventions
    (exercises both switch branches)."""
    gsf = Nc.GalaxyShapeFactorVarAdd.new(ellip_conv)

    g = Ncm.Complex.new()
    g.set(0.05, -0.02)
    e_obs = Ncm.Complex.new()
    e_obs.set(0.1, 0.03)

    lndet = gsf.lndet_jac(g, e_obs)

    assert lndet == lndet  # not NaN
    assert abs(lndet) < 1.0e10


def test_data_get_matches_data_set():
    """data_get() reads back exactly what data_set() wrote."""
    mset = _build_mset()
    gsf, data = _build_data(mset, Nc.GalaxyWLObsEllipConv.TRACE)

    gsf.data_set(data, 0.05, -0.02, 0.03, 0.005, -0.003, 0.05, Nc.WLEllipticityFrame.CELESTIAL)

    e1, e2, std_noise, c1, c2, m = gsf.data_get(data)

    assert (e1, e2, std_noise, c1, c2, m) == (0.05, -0.02, 0.03, 0.005, -0.003, 0.05)
    assert data.coord == Nc.WLEllipticityFrame.CELESTIAL


def test_data_get_radius():
    """NcGalaxyShapeFactorData.get_radius() returns the projected radius
    cached by prepare_data_array()."""
    mset = _build_mset()
    gsf, data = _build_data(mset, Nc.GalaxyWLObsEllipConv.TRACE)
    gsf.data_set(data, 0.05, -0.02, 0.03, 0.005, -0.003, 0.05, Nc.WLEllipticityFrame.CELESTIAL)

    gsf.prepare(mset)
    gsf.prepare_data_array(mset, [data], True, True)

    radius = data.get_radius()

    assert radius > 0.0

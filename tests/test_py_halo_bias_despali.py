#!/usr/bin/env python
#
# test_py_halo_bias_despali.py
#
# Wed Nov 16:32:00 12 2025
# Copyright  2025  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_py_halo_bias_despali.py
# Copyright (C) 2025 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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
# with this program.  If not, see <http://www.gnu.org/licenses/>.

"""Tests for Nc.HaloBiasDespali model."""

import pytest
import numpy as np

from numcosmo_py import Ncm, Nc

Ncm.cfg_init()


@pytest.fixture(name="mfunc", params=[Nc.MultiplicityFuncDespali])
def fixture_mass_function(request, dist, psf) -> Nc.HaloMassFunction:
    """Fixture for HaloMassFunction."""
    mulf = request.param.new()
    return Nc.HaloMassFunction.new(dist, psf, mulf)


@pytest.fixture(name="bias_despali")
def fixture_bias_despali(mfunc: Nc.HaloMassFunction) -> Nc.HaloBiasDespali:
    """Fixture for HaloBiasDespali."""
    return Nc.HaloBiasDespali.new(mfunc)


def test_halo_bias_despali_basic(bias_despali: Nc.HaloBiasDespali):
    """Basic instantiation and type checks."""
    assert isinstance(bias_despali, Nc.HaloBiasDespali)
    assert isinstance(bias_despali, Nc.HaloBias)


def test_halo_bias_despali_flags(bias_despali: Nc.HaloBiasDespali):
    """Test CMF and EO flags toggling."""
    # Defaults
    cmf_default = bias_despali.get_cmf()
    eo_default = bias_despali.get_eo()
    assert isinstance(cmf_default, bool)
    assert isinstance(eo_default, bool)

    # Toggle flags
    bias_despali.set_cmf(not cmf_default)
    bias_despali.set_eo(not eo_default)
    assert bias_despali.get_cmf() != cmf_default
    assert bias_despali.get_eo() != eo_default

    # Restore
    bias_despali.set_cmf(cmf_default)
    bias_despali.set_eo(eo_default)
    assert bias_despali.get_cmf() == cmf_default
    assert bias_despali.get_eo() == eo_default


def test_halo_bias_despali_delta_c_and_vir(
    bias_despali: Nc.HaloBiasDespali, cosmo: Nc.HICosmo
):
    """Test delta_c and delta_vir methods."""
    z_array = np.linspace(0.0, 2.0, 10)
    delta_c_values = [bias_despali.delta_c(cosmo, z) for z in z_array]
    delta_vir_values = [bias_despali.delta_vir(cosmo, z) for z in z_array]

    assert np.all(np.isfinite(delta_c_values))
    assert np.all(np.isfinite(delta_vir_values))
    assert np.all(np.array(delta_c_values) > 0.0)
    assert np.all(np.array(delta_vir_values) > 0.0)

    # sanity: they should vary smoothly with z
    assert all(np.diff(delta_c_values) > 0.0)
    assert all(np.diff(delta_vir_values) > 0.0)


def test_halo_bias_despali_new_full(mfunc: Nc.HaloMassFunction):
    """Test creation with explicit flags."""
    for eo, cmf in [(False, False), (True, False), (False, True), (True, True)]:
        bias_despali = Nc.HaloBiasDespali.new_full(mfunc, eo, cmf)
        assert bias_despali.get_eo() == eo
        assert bias_despali.get_cmf() == cmf


def test_halo_bias_despali_eval_virial(
    cosmo: Nc.HICosmo, dist: Nc.Distance, psf: Ncm.PowspecFilter
):
    """Test bias evaluation with virial mass definition."""
    mulf = Nc.MultiplicityFuncDespali.new()
    mulf.set_mdef(Nc.MultiplicityFuncMassDef.VIRIAL)
    mfunc = Nc.HaloMassFunction.new(dist, psf, mulf)

    # Test all flag combinations for virial
    for eo, cmf in [(False, False), (True, False), (False, True)]:
        bias_despali = Nc.HaloBiasDespali.new_full(mfunc, eo, cmf)
        bias_value = bias_despali.eval(cosmo, 0.5, 0.5)
        assert np.isfinite(bias_value)
        assert bias_value > 0.0


def test_halo_bias_despali_eval_mean(
    cosmo: Nc.HICosmo, dist: Nc.Distance, psf: Ncm.PowspecFilter
):
    """Test bias evaluation with mean mass definition."""
    mulf = Nc.MultiplicityFuncDespali.new()
    mulf.set_mdef(Nc.MultiplicityFuncMassDef.MEAN)
    mfunc = Nc.HaloMassFunction.new(dist, psf, mulf)

    # Test both eo flag values for mean
    for eo in [False, True]:
        bias_despali = Nc.HaloBiasDespali.new_full(mfunc, eo, False)
        bias_value = bias_despali.eval(cosmo, 0.5, 0.5)
        assert np.isfinite(bias_value)
        assert bias_value > 0.0


def test_halo_bias_despali_eval_critical(
    cosmo: Nc.HICosmo, dist: Nc.Distance, psf: Ncm.PowspecFilter
):
    """Test bias evaluation with critical mass definition."""
    mulf = Nc.MultiplicityFuncDespali.new()
    mulf.set_mdef(Nc.MultiplicityFuncMassDef.CRITICAL)
    mfunc = Nc.HaloMassFunction.new(dist, psf, mulf)

    # Test both eo flag values for critical
    for eo in [False, True]:
        bias_despali = Nc.HaloBiasDespali.new_full(mfunc, eo, False)
        bias_value = bias_despali.eval(cosmo, 0.5, 0.5)
        assert np.isfinite(bias_value)
        assert bias_value > 0.0


def test_halo_bias_despali_properties(mfunc: Nc.HaloMassFunction):
    """Test setting properties via GObject interface."""
    bias_despali = Nc.HaloBiasDespali.new(mfunc)

    # Test setting via properties
    bias_despali.props.eo = True
    bias_despali.props.cmf = True
    assert bias_despali.props.eo
    assert bias_despali.props.cmf

    bias_despali.props.eo = False
    bias_despali.props.cmf = False
    assert not bias_despali.props.eo
    assert not bias_despali.props.cmf


def test_serialization_deserialization(bias_despali: Nc.HaloBiasDespali):
    """Test serialization and deserialization."""
    for eo, cmf in [(False, False), (True, False), (False, True), (True, True)]:
        bias_despali.set_eo(eo)
        bias_despali.set_cmf(cmf)
        ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
        bias_despali_dup = ser.dup_obj(bias_despali)
        assert bias_despali_dup is not None
        assert isinstance(bias_despali_dup, Nc.HaloBiasDespali)
        assert bias_despali_dup is not bias_despali

        assert bias_despali_dup.get_eo() == bias_despali.get_eo()
        assert bias_despali_dup.get_cmf() == bias_despali.get_cmf()

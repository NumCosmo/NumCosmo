#!/usr/bin/env python
#
# test_py_halo_bias_tinker.py
#
# Wed Nov 16:32:00 12 2025
# Copyright  2025  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_py_halo_bias_tinker.py
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

"""Tests for Nc.HaloBiasTinker model."""

import pytest
import numpy as np

from numcosmo_py import Ncm, Nc

Ncm.cfg_init()


@pytest.fixture(name="mfunc", params=[Nc.MultiplicityFuncTinker])
def fixture_mass_function(request, dist, psf) -> Nc.HaloMassFunction:
    """Fixture for HaloMassFunction."""
    mulf = request.param.new()
    return Nc.HaloMassFunction.new(dist, psf, mulf)


@pytest.fixture(name="bias_tinker")
def fixture_bias_tinker(mfunc: Nc.HaloMassFunction) -> Nc.HaloBiasTinker:
    """Fixture for HaloBiasTinker."""
    return Nc.HaloBiasTinker.new(mfunc)


def test_halo_bias_tinker_basic(bias_tinker: Nc.HaloBiasTinker):
    """Basic instantiation and type checks."""
    assert isinstance(bias_tinker, Nc.HaloBiasTinker)
    assert isinstance(bias_tinker, Nc.HaloBias)


def test_halo_bias_tinker_parameters(bias_tinker: Nc.HaloBiasTinker):
    """Test parameter getters and setters."""
    # Test delta_c
    delta_c = bias_tinker.get_delta_c()
    assert isinstance(delta_c, float)
    assert delta_c >= 0.0

    bias_tinker.set_delta_c(1.686)
    assert bias_tinker.get_delta_c() == 1.686

    # Test B parameter
    B = bias_tinker.get_B()
    assert isinstance(B, float)
    assert B >= 0.0

    bias_tinker.set_B(0.183)
    assert bias_tinker.get_B() == 0.183

    # Test b parameter
    b = bias_tinker.get_b()
    assert isinstance(b, float)
    assert b >= 0.0

    bias_tinker.set_b(1.5)
    assert bias_tinker.get_b() == 1.5

    # Test c parameter
    c = bias_tinker.get_c()
    assert isinstance(c, float)
    assert c >= 0.0

    bias_tinker.set_c(2.4)
    assert bias_tinker.get_c() == 2.4


def test_halo_bias_tinker_new_full(mfunc: Nc.HaloMassFunction):
    """Test creation with explicit parameters."""
    delta_c, B, b, c = 1.686, 0.183, 1.5, 2.4
    bias_tinker = Nc.HaloBiasTinker.new_full(mfunc, delta_c, B, b, c)

    assert bias_tinker.get_delta_c() == delta_c
    assert bias_tinker.get_B() == B
    assert bias_tinker.get_b() == b
    assert bias_tinker.get_c() == c


def test_halo_bias_tinker_eval(bias_tinker: Nc.HaloBiasTinker, cosmo: Nc.HICosmo):
    """Test bias evaluation."""
    # Set parameters
    bias_tinker.set_delta_c(1.686)
    bias_tinker.set_B(0.183)
    bias_tinker.set_b(1.5)
    bias_tinker.set_c(2.4)

    # Test evaluation at different sigma and z values
    sigma_array = np.linspace(0.3, 1.5, 10)
    z_array = np.linspace(0.0, 2.0, 10)

    for sigma in sigma_array:
        bias_value = bias_tinker.eval(cosmo, sigma, 0.5)
        assert np.isfinite(bias_value)

    for z in z_array:
        bias_value = bias_tinker.eval(cosmo, 0.5, z)
        assert np.isfinite(bias_value)


def test_halo_bias_tinker_eval_different_mdef(
    cosmo: Nc.HICosmo, dist: Nc.Distance, psf: Ncm.PowspecFilter
):
    """Test bias evaluation with different mass definitions."""
    for mdef in [
        Nc.MultiplicityFuncMassDef.MEAN,
        Nc.MultiplicityFuncMassDef.CRITICAL,
    ]:
        mulf = Nc.MultiplicityFuncTinker.new()
        mulf.set_mdef(mdef)
        mfunc = Nc.HaloMassFunction.new(dist, psf, mulf)

        bias_tinker = Nc.HaloBiasTinker.new_full(mfunc, 1.686, 0.183, 1.5, 2.4)
        bias_value = bias_tinker.eval(cosmo, 0.5, 0.5)

        assert np.isfinite(bias_value)


def test_halo_bias_tinker_properties(mfunc: Nc.HaloMassFunction):
    """Test setting properties via GObject interface."""
    bias_tinker = Nc.HaloBiasTinker.new_full(mfunc, 1.686, 0.183, 1.5, 2.4)

    # Test reading via properties
    assert bias_tinker.props.critical_delta == 1.686
    assert bias_tinker.props.B == 0.183
    assert bias_tinker.props.b == 1.5
    assert bias_tinker.props.c == 2.4


def test_serialization_deserialization(bias_tinker: Nc.HaloBiasTinker):
    """Test serialization and deserialization."""
    test_params = [
        (1.686, 0.183, 1.5, 2.4),
        (1.5, 0.2, 1.6, 2.5),
        (1.8, 0.15, 1.4, 2.3),
    ]

    for delta_c, B, b, c in test_params:
        bias_tinker.set_delta_c(delta_c)
        bias_tinker.set_B(B)
        bias_tinker.set_b(b)
        bias_tinker.set_c(c)

        ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
        bias_tinker_dup = ser.dup_obj(bias_tinker)

        assert bias_tinker_dup is not None
        assert isinstance(bias_tinker_dup, Nc.HaloBiasTinker)
        assert bias_tinker_dup is not bias_tinker

        assert bias_tinker_dup.get_delta_c() == bias_tinker.get_delta_c()
        assert bias_tinker_dup.get_B() == bias_tinker.get_B()
        assert bias_tinker_dup.get_b() == bias_tinker.get_b()
        assert bias_tinker_dup.get_c() == bias_tinker.get_c()

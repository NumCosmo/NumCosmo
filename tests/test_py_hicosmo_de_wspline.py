#!/usr/bin/env python
#
# test_py_hicosmo_de_wspline.py
#
# Mon Sep 22 18:14:50 2025
# Copyright  2025  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_py_hicosmo_de_wspline.py
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

"""Unit tests for the NcmMSet class."""

import pytest
import numpy as np
from numpy.testing import assert_allclose
from numcosmo_py import Ncm, Nc

Ncm.cfg_init()


@pytest.fixture(name="cosmo_w", params=[3, 4, 6, 12])
def fixture_cosmo_w(request) -> Nc.HICosmoDEWSpline:
    """Returns a DEWSpline model with preset knots and values."""
    return Nc.HICosmoDEWSpline.new(nknots=request.param, z_f=2.0)


def test_wspline_knots(cosmo_w: Nc.HICosmoDEWSpline) -> None:
    """Test that knots are set and retrieved correctly."""

    alpha = cosmo_w.get_alpha()

    assert alpha.len() > 0

    alpha_array = alpha.dup_array()
    # Now we check if it is monotonically increasing
    assert np.all(np.diff(alpha_array) > 0)


def test_wspline_eval_w(cosmo_w: Nc.HICosmoDEWSpline) -> None:
    """Test evaluation of w(z) at knot points and in between."""
    alpha = cosmo_w.get_alpha()
    w_len = alpha.len()
    w_array = np.random.uniform(-1.5, -0.5, size=w_len).astype(float)
    z_test = np.expm1(alpha.dup_array())

    for z in z_test:
        w_test = cosmo_w.w_de(z)
        assert_allclose(w_test, -1.0)  # Default value is -1.0

    w_vec = Ncm.Vector.new_array(w_array)
    cosmo_w.props.w = w_vec

    for z, w in zip(z_test, w_array):
        w_test = cosmo_w.w_de(z)
        assert_allclose(w_test, w)


def test_wspline_mean_kappa(cosmo_w: Nc.HICosmoDEWSpline) -> None:
    """Test evaluation of mean kappa at knot points and in between."""
    alpha = cosmo_w.get_alpha()
    w_len = alpha.len()
    w_array = np.random.uniform(-1.5, -0.5, size=w_len).astype(float)

    kappa_mean = cosmo_w.mean_kappa()
    assert_allclose(kappa_mean, 0.0)  # Default value is 0.0

    w_vec = Ncm.Vector.new_array(w_array)
    cosmo_w.props.w = w_vec

    kappa_mean = cosmo_w.mean_kappa()
    assert isinstance(kappa_mean, float)
    assert not np.isnan(kappa_mean)
    assert not np.isinf(kappa_mean)
    assert kappa_mean != 0.0  # Should be different from default now


def test_wspline_mean_kappa_mset_func(cosmo_w: Nc.HICosmoDEWSpline) -> None:
    """Test evaluation of mean kappa at knot points and in between."""
    alpha = cosmo_w.get_alpha()
    w_len = alpha.len()
    w_array = np.random.uniform(-1.5, -0.5, size=w_len).astype(float)
    w_vec = Ncm.Vector.new_array(w_array)
    cosmo_w.props.w = w_vec

    kappa_mean = cosmo_w.mean_kappa()
    kappa_mean_func = Ncm.MSetFuncList.new("NcHICosmoDEWSpline:mean_kappa", None)

    mset = Ncm.MSet.new_array([cosmo_w])
    assert_allclose(kappa_mean, kappa_mean_func.eval0(mset))

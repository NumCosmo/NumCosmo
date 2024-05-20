#!/usr/bin/env python
#
# test_py_hipert_two_fluids.py
#
# Mon May 20 00:24:30 2024
# Copyright  2024  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_py_hipert_two_fluids.py
# Copyright (C) 2024 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

"""Tests on NcHIPertTwoFluids perturbations module."""

from itertools import product
import pytest
from numpy.testing import assert_allclose
import numpy as np

from numcosmo_py import Ncm, Nc

Ncm.cfg_init()


@pytest.fixture(name="two_fluids")
def fixture_two_fluids() -> Nc.HIPertTwoFluids:
    """Fixture for NcHIPertTwoFluids."""
    two_fluids = Nc.HIPertTwoFluids.new()

    return two_fluids


@pytest.fixture(name="cosmo_qgrw")
def fixture_cosmo_qgrw() -> Nc.HICosmo:
    """Fixture for NcCosmoQGrw."""
    cosmo = Nc.HICosmoQGRW()
    cosmo.props.w = 1.0e-5
    cosmo.props.Omegar = 1.0 * (1.0e-5)
    cosmo.props.Omegaw = 1.0 * (1.0 - 1.0e-5)
    cosmo.props.xb = 1.0e30

    return cosmo


def test_init(two_fluids):
    """Test NcHIPertTwoFluids initialization."""
    assert two_fluids is not None
    assert isinstance(two_fluids, Nc.HIPertTwoFluids)
    assert isinstance(two_fluids, Nc.HIPert)


def test_compute_full_spectrum(two_fluids, cosmo_qgrw):
    """Test NcHIPertTwoFluids compute_full_spectrum."""
    two_fluids.props.reltol = 1.0e-9

    def spec_params(Omegars=1.0e-5, w=1.0e-3, E0=1.0):
        cosmo_qgrw.props.w = w
        cosmo_qgrw.props.Omegar = E0 * Omegars
        cosmo_qgrw.props.Omegaw = E0 * (1.0 - Omegars)

        spec1 = two_fluids.compute_zeta_spectrum(
            cosmo_qgrw, 1, -cosmo_qgrw.abs_alpha(1.0e-14), -1.0, 1.0e-3, 1.0e8, 10
        )
        spec2 = two_fluids.compute_zeta_spectrum(
            cosmo_qgrw, 2, -cosmo_qgrw.abs_alpha(1.0e-14), -1.0, 1.0e-3, 1.0e8, 10
        )

        return spec1, spec2

    specs1 = []
    specs2 = []
    w_a = np.geomspace(1.0e-5, 1.0e-1, 6)
    for w in w_a:
        spec1, spec2 = spec_params(w=w)
        specs1.append(spec1)
        specs2.append(spec2)

    lnk_v = specs1[0].peek_xv()
    lnw_v = Ncm.Vector.new_array(np.log(w_a))
    zm = Ncm.Matrix.new(lnw_v.len(), lnk_v.len())

    for i, (spec1, spec2, w) in enumerate(zip(specs1, specs2, w_a)):
        lnPk_a1 = np.array(spec1.peek_yv().dup_array())
        lnPk_a2 = np.array(spec2.peek_yv().dup_array())

        lnPk0 = np.log(np.exp(lnPk_a1[0]) + np.exp(lnPk_a2[0]))
        lnPk = np.log(np.exp(lnPk_a1) + np.exp(lnPk_a2)) - lnPk0
        lnPk_v = Ncm.Vector.new_array(lnPk)

        zm.set_row(i, lnPk_v)

    s = Ncm.Spline2dBicubic.notaknot_new()
    s.set(lnk_v, lnw_v, zm, True)

    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
    default_calib_file = Ncm.cfg_get_data_filename("hiprim_2f_spline.bin", True)
    s_calib = ser.from_binfile(default_calib_file)

    for w, lnk in product(w_a, lnk_v.dup_array()):
        Pk0 = s.eval(lnk, np.log(w))
        Pk1 = s_calib.eval(lnk, np.log(w))
        assert_allclose(Pk0, Pk1, rtol=1.0e-3)


def test_evolve_array(two_fluids, cosmo_qgrw):
    """Test NcHIPertTwoFluids evolve_array."""
    two_fluids.props.reltol = 1.0e-9
    m = two_fluids.evolve_array(cosmo=cosmo_qgrw, alphaf=-1.0e-1)

    assert all(np.isfinite(m.dup_array()))

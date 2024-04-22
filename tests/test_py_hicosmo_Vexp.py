#!/usr/bin/env python
#
# test_py_hicosmo_Vexp.py
#
# Sun Apr 14 17:31:52 2023
# Copyright  2024  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_py_hicosmo_Vexp.py
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

"""Tests on NcHICosmoVexp cosmological model."""

from itertools import product
import pytest
from numpy.testing import assert_allclose
import numpy as np

from numcosmo_py import Ncm, Nc

Ncm.cfg_init()


@pytest.fixture(name="vexp", params=[1.0e-1, -1.0e-1], ids=["d>0", "d<0"])
def fixture_vexp(request) -> Nc.HICosmoVexp:
    """Fixture for NcHICosmoVexp."""
    vexp = Nc.HICosmoVexp()
    vexp.props.dphi = request.param  # pylint: disable=no-member

    return vexp


def test_init(vexp):
    """Test NcHICosmoVexp initialization."""
    assert vexp is not None
    assert isinstance(vexp, Nc.HICosmoVexp)
    assert isinstance(vexp, Nc.HICosmo)
    assert isinstance(vexp, Ncm.Model)


def test_init_new(vexp):
    """Test NcHICosmoVexp initialization with new."""
    assert vexp is not None
    assert isinstance(vexp, Nc.HICosmoVexp)
    assert isinstance(vexp, Nc.HICosmo)
    assert isinstance(vexp, Ncm.Model)


def test_hubble_negative_dphi(vexp):
    """Test NcHICosmoVexp.hubble."""
    assert vexp is not None

    # Negative d_phi means matching Omega_L in the expansion
    # phase.
    if vexp.props.dphi < 0.0:  # pylint: disable=no-member
        assert_allclose(vexp.Omega_t0(), 1.0)
        assert_allclose(vexp.E2(0.0), vexp.props.OmegaL)  # pylint: disable=no-member


def test_hubble_positive_dphi(vexp):
    """Test NcHICosmoVexp.hubble."""
    assert vexp is not None

    # Positive d_phi means matching Omega_c in the expansion
    # phase.
    if vexp.props.dphi > 0.0:  # pylint: disable=no-member
        assert_allclose(vexp.Omega_t0(), 1.0)
        assert_allclose(vexp.E2(0.0), vexp.props.Omegac)  # pylint: disable=no-member


def test_tau_min_max(vexp):
    """Test NcHICosmoVexp.tau_min and NcHICosmoVexp.tau_max."""
    assert vexp is not None

    tau_min = vexp.tau_min()
    tau_max = vexp.tau_max()

    assert tau_min < 0.0
    assert tau_max > 0.0
    assert tau_min < tau_max


def test_tau_q(vexp):
    """Test NcHICosmoVexp.tau_q."""
    assert vexp is not None

    tau_min = vexp.tau_min()
    tau_max = vexp.tau_max()

    # Time of the transition classical to quantum during the contraction phase.
    tau_qt_c = vexp.tau_qt_c()
    # Time of the transition classical to quantum during the expansion phase.
    tau_qt_e = vexp.tau_qt_e()

    assert tau_qt_c > tau_min
    assert tau_qt_c < tau_max
    assert tau_qt_e > tau_min
    assert tau_qt_e < tau_max

    assert tau_qt_c < 0.0
    assert tau_qt_e > 0.0


def test_xb(vexp):
    """Test NcHICosmoVexp.xb."""
    assert vexp is not None

    xbc = vexp.xbc()
    xbe = vexp.xbe()

    assert xbc > 0.0
    assert xbe > 0.0

    assert_allclose(vexp.tau_xe(xbe), 0.0, atol=1.0e-6)
    assert_allclose(vexp.tau_xc(xbc), 0.0, atol=1.0e-6)


def _compute_tau_array(vexp, abs_tau_min=1.0e-10):
    """Compute an array of tau values for testing.

    This includes a range from tau_min to -abs_tau_min and from abs_tau_min to tau_max.
    This is useful since we want to test the behavior of the model around the bounce.
    """
    tau_min = vexp.tau_min()
    tau_max = vexp.tau_max()

    tau_a = np.concatenate(
        (
            np.geomspace(tau_min, -abs_tau_min, 5000),
            np.geomspace(abs_tau_min, tau_max, 5000),
        )
    )
    return tau_a


def test_E_tau_negative_dphi(vexp):
    """Test NcHICosmoVexp.E_tau."""
    assert vexp is not None

    tau_a = _compute_tau_array(vexp)

    for tau in tau_a:
        assert np.isfinite(vexp.E_tau(tau))

    tau0c = vexp.tau_xc(1.0)
    tau0e = vexp.tau_xe(1.0)

    # Negative d_phi means matching Omega_L in the expansion phase.
    if vexp.props.dphi < 0.0:  # pylint: disable=no-member
        assert_allclose(
            vexp.E_tau(tau0c) ** 2, vexp.props.Omegac  # pylint: disable=no-member
        )
        assert_allclose(
            vexp.E_tau(tau0e) ** 2, vexp.props.OmegaL  # pylint: disable=no-member
        )


def test_E_tau_positive_dphi(vexp):
    """Test NcHICosmoVexp.E_tau."""
    assert vexp is not None

    tau_a = _compute_tau_array(vexp)

    for tau in tau_a:
        assert np.isfinite(vexp.E_tau(tau))

    tau0c = vexp.tau_xc(1.0)
    tau0e = vexp.tau_xe(1.0)

    # Positive d_phi means matching Omega_c in the expansion phase.
    if vexp.props.dphi > 0.0:  # pylint: disable=no-member
        assert_allclose(
            vexp.E_tau(tau0c) ** 2, vexp.props.OmegaL  # pylint: disable=no-member
        )
        assert_allclose(
            vexp.E_tau(tau0e) ** 2, vexp.props.Omegac  # pylint: disable=no-member
        )


def test_Ricci_scale(vexp):
    """Test NcHICosmoVexp.Ricci_scale."""
    assert vexp is not None

    tau_a = _compute_tau_array(vexp)

    for tau in tau_a:
        assert np.isfinite(vexp.Ricci_scale(tau))
        Ricci_in_Hubble = np.abs(
            vexp.Ricci_scale(tau) * vexp.E_tau(tau) / vexp.RH_planck()
        )
        # Ricci scale is close to the Hubble scale. This works for our parameters and
        # the range of tau, since the the Hubble parameter goes through zero at the
        # bounce if the tau goes close to zero the first assertion below will fail.
        assert Ricci_in_Hubble > 1.0e-12
        assert Ricci_in_Hubble < 1.0e1


def test_tau_x(vexp):
    """Test NcHICosmoVexp.tau_x."""
    assert vexp is not None

    xc_a = np.geomspace(1.0e-3, 1.03, 100)
    xe_a = np.geomspace(1.0e-3, 1.03, 100)
    tau_xc_a = [vexp.tau_xc(xc) for xc in xc_a]
    tau_xe_a = [vexp.tau_xe(xe) for xe in xe_a]

    assert_allclose([vexp.xc_tau(tau) for tau in tau_xc_a], xc_a, atol=1.0e-7)
    assert_allclose([vexp.xe_tau(tau) for tau in tau_xe_a], xe_a, atol=1.0e-7)


def test_serialize(vexp):
    """Test NcHICosmoVexp serialization."""
    assert vexp is not None

    vexp.props.dphi = -1.234e-1  # pylint: disable=no-member
    vexp.props.Omegac = 0.321  # pylint: disable=no-member
    vexp.props.OmegaL = 0.679  # pylint: disable=no-member

    vexp.props.glue_de = True  # pylint: disable=no-member
    vexp.props.xb = 1.0e30  # pylint: disable=no-member

    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)

    vexp2 = ser.dup_obj(vexp)

    assert vexp2 is not None
    assert isinstance(vexp2, Nc.HICosmoVexp)
    assert isinstance(vexp2, Nc.HICosmo)
    assert isinstance(vexp2, Ncm.Model)

    assert vexp2.props.dphi == vexp.props.dphi  # pylint: disable=no-member
    assert vexp2.props.Omegac == vexp.props.Omegac  # pylint: disable=no-member
    assert vexp2.props.OmegaL == vexp.props.OmegaL  # pylint: disable=no-member
    assert vexp2.props.glue_de == vexp.props.glue_de  # pylint: disable=no-member
    assert vexp2.props.xb == vexp.props.xb  # pylint: disable=no-member


def test_eval_at(vexp):
    """Evaluate NcHICosmoVexp at a given time."""
    assert vexp is not None

    tau_a = _compute_tau_array(vexp)

    for tau in tau_a:
        assert np.isfinite(vexp.xe_tau(tau))
        assert np.isfinite(vexp.xc_tau(tau))
        assert np.isfinite(vexp.alpha(tau))
        assert np.isfinite(vexp.phi(tau))
        assert np.isfinite(vexp.Ricci_scale(tau))
        assert np.all(np.isfinite(vexp.x_y(tau)))


def test_iadiab_eval_at(vexp):
    """Evaluate NcHICosmoVexp implementation of NcHIPertAdiab at a given time."""
    assert vexp is not None

    tau_a = _compute_tau_array(vexp)
    k_a = np.geomspace(1.0e-3, 1.0e3, 5)

    for tau, k in product(tau_a, k_a):
        assert np.isfinite(Nc.HIPertIAdiab.eval_F1(vexp, tau, k))
        assert np.isfinite(Nc.HIPertIAdiab.eval_m(vexp, tau, k))
        assert np.isfinite(Nc.HIPertIAdiab.eval_nu(vexp, tau, k))
        assert np.isfinite(Nc.HIPertIAdiab.eval_xi(vexp, tau, k))


def test_igw_eval_at(vexp):
    """Evaluate NcHICosmoVexp implementation of NcHIPertIGW at a given time."""
    assert vexp is not None

    tau_a = _compute_tau_array(vexp)
    k_a = np.geomspace(1.0e-3, 1.0e3, 5)

    for tau, k in product(tau_a, k_a):
        assert np.isfinite(Nc.HIPertIGW.eval_F1(vexp, tau, k))
        assert np.isfinite(Nc.HIPertIGW.eval_m(vexp, tau, k))
        assert np.isfinite(Nc.HIPertIGW.eval_nu(vexp, tau, k))
        assert np.isfinite(Nc.HIPertIGW.eval_xi(vexp, tau, k))


def test_iem_eval_at(vexp):
    """Evaluate NcHICosmoVexp implementation of NcHIPertIEM at a given time."""
    assert vexp is not None

    tau_a = _compute_tau_array(vexp)
    k_a = np.geomspace(1.0e-3, 1.0e3, 5)

    for tau, k in product(tau_a, k_a):
        assert np.isfinite(Nc.HIPertIEM.eval_F1(vexp, tau, k))
        assert np.isfinite(Nc.HIPertIEM.eval_m(vexp, tau, k))
        assert np.isfinite(Nc.HIPertIEM.eval_nu(vexp, tau, k))
        assert np.isfinite(Nc.HIPertIEM.eval_xi(vexp, tau, k))

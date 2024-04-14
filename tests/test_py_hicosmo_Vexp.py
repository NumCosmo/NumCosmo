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

from numpy.testing import assert_allclose
import numpy as np

from numcosmo_py import Ncm, Nc

Ncm.cfg_init()


def test_init():
    """Test NcHICosmoVexp initialization."""
    vexp = Nc.HICosmoVexp()
    assert vexp is not None
    assert isinstance(vexp, Nc.HICosmoVexp)
    assert isinstance(vexp, Nc.HICosmo)
    assert isinstance(vexp, Ncm.Model)


def test_hubble_negative_dphi():
    """Test NcHICosmoVexp.hubble."""
    vexp = Nc.HICosmoVexp()
    assert vexp is not None

    # Negative d_phi means matching Omega_L in the expansion
    # phase.
    assert vexp.props.dphi < 0.0  # pylint: disable=no-member
    assert_allclose(vexp.Omega_t0(), 1.0)
    assert_allclose(vexp.E2(0.0), vexp.props.OmegaL)  # pylint: disable=no-member


def test_hubble_positive_dphi():
    """Test NcHICosmoVexp.hubble."""
    vexp = Nc.HICosmoVexp()
    assert vexp is not None

    # Positive d_phi means matching Omega_c in the expansion
    # phase.
    vexp.props.dphi = 1.0e-1  # pylint: disable=no-member
    assert_allclose(vexp.Omega_t0(), 1.0)
    assert_allclose(vexp.E2(0.0), vexp.props.Omegac)  # pylint: disable=no-member


def test_tau_min_max():
    """Test NcHICosmoVexp.tau_min and NcHICosmoVexp.tau_max."""
    vexp = Nc.HICosmoVexp()
    assert vexp is not None

    tau_min = vexp.tau_min()
    tau_max = vexp.tau_max()

    assert tau_min < 0.0
    assert tau_max > 0.0
    assert tau_min < tau_max


def test_tau_q():
    """Test NcHICosmoVexp.tau_q."""
    vexp = Nc.HICosmoVexp()
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


def test_xb():
    """Test NcHICosmoVexp.xb."""
    vexp = Nc.HICosmoVexp()
    assert vexp is not None

    xbc = vexp.xbc()
    xbe = vexp.xbe()

    assert xbc > 0.0
    assert xbe > 0.0

    assert_allclose(vexp.tau_xe(xbe), 0.0, atol=1.0e-7)


def test_tau_x():
    """Test NcHICosmoVexp.tau_x."""
    vexp = Nc.HICosmoVexp()
    assert vexp is not None

    xc_a = np.geomspace(1.0e-3, 1.03, 100)
    xe_a = np.geomspace(1.0e-3, 1.03, 100)
    tau_xc_a = [vexp.tau_xc(xc) for xc in xc_a]
    tau_xe_a = [vexp.tau_xe(xe) for xe in xe_a]

    assert_allclose([vexp.x_tau(tau) for tau in tau_xc_a], xc_a, atol=1.0e-7)
    assert_allclose([vexp.x_tau(tau) for tau in tau_xe_a], xe_a, atol=1.0e-7)


if __name__ == "__main__":
    test_init()
    test_hubble_negative_dphi()
    test_hubble_positive_dphi()
    test_tau_min_max()
    test_tau_q()
    test_xb()
    test_tau_x()

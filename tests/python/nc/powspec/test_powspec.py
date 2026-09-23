#!/usr/bin/env python
#
# test_powspec.py
#
# Wed Feb 14 13:44:00 2024
# Copyright  2023  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_powspec.py
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

"""Unit tests for NumCosmo power-spectra without CCL dependency."""

import pytest

import numpy as np
from numpy.testing import assert_allclose

from numcosmo_py import Ncm, Nc

pytestmark = pytest.mark.powspec

Ncm.cfg_init()


def test_powspec_class(cosmo: Nc.HICosmo) -> None:
    """Test the power spectrum class."""
    ps_ml = Nc.PowspecMLCBE.new()
    cbe = ps_ml.peek_cbe()
    cbe.use_ppf(True)

    ps_ml.prepare(cosmo)

    ik_min = ps_ml.get_intern_k_min()
    ik_max = ps_ml.get_intern_k_max()
    kmin = ps_ml.get_kmin()
    kmax = ps_ml.get_kmax()

    assert ik_min > 0
    assert ik_max > ik_min

    k_a = np.geomspace(kmin, kmax, 100)

    ps_a = np.array([ps_ml.eval(cosmo, 0.0, k) for k in k_a])

    assert np.all(np.isfinite(ps_a))


def test_powspec_class_deriv_z(cosmo: Nc.HICosmo) -> None:
    """Test the power spectrum class."""
    ps_ml = Nc.PowspecMLCBE.new()
    cbe = ps_ml.peek_cbe()
    cbe.use_ppf(True)

    ps_ml.prepare(cosmo)

    ik_min = ps_ml.get_intern_k_min()
    ik_max = ps_ml.get_intern_k_max()
    kmin = ps_ml.get_kmin()
    kmax = ps_ml.get_kmax()

    assert kmin < ik_min
    assert kmax > ik_max
    assert ik_min > 0
    assert ik_max > ik_min

    k_a = np.geomspace(kmin, kmax, 100)

    diff = Ncm.Diff.new()

    for k in k_a:
        dps, _ = diff.rc_d1_1_to_1(0.2, lambda z, k0: ps_ml.eval(cosmo, z, k0), k)
        assert_allclose(dps, ps_ml.deriv_z(cosmo, 0.2, k), atol=0.0, rtol=1.0e-10)


def test_powspec_default_deriv_matches_finite_difference(cosmo: Nc.HICosmo) -> None:
    """A spectrum without closed-form derivatives inherits a finite-difference one.

    NcPowspecMNLHaloFit does not implement deriv_z or deriv_k, so it exercises
    the base class default, checked here against a central difference on its
    own eval().
    """
    ps_ml = Nc.PowspecMLTransfer.new(Nc.TransferFuncEH.new())
    ps_nl = Nc.PowspecMNLHaloFit.new(ps_ml, 3.0, 1.0e-6)
    ps_nl.prepare(cosmo)

    for z, k in ((0.0, 0.05), (0.7, 0.2), (2.0, 1.0e-3), (1.0, 3.0)):
        h_z = 1.0e-4 * (1.0 + z)
        fd_z = (
            (ps_nl.eval(cosmo, z + h_z, k) - ps_nl.eval(cosmo, z - h_z, k))
            / (2.0 * h_z)
            if z > h_z
            else (
                -3.0 * ps_nl.eval(cosmo, z, k)
                + 4.0 * ps_nl.eval(cosmo, z + h_z, k)
                - ps_nl.eval(cosmo, z + 2.0 * h_z, k)
            )
            / (2.0 * h_z)
        )
        h_k = 1.0e-4 * k
        fd_k = (ps_nl.eval(cosmo, z, k + h_k) - ps_nl.eval(cosmo, z, k - h_k)) / (
            2.0 * h_k
        )

        assert_allclose(ps_nl.deriv_z(cosmo, z, k), fd_z, rtol=1.0e-6)
        assert_allclose(ps_nl.deriv_k(cosmo, z, k), fd_k, rtol=1.0e-6)

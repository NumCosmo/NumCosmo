#!/usr/bin/env python
#
# test_py_sbessel_integrator_fftl.py
#
# Thu Jan 10 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_py_sbessel_integrator_fftl.py
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
# with this program.  If not, see <http://www.gnu.org/licenses/>.

"""Unit tests for FFTL spherical Bessel integrator."""

from typing import Callable
import pytest
import numpy as np
from numpy.testing import assert_allclose

from numcosmo_py import Ncm


class TestSBesselIntegratorFFTL:
    """Tests for NcmSBesselIntegratorFFTL."""

    @pytest.fixture
    def integrator(self) -> Ncm.SBesselIntegratorFFTL:
        """Create an FFTL integrator."""
        return Ncm.SBesselIntegratorFFTL.new(0, 10)

    def test_create(self, integrator: Ncm.SBesselIntegratorFFTL) -> None:
        """Test integrator creation."""
        assert integrator is not None
        assert integrator.get_lmin() == 0
        assert integrator.get_lmax() == 10

    def test_integrate_constant_function(
        self, integrator: Ncm.SBesselIntegratorFFTL
    ) -> None:
        """Test integration of f(x) = 1."""

        def constant_one(_x) -> float:
            return 1.0

        integrator.prepare()

        # Test for different multipoles
        for ell in [0, 1, 2, 3, 4, 5, 100]:
            result = integrator.integrate_ell(constant_one, 0.0, 2000.0, ell)
            assert np.isfinite(result)
            # Result should be non-zero for integration of j_ell from 0 to 20
            assert abs(result) > 1e-10

    def test_integrate_gaussian(self, integrator: Ncm.SBesselIntegratorFFTL) -> None:
        """Test integration of f(x) = 1."""

        def constant_one(x) -> float:
            return np.exp(-0.5 * ((x - 1.0) / 0.1) ** 2)

        integrator.prepare()

        print("Testing Gaussian function integration:")
        # Test for different multipoles
        for ell in range(0, 11):
            result = integrator.integrate_ell(constant_one, 0.2, 2.0, ell)
            assert np.isfinite(result)
            # Result should be non-zero for integration of j_ell from 0 to 20
            print(f"ell={ell}, result={result}")
            assert np.isfinite(result)

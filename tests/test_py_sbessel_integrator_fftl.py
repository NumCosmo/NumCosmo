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

import time
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
        """Test integration of Gaussian function."""
        center = 1000.0
        std = 10.0

        def gaussian(x) -> float:
            return np.exp(-0.5 * ((x - center) / std) ** 2)

        integrator.prepare()

        # Test for different multipoles
        for ell in range(0, 1100):
            result = integrator.integrate_ell(
                gaussian, center - 20.0 * std, center + 20.0 * std, ell
            )
            assert np.isfinite(result)

    def test_integrate_rational(self, integrator: Ncm.SBesselIntegratorFFTL) -> None:
        """Test integration of rational function."""
        center = 1000.0
        std = 10.0

        def rational(x) -> float:
            return x**2 / (1 + ((x - center) / std) ** 2) ** 2

        integrator.prepare()

        # Test for different multipoles
        for ell in range(0, 1100):
            result = integrator.integrate_ell(
                rational, center - 20.0 * std, center + 20.0 * std, ell
            )
            assert np.isfinite(result)

    def test_integrate_vs_integrate_ell(
        self, integrator: Ncm.SBesselIntegratorFFTL
    ) -> None:
        """Test that integrate method agrees with integrate_ell for each ell."""
        center = 1000.0
        std = 10.0

        def test_func(x) -> float:
            return x**2 / (1 + ((x - center) / std) ** 2) ** 2

        a = center - 20.0 * std
        b = center + 20.0 * std

        integrator.set_lmax(1000)
        integrator.set_lmin(0)

        integrator.prepare()

        # Get results from integrate_ell for each ell
        results_ell = []
        for ell in range(integrator.get_lmin(), integrator.get_lmax() + 1):
            result = integrator.integrate_ell(test_func, a, b, ell)
            results_ell.append(result)

        # Get results from integrate (all ells at once)
        n_ell = integrator.get_lmax() - integrator.get_lmin() + 1
        results_vec = Ncm.Vector.new(n_ell)
        integrator.integrate(test_func, a, b, results_vec)
        results_all = np.array([results_vec.get(i) for i in range(n_ell)])

        # Compare results
        assert_allclose(results_all, results_ell, rtol=1e-10, atol=0.0)

    def test_integrate_performance(self, integrator: Ncm.SBesselIntegratorFFTL) -> None:
        """Time comparison between integrate_ell loop and integrate method."""
        center = 1000.0
        std = 10.0

        def test_func(x) -> float:
            return x**2 / (1 + ((x - center) / std) ** 2) ** 2

        a = center - 20.0 * std
        b = center + 20.0 * std

        integrator.set_lmax(1000)
        integrator.set_lmin(0)

        integrator.prepare()

        n_ell = integrator.get_lmax() - integrator.get_lmin() + 1

        # Time integrate_ell approach
        start_time = time.time()
        results_ell = []
        for ell in range(integrator.get_lmin(), integrator.get_lmax() + 1):
            result = integrator.integrate_ell(test_func, a, b, ell)
            results_ell.append(result)
        time_integrate_ell = time.time() - start_time

        # Time integrate approach
        start_time = time.time()
        results_vec = Ncm.Vector.new(n_ell)
        integrator.integrate(test_func, a, b, results_vec)
        results_all = np.array([results_vec.get(i) for i in range(n_ell)])
        time_integrate = time.time() - start_time

        # Print timing results
        print(f"\nPerformance comparison for {n_ell} multipoles:")
        print(f"  integrate_ell (sequential): {time_integrate_ell:.4f} seconds")
        print(f"  integrate (all at once):    {time_integrate:.4f} seconds")
        print(f"  Speedup: {time_integrate_ell / time_integrate:.2f}x")

        # Verify results match
        assert_allclose(results_all, results_ell, rtol=1e-10, atol=0.0)

        # Assert that integrate is faster (should be significant)
        assert time_integrate < time_integrate_ell

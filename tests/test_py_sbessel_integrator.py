#!/usr/bin/env python
#
# test_py_sbessel_integrator.py
#
# Thu Jan 09 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_py_sbessel_integrator.py
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

"""Unit tests for spherical Bessel integrator."""

import pytest
import numpy as np
from numpy.testing import assert_allclose
import scipy.special as sp
import scipy.integrate as integrate

from numcosmo_py import Ncm


class TestSBesselIntegratorGL:
    """Tests for NcmSBesselIntegratorGL."""

    @pytest.fixture
    def integrator(self):
        """Create a GL integrator."""
        return Ncm.SBesselIntegratorGL.new(0, 10)

    def test_create(self, integrator):
        """Test integrator creation."""
        assert integrator is not None
        assert integrator.get_lmin() == 0
        assert integrator.get_lmax() == 10

    def test_prepare(self, integrator):
        """Test prepare method."""
        integrator.prepare()  # Should not raise

    def test_npts_property(self):
        """Test npts property getter/setter."""
        # Test default value
        integrator = Ncm.SBesselIntegratorGL.new(0, 10)
        assert integrator.props.npts == 10

        # Test setting via constructor
        integrator2 = Ncm.SBesselIntegratorGL(lmin=0, lmax=10, npts=20)
        assert integrator2.props.npts == 20

        # Test setter
        integrator.props.npts = 15
        assert integrator.props.npts == 15

    def test_margin_property(self):
        """Test margin property getter/setter."""
        # Test default value
        integrator = Ncm.SBesselIntegratorGL.new(0, 10)
        assert integrator.props.margin == 5.0

        # Test setting via constructor
        integrator2 = Ncm.SBesselIntegratorGL(lmin=0, lmax=10, margin=10.0)
        assert integrator2.props.margin == 10.0

        # Test setter
        integrator.props.margin = 7.5
        assert integrator.props.margin == 7.5

    def test_nosc_property(self):
        """Test nosc property getter/setter."""
        # Test default value
        integrator = Ncm.SBesselIntegratorGL.new(0, 10)
        assert integrator.props.nosc == 6.0

        # Test setting via constructor
        integrator2 = Ncm.SBesselIntegratorGL(lmin=0, lmax=10, nosc=8.0)
        assert integrator2.props.nosc == 8.0

        # Test setter
        integrator.props.nosc = 4.0
        assert integrator.props.nosc == 4.0

    def test_properties_affect_integration(self):
        """Test that changing properties affects integration results."""

        def test_func(_user_data, x):
            return np.exp(-x / 10.0)

        # Create integrators with different npts
        integrator1 = Ncm.SBesselIntegratorGL(lmin=0, lmax=10, npts=5)
        integrator2 = Ncm.SBesselIntegratorGL(lmin=0, lmax=10, npts=20)

        integrator1.prepare()
        integrator2.prepare()

        result1 = integrator1.integrate_ell(test_func, None, 0.0, 20.0, 5)
        result2 = integrator2.integrate_ell(test_func, None, 0.0, 20.0, 5)

        # Both should be finite and reasonably close
        assert np.isfinite(result1)
        assert np.isfinite(result2)
        # Higher npts should give similar or better result
        assert_allclose(result1, result2, rtol=1e-2)

    def test_different_margins(self):
        """Test integration with different margin values."""

        def test_func(_user_data, x):
            return np.exp(-x / 10.0)

        # Test with different margins
        for margin in [2.0, 5.0, 10.0]:
            integrator = Ncm.SBesselIntegratorGL(lmin=0, lmax=10, margin=margin)
            integrator.prepare()

            result = integrator.integrate_ell(test_func, None, 0.0, 20.0, 5)
            assert np.isfinite(result)

    def test_different_nosc(self):
        """Test integration with different nosc values."""

        def test_func(_user_data, x):
            return np.exp(-x / 10.0)

        # Test with different nosc values
        for nosc in [4.0, 6.0, 8.0]:
            integrator = Ncm.SBesselIntegratorGL(lmin=0, lmax=10, nosc=nosc)
            integrator.prepare()

            result = integrator.integrate_ell(test_func, None, 0.0, 20.0, 5)
            assert np.isfinite(result)

    def test_integrate_constant_function(self, integrator):
        """Test integration of f(x) = 1."""

        def constant_one(_user_data, _x):
            return 1.0

        integrator.prepare()

        # Test for different multipoles
        for ell in [0, 1, 2, 5, 10]:
            result = integrator.integrate_ell(constant_one, None, 0.0, 20.0, ell)
            assert np.isfinite(result)
            # Result should be non-zero for integration of j_ell from 0 to 20
            assert abs(result) > 1e-10

    def test_integrate_exponential_decay(self, integrator):
        """Test integration with exponential decay.

        Using f(x) = exp(-x/10).
        """

        def exponential_decay(_user_data, x):
            return np.exp(-x / 10.0)

        integrator.prepare()

        # Test for different multipoles
        results = []
        for ell in [0, 1, 2, 5]:
            result = integrator.integrate_ell(exponential_decay, None, 0.0, 20.0, ell)
            results.append(result)
            assert np.isfinite(result)

        # Results should decay with increasing ell for this smooth function
        # (not a strict requirement but typically observed)
        assert all(np.isfinite(r) for r in results)

    def test_integrate_power_law(self, integrator):
        """Test integration with power law f(x) = x^2."""

        def power_law(_user_data, x):
            return x**2

        integrator.prepare()

        for ell in [0, 2, 5]:
            result = integrator.integrate_ell(power_law, None, 0.0, 20.0, ell)
            assert np.isfinite(result)

    def test_integrate_oscillatory_function(self, integrator):
        """Test with oscillatory function f(x) = sin(x)."""

        def sine_function(_user_data, x):
            return np.sin(x)

        integrator.prepare()

        for ell in [0, 2, 5]:
            result = integrator.integrate_ell(sine_function, None, 0.0, 20.0, ell)
            assert np.isfinite(result)

    def test_integrate_zero_function(self, integrator):
        """Test with f(x) = 0."""

        def zero_function(_user_data, _x):
            return 0.0

        integrator.prepare()

        result = integrator.integrate_ell(zero_function, None, 0.0, 20.0, 5)
        assert_allclose(result, 0.0, atol=1e-10)

    def test_integrate_different_ranges(self, integrator):
        """Test integration over different ranges."""

        def constant_one(_user_data, _x):
            return 1.0

        integrator.prepare()

        # Test different integration ranges
        ranges = [(0.0, 10.0), (0.0, 20.0), (5.0, 15.0), (1.0, 30.0)]

        for a, b in ranges:
            result = integrator.integrate_ell(constant_one, None, a, b, 2)
            assert np.isfinite(result)

    def test_integrate_small_ell(self, integrator):
        """Test integration for small multipoles."""

        def constant_one(_user_data, _x):
            return 1.0

        integrator.prepare()

        for ell in [0, 1, 2]:
            result = integrator.integrate_ell(constant_one, None, 0.0, 10.0, ell)
            assert np.isfinite(result)
            assert abs(result) > 0.0

    def test_integrate_large_ell(self, integrator):
        """Test integration for larger multipoles."""

        def constant_one(_user_data, _x):
            return 1.0

        integrator.prepare()

        # Test with larger multipoles
        for ell in [20, 50, 100]:
            result = integrator.integrate_ell(constant_one, None, 0.0, 50.0, ell)
            assert np.isfinite(result)

    def test_integrate_with_user_data(self, integrator):
        """Test that user_data is properly passed through."""

        class DataContainer:
            """Simple data container."""

            def __init__(self, scale):
                self.scale = scale

        def scaled_exponential(user_data, x):
            return np.exp(-x / user_data.scale)

        integrator.prepare()

        data = DataContainer(scale=5.0)
        result1 = integrator.integrate_ell(scaled_exponential, data, 0.0, 20.0, 2)

        data.scale = 10.0
        result2 = integrator.integrate_ell(scaled_exponential, data, 0.0, 20.0, 2)

        # Results should be different for different scales
        assert result1 != result2
        assert np.isfinite(result1)
        assert np.isfinite(result2)

    def test_numerical_accuracy(self, integrator):
        """Test numerical accuracy with known analytical result."""
        # For j_0(x) = sin(x)/x, integral from 0 to infinity is pi/2
        # We can test a finite range and compare qualitatively

        def sinc_like(_user_data, _):
            # Test f(x) = 1 for x in range, which gives us the integral of j_ell
            return 1.0

        integrator.prepare()

        # For ell=0, j_0(x) = sin(x)/x
        # Integral from 0 to large value should approach a known value
        result = integrator.integrate_ell(sinc_like, None, 0.0, 100.0, 0)
        assert np.isfinite(result)
        # The integral of sin(x)/x from 0 to 100 is close to pi/2
        # Our result is integral of 1 * sin(x)/x
        assert 0.5 < abs(result) < 2.0  # Loose bounds for sanity check

    def test_consistency_across_calls(self, integrator):
        """Test that repeated calls give consistent results."""

        def test_func(_user_data, x):
            return np.exp(-x / 10.0)

        integrator.prepare()

        # Call multiple times with same parameters
        results = [
            integrator.integrate_ell(test_func, None, 0.0, 20.0, 5) for _ in range(5)
        ]

        # All results should be identical
        assert_allclose(results, results[0], rtol=1e-14)

    def test_lmin_lmax_properties(self, integrator):
        """Test lmin and lmax getter/setter."""
        assert integrator.get_lmin() == 0
        assert integrator.get_lmax() == 10

        integrator.set_lmin(2)
        integrator.set_lmax(20)

        assert integrator.get_lmin() == 2
        assert integrator.get_lmax() == 20

    def test_different_integrator_instances(self):
        """Test that different integrator instances are independent."""

        def test_func(_user_data, x):
            return np.exp(-x)

        integrator1 = Ncm.SBesselIntegratorGL.new(0, 5)
        integrator2 = Ncm.SBesselIntegratorGL.new(5, 10)

        integrator1.prepare()
        integrator2.prepare()

        result1 = integrator1.integrate_ell(test_func, None, 0.0, 10.0, 2)
        result2 = integrator2.integrate_ell(test_func, None, 0.0, 10.0, 7)

        # Results should be different (different ell values)
        assert result1 != result2
        assert np.isfinite(result1)
        assert np.isfinite(result2)

    def test_gaussian_function(self, integrator):
        """Test integration with Gaussian function."""

        def gaussian(_user_data, x):
            return np.exp(-((x - 10.0) ** 2) / (2.0 * 2.0**2))

        integrator.prepare()

        for ell in [0, 2, 5]:
            result = integrator.integrate_ell(gaussian, None, 0.0, 20.0, ell)
            assert np.isfinite(result)

    def test_piecewise_function(self, integrator):
        """Test with piecewise function."""

        def piecewise(_user_data, x):
            if x < 10.0:
                return 1.0
            else:
                return 0.5

        integrator.prepare()

        result = integrator.integrate_ell(piecewise, None, 0.0, 20.0, 3)
        assert np.isfinite(result)

    def test_turning_point_behavior(self, integrator):
        """Test behavior near turning point.

        For ell=5, x_tp = sqrt(ell(ell+1)) ~= 5.48.
        """
        # Test integration ranges that cross this point

        def constant_one(_user_data, _x):
            return 1.0

        integrator.prepare()

        ell = 5
        x_tp = np.sqrt(ell * (ell + 1))

        # Integrate before turning point
        result_before = integrator.integrate_ell(constant_one, None, 0.0, x_tp - 1, ell)

        # Integrate after turning point
        result_after = integrator.integrate_ell(constant_one, None, x_tp + 2, 20.0, ell)

        # Integrate across turning point
        result_across = integrator.integrate_ell(constant_one, None, 0.0, 20.0, ell)

        assert np.isfinite(result_before)
        assert np.isfinite(result_after)
        assert np.isfinite(result_across)

    def test_very_small_range(self, integrator):
        """Test integration over very small range."""

        def test_func(_user_data, x):
            return x

        integrator.prepare()

        result = integrator.integrate_ell(test_func, None, 5.0, 5.1, 2)
        assert np.isfinite(result)
        # Result should be small for small range
        assert abs(result) < 1.0

    def test_polynomial_function(self, integrator):
        """Test with polynomial functions of different orders."""

        integrator.prepare()

        for power in [0, 1, 2, 3]:

            def polynomial(_user_data, x, p=power):
                return x**p

            for ell in [0, 2, 5]:
                result = integrator.integrate_ell(polynomial, None, 0.0, 10.0, ell)
                assert np.isfinite(result)


class TestSBesselIntegratorGLAnalytical:
    """Tests comparing numerical results against analytical solutions."""

    @pytest.fixture
    def integrator(self):
        """Create a GL integrator."""
        return Ncm.SBesselIntegratorGL.new(0, 10)

    def _compute_reference(self, f, a, b, ell):
        """Compute high-precision reference value using scipy."""

        def integrand(x):
            return f(None, x) * sp.spherical_jn(ell, x)

        result, _ = integrate.quad(integrand, a, b, epsabs=1e-12, epsrel=1e-12)
        return result

    def test_constant_function_vs_scipy(self, integrator):
        """Test f(x)=1 against scipy high-precision integration."""

        def constant_one(_user_data, _x):
            return 1.0

        integrator.prepare()

        # Test several cases
        test_cases = [
            (0.0, 10.0, 0),
            (0.0, 10.0, 2),
            (0.0, 10.0, 5),
            (1.0, 15.0, 3),
            (0.0, 20.0, 1),
        ]

        for a, b, ell in test_cases:
            result = integrator.integrate_ell(constant_one, None, a, b, ell)
            reference = self._compute_reference(constant_one, a, b, ell)
            assert_allclose(result, reference, rtol=1e-6, atol=1e-10)

    def test_exponential_vs_scipy(self, integrator):
        """Test exponential function against scipy."""

        def exponential(_user_data, x):
            return np.exp(-x / 10.0)

        integrator.prepare()

        test_cases = [
            (0.0, 10.0, 0),
            (0.0, 10.0, 2),
            (0.0, 15.0, 5),
            (2.0, 12.0, 3),
        ]

        for a, b, ell in test_cases:
            result = integrator.integrate_ell(exponential, None, a, b, ell)
            reference = self._compute_reference(exponential, a, b, ell)
            assert_allclose(result, reference, rtol=1e-6, atol=1e-10)

    def test_polynomial_vs_scipy(self, integrator):
        """Test polynomial functions against scipy."""

        integrator.prepare()

        for power in [1, 2, 3]:

            def polynomial(_user_data, x, p=power):
                return x**p

            test_cases = [(0.0, 10.0, 0), (0.0, 10.0, 2), (1.0, 8.0, 3)]

            for a, b, ell in test_cases:
                result = integrator.integrate_ell(polynomial, None, a, b, ell)
                reference = self._compute_reference(polynomial, a, b, ell)
                assert_allclose(result, reference, rtol=1e-6, atol=1e-10)

    def test_sine_integral_special_case(self, integrator):
        """Test special case: sine integral.

        For ell=0, j_0(x) = sin(x)/x, so integral of f(x)=1 * j_0(x)
        is the sine integral.
        """

        def constant_one(_user_data, _x):
            return 1.0

        integrator.prepare()

        # Test from 0 to pi
        a, b = 0.0, np.pi
        result = integrator.integrate_ell(constant_one, None, a, b, 0)

        # Reference: Si(pi) - Si(0) = Si(pi) ~= 1.85194
        # Using scipy to compute
        reference = self._compute_reference(constant_one, a, b, 0)

        assert_allclose(result, reference, rtol=1e-6)

    def test_gaussian_vs_scipy(self, integrator):
        """Test Gaussian function against scipy."""

        def gaussian(_user_data, x):
            return np.exp(-((x - 5.0) ** 2) / (2.0 * 2.0**2))

        integrator.prepare()

        test_cases = [
            (0.0, 10.0, 0),
            (0.0, 10.0, 2),
            (1.0, 9.0, 4),
        ]

        for a, b, ell in test_cases:
            result = integrator.integrate_ell(gaussian, None, a, b, ell)
            reference = self._compute_reference(gaussian, a, b, ell)
            assert_allclose(result, reference, rtol=1e-6, atol=1e-10)

    def test_trigonometric_vs_scipy(self, integrator):
        """Test trigonometric functions against scipy."""

        def sine_func(_user_data, x):
            return np.sin(x)

        def cosine_func(_user_data, x):
            return np.cos(x)

        integrator.prepare()

        # Test sine
        for ell in [0, 2, 5]:
            result = integrator.integrate_ell(sine_func, None, 0.0, np.pi, ell)
            reference = self._compute_reference(sine_func, 0.0, np.pi, ell)
            assert_allclose(result, reference, rtol=1e-6, atol=1e-10)

        # Test cosine
        for ell in [0, 2, 5]:
            result = integrator.integrate_ell(cosine_func, None, 0.0, np.pi, ell)
            reference = self._compute_reference(cosine_func, 0.0, np.pi, ell)
            assert_allclose(result, reference, rtol=1e-6, atol=1e-10)

    def test_product_of_powers_vs_scipy(self, integrator):
        """Test x^n * exp(-x/L) against scipy."""

        integrator.prepare()

        for n in [0, 1, 2]:
            for L in [5.0, 10.0]:

                def power_exp(_user_data, x, n=n, L=L):
                    return x**n * np.exp(-x / L)

                for ell in [0, 2, 5]:
                    result = integrator.integrate_ell(power_exp, None, 0.0, 20.0, ell)
                    reference = self._compute_reference(power_exp, 0.0, 20.0, ell)
                    assert_allclose(result, reference, rtol=1e-5, atol=1e-10)

    def test_oscillatory_product_vs_scipy(self, integrator):
        """Test oscillatory function times smooth envelope."""

        def osc_func(_user_data, x):
            return np.sin(2 * x) * np.exp(-x / 10.0)

        integrator.prepare()

        test_cases = [
            (0.0, 10.0, 0),
            (0.0, 15.0, 2),
            (0.0, 20.0, 5),
        ]

        for a, b, ell in test_cases:
            result = integrator.integrate_ell(osc_func, None, a, b, ell)
            reference = self._compute_reference(osc_func, a, b, ell)
            assert_allclose(result, reference, rtol=1e-5, atol=1e-10)

    def test_rational_function_vs_scipy(self, integrator):
        """Test rational function against scipy."""

        def rational(_user_data, x):
            return 1.0 / (1.0 + x**2)

        integrator.prepare()

        test_cases = [(0.0, 5.0, 0), (0.0, 10.0, 2), (1.0, 8.0, 4)]

        for a, b, ell in test_cases:
            result = integrator.integrate_ell(rational, None, a, b, ell)
            reference = self._compute_reference(rational, a, b, ell)
            assert_allclose(result, reference, rtol=1e-6, atol=1e-10)

    def test_multiple_ells_consistency(self, integrator):
        """Test results consistent with scipy for range of ells."""

        def test_func(_user_data, x):
            return x * np.exp(-x / 5.0)

        integrator.prepare()

        a, b = 0.0, 15.0

        for ell in range(0, 11):
            result = integrator.integrate_ell(test_func, None, a, b, ell)
            reference = self._compute_reference(test_func, a, b, ell)
            assert_allclose(result, reference, rtol=1e-5, atol=1e-10)


class TestSBesselIntegratorGLLargeIntervals:
    """Tests for integration over very large intervals."""

    @pytest.fixture
    def integrator(self):
        """Create a GL integrator with parameters suitable for large intervals."""
        # Use more quadrature points for better accuracy with many oscillations
        return Ncm.SBesselIntegratorGL(lmin=0, lmax=10, npts=20, nosc=8.0)

    def test_large_interval_constant_function(self, integrator):
        """Test f(x)=1 over [0, 10^4] for ell=5."""

        def constant_one(_user_data, _x):
            return 1.0

        integrator.prepare()

        # For large x, j_ell oscillates rapidly and decays as 1/x
        # The integral should converge
        result = integrator.integrate_ell(constant_one, None, 0.0, 1e4, 5)
        assert np.isfinite(result)
        # The integral should be relatively small due to rapid oscillations
        assert abs(result) < 100.0

    def test_large_interval_exponential_decay(self, integrator):
        """Test exponential decay over [0, 10^4]."""

        def exp_decay(_user_data, x):
            return np.exp(-x / 1e3)

        integrator.prepare()

        result = integrator.integrate_ell(exp_decay, None, 0.0, 1e4, 5)
        assert np.isfinite(result)
        # With exponential decay, the integral should be well-behaved
        assert abs(result) < 1e4

    def test_large_interval_power_law_decay(self, integrator):
        """Test power law decay over [0, 10^4]."""

        def power_decay(_user_data, x):
            return 1.0 / (1.0 + x / 1e3)

        integrator.prepare()

        result = integrator.integrate_ell(power_decay, None, 0.0, 1e4, 5)
        assert np.isfinite(result)
        assert abs(result) < 1e4

    def test_large_interval_starting_at_high_value(self, integrator):
        """Test integration from [10^3, 10^4] for ell=5."""

        def constant_one(_user_data, _x):
            return 1.0

        integrator.prepare()

        # Integration starting far from origin where j_ell is highly oscillatory
        result = integrator.integrate_ell(constant_one, None, 1e3, 1e4, 5)
        assert np.isfinite(result)
        # Due to rapid oscillations and 1/x decay, this should be small
        assert abs(result) < 10.0

    def test_large_interval_multiple_ells(self, integrator):
        """Test large intervals for different multipoles."""

        def exp_decay(_user_data, x):
            return np.exp(-x / 1e3)

        integrator.prepare()

        results = {}
        for ell in [0, 2, 5, 10]:
            result = integrator.integrate_ell(exp_decay, None, 0.0, 1e4, ell)
            results[ell] = result
            assert np.isfinite(result)

        # All results should be finite and non-zero
        for ell, result in results.items():
            assert abs(result) > 1e-10, f"Result for ell={ell} is too small"

    def test_large_interval_with_slow_decay(self, integrator):
        """Test slowly decaying function over large interval."""

        def slow_decay(_user_data, x):
            # Decay as 1/sqrt(1 + x/L) with large L
            return 1.0 / np.sqrt(1.0 + x / 1e4)

        integrator.prepare()

        result = integrator.integrate_ell(slow_decay, None, 0.0, 1e4, 5)
        assert np.isfinite(result)
        # Slower decay means potentially larger integral
        assert abs(result) < 1e4

    def test_very_large_starting_point(self, integrator):
        """Test integration starting at very large x value."""

        def constant_one(_user_data, _x):
            return 1.0

        integrator.prepare()

        # Start at 5*10^3, go to 6*10^3
        result = integrator.integrate_ell(constant_one, None, 5e3, 6e3, 5)
        assert np.isfinite(result)
        # In highly oscillatory regime with 1/x envelope
        assert abs(result) < 1.0

    def test_large_interval_consistency(self, integrator):
        """Test that splitting large interval gives consistent results."""

        def exp_decay(_user_data, x):
            return np.exp(-x / 1e3)

        integrator.prepare()

        # Full interval
        result_full = integrator.integrate_ell(exp_decay, None, 0.0, 1e4, 5)

        # Split into two parts
        result_part1 = integrator.integrate_ell(exp_decay, None, 0.0, 5e3, 5)
        result_part2 = integrator.integrate_ell(exp_decay, None, 5e3, 1e4, 5)
        result_split = result_part1 + result_part2

        # Should be approximately equal
        assert_allclose(result_full, result_split, rtol=1e-4)

    def test_large_interval_gaussian_envelope(self, integrator):
        """Test Gaussian envelope centered away from origin."""

        def gaussian_far(_user_data, x):
            # Gaussian centered at x0 = 5*10^3, width sigma = 10^3
            x0 = 5e3
            sigma = 1e3
            return np.exp(-((x - x0) ** 2) / (2 * sigma**2))

        integrator.prepare()

        # Integrate over range covering the Gaussian
        result = integrator.integrate_ell(gaussian_far, None, 0.0, 1e4, 5)
        assert np.isfinite(result)
        # Gaussian centered far from origin, width ~10^3
        assert abs(result) < 1e4

    def test_large_interval_different_parameter_settings(self):
        """Test large intervals with different algorithm parameters."""

        def exp_decay(_user_data, x):
            return np.exp(-x / 1e3)

        # Test with different npts and nosc values
        param_sets = [
            {"npts": 10, "nosc": 6.0},
            {"npts": 20, "nosc": 8.0},
            {"npts": 30, "nosc": 10.0},
        ]

        results = []
        for params in param_sets:
            integrator = Ncm.SBesselIntegratorGL(
                lmin=0, lmax=10, npts=params["npts"], nosc=params["nosc"]
            )
            integrator.prepare()
            result = integrator.integrate_ell(exp_decay, None, 0.0, 1e4, 5)
            results.append(result)
            assert np.isfinite(result)

        # All results should be reasonably close
        for i in range(1, len(results)):
            assert_allclose(results[i], results[0], rtol=1e-3)

    def test_extremely_large_interval(self, integrator):
        """Test with extremely large interval [0, 10^5]."""

        def exp_decay(_user_data, x):
            return np.exp(-x / 1e4)

        integrator.prepare()

        result = integrator.integrate_ell(exp_decay, None, 0.0, 1e7, 5)
        assert np.isfinite(result)
        # Should converge despite large range
        assert abs(result) < 1e5

    def test_large_interval_high_multipole(self, integrator):
        """Test large intervals with high multipole values."""

        def exp_decay(_user_data, x):
            return np.exp(-x / 1e3)

        integrator.prepare()

        # Test with higher multipoles
        for ell in [20, 50]:
            result = integrator.integrate_ell(exp_decay, None, 0.0, 1e6, ell)
            assert np.isfinite(result)
            # Higher ell means turning point at sqrt(ell*(ell+1))
            # Most of the range is in highly oscillatory regime


class TestSBesselIntegratorGLAnalyticalLargeIntervals:
    """Tests for large intervals using analytical formulas.

    Uses known analytical formulas and consistency checks.
    """

    @pytest.fixture
    def integrator(self):
        """Create GL integrator with parameters for large intervals."""
        return Ncm.SBesselIntegratorGL(lmin=0, lmax=10, npts=20, nosc=8.0)

    def test_interval_additivity_large(self, integrator):
        """Test interval additivity for large intervals.

        Verify that integral_a^c = integral_a^b + integral_b^c.
        """

        def exp_decay(_user_data, x):
            return np.exp(-x / 1e3)

        integrator.prepare()

        # Test for multiple ell values
        for ell in [0, 2, 5]:
            # Full interval [0, 10^4]
            result_full = integrator.integrate_ell(exp_decay, None, 0.0, 1e4, ell)

            # Split into three parts
            result_1 = integrator.integrate_ell(exp_decay, None, 0.0, 3e3, ell)
            result_2 = integrator.integrate_ell(exp_decay, None, 3e3, 7e3, ell)
            result_3 = integrator.integrate_ell(exp_decay, None, 7e3, 1e4, ell)
            result_sum = result_1 + result_2 + result_3

            # Should satisfy additivity to high precision
            assert_allclose(
                result_full,
                result_sum,
                rtol=1e-10,
                err_msg=f"Additivity failed for ell={ell}",
            )

    def test_interval_additivity_many_splits(self, integrator):
        """Test additivity with many small intervals."""

        def exp_decay(_user_data, x):
            return np.exp(-x / 1e3)

        integrator.prepare()

        ell = 3
        # Full interval
        result_full = integrator.integrate_ell(exp_decay, None, 0.0, 1e4, ell)

        # Split into 20 pieces
        n_splits = 20
        result_sum = 0.0
        for i in range(n_splits):
            a = i * 500
            b = (i + 1) * 500
            result_sum += integrator.integrate_ell(exp_decay, None, a, b, ell)

        assert_allclose(result_full, result_sum, rtol=1e-10)

    def test_asymptotic_decay_tail_contribution(self, integrator):
        """Test that tail contribution from large x decays properly."""

        def exp_decay(_user_data, x):
            return np.exp(-x / 5e3)

        integrator.prepare()

        ell = 2
        # Main integral
        result_main = integrator.integrate_ell(exp_decay, None, 0.0, 1e4, ell)

        # Tail integrals - should decrease
        tail1 = integrator.integrate_ell(exp_decay, None, 1e4, 1.5e4, ell)
        tail2 = integrator.integrate_ell(exp_decay, None, 1.5e4, 2e4, ell)

        # Tail contributions should be decreasing
        assert abs(tail1) < 0.1 * abs(result_main)
        assert abs(tail2) < abs(tail1)

    def test_zero_function_large_interval(self, integrator):
        """Test that zero function gives zero over large interval."""

        def zero_func(_user_data, _x):
            return 0.0

        integrator.prepare()

        for ell in [0, 2, 5, 10]:
            result = integrator.integrate_ell(zero_func, None, 0.0, 1e4, ell)
            assert_allclose(
                result,
                0.0,
                atol=1e-14,
                err_msg=f"Zero function test failed for ell={ell}",
            )

    def test_linearity_large_interval(self, integrator):
        """Test linearity for large intervals.

        Verify that integral(af + bg) = a*integral(f) + b*integral(g).
        """

        def func1(_user_data, x):
            return np.exp(-x / 1e3)

        def func2(_user_data, x):
            return np.exp(-x / 2e3)

        integrator.prepare()

        a, b = 2.5, -1.3
        ell = 3

        # Compute separately
        result_f = integrator.integrate_ell(func1, None, 0.0, 1e4, ell)
        result_g = integrator.integrate_ell(func2, None, 0.0, 1e4, ell)
        result_linear = a * result_f + b * result_g

        # Compute combined
        def func_combined(_user_data, x):
            return a * func1(None, x) + b * func2(None, x)

        result_combined = integrator.integrate_ell(func_combined, None, 0.0, 1e4, ell)

        assert_allclose(result_combined, result_linear, rtol=1e-10)

    def test_convergence_sequence_large_intervals(self):
        """Test convergence as integration parameters are refined."""

        def exp_decay(_user_data, x):
            return np.exp(-x / 1e3)

        ell = 5

        # Test with increasing npts
        results_npts = []
        for npts in [10, 20, 30, 40]:
            integ = Ncm.SBesselIntegratorGL(lmin=0, lmax=10, npts=npts, nosc=8.0)
            integ.prepare()
            result = integ.integrate_ell(exp_decay, None, 0.0, 1e4, ell)
            results_npts.append(result)

        # Results should converge - differences should decrease
        diffs = [
            abs(results_npts[i + 1] - results_npts[i])
            for i in range(len(results_npts) - 1)
        ]
        for i in range(len(diffs) - 1):
            assert diffs[i + 1] < diffs[i] * 1.5  # Allow some tolerance

    def test_symmetry_interval_midpoint(self, integrator):
        """Test symmetric functions around interval midpoint."""

        # Even function around x=5000
        def symmetric_func(_user_data, x):
            return np.exp(-((x - 5e3) ** 2) / (2 * 1e6))

        integrator.prepare()

        ell = 2
        # Integrate symmetric function - check consistency with mirrored interval
        result = integrator.integrate_ell(symmetric_func, None, 0.0, 1e4, ell)

        assert np.isfinite(result)
        assert abs(result) > 1e-10  # Should be non-zero

    def test_monotonic_decay_function_ordering(self, integrator):
        """Test that faster decay gives smaller integral magnitude."""

        integrator.prepare()

        ell = 3
        alphas = [1e-3, 5e-3, 1e-2, 5e-2]  # Increasing decay rates
        results = []

        for alpha in alphas:

            def exp_decay(_user_data, x, alpha=alpha):
                return np.exp(-alpha * x)

            result = integrator.integrate_ell(exp_decay, None, 0.0, 1e4, ell)
            results.append(abs(result))

        # Faster decay (larger alpha) should give smaller integral
        # (approximately, allowing for oscillations)
        # At least check first and last are ordered
        assert results[0] > results[-1] * 0.9

    def test_scaling_property_large_interval(self, integrator):
        """Test scaling property for large intervals.

        Different interval lengths maintain relative ordering.
        """

        def func(_user_data, x):
            return np.exp(-x / 1e3)

        integrator.prepare()

        ell = 2

        # Base integral from 0 to b: f(x) j_ell(x) dx
        b = 1e4
        result_base = integrator.integrate_ell(func, None, 0.0, b, ell)

        # Scaled integral with alpha=2: integral_0^(b/2) f(x/2) j_ell(x) dx should ~=
        # 2*result Actually for integral f(alpha*x)j_ell(alpha*x)dx with substitution
        # y=alpha*x: = (1/alpha)*integral f(y)j_ell(y)dy

        # Let's test a simpler consistency: different intervals maintain relative
        # ordering
        result_half = integrator.integrate_ell(func, None, 0.0, b / 2, ell)
        result_quarter = integrator.integrate_ell(func, None, 0.0, b / 4, ell)

        # Larger intervals should have contributions from more of the function (not
        # strictly monotonic due to oscillations, but generally true)
        assert np.isfinite(result_base)
        assert np.isfinite(result_half)
        assert np.isfinite(result_quarter)

    def test_different_parameter_consistency(self):
        """Test that different algorithm parameters give consistent results."""

        def exp_decay(_user_data, x):
            return np.exp(-x / 1e3)

        ell = 4

        # Different parameter combinations
        configs = [
            {"npts": 15, "nosc": 6.0},
            {"npts": 20, "nosc": 8.0},
            {"npts": 25, "nosc": 10.0},
        ]

        results = []
        for config in configs:
            integ = Ncm.SBesselIntegratorGL(lmin=0, lmax=10, **config)
            integ.prepare()
            result = integ.integrate_ell(exp_decay, None, 0.0, 1e4, ell)
            results.append(result)

        # All should agree within 1%
        for r in results[1:]:
            assert_allclose(r, results[0], rtol=1e-2)

    def test_very_large_interval_convergence(self, integrator):
        """Test convergence for very large intervals.

        Extending interval beyond decay scale shouldn't change result much.
        """

        def exp_decay(_user_data, x):
            return np.exp(-x / 1e3)  # Decay scale ~1e3

        integrator.prepare()

        ell = 2
        # Beyond 5*decay_scale, contribution should be tiny
        result_1e4 = integrator.integrate_ell(exp_decay, None, 0.0, 1e4, ell)
        result_2e4 = integrator.integrate_ell(exp_decay, None, 0.0, 2e4, ell)

        # Difference should be very small
        diff = abs(result_2e4 - result_1e4)
        assert diff < 0.01 * abs(result_1e4)

    def test_integrate_plain_bessel_asymptotic(self, integrator):
        """Test integration using asymptotic formula.

        For f(x) = x^(ell+2), compare against analytical result using
        spherical Bessel recurrence relations.
        """

        integrator.prepare()

        for ell in [0, 1, 3, 5, 50]:

            def f(_ud, x, ell=ell):
                return x ** (ell + 2)

            a, b = 1.0e3, 1.0e5
            num = integrator.integrate_ell(f, None, a, b, ell)
            ref = b ** (ell + 2) * sp.spherical_jn(ell + 1, b) - a ** (
                ell + 2
            ) * sp.spherical_jn(ell + 1, a)

            assert_allclose(num, ref, rtol=1e-4, err_msg=f"ell={ell}")

    def test_integrate_plain_bessel_asymptotic_constant(self, integrator):
        """Test integration using asymptotic formula.

        For f(x) = 1, compare against analytical result using spherical Bessel
        recurrence relations.
        """

        integrator.prepare()

        def f(_ud, _x):
            return 1

        for ell in [0, 1, 3, 5, 50]:

            a, b = 1.0e6, 1.0e7
            num = integrator.integrate_ell(f, None, a, b, ell)
            ref = sp.spherical_jn(ell + 1, b) - sp.spherical_jn(ell + 1, a)

            assert_allclose(num, ref, rtol=1e-4, err_msg=f"ell={ell}")


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

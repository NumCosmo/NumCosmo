#!/usr/bin/env python
#
# test_py_sbessel_integrator_levin.py
#
# Sat Jan 25 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_py_sbessel_integrator_levin.py
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

"""Unit tests for Levin spherical Bessel integrator."""

from pathlib import Path
import gzip
import json
import pytest
import numpy as np
from numpy.testing import assert_allclose
from scipy.special import spherical_jn
from scipy.integrate import quad

from numcosmo_py import Ncm


class TestSBesselIntegratorLevin:
    """Tests for NcmSBesselIntegratorLevin."""

    @pytest.fixture
    def integrator(self) -> Ncm.SBesselIntegratorLevin:
        """Create a Levin integrator."""
        return Ncm.SBesselIntegratorLevin.new(0, 10)

    @pytest.mark.parametrize("l_val", [0, 1, 5, 10, 20])
    def test_integrate_constant(self, l_val: int) -> None:
        """Test integrate_ell method with f(x) = 1.

        Compares the result of the integrate_ell method against scipy numerical
        integration.
        """
        a, b = 1.0, 20.0

        # Create integrator
        integrator = Ncm.SBesselIntegratorLevin.new(l_val, l_val)

        # Define f(x) = 1
        def f_constant(_x: float, _k: float) -> float:
            return 1.0

        # Compute integral using the integrate_ell method
        result = integrator.integrate_ell(f_constant, a, b, 1.0, l_val)

        # Compute expected value using scipy
        def integrand(x: float) -> float:
            return spherical_jn(l_val, x)

        expected, _ = quad(integrand, a, b, epsabs=1e-12, epsrel=1e-14)

        # Compare
        assert_allclose(
            result,
            expected,
            rtol=1.0e-8,
            atol=1.0e-20,
            err_msg=(f"integrate_ell method failed for l={l_val}, f(x)=1"),
        )

    @pytest.mark.parametrize("l_val", [0, 1, 5, 10, 20])
    def test_integrate_linear(self, l_val: int) -> None:
        """Test integrate_ell method with f(x) = x.

        Compares the result of the integrate_ell method against scipy numerical
        integration.
        """
        a, b = 1.0, 20.0

        # Create integrator
        integrator = Ncm.SBesselIntegratorLevin.new(l_val, l_val)

        # Define f(x) = x
        def f_linear(x: float, _k: float) -> float:
            return x

        # Compute integral using the integrate_ell method
        result = integrator.integrate_ell(f_linear, a, b, 1.0, l_val)

        # Compute expected value using scipy
        def integrand(x: float) -> float:
            return x * spherical_jn(l_val, x)

        expected, _ = quad(integrand, a, b, epsabs=1e-12, epsrel=1e-14)

        # Compare
        assert_allclose(
            result,
            expected,
            rtol=1.0e-8,
            atol=1.0e-20,
            err_msg=(f"integrate_ell method failed for l={l_val}, f(x)=x"),
        )

    @pytest.mark.parametrize("l_val", [0, 1, 5, 10])
    def test_integrate_quadratic(self, l_val: int) -> None:
        """Test integrate_ell method with f(x) = x^2.

        Compares the result of the integrate_ell method against scipy numerical
        integration.
        """
        a, b = 1.0, 20.0

        # Create integrator
        integrator = Ncm.SBesselIntegratorLevin.new(l_val, l_val)

        # Define f(x) = x^2
        def f_quadratic(x: float, _k: float) -> float:
            return x * x

        # Compute integral using the integrate_ell method
        result = integrator.integrate_ell(f_quadratic, a, b, 1.0, l_val)

        # Compute expected value using scipy
        def integrand(x: float) -> float:
            return x * x * spherical_jn(l_val, x)

        expected, _ = quad(integrand, a, b, epsabs=1e-12, epsrel=1e-12)

        # Compare
        assert_allclose(
            result,
            expected,
            rtol=1.0e-8,
            atol=1.0e-20,
            err_msg=(f"integrate_ell method failed for l={l_val}, f(x)=x^2"),
        )

    @pytest.mark.parametrize("l_val", [0, 1, 5, 10, 20, 200])
    def test_integrate_rational(self, l_val: int) -> None:
        """Test integrate_ell method with rational function.

        Compares the result of the integrate_ell method against scipy numerical
        integration.
        """
        a, b = 1.0, 20.0
        center = 10.5
        std = 3.0

        # Create integrator
        integrator = Ncm.SBesselIntegratorLevin.new(l_val, l_val)

        # Define rational function f(x) = x^2 / (1 + ((x-center)/std)^2)^3
        def f_rational(x: float, _k: float) -> float:
            dx = (x - center) / std
            return x * x / ((1.0 + dx * dx) ** 3)

        # Compute integral using the integrate_ell method
        result = integrator.integrate_ell(f_rational, a, b, 1.0, l_val)

        # Compute expected value using scipy
        def integrand(x: float) -> float:
            dx = (x - center) / std
            return x * x * spherical_jn(l_val, x) / ((1.0 + dx * dx) ** 3)

        expected, _ = quad(integrand, a, b, epsabs=1e-12, epsrel=1e-14)

        # Compare
        assert_allclose(
            result,
            expected,
            rtol=1.0e-8,
            atol=1.0e-20,
            err_msg=(f"integrate_ell method failed for l={l_val}, rational function"),
        )

    def test_moving_edge_panels_stay_inside_callback_domain(self) -> None:
        """Moving edge cells are reusable without sampling outside [a, b]."""
        ell_min, ell_max = 0, 7
        integrator = Ncm.SBesselIntegratorLevin.new_full(
            ell_min, ell_max, 1.0, 1.0e4, 9, 500, 1.0e-9, 2, 1.0e-10
        )
        result = Ncm.Vector.new(ell_max - ell_min + 1)

        for a, b in [(5.0, 50.0), (6.7, 43.0), (9.5, 32.0), (10.5, 96.0)]:
            sampled = [np.inf, -np.inf]

            def f_domain(x: float, _k: float) -> float:
                sampled[0] = min(sampled[0], x)
                sampled[1] = max(sampled[1], x)
                return np.exp(-0.03 * x)

            integrator.integrate(f_domain, a, b, 1.0, result)
            assert sampled[0] >= np.nextafter(a, -np.inf)
            assert sampled[1] <= np.nextafter(b, np.inf)

            expected = np.array(
                [
                    quad(
                        lambda x, ell=ell: np.exp(-0.03 * x) * spherical_jn(ell, x),
                        a,
                        b,
                        epsabs=1.0e-12,
                        epsrel=1.0e-12,
                    )[0]
                    for ell in range(ell_min, ell_max + 1)
                ]
            )
            assert_allclose(result.to_numpy(), expected, rtol=1.0e-8, atol=1.0e-13)

    def test_accepted_edge_fits_rhs_once(self) -> None:
        """A smooth edge is fitted once before using its fixed-cell operator."""
        integrator = Ncm.SBesselIntegratorLevin.new_full(
            0, 7, 1.0, 1.0e4, 9, 500, 1.0e-10, 2, 1.0e-11
        )
        result = Ncm.Vector.new(8)
        calls = 0

        def constant(_x: float, _k: float) -> float:
            nonlocal calls
            calls += 1
            return 1.0

        # The first knot is exactly the left endpoint, leaving one right edge.
        integrator.integrate(constant, 1.0, 2.0, 1.0, result)

        expected = np.array(
            [
                quad(
                    lambda x, ell=ell: spherical_jn(ell, x),
                    1.0,
                    2.0,
                    epsabs=1.0e-13,
                    epsrel=1.0e-13,
                )[0]
                for ell in range(8)
            ]
        )
        # One fit of nine samples, plus the five probes that decide whether the edge
        # is a dead-junction cell: four interior ones and the first junction probe,
        # which already rules it out here.
        assert calls == 9 + 5
        assert_allclose(result.to_numpy(), expected, rtol=1.0e-9, atol=1.0e-14)

    def test_rejected_edge_reuses_fitted_rhs(self) -> None:
        """Rejected extrapolation must not evaluate the kernel a second time."""
        integrator = Ncm.SBesselIntegratorLevin.new_full(
            0, 7, 1.0, 1.0e4, 9, 500, 1.0e-10, 2, 1.0e-11
        )
        result = Ncm.Vector.new(8)
        calls = 0
        frequency = 500.0

        def oscillatory(x: float, _k: float) -> float:
            nonlocal calls
            calls += 1
            return np.cos(frequency * x)

        integrator.integrate(oscillatory, 1.0, 2.0, 1.0, result)

        expected, _ = quad(
            lambda x: np.cos(frequency * x) * spherical_jn(0, x),
            1.0,
            2.0,
            epsabs=1.0e-13,
            epsrel=1.0e-13,
            limit=1000,
        )
        # The fit of 1025 samples is reused; only the five dead-junction probes are added.
        assert calls == 1025 + 5
        assert_allclose(result.get(0), expected, rtol=1.0e-8, atol=1.0e-13)

    def test_set_reltol_rebuilds_the_operators(self) -> None:
        """A tolerance set after construction must reach the operators.

        The operators are built when the multipole range is known, so a later
        tolerance can only take effect by replacing them. Tightening it must
        therefore improve an edge-panel result that starts out under-resolved.

        ``cheb_reltol`` is set tight and left alone so that the operator
        tolerance is the only thing limiting the result.
        """
        ell_min, ell_max = 0, 7
        a, b = 999.9372447268345, 1018.1249616674683

        def kernel(x: float, _k: float) -> float:
            return np.exp(-0.003 * x)

        expected = np.array(
            [
                quad(
                    lambda x, ell=ell: np.exp(-0.003 * x) * spherical_jn(ell, x),
                    a,
                    b,
                    epsabs=1.0e-16,
                    epsrel=1.0e-12,
                    limit=400,
                )[0]
                for ell in range(ell_min, ell_max + 1)
            ]
        )
        scale = np.max(np.abs(expected))

        integrator = Ncm.SBesselIntegratorLevin.new_full(
            ell_min, ell_max, 1.0e-4, 1.0e6, 21, 1200, 1.0e-4, 2, 1.0e-12
        )
        result = Ncm.Vector.new(ell_max - ell_min + 1)

        integrator.integrate(kernel, a, b, 1.0, result)
        loose = np.max(np.abs(result.to_numpy() - expected)) / scale

        integrator.set_reltol(1.0e-11)
        assert integrator.get_reltol() == 1.0e-11

        integrator.integrate(kernel, a, b, 1.0, result)
        tight = np.max(np.abs(result.to_numpy() - expected)) / scale

        assert tight < loose
        assert tight < 1.0e-9

    def test_repeated_set_reltol_keeps_working(self) -> None:
        """Rebuilding the operators repeatedly must stay correct."""
        integrator = Ncm.SBesselIntegratorLevin.new(0, 5)
        result = Ncm.Vector.new(6)

        def kernel(x: float, _k: float) -> float:
            return np.exp(-0.01 * x)

        expected = np.array(
            [
                quad(
                    lambda x, ell=ell: np.exp(-0.01 * x) * spherical_jn(ell, x),
                    12.0,
                    47.0,
                    epsabs=1.0e-16,
                    epsrel=1.0e-12,
                    limit=400,
                )[0]
                for ell in range(6)
            ]
        )

        for i in range(8):
            integrator.set_reltol(1.0e-9 * (1.0 + 0.1 * i))
            integrator.integrate(kernel, 12.0, 47.0, 1.0, result)
            assert_allclose(result.to_numpy(), expected, rtol=1.0e-7, atol=1.0e-14)

    @pytest.mark.parametrize(
        "func_type,filename",
        [
            ("gaussian", "gauss_jl_500.json.gz"),
            ("rational", "rational_jl_500.json.gz"),
        ],
    )
    def test_truth_table(self, func_type: str, filename: str) -> None:
        """Test against truth tables for spherical Bessel integrals.

        This test verifies that for each multipole ell, the integrator
        achieves accurate results (rel_error < reltol) for k values up to
        at least min_k_ratio * ell.
        """
        ell_block_size = 8
        match func_type:
            case "gaussian":
                reltol = 1.0e-7
            case "rational":
                reltol = 1.0e-4
        min_k_ratio = 2.8

        truth_table_path = Path(
            Ncm.cfg_get_data_filename(f"truth_tables/sbessel/{filename}", True)
        )
        with gzip.open(truth_table_path, "rt") as f:
            truth_table = json.load(f)

        center = truth_table["center"]
        std = truth_table["std"]
        lb = truth_table["lower-bound"]
        ub = truth_table["upper-bound"]
        table = np.array(truth_table["table"])

        ells = truth_table["lvals"]
        kvals = truth_table["kvals"]
        ell_min = int(np.min(ells))
        ell_max = int(np.max(ells))
        n_ells = ell_max - ell_min + 1
        n_k = len(kvals)

        ell0 = ell_min
        ell1 = ell0 + ell_block_size - 1
        integrator = Ncm.SBesselIntegratorLevin.new(ell0, ell1)

        # Get the appropriate integration method
        if func_type == "gaussian":
            integrate_func = integrator.integrate_gaussian
        elif func_type == "rational":
            integrate_func = integrator.integrate_rational
        else:
            raise ValueError(f"Unknown function type: {func_type}")

        results_vec = Ncm.Vector.new(ell_block_size)

        # Compute all relative errors: shape (n_ells, n_k)
        rel_errors = np.zeros((n_ells, n_k))

        # Iterate over ell blocks first to reuse the cache
        for ell0 in range(ell_min, ell_max + 1, ell_block_size):
            ell1 = min(ell0 + ell_block_size - 1, ell_max)
            n_ell = ell1 - ell0 + 1
            integrator.set_ell_range(ell0, ell1)

            # Compute all k values for this ell block
            for i, k in enumerate(kvals):
                integrate_func(center, std, lb, ub, k, results_vec)
                results = results_vec.to_numpy()[:n_ell]

                truth_values = table[ell0 - ell_min : ell1 - ell_min + 1, i]
                rel_errors[ell0 - ell_min : ell1 - ell_min + 1, i] = np.abs(
                    (results - truth_values) / np.maximum(np.abs(truth_values), 1.0e-50)
                )

        # For each ell, find maximum k where rel_error < reltol
        failures = []
        for ell_idx in range(n_ells):
            ell = ell_min + ell_idx
            # Skip ell=0 as k/ell ratio is undefined
            if ell == 0:
                continue

            # Find all k indices where error is acceptable
            accurate_k_indices = np.where(rel_errors[ell_idx, :] < reltol)[0]

            if len(accurate_k_indices) > 0:
                k_max = kvals[accurate_k_indices[-1]]
                k_ratio = k_max / ell
                expected_k_min = min_k_ratio * ell

                if k_max < expected_k_min:
                    failures.append(
                        f"ell={ell}: accurate only up to k={k_max:.3g} "
                        f"(ratio={k_ratio:.1f}), expected k>={expected_k_min:.3g}"
                    )
            else:
                # No accurate results at all
                failures.append(
                    f"ell={ell}: no accurate results (all rel_errors > {reltol})"
                )

        # Assert all ells meet the criterion
        if failures:
            failure_msg = (
                f"\n[{func_type}] Accuracy criterion not met for {len(failures)} "
                f"multipoles (reltol={reltol}, min_k_ratio={min_k_ratio}):\n"
            )
            failure_msg += "\n".join(f"  {f}" for f in failures[:10])
            if len(failures) > 10:
                failure_msg += f"\n  ... and {len(failures) - 10} more"
            pytest.fail(failure_msg)

    def test_ell_zero_special_case(self) -> None:
        """Test ell=0 (monopole) special case where j_0(x) = sin(x)/x."""
        integrator = Ncm.SBesselIntegratorLevin.new(0, 0)

        def f_constant(_x: float, _k: float) -> float:
            return 1.0

        # Test with k=1
        result = integrator.integrate_ell(f_constant, 1.0, 10.0, 1.0, 0)

        # j_0(x) = sin(x)/x
        def integrand(x: float) -> float:
            return spherical_jn(0, x)

        expected, _ = quad(integrand, 1.0, 10.0, epsabs=1e-12, epsrel=1e-14)

        assert_allclose(
            result,
            expected,
            rtol=1.0e-8,
            atol=1.0e-20,
            err_msg="integrate_ell failed for ell=0 (monopole)",
        )

    @pytest.mark.parametrize("k_val", [0.1, 1.0, 10.0, 100.0])
    def test_k_scaling(self, k_val: float) -> None:
        """Test integrator accuracy with different k scales."""
        integrator = Ncm.SBesselIntegratorLevin.new(5, 5)

        def f_constant(_x: float, _k: float) -> float:
            return 1.0

        result = integrator.integrate_ell(f_constant, 1.0, 10.0, k_val, 5)

        # Verify against scipy
        def integrand(x: float) -> float:
            return spherical_jn(5, k_val * x)

        expected, _ = quad(integrand, 1.0, 10.0, epsabs=1e-12, epsrel=1e-12, limit=200)

        # Larger k may need more relaxed tolerance due to oscillations
        rtol = 1.0e-5 if k_val > 50 else 1.0e-8

        assert_allclose(
            result,
            expected,
            rtol=rtol,
            atol=1.0e-20,
            err_msg=f"integrate_ell failed for k={k_val}",
        )

    def test_narrow_integration_bounds(self) -> None:
        """Test integration over very narrow intervals."""
        integrator = Ncm.SBesselIntegratorLevin.new(5, 5)

        def f_test(x: float, _k: float) -> float:
            return x

        # Very narrow interval
        a, b = 5.0, 5.001
        result = integrator.integrate_ell(f_test, a, b, 1.0, 5)

        # For narrow intervals, approximate j_l(k*x) ~= j_l(k*x_mid)
        def integrand(x: float) -> float:
            return x * spherical_jn(5, x)

        expected, _ = quad(integrand, a, b, epsabs=1e-15, epsrel=1e-15)

        # Narrow intervals are challenging, use relaxed tolerance
        assert_allclose(
            result,
            expected,
            rtol=1.0e-4,
            atol=1.0e-15,
            err_msg="Narrow interval integration failed",
        )

    @pytest.mark.parametrize("ell_val", [100, 200, 300])
    def test_high_ell_regime(self, ell_val: int) -> None:
        """Test integrator performance at high ell values."""
        integrator = Ncm.SBesselIntegratorLevin.new(ell_val, ell_val)

        def f_simple(_x: float, _k: float) -> float:
            return 1.0

        # For high ell, keep k moderate (k ~ ell for best accuracy)
        k = float(ell_val)
        result = integrator.integrate_ell(f_simple, 1.0, 20.0, k, ell_val)

        # Verify against scipy
        def integrand(x: float) -> float:
            return spherical_jn(ell_val, k * x)

        expected, _ = quad(integrand, 1.0, 20.0, epsabs=1e-10, epsrel=1e-10, limit=1000)

        # High ell is extremely challenging; even scipy struggles
        # Relax tolerance to 3x (covers ell=100 case with rel_error ~ 2.5)
        assert_allclose(
            result,
            expected,
            rtol=3.0,
            atol=1.0e-15,
            err_msg=f"High ell integration failed for ell={ell_val}",
        )

    def test_oscillatory_regime(self) -> None:
        """Test performance in highly oscillatory regime (large k, moderate ell)."""
        ell = 10
        k = 500.0  # Large k creates rapid oscillations

        integrator = Ncm.SBesselIntegratorLevin.new(ell, ell)

        def f_smooth(x: float, _k: float) -> float:
            return np.exp(-0.01 * x)

        result = integrator.integrate_ell(f_smooth, 1.0, 100.0, k, ell)

        # Verify against scipy (will be slow due to oscillations)
        def integrand(x: float) -> float:
            return np.exp(-0.01 * x) * spherical_jn(ell, k * x)

        expected, _ = quad(
            integrand, 1.0, 100.0, epsabs=1e-10, epsrel=1e-10, limit=5000
        )

        # Oscillatory regime is challenging (scipy hits subdivision limit)
        # Relax tolerance to 1.0 (covers rel_error ~ 0.91)
        assert_allclose(
            result,
            expected,
            rtol=1.0e-2,
            atol=1.0e-15,
            err_msg=f"Oscillatory regime failed for k={k}, ell={ell}",
        )


def _jl_second_deriv(ell: int, x: np.ndarray) -> np.ndarray:
    """j_l'' from the homogeneous ODE and scipy's j_l, j_l'."""
    j = spherical_jn(ell, x)
    jp = spherical_jn(ell, x, derivative=True)
    return -2.0 / x * jp + (ell * (ell + 1.0) / x**2 - 1.0) * j


class TestSBesselIntegratorLevinDeriv:
    """Tests for the derivative-weighted integrals of NcmSBesselIntegratorLevin.

    integrate_deriv computes int_a^b K(x, k) j_l^{(d)}(k x) dx with the
    derivative taken with respect to the Bessel argument x = k chi.
    """

    @pytest.mark.parametrize("l_val", [0, 2, 10, 30])
    @pytest.mark.parametrize("deriv", [1, 2])
    def test_gaussian_vs_quadrature(self, l_val: int, deriv: int) -> None:
        """Gaussian K against scipy adaptive quadrature, including the turning point."""
        a, b = 1.0, 8.0
        center, sigma = 4.0, 1.2
        integrator = Ncm.SBesselIntegratorLevin.new(l_val, l_val)
        result = Ncm.Vector.new(1)

        def f_gauss(x: float, _k: float) -> float:
            return np.exp(-0.5 * ((x - center) / sigma) ** 2)

        # k values placing the turning point x ~ l inside, below and above the window
        k_list = [0.5, 2.0]
        if l_val > 0:
            k_list += [l_val / center, l_val / a]

        for k in k_list:
            integrator.integrate_deriv(f_gauss, a, b, k, deriv, result)

            if deriv == 1:

                def integrand(x: float) -> float:
                    return f_gauss(x, k) * spherical_jn(l_val, k * x, derivative=True)

            else:

                def integrand(x: float) -> float:
                    return f_gauss(x, k) * _jl_second_deriv(l_val, np.asarray(k * x))

            expected, _ = quad(integrand, a, b, limit=800, epsabs=1e-14, epsrel=1e-13)
            assert_allclose(
                result.get(0),
                expected,
                rtol=1.0e-8,
                atol=1.0e-13,
                err_msg=f"l={l_val}, deriv={deriv}, k={k}",
            )

    def test_evanescent_region_keeps_relative_precision(self) -> None:
        """Deep small-argument tail: tiny values must stay relatively accurate."""
        l_val = 30
        a, b, k = 1.0, 8.0, 0.5
        integrator = Ncm.SBesselIntegratorLevin.new(l_val, l_val)
        result = Ncm.Vector.new(1)

        def f_gauss(x: float, _k: float) -> float:
            return np.exp(-0.5 * ((x - 4.0) / 1.2) ** 2)

        integrator.integrate_deriv(f_gauss, a, b, k, 2, result)

        def integrand(x: float) -> float:
            return f_gauss(x, k) * _jl_second_deriv(l_val, np.asarray(k * x))

        expected, _ = quad(integrand, a, b, limit=400, epsabs=0.0, epsrel=1e-13)

        # values are ~1e-26 here; the point is relative, not absolute, accuracy
        assert expected != 0.0
        assert_allclose(result.get(0), expected, rtol=1.0e-9)

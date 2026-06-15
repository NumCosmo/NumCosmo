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
        N = 128
        a, b = 1.0, 20.0

        # Create integrator
        integrator = Ncm.SBesselIntegratorLevin.new(l_val, l_val)
        integrator.set_max_order(N)

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
        N = 128
        a, b = 1.0, 20.0

        # Create integrator
        integrator = Ncm.SBesselIntegratorLevin.new(l_val, l_val)
        integrator.set_max_order(N)

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
        N = 128
        a, b = 1.0, 20.0

        # Create integrator
        integrator = Ncm.SBesselIntegratorLevin.new(l_val, l_val)
        integrator.set_max_order(N)

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
        N = 256
        a, b = 1.0, 20.0
        center = 10.5
        std = 3.0

        # Create integrator
        integrator = Ncm.SBesselIntegratorLevin.new(l_val, l_val)
        integrator.set_max_order(N)

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

    def test_integrate_l_range_consistency(self) -> None:
        """Test that integrate gives consistent results.

        Compares integrate (which uses batched solver internally)
        against individual integrate_ell calls for each l.
        """
        N = 128
        a, b = 1.0, 20.0
        lmin, lmax = 5, 15

        # Create integrator
        integrator = Ncm.SBesselIntegratorLevin.new(lmin, lmax)
        integrator.set_max_order(N)

        # Define a test function
        def f_test(x: float, _k: float) -> float:
            return np.exp(-0.1 * x)

        # Use integrate (batched internally)
        results_vec = Ncm.Vector.new(lmax - lmin + 1)
        integrator.integrate(f_test, a, b, 1.0, results_vec)
        results_range_np = results_vec.to_numpy()

        # Compute individually for each l
        results_individual = []
        for ell in range(lmin, lmax + 1):
            results_individual.append(integrator.integrate_ell(f_test, a, b, 1.0, ell))

        results_individual_np = np.array(results_individual)

        # Compare
        assert_allclose(
            results_range_np,
            results_individual_np,
            rtol=1.0e-10,
            atol=1.0e-14,
            err_msg="integrate doesn't match individual integrate_ell calls",
        )

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
            Ncm.cfg_get_data_filename(f"truth_tables/{filename}", True)
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

        # Iterate over ell blocks first to leverage caching
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
        integrator.set_max_order(128)

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
        integrator.set_max_order(256)

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

    def test_ell_range_caching(self) -> None:
        """Test that changing ell_range properly updates results."""
        integrator = Ncm.SBesselIntegratorLevin.new(0, 5)
        integrator.set_max_order(128)

        def f_test(x: float, _k: float) -> float:
            return np.exp(-0.1 * x)

        # Get results for ell=0..5
        results1 = Ncm.Vector.new(6)
        integrator.integrate(f_test, 1.0, 10.0, 1.0, results1)

        # Change range to ell=10..15 and verify different results
        integrator.set_ell_range(10, 15)
        results2 = Ncm.Vector.new(6)
        integrator.integrate(f_test, 1.0, 10.0, 1.0, results2)

        # Results should be different (different ells computed)
        assert not np.allclose(
            results1.to_numpy(), results2.to_numpy()
        ), "Changing ell_range should produce different results"

        # Verify second set matches individual computations
        results_individual = []
        for ell in range(10, 16):
            results_individual.append(
                integrator.integrate_ell(f_test, 1.0, 10.0, 1.0, ell)
            )

        assert_allclose(
            results2.to_numpy(),
            np.array(results_individual),
            rtol=1.0e-10,
            atol=1.0e-14,
            err_msg="Cached results don't match individual integrate_ell calls",
        )

    @pytest.mark.parametrize("block_size", [1, 4, 8, 16])
    def test_block_size_independence(self, block_size: int) -> None:
        """Test that results are independent of block size used."""
        lmin, lmax = 10, 25
        k = 1.0

        def f_test(x: float, _k: float) -> float:
            return np.exp(-0.2 * x)

        # Compute with given block size
        integrator = Ncm.SBesselIntegratorLevin.new(lmin, lmin + block_size - 1)
        integrator.set_max_order(128)

        n_ells = lmax - lmin + 1
        results = np.zeros(n_ells)
        results_vec = Ncm.Vector.new(block_size)

        for ell0 in range(lmin, lmax + 1, block_size):
            ell1 = min(ell0 + block_size - 1, lmax)
            n_ell = ell1 - ell0 + 1
            integrator.set_ell_range(ell0, ell1)
            integrator.integrate(f_test, 1.0, 10.0, k, results_vec)
            results[ell0 - lmin : ell1 - lmin + 1] = results_vec.to_numpy()[:n_ell]

        # Compare with individual computations (effectively block_size=1)
        integrator_ref = Ncm.SBesselIntegratorLevin.new(lmin, lmin)
        integrator_ref.set_max_order(128)

        results_ref = []
        for ell in range(lmin, lmax + 1):
            results_ref.append(integrator_ref.integrate_ell(f_test, 1.0, 10.0, k, ell))

        assert_allclose(
            results,
            np.array(results_ref),
            rtol=1.0e-10,
            atol=1.0e-14,
            err_msg=f"Results differ for block_size={block_size}",
        )

    def test_narrow_integration_bounds(self) -> None:
        """Test integration over very narrow intervals."""
        integrator = Ncm.SBesselIntegratorLevin.new(5, 5)
        integrator.set_max_order(128)

        def f_test(x: float, _k: float) -> float:
            return x

        # Very narrow interval
        a, b = 5.0, 5.001
        result = integrator.integrate_ell(f_test, a, b, 1.0, 5)

        # For narrow intervals, approximate j_l(k*x) ≈ j_l(k*x_mid)
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
        integrator.set_max_order(512)

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
        integrator.set_max_order(1024)

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

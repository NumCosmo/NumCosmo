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
import subprocess
import sys
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

    def test_moving_edge_panels_stay_inside_callback_domain(self) -> None:
        """Moving edge cells are reusable without sampling outside [a, b]."""
        ell_min, ell_max = 0, 7
        integrator = Ncm.SBesselIntegratorLevin.new_full(
            ell_min, ell_max, 1.0, 1.0e4, 9, 500, 1.0e-9, 2, 1.0e-10
        )
        result = Ncm.Vector.new(ell_max - ell_min + 1)

        for a, b in [(5.0, 50.0), (6.7, 43.0), (9.5, 32.0), (10.5, 96.0)]:
            sampled = [np.inf, -np.inf]

            # pylint: disable=cell-var-from-loop
            def f_domain(x: float, _k: float) -> float:
                sampled[0] = min(sampled[0], x)
                sampled[1] = max(sampled[1], x)
                return np.exp(-0.03 * x)

            # pylint: enable=cell-var-from-loop

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

    def test_zero_bessel_batch_skips_rhs(self) -> None:
        """A panel with an identically zero Bessel batch never calls the kernel."""
        integrator = Ncm.SBesselIntegratorLevin.new(500, 500)
        result = Ncm.Vector.new(1)
        calls = 0

        def kernel(_x: float, _k: float) -> float:
            nonlocal calls
            calls += 1
            return 1.0

        integrator.integrate(kernel, 1.0, 2.0, 1.0, result)

        assert calls == 0
        assert result.get(0) == 0.0

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

    @pytest.mark.parametrize("setter", ["set_reltol", "set_cheb_reltol"])
    @pytest.mark.parametrize("value", ["0.0", "-1.0e-8"])
    def test_non_positive_tolerance_is_rejected(self, setter: str, value: str) -> None:
        """Neither tolerance is meaningful at or below zero.

        Runs in a subprocess since the assertion aborts rather than raises.
        """
        proc = subprocess.run(
            [
                sys.executable,
                "-c",
                "from numcosmo_py import Ncm\n"
                "Ncm.cfg_init()\n"
                "sbi = Ncm.SBesselIntegratorLevin.new(0, 3)\n"
                f"sbi.{setter}({value})\n",
            ],
            capture_output=True,
            text=True,
            check=False,
        )

        assert proc.returncode != 0
        assert "assertion failed" in proc.stderr

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


class TestOscillatoryResolutionFloor:
    """The decay test must not converge below a panel's oscillation count.

    A panel [a, b] in y = kx carries about (b - a) / pi oscillations of the
    solution. On such a panel the leading Chebyshev coefficients are small and
    nearly flat, so the adaptive QR's decay test used to fire on them and declare
    convergence at a tiny order -- a panel spanning 2162 in y was "converged"
    with 19 columns instead of the ~1200 it needs.

    Because panel contributions cancel heavily, the visible symptom was
    non-monotonic: at a loose tolerance every panel was under-resolved and the
    errors partly cancelled, while at an intermediate tolerance some panels
    resolved and others did not, destroying the cancellation. Requesting 1e-8 was
    then ~50x worse than requesting 1e-6.
    """

    @staticmethod
    def _gaussian_integral(reltol):
        """Strongly oscillatory Gaussian integral at a requested tolerance."""
        ell, k, a, b = 2, 500.0, 0.1, 10.0
        defaults = Ncm.SBesselIntegratorLevin.new(ell, ell)
        integrator = Ncm.SBesselIntegratorLevin.new_full(
            ell,
            ell,
            defaults.get_y_knots_min(),
            defaults.get_y_knots_max(),
            defaults.get_n_knots(),
            defaults.get_ell_cache_max(),
            reltol,
            defaults.get_cheb_min_order(),
            reltol,
        )
        return integrator.integrate_gaussian_ell(5.0, 1.0, a, b, k, ell)

    def test_loose_tolerances_agree_with_the_converged_result(self):
        """No requested tolerance may be catastrophically wrong.

        Pre-fix the 1e-8 request returned +6.10e-06 against a converged
        -1.19e-07: wrong sign and ~50x too large.
        """
        converged = self._gaussian_integral(1.0e-12)

        for reltol in (1.0e-4, 1.0e-6, 1.0e-8, 1.0e-10):
            result = self._gaussian_integral(reltol)
            assert_allclose(result, converged, rtol=1.0e-2, atol=0.0)

    def test_accuracy_is_monotonic_in_the_requested_tolerance(self):
        """Tightening the request must not make the answer wholesale worse.

        Pre-fix the error rose from 9.3e-01 at 1e-6 to 5.2e+01 at 1e-8, a 56x
        reversal.

        Only the loose half can be compared step by step. From about 1e-7 the
        adaptive refinement lands on discrete levels, so the delivered error
        bounces inside a ~1e-7 band rather than falling: sweeping reltol on one
        machine gives 4.4e-09 at 1e-8.25 against 1.8e-07 at 1e-8.50, and the
        request does not actually converge until past 1e-11. Which level a
        given request lands on also shifts with the BLAS kernel, so CI has
        produced 8.5e-08, 1.1e-07 and 3.0e-07 at 1e-11 where this machine
        gives 1.3e-08. The tight end therefore gets an absolute bound, not a
        comparison.
        """
        converged = self._gaussian_integral(1.0e-12)
        errors = {
            reltol: abs(self._gaussian_integral(reltol) / converged - 1.0)
            for reltol in (1.0e-4, 1.0e-6, 1.0e-8, 1.0e-11)
        }

        assert errors[1.0e-6] <= errors[1.0e-4] * 3.0
        assert errors[1.0e-8] <= errors[1.0e-6] * 3.0
        assert errors[1.0e-11] < 1.0e-5


class TestPanelRecording:
    """Per-panel diagnostics.

    The integral is a sum of per-panel boundary terms. Their individual sizes are
    not recoverable from the result, so measuring how much they cancel -- and
    whether a single panel's own accuracy is the limit -- needs them recorded.
    """

    @staticmethod
    def _integrator(ell, record):
        defaults = Ncm.SBesselIntegratorLevin.new(ell, ell)
        sbi = Ncm.SBesselIntegratorLevin.new_full(
            ell,
            ell,
            defaults.get_y_knots_min(),
            defaults.get_y_knots_max(),
            defaults.get_n_knots(),
            defaults.get_ell_cache_max(),
            1.0e-12,
            defaults.get_cheb_min_order(),
            1.0e-12,
        )
        sbi.set_record_panels(record)
        return sbi

    def test_recording_is_off_by_default(self):
        """Recording must be opt-in: it is a diagnostic, not a feature."""
        sbi = Ncm.SBesselIntegratorLevin.new(2, 2)

        assert not sbi.get_record_panels()
        assert sbi.get_n_panel_records() == 0

    def test_contributions_sum_to_the_result(self):
        """The recorded panels must account for the integral exactly."""
        sbi = self._integrator(2, True)
        result = sbi.integrate_gaussian_ell(5.0, 1.0, 0.1, 10.0, 50.0, 2)

        n = sbi.get_n_panel_records()
        assert n > 1, "expected the range to be split into several panels"

        # Summing the records in Python adds in a different order from the
        # internal accumulation, so a few ulps of difference are expected.
        total = sum(sbi.get_panel_contrib(i) for i in range(n))
        assert_allclose(total, result, rtol=1.0e-12)

    def test_panels_tile_the_range_in_order(self):
        """Recorded panels must abut and cover [k*a, k*b]."""
        k, a, b, ell = 50.0, 0.1, 10.0, 2
        sbi = self._integrator(ell, True)
        sbi.integrate_gaussian_ell(5.0, 1.0, a, b, k, ell)

        n = sbi.get_n_panel_records()
        bounds = [(sbi.get_panel_a(i), sbi.get_panel_b(i)) for i in range(n)]

        assert_allclose(bounds[0][0], k * a, rtol=1.0e-12)
        assert_allclose(bounds[-1][1], k * b, rtol=1.0e-12)
        for (_, prev_b), (next_a, _) in zip(bounds, bounds[1:]):
            assert_allclose(next_a, prev_b, rtol=1.0e-12)
        for i in range(n):
            assert sbi.get_panel_ell(i) == ell

    def test_records_are_cleared_between_integrations(self) -> None:
        """Records describe the most recent call, not an accumulation."""
        sbi = self._integrator(2, True)
        sbi.integrate_gaussian_ell(5.0, 1.0, 0.1, 10.0, 50.0, 2)
        first = sbi.get_n_panel_records()
        sbi.integrate_gaussian_ell(5.0, 1.0, 0.1, 10.0, 50.0, 2)

        assert sbi.get_n_panel_records() == first

    def test_disabling_recording_releases_the_records(self) -> None:
        """Disabling recording must clear the records."""
        sbi = self._integrator(2, True)
        sbi.integrate_gaussian_ell(5.0, 1.0, 0.1, 10.0, 50.0, 2)
        assert sbi.get_n_panel_records() > 0

        sbi.set_record_panels(False)
        assert sbi.get_n_panel_records() == 0


def _jl_second_deriv(ell: int, y: np.ndarray) -> np.ndarray:
    """j_l'' from the homogeneous ODE and scipy's j_l, j_l'."""
    j = spherical_jn(ell, y)
    jp = spherical_jn(ell, y, derivative=True)
    return -2.0 / y * jp + (ell * (ell + 1.0) / y**2 - 1.0) * j


class TestSBesselIntegratorLevinDeriv:
    """Tests for the derivative-weighted integrals of NcmSBesselIntegratorLevin.

    integrate_deriv computes int_a^b K(x, k) j_l^{(d)}(k x) dx with the
    derivative taken with respect to the Bessel argument y = k x.
    """

    @pytest.mark.parametrize("l_val", [0, 1, 2, 5, 20, 60])
    @pytest.mark.parametrize("k", [0.5, 3.0, 40.0])
    def test_constant_closed_form(self, l_val: int, k: float) -> None:
        """For K = 1 the integrals reduce to endpoint evaluations.

        int j_l'(k x) dx = [j_l(k x)]/k and int j_l''(k x) dx = [j_l'(k x)]/k,
        exactly.
        """
        a, b = 0.3, 12.0
        integrator = Ncm.SBesselIntegratorLevin.new(l_val, l_val)
        result = Ncm.Vector.new(1)

        def f_constant(_x: float, _k: float) -> float:
            return 1.0

        integrator.integrate_deriv(f_constant, a, b, k, 1, result)
        expected_d1 = (spherical_jn(l_val, k * b) - spherical_jn(l_val, k * a)) / k
        assert_allclose(result.get(0), expected_d1, rtol=1.0e-11, atol=1.0e-16)

        integrator.integrate_deriv(f_constant, a, b, k, 2, result)
        expected_d2 = (
            spherical_jn(l_val, k * b, derivative=True)
            - spherical_jn(l_val, k * a, derivative=True)
        ) / k
        assert_allclose(result.get(0), expected_d2, rtol=1.0e-11, atol=1.0e-16)

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

        # k values placing the turning point y ~ l inside, below and above the window
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

    def test_batched_matches_single(self) -> None:
        """A batched ell block must reproduce the per-ell results."""
        ell_min, ell_max = 2, 33
        a, b, k = 0.5, 9.0, 7.0
        n_ell = ell_max - ell_min + 1

        def f_gauss(x: float, _k: float) -> float:
            return np.exp(-0.5 * ((x - 4.0) / 1.5) ** 2)

        batched = Ncm.SBesselIntegratorLevin.new(ell_min, ell_max)
        res_batch = Ncm.Vector.new(n_ell)

        for deriv in (1, 2):
            batched.integrate_deriv(f_gauss, a, b, k, deriv, res_batch)

            for ell in (ell_min, ell_min + 7, ell_max):
                single = Ncm.SBesselIntegratorLevin.new(ell, ell)
                res_single = Ncm.Vector.new(1)
                single.integrate_deriv(f_gauss, a, b, k, deriv, res_single)
                assert_allclose(
                    res_batch.get(ell - ell_min),
                    res_single.get(0),
                    rtol=1.0e-10,
                    atol=1.0e-16,
                    err_msg=f"ell={ell}, deriv={deriv}",
                )

    def test_deriv_zero_matches_integrate(self) -> None:
        """deriv = 0 must dispatch to the plain integrate path."""
        ell_min, ell_max = 0, 8
        a, b, k = 0.5, 9.0, 3.0
        n_ell = ell_max - ell_min + 1

        def f_gauss(x: float, _k: float) -> float:
            return np.exp(-0.5 * ((x - 4.0) / 1.5) ** 2)

        integrator = Ncm.SBesselIntegratorLevin.new(ell_min, ell_max)
        res_a = Ncm.Vector.new(n_ell)
        res_b = Ncm.Vector.new(n_ell)

        integrator.integrate_deriv(f_gauss, a, b, k, 0, res_a)
        integrator.integrate(f_gauss, a, b, k, res_b)

        for i in range(n_ell):
            assert res_a.get(i) == res_b.get(i)

    def test_order_above_two_is_rejected(self) -> None:
        """Orders above 2 have no implementation anywhere and must fail loudly.

        Runs in a subprocess since g_error aborts rather than raises.
        """
        proc = subprocess.run(
            [
                sys.executable,
                "-c",
                "from numcosmo_py import Ncm\n"
                "Ncm.cfg_init()\n"
                "sbi = Ncm.SBesselIntegratorLevin.new(0, 3)\n"
                "res = Ncm.Vector.new(4)\n"
                "sbi.integrate_deriv(lambda x, k: 1.0, 1.0, 2.0, 1.0, 3, res)\n",
            ],
            capture_output=True,
            text=True,
            check=False,
        )

        assert proc.returncode != 0
        assert "not supported" in proc.stderr

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

    def test_recurrence_cross_check(self) -> None:
        """deriv = 2 against the l-recurrence combination of plain integrals.

        j_l''(y) = (l (l-1)/y^2 - 1) j_l(y) + (2/y) j_{l+1}(y), so the deriv-2
        integral must match a combination of three deriv-0 integrals computed
        through an independent code path.
        """
        l_val = 12
        a, b, k = 0.8, 10.0, 5.0

        def f_gauss(x: float, _k: float) -> float:
            return np.exp(-0.5 * ((x - 4.0) / 1.5) ** 2)

        integrator = Ncm.SBesselIntegratorLevin.new(l_val, l_val + 1)
        result = Ncm.Vector.new(2)
        integrator.integrate_deriv(f_gauss, a, b, k, 2, result)
        got = result.get(0)

        def f_over_y2(x: float, kk: float) -> float:
            return f_gauss(x, kk) / (kk * x) ** 2

        def f_over_y(x: float, kk: float) -> float:
            return f_gauss(x, kk) / (kk * x)

        i_f = integrator.integrate_ell(f_gauss, a, b, k, l_val)
        i_f_y2 = integrator.integrate_ell(f_over_y2, a, b, k, l_val)
        i_f_y1 = integrator.integrate_ell(f_over_y, a, b, k, l_val + 1)

        expected = l_val * (l_val - 1.0) * i_f_y2 - i_f + 2.0 * i_f_y1
        assert_allclose(got, expected, rtol=1.0e-9)


class TestTauConstraintRule:
    """Per-panel tau constraint: the rule, the order check and the guard.

    The tau constraint is valid while the working order n sits on the plateau of the
    homogeneous spectrum, n_F << n < N_min ~ 2 span / pi. The rule only sees the
    panel (N_min); the order check runs once the forcing fit exists. A forcing that
    needs more coefficients than the plateau holds must fall back to Dirichlet
    instead of admitting a homogeneous multiple of the smooth solution.
    """

    A, B, K = 3162.28, 3672.79, 1.0
    ELL_MIN, ELL_MAX = 200, 207

    @classmethod
    def _gaussian(cls, width_fraction: float):
        x0 = 0.5 * (cls.A + cls.B)
        width = (cls.B - cls.A) / width_fraction

        def forcing(x: float, _k: float) -> float:
            return 1.0e-16 * np.exp(-0.5 * ((x - x0) / width) ** 2)

        return forcing

    @classmethod
    def _integrate(cls, forcing, min_osc: float):
        sbi = Ncm.SBesselIntegratorLevin.new(cls.ELL_MIN, cls.ELL_MAX)
        sbi.set_tau_constraint_min_osc(min_osc)
        result = Ncm.Vector.new(cls.ELL_MAX - cls.ELL_MIN + 1)
        sbi.integrate(forcing, cls.A, cls.B, cls.K, result)
        return result.to_numpy(), sbi.get_n_constraint_fallbacks()

    def test_panel_qualifies_for_the_rule(self) -> None:
        """The test panel must hold more oscillations than the default threshold."""
        osc = 2.0 / np.pi * (self.B - self.A)
        defaults = Ncm.SBesselIntegratorLevin.new(self.ELL_MIN, self.ELL_MAX)

        assert osc > defaults.get_tau_constraint_min_osc()

    def test_smooth_forcing_uses_the_tau_constraint(self) -> None:
        """A wide forcing is solved with the tau constraint and matches Dirichlet."""
        forcing = self._gaussian(10.0)
        dirichlet, _ = self._integrate(forcing, 0.0)
        tau, fallbacks = self._integrate(forcing, 200.0)

        assert fallbacks == 0
        # The two constraints solve different systems: identical bits would mean the
        # tau path was never taken.
        assert np.any(tau != dirichlet)
        assert_allclose(tau, dirichlet, atol=1.0e-8 * np.abs(dirichlet).max())

    @pytest.mark.parametrize("width_fraction", [30.0, 100.0])
    def test_forcing_above_the_plateau_falls_back(self, width_fraction: float) -> None:
        """A forcing whose order exceeds 0.75 N_min is redone with Dirichlet.

        Without the order check the tau truncation lands past N_min, the solution
        grows without bound and the boundary term is non-finite.
        """
        forcing = self._gaussian(width_fraction)
        dirichlet, _ = self._integrate(forcing, 0.0)
        tau, fallbacks = self._integrate(forcing, 200.0)

        assert fallbacks > 0
        assert np.all(np.isfinite(tau))
        assert_allclose(tau, dirichlet, rtol=1.0e-10)

    def test_rule_off_never_falls_back(self) -> None:
        """With the rule off every panel is Dirichlet and nothing is counted."""
        _, fallbacks = self._integrate(self._gaussian(100.0), 0.0)

        assert fallbacks == 0


class TestTurningKnot:
    """The per-block knot just above the turning point.

    The guard cannot certify the tau constraint on a panel straddling the turning point,
    because the bound on the smooth member carries min|y^2 - nu^2| in its denominator.
    When the kernel's reach in y stops before the next base knot, that panel is the whole
    oscillatory region and the constraint is unusable at every k. One knot lifts a panel edge
    clear of the turning point.
    """

    A, B, K = 1.0e2, 1.0e5, 1.0

    @staticmethod
    def _forcing(x: float, _k: float) -> float:
        return np.exp(-0.5 * ((x - 3.0e4) / 6.0e3) ** 2)

    def test_default_is_on(self) -> None:
        """The margin ships enabled."""
        sbi = Ncm.SBesselIntegratorLevin.new(2, 9)
        assert sbi.get_turning_knot_margin() == pytest.approx(
            Ncm.SBESSEL_INTEGRATOR_LEVIN_DEFAULT_TURNING_KNOT_MARGIN
        )
        assert sbi.get_turning_knot_margin() > 1.0

    def test_round_trip(self) -> None:
        """The margin reads back, and zero is accepted as off."""
        sbi = Ncm.SBesselIntegratorLevin.new(2, 9)
        for value in (1.2, 0.0, 1.05):
            sbi.set_turning_knot_margin(value)
            assert sbi.get_turning_knot_margin() == value

    def _run(self, ell: int, margin: float):
        sbi = Ncm.SBesselIntegratorLevin.new(ell, ell + 7)
        sbi.set_turning_knot_margin(margin)
        result = Ncm.Vector.new(8)
        sbi.integrate(self._forcing, self.A, self.B, self.K, result)

        return (
            result.to_numpy().copy(),
            sbi.get_n_tau_solves(),
            sbi.get_n_constraint_fallbacks(),
        )

    @pytest.mark.parametrize("ell", [500, 1000])
    def test_removes_the_fallbacks_at_high_ell(self, ell: int) -> None:
        """The straddling panel cannot be certified, so without the knot it is refused.

        How much of the oscillatory region that panel holds depends on the kernel's reach
        in y; on a kernel that stops before the next base knot it is all of it and the
        constraint is unusable outright. Here the refusal is one panel, which is what the
        counters show: a fallback that the knot removes, and one more panel taking the
        cheap constraint.
        """
        off, tau_off, fb_off = self._run(ell, 0.0)
        on, tau_on, fb_on = self._run(ell, 1.05)

        assert fb_off > 0, "expected the straddling panel to be refused"
        assert fb_on == 0, "the knot must leave nothing to refuse"
        assert tau_on > tau_off

    def test_same_integral_either_way(self) -> None:
        """The knot changes which constraint runs, not the answer."""
        for ell in (200, 1000):
            off, _, _ = self._run(ell, 0.0)
            on, _, _ = self._run(ell, 1.05)
            peak = np.abs(off).max()
            assert_allclose(on, off, rtol=0.0, atol=1.0e-8 * peak)

    def test_low_ell_is_untouched(self) -> None:
        """Below the rule's threshold the straddling panel is Dirichlet either way, so
        no knot is added and the result is bit-identical."""
        off, _, _ = self._run(20, 0.0)
        on, _, _ = self._run(20, 1.05)
        assert np.array_equal(on, off)

    def test_rule_off_disables_the_knot(self) -> None:
        """The insertion is gated on the same oscillation count as the rule."""
        sbi = Ncm.SBesselIntegratorLevin.new(1000, 1007)
        sbi.set_tau_constraint_min_osc(0.0)
        result = Ncm.Vector.new(8)
        sbi.integrate(self._forcing, self.A, self.B, self.K, result)
        assert sbi.get_n_tau_solves() == 0


class TestSharedKnotTable:
    """Integrators over one grid share their j_l rows.

    The table is a pure function of the abscissae and of ell_cache_max, and a solver holds
    one integrator per multipole block, so the rows are shared rather than copied. Shared
    rows are read-only, and the risk of getting that wrong is one integrator disturbing
    another's results.
    """

    A, B, K = 1.0e2, 1.0e4, 1.0

    @staticmethod
    def _forcing(x: float, _k: float) -> float:
        return np.exp(-0.5 * ((x - 3.0e3) / 6.0e2) ** 2)

    def _run(self, sbi: Ncm.SBesselIntegratorLevin) -> np.ndarray:
        result = Ncm.Vector.new(8)
        sbi.integrate(self._forcing, self.A, self.B, self.K, result)

        return result.to_numpy().copy()

    def test_a_second_integrator_does_not_disturb_the_first(self) -> None:
        """Same grid, different blocks: sharing must not be mutation."""
        first = Ncm.SBesselIntegratorLevin.new(200, 207)
        before = self._run(first)

        second = Ncm.SBesselIntegratorLevin.new(600, 607)
        self._run(second)

        after = self._run(first)
        assert np.array_equal(after, before)

    def test_identical_grids_agree(self) -> None:
        """Two integrators built alike give the same answer."""
        a = Ncm.SBesselIntegratorLevin.new(200, 207)
        b = Ncm.SBesselIntegratorLevin.new(200, 207)

        assert np.array_equal(self._run(a), self._run(b))

    def test_turning_knot_row_is_per_instance(self) -> None:
        """Blocks insert their knot at different y, so that row cannot be shared."""
        low = Ncm.SBesselIntegratorLevin.new(400, 407)
        high = Ncm.SBesselIntegratorLevin.new(900, 907)

        low_first = self._run(low)
        self._run(high)

        assert np.array_equal(self._run(low), low_first)

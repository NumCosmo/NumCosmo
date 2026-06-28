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
from pathlib import Path
import gzip
import json
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
        ell_min, ell_max = integrator.get_ell_range()
        assert ell_min == 0
        assert ell_max == 10

    def test_integrate_constant_function(
        self, integrator: Ncm.SBesselIntegratorFFTL
    ) -> None:
        """Test integration of f(x) = 1."""

        def constant_one(_x, _k) -> float:
            return 1.0

        # Test for different multipoles
        for ell in [0, 1, 2, 3, 4, 5, 100]:
            result = integrator.integrate_ell(constant_one, 0.0, 2000.0, 1.0, ell)
            assert np.isfinite(result)
            # Result should be non-zero for integration of j_ell from 0 to 20
            assert abs(result) > 1e-10

    def test_integrate_vs_integrate_ell(
        self, integrator: Ncm.SBesselIntegratorFFTL
    ) -> None:
        """Test that integrate method agrees with integrate_ell for each ell."""
        center = 100.0
        std = 10.0

        def test_func(x: float, _k: float) -> float:
            return x**2 / (1 + ((x - center) / std) ** 2) ** 2

        a = center - 20.0 * std
        b = center + 20.0 * std

        integrator.set_ell_range(0, 600)

        # Get results from integrate_ell for each ell
        results_ell = []
        ell_min, ell_max = integrator.get_ell_range()
        for ell in range(ell_min, ell_max + 1):
            result = integrator.integrate_ell(test_func, a, b, 1.0, ell)
            results_ell.append(result)

        # Get results from integrate (all ells at once)
        ell_min, ell_max = integrator.get_ell_range()
        n_ell = ell_max - ell_min + 1
        results_vec = Ncm.Vector.new(n_ell)
        integrator.integrate(test_func, a, b, 1.0, results_vec)
        results_all = np.array([results_vec.get(i) for i in range(n_ell)])
        max_result = np.max(np.abs(results_all))

        # Compare results
        assert_allclose(
            results_all / max_result,
            results_ell / max_result,
            rtol=1e-5,
            atol=max_result * 1e-5,
        )

    def test_integrate_performance(self, integrator: Ncm.SBesselIntegratorFFTL) -> None:
        """Time comparison between integrate_ell loop and integrate method."""
        center = 1000.0
        std = 10.0

        def test_func(x: float, _k: float) -> float:
            return x**2 / (1 + ((x - center) / std) ** 2) ** 2

        a = center - 20.0 * std
        b = center + 20.0 * std

        integrator.set_ell_range(0, 200)

        ell_min, ell_max = integrator.get_ell_range()
        n_ell = ell_max - ell_min + 1

        # Time integrate_ell approach
        start_time = time.time()
        results_ell = []
        for ell in range(ell_min, ell_max + 1):
            result = integrator.integrate_ell(test_func, a, b, 1.0, ell)
            results_ell.append(result)
        time_integrate_ell = time.time() - start_time

        # Time integrate approach (run 10 times and average)
        times_integrate = []
        for _ in range(500):
            start_time = time.time()
            results_vec = Ncm.Vector.new(n_ell)
            integrator.integrate(test_func, a, b, 1.0, results_vec)
            results_all = np.array([results_vec.get(i) for i in range(n_ell)])
            times_integrate.append(time.time() - start_time)
        time_integrate = np.mean(times_integrate)

        # Print timing results
        print(f"\nPerformance comparison for {n_ell} multipoles:")
        print(f"  integrate_ell (sequential): {time_integrate_ell:.4f} seconds")
        print(f"  integrate (all at once):    {time_integrate:.4f} seconds")
        print(f"  Speedup: {time_integrate_ell / time_integrate:.2f}x")

        # Verify results match (correctness is what matters for tests)
        assert_allclose(results_all, results_ell, rtol=1e-7, atol=0.0)

    @pytest.mark.parametrize(
        "func_type,filename",
        [
            ("gaussian", "gauss_jl_500.json.gz"),
            ("rational", "rational_jl_500.json.gz"),
        ],
    )
    def test_truth_table(self, func_type: str, filename: str) -> None:
        """Test against truth tables for spherical Bessel integrals.

        This test verifies that for each multipole ell, the integrator achieves
        accurate results (rel_error < reltol) for k values up to at least min_k_ratio *
        ell. FFTL benefits from computing all ells at once in a single call.
        """
        match func_type:
            case "gaussian":
                reltol = 1.0e-13
            case "rational":
                reltol = 4.0e-2
        min_k_ratio = 1.2

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

        # Create integrator for all ells at once (FFTL's strength)
        integrator = Ncm.SBesselIntegratorFFTL.new(ell_min, ell_max)

        # Get the appropriate integration method
        if func_type == "gaussian":
            integrate_func = integrator.integrate_gaussian
        elif func_type == "rational":
            integrate_func = integrator.integrate_rational
        else:
            raise ValueError(f"Unknown function type: {func_type}")

        results_vec = Ncm.Vector.new(n_ells)

        # Compute all relative errors: shape (n_ells, n_k)
        rel_errors = np.zeros((n_ells, n_k))

        # Compute all k values (FFTL efficiently computes all ells at once)
        for i, k in enumerate(kvals):
            integrate_func(center, std, lb, ub, k, results_vec)
            results = results_vec.to_numpy()

            truth_values = table[:, i]
            rel_errors[:, i] = np.abs(
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

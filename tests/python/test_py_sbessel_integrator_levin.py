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

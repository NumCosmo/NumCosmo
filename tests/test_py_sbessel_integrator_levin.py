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

        Compares the result of the integrate_ell method against scipy numerical integration.
        """
        N = 128
        a, b = 1.0, 20.0

        # Create integrator
        integrator = Ncm.SBesselIntegratorLevin.new(l_val, l_val)
        integrator.set_max_order(N)

        # Define f(x) = 1
        def f_constant(_user_data: None, _x: float) -> float:
            return 1.0

        # Compute integral using the integrate_ell method
        result = integrator.integrate_ell(f_constant, a, b, l_val, None)

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

        Compares the result of the integrate_ell method against scipy numerical integration.
        """
        N = 128
        a, b = 1.0, 20.0

        # Create integrator
        integrator = Ncm.SBesselIntegratorLevin.new(l_val, l_val)
        integrator.set_max_order(N)

        # Define f(x) = x
        def f_linear(_user_data: None, x: float) -> float:
            return x

        # Compute integral using the integrate_ell method
        result = integrator.integrate_ell(f_linear, a, b, l_val, None)

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

        Compares the result of the integrate_ell method against scipy numerical integration.
        """
        N = 128
        a, b = 1.0, 20.0

        # Create integrator
        integrator = Ncm.SBesselIntegratorLevin.new(l_val, l_val)
        integrator.set_max_order(N)

        # Define f(x) = x^2
        def f_quadratic(_user_data: None, x: float) -> float:
            return x * x

        # Compute integral using the integrate_ell method
        result = integrator.integrate_ell(f_quadratic, a, b, l_val, None)

        # Compute expected value using scipy
        def integrand(x: float) -> float:
            return x * x * spherical_jn(l_val, x)

        expected, _ = quad(integrand, a, b, epsabs=1e-12, epsrel=1e-14)

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

        Compares the result of the integrate_ell method against scipy numerical integration.
        """
        N = 256
        a, b = 1.0, 20.0
        center = 10.5
        std = 3.0

        # Create integrator
        integrator = Ncm.SBesselIntegratorLevin.new(l_val, l_val)
        integrator.set_max_order(N)

        # Define rational function f(x) = x^2 / (1 + ((x-center)/std)^2)^3
        def f_rational(_user_data: None, x: float) -> float:
            dx = (x - center) / std
            return x * x / ((1.0 + dx * dx) ** 3)

        # Compute integral using the integrate_ell method
        result = integrator.integrate_ell(f_rational, a, b, l_val, None)

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
        integrator.prepare()

        # Define a test function
        def f_test(_user_data: None, x: float) -> float:
            return np.exp(-0.1 * x)

        # Use integrate (batched internally)
        results_vec = Ncm.Vector.new(lmax - lmin + 1)
        integrator.integrate(f_test, a, b, results_vec, None)
        results_range_np = results_vec.to_numpy()

        # Compute individually for each l
        results_individual = []
        for ell in range(lmin, lmax + 1):
            results_individual.append(integrator.integrate_ell(f_test, a, b, ell, None))

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
        """Test against truth tables for spherical Bessel integrals."""

        truth_table_path = Path(
            Ncm.cfg_get_data_filename(f"truth_tables/{filename}", True)
        )
        with gzip.open(truth_table_path, "rt") as f:
            truth_table = json.load(f)

        print(f"Loaded {func_type} truth table with {len(truth_table)} entries.")

        center = truth_table["center"]
        std = truth_table["std"]
        lb = truth_table["lower-bound"]
        ub = truth_table["upper-bound"]
        table = np.array(truth_table["table"])

        print(f"Preparing solver for {func_type} truth table...")
        print(f"  center = {center}")
        print(f"  std    = {std}")
        print(f"  lb     = {lb}")
        print(f"  ub     = {ub}")

        ells = truth_table["lvals"]
        ell_min = int(np.min(ells))
        ell_max = int(np.max(ells))
        integrator = Ncm.SBesselIntegratorLevin.new(ell_min, ell_max)
        integrator.prepare()

        # Get the appropriate integration method
        if func_type == "gaussian":
            integrate_func = integrator.integrate_gaussian
        elif func_type == "rational":
            integrate_func = integrator.integrate_rational
        else:
            raise ValueError(f"Unknown function type: {func_type}")

        N = 2**16  # Number of Chebyshev nodes
        print(f"Using N = {N} Chebyshev nodes")

        results_vec = Ncm.Vector.new(ell_max - ell_min + 1)
        print_rank = False
        print_ell: list[int] | None = [50]

        for i in range(1):
            print(f"Starting iteration {i}\r", end="", flush=True)
            for i, k in enumerate(truth_table["kvals"]):
                a = lb * k
                b = ub * k
                integrate_func(center, std, k, a, b, results_vec)
                results = results_vec.to_numpy()

                truth_values = table[: (ell_max - ell_min + 1), i]

                # Compute relative errors
                rel_errors = np.abs(
                    (results - truth_values) / np.maximum(np.abs(truth_values), 1.0e-50)
                )

                if print_ell is not None:
                    for ell in print_ell:
                        print(
                            f"[{func_type}] ell={ell:d}, k={k: 22.15g}, "
                            f"result={results[ell - ell_min]: 14.6e}, "
                            f"truth={truth_values[ell - ell_min]: 14.6e}, "
                            f"rel_error={rel_errors[ell - ell_min]: 4.2e}"
                        )

                if print_rank:
                    # Find best and worst agreement
                    best_idx = np.argmin(rel_errors)
                    worst_idx = np.argmax(rel_errors)
                    print(
                        f"\n[{func_type}] Testing k={k} with "
                        f"{ell_min}--{ell_max} multipoles:"
                    )
                    print(
                        f"  Best agreement:  ell={ell_min + best_idx}, "
                        f"rel_error={rel_errors[best_idx]:.2e}"
                    )
                    print(
                        f"  Worst agreement: ell={ell_min + worst_idx}, "
                        f"rel_error={rel_errors[worst_idx]:.2e}"
                    )

                    # Find 10 worst agreements
                    worst_10_indices = np.argsort(rel_errors)[-10:][::-1]
                    print("  Top 10 worst agreements:")
                    for rank, idx in enumerate(worst_10_indices, 1):
                        ell = ell_min + idx
                        print(
                            f"    {rank}. ell={ell}: result={results[idx]:.6e}, "
                            f"truth={truth_values[idx]:.6e}, "
                            f"rel_error={rel_errors[idx]:.2e}"
                        )
                    best_10_indices = np.argsort(rel_errors)[:10]
                    print("  Top 10 best agreements:")
                    for rank, idx in enumerate(best_10_indices, 1):
                        ell = ell_min + idx
                        print(
                            f"    {rank}. ell={ell}: result={results[idx]:.6e}, "
                            f"truth={truth_values[idx]:.6e}, "
                            f"rel_error={rel_errors[idx]:.2e}"
                        )

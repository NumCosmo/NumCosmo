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

from numcosmo_py import Ncm


class TestSBesselIntegratorLevin:
    """Tests for NcmSBesselIntegratorLevin."""

    @pytest.fixture
    def integrator(self) -> Ncm.SBesselIntegratorLevin:
        """Create a Levin integrator."""
        return Ncm.SBesselIntegratorLevin.new(0, 10)

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
        # ell_max = 200  # Limit for faster testing

        N = 2**16  # Number of Chebyshev nodes
        print(f"Using N = {N} Chebyshev nodes")

        print_rank = False
        print_ell: list[int] | None = [50]
        solver = Ncm.SBesselOdeSolver.new(0, lb, ub)

        for i in range(1):
            print(f"Starting iteration {i}\r", end="", flush=True)
            for i, k in enumerate(truth_table["kvals"]):
                # Get the appropriate integration method
                if func_type == "gaussian":
                    results_vec = solver.integrate_gaussian_l_range(
                        center, std, k, N, ell_min, ell_max
                    )
                elif func_type == "rational":
                    results_vec = solver.integrate_rational_l_range(
                        center, std, k, N, ell_min, ell_max
                    )
                else:
                    raise ValueError(f"Unknown function type: {func_type}")

                results = np.array(results_vec.dup_array())
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

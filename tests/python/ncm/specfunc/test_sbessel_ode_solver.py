#
# test_py_sbessel_ode_solver.py
#
# Fri Jan 31 00:00:00 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_py_sbessel_ode_solver.py
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

"""Tests for NcmSBesselOdeSolver.

The solver computes solutions to the modified spherical Bessel ODE:
    u''(x) + [x^2 - l(l+1)]u(x) = RHS

where u(x) = x*j_l(x) or u(x) = x*y_l(x) for the homogeneous case.
This formulation eliminates the first derivative term compared to the
standard spherical Bessel equation.
"""

import numpy as np
import pytest

from numpy.testing import assert_allclose
from scipy.special import spherical_jn, spherical_yn
from scipy.integrate import quad
from scipy.optimize import brentq

from numcosmo_py import Ncm


class TestSBesselStoredRotations:
    """Tests for stored rotation reuse: targeting the max_c_A/quiet_cols bug.

    The bug occurred when:
    1. Stored rotations were applied (updating max_c_A and quiet_cols)
    2. Extension was needed (more columns required)
    3. The continuation loop reset max_c_A=0.0 and quiet_cols=0
    4. This caused premature exit or incorrect convergence checks

    These tests ensure the fix properly carries state across stored/new rotations.
    """

    def test_stored_rotations_reuse_identical(self) -> None:
        """Test that reusing stored rotations gives identical results.

        Solves the same RHS twice - second solve should reuse stored rotations
        and produce identical solution.
        """
        solver = Ncm.SBesselOdeSolver.new()
        op = solver.create_operator(1.0, 10.0, 5, 5)

        # Create a non-trivial RHS
        rhs = np.zeros(50)
        rhs[10] = 1.0
        rhs[15] = 0.5
        rhs[20] = 0.25

        # First solve - builds rotations
        sol1, len1 = op.solve(rhs)
        sol1_np = np.array(sol1)
        n_cols_1 = op.get_n_cols()

        # Second solve - reuses rotations
        sol2, len2 = op.solve(rhs)
        sol2_np = np.array(sol2)
        n_cols_2 = op.get_n_cols()

        # Results should be identical
        assert len1 == len2, "Solution lengths should match"
        assert n_cols_1 == n_cols_2, "Number of columns should match"
        assert_allclose(
            sol1_np,
            sol2_np,
            rtol=1e-15,
            atol=1e-15,
            err_msg="Stored rotations should produce identical results",
        )

    def test_stored_rotations_extend_rhs(self) -> None:
        """Test extending stored rotations when RHS is longer.

        This specifically tests the bug scenario:
        1. Solve with short RHS (stores rotations)
        2. Solve with longer RHS (needs extension)
        3. Verify max_c_A and quiet_cols carry over correctly

        The bug would cause premature exit when max_c_A was reset to 0.0.
        """
        solver = Ncm.SBesselOdeSolver.new()
        op = solver.create_operator(1.0, 10.0, 10, 10)

        # First solve with short RHS
        rhs_short = np.zeros(40)
        rhs_short[10] = 1.0
        _, len_short = op.solve(rhs_short)
        n_cols_short = op.get_n_cols()

        print(f"Short RHS: n_cols={n_cols_short}, sol_len={len_short}")

        # Second solve with LONGER RHS - requires extension
        # This triggers the bug: stored rotations applied, then needs more columns
        rhs_long = np.zeros(80)
        rhs_long[10] = 1.0  # Same initial coefficients
        rhs_long[15] = 0.5  # Additional coefficients beyond stored
        rhs_long[50] = 0.1  # Even further out
        rhs_long[60] = 0.05

        # Solve from scratch (no stored rotations) - reference solution
        op_fresh = solver.create_operator(1.0, 10.0, 10, 10)
        sol_fresh, len_fresh = op_fresh.solve(rhs_long)
        sol_fresh_np = np.array(sol_fresh)
        n_cols_fresh = op_fresh.get_n_cols()

        # Solve with stored rotations + extension
        sol_extended, len_extended = op.solve(rhs_long)
        sol_extended_np = np.array(sol_extended)
        n_cols_extended = op.get_n_cols()

        # Should get same result as fresh solve
        assert (
            n_cols_extended == n_cols_fresh
        ), "Extended solve should converge at same column as fresh solve"
        assert len_extended == len_fresh, "Solution lengths should match"

        max_coeff = np.max(np.abs(sol_fresh_np))
        assert_allclose(
            sol_extended_np,
            sol_fresh_np,
            rtol=1e-12,
            atol=1e-12 * max_coeff,
            err_msg="Extended stored rotations should match fresh solve",
        )

        # Critical: verify solution is non-zero (bug would give premature zero exit)
        assert (
            np.max(np.abs(sol_extended_np)) > 1e-10
        ), "Solution should be non-zero (bug check: max_c_A=0 premature exit)"

    def test_stored_rotations_multiple_different_rhs(self) -> None:
        """Test solving multiple different RHS vectors with stored rotations.

        Each solve should reuse rotations correctly and give accurate results.
        """
        solver = Ncm.SBesselOdeSolver.new()
        op = solver.create_operator(1.0, 10.0, 8, 8)

        rhs_list = [
            (np.zeros(50), 10, 1.0),  # Single spike at index 10
            (np.zeros(50), 20, 0.8),  # Single spike at index 20
            (np.zeros(50), 15, 1.5),  # Single spike at index 15
            (np.zeros(60), 25, 1.2),  # Longer RHS, spike at 25
            (np.zeros(45), 12, 0.9),  # Shorter RHS, spike at 12
        ]

        solutions = []

        for idx, (rhs, spike_pos, spike_val) in enumerate(rhs_list):
            rhs[spike_pos] = spike_val

            # Solve with stored rotations
            sol, _ = op.solve(rhs)
            sol_np = np.array(sol)
            solutions.append(sol_np)

            # Verify against fresh solve
            op_fresh = solver.create_operator(1.0, 10.0, 8, 8)
            sol_fresh, _ = op_fresh.solve(rhs)
            sol_fresh_np = np.array(sol_fresh)

            max_coeff = np.max(np.abs(sol_fresh_np))
            assert_allclose(
                sol_np,
                sol_fresh_np,
                rtol=1e-12,
                atol=1e-12 * max_coeff,
                err_msg=f"Solve {idx+1}: stored rotations vs fresh",
            )

        # All solutions should be different (different RHS)
        for i in range(len(solutions) - 1):
            min_len = min(len(solutions[i]), len(solutions[i + 1]))
            diff = np.max(np.abs(solutions[i][:min_len] - solutions[i + 1][:min_len]))
            assert diff > 1e-8, f"Solutions {i} and {i+1} should differ"

    def test_stored_rotations_batched_extend(self) -> None:
        """Test batched stored rotations with RHS extension.

        This tests the batched version with multiple ell values,
        ensuring max_c_A and quiet_cols carry over correctly for each ell.
        """
        solver = Ncm.SBesselOdeSolver.new()

        # Test with multiple ell values
        lmin, lmax = 5, 8
        n_ell = lmax - lmin + 1
        op = solver.create_operator(1.0, 10.0, lmin, lmax)

        # First solve with short RHS
        rhs_short = np.zeros(40)
        rhs_short[12] = 1.0
        rhs_short[15] = 0.3

        sols_short, len_short = op.solve(rhs_short)
        _ = np.array(sols_short).reshape(n_ell, len_short)
        _ = op.get_n_cols()

        # Second solve with LONGER RHS - tests extension
        rhs_long = np.zeros(70)
        rhs_long[12] = 1.0
        rhs_long[15] = 0.3
        rhs_long[35] = 0.15  # Beyond stored rotations
        rhs_long[50] = 0.08

        # Fresh solve for reference
        op_fresh = solver.create_operator(1.0, 10.0, lmin, lmax)
        sols_fresh, len_fresh = op_fresh.solve(rhs_long)
        sols_fresh_np = np.array(sols_fresh).reshape(n_ell, len_fresh)
        n_cols_fresh = op_fresh.get_n_cols()

        # Extended solve (reuses stored rotations)
        sols_extended, len_extended = op.solve(rhs_long)
        sols_extended_np = np.array(sols_extended).reshape(n_ell, len_extended)
        n_cols_extended = op.get_n_cols()

        # Should converge at same column
        assert (
            n_cols_extended == n_cols_fresh
        ), f"Extended should match fresh: {n_cols_extended} vs {n_cols_fresh}"

        # Solutions should match for all ell values
        for i, ell in enumerate(range(lmin, lmax + 1)):
            max_coeff = np.max(np.abs(sols_fresh_np[i, :]))
            assert_allclose(
                sols_extended_np[i, :],
                sols_fresh_np[i, :],
                rtol=1e-12,
                atol=1e-12 * max_coeff,
                err_msg=f"ell={ell}: extended vs fresh with stored rotations",
            )

            # Critical: solutions should be non-zero
            assert (
                max_coeff > 1e-10
            ), f"ell={ell}: solution should be non-zero (bug check)"

    def test_stored_rotations_zero_rhs_edge_case(self) -> None:
        """Test that zero RHS portions don't cause max_c_A reset issues.

        When RHS has leading zeros followed by non-zero values,
        the algorithm should handle the zero detection correctly.
        """
        solver = Ncm.SBesselOdeSolver.new()
        op = solver.create_operator(1.0, 10.0, 5, 5)

        # First solve with non-zero RHS
        rhs1 = np.zeros(50)
        rhs1[15] = 1.0
        _, _ = op.solve(rhs1)
        _ = op.get_n_cols()

        # Second solve with more leading zeros, non-zero later
        # This tests that max_c_A from stored rotations carries over
        rhs2 = np.zeros(60)
        rhs2[40] = 0.5  # Non-zero far beyond stored rotations

        sol2, _ = op.solve(rhs2)
        sol2_np = np.array(sol2)
        _ = op.get_n_cols()

        # Solution should be non-zero despite leading zeros
        assert np.max(np.abs(sol2_np)) > 1e-10, (
            "Solution should be non-zero - bug would exit at 'if (max_c_A == 0.0)' "
            "if state wasn't carried over"
        )

        # Verify against fresh solve
        op_fresh = solver.create_operator(1.0, 10.0, 5, 5)
        sol_fresh, _ = op_fresh.solve(rhs2)
        sol_fresh_np = np.array(sol_fresh)

        max_coeff = np.max(np.abs(sol_fresh_np))
        assert_allclose(
            sol2_np,
            sol_fresh_np,
            rtol=1e-12,
            atol=1e-12 * max_coeff,
            err_msg="Extended solve with leading zeros should match fresh",
        )

    def test_stored_rotations_convergence_tracking(self) -> None:
        """Test that quiet_cols counter properly carries over.

        The quiet_cols counter tracks consecutive columns with small coefficients.
        It must carry over from stored rotations to new rotations for correct
        convergence detection.
        """
        solver = Ncm.SBesselOdeSolver.new()
        op = solver.create_operator(1.0, 15.0, 12, 12)

        # First solve - creates stored rotations
        rhs1 = np.zeros(60)
        rhs1[18] = 1.0
        rhs1[22] = 0.1

        _, _ = op.solve(rhs1)
        _ = op.get_n_cols()

        # Second solve - similar RHS but slightly longer
        # Should converge at similar column count
        rhs2 = np.zeros(70)
        rhs2[18] = 1.0
        rhs2[22] = 0.1
        rhs2[35] = 0.02  # Small coefficient that should quickly quiet

        sol2, _ = op.solve(rhs2)
        n_cols_2 = op.get_n_cols()

        # Verify against fresh solve - should converge at same point
        op_fresh = solver.create_operator(1.0, 15.0, 12, 12)
        sol_fresh, _ = op_fresh.solve(rhs2)
        n_cols_fresh = op_fresh.get_n_cols()

        # Convergence point should match (allow small tolerance due to numerical
        # accumulation)
        assert abs(n_cols_2 - n_cols_fresh) <= 2, (
            f"Extended solve should converge at approximately same column as fresh: "
            f"{n_cols_2} vs {n_cols_fresh} (diff={n_cols_2 - n_cols_fresh}, "
            f"quiet_cols tracking may have issues)"
        )

        # Solutions should match
        sol2_np = np.array(sol2)
        sol_fresh_np = np.array(sol_fresh)
        max_coeff = np.max(np.abs(sol_fresh_np))

        assert_allclose(
            sol2_np,
            sol_fresh_np,
            rtol=1e-12,
            atol=1e-12 * max_coeff,
            err_msg="Extended solve should match fresh (convergence tracking)",
        )

    def test_stored_rotations_batched_multiple_solves(self) -> None:
        """Test batched mode with multiple sequential solves and different RHS.

        Ensures max_c_A and quiet_cols arrays (one per ell) carry over correctly
        in batched mode.
        """
        solver = Ncm.SBesselOdeSolver.new()

        lmin, lmax = 3, 6
        n_ell = lmax - lmin + 1
        op = solver.create_operator(1.0, 12.0, lmin, lmax)

        # Solve with progressively longer RHS
        rhs_lengths = [40, 50, 65, 75]
        previous_n_cols = 0

        for rhs_len in rhs_lengths:
            rhs = np.zeros(rhs_len)
            rhs[15] = 1.0
            rhs[20] = 0.4
            if rhs_len > 40:
                rhs[35] = 0.15
            if rhs_len > 60:
                rhs[55] = 0.05

            # Solve with stored rotations
            sols, sol_len = op.solve(rhs)
            n_cols = op.get_n_cols()

            # Fresh solve for comparison
            op_fresh = solver.create_operator(1.0, 12.0, lmin, lmax)
            sols_fresh, sol_len_fresh = op_fresh.solve(rhs)
            sols_fresh_np = np.array(sols_fresh).reshape(n_ell, sol_len_fresh)

            # Verify each ell
            sols_np = np.array(sols).reshape(n_ell, sol_len)

            for i, ell in enumerate(range(lmin, lmax + 1)):
                max_coeff = np.max(np.abs(sols_fresh_np[i, :]))
                assert_allclose(
                    sols_np[i, :],
                    sols_fresh_np[i, :],
                    rtol=1e-12,
                    atol=1e-12 * max_coeff,
                    err_msg=f"RHS len={rhs_len}, ell={ell}: batched stored vs fresh",
                )

            # Verify n_cols increases or stays same (never decreases)
            assert (
                n_cols >= previous_n_cols
            ), f"n_cols should not decrease: {previous_n_cols} -> {n_cols}"
            previous_n_cols = n_cols

    def test_stored_rotations_high_ell_extension(self) -> None:
        """Test stored rotation extension with high ell values.

        High ell values require more columns to converge. This tests that
        extension works correctly when many columns beyond the stored ones
        are needed.
        """
        solver = Ncm.SBesselOdeSolver.new()
        op = solver.create_operator(1.0, 20.0, 50, 50)

        # First solve - short RHS
        rhs_short = np.zeros(80)
        rhs_short[25] = 1.0
        rhs_short[30] = 0.2

        _sol_short, _len_short = op.solve(rhs_short)
        n_cols_short = op.get_n_cols()

        # Second solve - much longer RHS requiring many more columns
        rhs_long = np.zeros(150)
        rhs_long[25] = 1.0
        rhs_long[30] = 0.2
        rhs_long[60] = 0.1
        rhs_long[90] = 0.05

        # Fresh reference
        op_fresh = solver.create_operator(1.0, 20.0, 50, 50)
        sol_fresh, _ = op_fresh.solve(rhs_long)
        sol_fresh_np = np.array(sol_fresh)
        n_cols_fresh = op_fresh.get_n_cols()

        # Extended solve
        sol_extended, _ = op.solve(rhs_long)
        sol_extended_np = np.array(sol_extended)
        n_cols_extended = op.get_n_cols()

        # Should match
        assert n_cols_extended == n_cols_fresh
        max_coeff = np.max(np.abs(sol_fresh_np))
        assert_allclose(
            sol_extended_np,
            sol_fresh_np,
            rtol=1e-11,
            atol=1e-11 * max_coeff,
            err_msg="High ell extended should match fresh",
        )

        # Verify many columns were added
        columns_added = n_cols_extended - n_cols_short
        assert columns_added >= 20, (
            f"Should have added many columns for high ell: " f"added {columns_added}"
        )


class TestSBesselTauConstraint:
    """The tau constraint against Dirichlet data.

    The boundary functional of the Levin reduction is invariant under adding
    homogeneous solutions to u, so both constraints must return the same panel
    integral. They differ in which member of the solution family is represented:
    Dirichlet data forces the oscillatory one, the tau constraint leaves the smooth
    one, which needs only the forcing's own order. The tau constraint is valid only
    while the homogeneous solutions are unrepresentable at that order; on a short or
    evanescent panel they are not, and the solve admits them. Both regimes are
    covered here so the guard in the integrator has a documented target.
    """

    ELL = 20

    @staticmethod
    def _forcing_rhs(a: float, b: float, order: int = 96) -> np.ndarray:
        """Endpoint data followed by the C^(2) coefficients of x F(x).

        F is a Gaussian bump centred on the panel with width a fixed fraction of
        it, so its Chebyshev order is the same on every panel.
        """
        m, h = 0.5 * (a + b), 0.5 * (b - a)
        spectral = Ncm.Spectral.new()

        def forcing(_user_data, t):
            x = m + h * t
            return x * np.exp(-0.5 * ((x - m) / (0.35 * h)) ** 2)

        cheb = np.array(
            spectral.compute_chebyshev_coeffs(forcing, -1.0, 1.0, order, None)
        )
        gegen = np.array(Ncm.Spectral.chebT_to_gegenbauer_alpha2(cheb))
        rhs = np.zeros(len(gegen) + 2)
        rhs[2:] = gegen

        return rhs

    @staticmethod
    def _boundary_functional(coeffs: np.ndarray, a: float, b: float, ell: int) -> float:
        """W(b) - W(a) with W = x j_l u' - (x j_l)' u, u from its Chebyshev series."""
        h = 0.5 * (b - a)
        u_a, u_b = np.polynomial.chebyshev.chebval([-1.0, 1.0], coeffs)
        du = (
            np.polynomial.chebyshev.chebval(
                [-1.0, 1.0], np.polynomial.chebyshev.chebder(coeffs)
            )
            / h
        )

        def term(x: float, u: float, dudy: float) -> float:
            jl = spherical_jn(ell, x)
            djl = spherical_jn(ell, x, derivative=True)
            return x * jl * dudy - (jl + x * djl) * u

        return term(b, u_b, du[1]) - term(a, u_a, du[0])

    def _solve(self, a: float, b: float, ell_min: int, ell_max: int, tau: bool):
        solver = Ncm.SBesselOdeSolver.new()
        solver.set_tolerance(1.0e-12)
        solver.set_default_constraint(
            Ncm.SBesselOdeConstraint.TAU if tau else Ncm.SBesselOdeConstraint.DIRICHLET
        )
        op = solver.create_operator(a, b, ell_min, ell_max)
        coeffs, n_cols = op.solve(self._forcing_rhs(a, b))

        return op, np.array(coeffs), n_cols

    def test_same_integral_far_fewer_columns_on_deep_panel(self) -> None:
        """Deep oscillatory panel: identical integral, ~400x fewer columns."""
        a, b = 1.0e4, 10.0**4.5
        op_d, c_d, n_d = self._solve(a, b, self.ELL, self.ELL, False)
        op_f, c_f, n_f = self._solve(a, b, self.ELL, self.ELL, True)

        i_d = self._boundary_functional(c_d, a, b, self.ELL)
        i_f = self._boundary_functional(c_f, a, b, self.ELL)

        assert_allclose(i_f, i_d, rtol=1.0e-9)
        # The tau solve stops at its forcing floor, 1.5 times the forcing order, which
        # is set by F alone; Dirichlet runs to the panel's oscillation count.
        n_forcing = len(self._forcing_rhs(a, b)) - 2
        assert (
            n_f <= 2 * n_forcing
        ), f"tau constraint used {n_f} columns for a forcing of order {n_forcing}"
        assert n_f * 20 < n_d, f"tau constraint used {n_f} columns, Dirichlet {n_d}"
        assert op_f.get_min_cols() == 0
        assert op_d.get_min_cols() > 0
        assert op_d.get_constraint() == Ncm.SBesselOdeConstraint.DIRICHLET
        assert op_f.get_constraint() == Ncm.SBesselOdeConstraint.TAU

    def test_tau_constraint_matches_quadrature_on_narrow_panel(self) -> None:
        """Narrow panel, where direct quadrature is affordable."""
        a, b = 100.0, 10.0**2.5
        m, h = 0.5 * (a + b), 0.5 * (b - a)
        _, c_f, _ = self._solve(a, b, self.ELL, self.ELL, True)

        ref, _ = quad(
            lambda x: np.exp(-0.5 * ((x - m) / (0.35 * h)) ** 2)
            * spherical_jn(self.ELL, x),
            a,
            b,
            limit=5000,
            epsabs=0.0,
            # A tighter request makes QUADPACK bisect into the oscillation and report
            # roundoff without moving the value: at 1e-10 it agrees with a composite
            # Gauss-Legendre reference to 1e-13, at 1e-13 only to 8e-12.
            epsrel=1.0e-10,
        )

        assert_allclose(
            self._boundary_functional(c_f, a, b, self.ELL), ref, rtol=1.0e-9
        )

    def test_batched_tau_constraint_matches_single(self) -> None:
        """A block shares one factorization; every member must match its own solve."""
        a, b = 1.0e4, 10.0**4.5
        lmin, lmax = self.ELL, self.ELL + 7
        n_ell = lmax - lmin + 1
        _, c_b, n_b = self._solve(a, b, lmin, lmax, True)
        c_b = c_b.reshape(n_ell, n_b)

        for i, ell in enumerate(range(lmin, lmax + 1)):
            _, c_s, n_s = self._solve(a, b, ell, ell, True)
            padded = np.zeros(max(n_b, n_s))
            padded[:n_s] = c_s
            assert_allclose(
                c_b[i],
                padded[:n_b],
                rtol=1.0e-9,
                atol=1.0e-13 * np.abs(c_s).max(),
                err_msg=f"batched tau constraint differs from single solve at ell={ell}",
            )

    def test_shallow_panel_admits_homogeneous_content(self) -> None:
        """Where N_min is small the tau solve is contaminated; where large it is not.

        This is the failure the integrator's guard detects: on [3.162, 10] at l=2 the
        panel holds ~4 oscillations against a working order of ~32, so the
        homogeneous solutions are representable and the coefficient norm explodes.
        """
        a, b, ell = 3.162, 10.0, 2
        _, c_d, _ = self._solve(a, b, ell, ell, False)
        _, c_f, _ = self._solve(a, b, ell, ell, True)
        ratio_shallow = np.abs(c_f).max() / np.abs(c_d).max()

        a, b = 100.0, 10.0**2.5
        _, c_d, _ = self._solve(a, b, ell, ell, False)
        _, c_f, _ = self._solve(a, b, ell, ell, True)
        ratio_deep = np.abs(c_f).max() / np.abs(c_d).max()

        assert (
            ratio_shallow > 1.0e6
        ), f"expected contamination, got ratio {ratio_shallow:.3e}"
        assert (
            ratio_deep < 10.0
        ), f"unexpected contamination on a deep panel: {ratio_deep:.3e}"

    @pytest.mark.parametrize("n_ell", [1, 8])
    def test_last_max_coeff_tracks_the_solution(self, n_ell: int) -> None:
        """The per-ell max |a_j| the guard reads equals the solution's own maximum."""
        a, b = 1.0e4, 10.0**4.5
        op, coeffs, n_cols = self._solve(a, b, self.ELL, self.ELL + n_ell - 1, True)
        coeffs = coeffs.reshape(n_ell, n_cols)

        for i in range(n_ell):
            assert op.get_last_max_coeff(i) == np.abs(coeffs[i]).max()

    def test_last_max_coeff_after_reconfigure_to_wider_block(self) -> None:
        """Regression: storage sized for one multipole, then reconfigured to eight.

        The tracking array is allocated with the accumulator arrays; a scalar solve
        allocates for one multipole, and a later reconfigure to a block must not
        write past it.
        """
        a, b = 1.0e4, 10.0**4.5
        solver = Ncm.SBesselOdeSolver.new()
        solver.set_tolerance(1.0e-12)
        solver.set_default_constraint(Ncm.SBesselOdeConstraint.TAU)
        rhs = self._forcing_rhs(a, b)

        op = solver.create_operator(a, b, self.ELL, self.ELL)
        op.solve(rhs)

        solver.reconfigure_operator(op, a, b, self.ELL, self.ELL + 7)
        coeffs, n_cols = op.solve(rhs)
        coeffs = np.array(coeffs).reshape(8, n_cols)

        for i in range(8):
            assert op.get_last_max_coeff(i) == np.abs(coeffs[i]).max()

    def test_pinned_constraint_same_integral(self) -> None:
        """Pinning two coefficients selects another member with the same integral."""
        a, b = 1.0e4, 10.0**4.5
        solver = Ncm.SBesselOdeSolver.new()
        solver.set_tolerance(1.0e-12)
        rhs = self._forcing_rhs(a, b)

        op = solver.create_operator(a, b, self.ELL, self.ELL)
        c_d, _ = op.solve(rhs)
        i_d = self._boundary_functional(np.array(c_d), a, b, self.ELL)
        floor = op.get_min_cols()

        op.set_pinned_constraint(120, 144)
        assert op.get_constraint() == Ncm.SBesselOdeConstraint.PINNED
        assert op.get_pins() == (120, 144)
        assert op.get_min_cols() == 145
        c_p, n_p = op.solve(rhs)
        assert n_p == 145
        assert_allclose(
            self._boundary_functional(np.array(c_p), a, b, self.ELL), i_d, rtol=1.0e-9
        )

        op.set_constraint(Ncm.SBesselOdeConstraint.DIRICHLET)
        assert op.get_constraint() == Ncm.SBesselOdeConstraint.DIRICHLET
        assert op.get_min_cols() == floor

    def test_constraint_setting_propagates_and_overrides(self) -> None:
        """Operators inherit the solver's constraint; a per-operator setter overrides it."""
        a, b = 1.0e4, 10.0**4.5
        solver = Ncm.SBesselOdeSolver.new()
        assert solver.get_default_constraint() == Ncm.SBesselOdeConstraint.DIRICHLET

        solver.set_default_constraint(Ncm.SBesselOdeConstraint.TAU)
        op = solver.create_operator(a, b, self.ELL, self.ELL)
        assert op.get_constraint() == Ncm.SBesselOdeConstraint.TAU
        assert op.get_min_cols() == 0

        op.set_constraint(Ncm.SBesselOdeConstraint.DIRICHLET)
        assert op.get_constraint() == Ncm.SBesselOdeConstraint.DIRICHLET
        assert op.get_min_cols() > 0

        op.set_min_cols(7)
        assert op.get_min_cols() == 7
        op.set_min_cols(-1)
        assert op.get_min_cols() > 7


class TestTauFloorFactorAndDiagnostics:
    """The tau working order, the roundoff bound it implies, and the dense matrix.

    The tau constraint has no oscillatory floor, and the decay test does not stop later than
    the floor in practice, so the floor factor sets the working order outright. It has to
    exceed one: at exactly the forcing's order the forcing itself is unresolved.
    """

    ELL = 20
    A, B = 1.0e4, 10.0**4.5

    def _rhs(self):
        """C^(2) coefficients of x F(x) for a bump on [A, B], boundary rows included."""
        solver = Ncm.SBesselOdeSolver.new()
        spectral = solver.peek_spectral()
        mid, half = 0.5 * (self.A + self.B), 0.5 * (self.B - self.A)
        cheb = np.array(
            spectral.compute_chebyshev_coeffs_adaptive(
                lambda _d, x: x * np.exp(-0.5 * ((x - mid) / (0.25 * half)) ** 2),
                self.A,
                self.B,
                3,
                1.0e-13,
                None,
            )[1]
        )
        n = len(cheb)
        padded = np.concatenate([cheb, np.zeros(8)])
        kk = np.arange(n)
        geg = (
            padded[:n] / (2 * (kk + 1))
            - (kk + 2) * padded[2 : n + 2] / ((kk + 1) * (kk + 3))
            + padded[4 : n + 4] / (2 * (kk + 3))
        )
        geg[0] += 0.5 * padded[0]

        return np.concatenate([[0.0, 0.0], geg]), n

    def test_default_factor(self) -> None:
        """The factor ships above one and reads back."""
        solver = Ncm.SBesselOdeSolver.new()
        assert solver.get_tau_floor_factor() > 1.0

        solver.set_tau_floor_factor(1.4)
        assert solver.get_tau_floor_factor() == 1.4

    def test_factor_sets_the_working_order(self) -> None:
        """The columns kept track the factor times the forcing order."""
        rhs, n_forcing = self._rhs()

        for factor in (1.1, 1.5, 2.0):
            solver = Ncm.SBesselOdeSolver.new()
            solver.set_tolerance(1.0e-12)
            solver.set_tau_floor_factor(factor)
            solver.set_default_constraint(Ncm.SBesselOdeConstraint.TAU)
            op = solver.create_operator(self.A, self.B, self.ELL, self.ELL)
            _, n_cols = op.solve(rhs)

            assert n_cols == int(np.ceil(factor * n_forcing))
            assert op.get_tau_constraint_order(len(rhs)) == n_cols

    def test_deriv_error_is_the_weighted_coefficient_sum(self) -> None:
        """The reported bound is sum_j j^2 |a_j| / h, on every solve path."""
        rhs, _ = self._rhs()
        half = 0.5 * (self.B - self.A)

        for constraint in (
            Ncm.SBesselOdeConstraint.DIRICHLET,
            Ncm.SBesselOdeConstraint.TAU,
        ):
            solver = Ncm.SBesselOdeSolver.new()
            solver.set_tolerance(1.0e-12)
            solver.set_default_constraint(constraint)
            op = solver.create_operator(self.A, self.B, self.ELL, self.ELL)
            coeffs, n_cols = op.solve(rhs)
            expected = np.sum(np.arange(n_cols) ** 2 * np.abs(np.array(coeffs))) / half

            assert_allclose(op.get_last_deriv_error(0), expected, rtol=1.0e-13)

            op.solve_values(rhs, self.A, self.B)
            assert_allclose(op.get_last_deriv_error(0), expected, rtol=1.0e-13)

            endpoints = np.array(op.solve_endpoints(rhs))
            assert_allclose(op.get_last_deriv_error(0), endpoints[2], rtol=1.0e-13)

    @pytest.mark.parametrize("n_ell", [1, 4, 8])
    def test_deriv_error_per_multipole(self, n_ell: int) -> None:
        """Each member of a block reports its own bound."""
        rhs, _ = self._rhs()
        half = 0.5 * (self.B - self.A)
        solver = Ncm.SBesselOdeSolver.new()
        solver.set_tolerance(1.0e-12)
        solver.set_default_constraint(Ncm.SBesselOdeConstraint.TAU)
        op = solver.create_operator(self.A, self.B, self.ELL, self.ELL + n_ell - 1)
        coeffs, n_cols = op.solve(rhs)
        block = np.array(coeffs).reshape(n_ell, n_cols)

        for i in range(n_ell):
            expected = np.sum(np.arange(n_cols) ** 2 * np.abs(block[i])) / half
            assert_allclose(op.get_last_deriv_error(i), expected, rtol=1.0e-13)

    def test_operator_matrix_follows_the_constraint(self) -> None:
        """The dense matrix carries the constraint rows the operator actually has.

        The solver-level accessor takes its constraint from the solver and so cannot express
        the pinned one; this is the entry point that can.
        """
        solver = Ncm.SBesselOdeSolver.new()
        solver.set_tolerance(1.0e-12)
        op = solver.create_operator(self.A, self.B, self.ELL, self.ELL)
        nrows = 40

        def rows01(matrix):
            data = np.array(matrix.dup_array()).reshape(matrix.nrows(), matrix.ncols())

            return np.count_nonzero(data[0]), np.count_nonzero(data[1])

        op.set_constraint(Ncm.SBesselOdeConstraint.DIRICHLET)
        assert rows01(op.get_matrix(nrows)) == (nrows, nrows)

        op.set_pinned_constraint(10, 20)
        assert rows01(op.get_matrix(nrows)) == (1, 1)

        op.set_constraint(Ncm.SBesselOdeConstraint.TAU)
        assert rows01(op.get_matrix(nrows)) == (0, 0)

    def test_operator_matrix_band_is_the_same_below_the_constraint(self) -> None:
        """Only the two constraint rows differ between constraints."""
        solver = Ncm.SBesselOdeSolver.new()
        op = solver.create_operator(self.A, self.B, self.ELL, self.ELL)

        op.set_constraint(Ncm.SBesselOdeConstraint.DIRICHLET)
        m_d = op.get_matrix(40)
        band_d = np.array(m_d.dup_array()).reshape(40, 40)[2:]

        op.set_constraint(Ncm.SBesselOdeConstraint.TAU)
        m_t = op.get_matrix(40)
        band_t = np.array(m_t.dup_array()).reshape(40, 40)[2:]

        assert_allclose(band_t, band_d, rtol=0.0, atol=0.0)


class TestDerivErrorBoundsTheFailure:
    """The reported bound has to bound the realized error, not merely correlate.

    Under the tau constraint the homogeneous content the truncation admits cancels exactly
    in the boundary functional, so what survives is the floating-point residue of that
    cancellation, and the functional weights coefficient j by j^2 / h. That is what
    get_last_deriv_error() reports, and the point of it is that it never under-reports.
    """

    ELL = 20
    X_A = 1.0e3

    def _panel(self, span: float):
        """A bump of fixed relative width, its right-hand side and a reference."""
        x_b = self.X_A + span
        mid, half = 0.5 * (self.X_A + x_b), 0.5 * span

        def bump(x):
            return np.exp(-0.5 * ((x - mid) / (0.25 * half)) ** 2)

        solver = Ncm.SBesselOdeSolver.new()
        spectral = solver.peek_spectral()
        cheb = np.array(
            spectral.compute_chebyshev_coeffs_adaptive(
                lambda _d, x: x * bump(x), self.X_A, x_b, 3, 1.0e-13, None
            )[1]
        )
        n = len(cheb)
        padded = np.concatenate([cheb, np.zeros(8)])
        kk = np.arange(n)
        geg = (
            padded[:n] / (2 * (kk + 1))
            - (kk + 2) * padded[2 : n + 2] / ((kk + 1) * (kk + 3))
            + padded[4 : n + 4] / (2 * (kk + 3))
        )
        geg[0] += 0.5 * padded[0]
        rhs = np.concatenate([[0.0, 0.0], geg])

        gx, gw = np.polynomial.legendre.leggauss(24)
        edges = np.linspace(self.X_A, x_b, max(int(np.ceil(2 * span)), 400) + 1)
        reference = sum(
            0.5
            * (e1 - e0)
            * np.sum(
                gw
                * bump(0.5 * (e0 + e1) + 0.5 * (e1 - e0) * gx)
                * spherical_jn(self.ELL, 0.5 * (e0 + e1) + 0.5 * (e1 - e0) * gx)
            )
            for e0, e1 in zip(edges[:-1], edges[1:])
        )

        return x_b, rhs, reference

    def _solve(self, x_b: float, rhs):
        solver = Ncm.SBesselOdeSolver.new()
        solver.set_tolerance(1.0e-12)
        solver.set_default_constraint(Ncm.SBesselOdeConstraint.TAU)
        op = solver.create_operator(self.X_A, x_b, self.ELL, self.ELL)
        coeffs = np.array(op.solve(rhs)[0])
        half = 0.5 * (x_b - self.X_A)
        u = np.polynomial.chebyshev.chebval([-1.0, 1.0], coeffs)
        du = (
            np.polynomial.chebyshev.chebval(
                [-1.0, 1.0], np.polynomial.chebyshev.chebder(coeffs)
            )
            / half
        )

        def term(x, u_val, du_val):
            jl = spherical_jn(self.ELL, x)
            djl = spherical_jn(self.ELL, x, derivative=True)

            return x * jl * du_val - (jl + x * djl) * u_val

        value = term(x_b, u[1], du[1]) - term(self.X_A, u[0], du[0])

        return value, op.get_last_deriv_error(0)

    @pytest.mark.parametrize("span", [30.0, 60.0, 100.0, 150.0, 250.0, 3000.0])
    def test_bound_is_never_below_the_error(self, span: float) -> None:
        """Across spans where the constraint is sound and where it fails outright."""
        x_b, rhs, reference = self._panel(span)
        value, deriv_error = self._solve(x_b, rhs)

        amplitude = max(
            abs(self.X_A * spherical_jn(self.ELL, self.X_A)),
            abs(x_b * spherical_jn(self.ELL, x_b)),
        )
        predicted = np.finfo(float).eps * deriv_error * amplitude / abs(reference)
        realized = abs(value - reference) / abs(reference)

        if realized > 1.0e-9:
            assert (
                predicted >= realized
            ), f"span {span}: bound {predicted:.2e} under the error {realized:.2e}"

    def test_it_separates_the_two_regimes(self) -> None:
        """A contaminated panel reports a bound orders above a sound one."""
        _, deriv_sound = self._solve(*self._panel(3000.0)[:2])
        _, deriv_broken = self._solve(*self._panel(30.0)[:2])

        assert deriv_broken > 1.0e6 * deriv_sound


class TestConjugatePoints:
    """Panels whose ends make the Dirichlet boundary matrix singular.

    A conjugate point is a zero of Phi = j(a) x(b) - j(b) x(a). There the two-point
    Dirichlet problem has no unique solution and the error grows as the machine epsilon
    times the boundary condition number. Neither of the other two constraints imposes an
    endpoint condition, so neither has that determinant to lose.
    """

    ELL = 2
    A0 = 200.0

    @classmethod
    def _phi(cls, span: float) -> float:
        b = cls.A0 + span

        return spherical_jn(cls.ELL, cls.A0) * spherical_yn(cls.ELL, b) - spherical_jn(
            cls.ELL, b
        ) * spherical_yn(cls.ELL, cls.A0)

    @classmethod
    def _conjugate_span(cls) -> float:
        """A span where Phi vanishes, found by root finding.

        Scanning for a maximum of the condition number is not enough: a scan of 2001
        points over 0.2 pi lands at cond ~ 1e4, where Dirichlet is still accurate to
        5e-13 and the effect looks absent.
        """
        return brentq(cls._phi, 7.9 * np.pi, 8.1 * np.pi, xtol=1.0e-13)

    @staticmethod
    def _boundary_cond(ell: int, a: float, b: float) -> float:
        mat = np.array(
            [
                [spherical_jn(ell, a), spherical_yn(ell, a)],
                [spherical_jn(ell, b), spherical_yn(ell, b)],
            ]
        )

        return np.linalg.cond(mat)

    def _rhs(self, a: float, b: float):
        solver = Ncm.SBesselOdeSolver.new()
        spectral = solver.peek_spectral()
        mid, half = 0.5 * (a + b), 0.5 * (b - a)
        cheb = np.array(
            spectral.compute_chebyshev_coeffs_adaptive(
                lambda _d, x: x * np.exp(-0.5 * ((x - mid) / (0.35 * half)) ** 2),
                a,
                b,
                3,
                1.0e-13,
                None,
            )[1]
        )
        n = len(cheb)
        padded = np.concatenate([cheb, np.zeros(8)])
        kk = np.arange(n)
        geg = (
            padded[:n] / (2 * (kk + 1))
            - (kk + 2) * padded[2 : n + 2] / ((kk + 1) * (kk + 3))
            + padded[4 : n + 4] / (2 * (kk + 3))
        )
        geg[0] += 0.5 * padded[0]

        return np.concatenate([[0.0, 0.0], geg])

    def _integral(self, a: float, b: float, constraint, pins=None) -> float:
        solver = Ncm.SBesselOdeSolver.new()
        solver.set_tolerance(1.0e-12)
        op = solver.create_operator(a, b, self.ELL, self.ELL)

        if pins is not None:
            op.set_pinned_constraint(*pins)
        else:
            op.set_constraint(constraint)

        coeffs = np.array(op.solve(self._rhs(a, b))[0])
        half = 0.5 * (b - a)
        u = np.polynomial.chebyshev.chebval([-1.0, 1.0], coeffs)
        du = (
            np.polynomial.chebyshev.chebval(
                [-1.0, 1.0], np.polynomial.chebyshev.chebder(coeffs)
            )
            / half
        )

        def term(x: float, u_val: float, du_val: float) -> float:
            jl = spherical_jn(self.ELL, x)
            djl = spherical_jn(self.ELL, x, derivative=True)

            return x * jl * du_val - (jl + x * djl) * u_val

        return term(b, u[1], du[1]) - term(a, u[0], du[0])

    def _reference(self, a: float, b: float) -> float:
        mid, half = 0.5 * (a + b), 0.5 * (b - a)
        value, _ = quad(
            lambda x: np.exp(-0.5 * ((x - mid) / (0.35 * half)) ** 2)
            * spherical_jn(self.ELL, x),
            a,
            b,
            limit=4000,
            epsabs=0.0,
            # See the note in TestSBesselTauConstraint: a tighter request buys nothing
            # on an oscillatory integrand and only reports roundoff.
            epsrel=1.0e-10,
        )

        return value

    def test_the_span_is_really_conjugate(self) -> None:
        """The root finder lands on a vanishing determinant, not merely a large one."""
        span = self._conjugate_span()

        assert abs(self._phi(span)) < 1.0e-17
        assert self._boundary_cond(self.ELL, self.A0, self.A0 + span) > 1.0e13

    def test_dirichlet_loses_digits_there(self) -> None:
        """Dirichlet error tracks the machine epsilon times the condition number."""
        a = self.A0
        b = a + self._conjugate_span()
        reference = self._reference(a, b)

        got = self._integral(a, b, Ncm.SBesselOdeConstraint.DIRICHLET)
        error = abs(got - reference) / abs(reference)

        assert (
            error > 1.0e-6
        ), f"expected the singular constraint to fail, got {error:.2e}"

    @pytest.mark.parametrize("pins", [(0, 1), (2, 3), (10, 11)])
    def test_pinned_is_immune(self, pins) -> None:
        """No endpoint condition means no determinant to degenerate."""
        a = self.A0
        b = a + self._conjugate_span()
        reference = self._reference(a, b)

        got = self._integral(a, b, None, pins=pins)

        assert_allclose(got, reference, rtol=1.0e-9)

    def test_away_from_the_conjugate_span_dirichlet_is_fine(self) -> None:
        """The failure is the determinant, not the panel."""
        a = self.A0
        b = a + self._conjugate_span() + 1.0e-2
        reference = self._reference(a, b)

        assert_allclose(
            self._integral(a, b, Ncm.SBesselOdeConstraint.DIRICHLET),
            reference,
            rtol=1.0e-9,
        )


class TestPhaseIncrementConditioning:
    """Which phase increments are safe to space knots by.

    Placing knots at the zeros of j_l would put both ends of every panel on a zero, so
    Phi = j(a) x(b) - j(b) x(a) vanishes identically and every panel is conjugate. The
    extrema are no better, being a further half period apart in the same sense. What
    decides it is the increment: a multiple of pi is singular, an odd multiple of pi/2 is
    where sin(dtheta) is one.
    """

    @staticmethod
    def _theta(ell: int, x: float) -> float:
        """WKB phase above the turning point."""
        nu = np.sqrt(ell * (ell + 1.0))

        return np.sqrt(x * x - nu * nu) - nu * np.arccos(nu / x)

    @classmethod
    def _x_at_phase(cls, ell: int, target: float, lo: float) -> float:
        return brentq(
            lambda x: cls._theta(ell, x) - target, lo + 1.0e-6, lo + 200.0, xtol=1.0e-12
        )

    @staticmethod
    def _cond(ell: int, a: float, b: float) -> float:
        mat = np.array(
            [
                [spherical_jn(ell, a), spherical_yn(ell, a)],
                [spherical_jn(ell, b), spherical_yn(ell, b)],
            ]
        )

        return np.linalg.cond(mat)

    @pytest.mark.parametrize("ell", [2, 20, 200])
    def test_multiples_of_pi_are_singular_and_half_odd_are_not(self, ell: int) -> None:
        """Spacing by pi is conjugate; spacing by pi/2 or 3pi/2 is well conditioned."""
        x0 = 3.0 * np.sqrt(ell * (ell + 1.0))
        theta0 = self._theta(ell, x0)

        worst_half_odd = 0.0
        best_multiple = np.inf
        for increment, is_multiple in (
            (np.pi, True),
            (2.0 * np.pi, True),
            (0.5 * np.pi, False),
            (1.5 * np.pi, False),
        ):
            x1 = self._x_at_phase(ell, theta0 + increment, x0)
            value = self._cond(ell, x0, x1)

            if is_multiple:
                best_multiple = min(best_multiple, value)
            else:
                worst_half_odd = max(worst_half_odd, value)

        assert worst_half_odd < 10.0, "odd multiples of pi/2 must be well conditioned"
        assert best_multiple > 1.0e2, "multiples of pi must be near conjugate"
        assert best_multiple > 50.0 * worst_half_odd

    def test_zeros_of_jl_are_the_worst_knots(self) -> None:
        """Both ends on a zero makes the determinant vanish identically."""
        ell = 20
        nu = np.sqrt(ell * (ell + 1.0))
        zeros, x, step = [], 3.0 * nu, 0.05
        previous = spherical_jn(ell, x)

        while len(zeros) < 3:
            nxt = x + step
            current = spherical_jn(ell, nxt)

            if previous * current < 0.0:
                zeros.append(
                    brentq(lambda t: spherical_jn(ell, t), x, nxt, xtol=1.0e-12)
                )

            x, previous = nxt, current

        for a, b in zip(zeros[:-1], zeros[1:]):
            assert self._cond(ell, a, b) > 1.0e12

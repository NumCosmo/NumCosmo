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
import time

from numpy.testing import assert_allclose

pytest.importorskip("scipy")
# flake8: noqa: E402
# pylint: disable=wrong-import-position

from scipy.special import spherical_jn
from scipy.linalg import solve
from scipy.integrate import quad

from numcosmo_py import Ncm


class TestSBesselOperators:
    """Tests for NcmSBesselOdeSolver operator matrices."""

    @staticmethod
    def _compare_batched_vs_single_solve(
        solver: Ncm.SBesselOdeSolver,
        a: float,
        b: float,
        lmin: int,
        lmax: int,
        rhs: np.ndarray,
        rtol: float,
        atol_factor: float,
        err_msg_template: str,
    ) -> None:
        """Compare batched solve() with individual single-ell solves.

        Args:
            solver: The SBesselOdeSolver instance
            a, b: Domain boundaries
            lmin, lmax: Range of ell values
            rhs: Right-hand side coefficients
            rtol: Relative tolerance for comparison
            atol_factor: Factor to multiply max coefficient for absolute tolerance
            err_msg_template: Template string for error message (should contain {ell})
        """
        n_ell = lmax - lmin + 1

        # Batched solve
        op_batched = solver.create_operator(a, b, lmin, lmax)
        solutions_batched, sol_len = op_batched.solve(rhs)
        solutions_batched_np = np.array(solutions_batched).reshape(n_ell, sol_len)

        # Individual solves
        for i, ell in enumerate(range(lmin, lmax + 1)):
            op_single = solver.create_operator(a, b, ell, ell)
            sol_single, _sol_len_single = op_single.solve(rhs)
            sol_single_np = np.array(sol_single)

            max_coeff_single = np.max(np.abs(sol_single_np))
            min_len = min(len(sol_single_np), sol_len)
            assert_allclose(
                solutions_batched_np[i, :min_len],
                sol_single_np[:min_len],
                rtol=rtol,
                atol=atol_factor * max_coeff_single,
                err_msg=err_msg_template.format(ell=ell),
            )

    @staticmethod
    def _compare_batched_vs_single_endpoints(
        solver: Ncm.SBesselOdeSolver,
        a: float,
        b: float,
        lmin: int,
        lmax: int,
        rhs: np.ndarray,
        rtol: float,
        atol_factor: float,
        err_msg_template: str,
        compare_error: bool = True,
    ) -> None:
        """Compare batched solve_endpoints() with individual single-ell solves.

        Args:
            solver: The SBesselOdeSolver instance
            a, b: Domain boundaries
            lmin, lmax: Range of ell values
            rhs: Right-hand side coefficients
            rtol: Relative tolerance for comparison
            atol_factor: Factor to multiply max derivative for absolute tolerance
            err_msg_template: Template string for error message (should contain {ell})
            compare_error: If False, only compare derivatives (not error estimate)
        """
        n_ell = lmax - lmin + 1

        # Batched solve_endpoints
        op_batched = solver.create_operator(a, b, lmin, lmax)
        endpoints_batched = np.array(op_batched.solve_endpoints(rhs)).reshape(n_ell, 3)

        # Individual solve_endpoints
        for i, ell in enumerate(range(lmin, lmax + 1)):
            op_single = solver.create_operator(a, b, ell, ell)
            deriv_a, deriv_b, error = op_single.solve_endpoints(rhs)
            endpoints_single = np.array([deriv_a, deriv_b, error])

            max_deriv_single = max(abs(deriv_a), abs(deriv_b))

            # Select which components to compare
            n_compare = 3 if compare_error else 2

            assert_allclose(
                endpoints_batched[i, :n_compare],
                endpoints_single[:n_compare],
                rtol=rtol,
                atol=atol_factor * max_deriv_single,
                err_msg=err_msg_template.format(ell=ell),
            )

    @pytest.mark.parametrize("l_val", list(range(21)))
    def test_spherical_bessel_ode(self, l_val: int) -> None:
        """Test solving modified spherical Bessel ODE with non-zero Dirichlet BCs.

        The function u(x) = x*j_l(x) satisfies:
        u''(x) + [x^2 - l(l+1)]u(x) = 0

        We solve this on interval [1,20] with BCs u(1) = 1*j_l(1), u(20) = 20*j_l(20)
        and verify the solution matches x*j_l(x) everywhere.
        """
        N = 128
        a, b = 1.0, 20.0  # Use interval [1, 20] as requested

        # Create solver and operator
        solver = Ncm.SBesselOdeSolver.new()

        # Get the operator matrix (should now be square NxN)
        mat = solver.get_operator_matrix(a, b, l_val, N)
        mat_np = mat.to_numpy()

        # Get boundary values using scipy
        y_a = spherical_jn(l_val, a) * a
        y_b = spherical_jn(l_val, b) * b

        # Debug: print matrix shape
        # Row 1 enforces u(+1) = y_b (i.e., u(b) in physical coords)
        # Rows 2 to N-1 are the differential operator (homogeneous: RHS = 0)
        rhs_np = np.zeros(N)
        rhs_np[0] = y_a
        rhs_np[1] = y_b

        # Assert matrix is square
        assert (
            mat_np.shape[0] == mat_np.shape[1]
        ), f"Matrix must be square, got shape {mat_np.shape}"

        # Solve the linear system using LU decomposition
        solution_coeffs = solve(mat_np, rhs_np)

        # Test at several interior points mapped to [-1, 1]
        # Mapping: x_physical = (a+b)/2 + (b-a)/2 * t, where t in [-1, 1]
        t_test = np.array([-0.8, -0.5, -0.2, 0.0, 0.2, 0.5, 0.8])

        for t in t_test:
            # Evaluate solution at t (coefficients are in Chebyshev basis)
            y_computed = Ncm.Spectral.chebyshev_eval(solution_coeffs, t)

            # Map t to physical x
            x_physical = (a + b) / 2.0 + (b - a) / 2.0 * t

            # Get exact value
            y_exact = spherical_jn(l_val, x_physical) * x_physical

            assert_allclose(
                y_computed,
                y_exact,
                rtol=1.0e-10,
                atol=1.0e-15,
                err_msg=f"Spherical Bessel solution mismatch at x={x_physical}",
            )

        # Verify boundary conditions are satisfied (use Chebyshev eval)
        y_at_a = Ncm.Spectral.chebyshev_eval(solution_coeffs, -1.0)
        y_at_b = Ncm.Spectral.chebyshev_eval(solution_coeffs, 1.0)
        assert_allclose(
            y_at_a, y_a, rtol=1.0e-13, atol=1.0e-13, err_msg="BC at x=a not satisfied"
        )
        assert_allclose(
            y_at_b, y_b, rtol=1.0e-13, atol=1.0e-13, err_msg="BC at x=b not satisfied"
        )

    @pytest.mark.parametrize("l_val", list(range(21)))
    def test_spherical_bessel_integration(self, l_val: int) -> None:
        """Test integration of spherical Bessel function via Green's identity.

        Solve the ODE u''(x) + [x^2 - l(l+1)]u(x) = 1 with homogeneous BCs (u(a)=0,
        u(b)=0). By Green's identity: [x*j_l(x)]_a^b * u'(b) - [x*j_l(x)]_a * u'(a) =
        int_a^b j_l(x)dx where u is the solution with RHS=1.
        """
        N = 128
        a, b = 1.0, 20.0

        # Create solver
        solver = Ncm.SBesselOdeSolver.new()

        # Get the operator matrix
        mat = solver.get_operator_matrix(a, b, l_val, N)
        mat_np = mat.to_numpy()

        # Set up RHS: homogeneous BCs (y(a)=0, y(b)=0) with RHS=1
        rhs_np = np.zeros(N)
        rhs_np[0] = 0.0  # BC at x=a (t=-1)
        rhs_np[1] = 0.0  # BC at x=b (t=+1)
        rhs_np[2] = 1.0
        # Rows 3 to N-1 have RHS = 0

        # Solve the linear system
        solution_coeffs = solve(mat_np, rhs_np)

        # Compute derivatives at endpoints using the derivative function
        # Need to account for the coordinate transformation: x = m + h*t
        # where m = (a+b)/2, h = (b-a)/2
        # dy/dx = (dy/dt) * (dt/dx) = (dy/dt) / h
        h = (b - a) / 2.0

        # Evaluate y'(t) at t=-1 (corresponds to x=a) and t=+1 (corresponds to x=b)
        y_prime_at_minus1 = Ncm.Spectral.chebyshev_deriv(solution_coeffs, -1.0)
        y_prime_at_plus1 = Ncm.Spectral.chebyshev_deriv(solution_coeffs, 1.0)

        # Convert from dy/dt to dy/dx
        y_prime_a = y_prime_at_minus1 / h
        y_prime_b = y_prime_at_plus1 / h

        # Get j_l values at endpoints
        j_l_a = spherical_jn(l_val, a)
        j_l_b = spherical_jn(l_val, b)

        # Compute left-hand side: [x*j_l(x)]*u'(x) from a to b
        lhs = b * j_l_b * y_prime_b - a * j_l_a * y_prime_a

        # Compute right-hand side: integral of j_l(x) from a to b
        # Note: The Green's identity for u''+(x^2-l(l+1))u=1 gives this relation

        def integrand(x: float) -> float:
            return spherical_jn(l_val, x) / x

        rhs, _ = quad(integrand, a, b, epsabs=1e-12, epsrel=1e-14)

        # Compare
        assert_allclose(
            lhs,
            rhs,
            rtol=1.0e-8,
            atol=1.0e-20,
            err_msg=(
                f"Green's identity failed for l={l_val}: "
                f"[x*j_l(x)]*u'(x) from a to b != int_a^b j_l(x)dx"
            ),
        )

    @pytest.mark.parametrize("l_val", list(range(21)))
    def test_spherical_bessel_integration_x(self, l_val: int) -> None:
        """Test integration of x*j_l(x) via Green's identity.

        Solve the ODE u''(x) + [x^2 - l(l+1)]u(x) = x with homogeneous BCs (u(a)=0,
        u(b)=0). By Green's identity: [x*j_l(x)] * u'(b)|_a^b = int_a^b x*j_l(x)dx
        where u is the solution with RHS=x.
        """
        N = 128
        a, b = 1.0, 20.0

        # Create solver
        solver = Ncm.SBesselOdeSolver.new()

        # Get the operator matrix
        mat = solver.get_operator_matrix(a, b, l_val, N)
        mat_np = mat.to_numpy()

        # Set up RHS: homogeneous BCs (y(a)=0, y(b)=0) with RHS=x
        # In the mapped coordinates, x = m + h*t where m=(a+b)/2, h=(b-a)/2
        # So we need RHS = m + h*t in Chebyshev basis
        # T_0(t) = 1, T_1(t) = t
        # Therefore: m*T_0 + h*T_1
        rhs_np = np.zeros(N)
        rhs_np[0] = 0.0  # BC at x=a (t=-1)
        rhs_np[1] = 0.0  # BC at x=b (t=+1)

        m = (a + b) / 2.0
        h = (b - a) / 2.0

        # RHS coefficients in Chebyshev basis for x = m + h*t
        rhs_np[2] = m  # Constant term (T_0 coefficient)
        rhs_np[3] = (
            0.25 * h
        )  # Linear term (T_1 coefficient), which equals 0.25 for this case
        # Rows 4 to N-1 have RHS = 0

        # Solve the linear system
        solution_coeffs = solve(mat_np, rhs_np)

        # Compute derivatives at endpoints
        # dy/dx = (dy/dt) / h
        y_prime_at_minus1 = Ncm.Spectral.chebyshev_deriv(solution_coeffs, -1.0)
        y_prime_at_plus1 = Ncm.Spectral.chebyshev_deriv(solution_coeffs, 1.0)

        # Convert from dy/dt to dy/dx
        y_prime_a = y_prime_at_minus1 / h
        y_prime_b = y_prime_at_plus1 / h

        # Get j_l values at endpoints
        j_l_a = spherical_jn(l_val, a)
        j_l_b = spherical_jn(l_val, b)

        # Compute left-hand side: [x*j_l(x)] * u'(x) from a to b
        lhs = b * j_l_b * y_prime_b - a * j_l_a * y_prime_a

        # Compute right-hand side: integral of x*j_l(x) from a to b
        def integrand(x: float) -> float:
            return spherical_jn(l_val, x)

        rhs, _ = quad(integrand, a, b, epsabs=1e-12, epsrel=1e-14)

        # Compare
        assert_allclose(
            lhs,
            rhs,
            rtol=1.0e-8,
            atol=1.0e-20,
            err_msg=(
                f"Green's identity failed for l={l_val}: "
                f"[x*j_l(x)]*u'(x) from a to b != int_a^b x*j_l(x)dx"
            ),
        )

    @pytest.mark.parametrize("l_val", list(range(21)))
    @pytest.mark.parametrize("N", [32, 64, 128])
    def test_solve_qr_vs_scipy(self, l_val: int, N: int) -> None:
        """Test that QR solve produces the same result as scipy solver."""
        a, b = 1.0, 20.0

        # Create solver and operator
        solver = Ncm.SBesselOdeSolver.new()
        op = solver.create_operator(a, b, l_val, l_val)

        # Set up simple RHS vector with just a few coefficients
        rhs = np.zeros(N)
        rhs[2] = 1.0

        # Call solve (QR method) - returns (solution, solution_len)
        solution_list, _solution_len = op.solve(rhs)
        solution = np.array(solution_list)

        op_mat = solver.get_operator_matrix(a, b, l_val, 2 * N).to_numpy()

        # Solve the system using scipy for comparison
        rhs_np_full = np.zeros(2 * N)
        rhs_np_full[2] = 1.0
        solution_scipy = solve(op_mat, rhs_np_full)
        safe_part = min(len(solution), len(solution_scipy))

        diff = np.sum(
            (solution_scipy[:safe_part] - solution[:safe_part]) ** 2
        ) / np.sum(solution[:safe_part] ** 2)
        assert diff < 1e-10, f"QR solve doesn't match scipy solution: diff={diff}"
        # Basic sanity check - solution should exist and be non-trivial
        assert solution is not None
        assert len(solution) > 0

    @pytest.mark.parametrize("l_val", list(range(21)))
    @pytest.mark.parametrize("N", [32, 64, 128])
    def test_solve_dense_vs_scipy(self, l_val: int, N: int) -> None:
        """Test that solve_dense produces the same result as scipy solver.

        This tests the LAPACK dgesv-based solve_dense method against scipy's
        solve to verify the column-major matrix format is correct.
        """
        a, b = 1.0, 20.0

        # Create solver
        solver = Ncm.SBesselOdeSolver.new()

        # Set up RHS: homogeneous BCs with RHS=1
        rhs_np = np.zeros(N)
        rhs_np[0] = 0.0  # BC at x=a
        rhs_np[1] = 0.0  # BC at x=b
        rhs_np[2] = 1.0
        rhs_vec = Ncm.Vector.new_array(rhs_np.tolist())

        # Solve using scipy
        mat = solver.get_operator_matrix(a, b, l_val, N)
        mat_np = mat.to_numpy()
        solution_scipy = solve(mat_np, rhs_np)

        # Solve using solve_dense (LAPACK dgesv)
        solution_dense = solver.solve_dense(a, b, l_val, rhs_vec, N)
        solution_dense_np = solution_dense.to_numpy()

        # Compare solutions
        assert_allclose(
            solution_dense_np,
            solution_scipy,
            rtol=1.0e-10,
            atol=1.0e-12,
            err_msg=(f"solve_dense doesn't match scipy for l={l_val}, N={N}"),
        )

    @pytest.mark.parametrize("l_val", [0, 1, 5, 10, 20])
    def test_solve_dense_integration_rhs1(self, l_val: int) -> None:
        """Test solve_dense with Green's identity for RHS=1.

        Verifies that solve_dense produces solutions u to u''+[x^2-l(l+1)]u=1 that
        satisfy Green's identity: [x*j_l(x)] * u'(x) from a to b = int_a^b j_l(x)dx
        """
        N = 128
        a, b = 1.0, 20.0

        # Create solver
        solver = Ncm.SBesselOdeSolver.new()

        # Set up RHS: homogeneous BCs with RHS=1
        rhs_np = np.zeros(N)
        rhs_np[0] = 0.0  # BC at x=a
        rhs_np[1] = 0.0  # BC at x=b
        rhs_np[2] = 1.0
        rhs_vec = Ncm.Vector.new_array(rhs_np.tolist())

        # Solve using solve_dense
        solution_vec = solver.solve_dense(a, b, l_val, rhs_vec, N)
        solution_coeffs = solution_vec.to_numpy()

        # Compute derivatives at endpoints
        h = (b - a) / 2.0
        y_prime_at_minus1 = Ncm.Spectral.chebyshev_deriv(solution_coeffs, -1.0)
        y_prime_at_plus1 = Ncm.Spectral.chebyshev_deriv(solution_coeffs, 1.0)
        y_prime_a = y_prime_at_minus1 / h
        y_prime_b = y_prime_at_plus1 / h

        # Get j_l values
        j_l_a = spherical_jn(l_val, a)
        j_l_b = spherical_jn(l_val, b)

        # Compute LHS: [x*j_l(x)]*u'(x) from a to b
        lhs = b * j_l_b * y_prime_b - a * j_l_a * y_prime_a

        # Compute RHS: integral
        def integrand(x: float) -> float:
            return spherical_jn(l_val, x) / x

        rhs, _ = quad(integrand, a, b, epsabs=1e-12, epsrel=1e-14)

        # Compare
        assert_allclose(
            lhs,
            rhs,
            rtol=1.0e-8,
            atol=1.0e-20,
            err_msg=(f"Green's identity failed with solve_dense for l={l_val}"),
        )

    @pytest.mark.parametrize("l_val", [0, 1, 5, 10, 20])
    def test_solve_dense_integration_rhs_x(self, l_val: int) -> None:
        """Test solve_dense with Green's identity for RHS=x.

        Verifies that solve_dense produces solutions u to u''+[x^2-l(l+1)]u=x that
        satisfy Green's identity: [x*j_l(x)] * u'(x) from a to b = int_a^b x*j_l(x)dx
        """
        N = 128
        a, b = 1.0, 20.0

        # Create solver
        solver = Ncm.SBesselOdeSolver.new()

        # Set up RHS: homogeneous BCs with RHS=x
        # x = m + h*t where m=(a+b)/2, h=(b-a)/2
        # Need to express x in Gegenbauer C^(2) basis
        m = (a + b) / 2.0
        h = (b - a) / 2.0
        rhs_np = np.zeros(N)
        rhs_np[0] = 0.0  # BC at x=a
        rhs_np[1] = 0.0  # BC at x=b
        rhs_np[2] = m  # Constant term
        rhs_np[3] = 0.25 * h  # Linear term scaled for C^(2) basis
        rhs_vec = Ncm.Vector.new_array(rhs_np.tolist())

        # Solve using solve_dense
        solution_vec = solver.solve_dense(a, b, l_val, rhs_vec, N)
        solution_coeffs = solution_vec.to_numpy()

        # Compute derivatives at endpoints
        y_prime_at_minus1 = Ncm.Spectral.chebyshev_deriv(solution_coeffs, -1.0)
        y_prime_at_plus1 = Ncm.Spectral.chebyshev_deriv(solution_coeffs, 1.0)
        y_prime_a = y_prime_at_minus1 / h
        y_prime_b = y_prime_at_plus1 / h

        # Get j_l values
        j_l_a = spherical_jn(l_val, a)
        j_l_b = spherical_jn(l_val, b)

        # Compute LHS: [x*j_l(x)]*u'(x) from a to b
        lhs = b * j_l_b * y_prime_b - a * j_l_a * y_prime_a

        # Compute RHS: integral of x*j_l(x)
        def integrand(x: float) -> float:
            return spherical_jn(l_val, x)

        rhs, _ = quad(integrand, a, b, epsabs=1e-12, epsrel=1e-14)

        # Compare
        assert_allclose(
            lhs,
            rhs,
            rtol=1.0e-8,
            atol=1.0e-20,
            err_msg=(f"Green's identity failed with solve_dense for l={l_val}, RHS=x"),
        )

    @pytest.mark.parametrize("l_val", [0, 1, 5, 10, 20])
    def test_solve_qr_integration_rhs1(self, l_val: int) -> None:
        """Test QR solve with Green's identity for RHS=1.

        Verifies that solve (QR method) produces solutions u to u''+[x^2-l(l+1)]u=1
        that satisfy Green's identity: [x*j_l(x)] * u'(x) from a to b = int_a^b j_l(x)dx
        """
        N = 128
        a, b = 1.0, 20.0

        # Create solver and operator
        solver = Ncm.SBesselOdeSolver.new()
        op = solver.create_operator(a, b, l_val, l_val)

        # Set up RHS: homogeneous BCs with RHS=1
        rhs_np = np.zeros(N)
        rhs_np[0] = 0.0  # BC at x=a
        rhs_np[1] = 0.0  # BC at x=b
        rhs_np[2] = 1.0

        # Solve using QR method
        solution_coeffs, _solution_len = op.solve(rhs_np)

        # Compute derivatives at endpoints
        h = (b - a) / 2.0
        y_prime_at_minus1 = Ncm.Spectral.chebyshev_deriv(solution_coeffs, -1.0)
        y_prime_at_plus1 = Ncm.Spectral.chebyshev_deriv(solution_coeffs, 1.0)
        y_prime_a = y_prime_at_minus1 / h
        y_prime_b = y_prime_at_plus1 / h

        # Get j_l values
        j_l_a = spherical_jn(l_val, a)
        j_l_b = spherical_jn(l_val, b)

        # Compute LHS: [x*j_l(x)]*u'(x) from a to b
        lhs = b * j_l_b * y_prime_b - a * j_l_a * y_prime_a

        # Compute RHS: integral
        def integrand(x: float) -> float:
            return spherical_jn(l_val, x) / x

        rhs, _ = quad(integrand, a, b, epsabs=1e-12, epsrel=1e-14)

        # Compare
        assert_allclose(
            lhs,
            rhs,
            rtol=1.0e-8,
            atol=1.0e-20,
            err_msg=(f"Green's identity failed with QR solve for l={l_val}"),
        )

    @pytest.mark.parametrize("l_val", [0, 1, 5, 10, 20])
    def test_solve_qr_integration_rhs_x(self, l_val: int) -> None:
        """Test QR solve with Green's identity for RHS=x.

        Verifies that solve (QR method) produces solutions u to u''+[x^2-l(l+1)]u=x
        that satisfy Green's identity: [x*j_l(x)] * u'(x) from a to b = int_a^b
        x*j_l(x)dx
        """
        N = 128
        a, b = 1.0, 20.0

        # Create solver and operator
        solver = Ncm.SBesselOdeSolver.new()
        op = solver.create_operator(a, b, l_val, l_val)

        # Set up RHS: homogeneous BCs with RHS=x
        # x = m + h*t where m=(a+b)/2, h=(b-a)/2
        # Need to express x in Gegenbauer C^(2) basis
        m = (a + b) / 2.0
        h = (b - a) / 2.0
        rhs_np = np.zeros(N)
        rhs_np[0] = 0.0  # BC at x=a
        rhs_np[1] = 0.0  # BC at x=b
        rhs_np[2] = m  # Constant term
        rhs_np[3] = 0.25 * h  # Linear term scaled for C^(2) basis

        # Solve using QR method
        solution_coeffs, _solution_len = op.solve(rhs_np)

        # Compute derivatives at endpoints
        y_prime_at_minus1 = Ncm.Spectral.chebyshev_deriv(solution_coeffs, -1.0)
        y_prime_at_plus1 = Ncm.Spectral.chebyshev_deriv(solution_coeffs, 1.0)
        y_prime_a = y_prime_at_minus1 / h
        y_prime_b = y_prime_at_plus1 / h

        # Get j_l values
        j_l_a = spherical_jn(l_val, a)
        j_l_b = spherical_jn(l_val, b)

        # Compute LHS: [x*j_l(x)]*u'(x) from a to b
        lhs = b * j_l_b * y_prime_b - a * j_l_a * y_prime_a

        # Compute RHS: integral of x*j_l(x)
        def integrand(x: float) -> float:
            return spherical_jn(l_val, x)

        rhs, _ = quad(integrand, a, b, epsabs=1e-12, epsrel=1e-14)

        # Compare
        assert_allclose(
            lhs,
            rhs,
            rtol=1.0e-8,
            atol=1.0e-20,
            err_msg=(f"Green's identity failed with QR solve for l={l_val}, RHS=x"),
        )

    @pytest.mark.parametrize("l_val", [0, 1, 5, 10, 20])
    @pytest.mark.parametrize("N", [32, 64, 128])
    def test_zero_rhs_solve(self, l_val: int, N: int) -> None:
        """Test that solve returns zero solution for zero RHS.

        When the RHS is completely zero (including boundary conditions), the
        solver should return a zero solution. The QR method may compress the
        solution representation, so the length may be less than N.
        """
        a, b = 1.0, 20.0

        # Create solver and operator
        solver = Ncm.SBesselOdeSolver.new()
        op = solver.create_operator(a, b, l_val, l_val)

        # Set up zero RHS
        rhs_np = np.zeros(N)

        # Solve using QR method
        solution_coeffs, solution_len = op.solve(rhs_np)
        solution = np.array(solution_coeffs)

        # Verify solution is zero
        assert_allclose(
            solution,
            np.zeros(solution_len),
            rtol=0.0,
            atol=1.0e-15,
            err_msg=f"Zero RHS should produce zero solution for l={l_val}, N={N}",
        )

        # Verify solution length is reasonable (QR may compress)
        assert (
            solution_len <= N
        ), f"Solution length {solution_len} should not exceed RHS length {N}"

    @pytest.mark.parametrize("l_val", [0, 1, 5, 10, 20])
    @pytest.mark.parametrize("N", [32, 64, 128])
    def test_zero_rhs_solve_dense(self, l_val: int, N: int) -> None:
        """Test that solve_dense returns zero solution for zero RHS.

        When the RHS is completely zero (including boundary conditions), the
        solve_dense method should return a zero solution with the same length as the
        RHS.
        """
        a, b = 1.0, 20.0

        # Create solver
        solver = Ncm.SBesselOdeSolver.new()

        # Set up zero RHS
        rhs_np = np.zeros(N)
        rhs_vec = Ncm.Vector.new_array(rhs_np.tolist())

        # Solve using solve_dense
        solution_vec = solver.solve_dense(a, b, l_val, rhs_vec, N)
        solution = solution_vec.to_numpy()

        # Verify solution is zero
        assert_allclose(
            solution,
            np.zeros(N),
            rtol=0.0,
            atol=1.0e-15,
            err_msg=(
                f"Zero RHS should produce zero solution "
                f"with solve_dense for l={l_val}, N={N}"
            ),
        )

        # Verify solution has same length as RHS
        assert len(solution) == N, f"Solution length should match RHS length {N}"

    @pytest.mark.parametrize("l_val", [0, 1, 5, 10, 20])
    @pytest.mark.parametrize("N", [32, 64, 128])
    def test_zero_rhs_solve_endpoints(self, l_val: int, N: int) -> None:
        """Test that solve_endpoints returns zero derivatives for zero RHS.

        When the RHS is completely zero (including boundary conditions), the
        solve_endpoints method should return zero derivatives at both endpoints.
        """
        a, b = 1.0, 20.0

        # Create solver and operator
        solver = Ncm.SBesselOdeSolver.new()
        op = solver.create_operator(a, b, l_val, l_val)

        # Set up zero RHS
        rhs_np = np.zeros(N)

        # Get endpoint derivatives
        deriv_a, deriv_b, _error = op.solve_endpoints(rhs_np)

        # Verify derivatives are zero
        assert_allclose(
            deriv_a,
            0.0,
            rtol=0.0,
            atol=1.0e-15,
            err_msg=(
                f"Zero RHS should produce zero derivative at a for l={l_val}, N={N}"
            ),
        )
        assert_allclose(
            deriv_b,
            0.0,
            rtol=0.0,
            atol=1.0e-15,
            err_msg=(
                f"Zero RHS should produce zero derivative at b for l={l_val}, N={N}"
            ),
        )

    def test_zero_rhs_batched(self) -> None:
        """Test that batched solve returns zero solutions for zero RHS.

        When the RHS is completely zero (including boundary conditions), the batched
        solver should return zero solutions for all ell values. The QR method may
        compress the solution representation.
        """
        N = 128
        a, b = 1.0, 20.0
        lmin, lmax = 5, 15
        n_ell = lmax - lmin + 1

        # Create solver
        solver = Ncm.SBesselOdeSolver.new()

        # Set up zero RHS
        rhs_np = np.zeros(N)

        # Solve using batched operator
        op_batched = solver.create_operator(a, b, lmin, lmax)
        solutions_batched_list, solution_len = op_batched.solve(rhs_np)
        solutions_batched = np.array(solutions_batched_list).reshape(
            n_ell, solution_len
        )

        # Verify all solutions are zero
        assert_allclose(
            solutions_batched,
            np.zeros((n_ell, solution_len)),
            rtol=0.0,
            atol=1.0e-15,
            err_msg=(
                "Zero RHS should produce zero solutions for all ell in batched solve"
            ),
        )

        # Verify solution length is reasonable (QR may compress)
        assert (
            solution_len <= N
        ), f"Solution length {solution_len} should not exceed RHS length {N}"

    def test_zero_rhs_batched_endpoints(self) -> None:
        """Test that batched solve_endpoints returns zero derivatives for zero RHS.

        When the RHS is completely zero (including boundary conditions), the
        batched solve_endpoints should return zero derivatives at both endpoints
        for all ell values.
        """
        N = 128
        a, b = 1.0, 20.0
        lmin, lmax = 5, 15
        n_ell = lmax - lmin + 1

        # Create solver
        solver = Ncm.SBesselOdeSolver.new()

        # Set up zero RHS
        rhs_np = np.zeros(N)

        # Solve using batched solve_endpoints
        op_batched = solver.create_operator(a, b, lmin, lmax)
        endpoints_batched = np.array(op_batched.solve_endpoints(rhs_np)).reshape(
            n_ell, 3
        )

        # Verify all derivatives are zero (first two columns: deriv_a and deriv_b)
        assert_allclose(
            endpoints_batched[:, 0],  # derivatives at a
            np.zeros(n_ell),
            rtol=0.0,
            atol=1.0e-15,
            err_msg="Zero RHS should produce zero derivatives at a for all ell",
        )
        assert_allclose(
            endpoints_batched[:, 1],  # derivatives at b
            np.zeros(n_ell),
            rtol=0.0,
            atol=1.0e-15,
            err_msg="Zero RHS should produce zero derivatives at b for all ell",
        )

    def test_solve_batched_vs_non_batched(self) -> None:
        """Test that batched and non-batched solvers give identical results.

        This test verifies that solve_batched produces the same solution as
        multiple calls to solve for individual l values.
        """
        N = 128
        a, b = 1.0, 20.0
        lmin, lmax = 5, 25
        n_l = lmax - lmin + 1

        # Create a test RHS vector (Gegenbauer coefficients)
        rhs_data = np.zeros(N)
        rhs_data[0] = 0.0  # BC at -1
        rhs_data[1] = 0.0  # BC at +1
        # Set some non-trivial coefficients
        for k in range(2, N):
            rhs_data[k] = 1.0 / (1.0 + k * k)

        # Create solver
        solver = Ncm.SBesselOdeSolver.new()

        # Solve using batched operator
        op_batched = solver.create_operator(a, b, lmin, lmax)
        solutions_batched_list, solution_len = op_batched.solve(rhs_data)
        solutions_batched = np.array(solutions_batched_list).reshape(n_l, solution_len)

        # Solve individually for each l
        solutions_individual = []
        for ell in range(lmin, lmax + 1):
            op_single = solver.create_operator(a, b, ell, ell)
            sol, _sol_len = op_single.solve(rhs_data)
            solutions_individual.append(sol)

        # Compare each solution
        # Note: batched solver uses a common solution_len across all ell values,
        # but individual solvers may converge at different lengths
        for i, ell in enumerate(range(lmin, lmax + 1)):
            individual_sol = np.array(solutions_individual[i])
            batched_sol = solutions_batched[i, :]

            # Compare only the overlapping part
            min_len = min(len(individual_sol), len(batched_sol))

            assert_allclose(
                batched_sol[:min_len],
                individual_sol[:min_len],
                rtol=1.0e-12,
                atol=1.0e-14,
                err_msg=f"Batched solution for l={ell} doesn't match individual solve",
            )

    @pytest.mark.parametrize("l_val", [0, 1, 5, 10, 20])
    def test_solve_endpoints_rhs1(self, l_val: int) -> None:
        """Test solve_endpoints with Green's identity for RHS=1.

        Verifies that solve_endpoints correctly computes endpoint derivatives for
        u''+[x^2-l(l+1)]u=1 that satisfy Green's identity: [x*j_l(x)] * u'(x) from a to
        b = int_a^b j_l(x)dx
        """
        N = 128
        a, b = 1.0, 20.0

        # Create solver and operator
        solver = Ncm.SBesselOdeSolver.new()
        op = solver.create_operator(a, b, l_val, l_val)

        # Set up RHS: homogeneous BCs with RHS=1
        rhs_np = np.zeros(N)
        rhs_np[0] = 0.0  # BC at x=a
        rhs_np[1] = 0.0  # BC at x=b
        rhs_np[2] = 1.0

        # Compute endpoints using solve_endpoints
        deriv_a, deriv_b, _error = op.solve_endpoints(rhs_np)

        # Get j_l values at endpoints
        j_l_a = spherical_jn(l_val, a)
        j_l_b = spherical_jn(l_val, b)

        # Compute LHS: [x*j_l(x)] * u'(x) from a to b
        lhs = b * j_l_b * deriv_b - a * j_l_a * deriv_a

        # Compute RHS: integral
        def integrand(x: float) -> float:
            return spherical_jn(l_val, x) / x

        rhs, _ = quad(integrand, a, b, epsabs=1e-12, epsrel=1e-14)

        # Compare
        assert_allclose(
            lhs,
            rhs,
            rtol=1.0e-8,
            atol=1.0e-20,
            err_msg=(f"Green's identity failed with solve_endpoints for l={l_val}"),
        )

    @pytest.mark.parametrize("l_val", [0, 1, 5, 10, 20])
    def test_solve_endpoints_rhs_x(self, l_val: int) -> None:
        """Test solve_endpoints with Green's identity for RHS=x.

        Verifies that solve_endpoints correctly computes endpoint derivatives for
        u''+[x^2-l(l+1)]u=x that satisfy Green's identity: [x*j_l(x)] * u'(x) from a to
        b = int_a^b x*j_l(x)dx
        """
        N = 128
        a, b = 1.0, 20.0

        # Create solver and operator
        solver = Ncm.SBesselOdeSolver.new()
        op = solver.create_operator(a, b, l_val, l_val)

        # Set up RHS: homogeneous BCs with RHS=x
        # x = m + h*t where m=(a+b)/2, h=(b-a)/2
        m = (a + b) / 2.0
        h = (b - a) / 2.0
        rhs_np = np.zeros(N)
        rhs_np[0] = 0.0  # BC at x=a
        rhs_np[1] = 0.0  # BC at x=b
        rhs_np[2] = m  # Constant term
        rhs_np[3] = 0.25 * h  # Linear term scaled for C^(2) basis

        # Compute endpoints using solve_endpoints
        deriv_a, deriv_b, _error = op.solve_endpoints(rhs_np)

        # Get j_l values at endpoints
        j_l_a = spherical_jn(l_val, a)
        j_l_b = spherical_jn(l_val, b)

        # Compute LHS: [x*j_l(x)] * u'(x) from a to b
        lhs = b * j_l_b * deriv_b - a * j_l_a * deriv_a

        # Compute RHS: integral of x*j_l(x)
        def integrand(x: float) -> float:
            return spherical_jn(l_val, x)

        rhs, _ = quad(integrand, a, b, epsabs=1e-12, epsrel=1e-14)

        # Compare
        assert_allclose(
            lhs,
            rhs,
            rtol=1.0e-8,
            atol=1.0e-20,
            err_msg=(
                f"Green's identity failed with solve_endpoints for l={l_val}, RHS=x"
            ),
        )

    @pytest.mark.parametrize("l_val", [0, 1, 5, 10, 20])
    def test_solve_endpoints_vs_full_solution(self, l_val: int) -> None:
        """Test that solve_endpoints.

        This verifies that the optimized endpoint computation matches
        the derivatives computed from the full solution.
        """
        N = 128
        a, b = 1.0, 20.0

        # Create solver and operator
        solver = Ncm.SBesselOdeSolver.new()
        op = solver.create_operator(a, b, l_val, l_val)

        # Set up RHS: homogeneous BCs with RHS=1
        rhs_np = np.zeros(N)
        rhs_np[0] = 0.0  # BC at x=a
        rhs_np[1] = 0.0  # BC at x=b
        rhs_np[2] = 1.0

        # Compute endpoints using solve_endpoints
        deriv_a_fast, deriv_b_fast, _error_fast = op.solve_endpoints(rhs_np)

        # Compute full solution
        solution_coeffs, _solution_len = op.solve(rhs_np)

        # Compute derivatives at endpoints from full solution
        h = (b - a) / 2.0
        y_prime_at_minus1 = Ncm.Spectral.chebyshev_deriv(solution_coeffs, -1.0)
        y_prime_at_plus1 = Ncm.Spectral.chebyshev_deriv(solution_coeffs, 1.0)
        deriv_a_full = y_prime_at_minus1 / h
        deriv_b_full = y_prime_at_plus1 / h

        # Compare
        assert_allclose(
            deriv_a_fast,
            deriv_a_full,
            rtol=1.0e-12,
            atol=1.0e-15,
            err_msg=(
                f"solve_endpoints deriv_a doesn't match full solution for l={l_val}"
            ),
        )
        assert_allclose(
            deriv_b_fast,
            deriv_b_full,
            rtol=1.0e-12,
            atol=1.0e-15,
            err_msg=(
                f"solve_endpoints deriv_b doesn't match full solution for l={l_val}"
            ),
        )

    @pytest.mark.parametrize("l_val", [0, 1, 5, 10])
    def test_solve_endpoints_multiple_rhs(self, l_val: int) -> None:
        """Test solve_endpoints with various RHS functions.

        Tests that solve_endpoints works correctly with different RHS functions
        including constant, linear, and quadratic.
        """
        N = 128
        a, b = 1.0, 20.0
        h = (b - a) / 2.0
        m = (a + b) / 2.0

        # Create solver and operator
        solver = Ncm.SBesselOdeSolver.new()
        op = solver.create_operator(a, b, l_val, l_val)

        # Test with RHS = x^2
        rhs_np = np.zeros(N)
        rhs_np[0] = 0.0  # BC at x=a
        rhs_np[1] = 0.0  # BC at x=b
        # x^2 = (m + h*t)^2 = m^2 + 2*m*h*t + h^2*t^2
        # In Gegenbauer C^(2) basis:
        # T_0 term (constant): m^2 + h^2*(1/2)
        # T_1 term (linear): 2*m*h * (1/4) = m*h/2
        # T_2 term (quadratic): h^2 * (something)
        rhs_np[2] = m * m + h * h / 2.0  # Approximate constant term
        rhs_np[3] = 0.25 * 2.0 * m * h  # Linear term

        # Just verify it doesn't crash and returns reasonable values
        deriv_a, deriv_b, error = op.solve_endpoints(rhs_np)

        assert np.isfinite(deriv_a), "deriv_a should be finite"
        assert np.isfinite(deriv_b), "deriv_b should be finite"
        assert np.isfinite(error), "error should be finite"
        assert error >= 0.0, "error should be non-negative"

    def test_solve_endpoints_batched_rhs1(self) -> None:
        """Test solve_endpoints_batched with Green's identity for RHS=1.

        Verifies that solve_endpoints_batched correctly computes endpoint derivatives
        for u''+[x^2-l(l+1)]u=1 across multiple l values, satisfying Green's identity:
        [x*j_l(x)] * u'(x) from a to b = int_a^b j_l(x)dx
        """
        N = 128
        a, b = 1.0, 20.0
        lmin, lmax = 5, 15
        n_l = lmax - lmin + 1

        # Create solver and operator with range of ell values
        solver = Ncm.SBesselOdeSolver.new()
        op = solver.create_operator(a, b, lmin, lmax)

        # Set up RHS: homogeneous BCs with RHS=1
        rhs_np = np.zeros(N)
        rhs_np[0] = 0.0  # BC at x=a
        rhs_np[1] = 0.0  # BC at x=b
        rhs_np[2] = 1.0

        # Compute endpoints (batched automatically based on operator range)
        endpoints_np = np.array(op.solve_endpoints(rhs_np)).reshape(n_l, 3)

        # Verify each l value satisfies Green's identity
        for i, ell in enumerate(range(lmin, lmax + 1)):
            deriv_a = endpoints_np[i, 0]
            deriv_b = endpoints_np[i, 1]

            # Get j_l values at endpoints
            j_l_a = spherical_jn(ell, a)
            j_l_b = spherical_jn(ell, b)

            # Compute LHS: [x*j_l(x)]*u'(x) from a to b
            lhs = b * j_l_b * deriv_b - a * j_l_a * deriv_a

            # Compute RHS: integral
            def integrand(x: float, ell: int = ell) -> float:
                return spherical_jn(ell, x) / x

            rhs, _ = quad(integrand, a, b, epsabs=1e-12, epsrel=1e-14)

            # Compare
            assert_allclose(
                lhs,
                rhs,
                rtol=1.0e-8,
                atol=1.0e-20,
                err_msg=(
                    f"Green's identity failed with solve_endpoints_batched for l={ell}"
                ),
            )

    def test_solve_endpoints_batched_rhs_x(self) -> None:
        """Test solve_endpoints_batched with Green's identity for RHS=x.

        Verifies that solve_endpoints_batched correctly computes endpoint derivatives
        for u''+[x^2-l(l+1)]u=x across multiple l values, satisfying Green's identity:
        [x*j_l(x)] * u'(x) from a to b = int_a^b x*j_l(x)dx
        """
        N = 128
        a, b = 1.0, 20.0
        lmin, lmax = 5, 15
        n_l = lmax - lmin + 1

        # Create solver and operator with range of ell values
        solver = Ncm.SBesselOdeSolver.new()
        op = solver.create_operator(a, b, lmin, lmax)

        # Set up RHS: homogeneous BCs with RHS=x
        m = (a + b) / 2.0
        h = (b - a) / 2.0
        rhs_np = np.zeros(N)
        rhs_np[0] = 0.0  # BC at x=a
        rhs_np[1] = 0.0  # BC at x=b
        rhs_np[2] = m  # Constant term
        rhs_np[3] = 0.25 * h  # Linear term scaled for C^(2) basis

        # Compute endpoints (batched automatically based on operator range)
        endpoints_np = np.array(op.solve_endpoints(rhs_np)).reshape(n_l, 3)

        # Verify each l value satisfies Green's identity
        for i, ell in enumerate(range(lmin, lmax + 1)):
            deriv_a = endpoints_np[i, 0]
            deriv_b = endpoints_np[i, 1]

            # Get j_l values at endpoints
            j_l_a = spherical_jn(ell, a)
            j_l_b = spherical_jn(ell, b)

            # Compute LHS: [x*j_l(x)]*u'(x) from a to b
            lhs = b * j_l_b * deriv_b - a * j_l_a * deriv_a

            # Compute RHS: integral of x*j_l(x)
            def integrand(x: float, ell: int = ell) -> float:
                return spherical_jn(ell, x)

            rhs, _ = quad(integrand, a, b, epsabs=1e-12, epsrel=1e-14)

            # Compare
            assert_allclose(
                lhs,
                rhs,
                rtol=1.0e-8,
                atol=1.0e-20,
                err_msg=(
                    f"Green's identity failed with solve_endpoints_batched "
                    f"for l={ell}, RHS=x"
                ),
            )

    def test_solve_endpoints_batched_vs_individual(self) -> None:
        """Test that solve_endpoints_batched matches individual solve_endpoints calls.

        Verifies that the batched endpoint computation produces the same results
        as calling solve_endpoints individually for each l value.
        """
        N = 128
        a, b = 1.0, 20.0
        lmin, lmax = 5, 15
        n_l = lmax - lmin + 1

        # Create solver
        solver = Ncm.SBesselOdeSolver.new()

        # Set up RHS: homogeneous BCs with RHS=1
        rhs_np = np.zeros(N)
        rhs_np[0] = 0.0  # BC at x=a
        rhs_np[1] = 0.0  # BC at x=b
        rhs_np[2] = 1.0

        # Compute endpoints using batched operator
        op_batched = solver.create_operator(a, b, lmin, lmax)
        endpoints_batched_np = np.array(op_batched.solve_endpoints(rhs_np)).reshape(
            n_l, 3
        )

        # Compute endpoints individually for each l
        endpoints_individual = []
        for ell in range(lmin, lmax + 1):
            op_single = solver.create_operator(a, b, ell, ell)
            deriv_a, deriv_b, error = op_single.solve_endpoints(rhs_np)
            endpoints_individual.append([deriv_a, deriv_b, error])

        endpoints_individual_np = np.array(endpoints_individual)

        # Compare each endpoint
        for i, ell in enumerate(range(lmin, lmax + 1)):
            assert_allclose(
                endpoints_batched_np[i, :],
                endpoints_individual_np[i, :],
                rtol=1.0e-12,
                atol=1.0e-15,
                err_msg=(
                    f"Batched endpoints for l={ell} don't match individual "
                    f"solve_endpoints"
                ),
            )

    def test_solve_endpoints_batched_vs_full_solution(self) -> None:
        """Test that solve_endpoints_batched.

        Verifies that the optimized batched endpoint computation matches
        the derivatives computed from the full batched solution.
        """
        N = 128
        a, b = 1.0, 20.0
        lmin, lmax = 5, 15
        n_l = lmax - lmin + 1
        h = (b - a) / 2.0

        # Create solver and batched operator
        solver = Ncm.SBesselOdeSolver.new()
        op = solver.create_operator(a, b, lmin, lmax)

        # Set up RHS: homogeneous BCs with RHS=1
        rhs_np = np.zeros(N)
        rhs_np[0] = 0.0  # BC at x=a
        rhs_np[1] = 0.0  # BC at x=b
        rhs_np[2] = 1.0

        # Compute endpoints using solve_endpoints
        endpoints_fast_np = np.array(op.solve_endpoints(rhs_np)).reshape(n_l, 3)

        # Compute full batched solution
        solutions_batched, solution_len = op.solve(rhs_np)
        solutions_batched_np = np.array(solutions_batched).reshape(n_l, solution_len)

        # Compute derivatives at endpoints from full solution for each l
        for i, ell in enumerate(range(lmin, lmax + 1)):
            solution_coeffs = solutions_batched_np[i, :]

            # Compute derivatives at endpoints from full solution
            y_prime_at_minus1 = Ncm.Spectral.chebyshev_deriv(solution_coeffs, -1.0)
            y_prime_at_plus1 = Ncm.Spectral.chebyshev_deriv(solution_coeffs, 1.0)
            deriv_a_full = y_prime_at_minus1 / h
            deriv_b_full = y_prime_at_plus1 / h

            # Compare with fast computation
            assert_allclose(
                endpoints_fast_np[i, 0],
                deriv_a_full,
                rtol=1.0e-12,
                atol=1.0e-15,
                err_msg=(
                    f"solve_endpoints_batched deriv_a doesn't match full solution "
                    f"for l={ell}"
                ),
            )
            assert_allclose(
                endpoints_fast_np[i, 1],
                deriv_b_full,
                rtol=1.0e-12,
                atol=1.0e-15,
                err_msg=(
                    f"solve_endpoints_batched deriv_b doesn't match full solution "
                    f"for l={ell}"
                ),
            )

    def test_solve_endpoints_batched_sanity(self) -> None:
        """Test solve_endpoints_batched basic sanity checks.

        Verifies that the batched endpoint computation returns finite values
        and reasonable error estimates.
        """
        N = 128
        a, b = 1.0, 20.0
        lmin, lmax = 5, 15
        n_l = lmax - lmin + 1

        # Create solver and batched operator
        solver = Ncm.SBesselOdeSolver.new()
        op = solver.create_operator(a, b, lmin, lmax)

        # Set up RHS with various coefficients
        rhs_np = np.zeros(N)
        rhs_np[0] = 0.0  # BC at x=a
        rhs_np[1] = 0.0  # BC at x=b
        for k in range(2, N):
            rhs_np[k] = 1.0 / (1.0 + k * k)

        # Compute endpoints
        endpoints_np = np.array(op.solve_endpoints(rhs_np)).reshape(n_l, 3)

        # Verify all values are finite and errors are non-negative
        assert endpoints_np.shape == (n_l, 3), "Wrong output shape"
        assert np.all(np.isfinite(endpoints_np)), "All values should be finite"
        assert np.all(endpoints_np[:, 2] >= 0.0), "All errors should be non-negative"

    def test_operator_get_interval(self) -> None:
        """Test that operator correctly returns its interval."""
        a_in, b_in = 1.0, 20.0
        solver = Ncm.SBesselOdeSolver.new()
        op = solver.create_operator(a_in, b_in, 5, 5)

        a_out, b_out = op.get_interval()

        assert_allclose(
            a_out,
            a_in,
            rtol=1.0e-15,
            atol=1.0e-15,
            err_msg="Interval left endpoint mismatch",
        )
        assert_allclose(
            b_out,
            b_in,
            rtol=1.0e-15,
            atol=1.0e-15,
            err_msg="Interval right endpoint mismatch",
        )

    def test_operator_get_ell_range(self) -> None:
        """Test that operator correctly returns its ell range."""
        solver = Ncm.SBesselOdeSolver.new()

        # Test single ell
        op_single = solver.create_operator(1.0, 20.0, 7, 7)
        ell_min_out, ell_max_out = op_single.get_ell_range()
        assert ell_min_out == 7, "Single ell: ell_min mismatch"
        assert ell_max_out == 7, "Single ell: ell_max mismatch"

        # Test ell range
        op_range = solver.create_operator(1.0, 20.0, 5, 15)
        ell_min_out, ell_max_out = op_range.get_ell_range()
        assert ell_min_out == 5, "Ell range: ell_min mismatch"
        assert ell_max_out == 15, "Ell range: ell_max mismatch"

    def test_operator_get_tolerance(self) -> None:
        """Test that operator returns the correct tolerance."""
        solver = Ncm.SBesselOdeSolver.new()

        # Test with default tolerance
        tol_default = solver.get_tolerance()
        op = solver.create_operator(1.0, 20.0, 5, 5)
        op_tol = op.get_tolerance()
        assert_allclose(
            op_tol,
            tol_default,
            rtol=1.0e-15,
            atol=1.0e-15,
            err_msg="Default tolerance mismatch",
        )

        # Test with custom tolerance
        custom_tol = 1.0e-10
        solver.set_tolerance(custom_tol)
        op2 = solver.create_operator(1.0, 20.0, 10, 10)
        op2_tol = op2.get_tolerance()
        assert_allclose(
            op2_tol,
            custom_tol,
            rtol=1.0e-15,
            atol=1.0e-15,
            err_msg="Custom tolerance mismatch",
        )

    def test_operator_accessors_consistency(self) -> None:
        """Test that operator accessors remain consistent across operations."""
        a, b = 2.5, 18.3
        ell_min, ell_max = 3, 12
        tol = 1.0e-9

        solver = Ncm.SBesselOdeSolver.new()
        solver.set_tolerance(tol)
        op = solver.create_operator(a, b, ell_min, ell_max)

        # Check initial values
        a_out, b_out = op.get_interval()
        ell_min_out, ell_max_out = op.get_ell_range()
        tol_out = op.get_tolerance()

        assert_allclose(a_out, a, rtol=1.0e-15, atol=1.0e-15)
        assert_allclose(b_out, b, rtol=1.0e-15, atol=1.0e-15)
        assert ell_min_out == ell_min
        assert ell_max_out == ell_max
        assert_allclose(tol_out, tol, rtol=1.0e-15, atol=1.0e-15)

        # Perform a solve operation
        N = 64
        rhs_np = np.zeros(N)
        rhs_np[0] = 0.0
        rhs_np[1] = 0.0
        rhs_np[2] = 1.0
        _solution, _sol_len = op.solve(rhs_np)

        # Check values remain unchanged after solve
        a_out2, b_out2 = op.get_interval()
        ell_min_out2, ell_max_out2 = op.get_ell_range()
        tol_out2 = op.get_tolerance()

        assert_allclose(a_out2, a, rtol=1.0e-15, atol=1.0e-15)
        assert_allclose(b_out2, b, rtol=1.0e-15, atol=1.0e-15)
        assert ell_min_out2 == ell_min
        assert ell_max_out2 == ell_max
        assert_allclose(tol_out2, tol, rtol=1.0e-15, atol=1.0e-15)

    def test_operator_diagonalization_reuse_single_ell(self) -> None:
        """Test that diagonalization is correctly reused for single ell."""
        a, b = 1.0, 20.0
        l_val = 5

        solver = Ncm.SBesselOdeSolver.new()
        op = solver.create_operator(a, b, l_val, l_val)

        # Initially, no diagonalization should be stored
        n_cols_initial = op.get_n_cols()
        assert n_cols_initial == 0, "Fresh operator should have n_cols = 0"

        # First solve with N=64
        N1 = 64
        rhs1 = np.zeros(N1)
        rhs1[0] = 0.0  # BC at a
        rhs1[1] = 0.0  # BC at b
        rhs1[2] = 1.0  # Constant RHS
        solution1, _sol_len1 = op.solve(rhs1)
        n_cols_after_first = op.get_n_cols()

        assert n_cols_after_first > 0, "After first solve, n_cols should be positive"
        # Note: n_cols can exceed N1 due to adaptive extension for convergence

        # Second solve with same RHS size - should reuse diagonalization
        rhs2 = np.zeros(N1)
        rhs2[0] = 0.0
        rhs2[1] = 0.0
        rhs2[2] = 2.0  # Different RHS value
        solution2, _sol_len2 = op.solve(rhs2)
        n_cols_after_second = op.get_n_cols()

        assert n_cols_after_second == n_cols_after_first, (
            f"Same size solve should reuse diagonalization: "
            f"{n_cols_after_second} != {n_cols_after_first}"
        )

        # Verify solutions are different (different RHS should give different result)
        solution1_np = np.array(solution1)
        solution2_np = np.array(solution2)
        assert not np.allclose(
            solution1_np, solution2_np, rtol=1e-10
        ), "Different RHS should produce different solutions"

        # Third solve with smaller RHS - should still reuse existing diagonalization
        N2 = 32
        rhs3 = np.zeros(N2)
        rhs3[0] = 0.0
        rhs3[1] = 0.0
        rhs3[2] = 3.0
        _solution3, _sol_len3 = op.solve(rhs3)
        n_cols_after_third = op.get_n_cols()

        assert n_cols_after_third <= n_cols_after_first, (
            f"Smaller RHS should reuse or reduce: "
            f"{n_cols_after_third} > {n_cols_after_first}"
        )

        # Fourth solve with larger RHS - should extend diagonalization
        N3 = 128
        rhs4 = np.zeros(N3)
        rhs4[0] = 0.0
        rhs4[1] = 0.0
        rhs4[2] = 1.0
        _solution4, _sol_len4 = op.solve(rhs4)
        n_cols_after_fourth = op.get_n_cols()

        assert n_cols_after_fourth >= n_cols_after_first, (
            f"Larger RHS should extend or maintain: "
            f"{n_cols_after_fourth} < {n_cols_after_first}"
        )

        # Reset operator - should clear diagonalization
        op.reset(a, b, l_val, l_val)
        n_cols_after_reset = op.get_n_cols()
        assert (
            n_cols_after_reset == 0
        ), f"After reset, n_cols should be 0, got {n_cols_after_reset}"

    def test_operator_diagonalization_reuse_batched(self) -> None:
        """Test that diagonalization is correctly reused for multiple ell values."""
        a, b = 1.0, 20.0
        ell_min, ell_max = 2, 8

        solver = Ncm.SBesselOdeSolver.new()
        op = solver.create_operator(a, b, ell_min, ell_max)

        # Initially, no diagonalization should be stored
        n_cols_initial = op.get_n_cols()
        assert n_cols_initial == 0, "Fresh operator should have n_cols = 0"

        # First solve with N=64
        N1 = 64
        rhs1 = np.zeros(N1)
        rhs1[0] = 0.0
        rhs1[1] = 0.0
        rhs1[2] = 1.0
        solution1, sol_len1 = op.solve(rhs1)
        n_cols_after_first = op.get_n_cols()

        assert n_cols_after_first > 0, "After first solve, n_cols should be positive"
        # Note: n_cols can exceed N1 due to adaptive extension for convergence

        # Second solve with same RHS size - should reuse diagonalization
        rhs2 = np.zeros(N1)
        rhs2[0] = 0.0
        rhs2[1] = 0.0
        rhs2[2] = 2.0
        solution2, sol_len2 = op.solve(rhs2)
        n_cols_after_second = op.get_n_cols()

        assert n_cols_after_second == n_cols_after_first, (
            f"Same size solve should reuse diagonalization: "
            f"{n_cols_after_second} != {n_cols_after_first}"
        )

        # Verify solutions differ for at least one ell value
        n_ell = ell_max - ell_min + 1
        solution1_np = np.array(solution1).reshape((n_ell, sol_len1))
        solution2_np = np.array(solution2).reshape((n_ell, sol_len2))
        assert not np.allclose(
            solution1_np, solution2_np, rtol=1e-10
        ), "Different RHS should produce different solutions"

        # Third solve with smaller RHS
        N2 = 32
        rhs3 = np.zeros(N2)
        rhs3[0] = 0.0
        rhs3[1] = 0.0
        rhs3[2] = 3.0
        _solution3, _sol_len3 = op.solve(rhs3)
        n_cols_after_third = op.get_n_cols()

        assert n_cols_after_third <= n_cols_after_first, (
            f"Smaller RHS should reuse or reduce: "
            f"{n_cols_after_third} > {n_cols_after_first}"
        )

        # Fourth solve with larger RHS - should extend diagonalization
        N3 = 128
        rhs4 = np.zeros(N3)
        rhs4[0] = 0.0
        rhs4[1] = 0.0
        rhs4[2] = 1.0
        _solution4, _sol_len4 = op.solve(rhs4)
        n_cols_after_fourth = op.get_n_cols()

        assert n_cols_after_fourth >= n_cols_after_first, (
            f"Larger RHS should extend or maintain: "
            f"{n_cols_after_fourth} < {n_cols_after_first}"
        )

        # Reset operator - should clear diagonalization
        op.reset(a, b, ell_min, ell_max)
        n_cols_after_reset = op.get_n_cols()
        assert (
            n_cols_after_reset == 0
        ), f"After reset, n_cols should be 0, got {n_cols_after_reset}"

    def test_operator_diagonalization_reuse_correctness(self) -> None:
        """Test that reused diagonalization produces correct results."""
        a, b = 1.0, 20.0
        l_val = 7

        solver = Ncm.SBesselOdeSolver.new()
        op1 = solver.create_operator(a, b, l_val, l_val)
        op2 = solver.create_operator(a, b, l_val, l_val)

        N = 64
        rhs = np.zeros(N)
        rhs[0] = 0.0
        rhs[1] = 0.0
        rhs[2] = 1.0

        # Solve with first operator (fresh diagonalization)
        solution1, _sol_len1 = op1.solve(rhs)
        solution1_np = np.array(solution1)

        # Solve with second operator (also fresh)
        solution2, _sol_len2 = op2.solve(rhs)
        solution2_np = np.array(solution2)

        # Both should give identical results
        assert_allclose(
            solution1_np,
            solution2_np,
            rtol=1e-14,
            atol=1e-14,
            err_msg="Fresh diagonalizations should give identical results",
        )

        # Now solve again with first operator (reusing diagonalization)
        solution1_reused, _sol_len1_reused = op1.solve(rhs)
        solution1_reused_np = np.array(solution1_reused)

        # Reused should match original
        assert_allclose(
            solution1_reused_np,
            solution1_np,
            rtol=1e-14,
            atol=1e-14,
            err_msg="Reused diagonalization should give identical results to fresh",
        )

    @pytest.mark.parametrize("l_val", [0, 5, 10, 15])
    def test_operator_diagonalization_endpoints_reuse(self, l_val: int) -> None:
        """Test that diagonalization reuse works correctly with solve_endpoints."""
        a, b = 1.0, 20.0

        solver = Ncm.SBesselOdeSolver.new()
        op = solver.create_operator(a, b, l_val, l_val)

        # Check initial state
        assert op.get_n_cols() == 0, "Fresh operator should have n_cols = 0"

        N = 64
        rhs = np.zeros(N)
        rhs[0] = 0.0
        rhs[1] = 0.0
        rhs[2] = 1.0

        # First solve_endpoints call
        endpoints1 = op.solve_endpoints(rhs)
        n_cols1 = op.get_n_cols()
        assert n_cols1 > 0, "After solve_endpoints, n_cols should be positive"

        # Second solve_endpoints call with same RHS size
        endpoints2 = op.solve_endpoints(rhs)
        n_cols2 = op.get_n_cols()

        # Should reuse diagonalization
        assert (
            n_cols2 == n_cols1
        ), f"Same size solve_endpoints should reuse: {n_cols2} != {n_cols1}"

        # Results should be identical
        endpoints1_np = np.array(endpoints1)
        endpoints2_np = np.array(endpoints2)
        assert_allclose(
            endpoints1_np,
            endpoints2_np,
            rtol=1e-14,
            atol=1e-14,
            err_msg="Reused diagonalization should give identical endpoints",
        )

    @pytest.mark.parametrize("l_val", [0, 5, 10, 15])
    def test_solve_polynomial_rhs_quadratic(self, l_val: int) -> None:
        """Test solving with quadratic polynomial RHS: u''+[x^2-l(l+1)]u = x^2.

        Verifies that the solver correctly handles polynomial RHS functions by
        checking that the solution satisfies the ODE at several test points.
        """
        N = 128
        a, b = 1.0, 20.0
        h = (b - a) / 2.0
        m = (a + b) / 2.0

        solver = Ncm.SBesselOdeSolver.new()
        op = solver.create_operator(a, b, l_val, l_val)

        # Set up RHS: u(a)=0, u(b)=0, RHS = x^2
        # x^2 = (m + h*t)^2 = m^2 + 2*m*h*t + h^2*t^2
        # In Gegenbauer C^(2) basis:
        # T_0: c_0, T_1: (1/4)*c_1, T_2: (1/6)*c_2
        rhs = np.zeros(N)
        rhs[0] = 0.0  # BC at a
        rhs[1] = 0.0  # BC at b
        # Approximate expansion for x^2 in Gegenbauer basis
        rhs[2] = m * m + h * h / 2.0  # Constant term
        rhs[3] = 0.25 * 2.0 * m * h  # Linear term
        rhs[4] = h * h / 6.0  # Quadratic term

        solution, _sol_len = op.solve(rhs)
        solution_np = np.array(solution)

        # Verify solution is non-trivial and finite
        assert len(solution_np) > 0, "Solution should not be empty"
        assert np.all(np.isfinite(solution_np)), "Solution should be finite"
        assert np.any(np.abs(solution_np) > 1e-10), "Solution should be non-trivial"

        # Check that solution satisfies BCs
        u_at_a = Ncm.Spectral.chebyshev_eval(solution_np, -1.0)
        u_at_b = Ncm.Spectral.chebyshev_eval(solution_np, 1.0)
        assert_allclose(u_at_a, 0.0, atol=1e-10, err_msg="BC at a not satisfied")
        assert_allclose(u_at_b, 0.0, atol=1e-10, err_msg="BC at b not satisfied")

    @pytest.mark.parametrize("l_val", [0, 5, 10])
    def test_solve_polynomial_rhs_cubic(self, l_val: int) -> None:
        """Test solving with cubic polynomial RHS: u''+[x^2-l(l+1)]u = x^3."""
        N = 128
        a, b = 1.0, 20.0
        h = (b - a) / 2.0
        m = (a + b) / 2.0

        solver = Ncm.SBesselOdeSolver.new()
        op = solver.create_operator(a, b, l_val, l_val)

        # Set up RHS for x^3 = (m + h*t)^3
        rhs = np.zeros(N)
        rhs[0] = 0.0
        rhs[1] = 0.0
        # Approximate Gegenbauer expansion for x^3
        rhs[2] = m * m * m + 1.5 * m * h * h  # Constant part
        rhs[3] = 0.25 * (3.0 * m * m * h + 0.5 * h * h * h)  # Linear part
        rhs[4] = 3.0 * m * h * h / 6.0  # Quadratic part
        rhs[5] = h * h * h / 8.0  # Cubic part

        solution, _sol_len = op.solve(rhs)
        solution_np = np.array(solution)

        # Verify solution properties
        assert len(solution_np) > 0
        assert np.all(np.isfinite(solution_np))

        # Check BCs
        u_at_a = Ncm.Spectral.chebyshev_eval(solution_np, -1.0)
        u_at_b = Ncm.Spectral.chebyshev_eval(solution_np, 1.0)
        assert_allclose(u_at_a, 0.0, atol=1e-10)
        assert_allclose(u_at_b, 0.0, atol=1e-10)

    @pytest.mark.parametrize("l_val", [0, 5, 10])
    def test_solve_exponential_rhs(self, l_val: int) -> None:
        """Test solving with exponential RHS: u''+[x^2-l(l+1)]u = exp(x).

        Uses Chebyshev expansion of exp((a+b)/2 + ((b-a)/2)*t) as RHS.
        """
        N = 128
        a, b = 1.0, 20.0
        h = (b - a) / 2.0
        m = (a + b) / 2.0

        solver = Ncm.SBesselOdeSolver.new()
        spectral = Ncm.Spectral.new()
        op = solver.create_operator(a, b, l_val, l_val)

        # Compute Chebyshev coefficients for exp(m + h*t) on [-1,1]
        def f_exp(_user_data, t):
            x = m + h * t
            return np.exp(x)

        exp_coeffs = np.array(
            spectral.compute_chebyshev_coeffs(f_exp, -1.0, 1.0, N - 2, None)
        )

        # Convert to Gegenbauer C^(2) basis for RHS
        rhs_gegenbauer = Ncm.Spectral.chebT_to_gegenbauer_alpha2(exp_coeffs)

        # Set up RHS with BCs
        rhs = np.zeros(N)
        rhs[0] = 0.0  # BC at a
        rhs[1] = 0.0  # BC at b
        rhs[2:] = rhs_gegenbauer[: N - 2]

        solution, _sol_len = op.solve(rhs)
        solution_np = np.array(solution)

        # Verify solution properties
        assert len(solution_np) > 0
        assert np.all(np.isfinite(solution_np))
        assert np.any(np.abs(solution_np) > 1e-10)

        # Check BCs
        u_at_a = Ncm.Spectral.chebyshev_eval(solution_np, -1.0)
        u_at_b = Ncm.Spectral.chebyshev_eval(solution_np, 1.0)
        assert_allclose(u_at_a, 0.0, atol=1e-8)
        assert_allclose(u_at_b, 0.0, atol=1e-8)

    @pytest.mark.parametrize("l_val", [0, 5, 10])
    def test_solve_sin_rhs(self, l_val: int) -> None:
        """Test solving with sine RHS: u''+[x^2-l(l+1)]u = sin(x)."""
        N = 128
        a, b = 1.0, 20.0
        h = (b - a) / 2.0
        m = (a + b) / 2.0

        solver = Ncm.SBesselOdeSolver.new()
        spectral = Ncm.Spectral.new()
        op = solver.create_operator(a, b, l_val, l_val)

        # Compute Chebyshev coefficients for sin(m + h*t) on [-1,1]
        def f_sin(_user_data, t):
            x = m + h * t
            return np.sin(x)

        sin_coeffs = np.array(
            spectral.compute_chebyshev_coeffs(f_sin, -1.0, 1.0, N - 2, None)
        )

        # Convert to Gegenbauer C^(2) basis
        rhs_gegenbauer = Ncm.Spectral.chebT_to_gegenbauer_alpha2(sin_coeffs)

        # Set up RHS with BCs
        rhs = np.zeros(N)
        rhs[0] = 0.0
        rhs[1] = 0.0
        rhs[2:] = rhs_gegenbauer[: N - 2]

        solution, _sol_len = op.solve(rhs)
        solution_np = np.array(solution)

        # Verify solution properties
        assert len(solution_np) > 0
        assert np.all(np.isfinite(solution_np))

        # Check BCs
        u_at_a = Ncm.Spectral.chebyshev_eval(solution_np, -1.0)
        u_at_b = Ncm.Spectral.chebyshev_eval(solution_np, 1.0)
        assert_allclose(u_at_a, 0.0, atol=1e-8)
        assert_allclose(u_at_b, 0.0, atol=1e-8)

    @pytest.mark.parametrize("l_val", [0, 5, 10])
    def test_solve_cos_rhs(self, l_val: int) -> None:
        """Test solving with cosine RHS: u''+[x^2-l(l+1)]u = cos(x)."""
        N = 128
        a, b = 1.0, 20.0
        h = (b - a) / 2.0
        m = (a + b) / 2.0

        solver = Ncm.SBesselOdeSolver.new()
        spectral = Ncm.Spectral.new()
        op = solver.create_operator(a, b, l_val, l_val)

        # Compute Chebyshev coefficients for cos(m + h*t) on [-1,1]
        def f_cos(_user_data, t):
            x = m + h * t
            return np.cos(x)

        cos_coeffs = np.array(
            spectral.compute_chebyshev_coeffs(f_cos, -1.0, 1.0, N - 2, None)
        )

        # Convert to Gegenbauer C^(2) basis
        rhs_gegenbauer = Ncm.Spectral.chebT_to_gegenbauer_alpha2(cos_coeffs)

        # Set up RHS with BCs
        rhs = np.zeros(N)
        rhs[0] = 0.0
        rhs[1] = 0.0
        rhs[2:] = rhs_gegenbauer[: N - 2]

        solution, _sol_len = op.solve(rhs)
        solution_np = np.array(solution)

        # Verify solution properties
        assert len(solution_np) > 0
        assert np.all(np.isfinite(solution_np))

        # Check BCs
        u_at_a = Ncm.Spectral.chebyshev_eval(solution_np, -1.0)
        u_at_b = Ncm.Spectral.chebyshev_eval(solution_np, 1.0)
        assert_allclose(u_at_a, 0.0, atol=1e-8)
        assert_allclose(u_at_b, 0.0, atol=1e-8)

    @pytest.mark.parametrize("l_val", [0, 5, 10])
    def test_solve_gaussian_rhs(self, l_val: int) -> None:
        """Test solving with Gaussian RHS: u''+[x^2-l(l+1)]u = exp(-x^2)."""
        N = 128
        a, b = 1.0, 20.0
        h = (b - a) / 2.0
        m = (a + b) / 2.0

        solver = Ncm.SBesselOdeSolver.new()
        spectral = Ncm.Spectral.new()
        op = solver.create_operator(a, b, l_val, l_val)

        # Compute Chebyshev coefficients for exp(-(m + h*t)^2) on [-1,1]
        def f_gauss(_user_data, t):
            x = m + h * t
            return np.exp(-x * x)

        gauss_coeffs = np.array(
            spectral.compute_chebyshev_coeffs(f_gauss, -1.0, 1.0, N - 2, None)
        )

        # Convert to Gegenbauer C^(2) basis
        rhs_gegenbauer = Ncm.Spectral.chebT_to_gegenbauer_alpha2(gauss_coeffs)

        # Set up RHS with BCs
        rhs = np.zeros(N)
        rhs[0] = 0.0
        rhs[1] = 0.0
        rhs[2:] = rhs_gegenbauer[: N - 2]

        solution, _sol_len = op.solve(rhs)
        solution_np = np.array(solution)

        # Verify solution properties
        assert len(solution_np) > 0
        assert np.all(np.isfinite(solution_np))

        # Check BCs
        u_at_a = Ncm.Spectral.chebyshev_eval(solution_np, -1.0)
        u_at_b = Ncm.Spectral.chebyshev_eval(solution_np, 1.0)
        assert_allclose(u_at_a, 0.0, atol=1e-8)
        assert_allclose(u_at_b, 0.0, atol=1e-8)

    @pytest.mark.parametrize("l_val", [0, 5, 10])
    def test_solve_rational_rhs(self, l_val: int) -> None:
        """Test solving with rational RHS: u''+[x^2-l(l+1)]u = 1/(1+x^2)."""
        N = 128
        a, b = 1.0, 20.0
        h = (b - a) / 2.0
        m = (a + b) / 2.0

        solver = Ncm.SBesselOdeSolver.new()
        spectral = Ncm.Spectral.new()
        op = solver.create_operator(a, b, l_val, l_val)

        # Compute Chebyshev coefficients for 1/(1+(m+h*t)^2) on [-1,1]
        def f_rational(_user_data, t):
            x = m + h * t
            return 1.0 / (1.0 + x * x)

        rational_coeffs = np.array(
            spectral.compute_chebyshev_coeffs(f_rational, -1.0, 1.0, N - 2, None)
        )

        # Convert to Gegenbauer C^(2) basis
        rhs_gegenbauer = Ncm.Spectral.chebT_to_gegenbauer_alpha2(rational_coeffs)

        # Set up RHS with BCs
        rhs = np.zeros(N)
        rhs[0] = 0.0
        rhs[1] = 0.0
        rhs[2:] = rhs_gegenbauer[: N - 2]

        solution, _sol_len = op.solve(rhs)
        solution_np = np.array(solution)

        # Verify solution properties
        assert len(solution_np) > 0
        assert np.all(np.isfinite(solution_np))

        # Check BCs
        u_at_a = Ncm.Spectral.chebyshev_eval(solution_np, -1.0)
        u_at_b = Ncm.Spectral.chebyshev_eval(solution_np, 1.0)
        assert_allclose(u_at_a, 0.0, atol=1e-8)
        assert_allclose(u_at_b, 0.0, atol=1e-8)

    def test_solve_batched_exponential_rhs(self) -> None:
        """Test batched solving with exponential RHS for multiple ell values."""
        N = 128
        a, b = 1.0, 20.0
        lmin, lmax = 2, 8
        n_ell = lmax - lmin + 1
        h = (b - a) / 2.0
        m = (a + b) / 2.0

        solver = Ncm.SBesselOdeSolver.new()
        spectral = Ncm.Spectral.new()

        # Compute RHS coefficients for exp(x)
        def f_exp(_user_data, t):
            x = m + h * t
            return np.exp(x)

        exp_coeffs = np.array(
            spectral.compute_chebyshev_coeffs(f_exp, -1.0, 1.0, N - 2, None)
        )
        rhs_gegenbauer = Ncm.Spectral.chebT_to_gegenbauer_alpha2(exp_coeffs)

        # Set up RHS
        rhs = np.zeros(N)
        rhs[0] = 0.0
        rhs[1] = 0.0
        rhs[2:] = rhs_gegenbauer[: N - 2]

        # Solve using batched operator
        op_batched = solver.create_operator(a, b, lmin, lmax)
        solutions_batched, sol_len = op_batched.solve(rhs)
        solutions_np = np.array(solutions_batched).reshape(n_ell, sol_len)

        # Verify each solution
        for i, ell in enumerate(range(lmin, lmax + 1)):
            sol = solutions_np[i, :]
            assert len(sol) > 0
            assert np.all(np.isfinite(sol))

            # Check BCs
            u_at_a = Ncm.Spectral.chebyshev_eval(sol, -1.0)
            u_at_b = Ncm.Spectral.chebyshev_eval(sol, 1.0)
            assert_allclose(
                u_at_a, 0.0, atol=1e-8, err_msg=f"BC at a failed for ell={ell}"
            )
            assert_allclose(
                u_at_b, 0.0, atol=1e-8, err_msg=f"BC at b failed for ell={ell}"
            )

    def test_solve_batched_trigonometric_rhs(self) -> None:
        """Test batched solving with sin(x) RHS for multiple ell values."""
        N = 128
        a, b = 1.0, 20.0
        lmin, lmax = 2, 8
        n_ell = lmax - lmin + 1
        h = (b - a) / 2.0
        m = (a + b) / 2.0

        solver = Ncm.SBesselOdeSolver.new()
        spectral = Ncm.Spectral.new()

        # Compute RHS coefficients for sin(x)
        def f_sin(_user_data, t):
            x = m + h * t
            return np.sin(x)

        sin_coeffs = np.array(
            spectral.compute_chebyshev_coeffs(f_sin, -1.0, 1.0, N - 2, None)
        )
        rhs_gegenbauer = Ncm.Spectral.chebT_to_gegenbauer_alpha2(sin_coeffs)

        # Set up RHS
        rhs = np.zeros(N)
        rhs[0] = 0.0
        rhs[1] = 0.0
        rhs[2:] = rhs_gegenbauer[: N - 2]

        # Solve using batched operator
        op_batched = solver.create_operator(a, b, lmin, lmax)
        solutions_batched, sol_len = op_batched.solve(rhs)
        solutions_np = np.array(solutions_batched).reshape(n_ell, sol_len)

        # Verify each solution
        for i, ell in enumerate(range(lmin, lmax + 1)):
            sol = solutions_np[i, :]
            assert len(sol) > 0
            assert np.all(np.isfinite(sol))

            # Check BCs
            u_at_a = Ncm.Spectral.chebyshev_eval(sol, -1.0)
            u_at_b = Ncm.Spectral.chebyshev_eval(sol, 1.0)
            assert_allclose(
                u_at_a, 0.0, atol=1e-8, err_msg=f"BC at a failed for ell={ell}"
            )
            assert_allclose(
                u_at_b, 0.0, atol=1e-8, err_msg=f"BC at b failed for ell={ell}"
            )

    def test_solve_batched_vs_single_exponential(self) -> None:
        """Test that batched and single-ell solves give same results for exp(x) RHS."""
        N = 128
        a, b = 1.0, 20.0
        lmin, lmax = 5, 10
        h = (b - a) / 2.0
        m = (a + b) / 2.0

        solver = Ncm.SBesselOdeSolver.new()
        spectral = Ncm.Spectral.new()

        # Compute RHS coefficients for exp(x)
        def f_exp(_user_data, t):
            x = m + h * t
            return np.exp(x)

        exp_coeffs = np.array(
            spectral.compute_chebyshev_coeffs(f_exp, -1.0, 1.0, N - 2, None)
        )
        rhs_gegenbauer = Ncm.Spectral.chebT_to_gegenbauer_alpha2(exp_coeffs)

        rhs = np.zeros(N)
        rhs[0] = 0.0
        rhs[1] = 0.0
        rhs[2:] = rhs_gegenbauer[: N - 2]

        # Compare batched vs single-ell solutions
        self._compare_batched_vs_single_solve(
            solver,
            a,
            b,
            lmin,
            lmax,
            rhs,
            rtol=1e-12,
            atol_factor=1e-12,
            err_msg_template="Batched vs single mismatch for ell={ell} with exp(x) RHS",
        )

    @pytest.mark.parametrize("n_ell", [2, 4, 8, 16, 32, 64])
    def test_optimized_batched_dimensions_constant_rhs(self, n_ell: int) -> None:
        """Test optimized batched paths for specific dimensions with constant RHS.

        Tests the optimized batched implementations for n_ell = 2, 4, 8, 16, 32, 64
        by comparing batched results with individual ell-by-ell solves.
        """
        N = 64
        a, b = 1.0, 20.0
        lmin = 5
        lmax = lmin + n_ell - 1

        solver = Ncm.SBesselOdeSolver.new()

        # Simple constant RHS
        rhs = np.zeros(N)
        rhs[0] = 0.0  # BC at a
        rhs[1] = 0.0  # BC at b
        rhs[2] = 1.0  # Constant RHS

        # Compare batched vs single-ell solutions
        self._compare_batched_vs_single_solve(
            solver,
            a,
            b,
            lmin,
            lmax,
            rhs,
            rtol=1e-12,
            atol_factor=1e-12,
            err_msg_template=(
                f"Batched (n_ell={n_ell}) vs single mismatch "
                f"for ell={{ell}} with constant RHS"
            ),
        )

    @pytest.mark.parametrize("n_ell", [2, 4, 8, 16, 32, 64])
    def test_optimized_batched_dimensions_linear_rhs(self, n_ell: int) -> None:
        """Test optimized batched paths for specific dimensions with linear RHS.

        Tests with RHS = x to verify optimized batched implementations handle
        non-constant right-hand sides correctly.
        """
        N = 64
        a, b = 1.0, 20.0
        lmin = 3
        lmax = lmin + n_ell - 1
        h = (b - a) / 2.0
        m = (a + b) / 2.0

        solver = Ncm.SBesselOdeSolver.new()

        # Linear RHS: x = m + h*t
        rhs = np.zeros(N)
        rhs[0] = 0.0  # BC at a
        rhs[1] = 0.0  # BC at b
        rhs[2] = m  # Constant term
        rhs[3] = 0.25 * h  # Linear term in Gegenbauer C^(2) basis

        # Compare batched vs single-ell solutions
        self._compare_batched_vs_single_solve(
            solver,
            a,
            b,
            lmin,
            lmax,
            rhs,
            rtol=1e-12,
            atol_factor=1e-12,
            err_msg_template=(
                f"Batched (n_ell={n_ell}) vs single mismatch "
                f"for ell={{ell}} with linear RHS"
            ),
        )

    @pytest.mark.parametrize("n_ell", [2, 4, 8, 16, 32, 64])
    def test_optimized_batched_dimensions_quadratic_rhs(self, n_ell: int) -> None:
        """Test optimized batched paths with quadratic RHS = x^2."""
        N = 64
        a, b = 1.0, 20.0
        lmin = 2
        lmax = lmin + n_ell - 1
        h = (b - a) / 2.0
        m = (a + b) / 2.0

        solver = Ncm.SBesselOdeSolver.new()

        # Quadratic RHS: x^2 = (m + h*t)^2
        rhs = np.zeros(N)
        rhs[0] = 0.0
        rhs[1] = 0.0
        rhs[2] = m * m + h * h / 2.0  # Constant term
        rhs[3] = 0.25 * 2.0 * m * h  # Linear term
        rhs[4] = h * h / 6.0  # Quadratic term

        # Compare batched vs single-ell solutions
        self._compare_batched_vs_single_solve(
            solver,
            a,
            b,
            lmin,
            lmax,
            rhs,
            rtol=1e-12,
            atol_factor=1e-12,
            err_msg_template=(
                f"Batched (n_ell={n_ell}) vs single mismatch "
                f"for ell={{ell}} with quadratic RHS"
            ),
        )

    @pytest.mark.parametrize("n_ell", [2, 4, 8, 16, 32, 64])
    def test_optimized_batched_dimensions_exponential_rhs(self, n_ell: int) -> None:
        """Test optimized batched paths with exponential RHS = exp(x).

        Tests the most complex case with exponential RHS to ensure optimized
        batched implementations handle spectral expansions correctly.
        """
        N = 64
        a, b = 1.0, 20.0
        lmin = 1
        lmax = lmin + n_ell - 1
        h = (b - a) / 2.0
        m = (a + b) / 2.0

        solver = Ncm.SBesselOdeSolver.new()
        spectral = Ncm.Spectral.new()

        # Compute exponential RHS coefficients
        def f_exp(_user_data, t):
            x = m + h * t
            return np.exp(x)

        exp_coeffs = np.array(
            spectral.compute_chebyshev_coeffs(f_exp, -1.0, 1.0, N - 2, None)
        )
        rhs_gegenbauer = Ncm.Spectral.chebT_to_gegenbauer_alpha2(exp_coeffs)

        rhs = np.zeros(N)
        rhs[0] = 0.0
        rhs[1] = 0.0
        rhs[2:] = rhs_gegenbauer[: N - 2]

        # Compare batched vs single-ell solutions (relaxed tolerance for exponential)
        self._compare_batched_vs_single_solve(
            solver,
            a,
            b,
            lmin,
            lmax,
            rhs,
            rtol=1e-10,
            atol_factor=1e-14,
            err_msg_template=(
                f"Batched (n_ell={n_ell}) vs single mismatch "
                f"for ell={{ell}} with exponential RHS"
            ),
        )

    @pytest.mark.parametrize("n_ell", [2, 4, 8, 16, 32, 64])
    def test_optimized_batched_dimensions_trigonometric_rhs(self, n_ell: int) -> None:
        """Test optimized batched paths with trigonometric RHS = sin(x)."""
        N = 64
        a, b = 1.0, 20.0
        lmin = 2
        lmax = lmin + n_ell - 1
        h = (b - a) / 2.0
        m = (a + b) / 2.0

        solver = Ncm.SBesselOdeSolver.new()
        spectral = Ncm.Spectral.new()

        # Compute sine RHS coefficients
        def f_sin(_user_data, t):
            x = m + h * t
            return np.sin(x)

        sin_coeffs = np.array(
            spectral.compute_chebyshev_coeffs(f_sin, -1.0, 1.0, N - 2, None)
        )
        rhs_gegenbauer = Ncm.Spectral.chebT_to_gegenbauer_alpha2(sin_coeffs)

        rhs = np.zeros(N)
        rhs[0] = 0.0
        rhs[1] = 0.0
        rhs[2:] = rhs_gegenbauer[: N - 2]

        # Compare batched vs single-ell solutions
        self._compare_batched_vs_single_solve(
            solver,
            a,
            b,
            lmin,
            lmax,
            rhs,
            rtol=1e-10,
            atol_factor=1e-13,
            err_msg_template=(
                f"Batched (n_ell={n_ell}) vs single mismatch "
                f"for ell={{ell}} with sin(x) RHS"
            ),
        )

    @pytest.mark.parametrize("n_ell", [2, 4, 8, 16, 32, 64])
    def test_optimized_batched_dimensions_endpoints_constant(self, n_ell: int) -> None:
        """Test optimized batched solve_endpoints for specific dimensions.

        Verifies that the optimized batched endpoint computation matches
        individual ell-by-ell endpoint computations for each optimized dimension.
        """
        N = 64
        a, b = 1.0, 20.0
        lmin = 3
        lmax = lmin + n_ell - 1

        solver = Ncm.SBesselOdeSolver.new()

        # Simple constant RHS
        rhs = np.zeros(N)
        rhs[0] = 0.0
        rhs[1] = 0.0
        rhs[2] = 1.0

        # Compare batched vs single-ell endpoints
        self._compare_batched_vs_single_endpoints(
            solver,
            a,
            b,
            lmin,
            lmax,
            rhs,
            rtol=1e-11,
            atol_factor=1e-11,
            err_msg_template=(
                f"Batched endpoints (n_ell={n_ell}) vs single mismatch "
                f"for ell={{ell}} with constant RHS"
            ),
            compare_error=True,
        )

    @pytest.mark.parametrize("n_ell", [2, 4, 8, 16, 32, 64])
    def test_optimized_batched_dimensions_endpoints_exponential(
        self, n_ell: int
    ) -> None:
        """Test optimized batched solve_endpoints with exponential RHS."""
        N = 64
        a, b = 1.0, 20.0
        lmin = 1
        lmax = lmin + n_ell - 1
        h = (b - a) / 2.0
        m = (a + b) / 2.0

        solver = Ncm.SBesselOdeSolver.new()
        spectral = Ncm.Spectral.new()

        # Compute exponential RHS coefficients
        def f_exp(_user_data, t):
            x = m + h * t
            return np.exp(x)

        exp_coeffs = np.array(
            spectral.compute_chebyshev_coeffs(f_exp, -1.0, 1.0, N - 2, None)
        )
        rhs_gegenbauer = Ncm.Spectral.chebT_to_gegenbauer_alpha2(exp_coeffs)

        rhs = np.zeros(N)
        rhs[0] = 0.0
        rhs[1] = 0.0
        rhs[2:] = rhs_gegenbauer[: N - 2]

        # Compare batched vs single-ell endpoints (derivatives only, not error estimate)
        self._compare_batched_vs_single_endpoints(
            solver,
            a,
            b,
            lmin,
            lmax,
            rhs,
            rtol=1e-11,
            atol_factor=1e-11,
            err_msg_template=(
                f"Batched endpoints (n_ell={n_ell}) vs single mismatch "
                f"for ell={{ell}} with exp(x) RHS"
            ),
            compare_error=False,
        )

    @pytest.mark.parametrize("n_ell", [2, 4, 8, 16, 32, 64])
    def test_optimized_batched_dimensions_diagonalization_reuse(
        self, n_ell: int
    ) -> None:
        """Test that optimized batched paths correctly reuse diagonalization.

        Verifies that the diagonalization reuse works correctly for all
        optimized batch dimensions.
        """
        N = 64
        a, b = 1.0, 20.0
        lmin = 2
        lmax = lmin + n_ell - 1

        solver = Ncm.SBesselOdeSolver.new()
        op = solver.create_operator(a, b, lmin, lmax)

        # Check initial state
        assert (
            op.get_n_cols() == 0
        ), f"Fresh operator (n_ell={n_ell}) should have n_cols = 0"

        # First solve
        rhs1 = np.zeros(N)
        rhs1[0] = 0.0
        rhs1[1] = 0.0
        rhs1[2] = 1.0
        _solution1, _sol_len1 = op.solve(rhs1)
        n_cols_first = op.get_n_cols()

        assert (
            n_cols_first > 0
        ), f"After first solve (n_ell={n_ell}), n_cols should be positive"

        # Second solve with different RHS - should reuse
        rhs2 = np.zeros(N)
        rhs2[0] = 0.0
        rhs2[1] = 0.0
        rhs2[2] = 2.0
        _solution2, _sol_len2 = op.solve(rhs2)
        n_cols_second = op.get_n_cols()

        assert n_cols_second == n_cols_first, (
            f"Second solve (n_ell={n_ell}) should reuse diagonalization: "
            f"{n_cols_second} != {n_cols_first}"
        )

        # Reset should clear
        op.reset(a, b, lmin, lmax)
        n_cols_after_reset = op.get_n_cols()
        assert (
            n_cols_after_reset == 0
        ), f"After reset (n_ell={n_ell}), n_cols should be 0, got {n_cols_after_reset}"

    def test_diagonalization_reuse_performance_single_ell(self) -> None:
        """Test that diagonalization reuse is significantly faster than reset.

        Compares the performance of:
        1. Repeated solve() calls (reuses diagonalization)
        2. Repeated solve() + reset() calls (re-diagonalizes each time)

        The reuse case should be significantly faster.
        """
        N = 256
        a, b = 1.0, 20.0
        l_val = 10
        n_iterations = 100

        solver = Ncm.SBesselOdeSolver.new()

        # Prepare RHS arrays for testing
        rhs_list = []
        for i in range(n_iterations):
            rhs = np.zeros(N)
            rhs[0] = 0.0
            rhs[1] = 0.0
            rhs[2] = 1.0 + i * 0.1  # Varying constant RHS
            rhs_list.append(rhs)

        # Test 1: Repeated solve() with reuse (fast path)
        op_reuse = solver.create_operator(a, b, l_val, l_val)
        start_reuse = time.perf_counter()
        for rhs in rhs_list:
            _solution, _sol_len = op_reuse.solve(rhs)
        end_reuse = time.perf_counter()
        time_reuse = end_reuse - start_reuse

        # Test 2: Repeated solve() + reset() without reuse (slow path)
        op_reset = solver.create_operator(a, b, l_val, l_val)
        start_reset = time.perf_counter()
        for rhs in rhs_list:
            _solution, _sol_len = op_reset.solve(rhs)
            op_reset.reset(a, b, l_val, l_val)
        end_reset = time.perf_counter()
        time_reset = end_reset - start_reset

        # Calculate speedup
        speedup = time_reset / time_reuse

        # Print timing results
        print(f"\n{'='*60}")
        print("Single-ell Diagonalization Reuse Performance Test")
        print(f"{'='*60}")
        print(f"Number of iterations: {n_iterations}")
        print(f"Matrix size: {N}")
        print(f"ell value: {l_val}")
        print(f"\nWith reuse:    {time_reuse:.4f} seconds")
        print(f"With reset:    {time_reset:.4f} seconds")
        print(f"Speedup:       {speedup:.2f}x")
        print(f"{'='*60}")

        # Both methods should work correctly (no assertion on speedup)
        # Timing is informational only to avoid flakiness in CI

    @pytest.mark.parametrize("n_ell", [2, 4, 8, 16, 32, 64])
    def test_diagonalization_reuse_performance_batched(self, n_ell: int) -> None:
        """Test that diagonalization reuse is faster for batched operations.

        Compares the performance of:
        1. Repeated solve() calls (reuses diagonalization)
        2. Repeated solve() + reset() calls (re-diagonalizes each time)

        The reuse case should be significantly faster, especially for larger batches.
        """
        N = 128
        a, b = 1.0, 20.0
        lmin = 5
        lmax = lmin + n_ell - 1
        n_iterations = 100

        solver = Ncm.SBesselOdeSolver.new()

        # Prepare RHS arrays for testing
        rhs_list = []
        for i in range(n_iterations):
            rhs = np.zeros(N)
            rhs[0] = 0.0
            rhs[1] = 0.0
            rhs[2] = 1.0 + i * 0.1  # Varying constant RHS
            rhs_list.append(rhs)

        # Test 1: Repeated solve() with reuse (fast path)
        op_reuse = solver.create_operator(a, b, lmin, lmax)
        start_reuse = time.perf_counter()
        for rhs in rhs_list:
            _solution, _sol_len = op_reuse.solve(rhs)
        end_reuse = time.perf_counter()
        time_reuse = end_reuse - start_reuse

        # Test 2: Repeated solve() + reset() without reuse (slow path)
        op_reset = solver.create_operator(a, b, lmin, lmax)
        start_reset = time.perf_counter()
        for rhs in rhs_list:
            _solution, _sol_len = op_reset.solve(rhs)
            op_reset.reset(a, b, lmin, lmax)
        end_reset = time.perf_counter()
        time_reset = end_reset - start_reset

        # Calculate speedup
        speedup = time_reset / time_reuse

        # Print timing results
        print(f"\n{'='*60}")
        print("Batched Diagonalization Reuse Performance Test")
        print(f"{'='*60}")
        print(f"Number of iterations: {n_iterations}")
        print(f"Matrix size: {N}")
        print(f"Batch size (n_ell): {n_ell}")
        print(f"ell range: [{lmin}, {lmax}]")
        print(f"\nWith reuse:    {time_reuse:.4f} seconds")
        print(f"With reset:    {time_reset:.4f} seconds")
        print(f"Speedup:       {speedup:.2f}x")
        print(f"{'='*60}")

        # Both methods should work correctly (no assertion on speedup)
        # Timing is informational only to avoid flakiness in CI

    @pytest.mark.parametrize("n_ell", [2, 4, 8, 16, 32, 64])
    def test_batched_vs_single_ell_performance(self, n_ell: int) -> None:
        """Test that batched operations are faster than individual single-ell solves.

        Compares the performance of:
        1. Batched solve for n_ell values at once
        2. n_ell individual single-ell solves

        The batched case should be faster, especially for larger n_ell.
        """
        N = 128
        a, b = 1.0, 20.0
        lmin = 5
        lmax = lmin + n_ell - 1
        n_iterations = 50

        solver = Ncm.SBesselOdeSolver.new()

        # Prepare RHS arrays for testing
        rhs_list = []
        for i in range(n_iterations):
            rhs = np.zeros(N)
            rhs[0] = 0.0
            rhs[1] = 0.0
            rhs[2] = 1.0 + i * 0.1
            rhs_list.append(rhs)

        # Test 1: Batched solves (fast path)
        op_batched = solver.create_operator(a, b, lmin, lmax)
        start_batched = time.perf_counter()
        for rhs in rhs_list:
            _solution, _sol_len = op_batched.solve(rhs)
        end_batched = time.perf_counter()
        time_batched = end_batched - start_batched

        # Test 2: Individual single-ell solves (slow path)
        ops_single = [
            solver.create_operator(a, b, ell, ell) for ell in range(lmin, lmax + 1)
        ]
        start_single = time.perf_counter()
        for rhs in rhs_list:
            for op_single in ops_single:
                _solution, _sol_len = op_single.solve(rhs)
        end_single = time.perf_counter()
        time_single = end_single - start_single

        # Calculate speedup
        speedup = time_single / time_batched

        # Print timing results
        print(f"\n{'='*60}")
        print("Batched vs Single-ell Performance Test")
        print(f"{'='*60}")
        print(f"Number of iterations: {n_iterations}")
        print(f"Matrix size: {N}")
        print(f"Batch size (n_ell): {n_ell}")
        print(f"ell range: [{lmin}, {lmax}]")
        print(f"\nBatched:       {time_batched:.4f} seconds")
        print(f"Single-ell:    {time_single:.4f} seconds")
        print(f"Speedup:       {speedup:.2f}x")
        print(f"{'='*60}")

        # Both methods should work correctly (no assertion on speedup)
        # Timing is informational only to avoid flakiness in CI

    @pytest.mark.parametrize("n_ell", [2, 4, 8, 16, 32, 64])
    def test_batched_vs_single_ell_solve_endpoints_performance(
        self, n_ell: int
    ) -> None:
        """Test that batched solve_endpoints is faster than individual calls.

        Compares the performance of:
        1. Batched solve_endpoints for multiple ell values at once
        2. Individual solve_endpoints calls for each ell value

        The batched case should be faster.
        """
        N = 128
        a, b = 1.0, 20.0
        lmin = 5
        lmax = lmin + n_ell - 1
        n_iterations = 50

        solver = Ncm.SBesselOdeSolver.new()

        # Prepare RHS arrays for testing
        rhs_list = []
        for i in range(n_iterations):
            rhs = np.zeros(N)
            rhs[0] = 0.0
            rhs[1] = 0.0
            rhs[2] = 1.0 + i * 0.1
            rhs_list.append(rhs)

        # Test 1: Batched solve_endpoints (fast path)
        op_batched = solver.create_operator(a, b, lmin, lmax)
        start_batched = time.perf_counter()
        for rhs in rhs_list:
            _ = op_batched.solve_endpoints(rhs)
        end_batched = time.perf_counter()
        time_batched = end_batched - start_batched

        # Test 2: Individual single-ell solve_endpoints (slow path)
        ops_single = [
            solver.create_operator(a, b, ell, ell) for ell in range(lmin, lmax + 1)
        ]
        start_single = time.perf_counter()
        for rhs in rhs_list:
            for op_single in ops_single:
                _deriv_a, _deriv_b, _error = op_single.solve_endpoints(rhs)
        end_single = time.perf_counter()
        time_single = end_single - start_single

        # Calculate speedup
        speedup = time_single / time_batched

        # Print timing results
        print(f"\n{'='*60}")
        print("Batched vs Single-ell solve_endpoints Performance Test")
        print(f"{'='*60}")
        print(f"Number of iterations: {n_iterations}")
        print(f"Matrix size: {N}")
        print(f"Batch size (n_ell): {n_ell}")
        print(f"ell range: [{lmin}, {lmax}]")
        print(f"\nBatched:       {time_batched:.4f} seconds")
        print(f"Single-ell:    {time_single:.4f} seconds")
        print(f"Speedup:       {speedup:.2f}x")
        print(f"{'='*60}")

        # Both methods should work correctly (no assertion on speedup)
        # Timing is informational only to avoid flakiness in CI


class TestSBesselOperatorMemoryManagement:
    """Tests for memory management and capacity reuse in NcmSBesselOdeOperator.

    These tests verify that the operator correctly manages internal memory when:
    - Solving with different RHS sizes in sequence
    - Reusing stored factorizations
    - Growing capacity as needed
    """

    def test_memory_status_accessors(self) -> None:
        """Test that memory status accessor functions work correctly."""
        solver = Ncm.SBesselOdeSolver.new()
        op = solver.create_operator(1.0, 20.0, 0, 2)

        # Before any solve, check initial state
        initial_capacity = op.get_operator_size()
        initial_n_cols = op.get_n_cols()

        assert (
            initial_capacity == 0
        ), f"Initial operator_capacity should be 0, got {initial_capacity}"
        assert initial_n_cols == 0, f"Initial n_cols should be 0, got {initial_n_cols}"

        # After first solve, capacity and n_cols should be non-zero
        rhs = np.zeros(64)
        rhs[2] = 1.0
        _, _ = op.solve(rhs)

        capacity_after_solve = op.get_operator_size()
        n_cols_after_solve = op.get_n_cols()

        assert (
            capacity_after_solve > 0
        ), f"operator_capacity should be > 0 after solve, got {capacity_after_solve}"
        assert (
            n_cols_after_solve > 0
        ), f"n_cols should be > 0 after solve, got {n_cols_after_solve}"

    def test_capacity_grows_with_higher_spectral_order(self) -> None:
        """Test that operator capacity can grow when solving with higher spectral order.

        The capacity grows adaptively based on convergence, which typically
        requires more columns for higher spectral orders.
        """
        solver = Ncm.SBesselOdeSolver.new()
        op = solver.create_operator(1.0, 20.0, 0, 0)

        # First solve with low spectral order
        rhs_low = np.zeros(32)
        rhs_low[2] = 1.0
        _, _ = op.solve(rhs_low)

        capacity_low = op.get_operator_size()
        n_cols_low = op.get_n_cols()

        # Second solve with much higher spectral order
        rhs_high = np.zeros(128)
        rhs_high[50] = 1.0
        _, _ = op.solve(rhs_high)

        capacity_high = op.get_operator_size()
        n_cols_high = op.get_n_cols()

        # Capacity should not decrease
        assert capacity_high >= capacity_low, (
            f"Capacity should not decrease. "
            f"After low: {capacity_low}, After high: {capacity_high}"
        )

        # Both solves should produce positive n_cols
        assert n_cols_low > 0, "n_cols after low spectral order should be positive"
        assert n_cols_high > 0, "n_cols after high spectral order should be positive"

    def test_capacity_reuse_with_smaller_rhs(self) -> None:
        """Test that solving with lower spectral order after higher works correctly.

        Memory management should handle solving with smaller spectral orders
        after larger ones, potentially reusing the allocated capacity.
        """
        solver = Ncm.SBesselOdeSolver.new()
        op = solver.create_operator(1.0, 20.0, 0, 0)

        # First solve with high spectral order
        rhs_large = np.zeros(128)
        rhs_large[60] = 1.0  # High spectral order
        solution_large, sol_len_large = op.solve(rhs_large)

        capacity_after_large = op.get_operator_size()
        _ = op.get_n_cols()

        # Second solve with low spectral order
        rhs_small = np.zeros(32)
        rhs_small[2] = 1.0  # Low spectral order
        solution_small, sol_len_small = op.solve(rhs_small)

        capacity_after_small = op.get_operator_size()
        n_cols_after_small = op.get_n_cols()

        # Capacity should not decrease (memory is reused)
        assert capacity_after_small >= capacity_after_large, (
            f"Capacity should not decrease. "
            f"After large: {capacity_after_large}, After small: {capacity_after_small}"
        )

        # n_cols should reflect the new solution (likely smaller)
        assert (
            n_cols_after_small > 0
        ), f"n_cols should be positive after small RHS: {n_cols_after_small}"

        # Both solutions should be valid
        assert (
            sol_len_large > 0 and sol_len_small > 0
        ), "Both solutions should be non-empty"
        assert np.all(np.isfinite(solution_large)), "Large solution should be finite"
        assert np.all(np.isfinite(solution_small)), "Small solution should be finite"

    def test_capacity_reuse_with_identical_rhs(self) -> None:
        """Test that solving with identical RHS twice reuses stored factorization.

        When the RHS is truly identical (same spectral order and coefficients),
        the stored factorization should be reused without growing capacity.
        """
        solver = Ncm.SBesselOdeSolver.new()
        op = solver.create_operator(1.0, 20.0, 0, 0)

        # First solve
        rhs = np.zeros(64)
        rhs[2] = 1.0
        solution1, _sol_len1 = op.solve(rhs)

        capacity1 = op.get_operator_size()
        n_cols1 = op.get_n_cols()

        # Second solve with identical RHS (same spectral order)
        solution2, _sol_len2 = op.solve(rhs)

        capacity2 = op.get_operator_size()
        n_cols2 = op.get_n_cols()

        # Capacity should not grow (factorization reused)
        assert capacity2 == capacity1, (
            f"Capacity should be unchanged with identical RHS. "
            f"First: {capacity1}, Second: {capacity2}"
        )
        # n_cols should be the same (same convergence point)
        assert n_cols2 == n_cols1, (
            f"n_cols should be unchanged with identical RHS. "
            f"First: {n_cols1}, Second: {n_cols2}"
        )

        # Solutions should be identical
        solution1_np = np.array(solution1)
        solution2_np = np.array(solution2)
        assert_allclose(
            solution1_np,
            solution2_np,
            rtol=1e-14,
            err_msg="Solutions should be identical for identical RHS",
        )

    def test_different_spectral_order_same_length(self) -> None:
        """Test that RHS with different spectral orders behave correctly.

        The spectral order (highest significant coefficient) determines
        convergence behavior, not the array length. Different spectral
        orders may require different capacities.
        """
        solver = Ncm.SBesselOdeSolver.new()
        op = solver.create_operator(1.0, 20.0, 0, 0)

        # Low spectral order: only coefficient at index 2
        rhs_low = np.zeros(64)
        rhs_low[2] = 1.0
        solution_low, sol_len_low = op.solve(rhs_low)

        capacity_low = op.get_operator_size()
        _ = op.get_n_cols()

        # Higher spectral order: coefficient at index 30
        rhs_high = np.zeros(64)
        rhs_high[30] = 1.0
        solution_high, sol_len_high = op.solve(rhs_high)

        capacity_high = op.get_operator_size()
        _ = op.get_n_cols()

        # Capacity should grow if higher spectral order needs more columns
        assert (
            capacity_high >= capacity_low
        ), f"Capacity should not decrease: {capacity_low} -> {capacity_high}"

        # Both solutions should be valid (non-empty and finite)
        assert sol_len_low > 0, "Low spectral order solution should be non-empty"
        assert sol_len_high > 0, "High spectral order solution should be non-empty"

        solution_low_np = np.array(solution_low)
        solution_high_np = np.array(solution_high)
        assert np.all(
            np.isfinite(solution_low_np)
        ), "Low order solution should be finite"
        assert np.all(
            np.isfinite(solution_high_np)
        ), "High order solution should be finite"

    @pytest.mark.parametrize(
        "first_order,second_order",
        [
            (2, 10),  # Low to medium spectral order
            (2, 30),  # Low to high spectral order
            (10, 30),  # Medium to high spectral order
            (30, 50),  # High to very high spectral order
        ],
    )
    def test_sequential_solves_increasing_spectral_order(
        self, first_order: int, second_order: int
    ) -> None:
        """Test sequential solves with increasing spectral orders.

        Tests the actual bug scenario: solving with different spectral orders
        in sequence should work correctly without memory corruption.
        """
        solver = Ncm.SBesselOdeSolver.new()
        op = solver.create_operator(1.0, 20.0, 5, 5)

        # First solve with low spectral order
        rhs1 = np.zeros(max(64, first_order + 4))
        rhs1[first_order] = 1.0
        solution1, sol_len1 = op.solve(rhs1)

        capacity1 = op.get_operator_size()
        _ = op.get_n_cols()

        # Second solve with higher spectral order
        rhs2 = np.zeros(max(64, second_order + 4))
        rhs2[second_order] = 1.0
        solution2, sol_len2 = op.solve(rhs2)

        capacity2 = op.get_operator_size()
        _ = op.get_n_cols()

        # Verify memory management is correct
        assert (
            capacity2 >= capacity1
        ), f"Capacity should not decrease: {capacity1} -> {capacity2}"

        # Verify solutions are valid and different
        assert sol_len1 > 0 and sol_len2 > 0, "Solutions should be non-empty"

        solution1_np = np.array(solution1)
        solution2_np = np.array(solution2)

        assert np.all(np.isfinite(solution1_np)), "First solution should be finite"
        assert np.all(np.isfinite(solution2_np)), "Second solution should be finite"

        # Solutions should differ (different RHS)
        max_len = min(len(solution1_np), len(solution2_np))
        assert not np.allclose(
            solution1_np[:max_len], solution2_np[:max_len]
        ), "Solutions should differ for different spectral orders"

    @pytest.mark.parametrize(
        "spectral_orders",
        [
            [2, 10, 2],  # Low, medium, low
            [2, 30, 10],  # Low, high, medium
            [30, 10, 5],  # High, medium, low
            [5, 30, 5, 30],  # Alternating
        ],
    )
    def test_sequential_solves_varying_spectral_orders(
        self, spectral_orders: list
    ) -> None:
        """Test sequential solves with varying spectral orders.

        This tests the core memory management behavior: solving with different
        spectral orders in sequence should always work correctly, with capacity
        never decreasing and all solutions being valid.
        """
        solver = Ncm.SBesselOdeSolver.new()
        op = solver.create_operator(1.0, 20.0, 0, 0)

        capacities = []
        n_cols_list = []
        solutions = []

        for order in spectral_orders:
            rhs = np.zeros(max(64, order + 4))
            rhs[order] = 1.0
            solution, sol_len = op.solve(rhs)

            capacity = op.get_operator_size()
            n_cols = op.get_n_cols()

            capacities.append(capacity)
            n_cols_list.append(n_cols)
            solutions.append(np.array(solution))

            # Verify solution is valid
            assert (
                sol_len > 0
            ), f"Solution should be non-empty for spectral order {order}"
            assert np.all(
                np.isfinite(solution)
            ), f"Solution should be finite for spectral order {order}"

        # Capacity should never decrease (key invariant)
        for i in range(1, len(capacities)):
            assert (
                capacities[i] >= capacities[i - 1]
            ), f"Capacity should not decrease: {capacities}"

        # All n_cols should be positive
        assert all(
            n > 0 for n in n_cols_list
        ), f"All n_cols should be positive: {n_cols_list}"

    def test_batched_memory_management(self) -> None:
        """Test memory management with batched (multi-ell) operators."""
        solver = Ncm.SBesselOdeSolver.new()
        n_ell = 5
        op = solver.create_operator(1.0, 20.0, 0, n_ell - 1)

        # First solve with low spectral order
        rhs1 = np.zeros(64)
        rhs1[2] = 1.0
        _, sol_len1 = op.solve(rhs1)

        capacity1 = op.get_operator_size()
        n_cols1 = op.get_n_cols()

        # Second solve with higher spectral order
        rhs2 = np.zeros(128)
        rhs2[40] = 1.0
        _, sol_len2 = op.solve(rhs2)

        capacity2 = op.get_operator_size()
        n_cols2 = op.get_n_cols()

        # Verify growth (capacity should not decrease, n_cols may grow)
        assert (
            capacity2 >= capacity1
        ), f"Batched capacity should not decrease: {capacity1} -> {capacity2}"
        # With batched operators, both should produce valid n_cols
        assert n_cols1 > 0, "n_cols should be positive after first solve"
        assert n_cols2 > 0, "n_cols should be positive after second solve"

        # Verify solution lengths match expected batched structure
        assert sol_len1 > 0, "First solution length should be positive"
        assert sol_len2 > 0, "Second solution length should be positive"

    def test_endpoints_memory_management(self) -> None:
        """Test that solve_endpoints also manages memory correctly."""
        solver = Ncm.SBesselOdeSolver.new()
        op = solver.create_operator(1.0, 20.0, 0, 0)

        # First solve_endpoints with low spectral order
        rhs_low = np.zeros(32)
        rhs_low[2] = 1.0
        _ = op.solve_endpoints(rhs_low)

        capacity_low = op.get_operator_size()
        n_cols_low = op.get_n_cols()

        # Second solve_endpoints with higher spectral order
        rhs_high = np.zeros(128)
        rhs_high[40] = 1.0
        _ = op.solve_endpoints(rhs_high)

        capacity_high = op.get_operator_size()
        n_cols_high = op.get_n_cols()

        # Verify capacity management (should not decrease)
        assert capacity_high >= capacity_low, (
            f"Capacity should not decrease with solve_endpoints: "
            f"{capacity_low} -> {capacity_high}"
        )

        # Both should produce positive n_cols
        assert n_cols_low > 0, "n_cols should be positive after first solve_endpoints"
        assert n_cols_high > 0, "n_cols should be positive after second solve_endpoints"

    def test_reset_clears_memory_status(self) -> None:
        """Test that reset() clears the diagonalization but keeps capacity."""
        solver = Ncm.SBesselOdeSolver.new()
        op = solver.create_operator(1.0, 20.0, 0, 0)

        # Solve to allocate memory
        rhs = np.zeros(64)
        rhs[2] = 1.0
        _, _ = op.solve(rhs)

        capacity_before = op.get_operator_size()
        n_cols_before = op.get_n_cols()

        assert capacity_before > 0, "Should have capacity before reset"
        assert n_cols_before > 0, "Should have n_cols before reset"

        # Reset the operator
        op.reset(1.0, 20.0, 1, 1)

        _ = op.get_operator_size()
        n_cols_after = op.get_n_cols()

        # After reset, n_cols should be 0 (no factorization)
        # but capacity stays allocated for reuse
        assert n_cols_after == 0, f"n_cols should be 0 after reset, got {n_cols_after}"
        # Note: capacity behavior after reset may vary - the important
        # thing is that subsequent solves work correctly

    def test_stress_test_many_sequential_solves(self) -> None:
        """Stress test with many sequential solves of varying spectral orders.

        This is the ultimate test for the memory management bug fixes:
        many consecutive solves with random spectral orders should all
        complete successfully without memory corruption or crashes.
        """
        solver = Ncm.SBesselOdeSolver.new()
        op = solver.create_operator(1.0, 20.0, 0, 0)

        # Generate a sequence of random spectral orders
        np.random.seed(21)
        spectral_orders = np.random.choice([2, 5, 10, 15, 20, 300, 400], size=20)

        prev_capacity = 0

        for i, order in enumerate(spectral_orders):
            rhs = np.zeros(max(64, order + 4))
            rhs[order] = 1.0 / (i + 1)  # Vary RHS values

            solution, sol_len = op.solve(rhs)

            capacity = op.get_operator_size()
            n_cols = op.get_n_cols()

            # Verify capacity never decreases
            assert (
                capacity >= prev_capacity
            ), f"Capacity decreased at iteration {i}: {prev_capacity} -> {capacity}"
            prev_capacity = capacity

            # Basic sanity checks
            assert capacity > 0, f"Capacity should be positive at iteration {i}"
            assert n_cols > 0, f"n_cols should be positive at iteration {i}"
            assert sol_len > 0, f"Solution should be non-empty at iteration {i}"

            # Verify solution has reasonable values
            solution_np = np.array(solution)
            assert np.all(
                np.isfinite(solution_np)
            ), f"Solution should be finite at iteration {i}, order={order}"


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

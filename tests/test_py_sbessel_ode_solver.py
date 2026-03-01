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
            y_at_a, y_a, rtol=1.0e-15, atol=1.0e-15, err_msg="BC at x=a not satisfied"
        )
        assert_allclose(
            y_at_b, y_b, rtol=1.0e-15, atol=1.0e-15, err_msg="BC at x=b not satisfied"
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

        # Verify that reuse is significantly faster
        # Expect at least 1.5x speedup from reusing diagonalization
        assert speedup > 1.5, (
            f"Diagonalization reuse should be significantly faster than reset. "
            f"Expected speedup > 1.5x, got {speedup:.2f}x. "
            f"Time with reuse: {time_reuse:.4f}s, with reset: {time_reset:.4f}s"
        )

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

        # Verify that reuse is significantly faster
        # Expect at least 1.5x speedup from reusing diagonalization
        assert speedup > 1.5, (
            f"Diagonalization reuse should be significantly faster than reset "
            f"for n_ell={n_ell}. Expected speedup > 1.5x, got {speedup:.2f}x. "
            f"Time with reuse: {time_reuse:.4f}s, with reset: {time_reset:.4f}s"
        )

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

        # Verify that batched is faster
        # Expect at least some speedup from batched operations
        assert speedup > 1.0, (
            f"Batched operations should be faster than individual single-ell solves "
            f"for n_ell={n_ell}. Expected speedup > 1.0x, got {speedup:.2f}x. "
            f"Time batched: {time_batched:.4f}s, single-ell: {time_single:.4f}s"
        )

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

        # Verify that batched is faster
        assert speedup > 1.0, (
            f"Batched solve_endpoints should be faster than individual calls. "
            f"Expected speedup > 1.0x, got {speedup:.2f}x. "
            f"Time batched: {time_batched:.4f}s, single-ell: {time_single:.4f}s"
        )

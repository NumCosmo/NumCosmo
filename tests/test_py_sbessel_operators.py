#
# test_py_sbessel_operators.py
#
# Fri Jan 31 00:00:00 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_py_sbessel_operators.py
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

"""Tests for NcmSBesselOdeSolver operator matrices."""

from pathlib import Path
import gzip
import json

import numpy as np
import pytest

from numpy.testing import assert_allclose
from scipy.special import spherical_jn
from scipy.linalg import solve
from scipy.integrate import quad

from numcosmo_py import Ncm


class TestSBesselOperators:
    """Tests for NcmSBesselOdeSolver operator matrices."""

    @pytest.mark.parametrize("l_val", list(range(21)))
    def test_spherical_bessel_ode(self, l_val: int) -> None:
        """Test solving spherical Bessel ODE with non-zero Dirichlet BCs.

        The spherical Bessel function j_l(x) satisfies:
        x^2*y'' + 2xy' + [x^2 - l(l+1)]y = 0

        We solve this on interval [1,20] with BCs y(1) = j_l(1), y(20) = j_l(20)
        and verify the solution matches j_l everywhere.
        """
        N = 128
        a, b = 1.0, 20.0  # Use interval [1, 20] as requested

        # Create solver with interval [a, b]
        solver = Ncm.SBesselOdeSolver.new(l_val, a, b)

        # Get the operator matrix (should now be square NxN)
        mat = solver.get_operator_matrix(N)
        mat_np = mat.to_numpy()

        # Get boundary values using scipy
        y_a = spherical_jn(l_val, a)
        y_b = spherical_jn(l_val, b)

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
        solution_coeffs_np = solve(mat_np, rhs_np)
        solution_coeffs = Ncm.Vector.new_array(solution_coeffs_np.tolist())

        # Test at several interior points mapped to [-1, 1]
        # Mapping: x_physical = (a+b)/2 + (b-a)/2 * t, where t in [-1, 1]
        t_test = np.array([-0.8, -0.5, -0.2, 0.0, 0.2, 0.5, 0.8])

        for t in t_test:
            # Evaluate solution at t (coefficients are in Chebyshev basis)
            y_computed = Ncm.Spectral.chebyshev_eval(solution_coeffs, t)

            # Map t to physical x
            x_physical = (a + b) / 2.0 + (b - a) / 2.0 * t

            # Get exact value
            y_exact = spherical_jn(l_val, x_physical)

            assert_allclose(
                y_computed,
                y_exact,
                rtol=1.0e-10,
                atol=1.0e-16,
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

        Solve the ODE with homogeneous BCs (y(a)=0, y(b)=0) and RHS = 1.
        By Green's identity: x^2*j_l(b)*y'(b) - x^2*j_l(a)*y'(a) = int_a^b j_l(x)dx
        where y is the solution with RHS=1.
        """
        N = 128
        a, b = 1.0, 20.0

        # Create solver with interval [a, b]
        solver = Ncm.SBesselOdeSolver.new(l_val, a, b)

        # Get the operator matrix
        mat = solver.get_operator_matrix(N)
        mat_np = mat.to_numpy()

        # Set up RHS: homogeneous BCs (y(a)=0, y(b)=0) with RHS=1
        rhs_np = np.zeros(N)
        rhs_np[0] = 0.0  # BC at x=a (t=-1)
        rhs_np[1] = 0.0  # BC at x=b (t=+1)
        rhs_np[2] = 1.0
        # Rows 3 to N-1 have RHS = 0

        # Solve the linear system
        solution_coeffs_np = solve(mat_np, rhs_np)
        solution_coeffs = Ncm.Vector.new_array(solution_coeffs_np.tolist())

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

        # Compute left-hand side: x^2*j_l(b)*y'(b) - x^2*j_l(a)*y'(a)
        lhs = b * b * j_l_b * y_prime_b - a * a * j_l_a * y_prime_a

        # Compute right-hand side: integral of j_l(x) from a to b
        # Use high-accuracy numerical integration

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
                f"x^2*j_l(b)*y'(b) - x^2*j_l(a)*y'(a) != int_a^b j_l(x)dx"
            ),
        )

    @pytest.mark.parametrize("l_val", list(range(21)))
    def test_spherical_bessel_integration_x(self, l_val: int) -> None:
        """Test integration of x*j_l(x) via Green's identity.

        Solve the ODE with homogeneous BCs (y(a)=0, y(b)=0) and RHS = x.
        By Green's identity: x^2*j_l(b)*y'(b) - x^2*j_l(a)*y'(a) = int_a^b x*j_l(x)dx
        where y is the solution with RHS=x.
        """
        N = 128
        a, b = 1.0, 20.0

        # Create solver with interval [a, b]
        solver = Ncm.SBesselOdeSolver.new(l_val, a, b)

        # Get the operator matrix
        mat = solver.get_operator_matrix(N)
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
        solution_coeffs_np = solve(mat_np, rhs_np)
        solution_coeffs = Ncm.Vector.new_array(solution_coeffs_np.tolist())

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

        # Compute left-hand side: x^2*j_l(b)*y'(b) - x^2*j_l(a)*y'(a)
        lhs = b * b * j_l_b * y_prime_b - a * a * j_l_a * y_prime_a

        # Compute right-hand side: integral of x*j_l(x) from a to b
        def integrand(x: float) -> float:
            return x * spherical_jn(l_val, x)

        rhs, _ = quad(integrand, a, b, epsabs=1e-12, epsrel=1e-14)

        # Compare
        assert_allclose(
            lhs,
            rhs,
            rtol=1.0e-8,
            atol=1.0e-20,
            err_msg=(
                f"Green's identity failed for l={l_val}: "
                f"x^2*j_l(b)*y'(b) - x^2*j_l(a)*y'(a) != int_a^b x*j_l(x)dx"
            ),
        )

    @pytest.mark.parametrize("l_val", list(range(21)))
    @pytest.mark.parametrize("N", [32, 64, 128])
    def test_solve_qr_vs_scipy(self, l_val: int, N: int) -> None:
        """Test that QR solve produces the same result as scipy solver."""
        a, b = 1.0, 20.0

        # Create solver with interval [a, b]
        solver = Ncm.SBesselOdeSolver.new(l_val, a, b)

        # Set up simple RHS vector with just a few coefficients
        rhs_np = np.zeros(N)
        rhs_np[2] = 1.0
        rhs_vec = Ncm.Vector.new_array(rhs_np.tolist())

        # Call solve (QR method) - just check it doesn't crash
        solution = solver.solve(rhs_vec)
        solution_np = solution.to_numpy()

        op_mat = solver.get_operator_matrix(N).to_numpy()

        # Solve the system using scipy for comparison
        rhs_np_full = np.zeros(N)
        rhs_np_full[2] = 1.0
        solution_scipy = solve(op_mat, rhs_np_full)
        safe_part = 20

        assert_allclose(
            solution_np[:safe_part],
            solution_scipy[:safe_part],
            rtol=1.0e-10,
            atol=0.0,
            err_msg="QR solve doesn't match scipy solution",
        )

        # Basic sanity check - solution should exist and be non-trivial
        assert solution is not None
        assert len(solution_np) > 0

    @pytest.mark.parametrize("l_val", list(range(21)))
    @pytest.mark.parametrize("N", [32, 64, 128])
    def test_solve_dense_vs_scipy(self, l_val: int, N: int) -> None:
        """Test that solve_dense produces the same result as scipy solver.

        This tests the LAPACK dgesv-based solve_dense method against scipy's
        solve to verify the column-major matrix format is correct.
        """
        a, b = 1.0, 20.0

        # Create solver with interval [a, b]
        solver = Ncm.SBesselOdeSolver.new(l_val, a, b)

        # Set up RHS: homogeneous BCs with RHS=1
        rhs_np = np.zeros(N)
        rhs_np[0] = 0.0  # BC at x=a
        rhs_np[1] = 0.0  # BC at x=b
        rhs_np[2] = 1.0
        rhs_vec = Ncm.Vector.new_array(rhs_np.tolist())

        # Solve using scipy
        mat = solver.get_operator_matrix(N)
        mat_np = mat.to_numpy()
        solution_scipy = solve(mat_np, rhs_np)

        # Solve using solve_dense (LAPACK dgesv)
        solution_dense = solver.solve_dense(rhs_vec, N)
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

        Verifies that solve_dense produces solutions that satisfy Green's
        identity: x^2*j_l(b)*y'(b) - x^2*j_l(a)*y'(a) = int_a^b j_l(x)dx
        """
        N = 128
        a, b = 1.0, 20.0

        # Create solver
        solver = Ncm.SBesselOdeSolver.new(l_val, a, b)

        # Set up RHS: homogeneous BCs with RHS=1
        rhs_np = np.zeros(N)
        rhs_np[0] = 0.0  # BC at x=a
        rhs_np[1] = 0.0  # BC at x=b
        rhs_np[2] = 1.0
        rhs_vec = Ncm.Vector.new_array(rhs_np.tolist())

        # Solve using solve_dense
        solution_coeffs = solver.solve_dense(rhs_vec, N)

        # Compute derivatives at endpoints
        h = (b - a) / 2.0
        y_prime_at_minus1 = Ncm.Spectral.chebyshev_deriv(solution_coeffs, -1.0)
        y_prime_at_plus1 = Ncm.Spectral.chebyshev_deriv(solution_coeffs, 1.0)
        y_prime_a = y_prime_at_minus1 / h
        y_prime_b = y_prime_at_plus1 / h

        # Get j_l values
        j_l_a = spherical_jn(l_val, a)
        j_l_b = spherical_jn(l_val, b)

        # Compute LHS: x^2*j_l(b)*y'(b) - x^2*j_l(a)*y'(a)
        lhs = b * b * j_l_b * y_prime_b - a * a * j_l_a * y_prime_a

        # Compute RHS: integral
        def integrand(x: float) -> float:
            return spherical_jn(l_val, x)

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

        Verifies that solve_dense produces solutions that satisfy Green's
        identity: x^2*j_l(b)*y'(b) - x^2*j_l(a)*y'(a) = int_a^b x*j_l(x)dx
        """
        N = 128
        a, b = 1.0, 20.0

        # Create solver
        solver = Ncm.SBesselOdeSolver.new(l_val, a, b)

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
        solution_coeffs = solver.solve_dense(rhs_vec, N)

        # Compute derivatives at endpoints
        y_prime_at_minus1 = Ncm.Spectral.chebyshev_deriv(solution_coeffs, -1.0)
        y_prime_at_plus1 = Ncm.Spectral.chebyshev_deriv(solution_coeffs, 1.0)
        y_prime_a = y_prime_at_minus1 / h
        y_prime_b = y_prime_at_plus1 / h

        # Get j_l values
        j_l_a = spherical_jn(l_val, a)
        j_l_b = spherical_jn(l_val, b)

        # Compute LHS: x^2*j_l(b)*y'(b) - x^2*j_l(a)*y'(a)
        lhs = b * b * j_l_b * y_prime_b - a * a * j_l_a * y_prime_a

        # Compute RHS: integral of x*j_l(x)
        def integrand(x: float) -> float:
            return x * spherical_jn(l_val, x)

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

        Verifies that solve (QR method) produces solutions that satisfy Green's
        identity: x^2*j_l(b)*y'(b) - x^2*j_l(a)*y'(a) = int_a^b j_l(x)dx
        """
        N = 128
        a, b = 1.0, 20.0

        # Create solver
        solver = Ncm.SBesselOdeSolver.new(l_val, a, b)

        # Set up RHS: homogeneous BCs with RHS=1
        rhs_np = np.zeros(N)
        rhs_np[0] = 0.0  # BC at x=a
        rhs_np[1] = 0.0  # BC at x=b
        rhs_np[2] = 1.0
        rhs_vec = Ncm.Vector.new_array(rhs_np.tolist())

        # Solve using QR method
        solution_coeffs = solver.solve(rhs_vec)

        # Compute derivatives at endpoints
        h = (b - a) / 2.0
        y_prime_at_minus1 = Ncm.Spectral.chebyshev_deriv(solution_coeffs, -1.0)
        y_prime_at_plus1 = Ncm.Spectral.chebyshev_deriv(solution_coeffs, 1.0)
        y_prime_a = y_prime_at_minus1 / h
        y_prime_b = y_prime_at_plus1 / h

        # Get j_l values
        j_l_a = spherical_jn(l_val, a)
        j_l_b = spherical_jn(l_val, b)

        # Compute LHS: x^2*j_l(b)*y'(b) - x^2*j_l(a)*y'(a)
        lhs = b * b * j_l_b * y_prime_b - a * a * j_l_a * y_prime_a

        # Compute RHS: integral
        def integrand(x: float) -> float:
            return spherical_jn(l_val, x)

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

        Verifies that solve (QR method) produces solutions that satisfy Green's
        identity: x^2*j_l(b)*y'(b) - x^2*j_l(a)*y'(a) = int_a^b x*j_l(x)dx
        """
        N = 128
        a, b = 1.0, 20.0

        # Create solver
        solver = Ncm.SBesselOdeSolver.new(l_val, a, b)

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

        # Solve using QR method
        solution_coeffs = solver.solve(rhs_vec)

        # Compute derivatives at endpoints
        y_prime_at_minus1 = Ncm.Spectral.chebyshev_deriv(solution_coeffs, -1.0)
        y_prime_at_plus1 = Ncm.Spectral.chebyshev_deriv(solution_coeffs, 1.0)
        y_prime_a = y_prime_at_minus1 / h
        y_prime_b = y_prime_at_plus1 / h

        # Get j_l values
        j_l_a = spherical_jn(l_val, a)
        j_l_b = spherical_jn(l_val, b)

        # Compute LHS: x^2*j_l(b)*y'(b) - x^2*j_l(a)*y'(a)
        lhs = b * b * j_l_b * y_prime_b - a * a * j_l_a * y_prime_a

        # Compute RHS: integral of x*j_l(x)
        def integrand(x: float) -> float:
            return x * spherical_jn(l_val, x)

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
    def test_integrate_constant(self, l_val: int) -> None:
        """Test integrate method with f(x) = 1.

        Compares the result of the integrate method against scipy numerical integration.
        """
        N = 128
        a, b = 1.0, 20.0

        # Create solver
        solver = Ncm.SBesselOdeSolver.new(l_val, a, b)

        # Define f(x) = 1
        def f_constant(_user_data: None, _x: float) -> float:
            return 1.0

        # Compute integral using the new integrate method
        result = solver.integrate(f_constant, N, None)

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
            err_msg=(f"integrate method failed for l={l_val}, f(x)=1"),
        )

    @pytest.mark.parametrize("l_val", [0, 1, 5, 10, 20])
    def test_integrate_linear(self, l_val: int) -> None:
        """Test integrate method with f(x) = x.

        Compares the result of the integrate method against scipy numerical integration.
        """
        N = 128
        a, b = 1.0, 20.0

        # Create solver
        solver = Ncm.SBesselOdeSolver.new(l_val, a, b)

        # Define f(x) = x
        def f_linear(_user_data: None, x: float) -> float:
            return x

        # Compute integral using the new integrate method
        result = solver.integrate(f_linear, N, None)

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
            err_msg=(f"integrate method failed for l={l_val}, f(x)=x"),
        )

    @pytest.mark.parametrize("l_val", [0, 1, 5, 10])
    def test_integrate_quadratic(self, l_val: int) -> None:
        """Test integrate method with f(x) = x^2.

        Compares the result of the integrate method against scipy numerical integration.
        """
        N = 128
        a, b = 1.0, 20.0

        # Create solver
        solver = Ncm.SBesselOdeSolver.new(l_val, a, b)

        # Define f(x) = x^2
        def f_quadratic(_user_data: None, x: float) -> float:
            return x * x

        # Compute integral using the new integrate method
        result = solver.integrate(f_quadratic, N, None)

        # Compute expected value using scipy
        def integrand(x: float) -> float:
            return x * x * spherical_jn(l_val, x)

        expected, _ = quad(integrand, a, b, epsabs=1e-11, epsrel=1e-14)

        # Compare
        assert_allclose(
            result,
            expected,
            rtol=1.0e-8,
            atol=1.0e-20,
            err_msg=(f"integrate method failed for l={l_val}, f(x)=x^2"),
        )

    @pytest.mark.parametrize("l_val", [0, 1, 5, 10, 20, 200])
    def test_integrate_rational(self, l_val: int) -> None:
        """Test integrate method with f(x) = x^2/(1+x^2)^4.

        Compares the result of the integrate method against scipy numerical integration.
        """
        N = 256
        a, b = 1.0, 20.0

        # Create solver
        solver = Ncm.SBesselOdeSolver.new(l_val, a, b)

        # Define f(x) = x^2/(1+x^2)^4
        def f_rational(_user_data: None, x: float) -> float:
            return 1.0e20 * x**2 / (1.0 + x**2) ** 4

        # Compute integral using the new integrate method
        result = solver.integrate(f_rational, N, None)

        # Compute expected value using scipy
        def integrand(x: float) -> float:
            return f_rational(None, x) * spherical_jn(l_val, x)

        expected, _ = quad(integrand, a, b, epsabs=1e-11, epsrel=1e-24)

        # Compare
        assert_allclose(
            result,
            expected,
            rtol=1.0e-7,
            atol=1.0e-18,
            err_msg=(f"integrate method failed for l={l_val}, f(x)=x^2/(1+x^2)^4"),
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
        rhs = Ncm.Vector.new(N)
        rhs_data = np.zeros(N)
        rhs_data[0] = 0.0  # BC at -1
        rhs_data[1] = 0.0  # BC at +1
        # Set some non-trivial coefficients
        for k in range(2, N):
            rhs_data[k] = 1.0 / (1.0 + k * k)
        rhs.set_array(rhs_data.tolist())

        # Create solver
        solver = Ncm.SBesselOdeSolver.new(lmin, a, b)

        # Solve using batched method
        solutions_batched = solver.solve_batched(rhs, lmin, n_l)
        solutions_batched_np = solutions_batched.to_numpy()

        # Solve individually for each l
        solutions_individual = []
        for ell in range(lmin, lmax + 1):
            solver.set_l(ell)
            sol = solver.solve(rhs)
            solutions_individual.append(sol.to_numpy())

        solutions_individual_np = np.array(solutions_individual)

        # Compare each solution
        for i, ell in enumerate(range(lmin, lmax + 1)):
            batched_sol = solutions_batched_np[i, :]
            individual_sol = solutions_individual_np[i, :]

            assert_allclose(
                batched_sol,
                individual_sol,
                rtol=1.0e-12,
                atol=1.0e-14,
                err_msg=f"Batched solution for l={ell} doesn't match individual solve",
            )

    def test_integrate_l_range_consistency(self) -> None:
        """Test that integrate_l_range gives consistent results.

        Compares integrate_l_range (which uses batched solver internally)
        against individual integrate calls for each l.
        """
        N = 128
        a, b = 1.0, 20.0
        lmin, lmax = 5, 15

        # Create solver
        solver = Ncm.SBesselOdeSolver.new(lmin, a, b)

        # Define a test function
        def f_test(x: float) -> float:
            return x * x / (1.0 + x * x)

        # Use integrate_l_range (batched internally)
        results_range = solver.integrate_l_range(f_test, N, lmin, lmax)
        results_range_np = results_range.to_numpy()

        # Compute individually for each l
        results_individual = []
        for ell in range(lmin, lmax + 1):
            solver.set_l(ell)
            result = solver.integrate(f_test, N)
            results_individual.append(result)

        results_individual_np = np.array(results_individual)

        # Compare
        assert_allclose(
            results_range_np,
            results_individual_np,
            rtol=1.0e-10,
            atol=1.0e-14,
            err_msg="integrate_l_range doesn't match individual integrate calls",
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

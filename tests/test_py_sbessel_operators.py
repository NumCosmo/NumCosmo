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

import numpy as np
import pytest
from numpy.testing import assert_allclose
from scipy.special import spherical_jn
from scipy.linalg import solve
from scipy.integrate import quad
from pathlib import Path
import gzip
import json

from numcosmo_py import Ncm


class TestSBesselOperators:
    """Tests for NcmSBesselOdeSolver operator matrices."""

    @staticmethod
    def get_analytical_chebyshev_coeffs(power: int, N: int) -> np.ndarray:
        """
        Get analytical Chebyshev coefficients for x^power.

        Uses the recurrence relation for x^n in terms of Chebyshev polynomials.
        """
        coeffs = np.zeros(N)

        if power == 0:
            coeffs[0] = 1.0
        elif power == 1:
            coeffs[1] = 1.0
        elif power == 2:
            coeffs[0] = 1.0 / 2.0
            coeffs[2] = 1.0 / 2.0
        elif power == 3:
            coeffs[1] = 3.0 / 4.0
            coeffs[3] = 1.0 / 4.0
        elif power == 4:
            coeffs[0] = 3.0 / 8.0
            coeffs[2] = 1.0 / 2.0
            coeffs[4] = 1.0 / 8.0
        elif power == 5:
            coeffs[1] = 5.0 / 8.0
            coeffs[3] = 5.0 / 16.0
            coeffs[5] = 1.0 / 16.0
        elif power == 6:
            coeffs[0] = 5.0 / 16.0
            coeffs[2] = 15.0 / 32.0
            coeffs[4] = 6.0 / 32.0
            coeffs[6] = 1.0 / 32.0

        return coeffs

    @staticmethod
    def matrix_to_numpy(ncm_mat: Ncm.Matrix) -> np.ndarray:
        """Convert NcmMatrix to numpy array."""
        nrows = ncm_mat.nrows()
        ncols = ncm_mat.ncols()
        return np.array(ncm_mat.dup_array()).reshape(nrows, ncols)

    @staticmethod
    def vector_to_numpy(ncm_vec: Ncm.Vector) -> np.ndarray:
        """Convert NcmVector to numpy array."""
        return np.array(ncm_vec.dup_array())

    def test_proj_matrix_identity(self) -> None:
        """Test that projection matrix correctly converts x^2 to Gegenbauer."""
        N = 32

        # Get Chebyshev coefficients for x^2
        cheb_coeffs = self.get_analytical_chebyshev_coeffs(2, N)
        cheb_vec = Ncm.Vector.new_array(cheb_coeffs.tolist())

        # Get projection matrix
        proj_mat = Ncm.SBesselOdeSolver.get_proj_matrix(N)
        proj_np = self.matrix_to_numpy(proj_mat)

        # Apply matrix: result should be Gegenbauer coeffs of x^2
        gegen_coeffs = proj_np @ cheb_coeffs

        # Verify by converting using the function
        gegen_vec_expected = Ncm.Vector.new(N)
        Ncm.SBesselOdeSolver.chebT_to_gegenbauer_lambda2(cheb_vec, gegen_vec_expected)
        gegen_expected = self.vector_to_numpy(gegen_vec_expected)

        assert_allclose(
            gegen_coeffs,
            gegen_expected,
            rtol=1.0e-12,
            atol=1.0e-14,
            err_msg="Projection matrix doesn't match conversion function",
        )

    def test_x_operator_on_x(self) -> None:
        """Test that x operator applied to x gives x^2."""
        N = 32

        # Get Chebyshev coefficients for x (input)
        cheb_x = self.get_analytical_chebyshev_coeffs(1, N)

        # Get x operator matrix
        x_mat = Ncm.SBesselOdeSolver.get_x_matrix(N)
        x_np = self.matrix_to_numpy(x_mat)

        # Apply: x * x = x^2
        gegen_result = x_np @ cheb_x

        # Expected: Gegenbauer coefficients of x^2
        cheb_x2 = self.get_analytical_chebyshev_coeffs(2, N)
        cheb_x2_vec = Ncm.Vector.new_array(cheb_x2.tolist())
        gegen_x2_vec = Ncm.Vector.new(N)
        Ncm.SBesselOdeSolver.chebT_to_gegenbauer_lambda2(cheb_x2_vec, gegen_x2_vec)
        gegen_expected = self.vector_to_numpy(gegen_x2_vec)

        assert_allclose(
            gegen_result,
            gegen_expected,
            rtol=1.0e-12,
            atol=1.0e-14,
            err_msg="x operator on x doesn't give x^2",
        )

    def test_x_operator_on_x2(self) -> None:
        """Test that x operator applied to x^2 gives x^3."""
        N = 32

        # Get Chebyshev coefficients for x^2 (input)
        cheb_x2 = self.get_analytical_chebyshev_coeffs(2, N)

        # Get x operator matrix
        x_mat = Ncm.SBesselOdeSolver.get_x_matrix(N)
        x_np = self.matrix_to_numpy(x_mat)

        # Apply: x * x^2 = x^3
        gegen_result = x_np @ cheb_x2

        # Expected: Gegenbauer coefficients of x^3
        cheb_x3 = self.get_analytical_chebyshev_coeffs(3, N)
        cheb_x3_vec = Ncm.Vector.new_array(cheb_x3.tolist())
        gegen_x3_vec = Ncm.Vector.new(N)
        Ncm.SBesselOdeSolver.chebT_to_gegenbauer_lambda2(cheb_x3_vec, gegen_x3_vec)
        gegen_expected = self.vector_to_numpy(gegen_x3_vec)

        assert_allclose(
            gegen_result,
            gegen_expected,
            rtol=1.0e-12,
            atol=1.0e-14,
            err_msg="x operator on x^2 doesn't give x^3",
        )

    def test_x2_operator_on_x(self) -> None:
        """Test that x^2 operator applied to x gives x^3."""
        N = 32

        # Get Chebyshev coefficients for x (input)
        cheb_x = self.get_analytical_chebyshev_coeffs(1, N)

        # Get x^2 operator matrix
        x2_mat = Ncm.SBesselOdeSolver.get_x2_matrix(N)
        x2_np = self.matrix_to_numpy(x2_mat)

        # Apply: x^2 * x = x^3
        gegen_result = x2_np @ cheb_x

        # Expected: Gegenbauer coefficients of x^3
        cheb_x3 = self.get_analytical_chebyshev_coeffs(3, N)
        cheb_x3_vec = Ncm.Vector.new_array(cheb_x3.tolist())
        gegen_x3_vec = Ncm.Vector.new(N)
        Ncm.SBesselOdeSolver.chebT_to_gegenbauer_lambda2(cheb_x3_vec, gegen_x3_vec)
        gegen_expected = self.vector_to_numpy(gegen_x3_vec)

        assert_allclose(
            gegen_result,
            gegen_expected,
            rtol=1.0e-12,
            atol=1.0e-14,
            err_msg="x^2 operator on x doesn't give x^3",
        )

    def test_x2_operator_on_x2(self) -> None:
        """Test that x^2 operator applied to x^2 gives x^4."""
        N = 32

        # Get Chebyshev coefficients for x^2 (input)
        cheb_x2 = self.get_analytical_chebyshev_coeffs(2, N)

        # Get x^2 operator matrix
        x2_mat = Ncm.SBesselOdeSolver.get_x2_matrix(N)
        x2_np = self.matrix_to_numpy(x2_mat)

        # Apply: x^2 * x^2 = x^4
        gegen_result = x2_np @ cheb_x2

        # Expected: Gegenbauer coefficients of x^4
        cheb_x4 = self.get_analytical_chebyshev_coeffs(4, N)
        cheb_x4_vec = Ncm.Vector.new_array(cheb_x4.tolist())
        gegen_x4_vec = Ncm.Vector.new(N)
        Ncm.SBesselOdeSolver.chebT_to_gegenbauer_lambda2(cheb_x4_vec, gegen_x4_vec)
        gegen_expected = self.vector_to_numpy(gegen_x4_vec)

        assert_allclose(
            gegen_result,
            gegen_expected,
            rtol=1.0e-12,
            atol=1.0e-14,
            err_msg="x^2 operator on x^2 doesn't give x^4",
        )

    def test_d_operator_on_x2(self) -> None:
        """Test that derivative operator applied to x^2 gives 2x."""
        N = 32

        # Get Chebyshev coefficients for x^2
        cheb_x2 = self.get_analytical_chebyshev_coeffs(2, N)

        # Get derivative operator matrix
        d_mat = Ncm.SBesselOdeSolver.get_d_matrix(N)
        d_np = self.matrix_to_numpy(d_mat)

        # Apply: d/dx(x^2) = 2x
        gegen_result = d_np @ cheb_x2

        # Expected: Gegenbauer coefficients of 2x
        cheb_2x = 2.0 * self.get_analytical_chebyshev_coeffs(1, N)
        cheb_2x_vec = Ncm.Vector.new_array(cheb_2x.tolist())
        gegen_2x_vec = Ncm.Vector.new(N)
        Ncm.SBesselOdeSolver.chebT_to_gegenbauer_lambda2(cheb_2x_vec, gegen_2x_vec)
        gegen_expected = self.vector_to_numpy(gegen_2x_vec)

        assert_allclose(
            gegen_result,
            gegen_expected,
            rtol=1.0e-12,
            atol=1.0e-14,
            err_msg="Derivative of x^2 doesn't give 2x",
        )

    def test_d_operator_on_x3(self) -> None:
        """Test that derivative operator applied to x^3 gives 3x^2."""
        N = 32

        # Get Chebyshev coefficients for x^3
        cheb_x3 = self.get_analytical_chebyshev_coeffs(3, N)

        # Get derivative operator matrix
        d_mat = Ncm.SBesselOdeSolver.get_d_matrix(N)
        d_np = self.matrix_to_numpy(d_mat)

        # Apply: d/dx(x^3) = 3x^2
        gegen_result = d_np @ cheb_x3

        # Expected: Gegenbauer coefficients of 3x^2
        cheb_3x2 = 3.0 * self.get_analytical_chebyshev_coeffs(2, N)
        cheb_3x2_vec = Ncm.Vector.new_array(cheb_3x2.tolist())
        gegen_3x2_vec = Ncm.Vector.new(N)
        Ncm.SBesselOdeSolver.chebT_to_gegenbauer_lambda2(cheb_3x2_vec, gegen_3x2_vec)
        gegen_expected = self.vector_to_numpy(gegen_3x2_vec)

        assert_allclose(
            gegen_result,
            gegen_expected,
            rtol=1.0e-12,
            atol=1.0e-14,
            err_msg="Derivative of x^3 doesn't give 3x^2",
        )

    def test_d2_operator_on_x3(self) -> None:
        """Test that second derivative operator applied to x^3 gives 6x."""
        N = 32

        # Get Chebyshev coefficients for x^3
        cheb_x3 = self.get_analytical_chebyshev_coeffs(3, N)

        # Get second derivative operator matrix
        d2_mat = Ncm.SBesselOdeSolver.get_d2_matrix(N)
        d2_np = self.matrix_to_numpy(d2_mat)

        # Apply: d^2/dx^2(x^3) = 6x
        gegen_result = d2_np @ cheb_x3

        # Expected: Gegenbauer coefficients of 6x
        cheb_6x = 6.0 * self.get_analytical_chebyshev_coeffs(1, N)
        cheb_6x_vec = Ncm.Vector.new_array(cheb_6x.tolist())
        gegen_6x_vec = Ncm.Vector.new(N)
        Ncm.SBesselOdeSolver.chebT_to_gegenbauer_lambda2(cheb_6x_vec, gegen_6x_vec)
        gegen_expected = self.vector_to_numpy(gegen_6x_vec)

        assert_allclose(
            gegen_result,
            gegen_expected,
            rtol=1.0e-12,
            atol=1.0e-14,
            err_msg="Second derivative of x^3 doesn't give 6x",
        )

    def test_d2_operator_on_x4(self) -> None:
        """Test that second derivative operator applied to x^4 gives 12x^2."""
        N = 32

        # Get Chebyshev coefficients for x^4
        cheb_x4 = self.get_analytical_chebyshev_coeffs(4, N)

        # Get second derivative operator matrix
        d2_mat = Ncm.SBesselOdeSolver.get_d2_matrix(N)
        d2_np = self.matrix_to_numpy(d2_mat)

        # Apply: d^2/dx^2(x^4) = 12x^2
        gegen_result = d2_np @ cheb_x4

        # Expected: Gegenbauer coefficients of 12x^2
        cheb_12x2 = 12.0 * self.get_analytical_chebyshev_coeffs(2, N)
        cheb_12x2_vec = Ncm.Vector.new_array(cheb_12x2.tolist())
        gegen_12x2_vec = Ncm.Vector.new(N)
        Ncm.SBesselOdeSolver.chebT_to_gegenbauer_lambda2(cheb_12x2_vec, gegen_12x2_vec)
        gegen_expected = self.vector_to_numpy(gegen_12x2_vec)

        assert_allclose(
            gegen_result,
            gegen_expected,
            rtol=1.0e-12,
            atol=1.0e-14,
            err_msg="Second derivative of x^4 doesn't give 12x^2",
        )

    def test_x_d_operator_on_x3(self) -> None:
        """Test that x*d/dx operator applied to x^3 gives 3x^3."""
        N = 32

        # Get Chebyshev coefficients for x^3
        cheb_x3 = self.get_analytical_chebyshev_coeffs(3, N)

        # Get x*d operator matrix
        x_d_mat = Ncm.SBesselOdeSolver.get_x_d_matrix(N)
        x_d_np = self.matrix_to_numpy(x_d_mat)

        # Apply: x * d/dx(x^3) = x * 3x^2 = 3x^3
        gegen_result = x_d_np @ cheb_x3

        # Expected: Gegenbauer coefficients of 3x^3
        cheb_3x3 = 3.0 * self.get_analytical_chebyshev_coeffs(3, N)
        cheb_3x3_vec = Ncm.Vector.new_array(cheb_3x3.tolist())
        gegen_3x3_vec = Ncm.Vector.new(N)
        Ncm.SBesselOdeSolver.chebT_to_gegenbauer_lambda2(cheb_3x3_vec, gegen_3x3_vec)
        gegen_expected = self.vector_to_numpy(gegen_3x3_vec)

        assert_allclose(
            gegen_result,
            gegen_expected,
            rtol=1.0e-11,
            atol=1.0e-13,
            err_msg="x*d operator on x^3 doesn't give 3x^3",
        )

    def test_x_d2_operator_on_x4(self) -> None:
        """Test that x*d^2/dx^2 operator applied to x^4 gives 12x^3."""
        N = 32

        # Get Chebyshev coefficients for x^4
        cheb_x4 = self.get_analytical_chebyshev_coeffs(4, N)

        # Get x*d^2 operator matrix
        x_d2_mat = Ncm.SBesselOdeSolver.get_x_d2_matrix(N)
        x_d2_np = self.matrix_to_numpy(x_d2_mat)

        # Apply: x * d^2/dx^2(x^4) = x * 12x^2 = 12x^3
        gegen_result = x_d2_np @ cheb_x4

        # Expected: Gegenbauer coefficients of 12x^3
        cheb_12x3 = 12.0 * self.get_analytical_chebyshev_coeffs(3, N)
        cheb_12x3_vec = Ncm.Vector.new_array(cheb_12x3.tolist())
        gegen_12x3_vec = Ncm.Vector.new(N)
        Ncm.SBesselOdeSolver.chebT_to_gegenbauer_lambda2(cheb_12x3_vec, gegen_12x3_vec)
        gegen_expected = self.vector_to_numpy(gegen_12x3_vec)

        assert_allclose(
            gegen_result,
            gegen_expected,
            rtol=1.0e-11,
            atol=1.0e-13,
            err_msg="x*d^2 operator on x^4 doesn't give 12x^3",
        )

    def test_x2_d2_operator_on_x4(self) -> None:
        """Test that x^2*d^2/dx^2 operator applied to x^4 gives 12x^4."""
        N = 32

        # Get Chebyshev coefficients for x^4
        cheb_x4 = self.get_analytical_chebyshev_coeffs(4, N)

        # Get x^2*d^2 operator matrix
        x2_d2_mat = Ncm.SBesselOdeSolver.get_x2_d2_matrix(N)
        x2_d2_np = self.matrix_to_numpy(x2_d2_mat)

        # Apply: x^2 * d^2/dx^2(x^4) = x^2 * 12x^2 = 12x^4
        gegen_result = x2_d2_np @ cheb_x4

        # Expected: Gegenbauer coefficients of 12x^4
        cheb_12x4 = 12.0 * self.get_analytical_chebyshev_coeffs(4, N)
        cheb_12x4_vec = Ncm.Vector.new_array(cheb_12x4.tolist())
        gegen_12x4_vec = Ncm.Vector.new(N)
        Ncm.SBesselOdeSolver.chebT_to_gegenbauer_lambda2(cheb_12x4_vec, gegen_12x4_vec)
        gegen_expected = self.vector_to_numpy(gegen_12x4_vec)

        assert_allclose(
            gegen_result,
            gegen_expected,
            rtol=1.0e-11,
            atol=1.0e-13,
            err_msg="x^2*d^2 operator on x^4 doesn't give 12x^4",
        )

    def test_d_operator_evaluation(self) -> None:
        """Test d/dx operator by evaluating result at points."""
        N = 32

        # Start with x^3
        cheb_x3 = self.get_analytical_chebyshev_coeffs(3, N)

        # Apply d operator: d/dx(x^3) = 3x^2
        d_mat = Ncm.SBesselOdeSolver.get_d_matrix(N)
        d_np = self.matrix_to_numpy(d_mat)
        gegen_result = d_np @ cheb_x3
        gegen_result_vec = Ncm.Vector.new_array(gegen_result.tolist())

        # Evaluate at test points and compare with 3*x^2
        test_points = np.array([0.0, 0.3, 0.5, 0.7, 0.9])
        for x in test_points:
            eval_result = Ncm.SBesselOdeSolver.gegenbauer_lambda2_eval(
                gegen_result_vec, x
            )
            expected = 3.0 * x**2
            assert_allclose(
                eval_result,
                expected,
                rtol=1.0e-12,
                atol=1.0e-14,
                err_msg=f"d/dx(x^3) evaluation at x={x} doesn't match 3x^2",
            )

    def test_d2_operator_evaluation(self) -> None:
        """Test d^2/dx^2 operator by evaluating result at points."""
        N = 32

        # Start with x^4
        cheb_x4 = self.get_analytical_chebyshev_coeffs(4, N)

        # Apply d^2 operator: d^2/dx^2(x^4) = 12x^2
        d2_mat = Ncm.SBesselOdeSolver.get_d2_matrix(N)
        d2_np = self.matrix_to_numpy(d2_mat)
        gegen_result = d2_np @ cheb_x4
        gegen_result_vec = Ncm.Vector.new_array(gegen_result.tolist())

        # Evaluate at test points and compare with 12*x^2
        test_points = np.array([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
        for x in test_points:
            eval_result = Ncm.SBesselOdeSolver.gegenbauer_lambda2_eval(
                gegen_result_vec, x
            )
            expected = 12.0 * x**2
            assert_allclose(
                eval_result,
                expected,
                rtol=1.0e-12,
                atol=1.0e-14,
                err_msg=f"d^2/dx^2(x^4) evaluation at x={x} doesn't match 12x^2",
            )

    def test_x_d_operator_evaluation(self) -> None:
        """Test x*d/dx operator by evaluating result at points."""
        N = 32

        # Start with x^4
        cheb_x4 = self.get_analytical_chebyshev_coeffs(4, N)

        # Apply x*d operator: x*d/dx(x^4) = x*4x^3 = 4x^4
        x_d_mat = Ncm.SBesselOdeSolver.get_x_d_matrix(N)
        x_d_np = self.matrix_to_numpy(x_d_mat)
        gegen_result = x_d_np @ cheb_x4
        gegen_result_vec = Ncm.Vector.new_array(gegen_result.tolist())

        # Evaluate at test points and compare with 4*x^4
        test_points = np.array([0.1, 0.3, 0.5, 0.7, 0.9])
        for x in test_points:
            eval_result = Ncm.SBesselOdeSolver.gegenbauer_lambda2_eval(
                gegen_result_vec, x
            )
            expected = 4.0 * x**4
            assert_allclose(
                eval_result,
                expected,
                rtol=1.0e-11,
                atol=1.0e-13,
                err_msg=f"x*d/dx(x^4) evaluation at x={x} doesn't match 4x^4",
            )

    def test_x2_d2_operator_evaluation(self) -> None:
        """Test x^2*d^2/dx^2 operator by evaluating result at points."""
        N = 32

        # Start with x^5
        cheb_x5 = self.get_analytical_chebyshev_coeffs(5, N)

        # Apply x^2*d^2 operator: x^2*d^2/dx^2(x^5) = x^2*20x^3 = 20x^5
        x2_d2_mat = Ncm.SBesselOdeSolver.get_x2_d2_matrix(N)
        x2_d2_np = self.matrix_to_numpy(x2_d2_mat)
        gegen_result = x2_d2_np @ cheb_x5
        gegen_result_vec = Ncm.Vector.new_array(gegen_result.tolist())

        # Evaluate at test points and compare with 20*x^5
        test_points = np.array([0.1, 0.3, 0.5, 0.7, 0.9])
        for x in test_points:
            eval_result = Ncm.SBesselOdeSolver.gegenbauer_lambda2_eval(
                gegen_result_vec, x
            )
            expected = 20.0 * x**5
            assert_allclose(
                eval_result,
                expected,
                rtol=1.0e-11,
                atol=1.0e-13,
                err_msg=f"x^2*d^2/dx^2(x^5) evaluation at x={x} doesn't match 20x^5",
            )

    def test_x_operator_evaluation(self) -> None:
        """Test x operator by evaluating result at points."""
        N = 32

        # Start with x^3
        cheb_x3 = self.get_analytical_chebyshev_coeffs(3, N)

        # Apply x operator: x*x^3 = x^4
        x_mat = Ncm.SBesselOdeSolver.get_x_matrix(N)
        x_np = self.matrix_to_numpy(x_mat)
        gegen_result = x_np @ cheb_x3
        gegen_result_vec = Ncm.Vector.new_array(gegen_result.tolist())

        # Evaluate at test points and compare with x^4
        test_points = np.array([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
        for x in test_points:
            eval_result = Ncm.SBesselOdeSolver.gegenbauer_lambda2_eval(
                gegen_result_vec, x
            )
            expected = x**4
            assert_allclose(
                eval_result,
                expected,
                rtol=1.0e-12,
                atol=1.0e-14,
                err_msg=f"x*x^3 evaluation at x={x} doesn't match x^4",
            )

    def test_x2_operator_evaluation(self) -> None:
        """Test x^2 operator by evaluating result at points."""
        N = 32

        # Start with x^3
        cheb_x3 = self.get_analytical_chebyshev_coeffs(3, N)

        # Apply x^2 operator: x^2*x^3 = x^5
        x2_mat = Ncm.SBesselOdeSolver.get_x2_matrix(N)
        x2_np = self.matrix_to_numpy(x2_mat)
        gegen_result = x2_np @ cheb_x3
        gegen_result_vec = Ncm.Vector.new_array(gegen_result.tolist())

        # Evaluate at test points and compare with x^5
        test_points = np.array([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
        for x in test_points:
            eval_result = Ncm.SBesselOdeSolver.gegenbauer_lambda2_eval(
                gegen_result_vec, x
            )
            expected = x**5
            assert_allclose(
                eval_result,
                expected,
                rtol=1.0e-12,
                atol=1.0e-14,
                err_msg=f"x^2*x^3 evaluation at x={x} doesn't match x^5",
            )

    def test_gaussian_operators(self) -> None:
        """Test operators on exp(-x^2), a non-polynomial function."""
        N = 64

        solver = Ncm.SBesselOdeSolver.new(0, -1.0, 1.0)

        # Define f(x) = exp(-x^2)
        def f_gaussian(_user_data: None, x: float) -> float:
            return np.exp(-(x**2))

        # Compute Chebyshev coefficients for exp(-x^2)
        cheb_vec = solver.compute_chebyshev_coeffs(f_gaussian, -1.0, 1.0, N, None)
        cheb_coeffs = self.vector_to_numpy(cheb_vec)

        # Get operator matrices
        proj_mat = Ncm.SBesselOdeSolver.get_proj_matrix(N)
        x_mat = Ncm.SBesselOdeSolver.get_x_matrix(N)
        x2_mat = Ncm.SBesselOdeSolver.get_x2_matrix(N)
        d_mat = Ncm.SBesselOdeSolver.get_d_matrix(N)
        x_d_mat = Ncm.SBesselOdeSolver.get_x_d_matrix(N)
        d2_mat = Ncm.SBesselOdeSolver.get_d2_matrix(N)
        x_d2_mat = Ncm.SBesselOdeSolver.get_x_d2_matrix(N)
        x2_d2_mat = Ncm.SBesselOdeSolver.get_x2_d2_matrix(N)

        # Convert to numpy
        proj_np = self.matrix_to_numpy(proj_mat)
        x_np = self.matrix_to_numpy(x_mat)
        x2_np = self.matrix_to_numpy(x2_mat)
        d_np = self.matrix_to_numpy(d_mat)
        x_d_np = self.matrix_to_numpy(x_d_mat)
        d2_np = self.matrix_to_numpy(d2_mat)
        x_d2_np = self.matrix_to_numpy(x_d2_mat)
        x2_d2_np = self.matrix_to_numpy(x2_d2_mat)

        # Test points for evaluation
        test_points = np.array([0.0, 0.3, 0.5, 0.7, 0.9])

        # Test 1: Projection gives same function
        gegen_proj = proj_np @ cheb_coeffs
        gegen_proj_vec = Ncm.Vector.new_array(gegen_proj.tolist())
        for x in test_points:
            result = Ncm.SBesselOdeSolver.gegenbauer_lambda2_eval(gegen_proj_vec, x)
            expected = np.exp(-(x**2))
            assert_allclose(
                result,
                expected,
                rtol=1.0e-8,
                atol=0.0,
                err_msg=f"Projection of exp(-x^2) at x={x} failed",
            )

        # Test 2: x operator: x*exp(-x^2)
        gegen_x = x_np @ cheb_coeffs
        gegen_x_vec = Ncm.Vector.new_array(gegen_x.tolist())
        for x in test_points:
            result = Ncm.SBesselOdeSolver.gegenbauer_lambda2_eval(gegen_x_vec, x)
            expected = x * np.exp(-(x**2))
            assert_allclose(
                result,
                expected,
                rtol=1.0e-8,
                atol=1.0e-15,
                err_msg=f"x*exp(-x^2) at x={x} failed",
            )

        # Test 3: x^2 operator: x^2*exp(-x^2)
        gegen_x2 = x2_np @ cheb_coeffs
        gegen_x2_vec = Ncm.Vector.new_array(gegen_x2.tolist())
        for x in test_points:
            result = Ncm.SBesselOdeSolver.gegenbauer_lambda2_eval(gegen_x2_vec, x)
            expected = x**2 * np.exp(-(x**2))
            assert_allclose(
                result,
                expected,
                rtol=1.0e-8,
                atol=1.0e-15,
                err_msg=f"x^2*exp(-x^2) at x={x} failed",
            )

        # Test 4: d operator: d/dx(exp(-x^2)) = -2x*exp(-x^2)
        gegen_d = d_np @ cheb_coeffs
        gegen_d_vec = Ncm.Vector.new_array(gegen_d.tolist())
        for x in test_points:
            result = Ncm.SBesselOdeSolver.gegenbauer_lambda2_eval(gegen_d_vec, x)
            expected = -2.0 * x * np.exp(-(x**2))
            assert_allclose(
                result,
                expected,
                rtol=1.0e-7,
                atol=1.0e-15,
                err_msg=f"d/dx(exp(-x^2)) at x={x} failed",
            )

        # Test 5: d^2 operator: d^2/dx^2(exp(-x^2)) = (4x^2 - 2)*exp(-x^2)
        gegen_d2 = d2_np @ cheb_coeffs
        gegen_d2_vec = Ncm.Vector.new_array(gegen_d2.tolist())
        for x in test_points:
            result = Ncm.SBesselOdeSolver.gegenbauer_lambda2_eval(gegen_d2_vec, x)
            expected = (4.0 * x**2 - 2.0) * np.exp(-(x**2))
            assert_allclose(
                result,
                expected,
                rtol=1.0e-6,
                atol=1.0e-15,
                err_msg=f"d^2/dx^2(exp(-x^2)) at x={x} failed",
            )

        # Test 6: x*d operator: x*d/dx(exp(-x^2)) = x*(-2x*exp(-x^2)) = -2x^2*exp(-x^2)
        gegen_x_d = x_d_np @ cheb_coeffs
        gegen_x_d_vec = Ncm.Vector.new_array(gegen_x_d.tolist())
        for x in test_points:
            result = Ncm.SBesselOdeSolver.gegenbauer_lambda2_eval(gegen_x_d_vec, x)
            expected = -2.0 * x**2 * np.exp(-(x**2))
            assert_allclose(
                result,
                expected,
                rtol=1.0e-7,
                atol=1.0e-9,
                err_msg=f"x*d/dx(exp(-x^2)) at x={x} failed",
            )

        # Test 7: x*d^2 operator: x*d^2/dx^2(exp(-x^2)) = x*(4x^2 - 2)*exp(-x^2)
        gegen_x_d2 = x_d2_np @ cheb_coeffs
        gegen_x_d2_vec = Ncm.Vector.new_array(gegen_x_d2.tolist())
        for x in test_points:
            result = Ncm.SBesselOdeSolver.gegenbauer_lambda2_eval(gegen_x_d2_vec, x)
            expected = x * (4.0 * x**2 - 2.0) * np.exp(-(x**2))
            assert_allclose(
                result,
                expected,
                rtol=1.0e-6,
                atol=1.0e-8,
                err_msg=f"x*d^2/dx^2(exp(-x^2)) at x={x} failed",
            )

        # Test 8: x^2*d^2 operator: x^2*d^2/dx^2(exp(-x^2)) = x^2*(4x^2 - 2)*exp(-x^2)
        gegen_x2_d2 = x2_d2_np @ cheb_coeffs
        gegen_x2_d2_vec = Ncm.Vector.new_array(gegen_x2_d2.tolist())
        for x in test_points:
            result = Ncm.SBesselOdeSolver.gegenbauer_lambda2_eval(gegen_x2_d2_vec, x)
            expected = x**2 * (4.0 * x**2 - 2.0) * np.exp(-(x**2))
            assert_allclose(
                result,
                expected,
                rtol=1.0e-6,
                atol=1.0e-8,
                err_msg=f"x^2*d^2/dx^2(exp(-x^2)) at x={x} failed",
            )

    def test_rational_operators(self) -> None:
        """Test operators on x^2/(1+x^2)^4, a rational function."""
        N = 64

        solver = Ncm.SBesselOdeSolver.new(0, -1.0, 1.0)

        # Define f(x) = x^2/(1+x^2)^4
        def f_rational(_user_data: None, x: float) -> float:
            return x**2 / (1.0 + x**2) ** 4

        # Analytical derivative: f'(x) = 2x(1-3x^2)/(1+x^2)^5
        def f_prime(x: float) -> float:
            return 2.0 * x * (1.0 - 3.0 * x**2) / (1.0 + x**2) ** 5

        # Analytical second derivative: f''(x) = (2 - 36x^2 + 42x^4)/(1+x^2)^6
        def f_double_prime(x: float) -> float:
            return (2.0 - 36.0 * x**2 + 42.0 * x**4) / (1.0 + x**2) ** 6

        # Compute Chebyshev coefficients
        cheb_vec = solver.compute_chebyshev_coeffs(f_rational, -1.0, 1.0, N, None)
        cheb_coeffs = self.vector_to_numpy(cheb_vec)

        # Get operator matrices
        proj_mat = Ncm.SBesselOdeSolver.get_proj_matrix(N)
        x_mat = Ncm.SBesselOdeSolver.get_x_matrix(N)
        x2_mat = Ncm.SBesselOdeSolver.get_x2_matrix(N)
        d_mat = Ncm.SBesselOdeSolver.get_d_matrix(N)
        x_d_mat = Ncm.SBesselOdeSolver.get_x_d_matrix(N)
        d2_mat = Ncm.SBesselOdeSolver.get_d2_matrix(N)
        x_d2_mat = Ncm.SBesselOdeSolver.get_x_d2_matrix(N)
        x2_d2_mat = Ncm.SBesselOdeSolver.get_x2_d2_matrix(N)

        # Convert to numpy
        proj_np = self.matrix_to_numpy(proj_mat)
        x_np = self.matrix_to_numpy(x_mat)
        x2_np = self.matrix_to_numpy(x2_mat)
        d_np = self.matrix_to_numpy(d_mat)
        x_d_np = self.matrix_to_numpy(x_d_mat)
        d2_np = self.matrix_to_numpy(d2_mat)
        x_d2_np = self.matrix_to_numpy(x_d2_mat)
        x2_d2_np = self.matrix_to_numpy(x2_d2_mat)

        # Test points for evaluation
        test_points = np.array([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])

        # Test 1: Projection gives same function
        gegen_proj = proj_np @ cheb_coeffs
        gegen_proj_vec = Ncm.Vector.new_array(gegen_proj.tolist())
        for x in test_points:
            result = Ncm.SBesselOdeSolver.gegenbauer_lambda2_eval(gegen_proj_vec, x)
            expected = x**2 / (1.0 + x**2) ** 4
            assert_allclose(
                result,
                expected,
                rtol=1.0e-8,
                atol=1.0e-10,
                err_msg=f"Projection of x^2/(1+x^2)^4 at x={x} failed",
            )

        # Test 2: x operator: x * x^2/(1+x^2)^4 = x^3/(1+x^2)^4
        gegen_x = x_np @ cheb_coeffs
        gegen_x_vec = Ncm.Vector.new_array(gegen_x.tolist())
        for x in test_points:
            result = Ncm.SBesselOdeSolver.gegenbauer_lambda2_eval(gegen_x_vec, x)
            expected = x**3 / (1.0 + x**2) ** 4
            assert_allclose(
                result,
                expected,
                rtol=1.0e-8,
                atol=1.0e-10,
                err_msg=f"x*f(x) at x={x} failed",
            )

        # Test 3: x^2 operator: x^2 * x^2/(1+x^2)^4 = x^4/(1+x^2)^4
        gegen_x2 = x2_np @ cheb_coeffs
        gegen_x2_vec = Ncm.Vector.new_array(gegen_x2.tolist())
        for x in test_points:
            result = Ncm.SBesselOdeSolver.gegenbauer_lambda2_eval(gegen_x2_vec, x)
            expected = x**4 / (1.0 + x**2) ** 4
            assert_allclose(
                result,
                expected,
                rtol=1.0e-8,
                atol=1.0e-10,
                err_msg=f"x^2*f(x) at x={x} failed",
            )

        # Test 4: d operator: d/dx of f(x)
        gegen_d = d_np @ cheb_coeffs
        gegen_d_vec = Ncm.Vector.new_array(gegen_d.tolist())
        for x in test_points:
            result = Ncm.SBesselOdeSolver.gegenbauer_lambda2_eval(gegen_d_vec, x)
            expected = f_prime(x)
            assert_allclose(
                result,
                expected,
                rtol=1.0e-7,
                atol=1.0e-9,
                err_msg=f"d/dx(f(x)) at x={x} failed",
            )

        # Test 5: d^2 operator: d^2/dx^2 of f(x)
        gegen_d2 = d2_np @ cheb_coeffs
        gegen_d2_vec = Ncm.Vector.new_array(gegen_d2.tolist())
        for x in test_points:
            result = Ncm.SBesselOdeSolver.gegenbauer_lambda2_eval(gegen_d2_vec, x)
            expected = f_double_prime(x)
            assert_allclose(
                result,
                expected,
                rtol=1.0e-6,
                atol=1.0e-8,
                err_msg=f"d^2/dx^2(f(x)) at x={x} failed",
            )

        # Test 6: x*d operator: x*f'(x)
        gegen_x_d = x_d_np @ cheb_coeffs
        gegen_x_d_vec = Ncm.Vector.new_array(gegen_x_d.tolist())
        for x in test_points:
            result = Ncm.SBesselOdeSolver.gegenbauer_lambda2_eval(gegen_x_d_vec, x)
            expected = x * f_prime(x)
            assert_allclose(
                result,
                expected,
                rtol=1.0e-7,
                atol=1.0e-9,
                err_msg=f"x*d/dx(f(x)) at x={x} failed",
            )

        # Test 7: x*d^2 operator: x*f''(x)
        gegen_x_d2 = x_d2_np @ cheb_coeffs
        gegen_x_d2_vec = Ncm.Vector.new_array(gegen_x_d2.tolist())
        for x in test_points:
            result = Ncm.SBesselOdeSolver.gegenbauer_lambda2_eval(gegen_x_d2_vec, x)
            expected = x * f_double_prime(x)
            assert_allclose(
                result,
                expected,
                rtol=1.0e-6,
                atol=1.0e-8,
                err_msg=f"x*d^2/dx^2(f(x)) at x={x} failed",
            )

        # Test 8: x^2*d^2 operator: x^2*f''(x)
        gegen_x2_d2 = x2_d2_np @ cheb_coeffs
        gegen_x2_d2_vec = Ncm.Vector.new_array(gegen_x2_d2.tolist())
        for x in test_points:
            result = Ncm.SBesselOdeSolver.gegenbauer_lambda2_eval(gegen_x2_d2_vec, x)
            expected = x**2 * f_double_prime(x)
            assert_allclose(
                result,
                expected,
                rtol=1.0e-6,
                atol=1.0e-8,
                err_msg=f"x^2*d^2/dx^2(f(x)) at x={x} failed",
            )

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
        mat_np = self.matrix_to_numpy(mat)

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
            y_computed = Ncm.SBesselOdeSolver.chebyshev_eval(solution_coeffs, t)

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
        y_at_a = Ncm.SBesselOdeSolver.chebyshev_eval(solution_coeffs, -1.0)
        y_at_b = Ncm.SBesselOdeSolver.chebyshev_eval(solution_coeffs, 1.0)
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
        mat_np = self.matrix_to_numpy(mat)

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
        y_prime_at_minus1 = Ncm.SBesselOdeSolver.chebyshev_deriv(solution_coeffs, -1.0)
        y_prime_at_plus1 = Ncm.SBesselOdeSolver.chebyshev_deriv(solution_coeffs, 1.0)

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
        mat_np = self.matrix_to_numpy(mat)

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
        y_prime_at_minus1 = Ncm.SBesselOdeSolver.chebyshev_deriv(solution_coeffs, -1.0)
        y_prime_at_plus1 = Ncm.SBesselOdeSolver.chebyshev_deriv(solution_coeffs, 1.0)

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
        solution_np = self.vector_to_numpy(solution)

        op_mat = self.matrix_to_numpy(solver.get_operator_matrix(N))

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
        mat_np = self.matrix_to_numpy(mat)
        solution_scipy = solve(mat_np, rhs_np)

        # Solve using solve_dense (LAPACK dgesv)
        solution_dense = solver.solve_dense(rhs_vec, N)
        solution_dense_np = self.vector_to_numpy(solution_dense)

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
        y_prime_at_minus1 = Ncm.SBesselOdeSolver.chebyshev_deriv(solution_coeffs, -1.0)
        y_prime_at_plus1 = Ncm.SBesselOdeSolver.chebyshev_deriv(solution_coeffs, 1.0)
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
        y_prime_at_minus1 = Ncm.SBesselOdeSolver.chebyshev_deriv(solution_coeffs, -1.0)
        y_prime_at_plus1 = Ncm.SBesselOdeSolver.chebyshev_deriv(solution_coeffs, 1.0)
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
        y_prime_at_minus1 = Ncm.SBesselOdeSolver.chebyshev_deriv(solution_coeffs, -1.0)
        y_prime_at_plus1 = Ncm.SBesselOdeSolver.chebyshev_deriv(solution_coeffs, 1.0)
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
        y_prime_at_minus1 = Ncm.SBesselOdeSolver.chebyshev_deriv(solution_coeffs, -1.0)
        y_prime_at_plus1 = Ncm.SBesselOdeSolver.chebyshev_deriv(solution_coeffs, 1.0)
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
        solutions_batched_np = self.matrix_to_numpy(solutions_batched)

        # Solve individually for each l
        solutions_individual = []
        for ell in range(lmin, lmax + 1):
            solver.set_l(ell)
            sol = solver.solve(rhs)
            solutions_individual.append(self.vector_to_numpy(sol))

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
        results_range_np = self.vector_to_numpy(results_range)

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
        ell_max = 200  # Limit for faster testing

        N = 2**14  # Number of Chebyshev nodes
        print(f"Using N = {N} Chebyshev nodes")

        print_rank = False
        print_ell: list[int] | None = [50]
        solver = Ncm.SBesselOdeSolver.new(0, lb, ub)

        for i in range(1):
            print(f"Starting iteration {i}\r", end="", flush=True)
            for i, k in enumerate(truth_table["kvals"]):
                # if (i != 50) and i < 100:
                #    continue

                # if k != 100.0:
                #    continue

                # Create solver for this ell

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
                            f"[{func_type}] ell={ell:5d}, k={k: 22.15g}, "
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

#!/usr/bin/env python
# Test for NcmSBesselOdeSolver operator matrices


import numpy as np
from numpy.testing import assert_allclose
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

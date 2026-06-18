#!/usr/bin/env python
#
# test_py_sf_sbessel.py
#
# Thu Feb 06 00:00:00 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_py_sf_sbessel.py
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

"""Tests on NcmSFSBessel class and spherical Bessel functions."""

import numpy as np
import pytest

from numpy.testing import assert_allclose
from scipy.special import spherical_jn

from numcosmo_py import Ncm

Ncm.cfg_init()


class TestSFSBessel:
    """Tests for spherical Bessel function evaluations."""

    @pytest.mark.parametrize("l_val", [0, 1, 2, 5, 10, 20, 50, 100])
    @pytest.mark.parametrize("x_val", [0.1, 1.0, 5.0, 10.0, 50.0, 100.0])
    def test_ncm_sf_sbessel(self, l_val: int, x_val: float) -> None:
        """Test single spherical Bessel function evaluation against scipy.

        Compares ncm_sf_sbessel(l, x) with scipy.special.spherical_jn(l, x)
        for various l and x values.
        """
        result = Ncm.sf_sbessel(l_val, x_val)
        expected = spherical_jn(l_val, x_val)

        assert_allclose(
            result,
            expected,
            rtol=1.0e-12,
            atol=1.0e-15,
            err_msg=f"Mismatch for j_{l_val}({x_val})",
        )

    @pytest.mark.parametrize("l_val", [0, 1, 2, 5, 10])
    def test_ncm_sf_sbessel_at_zero(self, l_val: int) -> None:
        """Test spherical Bessel function at x=0.

        j_0(0) = 1, j_l(0) = 0 for l > 0.
        """
        result = Ncm.sf_sbessel(l_val, 0.0)
        expected = 1.0 if l_val == 0 else 0.0

        assert_allclose(
            result,
            expected,
            rtol=1.0e-15,
            atol=1.0e-15,
            err_msg=f"Mismatch for j_{l_val}(0)",
        )

    @pytest.mark.parametrize("l_val", [0, 1, 2, 5, 10])
    def test_ncm_sf_sbessel_spline(self, l_val: int) -> None:
        """Test spline interpolation of spherical Bessel function.

        Creates a spline approximation over [0.1, 10.0] and verifies
        it matches the exact function at test points.
        """
        xi, xf = 0.1, 10.0
        reltol = 1.0e-10

        # Create spline
        spline = Ncm.sf_sbessel_spline(l_val, xi, xf, reltol)

        # Test at several points
        x_test = np.linspace(xi, xf, 20)
        for x in x_test:
            spline_val = spline.eval(x)
            exact_val = spherical_jn(l_val, x)

            assert_allclose(
                spline_val,
                exact_val,
                rtol=reltol * 10.0,  # Allow some margin
                atol=1.0e-15,
                err_msg=f"Spline mismatch for j_{l_val}({x})",
            )


class TestSFSBesselArray:
    """Tests for NcmSFSBesselArray class."""

    def test_array_creation_default(self) -> None:
        """Test creating array with default parameters."""
        sba = Ncm.SFSBesselArray.new()
        assert sba.get_lmax() == 10000
        assert sba.get_threshold() == 1.0e-100

    def test_array_creation_custom(self) -> None:
        """Test creating array with custom parameters."""
        lmax = 500
        threshold = 1.0e-50
        sba = Ncm.SFSBesselArray.new_full(lmax, threshold)

        assert sba.get_lmax() == lmax
        assert sba.get_threshold() == threshold

    @pytest.mark.parametrize("lmax", [10, 50, 100])
    @pytest.mark.parametrize("x_val", [1.0, 5.0, 10.0, 50.0])
    def test_array_eval(self, lmax: int, x_val: float) -> None:
        """Test array evaluation against scipy.

        Evaluates j_l(x) for l=0 to lmax and compares with scipy results.
        """
        sba = Ncm.SFSBesselArray.new_full(lmax, 1.0e-100)

        jl_x = sba.eval1(lmax, x_val)

        # Compare with scipy
        for ell in range(lmax + 1):
            expected = spherical_jn(ell, x_val)

            # If expected is very small, check if result is also small
            if abs(expected) < 1.0e-100:
                assert abs(jl_x[ell]) < 1.0e-50, f"Value too large for j_{ell}({x_val})"
            else:
                assert_allclose(
                    jl_x[ell],
                    expected,
                    rtol=1.0e-10,
                    atol=1.0e-15,
                    err_msg=f"Array mismatch for j_{ell}({x_val})",
                )

    def test_array_eval_at_zero(self) -> None:
        """Test that j_0(x) = 1 for x=0 and j_l(0) = 0 for l > 0."""
        sba = Ncm.SFSBesselArray.new()

        # Test at x=0
        jl_0 = sba.eval1(10, 0.0)

        assert_allclose(jl_0[0], 1.0, rtol=1.0e-15, atol=1.0e-15)
        for ell in range(1, 11):
            assert_allclose(
                jl_0[ell], 0.0, rtol=1.0e-15, atol=1.0e-15, err_msg=f"j_{ell}(0) != 0"
            )

    def test_array_cutoff(self) -> None:
        """Test automatic cutoff logic.

        For large l and small x, j_l(x) becomes negligible.
        Verify that the cutoff is applied correctly.
        """
        lmax = 200
        threshold = 1.0e-50
        sba = Ncm.SFSBesselArray.new_full(lmax, threshold)

        x_val = 10.0

        # Get cutoff
        cutoff = sba.eval_ell_cutoff(x_val)

        # Verify cutoff is reasonable
        assert cutoff > 0
        assert cutoff <= lmax

        # Evaluate with full lmax
        jl_x = sba.eval1(lmax, x_val)

        # Values beyond cutoff should be near zero (below threshold)
        for ell in range(cutoff, lmax + 1):
            assert (
                abs(jl_x[ell]) <= threshold * 10
            ), f"j_{ell}({x_val}) should be below threshold beyond cutoff"

    @pytest.mark.parametrize(
        "x_val,expected_min",
        [
            (1.0, 2),  # Small x -> small cutoff
            (10.0, 10),  # Medium x -> medium cutoff
            (50.0, 50),  # Large x -> large cutoff
        ],
    )
    def test_cutoff_scaling(self, x_val: float, expected_min: int) -> None:
        """Test that cutoff scales reasonably with x.

        Larger x should allow larger l before cutoff.
        """
        sba = Ncm.SFSBesselArray.new()
        cutoff = sba.eval_ell_cutoff(x_val)

        # Check cutoff is at least expected minimum and roughly proportional to x
        assert (
            cutoff >= expected_min
        ), f"Cutoff {cutoff} should be at least {expected_min} for x={x_val}"

    def test_array_small_x_taylor(self) -> None:
        """Test array evaluation for very small x using Taylor series.

        For x < 2*DBL_EPSILON^(1/4), the implementation uses Taylor series.
        """
        sba = Ncm.SFSBesselArray.new()
        x_val = 1.0e-3  # Very small x
        lmax = 5

        jl_x = sba.eval1(lmax, x_val)

        # Compare with scipy
        for ell in range(lmax + 1):
            expected = spherical_jn(ell, x_val)
            assert_allclose(
                jl_x[ell],
                expected,
                rtol=1.0e-8,
                atol=0.0,
                err_msg=f"Taylor series mismatch for j_{ell}({x_val})",
            )

    @pytest.mark.parametrize(
        "lmax,x_val",
        [(10, 50.0), (20, 100.0), (50, 500.0), (100, 1000.0), (200, 20000.0)],
    )
    def test_array_upward_recursion(self, lmax: int, x_val: float) -> None:
        """Test array evaluation using upward recursion.

        For x > lmax+1, the implementation uses upward recursion.
        """
        sba = Ncm.SFSBesselArray.new_full(lmax, 1.0e-100)

        jl_x = sba.eval1(lmax, x_val)

        # Compare with scipy
        for ell in range(lmax + 1):
            expected = spherical_jn(ell, x_val)
            assert_allclose(
                jl_x[ell],
                expected,
                rtol=1.0e-14,
                atol=0.0,
                err_msg=f"Upward recursion mismatch for j_{ell}({x_val})",
            )

    @pytest.mark.parametrize(
        "lmax,x_val",
        [(10, 10.0), (20, 15.0), (50, 30.0), (100, 50.0), (200, 100.0)],
    )
    def test_array_steed_barnett(self, lmax: int, x_val: float) -> None:
        """Test array evaluation using Steed/Barnett algorithm.

        For x <= lmax+1, the implementation uses the Steed/Barnett algorithm.
        """
        sba = Ncm.SFSBesselArray.new_full(lmax, 1.0e-100)

        jl_x = sba.eval1(lmax, x_val)

        # Compare with scipy
        for ell in range(lmax + 1):
            expected = spherical_jn(ell, x_val)
            assert_allclose(
                jl_x[ell],
                expected,
                rtol=1.0e-10,
                atol=0.0,
                err_msg=f"Steed/Barnett mismatch for j_{ell}({x_val})",
            )


class TestSFSBesselRecursionRelations:
    """Test recursion relations and mathematical properties of spherical Bessel
    functions."""

    @pytest.mark.parametrize("l_val", [1, 2, 5, 10])
    @pytest.mark.parametrize("x_val", [1.0, 5.0, 10.0])
    def test_recursion_relation(self, l_val: int, x_val: float) -> None:
        """Test three-term recursion relation.

        j_{l+1}(x) = (2l+1)/x * j_l(x) - j_{l-1}(x)
        """
        jl_minus = Ncm.sf_sbessel(l_val - 1, x_val)
        jl = Ncm.sf_sbessel(l_val, x_val)
        jl_plus = Ncm.sf_sbessel(l_val + 1, x_val)

        expected_jl_plus = (2.0 * l_val + 1.0) / x_val * jl - jl_minus

        assert_allclose(
            jl_plus,
            expected_jl_plus,
            rtol=1.0e-13,
            atol=0.0,
            err_msg=f"Recursion relation failed for l={l_val}, x={x_val}",
        )


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

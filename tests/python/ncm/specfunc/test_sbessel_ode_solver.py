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

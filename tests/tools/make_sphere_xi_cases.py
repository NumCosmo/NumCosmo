#!/usr/bin/env python3
#
# make_sphere_xi_cases.py
#
# Fri Oct 16 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# make_sphere_xi_cases.py
# Copyright (C) 2026 Sandro Dias Pinto Vitenti <vitenti@uel.br>
#
# numcosmo is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# numcosmo is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program.  If not, see <http://www.gnu.org/licenses/>.

"""Print the certified tables of tests/c/ncm/sphere/test_ncm_sphere_xi.c.

Each row is one entry of the NcmSphereXi kernel

    K_{b ell} = (2 ell + 1) / (4 pi) [G_ell(x_b) - G_ell(x_{b+1})] / (x_b - x_{b+1}),
    G_ell(x) = [P_{ell+1}(x) - P_{ell-1}(x)] / (2 ell + 1),  K_{b 0} = 1 / (4 pi),

with x_b = cos(theta_b), evaluated by mpmath with 40 digits at the edges as doubles and
rounded to the nearest double. The edges cover a bin starting at zero, arcminute bins, a
bin across pi / 2 and a bin ending at pi.

For the spectrum C_ell = 1 / ((ell + 1) (ell + 2)), ell = CHEB_LMIN..CHEB_LMAX (the doubles
the C test computes the same way), zero below CHEB_LMIN, with a_ell = (2 ell + 1) / (4 pi) C_ell:
the Chebyshev coefficients

    c_k = (2 - delta_k0) / pi sum_ell Lambda((ell - k) / 2) Lambda((ell + k) / 2) a_ell,
    Lambda(z) = Gamma(z + 1/2) / Gamma(z + 1), ell - k even,

at CHEB_KS, and xi(theta) = sum_ell a_ell P_ell(cos theta) at CHEB_THETAS (zero and pi are
the first and last angles of the grid of ncm_sphere_xi_eval_chebyshev_grid), with 40 digits
and rounded to double. The output replaces the four blocks between the BEGIN and END markers
of the test:

    python3 tests/tools/make_sphere_xi_cases.py
"""

from __future__ import annotations

import math

import mpmath

ARCMIN = math.pi / 180.0 / 60.0

EDGES = [
    0.0,
    1.0 * ARCMIN,
    2.5 * ARCMIN,
    3.0 * ARCMIN,
    10.0 * ARCMIN,
    60.0 * ARCMIN,
    250.0 * ARCMIN,
    1.5,
    1.7,
    3.0,
    math.pi,
]

ELLS = [0, 1, 2, 3, 10, 100, 1000, 3000]

CHEB_LMIN = 2
CHEB_LMAX = 300
CHEB_KS = [0, 1, 2, 3, 10, 100, 299, 300]
CHEB_THETAS = [0.0, ARCMIN, 1.0e-2, 0.5, 1.5, 3.0, math.pi]


def kernel(ell: int, theta_lo: float, theta_hi: float) -> mpmath.mpf:
    """Return K_{b ell} for the bin [theta_lo, theta_hi]."""
    if ell == 0:
        return 1 / (4 * mpmath.pi)

    x_lo = mpmath.cos(mpmath.mpf(theta_lo))
    x_hi = mpmath.cos(mpmath.mpf(theta_hi))

    def G(x: mpmath.mpf) -> mpmath.mpf:
        return (mpmath.legendre(ell + 1, x) - mpmath.legendre(ell - 1, x)) / (
            2 * ell + 1
        )

    return (2 * ell + 1) / (4 * mpmath.pi) * (G(x_lo) - G(x_hi)) / (x_lo - x_hi)


def cheb_a() -> list[mpmath.mpf]:
    """Return a_ell for the Chebyshev test spectrum, ell = 0..CHEB_LMAX."""
    a = [mpmath.mpf(0)] * (CHEB_LMAX + 1)
    for ell in range(CHEB_LMIN, CHEB_LMAX + 1):
        c_ell = 1.0 / ((ell + 1.0) * (ell + 2.0))
        a[ell] = (2 * ell + 1) * mpmath.mpf(c_ell) / (4 * mpmath.pi)
    return a


def cheb_coeff(a: list[mpmath.mpf], k: int) -> mpmath.mpf:
    """Return c_k of the Chebyshev test spectrum."""

    def lam(z: mpmath.mpf) -> mpmath.mpf:
        return mpmath.gamma(z + mpmath.mpf(1) / 2) / mpmath.gamma(z + 1)

    s = mpmath.mpf(0)
    for ell in range(k, CHEB_LMAX + 1, 2):
        s += lam(mpmath.mpf(ell - k) / 2) * lam(mpmath.mpf(ell + k) / 2) * a[ell]
    return (1 if k == 0 else 2) * s / mpmath.pi


def cheb_xi(a: list[mpmath.mpf], theta: float) -> mpmath.mpf:
    """Return xi(theta) of the Chebyshev test spectrum."""
    x = mpmath.cos(mpmath.mpf(theta))
    return sum(a[ell] * mpmath.legendre(ell, x) for ell in range(CHEB_LMAX + 1))


def main() -> None:
    """Print the edges, rows, Chebyshev coefficients and xi blocks."""
    mpmath.mp.dps = 40

    print("  /* BEGIN make_sphere_xi_cases.py edges */")
    for theta in EDGES:
        print(f"  {theta!r},")
    print("  /* END make_sphere_xi_cases.py edges */")
    print()
    print("  /* BEGIN make_sphere_xi_cases.py rows */")
    for b in range(len(EDGES) - 1):
        for ell in ELLS:
            value = float(kernel(ell, EDGES[b], EDGES[b + 1]))
            print(f"  {{{b}, {ell}, {value!r}}},")
    print("  /* END make_sphere_xi_cases.py rows */")
    print()

    a = cheb_a()
    print("  /* BEGIN make_sphere_xi_cases.py cheb */")
    for k in CHEB_KS:
        print(f"  {{{k}, {float(cheb_coeff(a, k))!r}}},")
    print("  /* END make_sphere_xi_cases.py cheb */")
    print()
    print("  /* BEGIN make_sphere_xi_cases.py xi */")
    for theta in CHEB_THETAS:
        print(f"  {{{theta!r}, {float(cheb_xi(a, theta))!r}}},")
    print("  /* END make_sphere_xi_cases.py xi */")


if __name__ == "__main__":
    main()

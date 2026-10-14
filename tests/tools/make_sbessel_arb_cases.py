#!/usr/bin/env python3
#
# make_sbessel_arb_cases.py
#
# Wed Oct 14 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# make_sbessel_arb_cases.py
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

"""Print a certified table of the C sbessel tests (tests/c/ncm/specfunc).

Each row is one integral

    I = int_a^b f(x) j_ell^(d)(k x) dx

certified by ``nc_xcor_kernel_analytic_arb --unnormalized`` to a relative radius of
1e-25 and rounded to the nearest double. The shapes and their two parameters p1, p2,
written the same way in the tool and in the C test:

- ``tophat``: f = 1;
- ``power``: f = x^n, n = p1;
- ``cos``: f = cos (omega x), omega = p1;
- ``exp``: f = exp (-x / s), s = p1 (the tool's ``power_exp`` with alpha 0, beta 1);
- ``exp2``: f = exp (-(x / s)^2), s = p1 (``power_exp`` with alpha 0, beta 2);
- ``sin``: f = sin (omega x), omega = p1;
- ``lorentz``: f = 1 / (1 + ((x - c) / w)^2), c = p1, w = p2;
- ``rational``: f = x^2 / (1 + ((x - c) / w)^2)^3, c = p1, w = p2;
- ``gauss``: f = exp (-((x - c) / w)^2 / 2), c = p1, w = p2.

A shape name ending in ``/x`` divides f by x (the tool's ``--over-x``).

Two tables, chosen by the first argument: ``levin`` for
test_ncm_sbessel_integrator_levin.c and ``ode`` for test_ncm_sbessel_ode_solver.c. The
output replaces the block between the BEGIN and END markers of that test:

    python3 tests/tools/make_sbessel_arb_cases.py levin

The tool is the binary named by XCOR_ARB_TOOL, else the meson-built one in
NUMCOSMO_BUILD_DIR (set by the build's numcosmo_export.sh).
"""

from __future__ import annotations

import os
import pathlib
import subprocess
import sys

import mpmath

TOOL = "nc_xcor_kernel_analytic_arb"

# name, shape, a, b, p1, p2, deriv, [(ell, [k, ...]), ...]
LEVIN_CASES = [
    (
        "constant",
        "tophat",
        1.0,
        20.0,
        0.0,
        0.0,
        0,
        [(ell, [1.0]) for ell in (0, 1, 5, 10, 20)],
    ),
    (
        "linear",
        "power",
        1.0,
        20.0,
        1.0,
        0.0,
        0,
        [(ell, [1.0]) for ell in (0, 1, 5, 10, 20)],
    ),
    (
        "quadratic",
        "power",
        1.0,
        20.0,
        2.0,
        0.0,
        0,
        [(ell, [1.0]) for ell in (0, 1, 5, 10)],
    ),
    (
        "rational",
        "rational",
        1.0,
        20.0,
        10.5,
        3.0,
        0,
        # ell 200 at k 20: the turning point x ~ 200.5 / k sits on the peak
        [(ell, [1.0]) for ell in (0, 1, 5, 10, 20)] + [(200, [20.0])],
    ),
    (
        "moving_edge_0",
        "exp",
        5.0,
        50.0,
        1.0 / 0.03,
        0.0,
        0,
        [(ell, [1.0]) for ell in range(8)],
    ),
    (
        "moving_edge_1",
        "exp",
        6.7,
        43.0,
        1.0 / 0.03,
        0.0,
        0,
        [(ell, [1.0]) for ell in range(8)],
    ),
    (
        "moving_edge_2",
        "exp",
        9.5,
        32.0,
        1.0 / 0.03,
        0.0,
        0,
        [(ell, [1.0]) for ell in range(8)],
    ),
    (
        "moving_edge_3",
        "exp",
        10.5,
        96.0,
        1.0 / 0.03,
        0.0,
        0,
        [(ell, [1.0]) for ell in range(8)],
    ),
    (
        "accepted_edge",
        "tophat",
        1.0,
        2.0,
        0.0,
        0.0,
        0,
        [(ell, [1.0]) for ell in range(8)],
    ),
    ("rejected_edge", "cos", 1.0, 2.0, 500.0, 0.0, 0, [(0, [1.0])]),
    (
        "set_reltol",
        "exp",
        100.0,
        1.0e4,
        1000.0,
        0.0,
        0,
        [(ell, [1.0]) for ell in range(8)],
    ),
    (
        "repeated_set_reltol",
        "exp",
        12.0,
        47.0,
        100.0,
        0.0,
        0,
        [(ell, [1.0]) for ell in range(6)],
    ),
    (
        # the rational truth table's shape at k = 10^0.5: one panel holds the whole
        # integral and the solution needs some 600 coefficients
        "reltol_bound",
        "rational",
        0.01,
        6.5,
        1.5,
        0.05,
        0,
        [(ell, [0.1 * 10.0**1.5]) for ell in range(8)],
    ),
    ("ell_zero", "tophat", 1.0, 10.0, 0.0, 0.0, 0, [(0, [1.0])]),
    ("k_scaling", "tophat", 1.0, 10.0, 0.0, 0.0, 0, [(5, [0.1, 1.0, 10.0, 100.0])]),
    ("narrow", "power", 5.0, 5.001, 1.0, 0.0, 0, [(5, [1.0])]),
    (
        "high_ell",
        "tophat",
        1.0,
        20.0,
        0.0,
        0.0,
        0,
        [(ell, [float(ell)]) for ell in (100, 200, 300)],
    ),
    ("oscillatory", "exp", 1.0, 100.0, 100.0, 0.0, 0, [(10, [500.0])]),
    ("evanescent", "gauss", 1.0, 8.0, 4.0, 1.2, 2, [(30, [0.5])]),
]


def _ode_bump(name: str, a: float, b: float, width: float, ell: int = 20) -> tuple:
    """A Gaussian, width ``width`` times the half-width, centred on [a, b]."""
    m, h = 0.5 * (a + b), 0.5 * (b - a)

    return (name, "gauss", a, b, m, width * h, 0, [(ell, [1.0])])


def _conjugate_span() -> float:
    """The root s in [7.9 pi, 8.1 pi] of j_2(200) y_2(200 + s) - j_2(200 + s) y_2(200).

    Computed with mpmath at 30 digits and rounded to double: the two-point Dirichlet
    problem at ell 2 on [200, 200 + s] has no unique solution.
    """
    mpmath.mp.dps = 30

    def jl(x):
        return mpmath.sqrt(mpmath.pi / (2 * x)) * mpmath.besselj(2.5, x)

    def yl(x):
        return mpmath.sqrt(mpmath.pi / (2 * x)) * mpmath.bessely(2.5, x)

    a = mpmath.mpf(200)
    root = mpmath.findroot(
        lambda s: jl(a) * yl(a + s) - jl(a + s) * yl(a),
        (7.9 * mpmath.pi, 8.1 * mpmath.pi),
        solver="anderson",
    )

    return float(root)


# int_1^20 j_l(x) / x dx and int_1^20 j_l(x) dx, ell 0..20: the Green identity of the
# ODE solver for the forcings 1 and x; then int_1^20 g(x) j_l(x) / x dx for the forcings
# g = x^2, x^3, exp (x), sin (x), cos (x), exp (-x^2) and 1 / (1 + x^2)
ODE_GREEN_ELLS = (0, 2, 3, 4, 5, 6, 7, 8, 10)
ODE_CASES = [
    (
        "green_inv_x",
        "power",
        1.0,
        20.0,
        -1.0,
        0.0,
        0,
        [(ell, [1.0]) for ell in range(21)],
    ),
    (
        "green_one",
        "tophat",
        1.0,
        20.0,
        0.0,
        0.0,
        0,
        [(ell, [1.0]) for ell in range(21)],
    ),
    (
        "rhs_quadratic",
        "power/x",
        1.0,
        20.0,
        2.0,
        0.0,
        0,
        [(ell, [1.0]) for ell in (0, 5, 10, 15)],
    ),
    (
        "rhs_cubic",
        "power/x",
        1.0,
        20.0,
        3.0,
        0.0,
        0,
        [(ell, [1.0]) for ell in (0, 5, 10)],
    ),
    (
        "rhs_exp",
        "exp/x",
        1.0,
        20.0,
        -1.0,
        0.0,
        0,
        [(ell, [1.0]) for ell in ODE_GREEN_ELLS],
    ),
    (
        "rhs_sin",
        "sin/x",
        1.0,
        20.0,
        1.0,
        0.0,
        0,
        [(ell, [1.0]) for ell in ODE_GREEN_ELLS],
    ),
    ("rhs_cos", "cos/x", 1.0, 20.0, 1.0, 0.0, 0, [(ell, [1.0]) for ell in (0, 5, 10)]),
    (
        "rhs_gauss",
        "exp2/x",
        1.0,
        20.0,
        1.0,
        0.0,
        0,
        [(ell, [1.0]) for ell in (0, 5, 10)],
    ),
    (
        "rhs_lorentz",
        "lorentz/x",
        1.0,
        20.0,
        0.0,
        1.0,
        0,
        [(ell, [1.0]) for ell in (0, 5, 10)],
    ),
    # int_a^b F j_20 dx, F a Gaussian centred on the panel: width 0.35 h on the narrow
    # tau panel [100, 10^2.5], 0.25 h on the panels [1e3, 1e3 + span] of the derivative
    # error bound (h the half-width)
    _ode_bump("tau_narrow", 100.0, 10.0**2.5, 0.35),
    *[
        _ode_bump(f"deriv_bound_{int(span)}", 1.0e3, 1.0e3 + span, 0.25)
        for span in (30.0, 60.0, 100.0, 150.0, 250.0, 3000.0)
    ],
    # int_a^b F j_2 dx on [200, 200 + s] with s the conjugate span (the Dirichlet
    # boundary matrix is singular there) and on [200, 200 + s + 0.01], F a Gaussian of
    # width 0.35 h
    _ode_bump("conjugate", 200.0, 200.0 + _conjugate_span(), 0.35, ell=2),
    _ode_bump("conjugate_away", 200.0, 200.0 + _conjugate_span() + 1.0e-2, 0.35, ell=2),
]

TABLES = {"levin": LEVIN_CASES, "ode": ODE_CASES}

C_SHAPE = {
    "tophat": "TEST_ARB_TOPHAT",
    "power": "TEST_ARB_POWER",
    "cos": "TEST_ARB_COS",
    "exp": "TEST_ARB_EXP",
    "rational": "TEST_ARB_RATIONAL",
    "gauss": "TEST_ARB_GAUSS",
    "exp2": "TEST_ARB_EXP2",
    "sin": "TEST_ARB_SIN",
    "lorentz": "TEST_ARB_LORENTZ",
}


def find_tool() -> pathlib.Path:
    """The Arb tool named by XCOR_ARB_TOOL, else the one in NUMCOSMO_BUILD_DIR."""
    override = os.environ.get("XCOR_ARB_TOOL")

    if override:
        return pathlib.Path(override)

    if "NUMCOSMO_BUILD_DIR" not in os.environ:
        sys.exit("NUMCOSMO_BUILD_DIR is not set; source the build's numcosmo_export.sh")

    tool = pathlib.Path(os.environ["NUMCOSMO_BUILD_DIR"]) / "tests" / "tools" / TOOL

    if not (tool.is_file() and os.access(tool, os.X_OK)):
        sys.exit(f"{tool} not found; build the {TOOL} target")

    return tool


def tool_args(shape: str, a: float, b: float, p1: float, p2: float) -> list[str]:
    """The tool's options for one shape, every double printed in full."""
    support = [f"--chi-lower={a!r}", f"--chi-upper={b!r}"]

    if shape == "tophat":
        return ["--shape=tophat", *support]
    if shape == "power":
        return ["--shape=power", f"--power-n={int(p1)}", *support]
    if shape.endswith("/x"):
        return [*tool_args(shape[:-2], a, b, p1, p2), "--over-x"]
    if shape == "cos":
        return ["--shape=cos", f"--omega={p1!r}", *support]
    if shape == "sin":
        return ["--shape=sin", f"--omega={p1!r}", *support]
    if shape == "exp2":
        return [
            "--shape=power_exp",
            "--alpha=0",
            "--beta=2",
            f"--chi-scale={p1!r}",
            *support,
        ]
    if shape == "lorentz":
        return [
            "--shape=lorentz",
            f"--chi-mean={p1!r}",
            f"--chi-scale={p2!r}",
            *support,
        ]
    if shape == "exp":
        return [
            "--shape=power_exp",
            "--alpha=0",
            "--beta=1",
            f"--chi-scale={p1!r}",
            *support,
        ]
    if shape == "rational":
        return [
            "--shape=rational",
            f"--chi-mean={p1!r}",
            f"--chi-scale={p2!r}",
            *support,
        ]
    if shape == "gauss":
        return [
            "--shape=gauss",
            f"--chi-mean={p1!r}",
            f"--chi-sigma={p2!r}",
            "--explicit-support",
            *support,
        ]

    raise ValueError(shape)


def certify(tool, shape, a, b, p1, p2, deriv, ell, ks) -> list[tuple[float, str]]:
    """(k, certified value as printed by the tool) for one multipole."""
    out = subprocess.run(
        [
            str(tool),
            *tool_args(shape, a, b, p1, p2),
            f"--ell={ell}",
            f"--bessel-deriv={deriv}",
            "--unnormalized",
        ],
        input=" ".join(repr(k) for k in ks),
        capture_output=True,
        text=True,
        check=True,
    ).stdout
    rows = [line.split("\t") for line in out.splitlines() if not line.startswith("#")]

    if len(rows) != len(ks):
        sys.exit(f"{shape} ell {ell}: expected {len(ks)} rows, got {len(rows)}")

    return [(float(r[2]), r[3]) for r in rows]


def main() -> None:
    """Print the C table named by the first argument."""
    if len(sys.argv) != 2 or sys.argv[1] not in TABLES:
        sys.exit(f"usage: {sys.argv[0]} {{{'|'.join(TABLES)}}}")

    table = sys.argv[1]
    tool = find_tool()
    print(f"  /* BEGIN make_sbessel_arb_cases.py {table} */")

    for name, shape, a, b, p1, p2, deriv, cells in TABLES[table]:
        for ell, ks in cells:
            base = shape.removesuffix("/x")
            over_x = "TRUE" if shape.endswith("/x") else "FALSE"
            for k, value in certify(tool, shape, a, b, p1, p2, deriv, ell, ks):
                print(
                    f'  {{"{name}", {C_SHAPE[base]}, {a!r}, {b!r}, {p1!r}, {p2!r}, '
                    f"{deriv}, {over_x}, {ell}, {k!r}, {float(value)!r}}},"
                )

    print(f"  /* END make_sbessel_arb_cases.py {table} */")


if __name__ == "__main__":
    main()

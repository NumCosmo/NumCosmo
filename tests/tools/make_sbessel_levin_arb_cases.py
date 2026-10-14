#!/usr/bin/env python3
#
# make_sbessel_levin_arb_cases.py
#
# Wed Oct 14 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# make_sbessel_levin_arb_cases.py
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

"""Print the certified table of the C Levin test (test_ncm_sbessel_integrator_levin.c).

Each row is one integral

    I = int_a^b f(x) j_ell^(d)(k x) dx

certified by ``nc_xcor_kernel_analytic_arb --unnormalized`` to a relative radius of
1e-25 and rounded to the nearest double. The shapes and their two parameters p1, p2,
written the same way in the tool and in the C test:

- ``tophat``: f = 1;
- ``power``: f = x^n, n = p1;
- ``cos``: f = cos (omega x), omega = p1;
- ``exp``: f = exp (-x / s), s = p1 (the tool's ``power_exp`` with alpha 0, beta 1);
- ``rational``: f = x^2 / (1 + ((x - c) / w)^2)^3, c = p1, w = p2;
- ``gauss``: f = exp (-((x - c) / w)^2 / 2), c = p1, w = p2.

The output replaces the block between the BEGIN and END markers of the C test:

    python3 tests/tools/make_sbessel_levin_arb_cases.py

The tool is the binary named by XCOR_ARB_TOOL, else the meson-built one in
NUMCOSMO_BUILD_DIR (set by the build's numcosmo_export.sh).
"""

from __future__ import annotations

import os
import pathlib
import subprocess
import sys

TOOL = "nc_xcor_kernel_analytic_arb"

# name, shape, a, b, p1, p2, deriv, [(ell, [k, ...]), ...]
CASES = [
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

C_SHAPE = {
    "tophat": "TEST_ARB_TOPHAT",
    "power": "TEST_ARB_POWER",
    "cos": "TEST_ARB_COS",
    "exp": "TEST_ARB_EXP",
    "rational": "TEST_ARB_RATIONAL",
    "gauss": "TEST_ARB_GAUSS",
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
    if shape == "cos":
        return ["--shape=cos", f"--omega={p1!r}", *support]
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
    """Print the C table."""
    tool = find_tool()
    print("  /* BEGIN make_sbessel_levin_arb_cases.py */")

    for name, shape, a, b, p1, p2, deriv, cells in CASES:
        for ell, ks in cells:
            for k, value in certify(tool, shape, a, b, p1, p2, deriv, ell, ks):
                print(
                    f'  {{"{name}", {C_SHAPE[shape]}, {a!r}, {b!r}, {p1!r}, {p2!r}, '
                    f"{deriv}, {ell}, {k!r}, {float(value)!r}}},"
                )

    print("  /* END make_sbessel_levin_arb_cases.py */")


if __name__ == "__main__":
    main()

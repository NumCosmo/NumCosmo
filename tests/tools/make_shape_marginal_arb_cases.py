#!/usr/bin/env python3
#
# make_shape_marginal_arb_cases.py
#
# Tue Nov 10 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# make_shape_marginal_arb_cases.py
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

"""Print the certified table of the shape-marginal tests (tests/c/nc/lss/galaxy).

Each row is one value

    ln P (eps_obs | g) = ln int_{|chi_I| < 1} P_2D (chi_I) N_2 (eps_obs - f_g (chi_I)) d^2 chi_I

certified by ``nc_galaxy_shape_marginal_arb`` (polar coordinates, its default target
radius) and rounded to the nearest double. The columns are the tool's input columns: the
shear convention (``trace`` or ``trace_det``), the population (``gauss`` with sigma = p1,
or ``beta`` with alpha = p1 and beta = p2), std_noise, g1, g2, e1, e2. The last two
columns are for the NcGalaxyShapeFactorQuad tests: whether its unit tier checks the row
(UNIT_ROWS) and whether its acceptance tier does (the first QUAD_ROWS rows). The
NcGalaxyShapeFactorSpectral test checks every row.

The output replaces the block between the BEGIN and END markers of
test_nc_galaxy_shape_marginal_arb_cases.h:

    python3 tests/tools/make_shape_marginal_arb_cases.py

The tool is the binary named by SHAPE_ARB_TOOL, else the meson-built one in
NUMCOSMO_BUILD_DIR (set by the build's numcosmo_export.sh). A point the tool does not
certify stops the script.
"""

from __future__ import annotations

import os
import pathlib
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

TOOL = "nc_galaxy_shape_marginal_arb"

C_CONV = {
    "trace": "NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE",
    "trace_det": "NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET",
}
C_POP = {"gauss": "TEST_SHAPE_ARB_POP_GAUSS", "beta": "TEST_SHAPE_ARB_POP_BETA"}

# conv, pop, p1, p2, std_noise, g1, g2, e1, e2
CASES = [
    ("trace_det", "gauss", 0.3, 0.0, 0.1, 0.3, 0.0, 0.0, 0.5),
    ("trace", "gauss", 0.3, 0.0, 0.01, 0.6, 0.0, 0.5, 0.0),
    ("trace", "gauss", 0.05, 0.0, 0.005, 0.1, 0.0, 0.05, 0.0),
    ("trace_det", "gauss", 0.05, 0.0, 0.03, 0.6, 0.0, -0.9, 0.0),
    ("trace_det", "gauss", 0.2, 0.0, 0.05, 0.2, -0.4, 0.3, 0.25),
    ("trace", "gauss", 0.3, 0.0, 1.0, 3.0, 0.0, 2.0, 0.0),
    ("trace_det", "gauss", 0.3, 0.0, 0.1, 1.5, 0.0, 1.2, 0.0),
    ("trace", "beta", 1.0, 1.0, 0.03, 0.6, 0.0, 0.0, 0.5),
    ("trace_det", "beta", 4.0, 1.0, 0.1, 0.6, 0.0, 0.9, 0.0),
    ("trace_det", "beta", 1.0, 1.0, 0.01, 0.3, 0.0, 0.0, 1.2),
    ("trace", "beta", 4.0, 1.0, 0.3, 1.5, 0.0, -0.5, 0.0),
    ("trace", "beta", 1.57, 1.65, 0.1, 0.3, 0.0, 0.5, 0.0),
    ("trace_det", "beta", 1.92, 1.98, 0.03, 0.6, 0.0, 0.0, 0.5),
    ("trace_det", "beta", 1.57, 1.65, 0.05, -0.1, 0.25, 0.2, -0.3),
    # Points where the default path was wrong (6.6), off by 7.5e-6, or did not return.
    ("trace_det", "gauss", 0.05, 0.0, 0.1, 0.005111126056639749, 0.0, 0.0, 0.9),
    ("trace_det", "beta", 1.57, 1.65, 0.1, 0.01, 0.0, 0.0, 0.02),
    ("trace", "beta", 4.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0),
    # Complex g at an angle to eps_obs. The first three are points where integration rules
    # with a window around the peak were wrong by 3.8e-4, 2.7e-7 and 4.8e-8.
    ("trace", "gauss", 0.05, 0.0, 0.03, 0.3, -0.219359538352843, 0.99, 0.0),
    (
        "trace",
        "gauss",
        0.05,
        0.0,
        0.005,
        0.3227416423252818,
        0.46019288007369547,
        0.99,
        0.0,
    ),
    (
        "trace",
        "gauss",
        0.4,
        0.0,
        0.01,
        -0.86591447813628053,
        0.24534081713806422,
        -0.93403982009271569,
        -0.32813048392547611,
    ),
    (
        "trace_det",
        "gauss",
        0.2,
        0.0,
        0.03,
        -0.36850938805797601,
        0.47349850148985267,
        -0.89999983214337287,
        0.00054967435876962577,
    ),
    (
        "trace",
        "beta",
        1.5,
        3.0,
        0.1,
        0.0020353798786394857,
        -0.0097906704954068141,
        -0.58556538052086571,
        -0.79825634049185845,
    ),
    (
        "trace_det",
        "beta",
        4.0,
        1.0,
        0.03,
        -0.53863697645386155,
        0.26432973271397636,
        0.80049145894708595,
        0.41135559331649529,
    ),
]

# The NcGalaxyShapeFactorQuad acceptance tier checks the first QUAD_ROWS rows.
QUAD_ROWS = 17


# Indices of CASES checked by the unit tier: six cases on which the default and log-polar
# paths take about 1 s in total (|g| < 1 and > 1, |eps_obs| < 1 and >= 1, Gauss and Beta,
# and two of the points where the default path had failed).
UNIT_ROWS = {1, 2, 6, 11, 14, 16}


def find_tool() -> pathlib.Path:
    """The Arb tool named by SHAPE_ARB_TOOL, else the one in NUMCOSMO_BUILD_DIR."""
    override = os.environ.get("SHAPE_ARB_TOOL")

    if override:
        return pathlib.Path(override)

    if "NUMCOSMO_BUILD_DIR" not in os.environ:
        sys.exit("NUMCOSMO_BUILD_DIR is not set; source the build's numcosmo_export.sh")

    tool = pathlib.Path(os.environ["NUMCOSMO_BUILD_DIR"]) / "tests" / "tools" / TOOL

    if not (tool.is_file() and os.access(tool, os.X_OK)):
        sys.exit(f"{tool} not found; build the {TOOL} target")

    return tool


def certify(tool: pathlib.Path, case: tuple) -> str:
    """ln P of one case as printed by the tool."""
    line = " ".join(c if isinstance(c, str) else repr(c) for c in case)
    out = subprocess.run(
        [str(tool), "--coords=polar"],
        input=line + "\n",
        capture_output=True,
        text=True,
        check=True,
    ).stdout
    rows = [r.split("\t") for r in out.splitlines() if not r.startswith("#")]

    if len(rows) != 1 or int(rows[0][11]) == 0:
        sys.exit(f"not certified: {line}")

    return rows[0][9]


def main() -> None:
    """Print the C table; SHAPE_ARB_JOBS sets the number of concurrent tool runs."""
    tool = find_tool()
    jobs = int(os.environ.get("SHAPE_ARB_JOBS", "1"))

    with ThreadPoolExecutor(max_workers=jobs) as pool:
        values = list(pool.map(lambda c: certify(tool, c), CASES))

    print("  /* BEGIN make_shape_marginal_arb_cases.py */")

    for i, ((conv, pop, p1, p2, sn, g1, g2, e1, e2), value) in enumerate(
        zip(CASES, values)
    ):
        print(
            f"  {{{C_CONV[conv]}, {C_POP[pop]}, {p1!r}, {p2!r}, {sn!r}, "
            f"{g1!r}, {g2!r}, {e1!r}, {e2!r}, {float(value)!r}, "
            f"{'TRUE' if i in UNIT_ROWS else 'FALSE'}, "
            f"{'TRUE' if i < QUAD_ROWS else 'FALSE'}}},"
        )

    print("  /* END make_shape_marginal_arb_cases.py */")


if __name__ == "__main__":
    main()

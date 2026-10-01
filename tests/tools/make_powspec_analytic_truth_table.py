#!/usr/bin/env python
#
# make_powspec_analytic_truth_table.py
#
# Mon Sep 28 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# make_powspec_analytic_truth_table.py
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

"""Generate the Arb truth table of the NcmPowspec integrals.

Runs ``ncm_powspec_analytic_arb --integrals`` (NcmPowspecAnalytic at its defaults,
BBKS shape and LCDM growth, integrated over ``k`` in [1e-6, 1e2] Mpc^-1) and writes
``data/truth_tables/powspec/ncm_powspec_analytic_integrals.bin``, an NcmObjDictStr:

``k_range``
    NcmVector ``k_lo, k_hi``.
``var``
    NcmMatrix, rows ``z, R, value, radius``: ncm_powspec_var_tophat_R().
``xi``
    NcmMatrix, rows ``z, r, value, radius``: ncm_powspec_corr3d().
``sproj``
    NcmMatrix, rows ``ell, z1, z2, xi1, xi2, value, radius``: ncm_powspec_sproj().

``radius`` is the certified Arb radius of ``value`` before rounding to double; each
value is the nearest double to the Arb midpoint. Lengths are in Mpc. The consumer is
``tests/c/ncm/powspec/test_ncm_powspec.c``.
"""

import argparse
import subprocess
from pathlib import Path

from numcosmo_py import Ncm

ROOT = Path(__file__).resolve().parents[2]


def _matrix(rows: list[list[float]]) -> Ncm.Matrix:
    return Ncm.Matrix.new_array([x for row in rows for x in row], len(rows[0]))


def main() -> None:
    """Write the table."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--tool",
        type=Path,
        default=ROOT / "Optimized" / "tests" / "tools" / "ncm_powspec_analytic_arb",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT
        / "data"
        / "truth_tables"
        / "powspec"
        / "ncm_powspec_analytic_integrals.bin",
    )
    args = parser.parse_args()

    k_lo, k_hi = 1.0e-6, 1.0e2
    out = subprocess.run(
        [
            str(args.tool),
            "--integrals",
            f"--k-lo={k_lo!r}",
            f"--k-hi={k_hi!r}",
            "--target-rel=1e-20",
        ],
        check=True,
        capture_output=True,
        text=True,
    ).stdout

    rows: dict[str, list[list[float]]] = {"var": [], "xi": [], "sproj": []}
    for line in out.splitlines():
        if line.startswith("#"):
            continue
        kind, ell, z1, z2, s1, s2, value, radius, _ = line.split("\t")
        if kind == "sproj":
            row = [ell, z1, z2, s1, s2, value, radius]
        else:
            row = [z1, s1, value, radius]
        rows[kind].append([float(x) for x in row])

    Ncm.cfg_init()
    table = Ncm.ObjDictStr.new()
    table.add("k_range", Ncm.Vector.new_array([k_lo, k_hi]))
    for kind, kind_rows in rows.items():
        table.add(kind, _matrix(kind_rows))

    args.output.parent.mkdir(parents=True, exist_ok=True)
    ser = Ncm.Serialize.new(Ncm.SerializeOpt.NONE)
    ser.dict_str_to_binfile(table, str(args.output))


if __name__ == "__main__":
    main()

#!/usr/bin/env python
#
# make_sbessel_truth_table.py
#
# Thu Sep 24 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# make_sbessel_truth_table.py
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

r"""Certify the sbessel truth tables in Arb.

The tables in ``data/truth_tables/sbessel/`` hold

.. math:: I_l(k) = \int_a^b f(x)\, j_l(k x)\, \mathrm{d}x

for the unnormalized integrands of ``ncm_sbessel_integrator_integrate_gaussian`` and
``_rational``, on every l of the table and 101 log-spaced k. This driver computes the
same cells with ``nc_xcor_kernel_analytic_arb`` (``--unnormalized``, and for the
Gaussian ``--explicit-support`` so the limits are exactly the table's doubles), one
process per l, and writes the generator's full-precision output: decimal midpoints and
certified radii.

Two kinds of cell are not certified, and the output says which:

- **Proven negligible.** Since :math:`|j_l(x)| \le x^l / (2l+1)!!` for real
  :math:`x \ge 0` (DLMF 10.14.4 at :math:`\nu = l + 1/2`),
  :math:`|I_l(k)| \le F\,(k b)^l / (2l+1)!!` with :math:`F` a bound on
  :math:`\int_a^b |f|`, which also bounds every cell of the table. Where that is
  below :data:`NEGLIGIBLE` times :math:`F` (and so below anything a double resolves
  against the table's values) or below the smallest normal double, the cell is
  skipped without running Arb and the bound is recorded. Certifying such a cell to a
  relative target costs precision without limit: one at 2.5e-308, l = 100, did not
  finish at 32768 bits.
- **Beyond the precision ceiling.** A cell Arb cannot certify to the target within
  ``--prec-max`` bits is retried alone at ``--prec-max-retry``; one that still fails is
  recorded as uncertified.

Usage::

    python3 tests/tools/make_sbessel_truth_table.py --shape gauss --tool TOOL \
        --output gauss_jl_500.arb.json.gz [--jobs 12] [--resume]
"""

from __future__ import annotations

import argparse
import concurrent.futures
import gzip
import json
import math
import pathlib
import subprocess
import sys
import time

HERE = pathlib.Path(__file__).resolve().parent
TABLES = HERE.parent.parent / "data" / "truth_tables" / "sbessel"
DBL_MIN = sys.float_info.min
NEGLIGIBLE = 1.0e-250

# The integrands, as the library and the committed tables define them.
SHAPES = {
    "gauss": {
        "table": "gauss_jl_500.json.gz",
        "args": [
            "--shape=gauss",
            "--chi-mean={center!r}",
            "--chi-sigma={std!r}",
            "--n-sigma=8",
            "--explicit-support",
            "--chi-lower={lower!r}",
            "--chi-upper={upper!r}",
        ],
        # int_a^b exp(-z^2/2) dx <= std sqrt(2 pi)
        "l1_bound": lambda t: t["std"] * math.sqrt(2.0 * math.pi),
    },
    "rational": {
        "table": "rational_jl_500.json.gz",
        "args": [
            "--shape=rational",
            "--chi-mean={center!r}",
            "--chi-scale={std!r}",
            "--chi-lower={lower!r}",
            "--chi-upper={upper!r}",
        ],
        # int_a^b x^2 (1 + z^2)^-3 dx <= b^2 std int (1 + z^2)^-3 dz = b^2 std 3 pi / 8
        "l1_bound": lambda t: t["upper-bound"] ** 2 * t["std"] * 3.0 * math.pi / 8.0,
    },
}

CONVENTION = (
    "table[l][i] = int_lower-bound^upper-bound f(x) j_l(kvals[i] x) dx, with f the "
    "unnormalized integrand named by shape (the one ncm_sbessel_integrator_integrate_"
    "gaussian or _rational integrates). Note the measure: dx, not d(kx). Values are "
    "full-precision decimal strings certified in Arb; radius is the certified absolute "
    "half-width, so |true - value| <= radius. A null value is a cell not certified: "
    "status says why, and bound is a proven bound on |true| for a negligible cell."
)


def log_double_factorial_odd(l: int) -> float:
    """Return log((2l+1)!!)."""
    return math.lgamma(2 * l + 2) - l * math.log(2.0) - math.lgamma(l + 1)


def underflow_bound(l: int, k: float, upper: float, l1: float) -> float:
    """Return a proven bound on |I_l(k)|, from |j_l(x)| <= x^l / (2l+1)!!."""
    # |j_l| <= 1 as well, so the integral never exceeds l1 whatever the power bound says.
    if l == 0:
        return l1

    log_bound = l * math.log(k * upper) - log_double_factorial_odd(l)

    if log_bound >= 0.0:
        return l1

    return l1 * math.exp(log_bound) if log_bound > -745.0 else 0.0


def run_tool(tool: str, args: list[str], ell: int, ks: list[float], prec_max: int):
    """Run the Arb tool for one l; return (values, radii, precs) or None on failure."""
    proc = subprocess.run(
        [tool, *args, f"--ell={ell}", f"--prec-max={prec_max}", "--unnormalized"],
        input="\n".join(repr(k) for k in ks) + "\n",
        capture_output=True,
        text=True,
        check=False,
    )
    rows = [
        line.split("\t")
        for line in proc.stdout.splitlines()
        if not line.startswith("#")
    ]

    if proc.returncode != 0 or len(rows) != len(ks):
        return None, proc.stderr.strip()

    return [(r[3], float(r[4]), int(r[5])) for r in rows], ""


def one_ell(tool, shape, table, ell, prec_max, prec_retry):
    """Certify every cell of one l."""
    spec = SHAPES[shape]
    fmt = {
        "center": table["center"],
        "std": table["std"],
        "lower": table["lower-bound"],
        "upper": table["upper-bound"],
    }
    args = [a.format(**fmt) for a in spec["args"]]
    l1 = spec["l1_bound"](table)
    ks = table["kvals"]
    cells: list[dict] = [{} for _ in ks]
    todo = []

    for i, k in enumerate(ks):
        bound = underflow_bound(ell, k, table["upper-bound"], l1)

        if bound < max(DBL_MIN, NEGLIGIBLE * l1):
            cells[i] = {
                "value": None,
                "radius": None,
                "status": "negligible",
                "bound": bound,
            }
        else:
            todo.append(i)

    t0 = time.monotonic()
    got, err = (
        run_tool(tool, args, ell, [ks[i] for i in todo], prec_max) if todo else ([], "")
    )

    if got is None:
        got = []

        for i in todo:
            one, err1 = run_tool(tool, args, ell, [ks[i]], prec_max)

            if one is None:
                one, err1 = run_tool(tool, args, ell, [ks[i]], prec_retry)

            got.append(one[0] if one is not None else (None, None, None, err1))

    for i, row in zip(todo, got):
        if row[0] is None:
            cells[i] = {
                "value": None,
                "radius": None,
                "status": "uncertified",
                "error": row[3],
            }
        else:
            cells[i] = {
                "value": row[0],
                "radius": row[1],
                "prec": row[2],
                "status": "certified",
            }

    print(
        f"  {shape} l={ell:<4d} {len(todo):3d} certified-or-tried, "
        f"{len(ks) - len(todo):3d} proven negligible, {time.monotonic() - t0:7.1f} s",
        file=sys.stderr,
        flush=True,
    )

    return cells


def main() -> int:
    """Certify one sbessel table and write the full-precision output."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--shape", choices=list(SHAPES), required=True)
    parser.add_argument("--tool", required=True, help="nc_xcor_kernel_analytic_arb")
    parser.add_argument("--output", type=pathlib.Path, required=True)
    parser.add_argument("--jobs", type=int, default=12)
    parser.add_argument("--prec-max", type=int, default=8192)
    parser.add_argument("--prec-max-retry", type=int, default=32768)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument(
        "--ells", default=None, help="comma-separated subset of l, for timing probes"
    )
    args = parser.parse_args()

    with gzip.open(TABLES / SHAPES[args.shape]["table"], "rt") as handle:
        table = json.load(handle)

    partial = args.output.with_suffix(".partial.json")
    done: dict[str, list[dict]] = {}

    if args.resume and partial.exists():
        done = json.loads(partial.read_text(encoding="utf-8"))

    ells = [int(l) for l in table["lvals"]]

    if args.ells is not None:
        ells = [int(l) for l in args.ells.split(",")]

    # most expensive first: the cost grows like l^2, so a high l started last is a tail
    tasks = sorted((l for l in ells if str(l) not in done), reverse=True)
    print(f"{len(tasks)} of {len(ells)} multipoles to certify", file=sys.stderr)

    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futures = {
            pool.submit(
                one_ell,
                args.tool,
                args.shape,
                table,
                l,
                args.prec_max,
                args.prec_max_retry,
            ): l
            for l in tasks
        }

        for future in concurrent.futures.as_completed(futures):
            done[str(futures[future])] = future.result()
            partial.write_text(json.dumps(done), encoding="utf-8")

    out = {
        "shape": args.shape,
        "center": table["center"],
        "std": table["std"],
        "lower-bound": table["lower-bound"],
        "upper-bound": table["upper-bound"],
        "kvals": table["kvals"],
        "lvals": ells,
        "cells": [done[str(l)] for l in ells],
        "convention": CONVENTION,
        "generator": "tests/tools/make_sbessel_truth_table.py",
        "target_rel": 1.0e-25,
        "negligible_below": NEGLIGIBLE,
    }

    with gzip.open(args.output, "wt") as handle:
        json.dump(out, handle)

    status = [c["status"] for row in out["cells"] for c in row]
    print(
        f"wrote {args.output}: {status.count('certified')} certified, "
        f"{status.count('negligible')} proven negligible, "
        f"{status.count('uncertified')} uncertified",
        file=sys.stderr,
    )

    return 1 if status.count("uncertified") else 0


if __name__ == "__main__":
    sys.exit(main())

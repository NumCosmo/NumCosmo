#!/usr/bin/env python
#
# make_xcor_pair_grid.py
#
# Tue Sep 8 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# make_xcor_pair_grid.py
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

"""Certify two windows on a *shared* grid of wavenumbers.

``make_xcor_window_truth_table.py`` lays each window's grid on its own turning
point, :math:`x_\\nu = k\\chi_\\mathrm{ref}/(\\ell + 1/2)`, which is right for
judging one window: it puts a certified value in each regime the solver
distinguishes. It is useless for judging a *pair*, because two windows at
different distances then carry certified values at different :math:`k`, and the
product that a cross spectrum integrates is never certified anywhere.

This tool fixes one grid and certifies both windows on it, so the product

.. math:: I^{(1)}_\\ell(k)\\, I^{(2)}_\\ell(k)

is a product of certified numbers rather than of one certified and one
interpolated. It is what lets the tail-against-tail figure state where a
far-separated pair's spectrum comes from, and how far a method is from the truth
*there* rather than at either window's peak.

The grid is logarithmic and deliberately spans both turning points, since the
point of a separated pair is that the two are far apart.

Emits ``data/truth_tables/xcor/xcor_pair_grid.json.gz``. Needs FLINT; reading it
does not.

Usage, and the command that reproduces the committed file::

    python3 tests/tools/make_xcor_pair_grid.py \\
        --pair gauss_near gauss_far --ells 10 50 --n-k 24
"""

from __future__ import annotations

import argparse
import concurrent.futures
import gzip
import json
import pathlib
import sys
import tempfile
import time

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "python" / "nc" / "xcor"))

import windows as W  # noqa: E402  pylint: disable=wrong-import-position

import make_xcor_window_truth_table as base  # noqa: E402  pylint: disable=wrong-import-position

OUTPUT = HERE.parent.parent / "data" / "truth_tables" / "xcor" / "xcor_pair_grid.json.gz"

CONVENTION = (
    "value[case][i_ell][i_k] = int W(chi) j_ell^(d)(kvals[i_ell][i_k] chi) dchi, "
    "with W the window named by case, normalized to unit integral over its "
    "truncated support. Both cases share the k grid at each multipole, so the "
    "product of two entries at the same index is a product of certified values. "
    "chi in Mpc, k in 1/Mpc. Values are full-precision decimal strings; radius "
    "is the certified absolute half-width, so |true - float(value)| <= radius."
)


def certify(exe: pathlib.Path, case: str, ell: int, ks, target_rel: float) -> dict:
    """One case at one multipole, on the supplied grid."""
    window = W.get(case)
    args = [
        str(exe),
        f"--ell={ell}",
        f"--target-rel={target_rel:g}",
        *W.arb_args(case),
    ]
    started = time.monotonic()
    out = base.subprocess.run(
        args,
        input=" ".join(repr(float(k)) for k in ks),
        capture_output=True,
        text=True,
        check=True,
    ).stdout
    vals, rads = [], []

    for line in out.splitlines():
        if line.startswith("#"):
            continue

        fields = line.split("\t")
        vals.append(fields[3])
        rads.append(float(fields[4]))

    if len(vals) != len(ks):
        sys.exit(f"{case} ell={ell}: {len(vals)} of {len(ks)} values")

    print(
        f"  {case:>18} ell={ell:<5} {len(ks)} k in "
        f"{time.monotonic() - started:6.1f} s",
        file=sys.stderr,
        flush=True,
    )

    return {"form": window.form, "deriv": window.deriv, "table": vals, "radius": rads}


def support_of(exe: pathlib.Path, case: str, ell: int, target_rel: float) -> tuple:
    """The window's support, read from the generator itself.

    Deliberately not routed through the truth-table generator's own helper: that
    one resolves a case through its ``CASES`` dict, which is keyed by *its*
    historical names, while everything here speaks the shared vocabulary of
    ``windows.py``. Reading the support from the tool's own output keeps one
    naming scheme in play.
    """
    out = base.subprocess.run(
        [
            str(exe),
            f"--ell={ell}",
            f"--target-rel={target_rel:g}",
            "--support-only",
            *W.arb_args(case),
        ],
        input="",
        capture_output=True,
        text=True,
        check=True,
    ).stdout

    for line in out.splitlines():
        if line.startswith("# shape="):
            lo, hi = line.split("support=[")[1].split("]")[0].split(",")

            return float(lo), float(hi)

    sys.exit(f"{case}: generator printed no support line")


def shared_grid(exe: pathlib.Path, pair, ell: int, n_k: int, target_rel: float):
    """A logarithmic grid spanning both windows' turning points, with margin.

    Each window's projection peaks near ``(ell + 1/2) / chi_ref``. A separated
    pair has two such peaks, and the product lives between and beyond them, so the
    grid brackets both rather than centring on either.
    """
    peaks = []

    for case in pair:
        lo, hi = support_of(exe, case, ell, target_rel)
        peaks.append((ell + 0.5) / (0.5 * (lo + hi)))

    lo, hi = min(peaks) / 4.0, max(peaks) * 8.0

    return [lo * (hi / lo) ** (i / (n_k - 1)) for i in range(n_k)]


def main() -> int:
    """Certify the pair on a shared grid and write the table."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pair", nargs=2, default=["gauss_near", "gauss_far"])
    parser.add_argument("--ells", nargs="*", type=int, default=[10, 50])
    parser.add_argument("--n-k", type=int, default=24)
    parser.add_argument("--target-rel", type=float, default=1.0e-25)
    parser.add_argument("--jobs", type=int, default=4)
    parser.add_argument("--output", type=pathlib.Path, default=OUTPUT)
    args = parser.parse_args()

    for case in args.pair:
        W.get(case)

    with tempfile.TemporaryDirectory() as tmp:
        exe = base.build_tool(pathlib.Path(tmp))
        grids = {
            ell: shared_grid(exe, args.pair, ell, args.n_k, args.target_rel)
            for ell in args.ells
        }
        cases: dict = {case: {} for case in args.pair}

        with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
            futures = {
                pool.submit(certify, exe, case, ell, grids[ell], args.target_rel): (
                    case,
                    ell,
                )
                for case in args.pair
                for ell in args.ells
            }

            for future in concurrent.futures.as_completed(futures):
                case, ell = futures[future]
                cases[case][ell] = future.result()

    payload = {
        "convention": CONVENTION,
        "generator": "tests/tools/make_xcor_pair_grid.py",
        "target_rel": args.target_rel,
        "pair": list(args.pair),
        "ells": list(args.ells),
        "kvals": [grids[ell] for ell in args.ells],
        "cases": {
            case: {
                "form": cases[case][args.ells[0]]["form"],
                "deriv": cases[case][args.ells[0]]["deriv"],
                "ctor": list(W.get(case).params),
                "table": [cases[case][ell]["table"] for ell in args.ells],
                "radius": [cases[case][ell]["radius"] for ell in args.ells],
            }
            for case in args.pair
        },
    }

    args.output.parent.mkdir(parents=True, exist_ok=True)

    with gzip.open(args.output, "wt") as handle:
        json.dump(payload, handle, indent=1, sort_keys=True)

    worst = max(
        radius / abs(float(value))
        for case in payload["cases"].values()
        for row_v, row_r in zip(case["table"], case["radius"])
        for value, radius in zip(row_v, row_r)
        if float(value) != 0.0
    )
    print(f"wrote {args.output} ({args.output.stat().st_size / 1024:.1f} KiB)")
    print(f"worst certified relative radius: {worst:.2e}")

    return 0


if __name__ == "__main__":
    sys.exit(main())

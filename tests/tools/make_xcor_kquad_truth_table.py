#!/usr/bin/env python3
#
# make_xcor_kquad_truth_table.py
#
# Fri Aug 28 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# make_xcor_kquad_truth_table.py
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

"""Generate the certified C_ell truth table by driving ``nc_xcor_kquad_arb``.

The pairs, windows and multipoles come from ``cases_k_integral.py``, so the
table covers exactly what the suite measures and cannot drift from it. Only
regeneration needs FLINT; the committed table is plain JSON and the test suite
reads it with no extra dependency.

    python3 tests/tools/make_xcor_kquad_truth_table.py \\
        --binary Optimized/tests/tools/nc_xcor_kquad_arb \\
        --out data/truth_tables/xcor/xcor_kquad.json.gz

Each entry carries its certified radius. A test that asserts against a value
whose accuracy is proved is doing something different from one that asserts
against a value someone believed.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import gzip
import json
import pathlib
import subprocess
import sys
import time

sys.path.insert(
    0, str(pathlib.Path(__file__).resolve().parents[1] / "python" / "nc" / "xcor")
)

import cases_k_integral as cases  # noqa: E402  pylint: disable=wrong-import-position

# The k range each pair is integrated over, before the tail is bounded rather
# than integrated. Below k_lo and above k_hi the integrand is bounded, not
# dropped -- see nc_xcor_kquad_arb.c. k_hi is per-pair because the cost of a
# panel grows with the phase k chi_max, and integrating far past where the
# windows have decayed is where all the expensive precision went.
K_LO = 1.0e-6
K_END = 1.0e4


def k_hi_for(pair: cases.PairSpec) -> float:
    """Where to stop integrating and start bounding.

    Past k ~ 30 / (smallest feature in chi) every window here has decayed, and
    the phase there is large enough that the 0F1 form of j_ell needs thousands
    of bits to confirm it. The bound costs nothing and says the same thing.
    """
    scales = []

    for name in (pair.kernel_a, pair.kernel_b):
        params = cases.KERNELS[name].params

        if "chi-sigma" in params:
            scales.append(float(params["chi-sigma"]))
        elif "chi-scale" in params:
            scales.append(float(params["chi-scale"]))
        elif "sigma" in params:
            scales.extend(float(v) for v in params["sigma"])
        else:
            # A hard edge has no smoothing scale at all: its I_ell decays only
            # as 1/k, so it needs a wider range than any smooth window.
            lower = float(params.get("chi-lower", 100.0))
            upper = float(params.get("chi-upper", 3000.0))
            scales.append(0.02 * (upper - lower))

    return min(30.0 / min(scales), 20.0)


def run_pair(
    binary: str, pair: cases.PairSpec, ells: list[int], target: float, timeout: float
) -> dict:
    """One pair, all multipoles, through the Arb generator."""
    spec_a = cases.KERNELS[pair.kernel_a]
    spec_b = cases.KERNELS[pair.kernel_b]
    k_hi = k_hi_for(pair)

    args = [binary]
    args += spec_a.arb_args("a")

    if not pair.isauto:
        args += spec_b.arb_args("b")

    args += [
        "--ells=" + ",".join(str(e) for e in ells),
        f"--k-lo={K_LO:.17g}",
        f"--k-hi={k_hi:.17g}",
        f"--k-end={K_END:.17g}",
        f"--target-rel={target:.17g}",
    ]

    start = time.perf_counter()

    try:
        out = subprocess.run(
            args, capture_output=True, text=True, timeout=timeout, check=False
        )
    except subprocess.TimeoutExpired:
        return {"case": pair.case, "failed": f"timed out after {timeout:.0f} s"}

    rows = [
        line.split("\t") for line in out.stdout.splitlines() if not line.startswith("#")
    ]

    if len(rows) != len(ells):
        return {
            "case": pair.case,
            "failed": f"{len(rows)} of {len(ells)} multipoles",
            "stderr": out.stderr[-500:],
        }

    return {
        "case": pair.case,
        "kernel_a": pair.kernel_a,
        "kernel_b": pair.kernel_b,
        "isauto": pair.isauto,
        "k_lo": K_LO,
        "k_hi": k_hi,
        "k_end": K_END,
        "ells": ells,
        "values": [r[1] for r in rows],
        "radius": [float(r[2]) for r in rows],
        "prec": [int(r[3]) for r in rows],
        "seconds": time.perf_counter() - start,
    }


def main() -> None:
    """Drive the generator over the case matrix."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", default="Optimized/tests/tools/nc_xcor_kquad_arb")
    parser.add_argument("--out", type=pathlib.Path, required=True)
    parser.add_argument("--ells", nargs="*", type=int, default=[2, 10, 50, 200])
    parser.add_argument("--cases", nargs="*", default=None)
    parser.add_argument("--target-rel", type=float, default=1.0e-10)
    parser.add_argument("--timeout", type=float, default=3600.0)
    parser.add_argument("--jobs", type=int, default=6)
    args = parser.parse_args()

    pairs = (
        cases.PAIRS
        if args.cases is None
        else [cases.PAIRS_BY_CASE[c] for c in args.cases]
    )

    entries: dict[str, dict] = {}
    failed: list[str] = []

    with concurrent.futures.ProcessPoolExecutor(max_workers=args.jobs) as pool:
        futures = {
            pool.submit(
                run_pair, args.binary, pair, args.ells, args.target_rel, args.timeout
            ): pair
            for pair in pairs
        }

        for future in concurrent.futures.as_completed(futures):
            entry = future.result()

            if "failed" in entry:
                failed.append(f"{entry['case']}: {entry['failed']}")
                print(f"  {entry['case']:4s} FAILED {entry['failed']}", flush=True)
                continue

            entries[entry["case"]] = entry
            print(
                f"  {entry['case']:4s} ok  {entry['seconds']:7.1f} s  "
                f"worst radius {max(entry['radius']):.2e}",
                flush=True,
            )

    payload = {
        "convention": (
            "C_ell = 2/pi INT dk k^2 P(k) I1_ell(k) I2_ell(k), k in 1/Mpc, "
            "chi in Mpc, each window normalized to unit integral over its "
            "truncated support, P the NcmPowspecAnalytic BBKS spectrum at z=0"
        ),
        "powspec": {
            "shape": "bbks",
            "amplitude": 2.08e7,
            "n_s": 0.9875,
            "k_eq": 0.10594,
        },
        "target_rel": args.target_rel,
        "ells": args.ells,
        "generator": "tests/tools/nc_xcor_kquad_arb.c",
        "cases": entries,
    }

    args.out.parent.mkdir(parents=True, exist_ok=True)

    with gzip.open(args.out, "wt") as handle:
        json.dump(payload, handle, indent=1)

    print(f"\nwrote {len(entries)} of {len(pairs)} cases to {args.out}")

    if failed:
        print("failed:")

        for line in failed:
            print(f"  {line}")


if __name__ == "__main__":
    main()

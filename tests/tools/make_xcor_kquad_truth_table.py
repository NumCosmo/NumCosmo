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

# The k range each pair is integrated over.
#
# Taken from the library's own closure range, not chosen. Both sides truncate:
# the library integrates over the intersection of the two closures' fitted
# domains, so a reference integrating further would differ from it by a tail
# neither side claims to compute, and the comparison would measure that tail
# instead of the quadrature.
#
# The truncation is worth reporting rather than hiding -- it is ~1e-8 relative
# on a 4 sigma Gaussian, far above the quadrature's own radius, and it is set
# by the window's n_sigma cut rather than its width: a hard edge makes I_ell
# fall off as 1/k instead of exponentially. Closing that gap would need k_hi a
# hundred times larger, where the phase k chi puts the 0F1 form of j_ell out of
# reach, so it is stated instead.
K_MARGIN = 1.5


def library_k_range(pair: cases.PairSpec, ell: int) -> tuple[float, float]:
    """The k interval the library's own closures span at ONE multipole, in 1/Mpc.

    Per multipole, emphatically. The closure's range tracks the Bessel turning
    point k ~ ell / chi, so it moves up by more than a decade between ell = 2
    and ell = 200. Taking the union over multipoles -- which this did at first
    -- integrates ell = 2 out to where only ell = 200 has support, which is
    exactly where the 0F1 form of j_ell needs thousands of bits and contributes
    nothing. It turned a 45 second case into one that had not finished in 100
    minutes.
    """
    from numcosmo_py import Nc, Ncm  # noqa: PLC0415

    Ncm.cfg_init()
    cosmo, dist, ps = cases.make_cosmo_bits()
    RH = Nc.HICosmo.RH_Mpc(cosmo)
    settings = cases.Settings(closure=Nc.XcorKernelClosure.CHEBYSHEV)

    kernel_a = cases.build_kernel(pair.kernel_a, cosmo, dist, ps, settings)
    kernel_b = (
        kernel_a
        if pair.isauto
        else cases.build_kernel(pair.kernel_b, cosmo, dist, ps, settings)
    )
    a_lo, a_hi = cases.build_integrand(kernel_a, cosmo, ell, ell, settings).get_range()

    if pair.isauto:
        b_lo, b_hi = a_lo, a_hi
    else:
        b_lo, b_hi = cases.build_integrand(
            kernel_b, cosmo, ell, ell, settings
        ).get_range()

    return min(a_lo, b_lo) / RH / K_MARGIN, max(a_hi, b_hi) / RH * K_MARGIN


def run_one(
    binary: str, pair: cases.PairSpec, ell: int, target: float, timeout: float
) -> dict:
    """One pair at one multipole: the unit of work, and of what survives.

    A whole pair was the unit at first, and a pair that stalls on its highest
    multipole then discards the multipoles that finished in under a minute.
    """
    spec_a = cases.KERNELS[pair.kernel_a]
    spec_b = cases.KERNELS[pair.kernel_b]
    k_lo, k_hi = library_k_range(pair, ell)

    args = [binary] + spec_a.arb_args("a")

    if not pair.isauto:
        args += spec_b.arb_args("b")

    args += [
        f"--ells={ell}",
        f"--k-lo={k_lo:.17g}",
        f"--k-hi={k_hi:.17g}",
        f"--target-rel={target:.17g}",
    ]

    start = time.perf_counter()

    try:
        out = subprocess.run(
            args, capture_output=True, text=True, timeout=timeout, check=False
        )
    except subprocess.TimeoutExpired:
        return {
            "case": pair.case,
            "ell": ell,
            "failed": f"timed out after {timeout:.0f} s",
            "k_lo": k_lo,
            "k_hi": k_hi,
        }

    rows = [
        line.split("\t") for line in out.stdout.splitlines() if not line.startswith("#")
    ]

    if len(rows) != 1:
        return {
            "case": pair.case,
            "ell": ell,
            "failed": "no row",
            "stderr": out.stderr[-400:],
        }

    value = float(rows[0][1])
    radius = float(rows[0][2])

    # Exiting cleanly is not the same as certifying. A run that ends with the
    # enclosure straddling zero has produced no significant digit, and recorded
    # as a value it would be a wrong entry in a table whose whole purpose is to
    # be trustworthy. Rejected as loudly as a timeout.
    if not value or radius > 1.0e-8 * abs(value):
        return {
            "case": pair.case,
            "ell": ell,
            "failed": f"radius {radius:.2e} against value {value:.4e}",
            "k_lo": k_lo,
            "k_hi": k_hi,
        }

    return {
        "case": pair.case,
        "ell": ell,
        "kernel_a": pair.kernel_a,
        "kernel_b": pair.kernel_b,
        "isauto": pair.isauto,
        "k_lo": k_lo,
        "k_hi": k_hi,
        "value": rows[0][1],
        "radius": radius,
        "prec": int(rows[0][3]),
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
    parser.add_argument(
        "--resume", action="store_true", help="skip entries already in --out"
    )
    args = parser.parse_args()

    pairs = (
        cases.PAIRS
        if args.cases is None
        else [cases.PAIRS_BY_CASE[c] for c in args.cases]
    )

    entries: dict[str, dict] = {}
    failed: list[str] = []

    # Resuming matters here: an entry is up to hours of arbitrary-precision
    # arithmetic, and the expensive ones are exactly the ones a first pass
    # fails to finish. Without this, extending the time limit means recomputing
    # everything that already succeeded.
    if args.resume and args.out.exists():
        with gzip.open(args.out, "rt") as handle:
            entries = json.load(handle)["cases"]

        print(f"resuming: {len(entries)} entries already certified", flush=True)

    tasks = [
        (pair, ell)
        for pair in pairs
        for ell in args.ells
        if f"{pair.case}_l{ell}" not in entries
    ]

    print(f"{len(tasks)} tasks over {len(pairs)} pairs", flush=True)

    with concurrent.futures.ProcessPoolExecutor(max_workers=args.jobs) as pool:
        futures = [
            pool.submit(run_one, args.binary, pair, ell, args.target_rel, args.timeout)
            for pair, ell in tasks
        ]

        for future in concurrent.futures.as_completed(futures):
            entry = future.result()
            key = f"{entry['case']}_l{entry['ell']}"

            if "failed" in entry:
                failed.append(f"{key}: {entry['failed']}")
                print(f"  {key:12s} FAILED {entry['failed']}", flush=True)
                continue

            entries[key] = entry
            print(
                f"  {key:12s} ok  {entry['seconds']:7.1f} s  prec {entry['prec']:5d}  "
                f"k<={entry['k_hi']:.3g}  radius {entry['radius']:.2e}",
                flush=True,
            )
            write_table(args.out, entries, args)

    write_table(args.out, entries, args)

    print(f"\nwrote {len(entries)} of {len(tasks)} entries to {args.out}")

    if failed:
        print("failed:")

        for line in failed:
            print(f"  {line}")


def write_table(out: pathlib.Path, entries: dict, args) -> None:
    """Serialize what is finished so far."""
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

    out.parent.mkdir(parents=True, exist_ok=True)

    with gzip.open(out, "wt") as handle:
        json.dump(payload, handle, indent=1)


if __name__ == "__main__":
    main()

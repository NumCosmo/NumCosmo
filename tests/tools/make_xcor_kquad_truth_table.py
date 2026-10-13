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

An entry one run cannot finish within ``--timeout`` is certified in parts of its
k range with ``--k-chunks N``, in parallel, and the parts are summed as balls;
a part that times out is re-run as ``--resplit M`` sub-parts, and a sub-part as
``--resplit-sub L`` pieces, each with ``--resume``. The finished parts live in
``<out>.chunks.json`` beside the table until the entry is merged, and a merged
entry records its parts under ``chunks``, so the split it was made with is in
the table.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import decimal
import gzip
import json
import pathlib
import subprocess
import sys
import time

sys.path.insert(
    0, str(pathlib.Path(__file__).resolve().parents[1] / "python" / "nc" / "xcor")
)

import cases_k_integral as cases  # noqa: E402

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

# Cases where K_MARGIN is not enough, as an absolute k_hi in 1/Mpc.
#
# The lensing window has a hard edge at chi_lower where it does not vanish, so
# its I_ell falls as 1/k and the C_ell integrand only as k^-3.5. Beyond the
# library's ell = 2 range times the margin (0.86/Mpc), the auto spectrum has
# 2.9e-7 of its own value still outstanding, against a target of 1e-10.
# Measured on the library's own closure built out to 79/Mpc: 1.8e-8 outstanding
# beyond 2/Mpc, 7.5e-10 beyond 5/Mpc, 5.0e-11 beyond 10/Mpc. The cross spectrum
# with a Gaussian needs no floor -- the Gaussian side cuts the product off, at
# 4e-10 beyond that same 0.86/Mpc. Cost grows as the square of the phase
# k_hi chi_max, which is why this is one entry and not a global bound.
K_HI_FLOOR = {"A4": 10.0}


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
    settings = cases.Settings()

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

    k_hi = max(a_hi, b_hi) / RH * K_MARGIN

    return min(a_lo, b_lo) / RH / K_MARGIN, max(k_hi, K_HI_FLOOR.get(pair.case, 0.0))


def library_cl(pair: cases.PairSpec, ell: int) -> float:
    """The library's C_ell of one pair at ONE multipole, at its default settings.

    The scale a chunked entry certifies against: every chunk is held to the target
    relative to the whole C_ell, as one run over the whole range would be.
    """
    from numcosmo_py import Nc, Ncm  # noqa: PLC0415

    Ncm.cfg_init()
    cosmo, dist, ps = cases.make_cosmo_bits()
    settings = cases.Settings()

    kernel_a = cases.build_kernel(pair.kernel_a, cosmo, dist, ps, settings)
    kernel_b = (
        None
        if pair.isauto
        else cases.build_kernel(pair.kernel_b, cosmo, dist, ps, settings)
    )
    xcor = Nc.Xcor.new(dist, ps, Nc.XcorMethod.KERNEL_EXACT)
    xcor.set_ell_batch_size(1)
    xcor.prepare(cosmo)

    vp = Ncm.Vector.new(1)
    xcor.compute(kernel_a, kernel_b, cosmo, ell, ell, vp)

    return vp.get(0)


def _generator_args(
    binary: str,
    pair: cases.PairSpec,
    ell: int,
    target: float,
    k_lo: float,
    k_hi: float,
    scale: float | None = None,
) -> list[str]:
    """The generator's command line for one pair at one multipole over [k_lo, k_hi]."""
    args = [binary] + cases.KERNELS[pair.kernel_a].arb_args("a")

    if not pair.isauto:
        args += cases.KERNELS[pair.kernel_b].arb_args("b")

    args += [
        f"--ells={ell}",
        f"--k-lo={k_lo:.17g}",
        f"--k-hi={k_hi:.17g}",
        f"--target-rel={target:.17g}",
    ]

    if scale is not None:
        args.append(f"--scale={scale:.17g}")

    return args


def run_chunk(
    binary: str,
    pair: cases.PairSpec,
    ell: int,
    target: float,
    timeout: float,
    k_lo: float,
    k_hi: float,
    scale: float,
    index: int | str,
) -> dict:
    """One part [k_lo, k_hi] of a pair's k range, certified against the whole C_ell.

    No radius check here: a part may be orders below the whole and need not carry a
    significant digit of its own; the check is made on the sum, in merge_chunks().
    """
    start = time.perf_counter()
    args = _generator_args(binary, pair, ell, target, k_lo, k_hi, scale)
    out = {"case": pair.case, "ell": ell, "index": index, "k_lo": k_lo, "k_hi": k_hi}

    try:
        run = subprocess.run(
            args, capture_output=True, text=True, timeout=timeout, check=False
        )
    except subprocess.TimeoutExpired:
        return out | {"failed": f"timed out after {timeout:.0f} s"}

    rows = [
        line.split("\t") for line in run.stdout.splitlines() if not line.startswith("#")
    ]

    if len(rows) != 1:
        return out | {"failed": "no row", "stderr": run.stderr[-400:]}

    return out | {
        "value": rows[0][1],
        "radius": float(rows[0][2]),
        "prec": int(rows[0][3]),
        "seconds": time.perf_counter() - start,
    }


def _sum_balls(chunks: list[dict]) -> tuple[decimal.Decimal, float]:
    """The sum of certified balls: values added in decimal, radii added.

    Each value is printed with only the digits its radius supports, so half a unit of its
    last printed digit is added to its radius: the sum of the printed values is then still
    inside the sum of the balls.
    """
    context = decimal.Context(prec=50)
    total = decimal.Decimal(0)
    radius = 0.0

    for chunk in chunks:
        value = decimal.Decimal(chunk["value"])
        total = context.add(total, value)
        radius += chunk["radius"] + 0.5 * 10.0 ** value.as_tuple().exponent

    return total, radius


def merge_subparts(index: int | str, subparts: list[dict]) -> dict:
    """Part @index of an entry from the certified sub-parts it was re-split into."""
    subparts = sorted(subparts, key=lambda c: int(str(c["index"]).split(".")[-1]))
    total, radius = _sum_balls(subparts)

    return {
        "case": subparts[0]["case"],
        "ell": subparts[0]["ell"],
        "index": index,
        "k_lo": subparts[0]["k_lo"],
        "k_hi": subparts[-1]["k_hi"],
        "value": f"{total:e}",
        "radius": radius,
        "prec": max(c["prec"] for c in subparts),
        "seconds": sum(c["seconds"] for c in subparts),
        "subparts": subparts,
    }


def merge_chunks(pair: cases.PairSpec, ell: int, chunks: list[dict]) -> dict:
    """The entry of a pair at one multipole from its certified parts, summed as
    _sum_balls() does, and held to the same certification as one run, in run_one().
    """
    chunks = sorted(chunks, key=lambda c: c["index"])
    total, radius = _sum_balls(chunks)

    entry = {
        "case": pair.case,
        "ell": ell,
        "kernel_a": pair.kernel_a,
        "kernel_b": pair.kernel_b,
        "isauto": pair.isauto,
        "k_lo": chunks[0]["k_lo"],
        "k_hi": chunks[-1]["k_hi"],
        "value": f"{total:e}",
        "radius": radius,
        "prec": max(c["prec"] for c in chunks),
        "seconds": sum(c["seconds"] for c in chunks),
        "chunks": chunks,
    }

    if not total or radius > 1.0e-8 * abs(float(total)):
        return entry | {
            "failed": f"radius {radius:.2e} against value {float(total):.4e}"
        }

    return entry


def run_one(
    binary: str, pair: cases.PairSpec, ell: int, target: float, timeout: float
) -> dict:
    """One pair at one multipole: the unit of work, and of what survives.

    A whole pair was the unit at first, and a pair that stalls on its highest
    multipole then discards the multipoles that finished in under a minute.
    """
    k_lo, k_hi = library_k_range(pair, ell)
    args = _generator_args(binary, pair, ell, target, k_lo, k_hi)

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
    parser.add_argument(
        "--k-chunks",
        type=int,
        default=1,
        help=(
            "split each entry's k range into this many parts, uniform in k (the "
            "phase k chi is linear in it), certified in parallel and summed; for "
            "entries one run cannot finish within --timeout"
        ),
    )
    parser.add_argument(
        "--resplit",
        type=int,
        default=1,
        help=(
            "with --k-chunks and --resume, certify each part still missing (one that "
            "timed out) as this many sub-parts of its own interval"
        ),
    )
    parser.add_argument(
        "--resplit-sub",
        type=int,
        default=1,
        help=(
            "with --resplit, certify each sub-part still missing as this many pieces "
            "of its own interval"
        ),
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

    if args.k_chunks > 1:
        run_chunked(args, tasks, entries, failed)
        report(args, entries, tasks, failed)

        return

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
    report(args, entries, tasks, failed)


def _entry_range(
    stored: dict, pair: cases.PairSpec, ell: int, n_chunks: int
) -> tuple[float, float, float]:
    """The k range an entry is split over and the scale its parts certify against.

    Fixed when the entry is first split and kept with its parts: the library's own range
    moves whenever the closures' domains do, and parts of two different splits do not sum
    to the integral. Parts saved before the range was kept give it back, since they tile the
    range uniformly.
    """
    if "_range" not in stored:
        done = [v for k, v in stored.items() if "." not in k and k.isdigit()]

        if done:
            p = done[0]
            width = (p["k_hi"] - p["k_lo"]) * n_chunks
            k_lo = p["k_lo"] - (p["k_hi"] - p["k_lo"]) * p["index"]
            k_range = (k_lo, k_lo + width)
        else:
            k_range = library_k_range(pair, ell)

        stored["_range"] = {
            "k_lo": k_range[0],
            "k_hi": k_range[1],
            "scale": abs(library_cl(pair, ell)),
        }

    r = stored["_range"]

    return r["k_lo"], r["k_hi"], r["scale"]


def _check_parts(key: str, stored: dict, edges: list[float], n_chunks: int) -> None:
    """Abort when a stored part does not sit on the edges it is summed between."""
    for name, part in stored.items():
        if name == "_range":
            continue

        i = int(name.split(".")[0])
        lo, hi = edges[i], edges[i + 1]

        if "." in name:
            continue

        if (
            abs(part["k_lo"] - lo) > 1.0e-12 * hi
            or abs(part["k_hi"] - hi) > 1.0e-12 * hi
        ):
            raise ValueError(
                f"{key}: part {name} covers [{part['k_lo']}, {part['k_hi']}], "
                f"the split puts it on [{lo}, {hi}]"
            )


def run_chunked(args, tasks: list, entries: dict, failed: list[str]) -> None:
    """Certify each entry as --k-chunks parts in parallel, and merge them as they complete.

    Finished parts are kept beside --out, so a rerun with --resume certifies only the
    parts still missing. With --resplit M a missing part is certified as M sub-parts of its
    own interval instead, labelled "<part>.<sub>", and summed into the part once all are in;
    with --resplit-sub N a missing sub-part is certified as N pieces "<part>.<sub>.<piece>",
    summed into the sub-part the same way.
    """
    parts_path = args.out.parent / (args.out.name + ".chunks.json")
    parts: dict[str, dict] = {}

    if args.resume and parts_path.exists():
        parts = json.loads(parts_path.read_text())

    jobs = []

    for pair, ell in tasks:
        key = f"{pair.case}_l{ell}"
        parts.setdefault(key, {})
        k_lo, k_hi, scale = _entry_range(parts[key], pair, ell, args.k_chunks)
        edges = [
            k_lo + (k_hi - k_lo) * i / args.k_chunks for i in range(args.k_chunks + 1)
        ]
        _check_parts(key, parts[key], edges, args.k_chunks)

        for i in range(args.k_chunks):
            if str(i) in parts[key]:
                continue

            if args.resplit <= 1:
                jobs.append((pair, ell, edges[i], edges[i + 1], scale, i))
                continue

            for j in range(args.resplit):
                if f"{i}.{j}" in parts[key]:
                    continue

                lo = edges[i] + (edges[i + 1] - edges[i]) * j / args.resplit
                hi = edges[i] + (edges[i + 1] - edges[i]) * (j + 1) / args.resplit

                if args.resplit_sub <= 1:
                    jobs.append((pair, ell, lo, hi, scale, f"{i}.{j}"))
                    continue

                for m in range(args.resplit_sub):
                    if f"{i}.{j}.{m}" not in parts[key]:
                        jobs.append(
                            (
                                pair,
                                ell,
                                lo + (hi - lo) * m / args.resplit_sub,
                                lo + (hi - lo) * (m + 1) / args.resplit_sub,
                                scale,
                                f"{i}.{j}.{m}",
                            )
                        )

    print(f"{len(jobs)} parts to certify", flush=True)

    with concurrent.futures.ProcessPoolExecutor(max_workers=args.jobs) as pool:
        futures = {
            pool.submit(
                run_chunk,
                args.binary,
                pair,
                ell,
                args.target_rel,
                args.timeout,
                lo,
                hi,
                scale,
                i,
            ): (pair, ell)
            for pair, ell, lo, hi, scale, i in jobs
        }

        for future in concurrent.futures.as_completed(futures):
            pair, ell = futures[future]
            chunk = future.result()
            key = f"{pair.case}_l{ell}"

            if "failed" in chunk:
                failed.append(f"{key} part {chunk['index']}: {chunk['failed']}")
                print(
                    f"  {key:12s} part {str(chunk['index']):>5s} FAILED {chunk['failed']}",
                    flush=True,
                )
                continue

            parts[key][str(chunk["index"])] = chunk
            print(
                f"  {key:12s} part {str(chunk['index']):>5s} ok  {chunk['seconds']:7.1f} s  "
                f"k in [{chunk['k_lo']:.3g}, {chunk['k_hi']:.3g}]",
                flush=True,
            )

            # A re-split sub-part is whole once all its pieces are in, and a re-split
            # part once all its sub-parts are.
            if isinstance(chunk["index"], str) and chunk["index"].count(".") == 2:
                sub = chunk["index"].rsplit(".", 1)[0]
                pieces = [
                    parts[key][f"{sub}.{m}"]
                    for m in range(args.resplit_sub)
                    if f"{sub}.{m}" in parts[key]
                ]

                if len(pieces) == args.resplit_sub:
                    parts[key][sub] = merge_subparts(sub, pieces)
                    for m in range(args.resplit_sub):
                        del parts[key][f"{sub}.{m}"]
                    print(
                        f"  {key:12s} part {sub:>5s} ok  in {args.resplit_sub} pieces",
                        flush=True,
                    )
                    chunk = parts[key][sub]

            if isinstance(chunk["index"], str) and chunk["index"].count(".") == 1:
                i = int(chunk["index"].split(".")[0])
                subs = [
                    parts[key][f"{i}.{j}"]
                    for j in range(args.resplit)
                    if f"{i}.{j}" in parts[key]
                ]

                if len(subs) == args.resplit:
                    parts[key][str(i)] = merge_subparts(i, subs)
                    for j in range(args.resplit):
                        del parts[key][f"{i}.{j}"]
                    print(
                        f"  {key:12s} part {i:5d} ok  in {args.resplit} sub-parts",
                        flush=True,
                    )

            parts_path.write_text(json.dumps(parts, indent=1))

            if all(str(i) in parts[key] for i in range(args.k_chunks)):
                entry = merge_chunks(
                    pair, ell, [parts[key][str(i)] for i in range(args.k_chunks)]
                )

                if "failed" in entry:
                    failed.append(f"{key}: {entry['failed']}")
                    print(f"  {key:12s} FAILED {entry['failed']}", flush=True)
                    continue

                entries[key] = entry
                print(
                    f"  {key:12s} ok  {entry['seconds']:7.1f} s in {args.k_chunks} parts  "
                    f"radius {entry['radius']:.2e}",
                    flush=True,
                )
                write_table(args.out, entries, args)

    write_table(args.out, entries, args)


def report(args, entries: dict, tasks: list, failed: list[str]) -> None:
    """What this run certified, against the whole table."""
    # `entries` is the whole table, resumed cells included, while `tasks` is
    # only this run's work: reporting one "of" the other read as "87 of 1".
    print(
        f"\nwrote {len(entries)} entries to {args.out}\n"
        f"  {len(tasks) - len(failed)} of {len(tasks)} cells certified this run"
    )

    if failed:
        print("failed:")

        for line in failed:
            print(f"  {line}")


def write_table(out: pathlib.Path, entries: dict, args) -> None:
    """Serialize what is finished so far.

    The ``ells`` header is derived from the entries rather than taken from the
    invocation: a partial rerun (one case, a subset of multipoles) must not
    relabel a table that still carries every other case.
    """
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
        "ells": sorted({entry["ell"] for entry in entries.values()}),
        "generator": "tests/tools/nc_xcor_kquad_arb.c",
        "cases": entries,
    }

    out.parent.mkdir(parents=True, exist_ok=True)

    with gzip.open(out, "wt") as handle:
        json.dump(payload, handle, indent=1)


if __name__ == "__main__":
    main()

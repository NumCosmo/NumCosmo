#!/usr/bin/env python3
"""Summarize or compare files from xcor_lsst_benchmark.py.

    python dev-notes/xcor_lsst_compare.py /tmp/xcor-before.json
    python dev-notes/xcor_lsst_compare.py /tmp/xcor-before.json /tmp/xcor-after.json \
        --skip-first

Timings are means/medians over repetitions of each multipole block. Ratios are
before/after (greater than one means faster). Numerical changes are normalized
by each reference closure or spectrum's maximum absolute value, so zeros and
sign changes do not create misleading pointwise relative errors. All repetitions
are checked numerically, including any omitted from the timing summary.
"""

import argparse
import json
from pathlib import Path
import statistics

import numpy as np


def groups(data: dict, skip_first: bool = False) -> dict:
    """Group repetitions of each ell block, optionally excluding repetition zero."""
    grouped = {}
    for block in data["blocks"]:
        if skip_first and block["repeat"] == 0:
            continue
        grouped.setdefault((block["lmin"], block["lmax"]), []).append(block)
    return grouped


def summary(blocks: list[dict], kernel_idx: int | None) -> tuple[float, float]:
    """Return mean and median closure or combined outer integration seconds."""
    values = [block["kernels"][kernel_idx]["closure_s"] if kernel_idx is not None
              else sum(spectrum["outer_s"] for spectrum in block["spectra"])
              for block in blocks]
    return statistics.mean(values), statistics.median(values)


def scaled_difference(before: list, after: list) -> tuple[float, bool]:
    """Compute a maximum difference relative to the reference peak."""
    lhs, rhs = np.asarray(before), np.asarray(after)
    if lhs.shape != rhs.shape or not np.isfinite(lhs).all() or not np.isfinite(rhs).all():
        raise ValueError("Incompatible shapes or nonfinite numerical results")
    exact = bool(np.array_equal(lhs, rhs))
    scale = float(np.max(np.abs(lhs)))
    change = float(np.max(np.abs(lhs - rhs)))
    return (change / scale if scale else (0.0 if change == 0 else float("inf"))), exact


def main() -> None:
    """Print timing summaries and, when requested, numerical regression changes."""
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("before", type=Path)
    parser.add_argument("after", type=Path, nargs="?")
    parser.add_argument("--skip-first", action="store_true",
                        help="Exclude repetition zero from timing statistics only.")
    args = parser.parse_args()
    before = json.loads(args.before.read_text())
    after = json.loads(args.after.read_text()) if args.after else None
    if after:
        for key in ("specs", "integrator_reltol", "integrator_cheb_reltol",
                    "closure_reltol", "closure_peak_epsilon"):
            if before[key] != after[key]:
                raise ValueError(f"Benchmark configurations differ: {key}")
    before_groups = groups(before, args.skip_first)
    after_groups = groups(after, args.skip_first) if after else {}
    if not before_groups:
        parser.error("No repetitions remain; omit --skip-first or benchmark more repeats")
    if after and before_groups.keys() != after_groups.keys():
        parser.error("Before/after multipole blocks differ")
    for ell_range, blocks in before_groups.items():
        print(f"ell={ell_range[0]}..{ell_range[1]} (seconds; mean / median)")
        labels = [kernel["label"] for kernel in blocks[0]["kernels"]] + ["Outer, all pairs"]
        for idx, label in enumerate(labels):
            kernel_idx = idx if idx < len(labels) - 1 else None
            bmean, bmedian = summary(blocks, kernel_idx)
            line = f"  {label}: {bmean:.6f} / {bmedian:.6f}"
            if after:
                amean, amedian = summary(after_groups[ell_range], kernel_idx)
                line += (f" -> {amean:.6f} / {amedian:.6f}"
                         f"  ratio={bmean / amean:.3f}x / {bmedian / amedian:.3f}x")
            print(line)
    print(f"Peak RSS before: {before['peak_rss_kib'] / 1024:.1f} MiB")
    if not after:
        return
    print(f"Peak RSS after:  {after['peak_rss_kib'] / 1024:.1f} MiB")
    reference_groups = groups(before)
    all_exact = True
    closure_change = 0.0
    spectrum_change = 0.0
    domains_equal = True
    for block in after["blocks"]:
        reference = reference_groups[(block["lmin"], block["lmax"])][0]
        if len(reference["kernels"]) != len(block["kernels"]):
            raise ValueError("Kernel counts differ")
        if len(reference["spectra"]) != len(block["spectra"]):
            raise ValueError("Spectrum counts differ")
        for lhs, rhs in zip(reference["kernels"], block["kernels"]):
            if lhs["k"] != rhs["k"]:
                raise ValueError("Sample grids differ; rerun after with --reference before.json")
            change, exact = scaled_difference(lhs["values"], rhs["values"])
            closure_change = max(closure_change, change)
            all_exact &= exact
            domains_equal &= all(lhs[key] == rhs[key]
                                 for key in ("range", "component_ranges", "panels"))
        for lhs, rhs in zip(reference["spectra"], block["spectra"]):
            if lhs["pair"] != rhs["pair"]:
                raise ValueError("Spectrum pairs differ")
            change, exact = scaled_difference(lhs["cl"], rhs["cl"])
            spectrum_change = max(spectrum_change, change)
            all_exact &= exact
    print(f"Largest closure change / reference peak: {closure_change:.6e}")
    print(f"Largest C_ell change / reference peak:   {spectrum_change:.6e}")
    print(f"All numerical samples identical: {all_exact}")
    print(f"All closure domains and panel counts identical: {domains_equal}")


if __name__ == "__main__":
    main()

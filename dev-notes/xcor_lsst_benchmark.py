#!/usr/bin/env python3
"""Time the LSST chi closures separately from the outer k integration.

Run from a NumCosmo build environment (source its numcosmo_export.sh first):

    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python dev-notes/xcor_lsst_benchmark.py \
        --starts 2,98,498 --repeat 7 --output /tmp/xcor-before.json
    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python dev-notes/xcor_lsst_benchmark.py \
        --starts 2,98,498 --repeat 7 --reference /tmp/xcor-before.json \
        --output /tmp/xcor-after.json
    python dev-notes/xcor_lsst_compare.py /tmp/xcor-before.json /tmp/xcor-after.json

The defaults match the LSST-Y1 bin-4 kernels in ``numcosmo xcor cls``:
full non-Limber, Chebyshev closures, eight multipoles per block, exact outer
quadrature, and unchanged library tolerances. Cosmology/kernel preparation and
integrator ell-range changes are timed separately from closure construction.
One integrator is shared between kernels and blocks, so this measures the
closure workload; it does not measure the solver's per-block integrator memory.

Reference files preserve k grids as well as closure samples and C_ell values.
The samples are a regression comparison, not an independent accuracy certificate.
Use a fixed CPU affinity and alternate before/after runs on an otherwise idle
machine for useful timings. Repeats retain the integrator's available caches;
the comparison helper can omit each block's first repetition with --skip-first.
"""

import argparse
import contextlib
import io
import json
import os
from pathlib import Path
import resource
import time

import numpy as np


def positive_int(value: str) -> int:
    """Parse a positive block size or repetition count."""
    number = int(value)
    if number < 1:
        raise argparse.ArgumentTypeError("must be positive")
    return number


def parse_args() -> argparse.Namespace:
    """Read the benchmark parameters without initializing NumCosmo."""
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--starts", default="2,98,498",
                        help="Comma-separated first multipoles of the blocks.")
    parser.add_argument("--size", type=positive_int, default=8)
    parser.add_argument("--repeat", type=positive_int, default=1)
    parser.add_argument("--output", type=Path, default=Path("/tmp/xcor-benchmark.json"))
    parser.add_argument("--reference", type=Path,
                        help="Evaluate closures on the saved reference k grids.")
    parser.add_argument("--kernels", default="number-counts,weak-lensing",
                        help="Comma-separated selection of the two LSST kernel types.")
    args = parser.parse_args()
    try:
        args.starts = [int(value) for value in args.starts.split(",")]
    except ValueError:
        parser.error("--starts must contain comma-separated integers")
    if any(value < 0 for value in args.starts):
        parser.error("--starts must be nonnegative")
    args.kernels = args.kernels.split(",")
    if not args.kernels or any(name not in ("number-counts", "weak-lensing")
                               for name in args.kernels):
        parser.error("--kernels must select number-counts and/or weak-lensing")
    return args


def main() -> None:
    """Construct, time, and save the requested kernel closures and spectra."""
    args = parse_args()
    # Importing the factory initializes NumCosmo, outside all measured sections.
    from numcosmo_py import Nc, Ncm
    from numcosmo_py.app.xcor.common import XcorKernelCommon

    specs = [f"{name} bin_idx=4 survey=LSST-Y1" for name in args.kernels]
    reference = json.loads(args.reference.read_text()) if args.reference else None
    reference_blocks = {}
    if reference:
        if reference["specs"] != specs:
            raise ValueError("Reference kernel specifications differ from this run")
        for block in reference["blocks"]:
            reference_blocks.setdefault((block["lmin"], block["lmax"]), block)
        for lmin in args.starts:
            if (lmin, lmin + args.size - 1) not in reference_blocks:
                raise ValueError(f"Reference does not contain the block starting at {lmin}")

    Ncm.cfg_init()
    start = time.perf_counter()
    # Reuse the actual CLI factories and suppress their descriptive output.
    with contextlib.redirect_stdout(io.StringIO()):
        common = XcorKernelCommon(kernel=specs)
    for _, kernel in common.kernels:
        kernel.set_l_limber(-1)
    xcor = Nc.Xcor.new(common.dist, common.ps_ml, Nc.XcorMethod.KERNEL_EXACT)
    xcor.set_closure_type(Nc.XcorKernelClosure.CHEBYSHEV)
    xcor.prepare(common.cosmo)
    result = {
        "setup_s": time.perf_counter() - start,
        "library_path": os.environ.get("LD_LIBRARY_PATH"),
        "ld_preload": os.environ.get("LD_PRELOAD"),
        "cpu_affinity": sorted(os.sched_getaffinity(0))
        if hasattr(os, "sched_getaffinity") else None,
        "specs": specs,
        "integrator_reltol": common.integrator.get_reltol(),
        "integrator_cheb_reltol": common.integrator.get_cheb_reltol(),
        "closure_reltol": [kernel.get_reltol() for _, kernel in common.kernels],
        "closure_peak_epsilon": [kernel.get_peak_epsilon() for _, kernel in common.kernels],
        "blocks": [],
    }

    def save() -> None:
        """Preserve partial progress after each completed closure and block."""
        result["peak_rss_kib"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        args.output.write_text(json.dumps(result, indent=2))

    for repeat in range(args.repeat):
        for lmin in args.starts:
            lmax = lmin + args.size - 1
            start = time.perf_counter()
            common.integrator.set_ell_range(lmin, lmax)
            block = {"lmin": lmin, "lmax": lmax, "repeat": repeat,
                     "ell_range_s": time.perf_counter() - start,
                     "kernels": [], "spectra": []}
            result["blocks"].append(block)
            closures = []
            for idx, (label, kernel) in enumerate(common.kernels):
                start = time.perf_counter()
                closure = kernel.get_eval_vectorized(
                    common.cosmo, lmin, lmax, Nc.XcorKernelClosure.CHEBYSHEV
                )
                elapsed = time.perf_counter() - start
                closures.append(closure)
                kmin, kmax = closure.get_range()
                if reference:
                    ref_block = reference_blocks[(lmin, lmax)]
                    ks = np.asarray(ref_block["kernels"][idx]["k"])
                else:
                    ks = np.geomspace(kmin, kmax, 257)
                values = [list(closure.eval_array(float(k))) for k in ks]
                block["kernels"].append({
                    "label": label, "closure_s": elapsed,
                    "range": [kmin, kmax], "panels": closure.get_n_panels(),
                    "component_ranges": [list(closure.get_range_comp(i))
                                         for i in range(closure.get_len())],
                    "k": ks.tolist(), "values": values,
                })
                save()
                print(f"ell={lmin}..{lmax} repeat={repeat} {label}: "
                      f"closure={elapsed:.6f}s panels={closure.get_n_panels()}", flush=True)
            for i, closure_i in enumerate(closures):
                for j in range(i, len(closures)):
                    vec = Ncm.Vector.new(args.size)
                    start = time.perf_counter()
                    xcor.integrate_block(closure_i, closures[j], lmin, lmax, i == j,
                                         Nc.XcorMethod.KERNEL_EXACT, vec, None)
                    elapsed = time.perf_counter() - start
                    block["spectra"].append({"pair": [i, j], "outer_s": elapsed,
                                             "cl": list(vec.dup_array())})
            save()
            print(f"  outer total={sum(s['outer_s'] for s in block['spectra']):.6f}s",
                  flush=True)
    print(f"Saved {args.output}", flush=True)


if __name__ == "__main__":
    main()

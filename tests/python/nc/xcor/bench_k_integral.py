#
# bench_k_integral.py
#
# Fri Aug 28 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# bench_k_integral.py
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

"""Sweep the outer k-integral over the case matrix: accuracy, cost, diagnostics.

Not a test. This is the driver the tolerances in ``test_k_integral.py`` are read
off, and the one a change to the quadrature or the closure is measured on. It
reports what the suite does not assert.

Usage::

    OMP_NUM_THREADS=1 python3 tests/python/nc/xcor/bench_k_integral.py \\
        --out sweep.json
    python3 tests/python/nc/xcor/bench_k_integral.py --baseline sweep.json

Every row carries the settings it was taken at. Threads are the caller's
business, but a method comparison belongs at ``OMP_NUM_THREADS=1``: the ell
block is where the OpenMP team lives, so a thread count silently changes what
is being compared.

Timing is reported in three numbers per block, never collapsed into one:

``build``
    both closures for the block, built through a *fresh* Levin integrator. A
    property of the closure, not of the method -- all three methods pay it.
    Timing a rebuild through an integrator that has already served this kernel
    reads 1.7x cheap, because it carries warm per-ell Bessel caches.
``compute``
    ``nc_xcor_compute()`` end to end, warm: it takes the integrator from the
    kernel, which this harness has already driven.
``quad``
    the outer quadrature alone, through ``Nc.Xcor.integrate_block()`` on the
    closures this harness already holds. Not defined for ``gsl``, which fits a
    closure per multipole and so has no block to be handed; that is what
    ``gsl_block`` is for.
``delta``
    ``compute`` minus the cheapest method's ``compute`` on the same pair,
    closure and block. Everything but the quadrature is common to the methods,
    so this is the part that is the method's, without the subtraction of two
    differently-warmed quantities that a ``compute - build`` would be. Where
    ``quad`` is defined it is the better number; ``delta`` stays the only one
    available for ``gsl``.
"""

from __future__ import annotations

import argparse
import json
import pathlib
import platform
import time
import typing

import numpy as np

from numcosmo_py import Nc, Ncm

try:
    import cases_k_integral as cases
except ImportError:  # running from the repository root
    import sys

    sys.path.insert(0, str(pathlib.Path(__file__).parent))
    import cases_k_integral as cases  # type: ignore[no-redef]


METHODS: typing.Final[dict[str, Nc.XcorMethod]] = {
    "exact": Nc.XcorMethod.KERNEL_EXACT,
    "cubature": Nc.XcorMethod.KERNEL_CUBATURE,
    "gsl": Nc.XcorMethod.KERNEL_GSL,
    "gsl_block": Nc.XcorMethod.KERNEL_GSL_BLOCK,
}

# KERNEL_GSL calls nc_xcor_kernel_get_eval once per multipole; every other
# method calls ..._get_eval_vectorized_full once per block. Each is measured
# against a reference built on the closure it actually integrates, or the row
# reports a closure difference under a method's name. gsl_block runs the same
# qagp as gsl over the block closure, which is what makes it comparable to
# exact and cubature at all.
PER_MULTIPOLE: typing.Final[frozenset[str]] = frozenset({"gsl"})

CLOSURES: typing.Final[dict[str, Nc.XcorKernelClosure]] = {
    "spline": Nc.XcorKernelClosure.SPLINE,
    "chebyshev": Nc.XcorKernelClosure.CHEBYSHEV,
}


def _closure_size(integrand: Nc.XcorKernelIntegrand) -> int:
    """Knots for a spline closure, coefficients summed over panels otherwise."""
    n_panels = integrand.get_n_panels()

    if n_panels == 0:
        knots = integrand.peek_knots()

        return 0 if knots is None else knots.len()

    total = 0

    for i in range(n_panels):
        coeffs, _, _ = integrand.peek_panel(i)
        total += Ncm.Matrix.ncols(coeffs)

    return total


def sweep_case(
    pair: cases.PairSpec,
    settings: cases.Settings,
    method_names: list[str],
    block_starts: list[int],
    repeats: int,
) -> list[dict]:
    """Every method on one pair, against the reference on frozen closures."""
    cosmo, dist, ps = cases.make_cosmo_bits()
    RH = Nc.HICosmo.RH_Mpc(cosmo)
    rows: list[dict] = []

    for lmin in block_starts:
        lmax = lmin + settings.ell_batch_size - 1
        ells = list(range(lmin, lmax + 1))

        kernel_a = cases.build_kernel(pair.kernel_a, cosmo, dist, ps, settings)
        kernel_b = (
            kernel_a
            if pair.isauto
            else cases.build_kernel(pair.kernel_b, cosmo, dist, ps, settings)
        )

        start = time.perf_counter()
        integrand_a = cases.build_integrand(kernel_a, cosmo, lmin, lmax, settings)
        integrand_b = (
            None
            if pair.isauto
            else cases.build_integrand(kernel_b, cosmo, lmin, lmax, settings)
        )
        build_time = time.perf_counter() - start

        reference = cases.reference_cl(RH, integrand_a, integrand_b)
        cancellation = cases.cancellation_ratio(integrand_a, integrand_b)

        references = {"block": reference.cl}

        if PER_MULTIPOLE.intersection(method_names):
            references["per_multipole"] = cases.per_multipole_reference(
                RH,
                kernel_a,
                None if pair.isauto else kernel_b,
                cosmo,
                lmin,
                lmax,
                settings,
            )

        common = {
            "case": pair.case,
            "regime": pair.regime,
            "kernel_a": pair.kernel_a,
            "kernel_b": pair.kernel_b,
            "isauto": pair.isauto,
            "lmin": lmin,
            "lmax": lmax,
            "closure": settings.closure.value_nick,
            "reltol": settings.reltol,
            "scaled_abstol": settings.scaled_abstol,
            "l_limber": settings.l_limber,
            "ell_batch_size": settings.ell_batch_size,
            "size_a": _closure_size(integrand_a),
            "size_b": None if integrand_b is None else _closure_size(integrand_b),
            "n_panels_a": integrand_a.get_n_panels(),
            "build_time": build_time,
            "reference_cells": reference.n_cells,
            "reference_max_order": reference.max_order,
            "reference_worst_move": reference.worst_cell_move,
            "k_min": reference.k_min,
            "k_max": reference.k_max,
        }

        compute_times: dict[str, float] = {}

        for method_name in method_names:
            xcor = Nc.Xcor.new(dist, ps, METHODS[method_name])
            xcor.set_closure_type(settings.closure)
            xcor.set_ell_batch_size(settings.ell_batch_size)
            xcor.set_reltol(settings.reltol)
            xcor.prepare(cosmo)

            vp = Ncm.Vector.new(len(ells))
            timings = []

            for _ in range(repeats):
                start = time.perf_counter()
                xcor.compute(kernel_a, kernel_b, cosmo, lmin, lmax, vp)
                timings.append(time.perf_counter() - start)

            # The quadrature on its own, over the closures already built above.
            # Only the block methods can be driven this way; gsl fits per
            # multipole and has nothing to be handed.
            quad_time = None

            if method_name not in PER_MULTIPOLE:
                vp_quad = Ncm.Vector.new(len(ells))
                quad_timings = []

                for _ in range(repeats):
                    start = time.perf_counter()
                    xcor.integrate_block(
                        integrand_a,
                        integrand_b,
                        lmin,
                        lmax,
                        pair.isauto,
                        METHODS[method_name],
                        vp_quad,
                        None,
                    )
                    quad_timings.append(time.perf_counter() - start)

                quad_time = float(np.median(quad_timings))

                # Same closures, same method: any difference here is the
                # entry point disagreeing with itself, not a tolerance.
                drift = np.abs(
                    np.array(vp_quad.dup_array()) - np.array(vp.dup_array())
                ).max()

                if drift > 0.0:
                    print(
                        f"  ! {pair.case} {method_name}: integrate_block "
                        f"differs from compute by {drift:.3e}",
                        flush=True,
                    )

            got = np.array(vp.dup_array())
            reference_kind = (
                "per_multipole" if method_name in PER_MULTIPOLE else "block"
            )
            truth = references[reference_kind]
            scale = np.abs(truth).max()
            relative = np.where(
                np.abs(truth) > 0.0,
                np.abs(got / np.where(truth == 0.0, 1.0, truth) - 1.0),
                np.abs(got),
            )
            compute_time = float(np.median(timings))
            compute_times[method_name] = compute_time

            for index, ell in enumerate(ells):
                rows.append(
                    {
                        **common,
                        "method": method_name,
                        "reference_kind": reference_kind,
                        "ell": ell,
                        "cl": float(got[index]),
                        "cl_reference": float(truth[index]),
                        "cl_block_reference": float(reference.cl[index]),
                        # Against the block's own peak, not against a value
                        # five orders below it: a relative tolerance on the
                        # smallest entry of a cancelling pair is a much harder
                        # request than the same number on an auto spectrum.
                        "rel_err": float(relative[index]),
                        "peak_err": (
                            float(abs(got[index] - truth[index]) / scale)
                            if scale > 0.0
                            else 0.0
                        ),
                        "cancellation": float(cancellation[index]),
                        "conditioned_err": (
                            float(relative[index] / cancellation[index])
                            if np.isfinite(cancellation[index])
                            else float("nan")
                        ),
                        "compute_time": compute_time,
                        "quad_time": quad_time,
                    }
                )

        cheapest = min(compute_times.values())

        for row in rows[-len(ells) * len(method_names) :]:
            row["compute_delta"] = row["compute_time"] - cheapest

    return rows


def summarize(rows: list[dict]) -> str:
    """Worst relative error and median cost per (case, closure, method)."""
    keys = sorted({(r["case"], r["closure"], r["method"]) for r in rows})
    lines = [
        f"{'case':5s} {'closure':10s} {'method':9s} "
        f"{'worst rel':>10s} {'worst peak':>10s} {'max C':>9s} "
        f"{'build s':>9s} {'quad s':>9s} {'delta s':>9s}"
    ]

    for case, closure, method in keys:
        selected = [
            r
            for r in rows
            if r["case"] == case and r["closure"] == closure and r["method"] == method
        ]
        worst_rel = max(r["rel_err"] for r in selected)
        worst_peak = max(r["peak_err"] for r in selected)
        max_c = max(r["cancellation"] for r in selected)
        build = float(np.median([r["build_time"] for r in selected]))
        delta = float(np.median([r["compute_delta"] for r in selected]))
        quads = [r["quad_time"] for r in selected if r["quad_time"] is not None]
        quad = f"{float(np.median(quads)):9.3f}" if quads else f"{'--':>9s}"
        lines.append(
            f"{case:5s} {closure:10s} {method:9s} "
            f"{worst_rel:10.2e} {worst_peak:10.2e} {max_c:9.2e} "
            f"{build:9.3f} {quad} {delta:9.3f}"
        )

    return "\n".join(lines)


def diff_against(rows: list[dict], baseline: list[dict]) -> str:
    """What moved since a baseline sweep, keyed by everything but the numbers."""

    def key(row: dict) -> tuple:
        return (row["case"], row["closure"], row["method"], row["ell"])

    old = {key(r): r for r in baseline}
    lines = ["moved rows (rel_err ratio, compute time ratio):"]

    for row in rows:
        previous = old.get(key(row))

        if previous is None:
            lines.append(f"  NEW      {key(row)}")
            continue

        err_ratio = (row["rel_err"] + 1.0e-300) / (previous["rel_err"] + 1.0e-300)
        time_ratio = (row["compute_time"] + 1.0e-12) / (
            previous["compute_time"] + 1.0e-12
        )

        if err_ratio > 2.0 or err_ratio < 0.5 or time_ratio > 1.5 or time_ratio < 0.67:
            lines.append(f"  {key(row)}  err x{err_ratio:.2f}  time x{time_ratio:.2f}")

    return "\n".join(lines) if len(lines) > 1 else "nothing moved"


def main() -> None:
    """Run the sweep."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=pathlib.Path, help="write the sweep here")
    parser.add_argument(
        "--baseline", type=pathlib.Path, help="compare against an earlier sweep"
    )
    parser.add_argument(
        "--cases", nargs="*", default=None, help="case ids, default all of them"
    )
    parser.add_argument(
        "--methods", nargs="*", default=list(METHODS), choices=list(METHODS)
    )
    parser.add_argument(
        "--closures", nargs="*", default=list(CLOSURES), choices=list(CLOSURES)
    )
    parser.add_argument("--ells", nargs="*", type=int, default=cases.ELLS)
    parser.add_argument("--reltol", type=float, default=1.0e-4)
    parser.add_argument("--scaled-abstol", type=float, default=1.0e-4)
    parser.add_argument("--ell-batch-size", type=int, default=8)
    parser.add_argument("--repeats", type=int, default=3)
    args = parser.parse_args()

    Ncm.cfg_init()

    selected = (
        cases.PAIRS
        if args.cases is None
        else [cases.PAIRS_BY_CASE[c] for c in args.cases]
    )
    rows: list[dict] = []

    for closure_name in args.closures:
        settings = cases.Settings(
            reltol=args.reltol,
            scaled_abstol=args.scaled_abstol,
            ell_batch_size=args.ell_batch_size,
            closure=CLOSURES[closure_name],
        )

        for pair in selected:
            start = time.perf_counter()
            rows.extend(
                sweep_case(pair, settings, args.methods, args.ells, args.repeats)
            )
            print(
                f"{pair.case} {closure_name}: {time.perf_counter() - start:.1f} s",
                flush=True,
            )

    print()
    print(summarize(rows))

    if args.out is not None:
        payload = {
            "platform": platform.platform(),
            "processor": platform.processor(),
            "rows": rows,
        }
        args.out.write_text(json.dumps(payload, indent=1))
        print(f"\nwrote {len(rows)} rows to {args.out}")

    if args.baseline is not None:
        print()
        print(diff_against(rows, json.loads(args.baseline.read_text())["rows"]))


if __name__ == "__main__":
    main()

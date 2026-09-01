#!/usr/bin/env python3
"""Benchmark NcmDiff first/second derivatives on a battery of test functions.

Measures, per (case, method): actual relative error against the analytic
derivative, the estimated error returned by NcmDiff, whether the estimate
covers the actual error (honesty), and the number of function evaluations.

Usage: python3 dev-notes/diff_benchmark.py [tag]
Writes dev-notes/diff_benchmark_<tag>.json and prints a summary table.
"""

import json
import math
import sys

import numpy as np

from numcosmo_py import Ncm

Ncm.cfg_init()


class Counter:
    """Wraps a scalar function, counting evaluations."""

    def __init__(self, f):
        self.f = f
        self.n = 0

    def __call__(self, x, *args):
        self.n += 1
        return self.f(x)


# Each case: name, f, df, d2f, x0
CASES = [
    ("exp@1", math.exp, math.exp, math.exp, 1.0),
    ("sin@1.3", math.sin, math.cos, lambda x: -math.sin(x), 1.3),
    (
        "sin100@1.3",
        lambda x: math.sin(100.0 * x),
        lambda x: 100.0 * math.cos(100.0 * x),
        lambda x: -1.0e4 * math.sin(100.0 * x),
        1.3,
    ),
    ("cube@2", lambda x: x**3, lambda x: 3.0 * x**2, lambda x: 6.0 * x, 2.0),
    ("log@1e-2", math.log, lambda x: 1.0 / x, lambda x: -1.0 / x**2, 1.0e-2),
    (
        "inv@1e-3",
        lambda x: 1.0 / x,
        lambda x: -1.0 / x**2,
        lambda x: 2.0 / x**3,
        1.0e-3,
    ),
    (
        "sqrt@1e-4",
        math.sqrt,
        lambda x: 0.5 / math.sqrt(x),
        lambda x: -0.25 / x**1.5,
        1.0e-4,
    ),
    (
        "gauss_n@1.02",
        lambda x: math.exp(-0.5 * ((x - 1.0) / 0.01) ** 2),
        lambda x: -((x - 1.0) / 1.0e-4) * math.exp(-0.5 * ((x - 1.0) / 0.01) ** 2),
        lambda x: (((x - 1.0) / 1.0e-4) ** 2 - 1.0e4)
        * math.exp(-0.5 * ((x - 1.0) / 0.01) ** 2),
        1.02,
    ),
    (
        "offset1e10@1.3",
        lambda x: 1.0e10 + math.sin(x),
        math.cos,
        lambda x: -math.sin(x),
        1.3,
    ),
    (
        "tanh50@1.001",
        lambda x: math.tanh(50.0 * (x - 1.0)),
        lambda x: 50.0 / math.cosh(50.0 * (x - 1.0)) ** 2,
        lambda x: -5000.0
        * math.tanh(50.0 * (x - 1.0))
        / math.cosh(50.0 * (x - 1.0)) ** 2,
        1.001,
    ),
    (
        "runge@0.3",
        lambda x: 1.0 / (1.0 + 25.0 * x**2),
        lambda x: -50.0 * x / (1.0 + 25.0 * x**2) ** 2,
        lambda x: (5000.0 * x**2 - 50.0 * (1.0 + 25.0 * x**2))
        / (1.0 + 25.0 * x**2) ** 3,
        0.3,
    ),
    ("exp@100", math.exp, math.exp, math.exp, 100.0),
    (
        "atan1000@1",
        lambda x: math.atan(1000.0 * (x - 1.0)),
        lambda x: 1000.0 / (1.0 + (1000.0 * (x - 1.0)) ** 2),
        lambda x: -2.0e6
        * (1000.0 * (x - 1.0))
        / (1.0 + (1000.0 * (x - 1.0)) ** 2) ** 2,
        1.0,
    ),
    ("sin@0", math.sin, math.cos, lambda x: -math.sin(x), 0.0),
    (
        "x7/2@0.5",
        lambda x: x**3.5,
        lambda x: 3.5 * x**2.5,
        lambda x: 8.75 * x**1.5,
        0.5,
    ),
]


def relerr(val, ref):
    """Relative error of val against ref (absolute if ref is 0)."""
    if ref == 0.0:
        return abs(val)
    return abs(val / ref - 1.0)


def run_method(diff, method, f, x0):
    """Run one NcmDiff method; returns (value, est_err, n_evals)."""
    cnt = Counter(f)
    if method == "rc_d1":
        val, err = diff.rc_d1_1_to_1(x0, cnt, None)
    elif method == "rf_d1":
        val, err = diff.rf_d1_1_to_1(x0, cnt, None)
    elif method == "rc_d2":
        val, err = diff.rc_d2_1_to_1(x0, cnt, None)
    else:
        raise ValueError(method)
    return val, err, cnt.n


def main():
    """Run the battery and print/save results."""
    tag = sys.argv[1] if len(sys.argv) > 1 else "baseline"
    extra_kwargs = {}
    for arg in sys.argv[2:]:
        key, val = arg.split("=")
        extra_kwargs[key] = json.loads(val)
    diff = Ncm.Diff(**extra_kwargs)

    results = []
    print(
        f"{'case':>16s} {'method':>7s} {'actual':>10s} {'est':>10s} "
        f"{'est/act':>9s} {'ok':>3s} {'nev':>5s}"
    )
    for name, f, df, d2f, x0 in CASES:
        for method in ("rc_d1", "rf_d1", "rc_d2"):
            ref = df(x0) if method != "rc_d2" else d2f(x0)
            val, err, nev = run_method(diff, method, f, x0)
            act = relerr(val, ref)
            est = abs(err / ref) if ref != 0.0 else abs(err)
            honest = est >= act
            ratio = est / act if act > 0.0 else float("inf")
            results.append(
                dict(
                    case=name,
                    method=method,
                    actual=act,
                    est=est,
                    honest=honest,
                    nev=nev,
                    val=val,
                    ref=ref,
                )
            )
            print(
                f"{name:>16s} {method:>7s} {act:10.2e} {est:10.2e} "
                f"{ratio:9.1e} {'Y' if honest else 'N':>3s} {nev:5d}"
            )

    # Aggregates per method
    print()
    for method in ("rc_d1", "rf_d1", "rc_d2"):
        rows = [r for r in results if r["method"] == method]
        acts = np.array([max(r["actual"], 1e-18) for r in rows])
        gm = math.exp(np.mean(np.log(acts)))
        wc = acts.max()
        nev = sum(r["nev"] for r in rows)
        dishonest = sum(1 for r in rows if not r["honest"])
        print(
            f"{method}: geo-mean act rel err {gm:.2e}  worst {wc:.2e}  "
            f"total evals {nev}  dishonest {dishonest}/{len(rows)}"
        )

    out = f"dev-notes/diff_benchmark_{tag}.json"
    with open(out, "w", encoding="utf-8") as fj:
        json.dump(results, fj, indent=1)
    print(f"\nsaved {out}")


if __name__ == "__main__":
    main()

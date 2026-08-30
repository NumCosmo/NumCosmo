#!/usr/bin/env python
#
# make_xcor_window_truth_table.py
#
# Wed Aug 27 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# make_xcor_window_truth_table.py
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

"""Generate the committed truth table of certified radial integrals.

Emits ``data/truth_tables/xcor/xcor_window_ilk.json.gz`` holding

.. math:: I_\\ell(k) = \\int W(\\chi)\\, j_\\ell(k\\chi)\\, \\mathrm{d}\\chi

for every analytic window, certified in Arb. This is the offline half of the
check: it needs FLINT, runs in minutes, and is re-run only when a window's
closed form or the case list changes. The test that consumes the table
(``tests/python/nc/xcor/test_xcor_window_truth_table.py``) needs neither FLINT
nor this script.

**Why the k grid is per-multipole.** :math:`I_\\ell(k)` peaks near
:math:`k\\chi \\simeq \\ell`, so a k grid shared across multipoles samples the
peak at one :math:`\\ell` and the dead tail at every other. The grid is laid
down in

.. math:: y = k \\chi_\\mathrm{ref} / (\\ell + 1/2)

instead, with :math:`\\chi_\\mathrm{ref}` the midpoint of the window's support,
which puts the peak at :math:`y \\simeq 1` for every multipole and every shape.

Usage::

    python tests/tools/make_xcor_window_truth_table.py [--ells=2,10,50,200]
"""

import argparse
import concurrent.futures
import gzip
import json
import pathlib
import shutil
import subprocess
import sys
import tempfile
import time

HERE = pathlib.Path(__file__).resolve().parent
SOURCE = HERE / "nc_xcor_kernel_analytic_arb.c"
OUTPUT = (
    HERE.parent.parent / "data" / "truth_tables" / "xcor" / "xcor_window_ilk.json.gz"
)

# The generator's arguments and the library's constructor arguments for each
# shape. `ctor` is written into the table so the test builds its kernel from
# the same numbers that were certified, and `arb` stays separate so a
# convention mismatch between the two sides still has somewhere to show up --
# the same split compare_xcor_window_arb.py makes, for the same reason.
CASES = {
    "gauss": {
        "arb": ["--chi-mean=1500", "--chi-sigma=300", "--n-sigma=4"],
        "ctor": [1500.0, 300.0, 4.0],
    },
    "tophat": {
        "arb": ["--chi-lower=500", "--chi-upper=2500"],
        "ctor": [500.0, 2500.0],
    },
    "tophat_smooth": {
        "arb": [
            "--chi-lower=1000",
            "--chi-upper=2000",
            "--chi-sigma=150",
            "--n-sigma=6",
        ],
        "ctor": [1000.0, 2000.0, 150.0, 6.0],
    },
    "student_t": {
        "arb": ["--chi-mean=1500", "--chi-scale=200", "--nu=2", "--n-scale=6"],
        "ctor": [1500.0, 200.0, 2.0, 6.0],
    },
    "power_exp": {
        "arb": [
            "--chi-scale=1200",
            "--alpha=2",
            "--beta=1.5",
            "--chi-lower=50",
            "--chi-upper=4000",
        ],
        "ctor": [1200.0, 2.0, 1.5, 50.0, 4000.0],
    },
    "lensing": {
        "arb": ["--chi-lower=50", "--chi-source-lower=2000", "--chi-source-upper=3000"],
        "ctor": [50.0, 2000.0, 3000.0],
    },
    "multi": {
        "arb": ["--mu=1000,1600", "--sigma=300,300", "--weight=1,0.6", "--n-sigma=4"],
        "ctor": [[1000.0, 1600.0], [300.0, 300.0], [1.0, 0.6], 4.0],
    },
}

# Below y = 1 the Bessel argument sits under the turning point and I_ell is
# exponentially small -- 1e-74 at ell = 200 -- which is a noise-floor case, not
# an accuracy case. Above it the integrand oscillates and the Levin ODE does
# the actual work, so the grid is denser there.
Y_GRID = [0.2, 0.5, 0.8, 1.0, 1.25, 1.6, 2.5, 4.0, 6.0]

CONVENTION = (
    "table[shape][i_ell][i_k] = int_{support[0]}^{support[1]} W(chi) "
    "j_ell(kvals[shape][i_ell][i_k] * chi) dchi, with W the normalized analytic "
    "window named by shape and constructed from ctor. chi is in Mpc and k in "
    "1/Mpc. Note the measure: dchi, not d(k chi). Values are full-precision "
    "decimal strings, not doubles; radius is the certified absolute half-width "
    "of the Arb ball, so |true - float(value)| <= radius."
)


def build_tool(workdir: pathlib.Path) -> pathlib.Path:
    """Compile the Arb reference generator, or exit with a clear reason."""
    if shutil.which("pkg-config") is None:
        sys.exit("pkg-config not found; cannot locate FLINT")

    if subprocess.run(
        ["pkg-config", "--exists", "flint"], check=False, capture_output=True
    ).returncode:
        sys.exit(
            "FLINT not found by pkg-config. Install it (conda-forge::libflint, "
            "Debian libflint-dev, Homebrew flint) to regenerate this table."
        )

    flags = subprocess.run(
        ["pkg-config", "--cflags", "--libs", "flint"],
        check=True,
        capture_output=True,
        text=True,
    ).stdout.split()
    exe = workdir / "xcor_window_arb"
    subprocess.run(
        ["cc", "-O2", "-o", str(exe), str(SOURCE), *flags, "-lm"], check=True
    )

    return exe


def run_tool(exe: pathlib.Path, shape: str, ell: int, ks, target_rel: float) -> tuple:
    """Return (support, [value strings], [radii]) for one shape and multipole."""
    out = subprocess.run(
        [
            str(exe),
            f"--shape={shape}",
            f"--ell={ell}",
            f"--target-rel={target_rel:g}",
            *CASES[shape]["arb"],
        ],
        input=" ".join(repr(float(k)) for k in ks),
        capture_output=True,
        text=True,
        check=True,
    ).stdout

    support = None
    vals, rads = [], []
    for line in out.splitlines():
        if line.startswith("# shape="):
            lo, hi = line.split("support=[")[1].split("]")[0].split(",")
            support = (float(lo), float(hi))
        elif not line.startswith("#"):
            f = line.split("\t")
            vals.append(f[3])
            rads.append(float(f[4]))

    if support is None or len(vals) != len(ks):
        sys.exit(f"generator produced {len(vals)} of {len(ks)} values for {shape}")

    return support, vals, rads


def one_shape(exe: pathlib.Path, shape: str, ells, target_rel: float) -> dict:
    """Certify every (ell, k) for one shape. Support fixes the k grid's origin."""
    support, _, _ = run_tool(exe, shape, ells[0], [1.0e-3], target_rel)
    chi_ref = 0.5 * (support[0] + support[1])

    kvals, table, radius = [], [], []
    for ell in ells:
        ks = [y * (ell + 0.5) / chi_ref for y in Y_GRID]
        _, vals, rads = run_tool(exe, shape, ell, ks, target_rel)
        kvals.append(ks)
        table.append(vals)
        radius.append(rads)

    return {
        "ctor": CASES[shape]["ctor"],
        "support": list(support),
        "chi_ref": chi_ref,
        "kvals": kvals,
        "table": table,
        "radius": radius,
    }


def main() -> int:
    """Certify every shape and write the compressed table."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ells", default="2,10,50,200")
    parser.add_argument("--target-rel", type=float, default=1.0e-25)
    parser.add_argument("--output", type=pathlib.Path, default=OUTPUT)
    args = parser.parse_args()

    ells = [int(e) for e in args.ells.split(",")]

    with tempfile.TemporaryDirectory() as tmp:
        exe = build_tool(pathlib.Path(tmp))

        # The shapes are independent single-threaded processes and their costs
        # differ by an order of magnitude, so wall time is the slowest shape
        # rather than their sum.
        started = time.monotonic()
        with concurrent.futures.ThreadPoolExecutor(max_workers=len(CASES)) as pool:
            futures = {
                shape: pool.submit(one_shape, exe, shape, ells, args.target_rel)
                for shape in CASES
            }
            shapes = {shape: f.result() for shape, f in futures.items()}
        elapsed = time.monotonic() - started

    payload = {
        "convention": CONVENTION,
        "generator": "tests/tools/make_xcor_window_truth_table.py",
        "target_rel": args.target_rel,
        "ells": ells,
        "y_grid": Y_GRID,
        "shapes": shapes,
    }

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(args.output, "wt") as f:
        json.dump(payload, f, indent=1, sort_keys=True)

    n = len(CASES) * len(ells) * len(Y_GRID)
    print(f"{n} certified values for {len(CASES)} shapes in {elapsed:.1f} s")
    print(f"wrote {args.output} ({args.output.stat().st_size / 1024:.1f} KiB)")

    worst = max(
        r / abs(float(v))
        for s in shapes.values()
        for row_v, row_r in zip(s["table"], s["radius"])
        for v, r in zip(row_v, row_r)
        if float(v) != 0.0
    )
    print(f"worst certified relative radius: {worst:.2e}")

    return 0


if __name__ == "__main__":
    sys.exit(main())

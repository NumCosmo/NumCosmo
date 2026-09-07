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

.. math:: x_\\nu = k \\chi_\\mathrm{ref} / (\\ell + 1/2)

instead, with :math:`\\chi_\\mathrm{ref}` the midpoint of the window's support,
which puts the peak at :math:`x_\\nu \\simeq 1` for every multipole and every shape.
The library writes the Bessel argument as :math:`x = k\\chi`, so this grid is that
argument at mid-support in units of the turning point.

**What the committed table covers.** It is the base check every change to the radial
machinery is measured against, so its grid spans the regimes the solver distinguishes
rather than a convenient corner of them:

===================  ===========================  =========================
:math:`x_\\nu`        regime                       what it exercises
===================  ===========================  =========================
0.2, 0.5             power law, below the         the sub-turning-point branch;
                     turning point                dropped above
                                                  ``POWER_LAW_ELL_MAX``, where
                                                  :math:`j_\\ell` underflows a double
0.8, 1.0, 1.25       the turning point            the runtime guard, the pinned
                                                  constraint
1.6, 2.5, 4, 6       moderately oscillatory       mixed panels
20, 60               deeply oscillatory           the tau constraint, which the
                                                  grid did not reach at all before
===================  ===========================  =========================

Multipoles run to 1000. Near the turning point the precision Arb needs grows steeply
with :math:`\\ell`, roughly a factor of six per decade, so the high multipoles are
certified to a looser target; that is what ``--target-rel`` takes per multipole. Every
target stays far below the :math:`10^{-8}` comparison the test makes, which
``test_table_is_certified_far_below_the_tolerance`` enforces.

Usage, and the command that reproduces the committed table::

    python tests/tools/make_xcor_window_truth_table.py \\
        --ells=2,10,50,200,500,1000 \\
        --target-rel=1e-25,1e-25,1e-25,1e-25,1e-18,1e-18

Two hours wall clock on twelve cores, and ell = 1000 alone is 98% of it: every shape
below ell = 500 finishes inside 25 s, while at ell = 1000 tophat takes 999 s and
multi 6828 s. The multipoles run concurrently, one process per (shape, ell), so the
wall clock is set by the slowest single cell rather than by the total.
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

# Below x_nu = 1 the Bessel argument sits under the turning point and I_ell is
# exponentially small -- 1e-74 at ell = 200 -- which is a noise-floor case, not
# an accuracy case. Above it the integrand oscillates and the Levin ODE does
# the actual work, so the grid is denser there.
X_NU_GRID = [0.2, 0.5, 0.8, 1.0, 1.25, 1.6, 2.5, 4.0, 6.0, 20.0, 60.0]

# Multipole above which the sub-turning-point points are dropped. There j_ell(x) is
# below the smallest double: at ell = 1000 and x_nu = 0.2 it is 1e-328 at mid-support
# and 5e-574 at the far edge, so the library's own answer is exactly 0.0 and comparing
# it against a certified value tests nothing. Below this multipole the same points are
# small but representable (1e-117 at ell = 200), and they do test the power-law branch.
POWER_LAW_ELL_MAX = 200

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
    # FLINT's headers trip -Wstringop-overread and -Warray-bounds under GCC 16, on its
    # own inline functions rather than on anything here; nothing in this compile is ours
    # to fix, so the diagnostics are filtered.
    quiet = ["-Wno-stringop-overread", "-Wno-array-bounds"]
    subprocess.run(
        ["cc", "-O2", "-o", str(exe), str(SOURCE), *flags, *quiet, "-lm"], check=True
    )

    return exe


def run_tool(
    exe: pathlib.Path, shape: str, ell: int, ks, target_rel: float, support_only=False
) -> tuple:
    """Return (support, [value strings], [radii]) for one shape and multipole."""
    out = subprocess.run(
        [
            str(exe),
            f"--shape={shape}",
            f"--ell={ell}",
            f"--target-rel={target_rel:g}",
            *(["--support-only"] if support_only else []),
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


def one_shape(exe: pathlib.Path, shape: str, ells, targets, x_nu_grid=None) -> dict:
    """Certify every (ell, k) for one shape. Support fixes the k grid's origin."""
    support, _, _ = run_tool(exe, shape, ells[0], [], targets[0], support_only=True)
    chi_ref = 0.5 * (support[0] + support[1])
    x_nu_grid = X_NU_GRID if x_nu_grid is None else x_nu_grid

    kvals, table, radius = [], [], []
    for ell, target_rel in zip(ells, targets):
        grid = [xnu for xnu in x_nu_grid if xnu >= 0.8 or ell <= POWER_LAW_ELL_MAX]
        ks = [xnu * (ell + 0.5) / chi_ref for xnu in grid]
        t0 = time.monotonic()
        _, vals, rads = run_tool(exe, shape, ell, ks, target_rel)
        # One line per (shape, ell): the only way to judge how long a run has left,
        # since the cost per multipole grows steeply and unevenly across shapes.
        print(
            f"  {shape:>14} ell={ell:<5} {len(ks)} k in {time.monotonic() - t0:7.1f} s",
            file=sys.stderr,
            flush=True,
        )
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
    parser.add_argument(
        "--x-nu-min",
        type=float,
        default=0.0,
        help=(
            "drop grid points below this x_nu = k chi_ref / (ell + 1/2), the Bessel "
            "argument at mid-support over the turning point. Deep in the "
            "power-law region the integral is hundreds of orders below anything double "
            "precision holds (j_1000 is 1e-328 at x_nu = 0.2), so certifying it to a "
            "relative target costs precision and tests nothing."
        ),
    )
    parser.add_argument(
        "--target-rel",
        default="1e-25",
        help=(
            "certified relative radius to reach, one value or one per multipole. Near "
            "the turning point the precision needed grows steeply with ell, so the high "
            "multipoles are certified less tightly: 1e-18 there is still ten orders "
            "below the comparison the test makes."
        ),
    )
    parser.add_argument("--output", type=pathlib.Path, default=OUTPUT)
    args = parser.parse_args()

    ells = [int(e) for e in args.ells.split(",")]
    targets = [float(t) for t in str(args.target_rel).split(",")]

    if len(targets) == 1:
        targets = targets * len(ells)
    elif len(targets) != len(ells):
        sys.exit("--target-rel takes one value or one per --ells entry")
    x_nu_grid = [xnu for xnu in X_NU_GRID if xnu >= args.x_nu_min]

    with tempfile.TemporaryDirectory() as tmp:
        exe = build_tool(pathlib.Path(tmp))

        # The shapes are independent single-threaded processes and their costs
        # differ by an order of magnitude, so wall time is the slowest shape
        # rather than their sum.
        started = time.monotonic()
        with concurrent.futures.ThreadPoolExecutor(max_workers=len(CASES)) as pool:
            futures = {
                shape: pool.submit(one_shape, exe, shape, ells, targets, x_nu_grid)
                for shape in CASES
            }
            shapes = {shape: f.result() for shape, f in futures.items()}
        elapsed = time.monotonic() - started

    payload = {
        "convention": CONVENTION,
        "generator": "tests/tools/make_xcor_window_truth_table.py",
        "target_rel": targets,
        "ells": ells,
        "x_nu_grid": x_nu_grid,
        "shapes": shapes,
    }

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(args.output, "wt") as f:
        json.dump(payload, f, indent=1, sort_keys=True)

    n = len(CASES) * len(ells) * len(X_NU_GRID)
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

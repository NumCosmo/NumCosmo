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

.. math:: I_\\ell(k) = \\int W(\\chi)\\, j_\\ell^{(d)}(k\\chi)\\, \\mathrm{d}\\chi

for every analytic window, certified in Arb. ``--bessel-deriv`` sets :math:`d`
and, with it, the output file: ``d = 1`` and ``d = 2`` write
``xcor_window_ilk_d1.json.gz`` and ``_d2.json.gz`` rather than overwriting the
``d = 0`` table. The derivative is with respect to the argument
:math:`k\\chi`; :math:`d = 2` is the weight a redshift-space distortion term
carries, so those two tables are what certifies the RSD path. This is the offline half of the
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

The derivative tables, which certify the Bessel-derivative components (RSD)::

    python tests/tools/make_xcor_window_truth_table.py --bessel-deriv=1
    python tests/tools/make_xcor_window_truth_table.py --bessel-deriv=2

Their default multipoles stop at 200 rather than 1000: a derivative order costs three
Bessel evaluations per point instead of one, and the component that carries
:math:`j_\\ell\\'\\'` is a redshift-space term, which is a low-multipole effect. Nothing
about the grid is otherwise different, so the two tables can be read alongside the
:math:`d = 0` one entry by entry.
"""

from __future__ import annotations


import argparse
import concurrent.futures
import gzip
import json
import os
import pathlib
import shutil
import subprocess
import sys
import tempfile
import time

HERE = pathlib.Path(__file__).resolve().parent
SOURCE = HERE / "nc_xcor_kernel_analytic_arb.c"
TABLE_DIR = HERE.parent.parent / "data" / "truth_tables" / "xcor"
OUTPUT = TABLE_DIR / "xcor_window_ilk.json.gz"


def default_output(bessel_deriv: int) -> pathlib.Path:
    """One file per derivative order, so a d > 0 run cannot clobber the d = 0 table.

    The d = 0 table runs to ell = 1000 and costs two hours; the derivative tables
    are a separate, cheaper grid. Keeping them apart means regenerating one never
    puts the other at risk, and each file validates on its own against
    ``test_table_is_certified_far_below_the_tolerance``.
    """
    if bessel_deriv == 0:
        return OUTPUT

    return TABLE_DIR / f"xcor_window_ilk_d{bessel_deriv}.json.gz"

# The generator's arguments and the library's constructor arguments for each
# case. `ctor` is written into the table so the test builds its kernel from
# the same numbers that were certified, and `arb` stays separate so a
# convention mismatch between the two sides still has somewhere to show up --
# the same split compare_xcor_window_arb.py makes, for the same reason.
#
# A case may carry `shape` when its key is not itself a generator shape name:
# several cases are the same closed form at different centres and widths, and
# the key names the case while `shape` names the form. Consumers read the shape
# back from the table rather than parsing the key.
CASES = {
    "gauss_mid": {
        "shape": "gauss",
        "arb": ["--chi-mean=1500", "--chi-sigma=300", "--n-sigma=4"],
        "ctor": [1500.0, 300.0, 4.0],
    },
    "tophat_mid": {
        "shape": "tophat",
        "arb": ["--chi-lower=500", "--chi-upper=2500"],
        "ctor": [500.0, 2500.0],
    },
    # Three further Gaussians, taken from the LSST-SRD bins rather than chosen:
    # the roster's localized shapes all sit at chi = 1500 with sigma/chi ~ 0.2,
    # while the SRD spans centres from 1.1 to 5.4 Gpc and sigma/chi from 0.04 to
    # 0.36. The solver's difficulty is set by the oscillation count across the
    # window and by sigma/chi, so covering one corner of that plane certifies
    # less than the shape count suggests.
    #
    # gauss_far    ~ Y10 lens bin 9  (chi 3790, sigma 163): the thinnest bin in
    #                relative terms, sigma/chi = 0.043, the regime section 10
    #                identifies as where the accuracy degrades.
    # gauss_broad  ~ a wide Gaussian at the scale of Y1 source bin 4 (chi 4520,
    #                sigma 715). It is the *width* that is taken from the survey,
    #                not the form: a source bin is skewed (see power_exp_srd4),
    #                and this case exists to pair with tophat_broad so that width
    #                and edge can be varied independently of skewness.
    # gauss_near   ~ Y10 lens bin 0  (chi 1095, sigma 182): the low-chi end.
    "gauss_far": {
        "shape": "gauss",
        "arb": ["--chi-mean=3800", "--chi-sigma=165", "--n-sigma=4"],
        "ctor": [3800.0, 165.0, 4.0],
    },
    "gauss_broad": {
        "shape": "gauss",
        "arb": ["--chi-mean=4520", "--chi-sigma=715", "--n-sigma=4"],
        "ctor": [4520.0, 715.0, 4.0],
    },
    "gauss_near": {
        "shape": "gauss",
        "arb": ["--chi-mean=1095", "--chi-sigma=182", "--n-sigma=4"],
        "ctor": [1095.0, 182.0, 4.0],
    },
    # Top-hats paired with each Gaussian above, at the same centre and the same
    # variance: a hard-edged window of variance sigma^2 has half-width sqrt(3)
    # sigma. Pairing them isolates the edge from the location and the width --
    # the transform of a hard edge decays as 1/k where the Gaussian's decays
    # exponentially, so a difference between a pair is the edge and nothing else.
    "tophat_far": {
        "arb": ["--chi-lower=3514", "--chi-upper=4086"],
        "ctor": [3514.0, 4086.0],
        "shape": "tophat",
    },
    "tophat_broad": {
        "arb": ["--chi-lower=3282", "--chi-upper=5758"],
        "ctor": [3282.0, 5758.0],
        "shape": "tophat",
    },
    "tophat_near": {
        "arb": ["--chi-lower=780", "--chi-upper=1410"],
        "ctor": [780.0, 1410.0],
        "shape": "tophat",
    },
    "tophat_smooth_mid": {
        "shape": "tophat_smooth",
        "arb": [
            "--chi-lower=1000",
            "--chi-upper=2000",
            "--chi-sigma=150",
            "--n-sigma=6",
        ],
        "ctor": [1000.0, 2000.0, 150.0, 6.0],
    },
    "student_t_mid": {
        "shape": "student_t",
        "arb": ["--chi-mean=1500", "--chi-scale=200", "--nu=2", "--n-scale=6"],
        "ctor": [1500.0, 200.0, 2.0, 6.0],
    },
    "power_exp_mid": {
        "shape": "power_exp",
        "arb": [
            "--chi-scale=1200",
            "--alpha=2",
            "--beta=1.5",
            "--chi-lower=50",
            "--chi-upper=4000",
        ],
        "ctor": [1200.0, 2.0, 1.5, 50.0, 4000.0],
    },
    # The high-redshift windows. CMB lensing and ISW both have support out to the
    # last-scattering surface, chi = 14146 Mpc at z = 1096, and are broad: the
    # lensing efficiency is weighted at chi_mean = 479 with sigma = 1456, the ISW
    # at 2155 with sigma = 1838 and its 5-95% range reaching 5850. No shape above
    # goes past 4000, so neither regime was emulated at all.
    #
    # lensing_cmb      ~ CMB lensing: the same efficiency integral with its sources
    #                    at the last-scattering surface instead of a galaxy bin.
    # power_exp_broad  ~ ISW: broad, skewed, weight from the observer out to z ~ 3.
    "lensing_cmb": {
        "shape": "lensing",
        "arb": [
            "--chi-lower=50",
            "--chi-source-lower=13900",
            "--chi-source-upper=14146",
        ],
        "ctor": [50.0, 13900.0, 14146.0],
    },
    # A skewed window at high chi, which is the regime a source bin occupies and
    # no other case here does: measured on the SRD dn/dz in chi, the Y10 lens bins
    # sit at skew 0.05 while Y1 source bin 4 reaches 0.82. The parameters are round
    # numbers at that location and width (mean 4800, sigma 850, skew 0.6) rather
    # than a fit to the bin. Reproducing a survey's window in detail is not the
    # point and would not survive the next data release; spanning the regimes in
    # forms Arb can certify exactly is.
    "power_exp_skew_far": {
        "shape": "power_exp",
        "arb": [
            "--chi-scale=2400",
            "--alpha=6",
            "--beta=2",
            "--chi-lower=3500",
            "--chi-upper=7200",
        ],
        "ctor": [2400.0, 6.0, 2.0, 3500.0, 7200.0],
    },
    "power_exp_isw": {
        "shape": "power_exp",
        "arb": [
            "--chi-scale=2500",
            "--alpha=1",
            "--beta=1.5",
            "--chi-lower=50",
            "--chi-upper=8000",
        ],
        "ctor": [2500.0, 1.0, 1.5, 50.0, 8000.0],
    },
    "lensing_gal": {
        "shape": "lensing",
        "arb": ["--chi-lower=50", "--chi-source-lower=2000", "--chi-source-upper=3000"],
        "ctor": [50.0, 2000.0, 3000.0],
    },
    "multi_mid": {
        "shape": "multi",
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


# Where meson leaves its test tools, whatever the build directory is called.
BUILD_DIRS = ("build", "builddir", "Optimized", "Debug", "_build")


def find_tool() -> pathlib.Path | None:
    """The meson-built generator, if this tree has been built.

    ``nc_xcor_kernel_analytic_arb`` is a meson target beside the other Arb tools
    (tests/tools/meson.build), so a built tree already has it, compiled by the
    same toolchain and flags as the rest of the project. Using it is not merely
    faster: compiling here again needs a working compiler *inside* whatever
    environment runs the driver, and on a Cray machine the default ``cc`` is a
    wrapper that fails on its own missing pkg-config modules.
    """
    override = os.environ.get("XCOR_ARB_TOOL")

    if override:
        return pathlib.Path(override)

    for build in BUILD_DIRS:
        candidate = HERE.parent.parent / build / "tests" / "tools" / SOURCE.stem

        if candidate.is_file() and os.access(candidate, os.X_OK):
            return candidate

    return None


def build_tool(workdir: pathlib.Path) -> pathlib.Path:
    """The Arb reference generator: the built one, or compiled here as a fallback.

    The fallback exists so the driver still works in a tree that was never built
    -- it is a standalone tool as much as a test one -- but a built tree should
    never reach it.
    """
    found = find_tool()

    if found is not None:
        print(f"using {found}", file=sys.stderr)

        return found

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
    # FLINT's headers trip -Wstringop-overread and -Warray-bounds under GCC 16, on
    # its own inline functions rather than on anything here; nothing in this
    # compile is ours to fix, so the diagnostics are filtered.
    quiet = ["-Wno-stringop-overread", "-Wno-array-bounds"]
    # $CC rather than a hardcoded `cc`, for the Cray case above.
    compiler = os.environ.get("CC", "cc")

    try:
        subprocess.run(
            [compiler, "-O2", "-o", str(exe), str(SOURCE), *flags, *quiet, "-lm"],
            check=True,
        )
    except subprocess.CalledProcessError:
        sys.exit(
            f"{compiler} failed to build the Arb generator, and no built one was "
            f"found in {'/'.join(BUILD_DIRS)}. Build the project, set "
            f"XCOR_ARB_TOOL to the binary, or set CC to a plain compiler."
        )

    return exe


def run_tool(
    exe: pathlib.Path,
    shape: str,
    ell: int,
    ks,
    target_rel: float,
    support_only=False,
    bessel_deriv: int = 0,
) -> tuple:
    """Return (support, [value strings], [radii]) for one shape and multipole."""
    out = subprocess.run(
        [
            str(exe),
            f"--shape={CASES[shape].get('shape', shape)}",
            f"--ell={ell}",
            f"--target-rel={target_rel:g}",
            f"--bessel-deriv={bessel_deriv}",
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


def shape_support(exe: pathlib.Path, shape: str, ell: int, target_rel: float) -> tuple:
    """The window's support and the reference distance its k grid hangs on.

    Both are properties of the window alone, so they are read at
    ``bessel_deriv = 0`` and are the same numbers every derivative table carries --
    which is what lets the tables be compared entry by entry across d.
    """
    support, _, _ = run_tool(exe, shape, ell, [], target_rel, support_only=True)

    return tuple(support), 0.5 * (support[0] + support[1])


def one_cell(
    exe: pathlib.Path,
    shape: str,
    ell: int,
    target_rel: float,
    x_nu_grid,
    bessel_deriv: int = 0,
) -> dict:
    """Certify one (shape, ell): the unit of work, and of loss.

    A cell at ell = 1000 can run for an hour, so it is also the unit the run
    records. Anything coarser means an interruption discards finished
    arbitrary-precision work -- which is not hypothetical: a reboot once cost six
    completed cells that were held in memory awaiting a whole pool.
    """
    support, chi_ref = shape_support(exe, shape, ell, target_rel)
    grid = [xnu for xnu in x_nu_grid if xnu >= 0.8 or ell <= POWER_LAW_ELL_MAX]
    ks = [xnu * (ell + 0.5) / chi_ref for xnu in grid]
    t0 = time.monotonic()
    _, vals, rads = run_tool(
        exe, shape, ell, ks, target_rel, bessel_deriv=bessel_deriv
    )
    elapsed = time.monotonic() - t0
    # One line per cell: the only way to judge how long a run has left, since the
    # cost per multipole grows steeply and unevenly across shapes.
    print(
        f"  {shape:>18} ell={ell:<5} {len(ks)} k in {elapsed:7.1f} s",
        file=sys.stderr,
        flush=True,
    )

    return {
        "shape_form": CASES[shape].get("shape", shape),
        "ctor": CASES[shape]["ctor"],
        "support": list(support),
        "chi_ref": chi_ref,
        "kvals": ks,
        "table": vals,
        "radius": rads,
        "seconds": elapsed,
    }


def partial_path(output: pathlib.Path) -> pathlib.Path:
    """Where finished cells are recorded before they can form a whole table.

    The committed table is rectangular -- every shape carries every multipole in
    the header -- so a shape cannot enter it until all of its cells exist. This
    file holds them meanwhile, and is what makes a run resumable.
    """
    return output.with_suffix(".partial.json")


def load_partial(output: pathlib.Path, resume: bool) -> dict:
    """Cells already certified, keyed ``shape|ell``.

    Seeded from the committed table too, so re-running with an extra multipole
    recomputes only the new cells rather than the whole grid.
    """
    cells: dict = {}

    if not resume:
        return cells

    path = partial_path(output)

    if path.exists():
        cells.update(json.loads(path.read_text())["cells"])

    if output.exists():
        with gzip.open(output, "rt") as handle:
            table = json.load(handle)

        for shape, entry in table["shapes"].items():
            for index, ell in enumerate(table["ells"]):
                key = f"{shape}|{ell}"
                cells.setdefault(
                    key,
                    {
                        "shape_form": entry.get("shape", shape),
                        "ctor": entry["ctor"],
                        "support": entry["support"],
                        "chi_ref": entry["chi_ref"],
                        "kvals": entry["kvals"][index],
                        "table": entry["table"][index],
                        "radius": entry["radius"][index],
                        "seconds": 0.0,
                    },
                )

    return cells


def write_partial(output: pathlib.Path, cells: dict) -> None:
    """Record the finished cells, atomically enough to survive a hard stop."""
    path = partial_path(output)
    tmp = path.with_suffix(".tmp")
    tmp.write_text(json.dumps({"cells": cells}, indent=1))
    tmp.replace(path)


def assemble(cells: dict, ells, targets, x_nu_grid, bessel_deriv: int) -> dict:
    """The rectangular payload, carrying every shape whose cells are all present."""
    shapes = {}

    for shape in CASES:
        keys = [f"{shape}|{ell}" for ell in ells]

        if not all(key in cells for key in keys):
            continue

        first = cells[keys[0]]
        shapes[shape] = {
            "shape": first["shape_form"],
            "ctor": first["ctor"],
            "support": first["support"],
            "chi_ref": first["chi_ref"],
            "kvals": [cells[key]["kvals"] for key in keys],
            "table": [cells[key]["table"] for key in keys],
            "radius": [cells[key]["radius"] for key in keys],
        }

    return {
        "convention": CONVENTION,
        "generator": "tests/tools/make_xcor_window_truth_table.py",
        "target_rel": targets,
        "ells": list(ells),
        "x_nu_grid": list(x_nu_grid),
        "bessel_deriv": bessel_deriv,
        "shapes": shapes,
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
    parser.add_argument(
        "--bessel-deriv",
        type=int,
        default=0,
        choices=(0, 1, 2),
        help=(
            "derivative order of the spherical Bessel weight, taken with respect to "
            "the argument k chi. 2 is what a redshift-space distortion term carries. "
            "Changes the default output file so a derivative run cannot overwrite the "
            "d = 0 table."
        ),
    )
    parser.add_argument(
        "--shapes",
        nargs="*",
        default=None,
        help=(
            "certify only these shapes. Adding a case to CASES otherwise means "
            "recomputing every existing one, including the ell = 1000 cells that "
            "cost hours; with this the new shapes go to their own file and are "
            "merged into the committed table per shape."
        ),
    )
    parser.add_argument(
        "--jobs",
        type=int,
        default=len(CASES),
        help=(
            "concurrent shapes, one single-threaded process each. Defaults to one per "
            "shape; lower it when sharing the machine with another generation lane."
        ),
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        help=(
            "skip cells already certified, read from the output table and from the "
            "partial store beside it. A cell is one (shape, multipole) and can run "
            "for an hour, so this is the difference between losing a run and losing "
            "a cell."
        ),
    )
    parser.add_argument("--output", type=pathlib.Path, default=None)
    args = parser.parse_args()

    output = default_output(args.bessel_deriv) if args.output is None else args.output

    ells = [int(e) for e in args.ells.split(",")]
    targets = [float(t) for t in str(args.target_rel).split(",")]

    if len(targets) == 1:
        targets = targets * len(ells)
    elif len(targets) != len(ells):
        sys.exit("--target-rel takes one value or one per --ells entry")
    x_nu_grid = [xnu for xnu in X_NU_GRID if xnu >= args.x_nu_min]

    selected = list(CASES) if args.shapes is None else args.shapes

    for shape in selected:
        if shape not in CASES:
            sys.exit(f"unknown shape {shape!r}; known: {', '.join(CASES)}")

    output.parent.mkdir(parents=True, exist_ok=True)
    cells = load_partial(output, args.resume)
    tasks = [
        (shape, ell, target)
        for shape in selected
        for ell, target in zip(ells, targets)
        if f"{shape}|{ell}" not in cells
    ]

    if cells:
        print(f"resuming: {len(cells)} cells already certified", flush=True)

    print(f"{len(tasks)} cells to certify", flush=True)

    with tempfile.TemporaryDirectory() as tmp:
        exe = build_tool(pathlib.Path(tmp))

        # Cells are independent single-threaded processes whose costs differ by
        # orders of magnitude, so the wall time is the slowest cell rather than
        # their sum -- and each is recorded as it lands, so an interruption costs
        # only what was in flight.
        started = time.monotonic()
        with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
            futures = {
                pool.submit(
                    one_cell, exe, shape, ell, target, x_nu_grid, args.bessel_deriv
                ): (shape, ell)
                for shape, ell, target in tasks
            }

            for future in concurrent.futures.as_completed(futures):
                shape, ell = futures[future]
                cells[f"{shape}|{ell}"] = future.result()
                write_partial(output, cells)

                payload = assemble(
                    cells, ells, targets, x_nu_grid, args.bessel_deriv
                )

                # The table is rewritten whenever a shape becomes complete, so a
                # run that is stopped early still leaves a usable, rectangular
                # table holding every shape that finished.
                if payload["shapes"]:
                    with gzip.open(output, "wt") as handle:
                        json.dump(payload, handle, indent=1, sort_keys=True)

        elapsed = time.monotonic() - started

    payload = assemble(cells, ells, targets, x_nu_grid, args.bessel_deriv)
    done = payload["shapes"]
    print(
        f"{len(cells)} cells certified at bessel-deriv {args.bessel_deriv} "
        f"({len(tasks)} this run, {elapsed:.1f} s)"
    )
    print(f"{len(done)} of {len(selected)} shapes complete over ells {ells}")

    if done:
        print(f"wrote {output} ({output.stat().st_size / 1024:.1f} KiB)")

    missing = [s for s in selected if s not in done]

    if missing:
        print(f"incomplete, still in {partial_path(output).name}: {missing}")

    # `done`, not a name from before the resumable rewrite: the summary reports on
    # what was written, and a run that completes no shape writes nothing to report.
    radii = [
        r / abs(float(v))
        for entry in done.values()
        for row_v, row_r in zip(entry["table"], entry["radius"])
        for v, r in zip(row_v, row_r)
        if float(v) != 0.0
    ]

    if radii:
        print(f"worst certified relative radius: {max(radii):.2e}")

    return 0


if __name__ == "__main__":
    sys.exit(main())

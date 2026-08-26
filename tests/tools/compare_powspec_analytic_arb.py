#!/usr/bin/env python
#
# compare_powspec_analytic_arb.py
#
# Tue Aug 26 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# compare_powspec_analytic_arb.py
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

"""Compare NcmPowspecAnalytic against certified Arb reference values.

Builds and runs ``ncm_powspec_analytic_arb.c``, which re-implements the same
closed form in Arb ball arithmetic, and reports the deviation of the library's
double-precision evaluation from the certified midpoint.

The reference carries a proved radius, so a disagreement is attributable: it is
the library's floating-point evaluation, never the reference's accuracy.

Usage::

    python tests/tools/compare_powspec_analytic_arb.py [--n-k=61] [--verbose]
"""

import argparse
import pathlib
import shutil
import subprocess
import sys
import tempfile
from decimal import Decimal

import numpy as np

from numcosmo_py import Ncm

HERE = pathlib.Path(__file__).resolve().parent
SOURCE = HERE / "ncm_powspec_analytic_arb.c"

SHAPES = {
    0: Ncm.PowspecAnalyticShape.POWER_LAW,
    1: Ncm.PowspecAnalyticShape.BBKS,
    2: Ncm.PowspecAnalyticShape.RATIONAL,
}
GROWTHS = {
    0: Ncm.PowspecAnalyticGrowth.NONE,
    1: Ncm.PowspecAnalyticGrowth.LCDM,
    2: Ncm.PowspecAnalyticGrowth.RATIONAL,
}
SHAPE_ARG = {0: "power_law", 1: "bbks", 2: "rational"}
GROWTH_ARG = {0: "none", 1: "lcdm", 2: "rational"}


def build_tool(workdir: pathlib.Path) -> pathlib.Path:
    """Compile the Arb reference generator, or exit with a clear reason."""
    if shutil.which("pkg-config") is None:
        sys.exit("pkg-config not found; cannot locate FLINT")

    probe = subprocess.run(
        ["pkg-config", "--exists", "flint"], check=False, capture_output=True
    )
    if probe.returncode != 0:
        sys.exit(
            "FLINT not found by pkg-config. Install it (conda-forge::libflint, "
            "Debian libflint-dev, Homebrew flint) to regenerate the reference."
        )

    flags = subprocess.run(
        ["pkg-config", "--cflags", "--libs", "flint"],
        check=True,
        capture_output=True,
        text=True,
    ).stdout.split()
    exe = workdir / "powspec_arb"
    subprocess.run(
        ["cc", "-O2", "-o", str(exe), str(SOURCE), *flags, "-lm"], check=True
    )
    return exe


def run_tool(exe: pathlib.Path, shape: int, growth: int, n_k: int, bao: float):
    """Run the generator for one shape/growth pair and yield its rows."""
    out = subprocess.run(
        [
            str(exe),
            f"--shape={SHAPE_ARG[shape]}",
            f"--growth={GROWTH_ARG[growth]}",
            f"--n-k={n_k}",
            f"--bao-amplitude={bao}",
        ],
        check=True,
        capture_output=True,
        text=True,
    ).stdout

    for line in out.splitlines():
        if line.startswith("#"):
            continue
        f = line.split("\t")
        # The value is printed with 30 significant digits, well past double, so
        # it is parsed as Decimal and only narrowed at the comparison.
        yield float(f[2]), float(f[3]), Decimal(f[4]), float(f[5]), int(f[6])


def main() -> int:
    """Compare every shape/growth combination and report the worst deviation."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n-k", type=int, default=61)
    parser.add_argument("--bao-amplitude", type=float, default=0.0)
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    Ncm.cfg_init()

    with tempfile.TemporaryDirectory() as tmp:
        exe = build_tool(pathlib.Path(tmp))

        print(
            f"{'shape':<10} {'growth':<9} {'n':>5} {'max |rel|':>11} "
            f"{'rms |rel|':>11} {'max radius':>11} {'prec':>5}"
        )
        worst = 0.0

        for shape in SHAPES:
            for growth in GROWTHS:
                psa = Ncm.PowspecAnalytic.new(SHAPES[shape], GROWTHS[growth])
                if args.bao_amplitude != 0.0:
                    psa.set_property("bao-amplitude", args.bao_amplitude)
                psa.prepare(None)

                rel, radii, precs = [], [], []
                for k, z, ref, radius, prec in run_tool(
                    exe, shape, growth, args.n_k, args.bao_amplitude
                ):
                    got = Decimal(repr(psa.eval(None, z, k)))
                    rel.append(abs(float((got - ref) / ref)))
                    radii.append(radius / abs(float(ref)) if ref != 0 else 0.0)
                    precs.append(prec)
                    if args.verbose:
                        print(f"    k={k:11.4e} z={z:6.2f} rel={rel[-1]:.3e}")

                rel = np.array(rel)
                worst = max(worst, rel.max())
                print(
                    f"{SHAPE_ARG[shape]:<10} {GROWTH_ARG[growth]:<9} {len(rel):>5} "
                    f"{rel.max():>11.3e} {np.sqrt((rel**2).mean()):>11.3e} "
                    f"{max(radii):>11.3e} {max(precs):>5}"
                )

    print(f"\nworst relative deviation over all combinations: {worst:.3e}")
    # Double precision is ~2.2e-16; a few ulps accumulated over the product of
    # transfer, growth and amplitude is the expected floor.
    if worst > 1.0e-13:
        print("FAIL: exceeds the double-precision floor", file=sys.stderr)
        return 1

    print("PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())

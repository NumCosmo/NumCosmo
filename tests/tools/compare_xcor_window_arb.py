#!/usr/bin/env python
#
# compare_xcor_window_arb.py
#
# Tue Aug 26 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# compare_xcor_window_arb.py
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

"""Compare the analytic xcor windows against certified Arb reference values.

Two things are checked, separately, because they fail for different reasons:

**The window itself.** ``nc_xcor_kernel_analytic_eval_W`` is compared pointwise
against the closed form re-implemented in Arb, including the normalization --
which the reference *measures* as a certified integral rather than copying the
library's formula.

**The radial integral.** ``I_ell(k) = int W(chi) j_ell(k chi) dchi`` computed by
``ncm_sbessel_integrator_integrate_ell`` on the library's own window, against
the same integral certified in Arb.

Keeping them apart matters: a disagreement in the first is a wrong window, a
disagreement in the second with the first clean is the integrator. The known
panel-alignment limit of ``NcmSBesselIntegratorLevin`` shows up only in the
second, and only where ``I_ell`` is far below its own peak.

Usage::

    python tests/tools/compare_xcor_window_arb.py [--shape=gauss] [--verbose]
"""

import argparse
import pathlib
import shutil
import subprocess
import sys
import tempfile
from decimal import Decimal, getcontext

import numpy as np

from numcosmo_py import Nc, Ncm

getcontext().prec = 45

HERE = pathlib.Path(__file__).resolve().parent
SOURCE = HERE / "nc_xcor_kernel_analytic_arb.c"

# One representative parameter set per shape: the arguments the Arb generator
# takes, and the constructor arguments the library takes. They are written out
# separately on purpose -- a shared helper that derived one from the other
# would hide exactly the convention mismatch this script exists to catch.
CASES = {
    "gauss": {
        "arb": ["--chi-mean=1500", "--chi-sigma=300", "--n-sigma=4"],
        "make": lambda d, ps, sbi: Nc.XcorKernelAnalyticGauss.new_full(
            d, ps, 1500.0, 300.0, 4.0, sbi
        ),
    },
    "tophat": {
        "arb": ["--chi-lower=500", "--chi-upper=2500"],
        "make": lambda d, ps, sbi: Nc.XcorKernelAnalyticTophat.new_full(
            d, ps, 500.0, 2500.0, sbi
        ),
    },
    "tophat_smooth": {
        "arb": [
            "--chi-lower=1000",
            "--chi-upper=2000",
            "--chi-sigma=150",
            "--n-sigma=6",
        ],
        "make": lambda d, ps, sbi: Nc.XcorKernelAnalyticTophatSmooth.new_full(
            d, ps, 1000.0, 2000.0, 150.0, 6.0, sbi
        ),
    },
    "student_t": {
        "arb": ["--chi-mean=1500", "--chi-scale=200", "--nu=2", "--n-scale=6"],
        "make": lambda d, ps, sbi: Nc.XcorKernelAnalyticStudentT.new_full(
            d, ps, 1500.0, 200.0, 2.0, 6.0, sbi
        ),
    },
    "power_exp": {
        "arb": [
            "--chi-scale=1200",
            "--alpha=2",
            "--beta=1.5",
            "--chi-lower=50",
            "--chi-upper=4000",
        ],
        "make": lambda d, ps, sbi: Nc.XcorKernelAnalyticPowerExp.new_full(
            d, ps, 1200.0, 2.0, 1.5, 50.0, 4000.0, sbi
        ),
    },
    "lensing": {
        "arb": [
            "--chi-lower=50",
            "--chi-source-lower=2000",
            "--chi-source-upper=3000",
        ],
        "make": lambda d, ps, sbi: Nc.XcorKernelAnalyticLensing.new_full(
            d, ps, 50.0, 2000.0, 3000.0, sbi
        ),
    },
    "multi": {
        "arb": [
            "--mu=1000,1600",
            "--sigma=300,300",
            "--weight=1,0.6",
            "--n-sigma=4",
        ],
        "make": lambda d, ps, sbi: Nc.XcorKernelAnalyticMulti.new_full(
            d,
            ps,
            Ncm.Vector.new_array([1000.0, 1600.0]),
            Ncm.Vector.new_array([300.0, 300.0]),
            Ncm.Vector.new_array([1.0, 0.6]),
            4.0,
            sbi,
        ),
    },
}


def build_tool(workdir: pathlib.Path) -> pathlib.Path:
    """Compile the Arb reference generator, or exit with a clear reason."""
    if shutil.which("pkg-config") is None:
        sys.exit("pkg-config not found; cannot locate FLINT")

    if subprocess.run(
        ["pkg-config", "--exists", "flint"], check=False, capture_output=True
    ).returncode:
        sys.exit(
            "FLINT not found by pkg-config. Install it (conda-forge::libflint, "
            "Debian libflint-dev, Homebrew flint) to run this comparison."
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


def run_tool(exe: pathlib.Path, shape: str, ell: int, ks) -> tuple:
    """Return (support, norm, {k: (value, radius)}) for one shape and ell."""
    out = subprocess.run(
        [str(exe), f"--shape={shape}", f"--ell={ell}", *CASES[shape]["arb"]],
        input=" ".join(repr(float(k)) for k in ks),
        capture_output=True,
        text=True,
        check=True,
    ).stdout

    support = norm = None
    vals = {}
    for line in out.splitlines():
        if line.startswith("# shape="):
            lo, hi = line.split("support=[")[1].split("]")[0].split(",")
            support = (float(lo), float(hi))
            norm = Decimal(line.split("norm=")[1].strip())
        elif not line.startswith("#"):
            f = line.split("\t")
            vals[float(f[2])] = (Decimal(f[3]), float(f[4]))

    return support, norm, vals


def run_window(exe: pathlib.Path, shape: str, chis) -> dict:
    """Return {chi: certified W(chi)} for one shape."""
    out = subprocess.run(
        [str(exe), f"--shape={shape}", "--window", *CASES[shape]["arb"]],
        input=" ".join(repr(float(c)) for c in chis),
        capture_output=True,
        text=True,
        check=True,
    ).stdout

    return {
        float(line.split("\t")[1]): Decimal(line.split("\t")[2])
        for line in out.splitlines()
        if not line.startswith("#")
    }


def main() -> int:
    """Compare every shape's window and radial integral against Arb."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--shape", default=None, choices=list(CASES))
    parser.add_argument("--ells", default="2,10,50")
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    Ncm.cfg_init()
    from numcosmo_py.cosmology import Cosmology  # noqa: PLC0415

    cosmology = Cosmology.default()
    cosmo, dist = cosmology.cosmo, cosmology.dist
    dist.prepare(cosmo)
    ps = Ncm.PowspecAnalytic.new(
        Ncm.PowspecAnalyticShape.BBKS, Ncm.PowspecAnalyticGrowth.LCDM
    )

    shapes = [args.shape] if args.shape else list(CASES)
    ells = [int(e) for e in args.ells.split(",")]
    ks = [1.0e-4, 1.0e-3, 3.0e-3, 1.0e-2, 3.0e-2, 0.1]

    worst_w = worst_i = worst_r = 0.0
    with tempfile.TemporaryDirectory() as tmp:
        exe = build_tool(pathlib.Path(tmp))

        print(
            f"{'shape':<14} {'dW/W_peak':>11} "
            f"{'dI/I_peak':>11} {'dI/I (>1e-8)':>11} {'max radius':>11}"
        )

        for shape in shapes:
            sbi = Ncm.SBesselIntegratorLevin.new(0, 8)
            kern = CASES[shape]["make"](dist, ps, sbi)
            kern.set_l_limber(-1)
            kern.prepare(cosmo)

            support, _, _ = run_tool(exe, shape, ells[0], ks[:1])
            lo, hi = support

            # 1. the window itself, pointwise, against the certified closed
            #    form -- normalization included, since the reference measures
            #    it rather than copying the library's formula.
            chis = np.linspace(lo + 1e-9 * (hi - lo), hi - 1e-9 * (hi - lo), 61)
            wref = run_window(exe, shape, chis)
            wpeak = max(abs(float(v)) for v in wref.values())
            # Measured against the window's own peak, not pointwise-relative.
            # A truncated shape is tiny near its edge and some closed forms
            # cancel there -- `lensing` loses every significant digit within an
            # ulp of chi = chi_source_upper, where the true value goes as
            # (b - chi)^2 / 2chi. Pointwise-relative that reads as a factor of
            # 90; against the peak it is 1e-15, which is what any integral of
            # the window actually sees.
            dw = (
                max(
                    abs(
                        float(Decimal(repr(kern.eval_W(float(chi)))) - wref[float(chi)])
                    )
                    for chi in chis
                )
                / wpeak
            )

            # 2. the radial integral, against certified Arb
            # The window returns a hard zero outside its support, so a node
            # landing a hair outside -- the y = k chi round trip makes that
            # unavoidable at the endpoints -- is a cliff the Chebyshev fit
            # cannot resolve, and the fit aborts on max-order. The library's own
            # integration path clamps for exactly this reason
            # (_nc_xcor_kernel_analytic_eval_W_comp_clamped), so the comparison
            # clamps too.
            def window_cb(_ud, chi, _k, _lo=lo, _hi=hi, _kern=kern):
                return _kern.eval_W(min(max(chi, _lo), _hi))

            di_worst = 0.0  # peak-relative
            di_rel = 0.0  # plain relative, where |I| is not negligible
            rad_max = 0.0
            for ell in ells:
                sbi = Ncm.SBesselIntegratorLevin.new(ell, ell)
                _, _, vals = run_tool(exe, shape, ell, ks)
                peak = max(abs(float(v)) for v, _ in vals.values())
                for k in ks:
                    ref, rad = vals[float(k)]
                    if ref == 0:
                        continue
                    got = Decimal(
                        repr(sbi.integrate_ell(window_cb, lo, hi, float(k), ell, None))
                    )
                    rel = abs(float(got / ref - 1))
                    rad_max = max(rad_max, rad / abs(float(ref)))
                    di_worst = max(di_worst, abs(float(got - ref)) / peak)

                    # A relative error on a value 100 orders below the peak
                    # says nothing -- both sides are zero for every purpose.
                    if abs(float(ref)) / peak > 1.0e-8:
                        di_rel = max(di_rel, rel)

                    if args.verbose:
                        print(f"    ell={ell:<4} k={k:<9.4g} rel={rel:.3e}")

            worst_w = max(worst_w, dw)
            worst_i = max(worst_i, di_worst)
            worst_r = max(worst_r, di_rel)
            print(
                f"{shape:<14} {dw:>11.2e} {di_worst:>11.2e} "
                f"{di_rel:>11.2e} {rad_max:>11.2e}"
            )

    print(f"\nworst window deviation, vs its peak      : {worst_w:.3e}")
    print(f"worst radial deviation, vs its peak      : {worst_i:.3e}")
    print(f"worst radial relative error, |I|/peak>1e-8: {worst_r:.3e}")

    return 0


if __name__ == "__main__":
    sys.exit(main())

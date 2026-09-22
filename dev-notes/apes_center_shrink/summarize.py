"""Summarise APES catalogs with the NumCosmo CLI diagnostics.

Two passes, as the convergence workflow requires: first with no burn-in so the
Constant Break statistic can suggest where burn-in ends, then again with that cut so
that tau is measured on the converged part only. tau is flagged as unreliable until
there are at least MIN_ITER iterations past the cut.
"""

import glob
import os
import re
import shutil
import subprocess
import sys

CLI = [shutil.which("numcosmo") or "numcosmo", "catalog", "analyze"]
MIN_ITER = 130


ENV = {**os.environ, "COLUMNS": "220"}  # keep rich from truncating the value column


def run(f, burnin):
    out = subprocess.run(
        CLI + [f, "--burnin", str(burnin)], capture_output=True, text=True, env=ENV
    ).stdout
    return "\n".join(l for l in out.splitlines() if "WARNING" not in l)


def columns(line):
    return [
        re.sub(r"\s+", " ", c).strip()
        for c in line.split("│")
        if re.sub(r"\s+", " ", c).strip()
    ]


def grab(out, label, idx=-1):
    for line in out.splitlines():
        if label in line:
            c = columns(line)
            if len(c) >= 2:
                return c[idx]
    return "?"


pattern = sys.argv[1] if len(sys.argv) > 1 else "cosmo_*.fits"
print(
    f"{'catalog':52s} {'chains':>6s} {'iters':>6s} {'CB':>5s} {'post':>5s} "
    f"{'tau':>9s} {'R-1':>8s} {'ESS':>7s} {'HW':>6s}"
)
for f in sorted(glob.glob(pattern)):
    # pass 1: no burn-in, ask the Constant Break where burn-in ends
    o0 = run(f, 0)
    nit = grab(o0, "Number of Iterations")
    cb = grab(o0, "Constant", 1)  # the suggested cut-off column
    try:
        cb_i = int(cb)
    except ValueError:
        cb_i = 0
    # pass 2: measure on the converged part only
    o1 = run(f, cb_i)
    nch = grab(o1, "Number of chains")
    post = grab(o1, "Number of Iterations")
    tau = grab(o1, "Autocorrelat")
    gr = grab(o1, "Gelman-Rubin")
    hw = grab(o1, "Heidelberger")
    ess = grab(o1, "Effective")
    try:
        flag = "" if int(post) >= MIN_ITER else " (short)"
    except ValueError:
        flag = ""
    print(
        f"{f[:52]:52s} {nch:>6s} {nit:>6s} {cb:>5s} {post:>5s} {tau:>9s} {gr:>8s} {ess:>7s} {hw:>6s}{flag}"
    )

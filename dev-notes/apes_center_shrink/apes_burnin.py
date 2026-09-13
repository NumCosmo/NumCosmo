"""Burn-in behaviour: shrinkage forces the proposal covariance to match the CURRENT
ensemble covariance, so a collapsed starting ensemble could expand more slowly than
with the (over-dispersed) plain KDE. Start the walkers in a tight ball and track how
the ensemble covariance grows toward the target."""
import sys
import numpy as np

from numcosmo_py import Ncm
from numcosmo_py.sampling.apes import APES
from numcosmo_py.interpolation.stats_dist import InterpolationKernel

Ncm.cfg_init()

d = int(sys.argv[1]); nwalkers = int(sys.argv[2]); niter = int(sys.argv[3])
kname = sys.argv[4]; h = float(sys.argv[5]); shrink = sys.argv[6] == "1"
lf = float(sys.argv[7]); interp = sys.argv[8] != "0"
start_scale = float(sys.argv[9]) if len(sys.argv) > 9 else 0.1
KMAP = {"cauchy": InterpolationKernel.CAUCHY, "st3": InterpolationKernel.ST3, "gauss": InterpolationKernel.GAUSS}

rng = np.random.default_rng(7)
ev = np.exp(rng.uniform(np.log(1e-2), 0.0, d))
Q, _ = np.linalg.qr(rng.normal(size=(d, d)))
S = Q @ np.diag(ev) @ Q.T
L = np.linalg.cholesky(S)
Sinv = np.linalg.inv(S)

apes = APES(nwalkers=nwalkers, ndim=d, model=None,
            log_prob=lambda t, _: -0.5 * np.asarray(t) @ Sinv @ np.asarray(t),
            interpolation_kernel=KMAP[kname], over_smooth=h, local_fraction=lf,
            center_shrink=shrink, use_interpolation=interp)

# tight, offset starting ball
x0 = start_scale * (rng.normal(size=(nwalkers, d)) @ L.T) + 2.0 * np.sqrt(np.diag(S))
mcat = apes.esmcmc.peek_catalog()
step = max(niter // 20, 1)
trace = []
for it in range(0, niter, step):
    apes.run_mcmc(x0, step)
    rows = np.array([mcat.peek_row(mcat.len() - nwalkers + i).dup_array() for i in range(nwalkers)])
    ens = rows[:, 1:]
    C = np.cov(ens.T)
    trace.append((it + step, np.trace(C @ Sinv) / d, np.mean(rows[:, 0])))
tag = f"{kname} h={h} shrink={int(shrink)} interp={int(interp)}"
print(f"# {tag} d={d} nwalkers={nwalkers} start_scale={start_scale}", flush=True)
for it, tr, m2 in trace:
    print(f"{tag},{it},{tr:.4f},{m2:.2f}", flush=True)
# iterations to reach 90 % of the target covariance trace
reached = [it for it, tr, _ in trace if tr > 0.9]
print(f"# {tag}: iters_to_90pct={reached[0] if reached else 'NOT REACHED'} final_trace={trace[-1][1]:.3f} final_mean_m2lnL={trace[-1][2]:.1f} (target {d})", flush=True)

"""End-to-end APES run on a correlated Gaussian: acceptance and tau, shrink on/off."""
import sys
import time

import numpy as np

from numcosmo_py import Ncm
from numcosmo_py.sampling.apes import APES
from numcosmo_py.interpolation.stats_dist import InterpolationKernel

Ncm.cfg_init()

d = int(sys.argv[1])
nwalkers = int(sys.argv[2])
niter = int(sys.argv[3])
kernel = {"cauchy": InterpolationKernel.CAUCHY, "st3": InterpolationKernel.ST3, "gauss": InterpolationKernel.GAUSS}[sys.argv[4]]
over_smooth = float(sys.argv[5])
center_shrink = sys.argv[6] == "1"
local_frac = float(sys.argv[7]) if len(sys.argv) > 7 else None
use_interp = (sys.argv[8] != "0") if len(sys.argv) > 8 else True

rng = np.random.default_rng(7)
ev = np.exp(rng.uniform(np.log(1e-2), 0.0, d))
Q, _ = np.linalg.qr(rng.normal(size=(d, d)))
S = Q @ np.diag(ev) @ Q.T
L = np.linalg.cholesky(S)
Sinv = np.linalg.inv(S)


def log_prob(theta, _args):
    theta = np.asarray(theta)
    return -0.5 * theta @ Sinv @ theta


apes = APES(
    nwalkers=nwalkers,
    ndim=d,
    model=None,
    log_prob=log_prob,
    interpolation_kernel=kernel,
    over_smooth=over_smooth,
    local_fraction=local_frac,
    center_shrink=center_shrink,
    use_interpolation=use_interp,
)
x0 = rng.normal(size=(nwalkers, d)) @ L.T  # stationary start
t0 = time.perf_counter()
apes.run_mcmc(x0, niter)
dt = time.perf_counter() - t0

acc = apes.esmcmc.get_accept_ratio()
mcat = apes.esmcmc.peek_catalog()
mcat.estimate_autocorrelation_tau(False)
tau = np.array(mcat.peek_autocorrelation_tau().dup_array())
print(
    f"d={d} nwalkers={nwalkers} niter={niter} kernel={sys.argv[4]} h={over_smooth} "
    f"lf={local_frac} shrink={int(center_shrink)} interp={int(use_interp)} acc={acc:.3f} "
    f"tau_mean={tau[1:].mean():.1f} tau_max={tau[1:].max():.1f} time={dt:.0f}s",
    flush=True,
)

"""Use the walker's OWN prepared proposal objects to score the acceptance, so the
comparison with the chain's measured acceptance cannot differ by a setup detail."""

import sys
import numpy as np
from numcosmo_py import Ncm
from numcosmo_py.sampling.apes import APES
from numcosmo_py.interpolation.stats_dist import InterpolationKernel

Ncm.cfg_init()
d = int(sys.argv[1])
nw = int(sys.argv[2])
niter = int(sys.argv[3])
kname = sys.argv[4]
h = float(sys.argv[5])
shrink = sys.argv[6] == "1"
lf = float(sys.argv[7])
interp = sys.argv[8] != "0"
KMAP = {
    "cauchy": InterpolationKernel.CAUCHY,
    "st3": InterpolationKernel.ST3,
    "gauss": InterpolationKernel.GAUSS,
}

rng = np.random.default_rng(7)
ev = np.exp(rng.uniform(np.log(1e-2), 0.0, d))
Q, _ = np.linalg.qr(rng.normal(size=(d, d)))
S = Q @ np.diag(ev) @ Q.T
L = np.linalg.cholesky(S)
Sinv = np.linalg.inv(S)
m2lnp = lambda X: np.einsum("ij,jk,ik->i", X, Sinv, X)

apes = APES(
    nwalkers=nw,
    ndim=d,
    model=None,
    log_prob=lambda t, _: -0.5 * np.asarray(t) @ Sinv @ np.asarray(t),
    interpolation_kernel=KMAP[kname],
    over_smooth=h,
    local_fraction=lf,
    center_shrink=shrink,
    use_interpolation=interp,
)
walker = apes.esmcmc.peek_walker()
x0 = rng.normal(size=(nw, d)) @ L.T
apes.run_mcmc(x0, niter)
chain_acc = apes.esmcmc.get_accept_ratio()
mcat = apes.esmcmc.peek_catalog()
rows = np.array([mcat.peek_row(mcat.len() - nw + i).dup_array() for i in range(nw)])
ens = rows[:, 1:]
half = nw // 2
sd0, sd1 = walker.peek_sds()
NRNG = Ncm.RNG.seeded_new(None, 5)


def score(sd, X):
    """expected acceptance of the independence step for current points X"""
    n = len(X)
    r_x = -0.5 * (
        np.array([sd.eval_m2lnp(Ncm.Vector.new_array(x.tolist())) for x in X])
        - m2lnp(X)
    )
    v = Ncm.Vector.new(d)
    Y = np.empty((n, d))
    for i in range(n):
        sd.sample(v, NRNG)
        Y[i] = v.dup_array()
    r_y = -0.5 * (
        np.array([sd.eval_m2lnp(Ncm.Vector.new_array(y.tolist())) for y in Y])
        - m2lnp(Y)
    )
    return np.exp(np.minimum(0.0, r_x - r_y)).mean()


# sd0 is built from the SECOND half and proposes for the FIRST half
print(f"d={d} {kname} h={h} shrink={int(shrink)} interp={int(interp)}", flush=True)
print(f"  chain acceptance            = {chain_acc:.3f}", flush=True)
print(
    f"  walker sd0, x = first half  = {score(sd0, ens[:half]):.3f}   (nk={sd0.get_n_kernels()})",
    flush=True,
)
print(
    f"  walker sd1, x = second half = {score(sd1, ens[half:]):.3f}   (nk={sd1.get_n_kernels()})",
    flush=True,
)
Xi = rng.normal(size=(half, d)) @ L.T
print(f"  walker sd0, x = iid target  = {score(sd0, Xi):.3f}", flush=True)

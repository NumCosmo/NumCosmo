"""Why does the chain accept less than the one-step benchmark predicts?

Runs an APES chain, then builds the SAME proposal twice: once from the chain's own
half-ensemble and once from an i.i.d. sample of the target of the same size, and
measures the expected independence-MH acceptance of both. The difference is the cost
of the ensemble not being an i.i.d. sample of the target.
"""
import sys
import numpy as np

from numcosmo_py import Ncm
from numcosmo_py.sampling.apes import APES
from numcosmo_py.interpolation.stats_dist import InterpolationKernel

Ncm.cfg_init()

d = int(sys.argv[1]); nwalkers = int(sys.argv[2]); niter = int(sys.argv[3])
kname = sys.argv[4]; h = float(sys.argv[5]); shrink = sys.argv[6] == "1"
lf = float(sys.argv[7]); interp = sys.argv[8] != "0"
KMAP = {"cauchy": (1.0, InterpolationKernel.CAUCHY), "st3": (3.0, InterpolationKernel.ST3),
        "gauss": (0.0, InterpolationKernel.GAUSS)}
nu, kenum = KMAP[kname]

rng = np.random.default_rng(7)
ev = np.exp(rng.uniform(np.log(1e-2), 0.0, d))
Q, _ = np.linalg.qr(rng.normal(size=(d, d)))
S = Q @ np.diag(ev) @ Q.T
L = np.linalg.cholesky(S)
Sinv = np.linalg.inv(S)
m2lnp = lambda X: np.einsum("ij,jk,ik->i", X, Sinv, X)

apes = APES(nwalkers=nwalkers, ndim=d, model=None, log_prob=lambda t, _: -0.5 * np.asarray(t) @ Sinv @ np.asarray(t),
            interpolation_kernel=kenum, over_smooth=h, local_fraction=lf,
            center_shrink=shrink, use_interpolation=interp)
x0 = rng.normal(size=(nwalkers, d)) @ L.T
apes.run_mcmc(x0, niter)
acc_chain = apes.esmcmc.get_accept_ratio()
mcat = apes.esmcmc.peek_catalog()
rows = np.array([mcat.peek_row(mcat.len() - nwalkers + i).dup_array() for i in range(nwalkers)])
ens = rows[:, 1:]            # row layout is [m2lnL, theta_0 .. theta_{d-1}]
ens_m2 = rows[:, 0]
assert ens.shape == (nwalkers, d), ens.shape
half = nwalkers // 2
NCMRNG = Ncm.RNG.seeded_new(None, 99)


def build_and_score(train, train_m2, label, x_points=None):
    kern = Ncm.StatsDistKernelGauss.new(d) if nu <= 0 else Ncm.StatsDistKernelST.new(d, nu)
    sd = Ncm.StatsDistVKDE.new(kern, Ncm.StatsDistCV.NONE)
    sd.set_over_smooth(h); sd.set_local_frac(lf); sd.set_center_shrink(shrink)
    for r in train:
        sd.add_obs(Ncm.Vector.new_array(r.tolist()))
    if interp:
        sd.prepare_interp(Ncm.Vector.new_array(train_m2.tolist()))
    else:
        sd.prepare()
    ntest = 800
    # x ~ pi (i.i.d.) unless the caller supplies the chain's own walker positions
    Xt = rng.normal(size=(ntest, d)) @ L.T if x_points is None else x_points[:ntest]
    ntest = len(Xt)
    r_test = -0.5 * (np.array([sd.eval_m2lnp(Ncm.Vector.new_array(x.tolist())) for x in Xt]) - m2lnp(Xt))
    v = Ncm.Vector.new(d)
    Y = np.empty((ntest, d))
    for i in range(ntest):
        sd.sample(v, NCMRNG)
        Y[i] = v.dup_array()
    r_Y = -0.5 * (np.array([sd.eval_m2lnp(Ncm.Vector.new_array(y.tolist())) for y in Y]) - m2lnp(Y))
    acc = np.exp(np.minimum(0.0, r_test - r_Y)).mean()
    # how degenerate is the training set?
    dist = np.linalg.norm(train[:, None, :] - train[None, :, :], axis=2) + np.eye(len(train)) * 1e9
    dmin = dist.min(axis=1)
    print(f"{label}: acc_1step={acc:.3f} std_r={np.std(r_test):.3f} "
          f"min_pair_dist={dmin.min():.3g} median_nn={np.median(dmin):.3g} "
          f"n_dup={(dmin < 1e-10).sum()}", flush=True)
    return acc


print(f"d={d} kernel={kname} h={h} shrink={int(shrink)} interp={int(interp)} chain_acc={acc_chain:.3f}", flush=True)
a_ens = build_and_score(ens[:half], ens_m2[:half], "train=ens  x=iid   ")
# the decisive comparison: same proposal, but x taken from the half the chain actually
# updates with it, instead of fresh i.i.d. target draws
a_real = build_and_score(ens[:half], ens_m2[:half], "train=ens  x=walkers", x_points=ens[half:])
Xiid = rng.normal(size=(half, d)) @ L.T
a_iid = build_and_score(Xiid, m2lnp(Xiid), "train=iid  x=iid   ")

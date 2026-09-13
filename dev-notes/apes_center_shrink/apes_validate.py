"""Correctness validation of NcmStatsDist centre shrinkage, emphasising high dimension.

Three independent checks per configuration:

1. exact  : reconstruct the mixture density in numpy from the public API
            (peek_center_array, peek_cov_decomp, peek_weights, get_lnnorm,
            get_href, get_center_shrink_factor) and compare with eval()/eval_m2lnp().
            Confirms the centre/bandwidth/normalisation bookkeeping to round-off.
2. ks     : project the draws on random directions u and compare with the EXACT
            one-dimensional law of the mixture along u (marginals of a Gaussian /
            Student-t kernel stay Gaussian / Student-t with the same nu). Tests that
            sample() draws from the density eval() returns -- an inconsistency here
            breaks detailed balance. Unlike a global importance-sampling check this
            stays sharp in high dimension.
3. cov    : covariance of draws from the mixture vs the sample covariance. With
            centre shrinkage the two must agree; without it the mixture is inflated
            by 1 + h^2 s^2.
"""

import argparse
import sys
import warnings

import numpy as np
from scipy import stats as sstats

from numcosmo_py import Ncm

Ncm.cfg_init()


def make_target(name, d, rng):
    if name == "gauss":
        ev = np.exp(rng.uniform(np.log(1e-2), 0.0, d))
        Q, _ = np.linalg.qr(rng.normal(size=(d, d)))
        S = Q @ np.diag(ev) @ Q.T
        L = np.linalg.cholesky(S)
        return lambda n: rng.normal(size=(n, d)) @ L.T, (
            lambda X: np.einsum("ij,jk,ik->i", X, np.linalg.inv(S), X)
        )
    if name == "funnel":
        def sample(n):
            v = 3.0 * rng.normal(size=n)
            x = rng.normal(size=(n, d - 1)) * np.exp(0.5 * v)[:, None]
            return np.concatenate([v[:, None], x], axis=1)

        def m2lnp(X):
            v = X[:, 0]
            return v**2 / 9.0 + np.sum(X[:, 1:] ** 2, axis=1) * np.exp(-v) + (d - 1) * v

        return sample, m2lnp
    raise ValueError(name)


def build(train, m2lnp, d, method, nu, h, local_frac, shrink, interp):
    kernel = Ncm.StatsDistKernelGauss.new(d) if nu <= 0 else Ncm.StatsDistKernelST.new(d, nu)
    if method == "kde":
        sd = Ncm.StatsDistKDE.new(kernel, Ncm.StatsDistCV.NONE)
    else:
        sd = Ncm.StatsDistVKDE.new(kernel, Ncm.StatsDistCV.NONE)
        sd.set_local_frac(local_frac)
    sd.set_over_smooth(h)
    sd.set_center_shrink(shrink)
    for row in train:
        sd.add_obs(Ncm.Vector.new_array(row.tolist()))
    if interp:
        sd.prepare_interp(Ncm.Vector.new_array(m2lnp.tolist()))
    else:
        sd.prepare()
    return sd, kernel


def reconstruct(sd, nu, d, X):
    """Mixture density at rows of X, computed in numpy from the public API."""
    nk = sd.get_n_kernels()
    w = np.array(sd.peek_weights().dup_array())
    a = sd.get_center_shrink_factor()
    href = a * sd.get_href()
    C = np.array([v.dup_array() for v in sd.peek_center_array()])
    lnnorm = np.array([sd.get_lnnorm(i) for i in range(nk)])
    # cov_decomp is upper triangular U with C_i = U^T U; the strict lower part is
    # scratch and must be dropped. chi2 = |U^-T (x - c)|^2 / href^2.
    UinvT = [np.linalg.inv(np.triu(np.array(sd.peek_cov_decomp(i).dup_array()).reshape(d, d))).T
             for i in range(nk)]
    out = np.empty(len(X))
    for j, x in enumerate(X):
        chi2 = np.empty(nk)
        for i in range(nk):
            y = UinvT[i] @ (x - C[i])
            chi2[i] = y @ y / href**2
        K = np.exp(-0.5 * chi2) if nu <= 0 else np.exp(-0.5 * (nu + d) * np.log1p(chi2 / nu))
        out[j] = np.sum(w * K * np.exp(-lnnorm))
    return out, href, a


def ks_projection(sd, nu, d, Y, ndir, rng):
    """Max KS distance between the projected draws and the exact projected mixture."""
    nk = sd.get_n_kernels()
    w = np.array(sd.peek_weights().dup_array())
    w = w / w.sum()
    a = sd.get_center_shrink_factor()
    href = a * sd.get_href()
    C = np.array([v.dup_array() for v in sd.peek_center_array()])
    U = [np.triu(np.array(sd.peek_cov_decomp(i).dup_array()).reshape(d, d)) for i in range(nk)]
    worst = 0.0
    for _ in range(ndir):
        u = rng.normal(size=d)
        u /= np.linalg.norm(u)
        # scale of kernel i along u: href^2 u^T C_i u with C_i = U_i^T U_i
        sig = href * np.array([np.linalg.norm(U[i] @ u) for i in range(nk)])
        mu_i = C @ u
        z = np.sort(Y @ u)
        if nu <= 0:
            cdf = np.array([np.sum(w * sstats.norm.cdf(zz, loc=mu_i, scale=sig)) for zz in z])
        else:
            cdf = np.array([np.sum(w * sstats.t.cdf(zz, nu, loc=mu_i, scale=sig)) for zz in z])
        n = len(z)
        emp_hi = np.arange(1, n + 1) / n
        emp_lo = np.arange(0, n) / n
        worst = max(worst, np.max(np.maximum(np.abs(cdf - emp_hi), np.abs(cdf - emp_lo))))
    return worst


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--target", default="gauss")
    p.add_argument("--d", type=int, nargs="+", default=[10])
    p.add_argument("--nmult", type=int, default=40, help="training points = nmult * d")
    p.add_argument("--method", nargs="+", default=["vkde"])
    p.add_argument("--nu", type=float, nargs="+", default=[0.0, 3.0])
    p.add_argument("--h", type=float, nargs="+", default=[2.0])
    p.add_argument("--local-frac", type=float, default=0.2)
    p.add_argument("--interp", type=int, nargs="+", default=[0, 1])
    p.add_argument("--nexact", type=int, default=4)
    p.add_argument("--nks", type=int, default=4000, help="draws used for the KS projection check")
    p.add_argument("--ndir", type=int, default=3, help="random projection directions")
    p.add_argument("--ncov", type=int, default=40000, help="draws used for the covariance check")
    p.add_argument("--seed", type=int, default=11)
    a = p.parse_args()

    print("target,d,n,method,nu,h,interp,shrink,kappa,a,a_exp,href,max_rel_eval,max_rel_m2lnp,ks,ks_crit,cov_trace,cov_frob,warn", flush=True)
    for d in a.d:
        rng = np.random.default_rng(a.seed + d)
        sampler, m2lnp_f = make_target(a.target, d, rng)
        n = a.nmult * d
        train = sampler(n)
        m2 = m2lnp_f(train)
        Sig = np.cov(train.T)
        mu = train.mean(axis=0)
        Siginv = np.linalg.inv(Sig)
        # narrow reference density for the normalisation check (lighter tails than p)
        Lf = np.linalg.cholesky(0.5 * Sig)
        logdetf = 2.0 * np.sum(np.log(np.diag(Lf)))
        Finv = np.linalg.inv(0.5 * Sig)

        def logf(X):
            D = X - mu
            return -0.5 * (np.einsum("ij,jk,ik->i", D, Finv, D) + logdetf + d * np.log(2 * np.pi))

        for method in a.method:
            for nu in a.nu:
                for h in a.h:
                    for interp in a.interp:
                        for shrink in (0, 1):
                            with warnings.catch_warnings(record=True) as wl:
                                try:
                                    sd, _k = build(train, m2, d, method, nu, h, a.local_frac, bool(shrink), bool(interp))
                                except Exception as e:  # noqa
                                    print(f"# FAILED d={d} {method} nu={nu} h={h} interp={interp} shrink={shrink}: {e}", file=sys.stderr, flush=True)
                                    continue
                            # 1. exact reconstruction
                            Xt = sampler(a.nexact)
                            ref, href, af = reconstruct(sd, nu, d, Xt)
                            got = np.array([sd.eval(Ncm.Vector.new_array(x.tolist())) for x in Xt])
                            got2 = np.array([np.exp(-0.5 * sd.eval_m2lnp(Ncm.Vector.new_array(x.tolist()))) for x in Xt])
                            rel1 = np.max(np.abs(got / ref - 1.0))
                            rel2 = np.max(np.abs(got2 / ref - 1.0))
                            # 2. normalisation / sample-eval consistency
                            v = Ncm.Vector.new(d)
                            Y = np.empty((a.ncov, d))
                            for i in range(a.ncov):
                                sd.sample(v, RNG)
                                Y[i] = v.dup_array()
                            ks = ks_projection(sd, nu, d, Y[: a.nks], a.ndir, rng)
                            ks_crit = 1.63 / np.sqrt(a.nks)  # 1% two-sided KS threshold
                            kappa = 1.0 if nu <= 0 else (nu / (nu - 2.0) if nu > 2.0 else float("inf"))
                            s2 = (1.0 / af**2 - 1.0) / (kappa * sd.get_href() ** 2) if shrink and np.isfinite(kappa) else float("nan")
                            a_exp = af
                            # 3. covariance identity
                            Cest = np.cov(Y.T)
                            ctr = np.trace(Cest @ Siginv) / d
                            cfr = np.linalg.norm(Cest - Sig) / np.linalg.norm(Sig)
                            print(f"{a.target},{d},{n},{method},{nu},{h},{interp},{shrink},{kappa:.3f},{af:.4f},{a_exp:.4f},{href:.4f},"
                                  f"{rel1:.2e},{rel2:.2e},{ks:.4f},{ks_crit:.4f},{ctr:.4f},{cfr:.4f},{len(wl)}", flush=True)


RNG = Ncm.RNG.seeded_new(None, 1234)

if __name__ == "__main__":
    main()

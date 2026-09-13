"""End-to-end APES on hard targets, with target-specific diagnostics.

Targets
  gauss     correlated Gaussian (reference)
  funnel    Neal funnel: v ~ N(0, 3^2), x_i | v ~ N(0, e^v)
  gaussmix  two well-separated Gaussian modes (the case where matching the GLOBAL
            covariance is the wrong criterion for centre shrinkage)
  banana    Haario twisted Gaussian, curved degeneracy in d dimensions
  studentt  multivariate Student-t target, nu = 3 (heavy tails)
  rosenbrock  the 2d Rosenbrock of the APES paper
"""
import sys
import time

import numpy as np

from numcosmo_py import Ncm
from numcosmo_py.sampling.apes import APES
from numcosmo_py.interpolation.stats_dist import InterpolationKernel

Ncm.cfg_init()

target = sys.argv[1]
d = int(sys.argv[2]); nw = int(sys.argv[3]); niter = int(sys.argv[4])
kname = sys.argv[5]; h = float(sys.argv[6]); shrink = sys.argv[7] == "1"
lf = float(sys.argv[8]) if sys.argv[8] != "-" else None
interp = sys.argv[9] != "0"
KMAP = {"cauchy": InterpolationKernel.CAUCHY, "st3": InterpolationKernel.ST3, "gauss": InterpolationKernel.GAUSS}

rng = np.random.default_rng(7)
SEP = 8.0  # mode separation for gaussmix, in units of the marginal sigma


def make(target, d):
    if target == "gauss":
        ev = np.exp(rng.uniform(np.log(1e-2), 0.0, d))
        Q, _ = np.linalg.qr(rng.normal(size=(d, d)))
        S = Q @ np.diag(ev) @ Q.T
        L = np.linalg.cholesky(S); Si = np.linalg.inv(S)
        return (lambda n: rng.normal(size=(n, d)) @ L.T,
                lambda t: -0.5 * t @ Si @ t)
    if target == "funnel":
        def samp(n):
            v = 3.0 * rng.normal(size=n)
            return np.concatenate([v[:, None], rng.normal(size=(n, d - 1)) * np.exp(0.5 * v)[:, None]], axis=1)

        def lp(t):
            v = t[0]
            return -0.5 * (v * v / 9.0 + np.dot(t[1:], t[1:]) * np.exp(-v) + (d - 1) * v)
        return samp, lp
    if target == "gaussmix":
        mu = np.zeros(d); mu[0] = SEP / 2.0

        def samp(n):
            s = rng.integers(0, 2, n) * 2 - 1
            return rng.normal(size=(n, d)) + s[:, None] * mu

        def lp(t):
            a = -0.5 * np.dot(t - mu, t - mu)
            b = -0.5 * np.dot(t + mu, t + mu)
            m = max(a, b)
            return m + np.log(np.exp(a - m) + np.exp(b - m)) - np.log(2.0)
        return samp, lp
    if target == "banana":
        b = 0.03
        var1 = 100.0

        def samp(n):
            z = rng.normal(size=(n, d)); z[:, 0] *= np.sqrt(var1)
            x = z.copy(); x[:, 1] = z[:, 1] + b * (z[:, 0] ** 2 - var1)
            return x

        def lp(t):
            y1 = t[0]
            y2 = t[1] - b * (y1 * y1 - var1)
            return -0.5 * (y1 * y1 / var1 + y2 * y2 + np.dot(t[2:], t[2:]))
        return samp, lp
    if target == "studentt":
        nu_t = 3.0
        ev = np.exp(rng.uniform(np.log(1e-2), 0.0, d))
        Q, _ = np.linalg.qr(rng.normal(size=(d, d)))
        S = Q @ np.diag(ev) @ Q.T
        L = np.linalg.cholesky(S); Si = np.linalg.inv(S)

        def samp(n):
            g = rng.normal(size=(n, d)) @ L.T
            w = np.sqrt(nu_t / rng.chisquare(nu_t, n))
            return g * w[:, None]

        def lp(t):
            return -0.5 * (nu_t + d) * np.log1p(t @ Si @ t / nu_t)
        return samp, lp
    if target == "rosenbrock":
        def samp(n):
            x1 = 1.0 + np.sqrt(10.0) * rng.normal(size=n)
            return np.stack([x1, x1**2 + np.sqrt(0.1) * rng.normal(size=n)], axis=1)

        def lp(t):
            return -0.05 * (100.0 * (t[1] - t[0] ** 2) ** 2 + (1.0 - t[0]) ** 2)
        return samp, lp
    raise ValueError(target)


sampler, log_prob = make(target, d)
apes = APES(nwalkers=nw, ndim=d, model=None, log_prob=lambda t, _: log_prob(np.asarray(t)),
            interpolation_kernel=KMAP[kname], over_smooth=h, local_fraction=lf,
            center_shrink=shrink, use_interpolation=interp)
x0 = sampler(nw)
t0 = time.perf_counter()
apes.run_mcmc(x0, niter)
dt = time.perf_counter() - t0
acc = apes.esmcmc.get_accept_ratio()
sd0, _sd1 = apes.esmcmc.peek_walker().peek_sds()
a_shrink = sd0.get_center_shrink_factor()
href_applied = sd0.get_href()
mcat = apes.esmcmc.peek_catalog()
mcat.estimate_autocorrelation_tau(False)
tau = np.array(mcat.peek_autocorrelation_tau().dup_array())[1:]
extra = ""
rows = np.array([mcat.peek_row(mcat.len() - nw + i).dup_array() for i in range(nw)])
ens = rows[:, 1:]
if target == "gaussmix":
    # occupancy, mixing between modes, and whether the modes are still RESOLVED:
    # for the true target |x0| ~ 4 +- 1 and only 2.3 % of the mass has |x0| < 2
    frac = np.mean(ens[:, 0] > 0.0)
    first = np.array([mcat.peek_row(i).dup_array() for i in range(nw)])[:, 1]
    hops = np.mean(np.sign(first) != np.sign(ens[:, 0]))
    gap = np.mean(np.abs(ens[:, 0]) < 2.0)
    ref = sampler(20000)
    extra = (f" mode_frac={frac:.3f} changed_mode={hops:.3f} tau_x0={tau[0]:.1f}"
             f" mean|x0|={np.mean(np.abs(ens[:,0])):.2f}(true {np.mean(np.abs(ref[:,0])):.2f})"
             f" frac_in_gap={gap:.3f}(true {np.mean(np.abs(ref[:,0])<2.0):.3f})")
elif target == "funnel":
    ref = sampler(20000)
    extra = (f" tau_v={tau[0]:.1f} v_mean={ens[:,0].mean():.2f}(true {ref[:,0].mean():.2f})"
             f" v_sd={ens[:,0].std():.2f}(true {ref[:,0].std():.2f})"
             f" v_range=[{ens[:,0].min():.1f},{ens[:,0].max():.1f}]")
print(f"{target} d={d} nw={nw} it={niter} {kname} h={h} lf={lf} shrink={int(shrink)} interp={int(interp)} "
      f"a={a_shrink:.3f} href={href_applied:.3f} acc={acc:.3f} tau_mean={np.mean(tau):.1f} tau_max={np.max(tau):.1f} time={dt:.0f}s{extra}", flush=True)

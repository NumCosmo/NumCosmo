"""Standalone benchmark of the APES approximation quality.

Builds the NcmStatsDistVKDE interpolant from an i.i.d. sample of the target
(the converged half-ensemble) and measures, on independent test points, the
scatter of ln(pi_tilde/pi) and the expected acceptance of the independence
MH step  alpha = min(1, pi(y) pt(x) / (pi(x) pt(y))),  x ~ pi, y ~ pi_tilde.
"""

import argparse
import sys
import time

import numpy as np

from numcosmo_py import Ncm

Ncm.cfg_init()


# ----------------------------------------------------------------------------- targets
class Gauss:
    def __init__(self, d, rng):
        self.d = d
        # eigenvalues log-uniform in [1e-2, 1], random rotation
        ev = np.exp(rng.uniform(np.log(1e-2), 0.0, d))
        Q, _ = np.linalg.qr(rng.normal(size=(d, d)))
        self.S = Q @ np.diag(ev) @ Q.T
        self.L = np.linalg.cholesky(self.S)
        self.Sinv = np.linalg.inv(self.S)

    def sample(self, n, rng):
        return rng.normal(size=(n, self.d)) @ self.L.T

    def m2lnp(self, X):
        return np.einsum("ij,jk,ik->i", X, self.Sinv, X)


class Rosenbrock:
    d = 2

    def sample(self, n, rng):
        x1 = 1.0 + np.sqrt(10.0) * rng.normal(size=n)
        x2 = x1**2 + np.sqrt(0.1) * rng.normal(size=n)
        return np.stack([x1, x2], axis=1)

    def m2lnp(self, X):
        x1, x2 = X[:, 0], X[:, 1]
        return 0.1 * (100.0 * (x2 - x1**2) ** 2 + (1.0 - x1) ** 2)


class Funnel:
    """v ~ N(0, 3^2), x_i | v ~ N(0, e^v), i = 1..d-1."""

    def __init__(self, d):
        self.d = d

    def sample(self, n, rng):
        v = 3.0 * rng.normal(size=n)
        x = rng.normal(size=(n, self.d - 1)) * np.exp(0.5 * v)[:, None]
        return np.concatenate([v[:, None], x], axis=1)

    def m2lnp(self, X):
        v = X[:, 0]
        x = X[:, 1:]
        return v**2 / 9.0 + np.sum(x**2, axis=1) * np.exp(-v) + (self.d - 1) * v


def make_target(name, d, rng):
    if name == "gauss":
        return Gauss(d, rng)
    if name == "rosen":
        return Rosenbrock()
    if name == "funnel":
        return Funnel(d)
    raise ValueError(name)


# ----------------------------------------------------------------------------- build
CV = {
    "none": Ncm.StatsDistCV.NONE,
    "split": Ncm.StatsDistCV.SPLIT,
    "split_m2lnp": Ncm.StatsDistCV.SPLIT_M2LNP,
}


def rot_h(nu, n, d):
    if nu <= 0:
        return (4.0 / (n * (d + 2.0))) ** (1.0 / (d + 4.0))
    nu = max(nu, 3.0)
    return (
        16.0
        * (nu - 2) ** 2
        * (1 + d + nu)
        * (3 + d + nu)
        / ((2 + d) * (d + nu) * (2 + d + nu) * (d + 2 * nu) * (2 + d + 2 * nu) * n)
    ) ** (1.0 / (d + 4.0))


def local_scale2(sd, train):
    """mean_i tr(C_i Sigma^-1)/d : kernel covariance scale relative to the sample covariance."""
    d = train.shape[1]
    Sinv = np.linalg.inv(np.cov(train.T))
    acc = 0.0
    nk = sd.get_n_kernels()
    for i in range(nk):
        U = np.array(sd.peek_cov_decomp(i).dup_array()).reshape(d, d)
        U = np.triu(U)
        C = U.T @ U
        acc += np.trace(C @ Sinv) / d
    return acc / nk


def build(
    train,
    m2lnp_train,
    nu,
    h,
    local_frac,
    cv,
    split_frac,
    shrink,
    interp,
    kde=False,
    cshrink=False,
    target=None,
    native=False,
):
    d = train.shape[1]
    if cshrink:
        if kde:
            s2 = 1.0
        else:
            sd0, keep0, _ = build(
                train, m2lnp_train, nu, h, local_frac, "none", split_frac, shrink, False
            )
            s2 = local_scale2(sd0, train)
        a = 1.0 / np.sqrt(1.0 + h * h * s2)
        mu = train.mean(axis=0)
        train = mu + a * (train - mu)
        m2lnp_train = target.m2lnp(train)
    if nu > 0:
        kernel = Ncm.StatsDistKernelST.new(d, nu)
    else:
        kernel = Ncm.StatsDistKernelGauss.new(d)
    if kde:
        sd = Ncm.StatsDistKDE.new(kernel, CV[cv])
        sd.set_over_smooth(h / rot_h(nu, train.shape[0], d))
    else:
        sd = Ncm.StatsDistVKDE.new(kernel, CV[cv])
        sd.set_over_smooth(h)
        sd.set_local_frac(local_frac)
    sd.set_split_frac(split_frac)
    if shrink != 0.0:
        raise SystemExit(
            "NcmStatsDist:shrink (weight floor) was removed on 2026-09-16; only --shrink 0 is valid"
        )
    if native:
        sd.set_center_shrink(True)
    keep = []
    for row in train:
        v = Ncm.Vector.new_array(row.tolist())
        keep.append(v)
        sd.add_obs(v)
    t0 = time.perf_counter()
    if interp:
        sd.prepare_interp(Ncm.Vector.new_array(m2lnp_train.tolist()))
    else:
        sd.prepare()
    dt = time.perf_counter() - t0
    return sd, keep, dt


def evaluate(sd, target, test, m2lnp_test, nprop, ncm_rng):
    d = test.shape[1]
    m2lnq_test = np.array(
        [sd.eval_m2lnp(Ncm.Vector.new_array(row.tolist())) for row in test]
    )
    r_test = -0.5 * (m2lnq_test - m2lnp_test)  # ln(pt/pi) + const
    v = Ncm.Vector.new(d)
    Y = np.empty((nprop, d))
    for i in range(nprop):
        sd.sample(v, ncm_rng)
        Y[i] = v.dup_array()
    m2lnq_Y = np.array([sd.eval_m2lnp(Ncm.Vector.new_array(row.tolist())) for row in Y])
    r_Y = -0.5 * (m2lnq_Y - target.m2lnp(Y))
    n = min(len(r_test), nprop)
    log_alpha = np.minimum(0.0, r_test[:n] - r_Y[:n])
    acc = np.exp(log_alpha).mean()
    w = np.array(sd.peek_weights().dup_array())
    nz = np.mean(w > 1e-3 / len(w))
    return dict(
        std_r=float(np.std(r_test)),
        std_rY=float(np.std(r_Y)),
        acc=float(acc),
        nz=float(nz),
        href=float(sd.get_href()),
        os=float(sd.get_over_smooth()),
        a=(
            float(sd.get_center_shrink_factor())
            if hasattr(sd, "get_center_shrink_factor")
            else 1.0
        ),
        nk=int(sd.get_n_kernels()),
    )


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--target", default="gauss")
    p.add_argument("--d", type=int, default=10)
    p.add_argument("--n", type=int, default=500, help="training points (half ensemble)")
    p.add_argument("--ntest", type=int, default=1000)
    p.add_argument("--nu", type=float, nargs="+", default=[1.0])
    p.add_argument("--h", type=float, nargs="+", default=[0.2])
    p.add_argument("--local-frac", type=float, nargs="+", default=[0.05])
    p.add_argument("--cv", default="none")
    p.add_argument("--split-frac", type=float, nargs="+", default=[0.5])
    p.add_argument("--shrink", type=float, nargs="+", default=[0.0])
    p.add_argument("--no-interp", action="store_true")
    p.add_argument(
        "--kde",
        action="store_true",
        help="global-covariance KDE instead of VKDE; h in units of the sample covariance",
    )
    p.add_argument(
        "--cshrink",
        action="store_true",
        help="shrink kernel centres toward the mean, a = 1/sqrt(1+h^2 s^2) (python emulation)",
    )
    p.add_argument(
        "--native-cshrink",
        action="store_true",
        help="use NcmStatsDist center-shrink property",
    )
    p.add_argument("--seed", type=int, default=1)
    a = p.parse_args()

    rng = np.random.default_rng(a.seed)
    ncm_rng = Ncm.RNG.seeded_new(None, a.seed)
    target = make_target(a.target, a.d, rng)
    d = target.d
    train = target.sample(a.n, rng)
    test = target.sample(a.ntest, rng)
    m2_train = target.m2lnp(train)
    m2_test = target.m2lnp(test)

    print(
        "target,d,n,kde,cshrink,interp,cv,split,lf,nu,h,shrink,nk,os,href,a,nz,std_r,std_rY,acc,t_prep",
        flush=True,
    )
    for lf in a.local_frac:
        for nu in a.nu:
            for h in a.h:
                for sf in a.split_frac:
                    for sh in a.shrink:
                        try:
                            sd, keep, dt = build(
                                train,
                                m2_train,
                                nu,
                                h,
                                lf,
                                a.cv,
                                sf,
                                sh,
                                not a.no_interp,
                                kde=a.kde,
                                cshrink=a.cshrink,
                                target=target,
                                native=a.native_cshrink,
                            )
                            r = evaluate(sd, target, test, m2_test, a.ntest, ncm_rng)
                        except Exception as e:  # noqa
                            print(
                                f"# FAILED nu={nu} h={h} lf={lf} sf={sf}: {e}",
                                file=sys.stderr,
                                flush=True,
                            )
                            continue
                        print(
                            f"{a.target},{d},{a.n},{int(a.kde)},{int(a.cshrink) + 2 * int(a.native_cshrink)},{int(not a.no_interp)},{a.cv},{sf},{lf},{nu},{h},{sh},"
                            f"{r['nk']},{r['os']:.4g},{r['href']:.4g},{r['a']:.4f},{r['nz']:.3f},{r['std_r']:.4f},{r['std_rY']:.4f},"
                            f"{r['acc']:.4f},{dt:.2f}",
                            flush=True,
                        )


if __name__ == "__main__":
    main()

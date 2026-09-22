"""Set up Planck 2018 TT + LCDM and time one likelihood evaluation."""
import time
import numpy as np
from numcosmo_py import Nc, Ncm
from numcosmo_py.experiments.planck18 import generate_planck18_tt

Ncm.cfg_init()
exp, _oa = generate_planck18_tt()
print("experiment keys:", exp.keys())
likelihood = exp.get("likelihood")
mset = exp.get("model-set")
mset.prepare_fparam_map()
n = mset.fparams_len()
print(f"free parameters: {n}")
for i in range(n):
    print(f"  {i:2d} {mset.fparam_full_name(i):28s} "
          f"[{mset.fparam_get_lower_bound(i):.4g}, {mset.fparam_get_upper_bound(i):.4g}] "
          f"scale={mset.fparam_get_scale(i):.4g}")
fit = Ncm.Fit.factory(Ncm.FitType.NLOPT, "ln-neldermead", likelihood, mset,
                      Ncm.FitGradType.NUMDIFF_FORWARD)
t0 = time.perf_counter(); v0 = fit.m2lnL_val(); t1 = time.perf_counter()
print(f"first m2lnL = {v0:.4f} in {t1-t0:.2f}s")
ts = []
v = Ncm.Vector.new(n)
mset.fparams_get_vector(v)
p0 = np.array(v.dup_array())
rng = np.random.default_rng(0)
for i in range(4):
    p = p0 + 1e-3 * rng.normal(size=n) * np.array([mset.fparam_get_scale(j) for j in range(n)])
    mset.fparams_set_vector(Ncm.Vector.new_array(p.tolist()))
    t0 = time.perf_counter(); v = fit.m2lnL_val(); t1 = time.perf_counter()
    ts.append(t1 - t0)
    print(f"  perturbed m2lnL = {v:.4f} in {t1-t0:.2f}s")
print(f"median eval time: {np.median(ts):.2f}s")

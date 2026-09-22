"""Phase A: adapt an ensemble on the real Planck 2018 TT + LCDM posterior.

Runs the APES sampler in its current default configuration from the prior sampler.
Produces (a) the baseline acceptance including burn-in and (b) a catalog that phase B
reuses to seed every candidate configuration, so the comparison is controlled.
"""
import sys
import time

import numpy as np

from numcosmo_py import Ncm
from numcosmo_py.experiments.planck18 import generate_planck18_tt
from numcosmo_py.sampling.esmcmc import create_esmcmc, WalkerTypes
from numcosmo_py.interpolation.stats_dist import InterpolationMethod, InterpolationKernel

Ncm.cfg_init()

nw = int(sys.argv[1]); niter = int(sys.argv[2]); lf = float(sys.argv[3])

exp, _oa = generate_planck18_tt()
likelihood = exp.get("likelihood")
mset = exp.get("model-set")
mset.prepare_fparam_map()
print(f"# Planck18 TT + LCDM, {mset.fparams_len()} free parameters, "
      f"{nw} walkers, local_frac {lf}", flush=True)

esmcmc = create_esmcmc(
    likelihood, mset, "planckA",
    verbose=False, sampler=WalkerTypes.APES,
    interpolation_method=InterpolationMethod.VKDE,
    interpolation_kernel=InterpolationKernel.CAUCHY,   # current default
    nwalkers=nw, use_threads=True, use_apes_threads=False,
    over_smooth=1.1, local_fraction=lf,
    use_apes_interpolation=True, use_apes_center_shrink=False,
    init_sampling_scale=1.0,
)
esmcmc.start_run()
t0 = time.perf_counter()
for it in range(niter):
    esmcmc.run(it + 1)
    acc = esmcmc.get_accept_ratio()
    lacc = esmcmc.get_accept_ratio_last_update()
    el = time.perf_counter() - t0
    print(f"iter {it+1:3d}  acc_cumulative={acc:.4f}  acc_last={lacc:.4f}  "
          f"elapsed={el/60:.1f} min  ({(it+1)*nw/el:.2f} evals/s)", flush=True)
esmcmc.end_run()
print(f"# catalog: {esmcmc.peek_catalog().peek_filename()}", flush=True)

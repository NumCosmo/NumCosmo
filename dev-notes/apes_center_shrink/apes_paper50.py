"""The APES paper's own hardest case: the truncated Gaussian in d dimensions.

Uses the model and data objects from numcosmo_py.experiments.gauss_constraint, and
starts from the prior sampler as the paper does, so burn-in is included. Reports the
acceptance and the autocorrelation time after discarding the first half of the chain.
"""
import sys
import time

import numpy as np

from numcosmo_py import Ncm
from numcosmo_py.experiments.gauss_constraint import create_mset, create_data_object
from numcosmo_py.sampling.esmcmc import create_esmcmc, WalkerTypes
from numcosmo_py.interpolation.stats_dist import InterpolationMethod, InterpolationKernel

Ncm.cfg_init()

dim = int(sys.argv[1]); nw = int(sys.argv[2]); niter = int(sys.argv[3])
kname = sys.argv[4]; h = float(sys.argv[5]); shrink = sys.argv[6] == "1"
lf = float(sys.argv[7]); interp = sys.argv[8] != "0"
tag = sys.argv[9]
KMAP = {"cauchy": InterpolationKernel.CAUCHY, "st3": InterpolationKernel.ST3, "gauss": InterpolationKernel.GAUSS}

rng = Ncm.RNG.seeded_new(None, 0)
mset, _mgc = create_mset(dim)
dgc = create_data_object(mset, dim, rng)
likelihood = Ncm.Likelihood.new(Ncm.Dataset.new_array([dgc]))

esmcmc = create_esmcmc(
    likelihood, mset, f"paper{dim}d_{tag}",
    verbose=False, sampler=WalkerTypes.APES,
    interpolation_method=InterpolationMethod.VKDE,
    interpolation_kernel=KMAP[kname],
    nwalkers=nw, use_threads=True, use_apes_threads=True,
    over_smooth=h, local_fraction=lf,
    use_apes_interpolation=interp,
    use_apes_center_shrink=shrink,
    init_sampling_scale=1.0,
)
t0 = time.perf_counter()
esmcmc.start_run()
esmcmc.run(niter)
esmcmc.end_run()
dt = time.perf_counter() - t0
acc = esmcmc.get_accept_ratio()
fname = esmcmc.peek_catalog().peek_filename()
sd0, _ = esmcmc.peek_walker().peek_sds()
a_sh = sd0.get_center_shrink_factor()

# autocorrelation after discarding the first half as burn-in
burn = (niter // 2) * nw
mcat = Ncm.MSetCatalog.new_from_file_ro(fname, burn)
mcat.estimate_autocorrelation_tau(False)
tau = np.array(mcat.peek_autocorrelation_tau().dup_array())[1:]
print(f"paper-truncgauss d={dim} nw={nw} it={niter} {kname} h={h} lf={lf} shrink={int(shrink)} "
      f"interp={int(interp)} a={a_sh:.3f} acc={acc:.4f} tau_mean={np.mean(tau):.1f} "
      f"tau_max={np.max(tau):.1f} time={dt:.0f}s", flush=True)

"""Phase B: seed a fresh chain from the phase-A ensemble and measure the acceptance
of one candidate configuration. Every configuration starts from the same catalog."""
import sys
import time

from numcosmo_py import Ncm
from numcosmo_py.experiments.planck18 import generate_planck18_tt
from numcosmo_py.sampling.esmcmc import create_esmcmc, WalkerTypes
from numcosmo_py.interpolation.stats_dist import InterpolationMethod, InterpolationKernel

Ncm.cfg_init()

cat_file = sys.argv[1]
nw = int(sys.argv[2]); niter = int(sys.argv[3]); lf = float(sys.argv[4])
kname = sys.argv[5]; h = float(sys.argv[6]); shrink = sys.argv[7] == "1"
tag = sys.argv[8]
KMAP = {"cauchy": InterpolationKernel.CAUCHY, "st3": InterpolationKernel.ST3, "gauss": InterpolationKernel.GAUSS}

exp, _oa = generate_planck18_tt()
likelihood = exp.get("likelihood")
mset = exp.get("model-set")
mset.prepare_fparam_map()
start_mcat = Ncm.MSetCatalog.new_from_file_ro(cat_file, 0)

esmcmc = create_esmcmc(
    likelihood, mset, f"planckB_{tag}",
    verbose=False, sampler=WalkerTypes.APES,
    interpolation_method=InterpolationMethod.VKDE,
    interpolation_kernel=KMAP[kname],
    nwalkers=nw, use_threads=True, use_apes_threads=False,
    over_smooth=h, local_fraction=lf,
    use_apes_interpolation=True, use_apes_center_shrink=shrink,
    start_mcat=start_mcat,
)
esmcmc.start_run()
t0 = time.perf_counter()
esmcmc.run(niter)
esmcmc.end_run()
dt = time.perf_counter() - t0
sd0, _ = esmcmc.peek_walker().peek_sds()
print(f"PLANCK18-TT d=21 nw={nw} it={niter} {kname} h={h} lf={lf} shrink={int(shrink)} "
      f"a={sd0.get_center_shrink_factor():.3f} acc={esmcmc.get_accept_ratio():.4f} "
      f"time={dt/60:.1f}min", flush=True)

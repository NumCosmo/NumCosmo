"""APES on real but cheap cosmological likelihoods, built by combining NcmData objects.

Datasets
  snia        JLA SNIa alone (alpha, beta, M1, M2 free)
  snia_bao    JLA + BAO (SDSS DR12 + DESI DR2)
  full        JLA + BAO + H(z) (Moresco 2016)
Cosmology is NcHICosmoDEXcdm with Omega_x traded for Omega_k.
"""
import sys
import time

import numpy as np

from numcosmo_py import Nc, Ncm
from numcosmo_py.sampling.esmcmc import create_esmcmc, WalkerTypes
from numcosmo_py.interpolation.stats_dist import InterpolationMethod, InterpolationKernel

Ncm.cfg_init()

dataset = sys.argv[1]
nw = int(sys.argv[2]); niter = int(sys.argv[3])
kname = sys.argv[4]; h = float(sys.argv[5]); shrink = sys.argv[6] == "1"
lf = float(sys.argv[7]); interp = sys.argv[8] != "0"
tag = sys.argv[9] if len(sys.argv) > 9 else "run"
KMAP = {"cauchy": InterpolationKernel.CAUCHY, "st3": InterpolationKernel.ST3, "gauss": InterpolationKernel.GAUSS}


def build(dataset):
    """dataset is one of snia / snia_bao / full / wsplineN (N = number of w knots)."""
    if dataset.startswith("wspline"):
        nknots = int(dataset[len("wspline"):])
        cosmo = Nc.HICosmoDEWSpline.new(nknots, 2.0)
        cosmo.param_set_by_name("H0", 70.0)
        cosmo.param_set_by_name("Omegab", 0.05)
        cosmo.param_set_by_name("Omegac", 0.25)
        free = ["H0", "Omegac", "Omegab"]
        for k in range(nknots):
            cosmo.param_set_by_name(f"w_{k}", -1.0)
            free.append(f"w_{k}")
        for p in free:
            cosmo.param_set_ftype(cosmo.param_index_from_name(p)[1], Ncm.ParamType.FREE)
        dataset = "full"
    else:
        cosmo = Nc.HICosmoDEXcdm()
        cosmo.omega_x2omega_k()
        cosmo.param_set_by_name("H0", 70.0)
        cosmo.param_set_by_name("Omegab", 0.05)
        cosmo.param_set_by_name("Omegac", 0.25)
        cosmo.param_set_by_name("Omegak", 0.0)
        cosmo.param_set_by_name("w", -1.0)
        for p in ("H0", "Omegac", "Omegak", "w"):
            cosmo.param_set_ftype(cosmo.param_index_from_name(p)[1], Ncm.ParamType.FREE)

    dist = Nc.Distance(zf=3.0)
    snia_id = Nc.DataSNIAId.COV_JLA_SNLS3_SDSS_SYS_STAT
    snia_model = Nc.SNIADistCov.new_by_id(dist, snia_id)
    for p in ("alpha", "beta", "M1", "M2"):
        snia_model.param_set_ftype(snia_model.param_index_from_name(p)[1], Ncm.ParamType.FREE)
    snia = Nc.DataSNIACov.new_from_cat_id(snia_id, False)

    dset = Ncm.Dataset()
    dset.append_data(snia)
    if dataset in ("snia_bao", "full"):
        for bid in (Nc.DataBaoId.DTR_DHR_SDSS_DR12_2016_DR16_COMPATIBLE,
                    Nc.DataBaoId.DVR_DTDH_DESI_DR2_2025):
            dset.append_data(Nc.data_bao_create(dist, bid))
    if dataset == "full":
        dset.append_data(Nc.DataHubble.new_from_id(Nc.DataHubbleId.MORESCO2016_DR9_BC03))

    mset = Ncm.MSet()
    mset.set(cosmo)
    mset.set(snia_model)
    mset.prepare_fparam_map()
    return Ncm.Likelihood(dataset=dset), mset


likelihood, mset = build(dataset)
n = mset.fparams_len()
names = [mset.fparam_full_name(i) for i in range(n)]
print(f"# {dataset}: {n} free parameters: {', '.join(names)}", flush=True)

esmcmc = create_esmcmc(
    likelihood, mset, f"cosmo_{dataset}_{tag}",
    verbose=False, sampler=WalkerTypes.APES,
    interpolation_method=InterpolationMethod.VKDE,
    interpolation_kernel=KMAP[kname],
    nwalkers=nw, use_threads=True, use_apes_threads=False,
    over_smooth=h, local_fraction=lf,
    use_apes_interpolation=interp, use_apes_center_shrink=shrink,
    init_sampling_scale=1.0,
)
t0 = time.perf_counter()
esmcmc.start_run()
esmcmc.run(niter)
esmcmc.end_run()
dt = time.perf_counter() - t0
acc = esmcmc.get_accept_ratio()
sd0, _ = esmcmc.peek_walker().peek_sds()
fname = esmcmc.peek_catalog().peek_filename()
mcat = Ncm.MSetCatalog.new_from_file_ro(fname, (niter // 2) * nw)
mcat.estimate_autocorrelation_tau(False)
tau = np.array(mcat.peek_autocorrelation_tau().dup_array())[1:]
print(f"{dataset} d={n} nw={nw} it={niter} {kname} h={h} lf={lf} shrink={int(shrink)} "
      f"interp={int(interp)} a={sd0.get_center_shrink_factor():.3f} acc={acc:.4f} "
      f"tau_mean={np.mean(tau):.1f} tau_max={np.max(tau):.1f} time={dt:.0f}s", flush=True)

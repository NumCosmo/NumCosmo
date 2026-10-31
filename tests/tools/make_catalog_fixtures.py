#!/usr/bin/env python
"""Write small catalog files covering every state the catalog readers must accept.

The files go to data/truth_tables/catalog/ and are read by tests/c/ncm/fit/
test_ncm_mset_catalog_fixtures.c, which asserts the state each one carries.

Four files are real runs of the current library on a two-parameter NcmModelMVND
with an NcmDataGaussCovMVND likelihood, one per catalog producer, with the run
type and sampler strings that producer writes today:

  esmcmc.fits        NcmFitESMCMC, 6 walkers, 3 ensembles: RTYPE "Ensemble
                     Sampler MCMC", INITSMP, SAMPLER, SAMPOPT, MARKID
  mc.fits            NcmFitMC, RTYPE "from-model", 10 rows
  mcbs.fits          NcmFitMCBS, RTYPE "bootstrap-mean", 5 rows
  mcmc.fits          NcmFitMCMC, RTYPE "Multivariate Gaussian Sampler", 10 rows

Two are written through NcmMSetCatalog directly:

  weighted.fits      weighted rows, one chain, RTYPE "Fixture"
  functions.fits     a `functions' key in HDU0 holding one NcmMSetFuncList

The rest are header edits of esmcmc.fits, each a state the library has written
or must tolerate:

  esmcmc_no_rng.fits        no RNG_ALGO, RNG_SEED, RNG_INIS, RNG_STAT
  esmcmc_no_rng_algo.fits   RNG_SEED, RNG_INIS, RNG_STAT without RNG_ALGO
                            (no library version writes this; the readers
                            fall back to the default algorithm)
  esmcmc_no_m2lnp_id.fits   no M2LNP_ID
  esmcmc_first_id.fits      FIRST_ID 6 (a multiple of the 6 chains), MARKID shifted by
                            the same amount
  esmcmc_no_markid.fits     no MARKID
  esmcmc_header_only.fits   the table header with no rows

The three HDU0 containers (vardict, bare object, `.mset' sidecar) are the
old_reparam_* files of data/truth_tables/bbn/.
"""

import os
import shutil
import tempfile

import numpy as np
from astropy.io import fits

from numcosmo_py import Nc, Ncm

Ncm.cfg_init()

OUT = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "..",
    "..",
    "data",
    "truth_tables",
    "catalog",
)

SEED = 20261031
DIM = 2
NWALKERS = 6
NENSEMBLES = 3
NROWS_MC = 10
NROWS_MCBS = 5
NROWS_MCMC = 10
NROWS_WEIGHTED = 10
FIRST_ID = 6

RNG_KEYS = ("RNG_ALGO", "RNG_SEED", "RNG_INIS", "RNG_STAT")


def build_fit(rng):
    """A fit of a two-parameter MVND model to MVND data, as the catalog tests use."""
    data_mvnd = Ncm.DataGaussCovMVND.new_full(DIM, 1.0e-1, 1.0, 50.0, -1.0, 1.0, rng)
    model = Ncm.ModelMVND.new(DIM)
    mset = Ncm.MSet.new_array([model])
    mset.param_set_all_ftype(Ncm.ParamType.FREE)
    mset.prepare_fparam_map()

    dset = Ncm.Dataset.new_array([data_mvnd])
    lh = Ncm.Likelihood.new(dset)
    fit = Ncm.Fit.factory(
        Ncm.FitType.NLOPT,
        "ln-neldermead",
        lh,
        mset,
        Ncm.FitGradType.NUMDIFF_FORWARD,
    )

    return fit, data_mvnd


def write_esmcmc(path):
    """Three ensembles of the stretch walker from a Gaussian initial sampler."""
    rng = Ncm.RNG.seeded_new(None, SEED)
    fit, _ = build_fit(rng)
    mset = fit.peek_mset()

    init_sampler = Ncm.MSetTransKernGauss.new(0)
    init_sampler.set_mset(mset)
    init_sampler.set_prior_from_mset()
    init_sampler.set_cov_from_rescale(1.0e-1)

    walker = Ncm.FitESMCMCWalkerStretch.new(NWALKERS, mset.fparams_len())
    esmcmc = Ncm.FitESMCMC.new(fit, NWALKERS, init_sampler, walker, Ncm.FitRunMsgs.NONE)
    esmcmc.set_rng(rng)
    esmcmc.set_data_file(path)
    esmcmc.start_run()
    esmcmc.run(NENSEMBLES)
    esmcmc.end_run()


def write_mc(path):
    """Monte Carlo resampling from the fiducial model."""
    rng = Ncm.RNG.seeded_new(None, SEED)
    fit, _ = build_fit(rng)

    mc = Ncm.FitMC.new(fit, Ncm.FitMCResampleType.FROM_MODEL, Ncm.FitRunMsgs.NONE)
    mc.set_rng(rng)
    mc.set_data_file(path)
    mc.start_run()
    mc.run(NROWS_MC)
    mc.end_run()


def write_mcbs(path):
    """Bootstrap means; the per-realization bootstrap catalogs are written next
    to the file and not kept."""
    rng = Ncm.RNG.seeded_new(None, SEED)
    fit, _ = build_fit(rng)

    mcbs = Ncm.FitMCBS.new(fit)
    mcbs.set_rng(rng)
    mcbs.set_filename(path)
    mcbs.run(
        None,
        0,
        NROWS_MCBS,
        10,
        Ncm.FitMCResampleType.BOOTSTRAP_NOMIX,
        Ncm.FitRunMsgs.NONE,
        False,
    )
    # ncm_fit_mcbs_set_filename detaches the file again through
    # ncm_mset_catalog_reset, so the rows reach the file only by setting it after
    # the run (ncm_mset_catalog_set_file syncs the rows it holds).
    mcat = mcbs.get_catalog()
    mcat.set_file(path)
    mcat.sync(True)


def write_mcmc(path):
    """Metropolis-Hastings with a Gaussian transition kernel."""
    rng = Ncm.RNG.seeded_new(None, SEED)
    fit, data_mvnd = build_fit(rng)

    tkern = Ncm.MSetTransKernGauss.new(0)
    mcmc = Ncm.FitMCMC.new(fit, tkern, Ncm.FitRunMsgs.NONE)
    tkern.set_cov(data_mvnd.peek_cov())
    mcmc.set_rng(rng)
    mcmc.set_data_file(path)
    mcmc.start_run()
    mcmc.run(NROWS_MCMC)
    mcmc.end_run()


def build_catalog(nadd_vals, names, symbols, weighted):
    """A catalog on the MVND model set, built the way the old API builds one."""
    rng = Ncm.RNG.seeded_new(None, SEED)
    fit, _ = build_fit(rng)
    mcat = Ncm.MSetCatalog(
        mset=fit.peek_mset(),
        nadd_vals=nadd_vals,
        nadd_val_names=names,
        nadd_val_symbols=symbols,
        nchains=1,
        m2lnp_var=0,
        weighted=weighted,
    )
    mcat.set_run_type("Fixture")

    return mcat


def write_weighted(path):
    """Weighted rows added as vectors: m2lnL, the two parameters and the weight."""
    mcat = build_catalog(1, ["m2lnL"], ["-2\\ln(L)"], True)
    mcat.set_file(path)

    rng = np.random.default_rng(SEED)

    for i in range(NROWS_WEIGHTED):
        x = rng.normal(size=DIM)
        row = [float(np.dot(x, x))] + [float(v) for v in x] + [1.0 / (i + 1)]
        mcat.add_from_vector(Ncm.Vector.new_array(row))

    mcat.sync(True)


def write_functions(path):
    """A `functions' key in HDU0: one registered NcmMSetFuncList object."""
    mcat = build_catalog(2, ["m2lnL", "H0"], ["-2\\ln(L)", "H_0"], False)
    functions = Ncm.ObjArray.new()
    functions.add(Ncm.MSetFuncList.new("NcHICosmo:H0", None))
    mcat.set_functions_array(functions)
    mcat.set_file(path)

    rng = np.random.default_rng(SEED)

    for _ in range(NROWS_WEIGHTED):
        x = rng.normal(size=DIM)
        row = [float(np.dot(x, x)), 70.0] + [float(v) for v in x]
        mcat.add_from_vector(Ncm.Vector.new_array(row))

    mcat.sync(True)


def edit_header(src, dst, remove=(), update=None, drop_rows=False):
    """Copy @src to @dst without the keys in @remove and with the pairs in @update."""
    with fits.open(src) as hdul:
        table = hdul[1]

        for key in remove:
            del table.header[key]

        for key, value in (update or {}).items():
            table.header[key] = value

        if drop_rows:
            table = fits.BinTableHDU(data=table.data[:0], header=table.header)

        fits.HDUList([hdul[0].copy(), table]).writeto(dst, overwrite=True)


def main():
    """Write every fixture."""
    os.makedirs(OUT, exist_ok=True)

    with tempfile.TemporaryDirectory() as tmp:
        # MCBS writes its bootstrap catalogs next to its file.
        mcbs_tmp = os.path.join(tmp, "mcbs.fits")
        write_mcbs(mcbs_tmp)
        shutil.copyfile(mcbs_tmp, os.path.join(OUT, "mcbs.fits"))

    esmcmc = os.path.join(OUT, "esmcmc.fits")

    for path in (
        "esmcmc.fits",
        "mc.fits",
        "mcmc.fits",
        "weighted.fits",
        "functions.fits",
    ):
        full = os.path.join(OUT, path)

        if os.path.exists(full):
            os.remove(full)

    write_esmcmc(esmcmc)
    write_mc(os.path.join(OUT, "mc.fits"))
    write_mcmc(os.path.join(OUT, "mcmc.fits"))
    write_weighted(os.path.join(OUT, "weighted.fits"))
    write_functions(os.path.join(OUT, "functions.fits"))

    with fits.open(esmcmc) as hdul:
        markid = int(hdul[1].header["MARKID"])

    edit_header(esmcmc, os.path.join(OUT, "esmcmc_no_rng.fits"), remove=RNG_KEYS)
    edit_header(
        esmcmc, os.path.join(OUT, "esmcmc_no_rng_algo.fits"), remove=("RNG_ALGO",)
    )
    edit_header(
        esmcmc, os.path.join(OUT, "esmcmc_no_m2lnp_id.fits"), remove=("M2LNP_ID",)
    )
    edit_header(
        esmcmc,
        os.path.join(OUT, "esmcmc_first_id.fits"),
        update={"FIRST_ID": FIRST_ID, "MARKID": markid + FIRST_ID},
    )
    edit_header(esmcmc, os.path.join(OUT, "esmcmc_no_markid.fits"), remove=("MARKID",))
    edit_header(esmcmc, os.path.join(OUT, "esmcmc_header_only.fits"), drop_rows=True)

    for name in sorted(os.listdir(OUT)):
        print(f"{name:28s} {os.path.getsize(os.path.join(OUT, name)):6d} bytes")


if __name__ == "__main__":
    main()

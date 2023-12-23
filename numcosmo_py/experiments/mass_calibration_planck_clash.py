#!/usr/bin/env python
#
# mass_calibration_planck_clash.py
#
# Wed May 3 10:21:00 2015
# Copyright  2015  Mariana Penna-Lima
# <pennalima@gmail.com>
#
# mass_calibration_planck_clash.py
# Copyright (C) 2015 Mariana Penna-Lima <pennalima@gmail.com>
#
# numcosmo is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# numcosmo is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program.  If not, see <http://www.gnu.org/licenses/>.

"""Old Planck-CLASH mass calibration experiment."""

from typing import Tuple

from numcosmo_py import Ncm, Nc

Ncm.cfg_init()


def create_cosmology() -> Nc.HICosmo:
    """Create a cosmology."""

    cosmo = Nc.HICosmoDEXcdm()
    prim = Nc.HIPrimPowerLaw.new()
    reion = Nc.HIReionCamb.new()

    cosmo.add_submodel(reion)
    cosmo.add_submodel(prim)

    # Cosmological parameters
    cosmo.props.H0 = 70.0
    cosmo.props.Omegab = 0.049
    cosmo.props.Omegac = 0.251
    cosmo.props.Omegax = 0.7
    cosmo.props.Tgamma0 = 2.72
    prim.props.n_SA = 0.967
    prim.props.ln10e10ASA = 3.064  # cosmo.props.sigma8  = 0.816
    cosmo.props.w = -1.0
    cosmo.omega_x2omega_k()
    cosmo.param_set_by_name("Omegak", 0.0)

    return cosmo


def create_cluster_abundance(
    Tinker_lin_interp=True,
) -> Tuple[Nc.ClusterAbundance, Nc.ClusterRedshiftNodist, Nc.ClusterMassPlCL]:
    """Create a cluster abundance."""
    # Cosmological distances
    dist = Nc.Distance.new(4.0)
    # Transfer function
    tf = Nc.TransferFuncEH()
    # Linear matter power spectrum
    ps_ml = Nc.PowspecMLTransfer.new(tf)
    psf = Ncm.PowspecFilter.new(ps_ml, Ncm.PowspecFilterType.TOPHAT)
    mulf = Nc.MultiplicityFuncTinker.new()
    mulf.set_linear_interp(Tinker_lin_interp)
    mulf.set_mdef(Nc.MultiplicityFuncMassDef.MEAN)
    mulf.set_Delta(500.0)
    mf = Nc.HaloMassFunction.new(dist, psf, mulf)
    cad = Nc.ClusterAbundance.new(mf, None)
    clusterz = Nc.ClusterRedshiftNodist(z_min=0.0, z_max=2.0)
    # Planck - CLASH cluster mass distribution, pivot mass M0
    clusterm = Nc.ClusterMassPlCL(M0=5.7e14)

    for i in range(7):
        clusterm.param_set_ftype(i, Ncm.ParamType.FREE)

    clusterm.param_set_by_name("Asz", 1.00)
    clusterm.param_set_by_name("Bsz", 0.25)
    clusterm.param_set_by_name("sigma_sz", 0.12)
    clusterm.param_set_by_name("Al", 1.00)
    clusterm.param_set_by_name("Bl", 0.0)
    clusterm.param_set_by_name("sigma_l", 0.27)
    clusterm.param_set_by_name("cor", 0.0)

    return cad, clusterz, clusterm


def run_example():
    """Run the example."""

    Tinker_lin_interp = True

    NT = 3  # Number of threads
    NClusters = 21  # Number of clusters
    NWalkers = 100  # Number of walkers / chains

    cosmo = create_cosmology()
    cad, clusterz, clusterm = create_cluster_abundance(Tinker_lin_interp)
    cad.prepare(cosmo, clusterz, clusterm)

    cpc = Nc.ClusterPseudoCounts.new(NClusters)

    mset = Ncm.MSet.empty_new()
    mset.set(cosmo)
    mset.set(clusterm)
    mset.set(clusterz)
    mset.set(cpc)

    for i in range(4):
        cpc.param_set_ftype(i, Ncm.ParamType.FREE)

    cpc.param_set_by_name("lnMCut", 33.0)
    cpc.param_set_by_name("sigma_Mcut", 0.10)
    cpc.param_set_by_name("zmin", 0.188)
    cpc.param_set_by_name("Deltaz", 0.70214)

    plclash = Nc.DataClusterPseudoCounts.new_from_file(
        "nc_data_cluster_planck_clash.obj"
    )
    plclash.set_cad(cad)

    dset = Ncm.Dataset.new()

    dset.append_data(plclash)
    lh = Ncm.Likelihood(dataset=dset)

    # Gaussian prior on the lensing bias, b_l = 0 \pm 0.08
    lh.priors_add_gauss_param(clusterm.id(), 4, 0.0, 0.08)

    algorithm = "ln-neldermead"
    fit = Ncm.Fit.factory(
        Ncm.FitType.NLOPT, algorithm, lh, mset, Ncm.FitGradType.NUMDIFF_CENTRAL
    )

    Ncm.func_eval_set_max_threads(NT)
    Ncm.func_eval_log_pool_stats()

    init_sampler = Ncm.MSetTransKernGauss.new(0)
    stretch = Ncm.FitESMCMCWalkerStretch.new(NWalkers, mset.fparams_len())
    esmcmc = Ncm.FitESMCMC.new(
        fit, NWalkers, init_sampler, stretch, Ncm.FitRunMsgs.FULL
    )

    init_sampler.set_mset(mset)
    init_sampler.set_prior_from_mset()

    esmcmc.set_nthreads(NT)
    init_sampler.set_cov_from_scale()

    esmcmc.set_data_file("test.fits")

    esmcmc.start_run()
    esmcmc.run(1000)  # Number of points to be computed in each chain
    esmcmc.end_run()

    esmcmc.mean_covar()
    fit.log_covar()


if __name__ == "__main__":
    run_example()

#!/usr/bin/env python
#
# example_wl_likelihood.py
#
# Mon May 22 16:00:00 2023
# Copyright  2023  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# example_wl_likelihood.py
# Copyright (C) 2023 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

"""Example using weak lensing likelihood."""

import numpy as np
from numcosmo_py import Nc, Ncm

#
#  Initializing the library objects, this must be called before
#  any other library function.
#
Ncm.cfg_init()


def test_wl_likelihood() -> None:
    """Example using weak lensing likelihood."""
    seed = 1235
    H0 = 70.0
    Omegab = 0.045
    Omegac = 0.255
    Omegak = 0.0
    w = -1.0
    cluster_ra, cluster_dec = 12.34, -55.123
    cluster_z = 0.2
    cluster_c = 4.0
    cluster_mass = 1.0e14
    n_galaxies = 10000
    sigma0 = 0.03
    galaxies_ang_dist = 0.15
    galaxy_shape_sigma = 3.0e-1
    galaxy_shape_std_noise = 1.0e-1
    min_r, max_r = 0.3, 3.0

    dist = Nc.Distance.new(5.0)

    cosmo = Nc.HICosmoDEXcdm()
    cosmo.omega_x2omega_k()

    halo_mass_summary = Nc.HaloCMParam.new(Nc.HaloMassSummaryMassDef.CRITICAL, 200.0)
    density_profile = Nc.HaloDensityProfileNFW.new(halo_mass_summary)
    surface_mass_density = Nc.WLSurfaceMassDensity.new(dist)
    halo_position = Nc.HaloPosition.new(dist)

    halo_mass_summary["cDelta"] = cluster_c
    halo_mass_summary["log10MDelta"] = np.log10(cluster_mass)

    galaxies_ra0 = cluster_ra - galaxies_ang_dist
    galaxies_ra1 = cluster_ra + galaxies_ang_dist
    galaxies_dec0 = cluster_dec - galaxies_ang_dist
    galaxies_dec1 = cluster_dec + galaxies_ang_dist

    position_factor = Nc.GalaxyPositionFactorFlat.new(
        galaxies_ra0, galaxies_ra1, galaxies_dec0, galaxies_dec1
    )
    redshift_factor = Nc.GalaxyRedshiftFactorComposed.new(0.01, 6.0)
    redshift_pop = Nc.GalaxyRedshiftPopLSSTSRD.new_y1_source()
    redshift_obs_sel = Nc.GalaxyRedshiftObsGauss.new()

    ellip_conv = Nc.GalaxyWLObsEllipConv.TRACE
    ellip_frame = Nc.WLEllipticityFrame.CARTESIAN

    shape_factor = Nc.GalaxyShapeFactorVarAdd.new(ellip_conv)
    shape_pop = Nc.GalaxyShapePopGauss.new()

    cosmo["H0"] = H0
    cosmo["Omegab"] = Omegab
    cosmo["Omegac"] = Omegac
    cosmo["Omegak"] = Omegak
    cosmo["w"] = w
    halo_position["ra"] = cluster_ra
    halo_position["dec"] = cluster_dec
    halo_position["z"] = cluster_z
    shape_pop["sigma"] = galaxy_shape_sigma

    dist.prepare(cosmo)
    halo_position.prepare(cosmo)
    surface_mass_density.prepare(cosmo)

    halo_position.param_set_desc("ra", {"fit": True})
    halo_position.param_set_desc("dec", {"fit": True})
    halo_mass_summary.param_set_desc("cDelta", {"fit": False})
    halo_mass_summary.param_set_desc("log10MDelta", {"fit": True})

    mset = Ncm.MSet.new_array(
        [
            cosmo,
            density_profile,
            surface_mass_density,
            halo_position,
            redshift_pop,
            redshift_obs_sel,
            shape_pop,
        ]
    )
    mset.prepare_fparam_map()

    rng = Ncm.RNG.seeded_new("mt19937", seed)

    pos_data = Nc.GalaxyPositionFactorData.new(position_factor, mset)
    z_data = Nc.GalaxyRedshiftFactorData.new(redshift_factor, mset)
    s_data = Nc.GalaxyShapeFactorData.new(shape_factor, mset, pos_data, z_data)
    cols = Nc.GalaxyShapeFactorData.required_columns(s_data)

    obs = Nc.GalaxyWLObs.new(ellip_conv, ellip_frame, n_galaxies, cols)

    for i in range(n_galaxies):
        for c in cols:
            obs.set(c, i, 0.0)

        c1 = rng.gaussian_gen(0.0, 0.01)
        c2 = rng.gaussian_gen(0.0, 0.01)
        m = np.exp(rng.gaussian_gen(0.0, 0.08))

        obs.set(Nc.GALAXY_REDSHIFT_OBS_GAUSS_COL_SIGMA0, i, sigma0)
        obs.set("std_noise", i, galaxy_shape_std_noise)
        obs.set("c1", i, c1)
        obs.set("c2", i, c2)
        obs.set("m", i, m)

    cluster = Nc.DataClusterWLFactor.new(position_factor, redshift_factor, shape_factor)
    cluster.set_cut(min_r / cosmo.h(), max_r / cosmo.h())
    cluster.set_obs(obs)
    cluster.resample(mset, rng)

    dset = Ncm.Dataset()
    dset.append_data(cluster)

    lh = Ncm.Likelihood(dataset=dset)
    fit = Ncm.Fit.factory(
        Ncm.FitType.NLOPT, "ln-neldermead", lh, mset, Ncm.FitGradType.NUMDIFF_FORWARD
    )

    fit.run(Ncm.FitRunMsgs.SIMPLE)
    fit.log_info()

    # fit.obs_fisher()
    # fit.log_covar()

    init_sampler = Ncm.MSetTransKernGauss.new(0)
    init_sampler.set_mset(mset)
    init_sampler.set_prior_from_mset()
    # init_sampler.set_cov(fit.get_covar())
    init_sampler.set_cov_from_rescale(1.0)

    nwalkers = 300
    nthreads = 6
    Ncm.cfg_set_openmp_nthreads(nthreads)
    apes = Ncm.FitESMCMCWalkerAPES.new(nwalkers, mset.fparams_len())
    apes.set_use_threads(True)
    esmcmc = Ncm.FitESMCMC.new(fit, nwalkers, init_sampler, apes, Ncm.FitRunMsgs.SIMPLE)
    esmcmc.set_nthreads(nthreads)

    esmcmc.set_data_file("example_wl_likelihood.fits")

    esmcmc.start_run()
    esmcmc.run(100)
    esmcmc.end_run()

    esmcmc.mean_covar()
    fit.log_covar()


if __name__ == "__main__":
    test_wl_likelihood()

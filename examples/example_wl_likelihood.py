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
    use_photoz = False
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
    sigma_z = 0.03
    galaxy_true_z_min = 0.2
    galaxy_true_z_max = 1.2
    galaxies_ang_dist = 0.15
    galaxy_shape_e_rms = 3.0e-1
    galaxy_shape_e_sigma = 1.0e-1
    min_r, max_r = 0.3, 3.0

    dist = Nc.Distance.new(5.0)

    cosmo = Nc.HICosmoDEXcdm()
    cosmo.omega_x2omega_k()

    density_profile = Nc.HaloDensityProfileNFW.new(
        Nc.HaloDensityProfileMassDef.CRITICAL, 200.0
    )
    surface_mass_density = Nc.WLSurfaceMassDensity.new(dist)
    halo_position = Nc.HaloPosition.new(dist)

    density_profile["cDelta"] = cluster_c
    density_profile["log10MDelta"] = np.log10(cluster_mass)

    galaxies_ra0 = cluster_ra - galaxies_ang_dist
    galaxies_ra1 = cluster_ra + galaxies_ang_dist
    galaxies_dec0 = cluster_dec - galaxies_ang_dist
    galaxies_dec1 = cluster_dec + galaxies_ang_dist

    galaxy_position = Nc.GalaxySDPositionFlat.new(
        galaxies_ra0, galaxies_ra1, galaxies_dec0, galaxies_dec1
    )
    galaxy_redshift_true = Nc.GalaxySDTrueRedshiftLSSTSRD.new(
        galaxy_true_z_min, galaxy_true_z_max
    )
    if use_photoz:
        galaxy_redshift = Nc.GalaxySDObsRedshiftSpec.new(galaxy_redshift_true)
    else:
        galaxy_redshift = Nc.GalaxySDObsRedshiftGauss.new(galaxy_redshift_true)

    galaxy_shape = Nc.GalaxySDShapeGauss.new()

    cosmo["H0"] = H0
    cosmo["Omegab"] = Omegab
    cosmo["Omegac"] = Omegac
    cosmo["Omegak"] = Omegak
    cosmo["w"] = w
    halo_position["ra"] = cluster_ra
    halo_position["dec"] = cluster_dec
    halo_position["z"] = cluster_z
    galaxy_shape["e-rms"] = galaxy_shape_e_rms

    dist.prepare(cosmo)
    halo_position.prepare(cosmo)
    surface_mass_density.prepare(cosmo)

    cluster = Nc.DataClusterWL.new()
    cluster.set_cut(min_r, max_r)

    halo_position.param_set_desc("ra", {"fit": True})
    halo_position.param_set_desc("dec", {"fit": True})
    density_profile.param_set_desc("cDelta", {"fit": False})
    density_profile.param_set_desc("log10MDelta", {"fit": True})

    mset = Ncm.MSet.new_array(
        [
            cosmo,
            density_profile,
            surface_mass_density,
            halo_position,
            galaxy_redshift,
            galaxy_position,
            galaxy_shape,
        ]
    )

    rng = Ncm.RNG.seeded_new("mt19937", seed)

    z_data = Nc.GalaxySDObsRedshiftData.new(galaxy_redshift)
    p_data = Nc.GalaxySDPositionData.new(galaxy_position, z_data)
    s_data = Nc.GalaxySDShapeData.new(galaxy_shape, p_data)

    obs = Nc.GalaxyWLObs.new(
        Nc.GalaxyWLObsCoord.EUCLIDEAN, n_galaxies, list(s_data.required_columns())
    )

    for i in range(n_galaxies):
        if use_photoz:
            galaxy_redshift.gen(mset, z_data, rng)
        else:
            galaxy_redshift.gen(mset, z_data, sigma_z, rng)  # type: ignore

        galaxy_position.gen(mset, p_data, rng)
        galaxy_shape.gen(
            mset,
            s_data,
            galaxy_shape_e_sigma,
            galaxy_shape_e_sigma,
            Nc.GalaxyWLObsCoord.EUCLIDEAN,
            rng,
        )
        s_data.write_row(obs, i)
    cluster.set_obs(obs)

    dset = Ncm.Dataset()
    dset.append_data(cluster)

    lh = Ncm.Likelihood(dataset=dset)
    fit = Ncm.Fit.factory(
        Ncm.FitType.NLOPT, "ln-neldermead", lh, mset, Ncm.FitGradType.NUMDIFF_FORWARD
    )

    fit.run(Ncm.FitRunMsgs.SIMPLE)
    fit.log_info()

    fit.obs_fisher()
    fit.log_covar()

    init_sampler = Ncm.MSetTransKernGauss.new(0)
    init_sampler.set_mset(mset)
    init_sampler.set_prior_from_mset()
    init_sampler.set_cov(fit.get_covar())

    nwalkers = 300
    apes = Ncm.FitESMCMCWalkerAPES.new(nwalkers, mset.fparams_len())
    esmcmc = Ncm.FitESMCMC.new(fit, nwalkers, init_sampler, apes, Ncm.FitRunMsgs.SIMPLE)

    esmcmc.set_data_file("example_wl_likelihood.fits")

    esmcmc.start_run()
    esmcmc.run(100)
    esmcmc.end_run()

    esmcmc.mean_covar()
    fit.log_covar()


if __name__ == "__main__":
    test_wl_likelihood()

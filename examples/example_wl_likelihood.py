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

from numcosmo_py import Nc, Ncm

#
#  Initializing the library objects, this must be called before
#  any other library function.
#
Ncm.cfg_init()


def test_wl_likelihood() -> None:
    """Example using weak lensing likelihood."""

    dist = Nc.Distance.new(5.0)

    cosmo = Nc.HICosmoDEXcdm()
    cosmo.omega_x2omega_k()

    cosmo.props.H0 = 70.0
    cosmo.props.Omegab = 0.045
    cosmo.props.Omegac = 0.255
    cosmo.props.Omegax = 0.7
    cosmo.param_set_by_name("Omegak", 0.0)

    nfw = Nc.HaloDensityProfile.new_from_name(
        "NcHaloDensityProfileNFW{'Delta':<200.0>}"
    )
    # 4 as Douglas. In LCDM c = 5 corresponds to cluster masses.
    # (see Lokas and G. Mamon, astro-ph/0002395)
    nfw.param_set_by_name("cDelta", 5.0)
    nfw.param_set_by_name("MDelta", 1.0e15)

    smd = Nc.WLSurfaceMassDensity.new(dist)
    rs = Nc.ReducedShearClusterMass.new()

    mset = Ncm.MSet.new_array([cosmo, nfw, smd, rs])
    mset.param_set_ftype(
        Nc.HaloDensityProfile.id(),
        Nc.HaloDensityProfileSParams.LOG10M_DELTA,
        Ncm.ParamType.FREE,
    )
    mset.param_set_ftype(
        Nc.HaloDensityProfile.id(),
        Nc.HaloDensityProfileSParams.C_DELTA,
        Ncm.ParamType.FREE,
    )

    d1 = Nc.DataReducedShearClusterMass.new(dist)
    d1.load_hdf5("cat07_sim_leftra.hdf5", ord("i"), 0.3, 0.1, 0.1)

    ser = Ncm.Serialize.new(Ncm.SerializeOpt.NONE)
    ser.to_file(d1, "test.obj")

    dset = Ncm.Dataset()
    dset.append_data(d1)

    lh = Ncm.Likelihood(dataset=dset)
    fit = Ncm.Fit.factory(
        Ncm.FitType.NLOPT, "ln-neldermead", lh, mset, Ncm.FitGradType.NUMDIFF_FORWARD
    )

    fit.run(Ncm.FitRunMsgs.SIMPLE)

    fit.log_info()

    fit.numdiff_m2lnL_covar()

    fit.log_covar()

    init_sampler = Ncm.MSetTransKernGauss.new(0)
    init_sampler.set_mset(mset)
    init_sampler.set_prior_from_mset()
    init_sampler.set_cov_from_rescale(1.0)

    nwalkers = 50
    stretch = Ncm.FitESMCMCWalkerStretch.new(nwalkers, mset.fparams_len())

    stretch.set_scale(2.5)
    stretch.set_box_mset(mset)

    esmcmc = Ncm.FitESMCMC.new(
        fit, nwalkers, init_sampler, stretch, Ncm.FitRunMsgs.SIMPLE
    )

    esmcmc.set_auto_trim(True)
    esmcmc.set_auto_trim_div(100)
    esmcmc.set_max_runs_time(2.0 * 60.0)

    esmcmc.set_data_file("example_iii.fits")

    esmcmc.start_run()
    esmcmc.run_lre(10, 1.0e-3)
    esmcmc.end_run()

    esmcmc.mean_covar()
    fit.log_covar()


if __name__ == "__main__":
    test_wl_likelihood()

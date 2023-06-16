#!/usr/bin/env python
#
# example_rg_snia_bao.py
#
# Mon May 22 16:00:00 2023
# Copyright  2023  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# example_rg_snia_bao.py
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

"""Example computing confidence regions."""

import matplotlib.pyplot as plt

from numcosmo_py import Nc, Ncm

#
#  Initializing the library objects, this must be called before
#  any other library function.
#
Ncm.cfg_init()


def test_rg_snia_bao() -> None:
    """Example computing confidence regions."""

    #
    #  New homogeneous and isotropic cosmological model NcHICosmoDEXcdm
    #
    cosmo = Nc.HICosmoDEXcdm()

    #
    #  Setting values for the cosmological model, those not set stay in the
    #  default values. Remeber to use the _orig_ version to set the original
    #  parameters in case when a reparametrization is used.
    #

    #
    # OO-like
    #
    cosmo.props.H0 = 70.0
    cosmo.props.Omegab = 0.05
    cosmo.props.Omegac = 0.25
    cosmo.props.Omegax = 0.70
    cosmo.props.Tgamma0 = 2.72
    cosmo.props.w = -1.0

    #
    #  A new Modelset with cosmo as the HICosmo model to be used.
    #
    mset = Ncm.MSet()
    mset.set(cosmo)

    #
    #  Setting parameters Omega_c, Omega_x and w to be fitted (and change
    #  parameter Omega_x -> Omega_k).
    #
    # cosmo.de_omega_x2omega_k ()
    cosmo.props.Omegac_fit = True
    cosmo.props.Omegax_fit = True
    cosmo.props.w_fit = True

    #
    #  A new Distance object optimized to redshift 2.5.
    #
    dist = Nc.Distance(zf=2.5)

    #
    #  A new Data object from distance modulus catalogs. Type Ia supernovae (snia).
    #  A new Data object from BAO catalogs.
    #
    snia = Nc.DataDistMu.new_from_id(dist, Nc.DataSNIAId.SIMPLE_UNION2_1)
    bao1 = Nc.data_bao_create(dist, Nc.DataBaoId.RDV_BEUTLER2011)
    bao2 = Nc.data_bao_create(dist, Nc.DataBaoId.EMPIRICAL_FIT_ROSS2015)
    bao3 = Nc.data_bao_create(dist, Nc.DataBaoId.RDV_BOSS_QSO_ATA2017)

    #
    #  A new Dataset with snia and bao set.
    #
    dset = Ncm.Dataset()
    dset.append_data(snia)
    dset.append_data(bao1)
    dset.append_data(bao2)
    dset.append_data(bao3)

    #
    #  Creating a Likelihood from the Dataset.
    #
    lh = Ncm.Likelihood(dataset=dset)

    #
    #  Creating a Fit object of type NLOPT using the fitting algorithm ln-neldermead to
    #  fit the Modelset mset using the Likelihood lh and using a numerical differentiation
    #  algorithm (NUMDIFF_FORWARD) to obtain the gradient (if needed).
    #
    fit = Ncm.Fit.new(
        Ncm.FitType.NLOPT, "ln-neldermead", lh, mset, Ncm.FitGradType.NUMDIFF_FORWARD
    )

    #
    #  Running the fitter printing messages.
    #
    fit.run(Ncm.FitRunMsgs.SIMPLE)

    #
    #  Printing fitting informations.
    #
    fit.log_info()

    #
    #  Calculating the parameters covariance using numerical differentiation.
    #
    fit.numdiff_m2lnL_covar()

    #
    #  Printing the covariance matrix.
    #
    fit.log_covar()

    #
    #  Creating a new Likelihood ratio test object.
    #  First we create two PIndex indicating which parameter
    #    we are going to study.
    #
    p1 = Ncm.MSetPIndex.new(cosmo.id(), Nc.HICosmoDESParams.OMEGA_C)
    p2 = Ncm.MSetPIndex.new(cosmo.id(), Nc.HICosmoDEXCDMSParams.W)

    lhr2d = Ncm.LHRatio2d.new(fit, p1, p2, 1.0e-3)

    #
    #  Calculating the confidence region using the Likelihood ratio test.
    #  Also calculate using the Fisher matrix approach.
    #
    cr_rg = lhr2d.conf_region(0.6827, 300.0, Ncm.FitRunMsgs.SIMPLE)
    fisher_rg = lhr2d.fisher_border(0.6827, 300.0, Ncm.FitRunMsgs.SIMPLE)

    cr_p1array = cr_rg.p1.dup_array()
    cr_p2array = cr_rg.p2.dup_array()

    fisher_p1array = fisher_rg.p1.dup_array()
    fisher_p2array = fisher_rg.p2.dup_array()

    #
    #  Ploting the confidence regions obtained from both methods.
    #

    plt.title(f"Confidence regions ({cr_rg.clevel * 100.0:.2f})")
    plt.plot(cr_p1array, cr_p2array, "r", label="Likelihood Ratio")
    plt.plot(fisher_p1array, fisher_p2array, "b-", label="Fisher Matrix")

    plt.xlabel(r"$\Omega_c$")
    plt.ylabel(r"$w$")

    plt.legend(loc=4)

    plt.savefig("snia_bao_rg_omegac_w.svg")


if __name__ == "__main__":
    test_rg_snia_bao()

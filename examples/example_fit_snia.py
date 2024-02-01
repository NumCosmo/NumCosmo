#!/usr/bin/env python
#
# example_fit_snia.py
#
# Mon May 22 16:00:00 2023
# Copyright  2023  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# example_fit_snia.py
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

"""Example fitting SNIa data."""

from numcosmo_py import Nc, Ncm

#
#  Initializing the library objects, this must be called before
#  any other library function.
#
Ncm.cfg_init()


def test_fit_snia() -> None:
    """Example fitting SNIa data."""

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
    #  Creating a new Modelset and set cosmo as the HICosmo model to be used.
    #
    mset = Ncm.MSet()
    mset.set(cosmo)

    #
    #  Setting parameters Omega_c and w to be fitted.
    #
    cosmo.props.Omegac_fit = True
    cosmo.props.w_fit = True

    #
    #  Creating a new Distance object optimized to redshift 2.
    #
    dist = Nc.Distance(zf=2.0)

    #
    #  Creating a new Data object from distance modulus catalogs.
    #
    snia = Nc.DataDistMu.new_from_id(dist, Nc.DataSNIAId.SIMPLE_UNION2_1)

    #
    #  Creating a new Dataset and add snia to it.
    #
    dset = Ncm.Dataset()
    dset.append_data(snia)

    #
    #  Creating a Likelihood from the Dataset.
    #
    lh = Ncm.Likelihood(dataset=dset)

    #
    # Creating a Fit object of type NLOPT using the fitting algorithm ln-neldermead to
    # fit the Modelset mset using the Likelihood lh and using a numerical
    # differentiation algorithm (NUMDIFF_FORWARD) to obtain the gradient (if needed).
    #
    fit = Ncm.Fit.factory(
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


if __name__ == "__main__":
    test_fit_snia()

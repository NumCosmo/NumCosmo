#!/usr/bin/env python
#
# example_mcmc.py
#
# Mon May 22 16:00:00 2023
# Copyright  2023  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# example_mcmc.py
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

"""Example of running a MCMC analysis on a SLine model."""

import sys
import os.path
from typing import Union

from py_sline_model import PySLineModel
from py_sline_data import PySLineData
from py_sline_gauss import PySLineGauss

from numcosmo_py import Ncm


#
#  Initializing the library objects, this must be called before
#  any other library function.
#
Ncm.cfg_init()


def run_mcmc() -> None:
    """Example of running a MCMC analysis on a SLine model."""

    #
    # Instantiating a new SLine model object and setting
    # some values for its parameters.
    #
    slm = PySLineModel()
    slm.props.alpha = 0.9
    slm.props.a = 0.1

    #
    # New Model set object including slm with parameters
    # set as free.
    #
    mset = Ncm.MSet.empty_new()
    mset.set(slm)
    mset.param_set_all_ftype(Ncm.ParamType.FREE)
    mset.prepare_fparam_map()

    #
    # Creating a new Serialization object, and load
    # the data file.
    #
    sld: Union[PySLineData, PySLineGauss]
    data_file = "example_data.obj"
    ser = Ncm.Serialize.new(Ncm.SerializeOpt.NONE)
    if not os.path.exists(data_file):
        print("data file does not exists, run example_create_data.py first.")
        sys.exit(-1)
    else:
        data = ser.from_binfile(data_file)
        assert isinstance(data, (PySLineData, PySLineGauss))
        sld = data

    #
    # New data set object with sld added.
    #
    dset = Ncm.Dataset.new()
    dset.append_data(sld)

    #
    # New likelihood object using dset.
    #
    lh = Ncm.Likelihood.new(dset)

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
    # New Gaussian transition kernel to be used in MCMC algorithm.
    # It was created with size 0 (number of parameters), but once
    # added to the MCMC object the correct size is assigned.
    #
    gtkern = Ncm.MSetTransKernGauss.new(0)
    mcmc = Ncm.FitMCMC.new(fit, gtkern, Ncm.FitRunMsgs.SIMPLE)

    #
    # Getting the Fisher matrix calculated above, scaling it by
    # multiplying by 2 and setting it as the covariance matrix
    # of the transition kernel.
    #
    cov = fit.fstate.covar.dup()
    cov.scale(2.0)
    gtkern.set_cov(cov)

    #
    # Using `example_mcmc_out.fits' as the catalog file, if there
    # is already data in it, the sampler continues from where it stopped.
    #
    mcmc.set_data_file("example_mcmc_out.fits")

    #
    # Running the mcmc, it will first calculate 1000 points, after that
    # it will estimate the error in the parameters mean. Using the current
    # errors the algorithm tries to calculated how many extra steps are
    # necessary to obtain the required error `10^-3' in every parameters,
    # and it will run such extra steps. It will repeat this procedure
    # until it attains the required error in every parameter.
    #
    #
    mcmc.start_run()
    mcmc.run_lre(1000, 1.0e-3)
    mcmc.end_run()

    #
    # Calculates the parameter means and covariance and set it into
    # the fit object and then print.
    #
    mcmc.mean_covar()
    fit.log_covar()


if __name__ == "__main__":
    run_mcmc()

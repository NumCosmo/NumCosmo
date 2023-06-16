#!/usr/bin/env python
#
# example_mc.py
#
# Mon May 22 16:00:00 2023
# Copyright  2023  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# example_mc.py
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

"""Example of running a Monte Carlo analysis."""

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


def run_mc() -> None:
    """Example of running a Monte Carlo analysis.""" ""

    #
    # Instantiating a new SLine model object and setting
    # some values for its parameters.
    #
    slm = PySLineModel()
    slm.props.alpha = 0.9
    slm.props.a = 0.2

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
    #  Printing fitting information.
    #
    fit.log_info()

    #
    #  Calculating the parameters' covariance using numerical differentiation.
    #
    fit.numdiff_m2lnL_covar()

    #
    #  Printing the covariance matrix.
    #
    fit.log_covar()

    #
    # Creates a new Monte Carlo object, using resample from model method.
    # Since no fiducial mset was specified it will use the mset from fit,
    # i.e., it will use the best-fit found above to resample during the MC
    # run.
    #
    mc = Ncm.FitMC.new(fit, Ncm.FitMCResampleType.FROM_MODEL, Ncm.FitRunMsgs.SIMPLE)
    mc.set_nthreads(2)

    #
    # Using `example_mcmc_out.fits' as the catalog file, if there
    # is already data in it, the sampler continues from where it stopped.
    #
    mc.set_data_file("example_mc_out.fits")

    #
    # Running the mcmc, it will first calculate 1000 points, after that
    # it will estimate the error in the parameters mean. Using the current
    # errors the algorithm tries to calculate how many extra steps are
    # necessary to obtain the required error `10^-3' in every parameter,
    # and it will run such extra steps. It will repeat this procedure
    # until it attains the required error in every parameter.
    #
    #
    mc.start_run()
    mc.run_lre(20000, 1.0e-3)
    mc.end_run()

    #
    # Calculates the parameter means and covariance and set it into
    # the fit object and then print.
    #
    mc.mean_covar()
    fit.log_covar()


if __name__ == "__main__":
    run_mc()

#!/usr/bin/env python
#
# example_esmcmc.py
#
# Mon May 22 16:00:00 2023
# Copyright  2023  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# example_esmcmc.py
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
from py_sline_mfunc import PyTestFunc

from numcosmo_py import Ncm

#
#  Initializing the library objects, this must be called before
#  any other library function.
#
Ncm.cfg_init()


def run_esmcmc() -> None:
    """Example of running a MCMC analysis on a SLine model."""

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
    # Printing fitting informations.
    #
    fit.log_info()

    #
    # Setting single thread calculation.
    #
    Ncm.func_eval_set_max_threads(0)
    Ncm.func_eval_log_pool_stats()

    #
    # New Gaussian prior to provide the initial points for the chain.
    # It was created with size 0 (number of parameters), but once
    # initialized with mset the correct size is assigned.
    #
    # The initial sampler will use a diagonal covariance with the
    # diagonal terms being the parameters scale set by each model.
    #
    init_sampler = Ncm.MSetTransKernGauss.new(0)
    init_sampler.set_mset(mset)
    init_sampler.set_prior_from_mset()
    init_sampler.set_cov_from_rescale(1.0)

    #
    # Creates the ESMCMC walker object, this object is responsible
    # for moving the walkers in each interation, the stretch move
    # is affine invariant and therefore gives good results even for
    # very correlated parametric space.
    #
    nwalkers = 300
    apes = Ncm.FitESMCMCWalkerAPES.new(nwalkers, mset.fparams_len())

    #
    # The methods below set the walk scale, which controls the size of the
    # step done between two walkers and circumscribe the walkers inside
    # the box defined by the parameters inside the mset object.
    #
    # stretch.set_scale(3.0)
    # stretch.set_box_mset(mset)

    #
    # Additional functions of mset to be computed during the sampling
    # process.
    #
    mfunc_oa = Ncm.ObjArray.new()

    tf = PyTestFunc()
    mfunc_oa.add(tf)

    #
    # Initialize the ESMCMC object using the objects above. It will
    # use 50 walkers, i.e., each point in the MCMC chain contains
    # 50 points in the parametric space. Each step uses the last point
    # in the chain (the last 50 parametric points) to calculate the
    # proposal points.
    #
    esmcmc = Ncm.FitESMCMC.new_funcs_array(
        fit, nwalkers, init_sampler, apes, Ncm.FitRunMsgs.SIMPLE, mfunc_oa
    )

    #
    # These methods enable the auto-trim options on ESMCMC. This option
    # makes the sampler check the chains' health and trim any unnecessary
    # burn-in part. We set the number of divisions to 100 so we test the
    # chains in blocks of n/100. The last method asserts that each 2min
    # the catalog will be checked.
    #
    # esmcmc.set_auto_trim(True)
    # esmcmc.set_auto_trim_div(100)
    esmcmc.set_max_runs_time(2.0 * 60.0)
    esmcmc.set_nthreads(2)

    #
    # Using `example_esmcmc_out.fits' as the catalog file, if there
    # is already data in it, the sampler continues from where it stopped.
    #
    esmcmc.set_data_file("example_esmcmc_out.fits")

    #
    # Running the esmcmc, it will first calculate 1000 points, after that
    # it will estimate the error in the parameters mean. Using the current
    # errors the algorithm tries to calculate how many extra steps are
    # necessary to obtain the required error `10^-3' in every parameter,
    # and it will run such extra steps. It will repeat this procedure
    # until it attains the required error in every parameter.
    #
    #
    esmcmc.start_run()
    esmcmc.run_lre(10, 1.0e-3)
    esmcmc.end_run()

    #
    # Calculates the parameter means and covariance and set it into
    # the fit object and then print.
    #
    esmcmc.mean_covar()
    fit.log_covar()


if __name__ == "__main__":
    run_esmcmc()

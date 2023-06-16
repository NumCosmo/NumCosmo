#!/usr/bin/env python
#
# example_fit.py
#
# Mon May 22 16:00:00 2023
# Copyright  2023  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# example_fit.py
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

"""Example of fitting a SLine model."""

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


def run_fit() -> None:
    """Example of fitting a SLine model."""
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
    # New data set object with sld added and the likelihood using it.
    #
    dset = Ncm.Dataset.new()
    dset.append_data(sld)
    lh = Ncm.Likelihood.new(dset)

    #
    #  Creating a Fit object of type NLOPT using the fitting algorithm ln-neldermead to
    #  fit the model set mset using the Likelihood lh and using a numerical differentiation
    #  algorithm (NUMDIFF_FORWARD) to obtain the gradient (if needed).
    #
    fit = Ncm.Fit.new(
        Ncm.FitType.NLOPT, "ln-neldermead", lh, mset, Ncm.FitGradType.NUMDIFF_FORWARD
    )

    #
    #  Running the fitter printing messages.
    #
    if not (len(sys.argv) > 1 and sys.argv[1] == "n"):
        fit.run(Ncm.FitRunMsgs.FULL)

    #
    #  Printing fitting information.
    #
    fit.log_info()

    #
    # Calculating the parameters' covariance using numerical differentiation
    # (observed Fisher matrix).
    #
    fit.obs_fisher()
    fit.log_covar()

    #
    # Calculating the parameters' covariance using numerical differentiation
    # (expected Fisher matrix).
    #
    fit.fisher()
    fit.log_covar()


if __name__ == "__main__":
    run_fit()

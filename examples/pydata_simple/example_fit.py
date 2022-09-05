#!/usr/bin/env python

import sys
import math
import numpy as np
import matplotlib.pyplot as plt
import os.path

try:
    import gi

    gi.require_version('NumCosmo', '1.0')
    gi.require_version('NumCosmoMath', '1.0')
except:
    pass

from gi.repository import GObject
from gi.repository import NumCosmo as Nc
from gi.repository import NumCosmoMath as Ncm

from py_sline_model import PySLineModel
from py_sline_data import PySLineData
from py_sline_gauss import PySLineGauss

#
#  Initializing the library objects, this must be called before
#  any other library function.
#
Ncm.cfg_init()

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
sld = None
data_file = "example_data.obj"
ser = Ncm.Serialize.new(0)
if not os.path.exists(data_file):
    print("data file does not exists, run example_create_data.py first.")
    sys.exit(-1)
else:
    sld = ser.from_binfile(data_file)

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
fit = Ncm.Fit.new(Ncm.FitType.NLOPT, "ln-neldermead", lh, mset,
                  Ncm.FitGradType.NUMDIFF_FORWARD)

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
#  Calculating the parameters' covariance using numerical differentiation (observed Fisher matrix).
#
fit.obs_fisher()
fit.log_covar()

#
#  Calculating the parameters' covariance using numerical differentiation (expected Fisher matrix).
#
fit.fisher()
fit.log_covar()

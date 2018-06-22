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
Ncm.cfg_init ()

#
# Instantiating a new SLine model object and setting
# some values for its parameters.
#
slm = PySLineModel ()
slm.props.alpha = 0.9
slm.props.a     = 0.2

#
# New Model set object including slm with parameters
# set as free.
#
mset = Ncm.MSet.empty_new ()
mset.set (slm)
mset.param_set_all_ftype (Ncm.ParamType.FREE)
mset.prepare_fparam_map ()

#
# Creating a new Serialization object, and load
# the data file.
#
sld = None
data_file = "example_data.obj"
ser = Ncm.Serialize.new (0)
if not os.path.exists (data_file):
  print ("data file does not exists, run example_create_data.py first.")
  sys.exit (-1)
else:
  sld = ser.from_binfile (data_file)

#
# New data set object with sld added.
#
dset = Ncm.Dataset.new ()
dset.append_data (sld)

#
# New likelihood object using dset.
#
lh = Ncm.Likelihood.new (dset)

#
#  Creating a Fit object of type NLOPT using the fitting algorithm ln-neldermead to
#  fit the Modelset mset using the Likelihood lh and using a numerical differentiation
#  algorithm (NUMDIFF_FORWARD) to obtain the gradient (if needed).
#
fit = Ncm.Fit.new (Ncm.FitType.NLOPT, "ln-neldermead", lh, mset, Ncm.FitGradType.NUMDIFF_FORWARD)

#
#  Running the fitter printing messages.
#
fit.run (Ncm.FitRunMsgs.SIMPLE)

#
#  Printing fitting informations.
#
fit.log_info ()

#
#  Calculating the parameters covariance using numerical differentiation.
#
fit.numdiff_m2lnL_covar ()

#
#  Printing the covariance matrix.
# 
fit.log_covar ()

#
# Creates a new Monte Carlo object, using resample from model method.
# Since no fiducial mset was specified it will use the mset from fit,
# i.e., it will use the best-fit found above to resample during the MC
# run.
#
mc = Ncm.FitMC.new (fit, Ncm.FitMCResampleType.FROM_MODEL, 
                    Ncm.FitRunMsgs.SIMPLE)
mc.set_nthreads (2)

#
# Using `example_mcmc_out.fits' as the catalog file, if there
# is already data in it, the sampler continues from where it stopped.
#
mc.set_data_file ("example_mc_out.fits")

#
# Running the mcmc, it will first calculate 1000 points, after that
# it will estimate the error in the parameters mean. Using the current
# errors the algorithm tries to calculated how many extra steps are 
# necessary to obtain the required error `10^-3' in every parameters,
# and it will run such extra steps. It will repeat this procedure
# until it attains the required error in every parameter.
# 
#
mc.start_run ()
mc.run_lre (1000, 1.0e-3)
mc.end_run ()

#
# Calculates the parameter means and covariance and set it into 
# the fit object and then print.
# 
mc.mean_covar ()
fit.log_covar ()

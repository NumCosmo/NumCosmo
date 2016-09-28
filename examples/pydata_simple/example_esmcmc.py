#!/usr/bin/python2

import gi
import math
import numpy as np
import matplotlib.pyplot as plt
import os.path

gi.require_version('NumCosmo', '1.0')
gi.require_version('NumCosmoMath', '1.0')

from gi.repository import GObject
from gi.repository import NumCosmo as Nc
from gi.repository import NumCosmoMath as Ncm

from py_sline_model import PySLineModel
from py_sline_data import PySLineData

#
# Instantiating a new empty SLine data object.
#
sld = PySLineData ()

#
# Instantiating a new SLine model object and setting
# some values for its parameters.
#
slm = PySLineModel ()
slm.props.m = 0.9
slm.props.b = 0.1

#
# New Model set object including slm with parameters
# set as free.
#
mset = Ncm.MSet.empty_new ()
mset.set (slm)
mset.param_set_all_ftype (Ncm.ParamType.FREE)
mset.prepare_fparam_map ()

#
# Creating a new Serialization object, with a data
# file does not exists, generate a new sample using
# mset as fiducial model and save it to data_file.
# 
# If data_file already exists, reload sld from it.
#
ser = Ncm.Serialize.new (0)
data_file = "example_data.obj"
if not os.path.exists (data_file):
  rng = Ncm.RNG.new ()
  sld.resample (mset, rng)
  ser.to_file (sld, data_file)
else:
  sld = ser.from_file (data_file)

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
# Printing fitting informations.
#
fit.log_info ()

#
# Setting single thread calculation.
#
Ncm.func_eval_set_max_threads (1)
Ncm.func_eval_log_pool_stats ()

#
# New Gaussian prior to provide the initial points for the chain.
# It was created with size 0 (number of parameters), but once 
# initialized with mset the correct size is assigned. 
#
# The initial sampler will use a diagonal covariance with the
# diagonal terms being the parameters scale set by each model.
#
init_sampler = Ncm.MSetTransKernGauss.new (0)
init_sampler.set_mset (mset)
init_sampler.set_prior_from_mset ()
init_sampler.set_cov_from_rescale (1.0)

#
# Creates the ESMCMC walker object, this object is responsible
# for moving the walkers in each interation, the stretch move
# is affine invariant and therefore gives good results even for
# very correlated parametric space.
# 
nwalkers = 50
stretch = Ncm.FitESMCMCWalkerStretch.new (nwalkers, mset.fparams_len ())

#
# Initialize the ESMCMC object using the objects above. It will
# use 50 walkers, i.e., each point in the MCMC chain contains
# 50 points in the parametric space. Each step uses the last point
# in the chain (the last 50 parametric points) to calculate the
# proposal points.
#
esmcmc  = Ncm.FitESMCMC.new (fit, nwalkers, init_sampler, stretch, Ncm.FitRunMsgs.SIMPLE)

#
# Using `example_esmcmc_out.fits' as the catalog file, if there
# is already data in it, the sampler continues from where it stopped.
# 
esmcmc.set_data_file ("example_esmcmc_out.fits")

#
# Running the esmcmc, it will first calculate 1000 points, after that
# it will estimate the error in the parameters mean. Using the current
# errors the algorithm tries to calculated how many extra steps are 
# necessary to obtain the required error `10^-3' in every parameters,
# and it will run such extra steps. It will repeat this procedure
# until it attains the required error in every parameter.
# 
#
esmcmc.start_run ()
esmcmc.run_lre (1000, 1.0e-3)
esmcmc.end_run ()

#
# Calculates the parameter means and covariance and set it into 
# the fit object and then print.
# 
esmcmc.mean_covar ()
fit.log_covar ()

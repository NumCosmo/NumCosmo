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

#
#  Initializing the library objects, this must be called before
#  any other library function.
#
Ncm.cfg_init ()

dim = 10
#
# Instantiating a new SLine model object and setting
# some values for its parameters.
#
mrb = Ncm.ModelFunnel.new (dim - 1)

#
# New Model set object including slm with parameters
# set as free.
#
mset = Ncm.MSet.empty_new ()
mset.set (mrb)
mset.param_set_all_ftype (Ncm.ParamType.FREE)
mset.prepare_fparam_map ()

#
# Creating a new Serialization object, and load
# the data file.
#
drb = Ncm.DataFunnel.new ()

#
# New data set object with sld added.
#
dset = Ncm.Dataset.new ()
dset.append_data (drb)

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

#fit.run (Ncm.FitRunMsgs.SIMPLE)

#
# Printing fitting informations.
#
fit.log_info ()

#
# Setting single thread calculation.
#
Ncm.func_eval_set_max_threads (0)
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
sampler = 'apes'
#sampler  = 'stretch'
nwalkers = int (math.ceil (3000 * 2))
ssize    = 20000000

if sampler == 'apes':
  walker = Ncm.FitESMCMCWalkerAPES.new (nwalkers, mset.fparams_len ())
  walker.set_over_smooth (0.3)
  #sd0, sd1 = walker.peek_sds ()
  #sd0.set_local_frac (0.03)
  #sd1.set_local_frac (0.03)
elif sampler == "stretch":
  walker = Ncm.FitESMCMCWalkerStretch.new (nwalkers, mset.fparams_len ())

#
# The methods below set the walk scale, which controls the size of the
# step done between two walkers and circumscribe the walkers inside
# the box defined by the parameters inside the mset object.
#
#walker.set_scale (3.0)
#walker.set_box_mset (mset)

#
# Initialize the ESMCMC object using the objects above. It will
# use 50 walkers, i.e., each point in the MCMC chain contains
# 50 points in the parametric space. Each step uses the last point
# in the chain (the last 50 parametric points) to calculate the
# proposal points.
#
esmcmc  = Ncm.FitESMCMC.new (fit, nwalkers, init_sampler, walker, Ncm.FitRunMsgs.SIMPLE)

#
# These methods enable the auto-trim options on ESMCMC. This option 
# makes the sampler check the chains' health and trim any unnecessary 
# burn-in part. We set the number of divisions to 100 so we test the
# chains in blocks of n/100. The last method asserts that each 2min
# the catalog will be checked.
#
#esmcmc.set_auto_trim (True)
#esmcmc.set_auto_trim_div (100)
#esmcmc.set_max_runs_time (2.0 * 60.0)
#esmcmc.set_nthreads (4)
esmcmc.set_data_file ("example_funnel_%d_%s_st_%d.fits" % (dim, sampler, nwalkers))

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
esmcmc.run (ssize / nwalkers)
esmcmc.end_run ()

#
# Calculates the parameter means and covariance and set it into 
# the fit object and then print.
# 
esmcmc.mean_covar ()
fit.log_covar ()

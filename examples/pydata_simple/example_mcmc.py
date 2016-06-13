#!/usr/bin/python2

import gi
import math
import numpy as np
import matplotlib.pyplot as plt

gi.require_version('NumCosmo', '1.0')
gi.require_version('NumCosmoMath', '1.0')

from gi.repository import GObject
from gi.repository import NumCosmo as Nc
from gi.repository import NumCosmoMath as Ncm

from py_sline_model import PySLineModel
from py_sline_data import PySLineData

#
#
#
sld = PySLineData ()

#
#
#
slm = PySLineModel ()
slm.props.m = 0.9
slm.props.b = 0.1

#
#
#
mset = Ncm.MSet.empty_new ()
mset.set (slm)
mset.param_set_all_ftype (Ncm.ParamType.FREE)

#
#
#
sld.simulate_data (slm)

#
#
#
dset = Ncm.Dataset.new ()
dset.append_data (sld)

#
#
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
#
#
init_sampler = Ncm.MSetTransKernGauss.new (0)

#
#
#
gtkern = Ncm.MSetTransKernGauss.new (0)
mcmc   = Ncm.FitMCMC.new (fit, gtkern, Ncm.FitRunMsgs.SIMPLE)
cov    = fit.fstate.covar.dup ()

mcmc.set_data_file ("example_mcmc_out.fits")

cov.scale (2.0)
gtkern.set_cov (cov)

init_sampler.set_mset (mset)
init_sampler.set_prior_from_mset ()

mcmc.start_run ()
mcmc.run_lre (1000, 1.0e-3)
mcmc.end_run ()

mcmc.mean_covar ()
fit.log_covar ()

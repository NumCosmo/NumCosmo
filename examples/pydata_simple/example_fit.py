#!/usr/bin/python2

import sys
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
from py_sline_gauss import PySLineGauss

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
# Creating a new Serialization object, and load
# the data file.
#
sld = None
data_file = "example_data.obj"
ser = Ncm.Serialize.new (0)
if not os.path.exists (data_file):
  print "data file does not exists, run example_create_data.py first."
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
#  fit the model set mset using the Likelihood lh and using a numerical differentiation
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

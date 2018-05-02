#!/usr/bin/env python

try:
  import gi
  gi.require_version('NumCosmo', '1.0')
  gi.require_version('NumCosmoMath', '1.0')
except:
  pass

import timeit
import sys
import time
import math
import numpy as np
from gi.repository import NumCosmo as Nc
from gi.repository import NumCosmoMath as Ncm


#
#  Initializing the library objects, this must be called before 
#  any other library function.
#
sys.argv = Ncm.cfg_init_full (sys.argv)


cosmo = Nc.HICosmo.new_from_name (Nc.HICosmo, "NcHICosmoDEXcdm")

cosmo.props.H0      = 70.0
cosmo.props.Omegab  = 0.05
cosmo.props.Omegac  = 0.25
cosmo.props.Omegax  = 0.70
cosmo.props.Tgamma0 = 2.72
cosmo.props.w       = -1.0

#
#  Creating a new Modelset and set cosmo as the HICosmo model to be used.
#
mset = Ncm.MSet ()
mset.set (cosmo)

#
#  Setting parameters Omega_c and w to be fitted.
#
cosmo.props.Omegac_fit = True
cosmo.props.w_fit = True

#
#  Creating a new Distance object optimized to redshift 2.
#
dist = Nc.Distance (zf = 2.0)

#
#  Creating a new Data object from distance modulus catalogs.
#
if False:
  data_snia = Nc.DataDistMu.new_from_id (dist, Nc.DataSNIAId.SIMPLE_UNION2_1)
else:
  snia = Nc.SNIADistCov.new (dist, 4)
  mset.set (snia)
  mset.param_set_mid_ftype (Nc.SNIADistCov.id (), Ncm.ParamType.FIXED)
  mset.param_set_ftype (Nc.SNIADistCov.id (), Nc.SNIADistCovSParams.ALPHA, Ncm.ParamType.FREE)
  mset.param_set_ftype (Nc.SNIADistCov.id (), Nc.SNIADistCovSParams.BETA,  Ncm.ParamType.FREE)
  #mset.param_set_ftype (Nc.SNIADistCov.id (), Nc.SNIADistCovSParams.M1,    Ncm.ParamType.FREE)
  #mset.param_set_ftype (Nc.SNIADistCov.id (), Nc.SNIADistCovSParams.M2,    Ncm.ParamType.FREE)
  data_snia = Nc.DataSNIACov.new (False)
  Nc.data_snia_load_cat (data_snia, Nc.DataSNIAId.COV_JLA_SNLS3_SDSS_SYS_STAT_CMPL)

#
#  Creating a new Dataset and add snia to it.
#
dset = Ncm.Dataset ()
dset.append_data (data_snia)

#
#  Creating a Likelihood from the Dataset.
#
lh = Ncm.Likelihood (dataset = dset)

#
#  Creating a Fit object of type NLOPT using the fitting algorithm ln-neldermead to
#  fit the Modelset mset using the Likelihood lh and using a numerical differentiation
#  algorithm (NUMDIFF_FORWARD) to obtain the gradient (if needed).
#
fit = Ncm.Fit.new (Ncm.FitType.NLOPT, "ln-neldermead", lh, mset, Ncm.FitGradType.NUMDIFF_FORWARD)

mj  = Ncm.MPIJobFit.new (fit, None)

ser = Ncm.Serialize.new (0)
mj.init_all_slaves (ser)

#time.sleep (10000)

rng = Ncm.RNG.new ()
rng.set_random_seed (False)

param_a = []
m2lnL_a = []
for i in range (24):
  w_i  = rng.gaussian_gen (-1.0,  0.05)
  Oc_i = rng.gaussian_gen ( 0.25, 0.01)  
  a_i  = rng.gaussian_gen ( 0.141, 0.01)
  b_i  = rng.gaussian_gen ( 3.101, 0.1)
  param_a.append (Ncm.Vector.new_array ([Oc_i, w_i, a_i, b_i]))
  m2lnL_a.append (Ncm.Vector.new (1 + 4))

#mj.run_array (param_a, m2lnL_a)

t1 = timeit.Timer('mj.run_array (param_a, m2lnL_a)', "from __main__ import mj, param_a, m2lnL_a")
print ("Timming: ", t1.timeit (1))

for param, m2lnL in zip (param_a, m2lnL_a):
  fit.params_set_vector (param)
  print ("% -22.15g % -22.15g % -22.15g % -22.15g % -22.15g" % (m2lnL.get (0), m2lnL.get (1), m2lnL.get (2), m2lnL.get (3), m2lnL.get (4)))


mj.free_all_slaves ()




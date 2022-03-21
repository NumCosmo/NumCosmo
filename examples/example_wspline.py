#!/usr/bin/env python

try:
  import gi
  gi.require_version('NumCosmo', '1.0')
  gi.require_version('NumCosmoMath', '1.0')
except:
  pass

import math
import numpy as np
from gi.repository import NumCosmo as Nc
from gi.repository import NumCosmoMath as Ncm

#
#  Initializing the library objects, this must be called before 
#  any other library function.
#
Ncm.cfg_init ()

#
#  New homogeneous and isotropic cosmological model NcHICosmoDEXcdm 
#  with one massive neutrino.
#
zf = 3.0
cosmo  = Nc.HICosmoDEWSpline.new (12, zf)
cosmo2 = Nc.HICosmoDEXcdm.new ()
zf = 3000.0

#
#  New cosmological distance objects optimizied to perform calculations
#  up to redshift 2.0.
#
dist = Nc.Distance.new (zf)

#
#  Setting values for the cosmological model, those not set stay in the
#  default values. Remember to use the _orig_ version to set the original
#  parameters when a reparametrization is used.
#

#
# C-like
#
cosmo.orig_param_set (Nc.HICosmoDESParams.H0,       70.00)
cosmo.orig_param_set (Nc.HICosmoDESParams.OMEGA_C,   0.25)
cosmo.orig_param_set (Nc.HICosmoDESParams.OMEGA_X,   0.70)
cosmo.orig_param_set (Nc.HICosmoDESParams.T_GAMMA0,  2.72)
cosmo.orig_param_set (Nc.HICosmoDESParams.OMEGA_B,   0.05)

cosmo2.orig_param_set (Nc.HICosmoDESParams.H0,       70.00)
cosmo2.orig_param_set (Nc.HICosmoDESParams.OMEGA_C,   0.25)
cosmo2.orig_param_set (Nc.HICosmoDESParams.OMEGA_X,   0.70)
cosmo2.orig_param_set (Nc.HICosmoDESParams.T_GAMMA0,  2.72)
cosmo2.orig_param_set (Nc.HICosmoDESParams.OMEGA_B,   0.05)
cosmo2.orig_param_set (Nc.HICosmoDEXCDMSParams.W,    -1.0)

alpha_a = np.array (cosmo.get_alpha ().dup_array ())
#w_a = -1.0 + 1.0 * alpha_a / alpha_a[-1]
w_a = [-1.0] * 12

cosmo.orig_vparam_set_vector (Nc.HICosmoDEWSplineVParams.W, Ncm.Vector.new_array (w_a))
#cosmo.orig_vparam_set_vector (Nc.HICosmoDEWSplineVParams.W, Ncm.Vector.new_array (w_a))

#
#  Printing the parameters used.
#
print ("# Model parameters: ")
cosmo.params_log_all ()

dist.prepare (cosmo)

#
#  Printing some distances up to redshift 1.0.
#

N = 20
for i in range (0, N):
  z  = zf / (N - 1.0) * i
  
  w = cosmo.w_de (z)
  E2               = cosmo.E2 (z)
  E2Omega_de       = cosmo.E2Omega_de (z)
  dE2Omega_de_dz   = cosmo.dE2Omega_de_dz (z)
  d2E2Omega_de_dz2 = cosmo.d2E2Omega_de_dz2 (z)

  print ("% 22.15f w % 22.15g E2 % 22.15g E2Omega_de % 22.15g dE2Omega_de_dz % 22.15g d2E2Omega_de_dz2 % 22.15g" % (z, w, E2, E2Omega_de, dE2Omega_de_dz, d2E2Omega_de_dz2))

  w = cosmo2.w_de (z)
  E2               = cosmo2.E2 (z)
  E2Omega_de       = cosmo2.E2Omega_de (z)
  dE2Omega_de_dz   = cosmo2.dE2Omega_de_dz (z)
  d2E2Omega_de_dz2 = cosmo2.d2E2Omega_de_dz2 (z)
  
  print ("% 22.15f w % 22.15g E2 % 22.15g E2Omega_de % 22.15g dE2Omega_de_dz % 22.15g d2E2Omega_de_dz2 % 22.15g" % (z, w, E2, E2Omega_de, dE2Omega_de_dz, d2E2Omega_de_dz2))



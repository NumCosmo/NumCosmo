#!/usr/bin/env python

try:
  import gi
  gi.require_version('NumCosmo', '1.0')
  gi.require_version('NumCosmoMath', '1.0')
except:
  pass

import math
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
cosmo = Nc.HICosmo.new_from_name (Nc.HICosmo, "NcHICosmoDEXcdm{'massnu-length':<1>}")
cosmo.set_reparam (Nc.HICosmoDEReparamCMB.new (cosmo.len ()))

#
#  New cosmological distance objects optimizied to perform calculations
#  up to redshift 2.0.
#
dist = Nc.Distance.new (2.0)

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
cosmo.orig_param_set (Nc.HICosmoDEXCDMSParams.W,    -1.10)

cosmo.orig_vparam_set (Nc.HICosmoDEVParams.M, 0, 0.06)

#
# OO-like
#
cosmo.props.H0      = 70.00
cosmo.props.Omegab  =  0.04
cosmo.props.Omegac  =  0.25
cosmo.props.Omegax  =  0.70
cosmo.props.Tgamma0 =  2.72
cosmo.props.w       = -1.10

massnu_v = Ncm.Vector.new_array ([0.06])
cosmo.props.massnu  = massnu_v

#
#  Printing the parameters used.
#
print ("# Model parameters: ")
cosmo.params_log_all ()

dist.prepare (cosmo)

#
#  Printing some distances up to redshift 1.0.
#
N      = 20
RH_Mpc = cosmo.RH_Mpc ()

for i in range (0, N):
  z  = 1.0 / (N - 1.0) * i
  Dc = dist.comoving (cosmo, z)
  dc = RH_Mpc * Dc
  
  print ("% 10.8f % 22.15g [c/H0] % 22.15g [Mpc]" % (z, Dc, dc))



#!/usr/bin/env python

import math
import gi
gi.require_version('NumCosmo', '1.0')
gi.require_version('NumCosmoMath', '1.0')

from gi.repository import NumCosmo as Nc
from gi.repository import NumCosmoMath as Ncm

#
#  Initializing the library objects, this must be called before 
#  any other library function.
#
Ncm.cfg_init ()

#
#  New homogeneous and isotropic cosmological model NcHICosmoDEXcdm 
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
cosmo.orig_param_set (Nc.HICosmoDESParams.ENNU,      2.0328)
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
cosmo.props.ENnu    =  2.0328
cosmo.props.w       = -1.10

massnu_v = Ncm.Vector.new_array ([0.06])
cosmo.props.massnu  = massnu_v

cosmo.omega_x2omega_k ()
cosmo.param_set_by_name ("Omegak", 0.0)

#
#  Printing the parameters used.
#
print ("# Model parameters: ", end='', flush = True) 
cosmo.params_log_all ()

#
#  New CLASS backend object.
#
cbe = Nc.CBE.new ()
cbe.set_calc_transfer (True)

#
# Submodels necessary for CLASS
#
reion = Nc.HIReionCamb.new ()
prim  = Nc.HIPrimPowerLaw.new ()
cosmo.add_submodel (reion)
cosmo.add_submodel (prim)

#
# Preparing CLASS backend
#
cbe.prepare (cosmo)
dist.prepare (cosmo)

print ("# theta100CMB % 22.15e" % (dist.theta100CMB (cosmo)))
print ("# zt          % 22.15e" % (cosmo.zt (5.0)))
print ("# Omega_mnu0  % 22.15e" % (cosmo.Omega_mnu0 ()))
print ("# Press_mnu0  % 22.15e" % (cosmo.Press_mnu0 ()))
print ("# Omega_k0    % 22.15e" % (cosmo.Omega_k0 ()))

ztest = 1.0e4
print ("# E2Omega_mnu   (% 22.15g) % 22.15e" % (ztest, cosmo.E2Omega_mnu (ztest)))
print ("# E2Press_mnu   (% 22.15g) % 22.15e" % (ztest, cosmo.E2Press_mnu (ztest)))
print ("# E2Omega_mnu_d (% 22.15g) % 22.15e" % (ztest, cosmo.E2Omega_mnu (ztest) - 3.0 * cosmo.E2Press_mnu (ztest)))
print ("# E2Omega_b     (% 22.15g) % 22.15e" % (ztest, cosmo.E2Omega_b (ztest)))
print ("# E2Omega_c     (% 22.15g) % 22.15e" % (ztest, cosmo.E2Omega_c (ztest)))
print ("# Neff           % 22.15e" % (cosmo.Neff ()))

#
# Printing comparison between CLASS and NumCosmo background
#
cbe.compare_bg (cosmo, True)


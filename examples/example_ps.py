#!/usr/bin/python2

import gi
gi.require_version('NumCosmo', '1.0')
gi.require_version('NumCosmoMath', '1.0')

import sys
import time
from math import *
import numpy as np
from gi.repository import GObject
import matplotlib
import matplotlib.pyplot as plt
from gi.repository import NumCosmo as Nc
from gi.repository import NumCosmoMath as Ncm

matplotlib.rcParams.update({'font.size': 11})

#
#  Initializing the library objects, this must be called before 
#  any other library function.
#
Ncm.cfg_init ()

#
#  New homogeneous and isotropic cosmological model NcHICosmoDEXcdm 
#
cosmo = Nc.HICosmo.new_from_name (Nc.HICosmo, "NcHICosmoDEXcdm{'massnu-length':<1>, 'Tnu-length':<1>}")
cosmo.omega_x2omega_k ()
cosmo.param_set_by_name ("Omegak", 0.0)
cosmo.param_set_by_name ("w", -1.0)
cosmo.param_set_by_name ("Omegab", 0.04909244421)
cosmo.param_set_by_name ("Omegac", 0.26580755578)

reion = Nc.HIReionCamb.new ()
prim  = Nc.HIPrimPowerLaw.new ()

cosmo.param_set_by_name ("H0", 67.31)

prim.param_set_by_name ("n_SA", 0.9658)
prim.param_set_by_name ("ln10e10ASA", 3.0904)

reion.param_set_by_name ("z_re", 9.9999)

cosmo.add_submodel (reion)
cosmo.add_submodel (prim)

ps_cbe  = Nc.PowspecMLCBE.new ()
ps_eh   = Nc.PowspecMLTransfer.new (Nc.TransferFuncEH.new ())

z_min = 0.0
z_max = 2.0
zdiv  = 0.49999999999

k_min = 1.0e-5
k_max = 1.0e3

nk = 2000
nR = 2000
Rh8 = 8.0 / cosmo.h ()

ps_cbe.set_kmin (k_min)
ps_eh.set_kmin (k_min)

ps_cbe.set_kmax (k_max)
ps_eh.set_kmax (k_max)

ps_cbe.require_zi (z_min)
ps_cbe.require_zf (z_max)

ps_eh.require_zi (z_min)
ps_eh.require_zf (z_max)

ps_eh.prepare (cosmo)
ps_cbe.prepare (cosmo)

for z in np.arange (z_min, z_max, (z_max - z_min) * zdiv):
  k_a = []
  Pk_eh_a = []
  Pk_cbe_a = []
  for lnk in np.arange (log (ps_cbe.props.kmin), log (ps_cbe.props.kmax), log (k_max / k_min) / nk):
    k = exp (lnk)
    k2 = k * k
    k3 = k2 * k
    k_a.append (k)
    Pk_cbe_a.append (k3 * ps_cbe.eval (cosmo, z, k))
    Pk_eh_a.append (k3 * ps_eh.eval (cosmo, z, k))
    #print k, ps_eh.eval (cosmo, z, k) / ps_cbe.eval (cosmo, z, k)
  
  plt.plot (k_a, Pk_cbe_a, label = r'CLASS $z = %.2f$' % (z))
  plt.plot (k_a, Pk_eh_a, label  = r'EH    $z = %.2f$' % (z))

plt.legend (loc="lower right")
plt.xscale('log')
plt.yscale('log')
plt.savefig ("ps_cbe_eh.png")
plt.clf ()

for z in np.arange (z_min, z_max, (z_max - z_min) * zdiv * 0.5):
  k_a = []
  Pk_eh_a = []
  Pk_cbe_a = []
  for lnk in np.arange (log (ps_cbe.props.kmin), log (ps_cbe.props.kmax), log (k_max / k_min) / nk):
    k = exp (lnk)
    k2 = k * k
    k3 = k2 * k
    k_a.append (k)
    Pk_cbe_a.append (k3 * ps_cbe.eval (cosmo, z, k))
    Pk_eh_a.append (k3 * ps_eh.eval (cosmo, z, k))
    #print k, ps_eh.eval (cosmo, z, k) / ps_cbe.eval (cosmo, z, k)
  
  plt.plot (k_a, np.abs (1.0 - np.array (Pk_eh_a) / np.array (Pk_cbe_a)), label = r'CLASS EH cmp $z = %.2f$' % (z))

plt.legend (loc="lower right")
plt.xscale('log')
plt.yscale('log')
plt.savefig ("ps_diff_cbe_eh.png")
plt.clf ()

#
# Filtering 
#
psf_cbe = Ncm.PowspecFilter.new (ps_cbe, Ncm.PowspecFilterType.TOPHAT)
psf_eh  = Ncm.PowspecFilter.new (ps_eh, Ncm.PowspecFilterType.TOPHAT)

psf_cbe.set_best_lnr0 ()
psf_eh.set_best_lnr0 ()

psf_cbe.prepare (cosmo)
psf_eh.prepare (cosmo)


print "# CBE sigma8 = % 20.15g, EH sigma8 = % 20.15g" % (psf_cbe.eval_sigma (0.0, Rh8), psf_eh.eval_sigma (0.0, Rh8))
print "# kmin % 20.15g kmax % 20.15g" % (ps_cbe.props.kmin, ps_cbe.props.kmax)
print "# Rmin % 20.15g Rmax % 20.15g" % (psf_cbe.get_r_min (), psf_cbe.get_r_max ())

lnRmin = log (psf_cbe.get_r_min ())
lnRmax = log (psf_cbe.get_r_max ())

for z in np.arange (z_min, z_max, (z_max - z_min) * zdiv):
  Rh_a = []
  sigma_eh_a  = []
  sigma_cbe_a = []
  for lnR in np.arange (lnRmin, lnRmax, (lnRmax - lnRmin) / nR):
    R  = exp (lnR)
    Rh = R * cosmo.h ()
    
    Rh_a.append (Rh)
    sigma_cbe_a.append (psf_cbe.eval_sigma_lnr (z, lnR))
    sigma_eh_a.append (psf_eh.eval_sigma_lnr (z, lnR))

  plt.plot (Rh_a, sigma_cbe_a, label = r'$\sigma$ CLASS $z = %.2f$' % (z))
  plt.plot (Rh_a, sigma_cbe_a, label = r'$\sigma$ EH    $z = %.2f$' % (z))

plt.legend (loc="lower left")
plt.xscale('log')
plt.yscale('log')
plt.savefig ("ps_var_cbe_eh.png")
plt.clf ()

for z in np.arange (z_min, z_max, (z_max - z_min) * zdiv):
  Rh_a = []
  dvar_eh_a   = []
  dvar_cbe_a  = []
  for lnR in np.arange (lnRmin, lnRmax, (lnRmax - lnRmin) / nR):
    R  = exp (lnR)
    Rh = R * cosmo.h ()
    
    Rh_a.append (Rh)
    dvar_cbe_a.append (psf_cbe.eval_dlnvar_dlnr (z, lnR))
    dvar_eh_a.append (psf_eh.eval_dlnvar_dlnr (z, lnR))

  plt.plot (Rh_a, dvar_cbe_a, label  = r'$\frac{\mathrm{d}\ln\sigma^2}{\mathrm{d}\ln{}R}$ CLASS $z = %.2f$' % (z))
  plt.plot (Rh_a, dvar_eh_a, label   = r'$\frac{\mathrm{d}\ln\sigma^2}{\mathrm{d}\ln{}R}$ EH    $z = %.2f$' % (z))

plt.legend (loc="lower left")
plt.xscale('log')
plt.savefig ("ps_dvar_cbe_eh.png")
plt.clf ()

zmaxnl = 10.0
z_max  = zmaxnl
pshf   = Nc.PowspecMNLHaloFit.new (ps_cbe, zmaxnl, 1.0e-5)

pshf.set_kmin (k_min)
pshf.set_kmax (k_max)
pshf.require_zi (z_min)
pshf.require_zf (z_max)

pshf.prepare (cosmo)

for z in np.arange (z_min, z_max, (z_max - z_min) * zdiv):
  k_a = []
  Pk_cbe_a = []
  Pknl_cbe_a = []
  for lnk in np.arange (log (ps_cbe.props.kmin), log (ps_cbe.props.kmax), log (k_max / k_min) / nk):
    k = exp (lnk)
    k2 = k * k
    k3 = k2 * k
    k_a.append (k)
    Pk_cbe_a.append (k3 * ps_cbe.eval (cosmo, z, k))
    Pknl_cbe_a.append (k3 * pshf.eval (cosmo, z, k))
  
  plt.plot (k_a, Pk_cbe_a, label = r'CLASS $z = %.2f$' % (z))
  plt.plot (k_a, Pknl_cbe_a, label  = r'CLASS+HaloFit $z = %.2f$' % (z))

plt.legend (loc="lower right")
plt.xscale('log')
plt.yscale('log')
plt.savefig ("ps_cbe_halofit.png")
plt.clf ()

psf_cbenl = Ncm.PowspecFilter.new (pshf, Ncm.PowspecFilterType.TOPHAT)
psf_cbenl.set_best_lnr0 ()
psf_cbenl.prepare (cosmo)

print "# CBE+HaloFit sigma8 = % 20.15g" % (psf_cbenl.eval_sigma (0.0, Rh8))


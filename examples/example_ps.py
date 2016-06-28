#!/usr/bin/python2

import gi
gi.require_version('NumCosmo', '1.0')
gi.require_version('NumCosmoMath', '1.0')

import sys
import time
from math import *
import numpy as np
from gi.repository import GObject
import matplotlib.pyplot as plt
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
cosmo = Nc.HICosmo.new_from_name (Nc.HICosmo, "NcHICosmoDEXcdm")
cosmo.omega_x2omega_k ()
cosmo.param_set_by_name ("Omegak", 0.0)
cosmo.param_set_by_name ("w", -1.0)
cosmo.param_set_by_name ("Omegab", 0.04909244421)
cosmo.param_set_by_name ("Omegac", 0.26580755578)

reion = Nc.HIReionCamb.new ()
prim  = Nc.HIPrimPowerLaw.new ()

cosmo.param_set_by_name ("H0", 67.31)
cosmo.param_set_by_name ("ns", 0.9658)

prim.param_set_by_name ("n_SA", 0.9658)
prim.param_set_by_name ("ln10e10ASA", 3.0904)

reion.param_set_by_name ("z_re", 9.999)

cosmo.add_submodel (reion)
cosmo.add_submodel (prim)

ps_cbe  = Nc.PowspecMLCBE.new ()
ps_eh   = Nc.PowspecMLTransfer.new (Nc.TransferFuncEH.new ())
ps_bbks = Nc.PowspecMLTransfer.new (Nc.TransferFuncBBKS.new ())

ps_cbe.set_kmax (1.0e1)
ps_eh.set_kmax (1.0e1)
ps_bbks.set_kmax (1.0e1)

ps_eh.prepare (cosmo)
ps_bbks.prepare (cosmo)
ps_cbe.prepare (cosmo)

for z in np.arange (0.0, 1.1, 0.1):
  k_a = []
  Pk_a = []
  Pk_bbks_a = []
  Pk_cbe_a = []
  for lnk in np.arange (log (ps_cbe.props.kmin), log (ps_cbe.props.kmax), log (1.0e2) / 100.0):
    k = exp (lnk)
    k2 = k * k
    k3 = k2 * k
    k_a.append (k)
    Pk_a.append (k3 * ps_eh.eval (cosmo, z, k))
    Pk_bbks_a.append (k3 * ps_bbks.eval (cosmo, z, k))
    Pk_cbe_a.append (k3 * ps_cbe.eval (cosmo, z, k))
    #print z, k, ps_eh.eval (cosmo, z, k), fac * ps_cbe.eval (cosmo, z, k)
  
  plt.plot (k_a, Pk_a, label = "EH")
  plt.plot (k_a, Pk_bbks_a, label = "BBKS")
  plt.plot (k_a, Pk_cbe_a, label = "CLASS")
  plt.legend (loc="lower right")
  plt.xscale('log')
  plt.yscale('log')
  #plt.show()

#
# Filtering 
#
psf = Ncm.PowspecFilter.new (ps_cbe, Ncm.PowspecFilterType.TOPHAT)

psf.prepare (cosmo)

wp =  Nc.Window.new_from_name ("NcWindowTophat")
tf = Nc.TransferFunc.new_from_name ("NcTransferFuncEH")
vp = Nc.MatterVar.new (Nc.MatterVarStrategy.FFT, wp, tf)
vp.prepare (cosmo)

np = 2000
divfac = 1.0 / (np - 1.0)

factor = psf.eval_sigma (0.0, 8.0 / cosmo.h ()) / sqrt (vp.var0 (cosmo, log (8.0)))

print "# sigma8 == % 20.15g" % (psf.eval_sigma (0.0, 8.0 / cosmo.h ()))
print "# sigma8 == % 20.15g" % (factor * sqrt (vp.var0 (cosmo, log (8.0))))

print "# kmin % 20.15g kmax % 20.15g" % (ps_cbe.props.kmin, ps_cbe.props.kmax)
print "# Rmin % 20.15g Rmax % 20.15g" % (psf.get_r_min (), psf.get_r_max ())

lnRmin = log (psf.get_r_min () * cosmo.h ())
lnRmax = log (psf.get_r_max () * cosmo.h ())

print "# lnRmin % 20.15g lnRmax % 20.15g" % (exp (lnRmin), exp (lnRmax))
z = 0.0

for i in range (0, np):
  lnRh = lnRmin +  (lnRmax - lnRmin) * divfac * i
  Rh   = exp (lnRh)
  R    = Rh / cosmo.h ()
  lnR  = log (R)
  
  sigma2    = vp.var0 (cosmo, lnRh)
  dlnsigma2 = vp.dlnvar0_dlnR (cosmo, lnRh)
  
  sigma_n     = psf.eval_sigma (z, R)
  dlnsigma2_n = psf.eval_dlnvar_dlnr (z, lnR)

  print lnRh, Rh, factor * sqrt (sigma2), sigma_n, dlnsigma2_n

print "\n"


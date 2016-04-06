#!/usr/bin/python2

import gi
gi.require_version('NumCosmo', '1.0')
gi.require_version('NumCosmoMath', '1.0')

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
#cosmo.param_set_by_name ("Omegac", 0.2)

reion = Nc.HIReionCamb.new ()
prim  = Nc.HIPrimPowerLaw.new ()

prim.param_set_by_name ("n_SA", 0.5)

cosmo.add_submodel (reion)
cosmo.add_submodel (prim)

ps_cbe  = Nc.PowspecMLCBE.new ()
ps_eh   = Nc.PowspecMLTransfer.new (Nc.TransferFuncEH.new ())
ps_bbks = Nc.PowspecMLTransfer.new (Nc.TransferFuncBBKS.new ())

ps_eh.prepare (cosmo)
ps_bbks.prepare (cosmo)
ps_cbe.prepare (cosmo)

#fac_eh = ps_cbe.eval (cosmo, 0.0, 0.001) / ps_eh.eval (cosmo, 0.0, 0.001) 
#fac_bbks = ps_cbe.eval (cosmo, 0.0, 0.001) / ps_bbks.eval (cosmo, 0.0, 0.001)

for z in np.arange (0.0, 0.1, 0.1):
  k_a = []
  Pk_a = []
  Pk_bbks_a = []
  Pk_cbe_a = []
  for lnk in np.arange (log (1.0e-3), log (1.0e-1), log (1.0e2) / 100.0):
    k = exp (lnk)
    k_a.append (k)
    Pk_a.append (ps_eh.eval (cosmo, z, k))
    Pk_bbks_a.append (ps_bbks.eval (cosmo, z, k))
    Pk_cbe_a.append (ps_cbe.eval (cosmo, z, k))
#    print z, k, ps_eh.eval (cosmo, z, k), fac * ps_cbe.eval (cosmo, z, k)
  
  plt.plot (k_a, Pk_a, label = "EH")
  plt.plot (k_a, Pk_bbks_a, label = "BBKS")
  plt.plot (k_a, Pk_cbe_a, label = "CLASS")
  plt.legend (loc="upper right")
  plt.yscale('log')
  plt.show()



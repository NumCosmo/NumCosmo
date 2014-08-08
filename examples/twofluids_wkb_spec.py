#!/usr/bin/python2

from gi.repository import NumCosmo as Nc
from gi.repository import NumCosmoMath as Ncm
import matplotlib.pyplot as plt
import time
import math

#
#  Initializing the library objects, this must be called before 
#  any other library function.
#
Ncm.cfg_init ()

#
#  New homogeneous and isotropic cosmological model NcHICosmoQGRW
#
cosmo = Nc.HICosmo.new_from_name (Nc.HICosmo, "NcHICosmoQGRW")

w        = 1.0e-8
prec     = 1.0e-9
k_min    = 1.0e-1
k_max    = 1.0e4
wkb_prec = prec

cosmo.props.w      = w
cosmo.props.Omegar = 1.0e-5
cosmo.props.Omegaw = 1.0 - 1.0e-5

pertZ = Nc.HIPertTwoFluids.new ()
pertQ = Nc.HIPertTwoFluids.new ()

pertZ.props.reltol = prec
pertQ.props.reltol = prec

Delta_zeta_A = []
Delta_zeta_B = []
Delta_Q_A = []
Delta_Q_B = []

k_a = []

alpha_minf      = -cosmo.abs_alpha (1.0e-27)
alpha_end       = cosmo.abs_alpha (1.0e25)  

for i in range (101):
  mode_k = k_min * math.exp (math.log (k_max / k_min) * i / 100.0)
  
  k_a.append (mode_k)

  pertZ.set_mode_k (mode_k);
  pertQ.set_mode_k (mode_k);

  alpha_zeta_wkb  = pertZ.wkb_zeta_maxtime (cosmo, alpha_minf, -alpha_end)
  alpha_S_wkb     = pertZ.wkb_S_maxtime (cosmo, alpha_minf, -alpha_end)   
  alpha_zeta_prec = pertZ.wkb_zeta_maxtime_prec (cosmo, wkb_prec, alpha_minf, -alpha_end)
  alpha_S_prec    = pertZ.wkb_S_maxtime_prec (cosmo, wkb_prec, alpha_minf, -alpha_end)   
  alphai          = -cosmo.abs_alpha (cosmo.x_alpha (min (alpha_S_prec, alpha_zeta_prec)) * 1.0e-1)
  
  print "# Mode[%d] k = %f, prec = %e" % (i, mode_k, prec)
  print "# Maxtime wkb zeta: %f %e" % (alpha_zeta_wkb, cosmo.x_alpha (alpha_zeta_wkb))
  print "# Maxtime wkb Q:    %f %e" % (alpha_S_wkb, cosmo.x_alpha (alpha_S_wkb))
  print "# Prec    wkb zeta: %f %e" % (alpha_zeta_prec, cosmo.x_alpha (alpha_zeta_prec))
  print "# Prec    wkb Q:    %f %e" % (alpha_S_prec, cosmo.x_alpha (alpha_S_prec))
  print "# Inital time:      %f %e" % (alphai, cosmo.x_alpha (alphai))

  pertZ.set_stiff_solver (True)
  pertQ.set_stiff_solver (True)

  print "# Preparing zeta"
  pertZ.prepare_patched_zeta (cosmo, wkb_prec, alphai, -alpha_end)
  print "# Preparing Q"
  pertZ.prepare_patched_S (cosmo, wkb_prec, alphai, -alpha_end)
  
  print "# Preparing zeta"
  pertQ.prepare_patched_zeta (cosmo, wkb_prec, alphai, -alpha_end)
  print "# Preparing Q"
  pertQ.prepare_patched_S (cosmo, wkb_prec, alphai, -alpha_end)

  alphai = -cosmo.abs_alpha (1.0e-18)
  
  print "# Setting inital conditions"
  pertZ.set_init_cond_patched_zeta (cosmo, alphai)
  pertQ.set_init_cond_patched_Q (cosmo, alphai)

  varsZ = [1.234] * 8
  varsQ = [1.234] * 8
  varsZ_wkb = [1.234] * 8
  varsQ_wkb = [1.234] * 8

  pertZ.evolve (cosmo, alpha_end)
  pertQ.evolve (cosmo, alpha_end)
  (alphas, varsZ) = pertZ.get_values (varsZ)
  (alphas, varsQ) = pertQ.get_values (varsQ)

  Delta_zeta_A.append (mode_k**3 * math.hypot (varsZ[0], varsZ[1])**2)
  Delta_zeta_B.append (mode_k**3 * math.hypot (varsQ[0], varsQ[1])**2)
  Delta_Q_A.append (mode_k**3 * math.hypot (varsZ[4], varsZ[5])**2)
  Delta_Q_B.append (mode_k**3 * math.hypot (varsQ[4], varsQ[5])**2)

#
# Plotting results
#

plt.title (r"Mode $\vert\zeta\vert^2$,  k = " + str (mode_k) + ",  w = " + str (w))
plt.xscale('log')
plt.yscale('log')
plt.plot (k_a, Delta_zeta_A, 'r', label=r"$\Delta_{\zeta_A}$")
plt.plot (k_a, Delta_zeta_B, 'b', label=r"$\Delta_{\zeta_B}$")
plt.xlabel (r'$k$')
plt.ylabel (r'$\Delta_\zeta$')
plt.legend (loc='best')

plt.savefig ("pspec_zeta_w_" + str (w) + ".pdf")



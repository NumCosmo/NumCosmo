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

w      = 1.0e-4
prec   = 1.0e-11
k_min  = 1.0e-1
k_max  = 1.0e4

cosmo.props.w      = w
cosmo.props.Omegar = 1.0e-5
cosmo.props.Omegaw = 1.0 - 1.0e-5

pertZ = Nc.HIPertTwoFluids.new ()
pertQ = Nc.HIPertTwoFluids.new ()

pertZ.props.reltol = prec
pertQ.props.reltol  = prec

Delta_zeta_A = []
Delta_zeta_B = []
Delta_Q_A = []
Delta_Q_B = []

k_a = []

for i in range (1001):
  mode_k = k_min * math.exp (math.log (k_max / k_min) * i / 1000.0)
  
  k_a.append (mode_k)

  pertZ.set_mode_k (mode_k);
  pertQ.set_mode_k (mode_k);

  alphaZt = pertZ.wkb_zeta_maxtime (cosmo, -cosmo.abs_alpha (1.0e-30), -cosmo.abs_alpha (1.0e20))
  alphaQt = pertZ.wkb_S_maxtime (cosmo, -cosmo.abs_alpha (1.0e-30), -cosmo.abs_alpha (1.0e20))

  xi = min (cosmo.x_alpha (min (alphaZt, alphaQt)) * 1.0e-3, 1.0e-3)

  alphai = -cosmo.abs_alpha (xi)
  alphaf = cosmo.abs_alpha (1.0e20)

  print "# Mode[%d] k = %f" % (i, mode_k)
  print "# Maxtime zeta", alphaZt, cosmo.x_alpha (alphaZt)
  print "# Maxtime Q   ", alphaQt, cosmo.x_alpha (alphaQt)
  print "# Initial time", alphai, cosmo.x_alpha (alphai)

  #pert.set_stiff_solver (True)
  pertZ.prepare_wkb_zeta (cosmo, alphai, alphaZt)
  pertZ.prepare_wkb_S (cosmo, alphai, alphaQt)

  pertQ.prepare_wkb_zeta (cosmo, alphai, alphaZt)
  pertQ.prepare_wkb_S (cosmo, alphai, alphaQt)

  pertZ.set_init_cond_wkb_zeta (cosmo, alphai)
  pertQ.set_init_cond_wkb_Q (cosmo, alphai)

  varsZ = [1.234] * 8
  varsQ = [1.234] * 8
  varsZ_wkb = [1.234] * 8
  varsQ_wkb = [1.234] * 8

  pertZ.evolve (cosmo, alphaf)
  pertQ.evolve (cosmo, alphaf)
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



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
prec   = 1.0e-9
mode_k = 1.0

cosmo.props.w      = w
cosmo.props.Omegar = 1.0e-5
cosmo.props.Omegaw = 1.0 - 1.0e-5

pertZ = Nc.HIPertTwoFluids.new ()
pertQ = Nc.HIPertTwoFluids.new ()

pertZ.props.reltol = prec
pertZ.set_mode_k (mode_k);

pertQ.props.reltol  = prec
pertQ.set_mode_k (mode_k);

wkb_prec = prec

alpha_minf      = -cosmo.abs_alpha (1.0e-30)
alpha_end       = cosmo.abs_alpha (1.0e25)
alpha_zeta_wkb  = pertZ.wkb_zeta_maxtime (cosmo, alpha_minf, -alpha_end)
alpha_S_wkb     = pertZ.wkb_S_maxtime (cosmo, alpha_minf, -alpha_end)
alpha_zeta_prec = pertZ.wkb_zeta_maxtime_prec (cosmo, wkb_prec, alpha_minf, -alpha_end)
alpha_S_prec    = pertZ.wkb_S_maxtime_prec (cosmo, wkb_prec, alpha_minf, -alpha_end)
alphai          = -cosmo.abs_alpha (cosmo.x_alpha (min (alpha_S_prec, alpha_zeta_prec)) * 1.0e-1)

print "# Maxtime wkb zeta: %f %e" % (alpha_zeta_wkb, cosmo.x_alpha (alpha_zeta_wkb))
print "# Maxtime wkb Q:    %f %e" % (alpha_S_wkb, cosmo.x_alpha (alpha_S_wkb))
print "# Prec    wkb zeta: %f %e" % (alpha_zeta_prec, cosmo.x_alpha (alpha_zeta_prec))
print "# Prec    wkb Q:    %f %e" % (alpha_S_prec, cosmo.x_alpha (alpha_S_prec))
print "# Inital time:      %f %e" % (alphai, cosmo.x_alpha (alphai))


print "# Preparing zeta"
pertZ.prepare_patched_zeta (cosmo, wkb_prec, alphai, -alpha_end)
print "# Preparing Q"
pertZ.prepare_patched_S (cosmo, wkb_prec, alphai, -alpha_end)

print "# Preparing zeta"
pertQ.prepare_patched_zeta (cosmo, wkb_prec, alphai, -alpha_end)
print "# Preparing Q"
pertQ.prepare_patched_S (cosmo, wkb_prec, alphai, -alpha_end)

pertZ.set_stiff_solver (True)
pertQ.set_stiff_solver (True)

#alphai = alpha_S_prec
#alphai = -cosmo.abs_alpha (1.0e-20)

print "# Setting inital conditions"
pertZ.set_init_cond_patched_zeta (cosmo, alphai)
pertQ.set_init_cond_patched_Q (cosmo, alphai)

varsZ = [1.234] * 8
varsQ = [1.234] * 8
varsZ_wkb = [1.234] * 8
varsQ_wkb = [1.234] * 8
 
alpha_a = []

Re_zeta_A = []
Im_zeta_A = []
Re_zeta_B = []
Im_zeta_B = []
zeta_A = []
zeta_B = []

Re_Q_A = []
Im_Q_A = []
Re_Q_B = []
Im_Q_B = []
Q_A = []
Q_B = []

for i in range (10000):
  alpha = alphai + (alpha_end - alphai) / 10000.0 * (i + 1)
  pertQ.evolve (cosmo, alpha)
  pertZ.evolve (cosmo, alpha)
  (alphas, varsZ) = pertZ.get_values (varsZ)
  (alphas, varsQ) = pertQ.get_values (varsQ)
  
  (varsZ_wkb[0], varsZ_wkb[1], varsZ_wkb[2], varsZ_wkb[3]) = pertZ.patched_zeta_Pzeta (cosmo, alphas)
  (varsQ_wkb[4], varsQ_wkb[5], varsQ_wkb[6], varsQ_wkb[7]) = pertQ.patched_Q_PQ (cosmo, alphas)
  
#  varsZ_wkb = pertZ.wkb_full_zeta (cosmo, alphas, varsZ_wkb)
#  varsQ_wkb = pertQ.wkb_full_Q (cosmo, alphas, varsQ_wkb)
  
  alpha_a.append (alphas)
  
  main_zeta = math.hypot (varsZ[0], varsZ[1])
  zeta_wkb  = math.hypot (varsZ_wkb[0], varsZ_wkb[1])
  main_Q    = math.hypot (varsQ[4], varsQ[5])
  Q_wkb     = math.hypot (varsQ_wkb[4], varsQ_wkb[5])

  print alphas, cosmo.x_alpha (alphas), main_zeta, zeta_wkb, main_Q, Q_wkb, math.fabs ((zeta_wkb - main_zeta) / main_zeta), math.fabs ((Q_wkb - main_Q) / main_Q)

  Re_zeta_A.append (varsZ[0])
  Im_zeta_A.append (varsZ[1])
  Re_zeta_B.append (varsQ[0])
  Im_zeta_B.append (varsQ[1])
  
  zeta_A.append (math.hypot (varsZ[0], varsZ[1])**2)
  zeta_B.append (math.hypot (varsQ[0], varsQ[1])**2)

  Re_Q_A.append (varsZ[4])
  Im_Q_A.append (varsZ[5])
  Re_Q_B.append (varsQ[4])
  Im_Q_B.append (varsQ[5])

  Q_A.append (math.hypot (varsZ[4], varsZ[5])**2)
  Q_B.append (math.hypot (varsQ[4], varsQ[5])**2)



#
#  Ploting ionization history.
#

plt.title (r"Mode $\vert\zeta\vert^2$,  k = " + str (mode_k) + ",  w = " + str (w))
plt.yscale('log')
plt.plot (alpha_a, zeta_A, 'r', label=r"$\vert\zeta_A\vert^2$")
plt.plot (alpha_a, zeta_B, 'b', label=r"$\vert\zeta_B\vert^2$")
plt.xlabel(r'$\alpha$')
plt.ylabel(r'$\zeta$')
plt.legend(loc=2)

plt.savefig ("mode_zeta_k_" + str (mode_k) + "_w_" + str (w) + ".pdf")

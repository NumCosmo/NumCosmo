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
mode_k = 1.0e-1

cosmo.props.w      = w
cosmo.props.Omegar = 1.0e-5
cosmo.props.Omegaw = 1.0 - 1.0e-5

pertZ = Nc.HIPertTwoFluids.new ()
pertQ = Nc.HIPertTwoFluids.new ()


pertZ.props.reltol = prec
pertZ.set_mode_k (mode_k);

pertQ.props.reltol  = prec
pertQ.set_mode_k (mode_k);

alphaZt = pertZ.wkb_zeta_maxtime (cosmo, -cosmo.abs_alpha (1.0e-30), -cosmo.abs_alpha (1.0e20))
alphaQt = pertZ.wkb_S_maxtime (cosmo, -cosmo.abs_alpha (1.0e-30), -cosmo.abs_alpha (1.0e20))

alphai = -cosmo.abs_alpha (cosmo.x_alpha (min (alphaZt, alphaQt)) * 1.0e-3)
alphaf = cosmo.abs_alpha (1.0e20)
print "# Maxtime zeta", alphaZt, cosmo.x_alpha (alphaZt)
print "# Maxtime Q   ", alphaQt, cosmo.x_alpha (alphaQt)
print "# Initial time", alphai, cosmo.x_alpha (alphai)

pertZ.prepare_wkb_zeta (cosmo, alphai, alphaZt)
pertZ.prepare_wkb_S (cosmo, alphai, alphaQt)

pertQ.prepare_wkb_zeta (cosmo, alphai, alphaZt)
pertQ.prepare_wkb_S (cosmo, alphai, alphaQt)

#pert.set_stiff_solver (True)

pertZ.set_init_cond_wkb_zeta (cosmo, alphai)
pertQ.set_init_cond_wkb_Q (cosmo, alphai)

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

for i in range (100000):
  alpha = alphai + (alphaf - alphai) / 100000.0 * (i + 1)
  pertZ.evolve (cosmo, alpha)
  pertQ.evolve (cosmo, alpha)
  (alphas, varsZ) = pertZ.get_values (varsZ)
  (alphas, varsQ) = pertQ.get_values (varsQ)
  varsZ_wkb = pertZ.wkb_full_zeta (cosmo, alphas, varsZ_wkb)
  varsQ_wkb = pertQ.wkb_full_Q (cosmo, alphas, varsQ_wkb)
  
  alpha_a.append (alphas)

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

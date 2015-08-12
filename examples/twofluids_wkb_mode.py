#!/usr/bin/python2

from gi.repository import GObject
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

w      = 1.0e-16
prec   = 1.0e-10
mode_k = 1.0e0

cosmo.props.w      = w
cosmo.props.Omegar = 1.0e-5
cosmo.props.Omegaw = 1.0 - 1.0e-5
cosmo.props.xb     = 1.e30

pert = Nc.HIPertTwoFluids.new ()

pert.props.reltol = prec
pert.set_mode_k (mode_k);

wkb_prec = prec

alpha_minf      = -cosmo.abs_alpha (1.0e-40)
alpha_end       = cosmo.abs_alpha (1.0e25)
alpha_zeta_wkb  = pert.wkb_zeta_maxtime (cosmo, alpha_minf, -alpha_end)
alpha_S_wkb     = pert.wkb_S_maxtime (cosmo, alpha_minf, -alpha_end)
alpha_zeta_prec = pert.wkb_zeta_maxtime_prec (cosmo, Nc.HIPertWKBCmp.POTENTIAL, wkb_prec, alpha_minf, -alpha_end)
alpha_S_prec    = pert.wkb_S_maxtime_prec (cosmo, Nc.HIPertWKBCmp.POTENTIAL, wkb_prec, alpha_minf, -alpha_end)
alphai          = -cosmo.abs_alpha (cosmo.x_alpha (min (alpha_S_prec, alpha_zeta_prec)) * 1.0)

print "# Maxtime wkb zeta: %f %e" % (alpha_zeta_wkb, cosmo.x_alpha (alpha_zeta_wkb))
print "# Maxtime wkb Q:    %f %e" % (alpha_S_wkb, cosmo.x_alpha (alpha_S_wkb))
print "# Prec    wkb zeta: %f %e" % (alpha_zeta_prec, cosmo.x_alpha (alpha_zeta_prec))
print "# Prec    wkb Q:    %f %e" % (alpha_S_prec, cosmo.x_alpha (alpha_S_prec))
alpha_zeta_prec_a2 = pert.wkb_zeta_maxtime_prec (cosmo, Nc.HIPertWKBCmp.ALPHA2, wkb_prec**2, alpha_minf, -alpha_end)
alpha_S_prec_a2    = pert.wkb_S_maxtime_prec (cosmo, Nc.HIPertWKBCmp.ALPHA2, wkb_prec**2, alpha_minf, -alpha_end)
print "# Prec a2 wkb zeta: %f %e" % (alpha_zeta_prec_a2, cosmo.x_alpha (alpha_zeta_prec_a2))
print "# Prec a2 wkb Q:    %f %e" % (alpha_S_prec_a2, cosmo.x_alpha (alpha_S_prec_a2))
alpha_zeta_prec_a2 = pert.wkb_zeta_maxtime_prec (cosmo, Nc.HIPertWKBCmp.ALPHA2, wkb_prec**1, alpha_minf, -alpha_end)
alpha_S_prec_a2    = pert.wkb_S_maxtime_prec (cosmo, Nc.HIPertWKBCmp.ALPHA2, wkb_prec**1, alpha_minf, -alpha_end)
print "# Prec a2 wkb zeta: %f %e" % (alpha_zeta_prec_a2, cosmo.x_alpha (alpha_zeta_prec_a2))
print "# Prec a2 wkb Q:    %f %e" % (alpha_S_prec_a2, cosmo.x_alpha (alpha_S_prec_a2))
print "# Inital time:      %f %e" % (alphai, cosmo.x_alpha (alphai))

print "# Preparing zeta"
#pert.prepare_wkb_zeta (cosmo, wkb_prec, alpha_minf, -alpha_end)
print "# Preparing Q"
#pert.prepare_wkb_S (cosmo, wkb_prec, alpha_minf, -alpha_end)

pert.set_stiff_solver (True)

vars = [1.234] * 8
vars_wkb = [1.234] * 8
 
alpha_A = []
alpha_B = []

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
  eom_full = pert.eom_full (cosmo, alpha)
  
  nup = eom_full.nu_plus
  num = eom_full.nu_minus
  laq = eom_full.lambda_s
  laz = eom_full.lambda_zeta
  Uplus  = eom_full.Uplus
  Uminus = eom_full.Uminus
  Wplus  = eom_full.Wplus
  Wminus = eom_full.Wminus
  U11    = eom_full.US[0]
  U22    = eom_full.US[1]
  U12    = eom_full.US[2]
  
#  qpz = math.sqrt (nuz / nup) * laq
#  qpq = math.sqrt (nuq / nup) * laz
#  qmz = math.sqrt (nuz / num) * laz
#  qmq = math.sqrt (nuq / num) * laq
  
#  Ppz = math.sqrt (nup / nuz) * laq  
#  Ppq = math.sqrt (nup / nuq) * laz
#  Pmz = math.sqrt (num / nuz) * laz  
#  Pmq = math.sqrt (num / nuq) * laq

  #print ("% 20.15g % 20.15g % 20.15g % 20.15g % 20.15g % 20.15g" % (alpha, pert.nuA (cosmo, alpha), pert.nuB (cosmo, alpha), pert.YAB (cosmo, alpha), pert.LA (cosmo, alpha), pert.LB (cosmo, alpha)))
  #print ("% 20.15f % 20.15e % 20.15e % 20.15e % 20.15e % 20.15e % 20.15e % 20.15e % 20.15e" % (alpha, nuz, nuq, y, nu1, nup, num, laq, laz))
  #print ("% 20.15f % 20.15e % 20.15e % 20.15e % 20.15e % 20.15e % 20.15e % 20.15e % 20.15e" % (alpha, qpz, qpq, qmz, qmq, Ppz, Ppq, Pmz, Pmq))
  print ("% 20.15f % 20.15e % 20.15e % 20.15e % 20.15e % 20.15e % 20.15e % 20.15e % 20.15e % 20.15e % 20.15e % 20.15e" % (alpha, nup, num, laq, laz, Uplus, Uminus, Wplus, Wminus, U11, U22, U12))


exit ()


print "# Setting inital conditions (zeta)"
pert.set_init_cond_wkb_zeta (cosmo, alphai)

for i in range (10000):
  alpha = alphai + (alpha_end - alphai) / 10000.0 * (i + 1)
  pert.evolve (cosmo, alpha)

  (alphas, vars) = pert.get_values (vars)  
  (vars_wkb[0], vars_wkb[1], vars_wkb[2], vars_wkb[3]) = pert.patched_zeta_Pzeta (cosmo, alphas)
  vars_wkb = pert.wkb_full_zeta (cosmo, alphas, vars_wkb)
  
  alpha_A.append (alphas)
  
  main_zeta = math.hypot (vars[0], vars[1])
  zeta_wkb  = math.hypot (vars_wkb[0], vars_wkb[1])
  
  sub_Q = math.hypot (vars[4], vars[5])
  Q_wkb = math.hypot (vars_wkb[4], vars_wkb[5])

  print alphas, cosmo.x_alpha (alphas), main_zeta, zeta_wkb, math.fabs ((zeta_wkb - main_zeta) / main_zeta), sub_Q, Q_wkb, math.fabs ((Q_wkb - sub_Q) / sub_Q)

  Re_zeta_A.append (vars[0])
  Im_zeta_A.append (vars[1])
  
  zeta_A.append (math.hypot (vars[0], vars[1])**2)

  Re_Q_A.append (vars[4])
  Im_Q_A.append (vars[5])

  Q_A.append (math.hypot (vars[4], vars[5])**2)

print "# Setting inital conditions (Q)"
pert.set_init_cond_patched_Q (cosmo, alphai)

for i in range (10000):
  alpha = alphai + (alpha_end - alphai) / 10000.0 * (i + 1)
  pert.evolve (cosmo, alpha)

  (alphas, vars) = pert.get_values (vars)  
  (vars_wkb[4], vars_wkb[5], vars_wkb[6], vars_wkb[7]) = pert.patched_Q_PQ (cosmo, alphas)
  vars_wkb = pert.wkb_full_Q (cosmo, alphas, vars_wkb)
  
  alpha_B.append (alphas)
  
  main_Q = math.hypot (vars[4], vars[5])
  Q_wkb  = math.hypot (vars_wkb[4], vars_wkb[5])

  sub_zeta = math.hypot (vars[0], vars[1])
  zeta_wkb = math.hypot (vars_wkb[0], vars_wkb[1])

  print alphas, cosmo.x_alpha (alphas), main_Q, Q_wkb, math.fabs ((Q_wkb - main_Q) / main_Q), sub_zeta, zeta_wkb, math.fabs ((zeta_wkb - sub_zeta) / sub_zeta)

  Re_zeta_B.append (vars[0])
  Im_zeta_B.append (vars[1])
  
  zeta_B.append (math.hypot (vars[0], vars[1])**2)

  Re_Q_B.append (vars[4])
  Im_Q_B.append (vars[5])

  Q_B.append (math.hypot (vars[4], vars[5])**2)


#
#  
#

plt.title (r"Mode $\vert\zeta\vert^2$,  k = " + str (mode_k) + ",  w = " + str (w))
plt.yscale('log')
plt.plot (alpha_A, zeta_A, 'r', label=r"$\vert\zeta_A\vert^2$")
plt.plot (alpha_B, zeta_B, 'b', label=r"$\vert\zeta_B\vert^2$")
plt.xlabel(r'$\alpha$')
plt.ylabel(r'$\zeta$')
plt.legend(loc=2)

plt.savefig ("mode_zeta_k_" + str (mode_k) + "_w_" + str (w) + ".pdf")

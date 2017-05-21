#!/usr/bin/python2

try:
  import gi
  gi.require_version('NumCosmo', '1.0')
  gi.require_version('NumCosmoMath', '1.0')
except:
  pass

import math
import matplotlib
import matplotlib.pyplot as plt
from gi.repository import NumCosmo as Nc
from gi.repository import NumCosmoMath as Ncm

import numpy as np

#
#  Initializing the library objects, this must be called before 
#  any other library function.
#
Ncm.cfg_init ()

adiab = Nc.HIPertAdiab.new ()
gw    = Nc.HIPertGW.new ()
Vexp  = Nc.HICosmoVexp.new ()

Vexp.props.alphab   = +1.0e-30
Vexp.props.sigmaphi = +6.0e-1
Vexp.props.dphi     = -3.0e-1
Vexp.props.xb       = 1.0e37
Vexp.props.OmegaL   = 1.0
Vexp.props.Omegac   = 1.0
Vexp.props.H0       = 67.8

print "# (t0, t1) = (% 21.15e, % 21.15e)" % (Vexp.tau_min (), Vexp.tau_max ())

k = 1.0e0

adiab.set_ti (Vexp.tau_min ())
adiab.set_tf (Vexp.tau_max ())
adiab.set_k (k)
adiab.set_reltol (1.0e-13)

gw.set_ti (Vexp.tau_min ())
gw.set_tf (Vexp.tau_max ())
gw.set_k (k)
gw.set_reltol (1.0e-13)

gw.prepare (Vexp)

adiab.prepare (Vexp)

(t0, t1) = adiab.get_t0_t1 (Vexp)

print "# (t0, t1) = (% 21.15e, % 21.15e)" % (t0, t1)

t_a = np.linspace (t0, t1, 1000000)

Delta_zeta  = []
Delta_Pzeta = []

zeta_a_a = []
zeta_b_a = []

h_a_a = []
h_b_a = []

Delta_h  = []
Delta_Ph = []

lp2_RH2 = Vexp.RH_planck ()**(-2.0)
lp_RH   = math.sqrt (lp2_RH2)

for t in t_a:
  (Delta_zeta_v, Delta_Pzeta_v) = adiab.eval_Delta (Vexp, t) 
  (Delta_h_v, Delta_Ph_v)       = gw.eval_Delta (Vexp, t) 

  (zeta_a_v, zeta_b_v, Pzeta_a_v, Pzeta_b_v)  = adiab.eval_QV (Vexp, t) 
  (h_a_v, h_b_v, Ph_a_v, Ph_b_v)              = gw.eval_QV (Vexp, t) 

  Delta_zeta.append (lp2_RH2 * Delta_zeta_v)
  Delta_h.append (lp2_RH2 * Delta_h_v)
  
  zeta_a_a.append (lp_RH * zeta_a_v)
  h_a_a.append (lp_RH * h_a_v)

mylw = 1

#plt.plot (eta_a, RLnI,   lw=mylw, label = r'$I$')
#plt.plot (eta_a, ILnI,   lw=mylw, label = r'$J$')

#plt.plot (eta_a, RQ,   lw=mylw, label = r'$\theta$')
#plt.plot (eta_a, IQ,   lw=mylw, label = r'$\psi$')

plt.plot (t_a, Delta_zeta, lw=mylw, label = r'$\Delta_{\zeta}$')
plt.plot (t_a, Delta_h,    lw=mylw, label = r'$\Delta_{h}$')

plt.grid (b=True, which='both', linestyle=':', color='0.75', linewidth=0.5)
leg = plt.legend (loc="upper left")

#plt.xscale('symlog', linthreshx=1.0e-5)
#plt.yscale('symlog', linthreshy=1.0e-20)
plt.yscale('log', linthreshy=1.0e-20)

plt.show ()
plt.clf ()


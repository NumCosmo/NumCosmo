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
import scipy.integrate as integrate

#
#  Initializing the library objects, this must be called before 
#  any other library function.
#
Ncm.cfg_init ()

adiab = Nc.HIPertAdiab.new ()
gw    = Nc.HIPertGW.new ()
Vexp  = Nc.HICosmoVexp.new ()

Vexp.props.alphab   = +1.0e-30
Vexp.props.sigmaphi = +5.0e-1
Vexp.props.dphi     = -5.0e-2
Vexp.props.xb       = 1.0e38
Vexp.props.OmegaL   = 1.0
Vexp.props.Omegac   = 1.0
Vexp.props.H0       = 67.8

k = 1.0e-2
tc = Vexp.tau_xe (1.0e15)

adiab.set_ti (Vexp.tau_min ())
adiab.set_tf (tc)
adiab.set_k (k)
adiab.set_reltol (1.0e-14)

gw.set_ti (Vexp.tau_min ())
gw.set_tf (tc)
gw.set_k (k)
gw.set_reltol (1.0e-14)

#If = lambda x: gw.eval_nu (Vexp, x, k) / gw.eval_mnu (Vexp, x, k)
#print integrate.quad (If, -2.0, 2.0)
#for n in np.logspace (-40, -1, 10000):
#  print n, If (n), If (-n)
#exit ()

print "# Preparing ADIAB"
adiab.prepare (Vexp)

print "# Preparing GW"
gw.prepare (Vexp)

#(t0, t1) = adiab.get_t0_t1 (Vexp)
(t0, t1) = gw.get_t0_t1 (Vexp)

print "# BACKG (t0, t1) = (% 21.15e, % 21.15e)" % (Vexp.tau_min (), Vexp.tau_max ())
print "# ADIAB (t0, t1) = (% 21.15e, % 21.15e)" % (t0, t1)

(Delta_zeta_c, Delta_Pzeta_c) = adiab.eval_Delta (Vexp, tc)
(Delta_h_c, Delta_Ph_c)       = gw.eval_Delta (Vexp, tc)

print "# Time of x_e = 10^15:     tau_c     = % 21.15f" % tc
print "# Power spectrum at tau_c: PS_ADIAB  = % 21.15e" % Delta_zeta_c
print "# Power spectrum at tau_c: PS_GW     = % 21.15e" % Delta_h_c
print "# Power spectrum at tau_c: r         = % 21.15e" % (2.0 * Delta_h_c / Delta_zeta_c)

t_a = np.linspace (t0, t1, 1000000)

Delta_zeta  = []
Delta_Pzeta = []

Delta_h  = []
Delta_Ph = []

zeta_a_a  = []
zeta_b_a  = []
Pzeta_a_a = []
Pzeta_b_a = []

h_a_a  = []
h_b_a  = []
Ph_a_a = []
Ph_b_a = []

nu_a     = []
m_a      = []
mnu_a    = []
dlnmnu_a = []

epsilon_a    = []
gamma_a      = []
sin_thetab_a = []
cos_thetab_a = []

for t in t_a:
  (Delta_zeta_v, Delta_Pzeta_v) = adiab.eval_Delta (Vexp, t) 
  (Delta_h_v, Delta_Ph_v)       = gw.eval_Delta (Vexp, t) 

  (zeta_a_v, zeta_b_v, Pzeta_a_v, Pzeta_b_v)  = adiab.eval_QV (Vexp, t) 
  (h_a_v, h_b_v, Ph_a_v, Ph_b_v)              = gw.eval_QV (Vexp, t) 

  (epsilon_v, gamma_v, sin_thetab_v, cos_thetab_v) = gw.eval_AA (Vexp, t) 

  epsilon_a.append (epsilon_v)
  gamma_a.append (gamma_v)
  sin_thetab_a.append (sin_thetab_v)
  cos_thetab_a.append (cos_thetab_v)

  nu_adiab     = adiab.eval_nu (Vexp, t, k)
  mnu_adiab    = adiab.eval_mnu (Vexp, t, k)
  m_adiab      = mnu_adiab / nu_adiab
  dlnmnu_adiab = math.fabs (adiab.eval_dlnmnu (Vexp, t, k))

  nu_a.append (nu_adiab)
  m_a.append (m_adiab)
  mnu_a.append (mnu_adiab)
  dlnmnu_a.append (dlnmnu_adiab)

  Delta_zeta.append (Delta_zeta_v)
  Delta_h.append (Delta_h_v)
  
  zeta_a_a.append (zeta_a_v)
  zeta_b_a.append (zeta_b_v)

  Pzeta_a_a.append (Pzeta_a_v)
  Pzeta_b_a.append (Pzeta_b_v)

  h_a_a.append (h_a_v)
  h_b_a.append (h_b_v)

  Ph_a_a.append (Ph_a_v)
  Ph_b_a.append (Ph_b_v)

mylw = 1

h_a_a  = np.array (h_a_a)
h_b_a  = np.array (h_b_a)

Ph_a_a = np.array (Ph_a_a)
Ph_b_a = np.array (Ph_b_a)

zeta_a_a  = np.array (zeta_a_a)
zeta_b_a  = np.array (zeta_b_a)

Pzeta_a_a = np.array (Pzeta_a_a)
Pzeta_b_a = np.array (Pzeta_b_a)

#plt.plot (t_a, Delta_zeta, lw=mylw, label = r'$\Delta_{\zeta}$')
#plt.plot (t_a, Delta_h,    lw=mylw, label = r'$\Delta_{h}$')

plt.plot (t_a, zeta_a_a, lw=mylw, label = r'$\tilde{\zeta}^a$')
plt.plot (t_a, zeta_b_a, lw=mylw, label = r'$\tilde{\zeta}^b$')
plt.plot (t_a, h_a_a,    lw=mylw, label = r'$\tilde{h}^a$')
plt.plot (t_a, h_b_a,    lw=mylw, label = r'$\tilde{h}^b$')

#plt.plot (t_a, epsilon_a,    lw=mylw, label = r'$\epsilon$')
#plt.plot (t_a, gamma_a,      lw=mylw, label = r'$\gamma$')
#plt.plot (t_a, sin_thetab_a, lw=mylw, label = r'$\sin(\theta_b)$')
#plt.plot (t_a, cos_thetab_a, lw=mylw, label = r'$\cos(\theta_b)$')

#plt.plot (t_a, nu_a,     lw=mylw, label = r'$\nu_h$')
#plt.plot (t_a, m_a,      lw=mylw, label = r'$m_h$')
#plt.plot (t_a, mnu_a,    lw=mylw, label = r'$m_h\nu_h$')
#plt.plot (t_a, 1.0 / (np.array (mnu_a)),    lw=mylw, label = r'$(m_h\nu_h)^{-1}$')
#plt.plot (t_a, dlnmnu_a, lw=mylw, label = r'$d\ln(m_h\nu_h)$')

#plt.plot (t_a, (Ph_a_a * h_b_a), lw=mylw, label = r't_1')
#plt.plot (t_a, (Ph_b_a * h_a_a), lw=mylw, label = r't_2')

plt.grid (b=True, which='both', linestyle=':', color='0.75', linewidth=0.5)
leg = plt.legend (loc="upper left")

#plt.xscale('symlog', linthreshx=1.0e-5)
plt.yscale('symlog', linthreshy=1.0e-120)
#plt.yscale('log', linthreshy=1.0e-20)

plt.show ()
plt.clf ()


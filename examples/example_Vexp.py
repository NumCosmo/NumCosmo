#!/usr/bin/env python

try:
  import gi
  gi.require_version('NumCosmo', '1.0')
  gi.require_version('NumCosmoMath', '1.0')
except:
  pass

import math
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from gi.repository import NumCosmo as Nc
from gi.repository import NumCosmoMath as Ncm

#
#  Initializing the library objects, this must be called before 
#  any other library function.
#
Ncm.cfg_init ()

#
#  New homogeneous and isotropic cosmological model: NcHICosmoVexp 
#
cosmo1 = Nc.HICosmo.new_from_name (Nc.HICosmo, "NcHICosmoVexp")
cosmo2 = Nc.HICosmo.new_from_name (Nc.HICosmo, "NcHICosmoVexp")

cosmo1.props.alphab   = +1.0e-1
cosmo1.props.sigmaphi = +8.0e-1
cosmo1.props.dphi     = +5.0e-1
#cosmo1.props.xb       = 4.8e37
cosmo1.props.xb       = 1.0e37
cosmo1.props.OmegaL   = 1.0e0
cosmo1.props.Omegac   = 1.0
cosmo1.props.H0       = 67.8

cosmo2.props.alphab   = +1.0e-1
cosmo2.props.sigmaphi = +8.0e-1
cosmo2.props.dphi     = +5.0e-1
#cosmo2.props.xb       = 4.8e37
cosmo2.props.xb       = 1.0e37
cosmo2.props.OmegaL   = 1.0e0
cosmo2.props.Omegac   = 1.0
cosmo2.props.H0       = 67.8

tau_min = max (cosmo1.tau_min (), cosmo1.tau_xc (1.0e-5))
tau_max = min (cosmo1.tau_max (), -cosmo1.tau_xc (1.0e-5))
k       = 1.0e0
xb      = cosmo1.xbe ()

print ("# tau interval:  (% 22.15g, % 22.15g)" % (tau_min, tau_max))

Nn    = 100000
LS    = 1.0e-30
tau_a = np.concatenate ((np.geomspace (tau_min, -LS, Nn/2), np.geomspace (LS, tau_max, Nn/2)), axis = 0)

alpha_a     = []

nu1_a        = []
mnu1_a       = []
mnu_gw1_a    = []
dlnmnu1_a    = []
dlnmnu_gw1_a = []
Eaa01_a      = []
x1_a         = []
y1_a         = []

nu2_a        = []
mnu2_a       = []
mnu_gw2_a    = []
dlnmnu2_a    = []
dlnmnu_gw2_a = []
Eaa02_a      = []
x2_a         = []
y2_a         = []

for tau in tau_a:
  nu1, dlnmnu1 = cosmo1.eval_system (tau, k)
  mnu1         = cosmo1.eval_mnu (tau, k)
  mnu_gw1      = Nc.HIPertIGW.eval_mnu (cosmo1, tau, k)

  nu2, dlnmnu2 = cosmo2.eval_system (tau, k)
  mnu2         = cosmo2.eval_mnu (tau, k)
  mnu_gw2      = Nc.HIPertIGW.eval_mnu (cosmo2, tau, k)

  nu_gw1, dlnmnu_gw1 = Nc.HIPertIGW.eval_system (cosmo1, tau, k)  
  nu_gw2, dlnmnu_gw2 = Nc.HIPertIGW.eval_system (cosmo2, tau, k)  
  
  Eaa01 = k * tau / nu1;
  Eaa02 = k * tau / nu2;
  
  (x1, y1) = cosmo1.x_y (tau)
  (x2, y2) = cosmo2.x_y (tau)

  alpha = 0
  if tau > 0.0:
    alpha = 0.5 * tau * tau / math.log (10.0)
  else:
    alpha = - 0.5 * tau * tau / math.log (10.0)

  alpha_a.append (alpha)

  nu1_a.append (nu1)
  mnu1_a.append (nu1 / mnu1)
  mnu_gw1_a.append (nu1 / mnu_gw1)
  
  dlnmnu1_a.append (math.fabs (k * dlnmnu1 / nu1))
  dlnmnu_gw1_a.append (math.fabs (k * dlnmnu_gw1 / nu1))
  Eaa01_a.append (math.fabs (Eaa01) * math.exp (-0.5 * tau * tau) * xb)
  
  x1_a.append (x1)
  y1_a.append (y1)

  nu2_a.append (nu2)
  mnu2_a.append (nu2 / mnu2)
  mnu_gw2_a.append (nu2 / mnu_gw2)
  
  dlnmnu2_a.append (math.fabs (k * dlnmnu2 / nu2))
  dlnmnu_gw2_a.append (math.fabs (k * dlnmnu_gw2 / nu2))
  Eaa02_a.append (math.fabs (Eaa02) * math.exp (-0.5 * tau * tau) * xb)
  
  x2_a.append (x2)
  y2_a.append (y2)

mylw = 1

print ("# tau classical: (% 22.15g, % 22.15g)" % (cosmo1.tau_qt_c (), cosmo1.tau_qt_e ()))

#plt.plot (alpha_a, np.array (x1_a)**(-2),     lw=mylw, label = r'$x^{-2}_1$')
#plt.plot (alpha_a, np.array (x2_a)**(-2),     lw=mylw, label = r'$x^{-2}_2$')
#plt.plot (alpha_a, nu_a,     lw=mylw, label = r'$\nu$')
#plt.plot (alpha_a, mnu1_a,       lw=mylw, label = r'$m_1\nu_1$')
#plt.plot (alpha_a, mnu2_a,       lw=mylw, label = r'$m_2\nu_2$')
#plt.plot (alpha_a, mnu_gw_a,    lw=mylw, label = r'$m_\mathrm{gw}\nu_\mathrm{gw}$')
plt.plot (alpha_a, dlnmnu1_a, lw=mylw, label = r'$\mathrm{d}\ln(m\nu)/\mathrm{d}\tau$')
#plt.plot (alpha_a, dlnmnu_gw1_a, lw=mylw, label = r'$\mathrm{d}\ln(m\nu)_\mathrm{gw}/\mathrm{d}\tau$')
plt.plot (alpha_a, dlnmnu2_a, lw=mylw, label = r'$\mathrm{d}\ln(m\nu)/\mathrm{d}\tau$')
#plt.plot (alpha_a, dlnmnu_gw2_a, lw=mylw, label = r'$\mathrm{d}\ln(m\nu)_\mathrm{gw}/\mathrm{d}\tau$')
plt.plot (alpha_a, Eaa01_a, lw=mylw, label = r'$E_1$')
plt.plot (alpha_a, Eaa02_a, lw=mylw, label = r'$E_2$')

plt.grid (b=True, which='both', linestyle=':', color='0.75', linewidth=0.5)
leg = plt.legend (loc="best")

#plt.xscale('symlog', linthreshx=1.0e-1)
#plt.xscale('log', linthreshx=1.0e-5)
plt.yscale('log', subsy = [2, 4, 6, 8])
#plt.yscale('symlog', linthreshy=1.0e-2, linscaley = 1.0, subsy = [2])
#plt.ylim (1.0e-10, 1.0e12)

plt.show ()
plt.clf ()

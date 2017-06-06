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
import numpy as np
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
cosmo = Nc.HICosmo.new_from_name (Nc.HICosmo, "NcHICosmoVexp")

if False: # Boa!
  cosmo.props.alphab   = +11.0e-2
  cosmo.props.sigmaphi = +6.8e0
  cosmo.props.dphi     = -5.0e-2
  cosmo.props.xb       = 4.8e37
  cosmo.props.OmegaL   = 1.0
  cosmo.props.Omegac   = 1.0
  cosmo.props.H0       = 67.8
else:
  cosmo.props.alphab   = +1.0e-20
  cosmo.props.sigmaphi = +8.0e-1
  cosmo.props.dphi     = -5.0e-1
  cosmo.props.xb       = 2.0e38
  cosmo.props.OmegaL   = 1.0
  cosmo.props.Omegac   = 1.0
  cosmo.props.H0       = 67.8

#
#  New cosmological distance objects optimizied to perform calculations
#  up to redshift 2.0.
#
dist = Nc.Distance.new (2.0)

tau_min = cosmo.tau_min ()
tau_max = cosmo.tau_max ()
k = 1.0e0
xb = cosmo.xbe ()

print "# ", tau_min, tau_max

tau_a = np.linspace (tau_min, tau_max, 100000)

alpha_a     = []
nu_a        = []
mnu_a       = []
dlnmnu_a    = []
dlnmnu_gw_a = []
Eaa0_a      = []

for tau in tau_a:
  nu, dlnmnu = cosmo.eval_system (tau, k)
  mnu        = cosmo.eval_mnu (tau, k)

  nu_gw, dlnmnu_gw = Nc.HIPertIGW.eval_system (cosmo, tau, k)  
  
  Eaa0 = k * tau / nu;
  
  alpha = 0
  if tau > 0.0:
    alpha = 0.5 * tau * tau / math.log (10.0)
  else:
    alpha = - 0.5 * tau * tau / math.log (10.0)

  alpha_a.append (alpha)
  nu_a.append (nu)
  mnu_a.append (mnu)
  dlnmnu_a.append (math.fabs (k * dlnmnu / nu))
  dlnmnu_gw_a.append (math.fabs (k * dlnmnu_gw / nu))
  Eaa0_a.append (math.fabs (Eaa0) * math.exp (-0.5 * tau * tau) * xb)

mylw = 1

print "# tau_qt_c % 21.15g % 21.15g" % (cosmo.tau_qt_c (), cosmo.tau_qt_e ())

#plt.plot (alpha_a, nu_a,     lw=mylw, label = r'$\nu$')
plt.plot (alpha_a, mnu_a,    lw=mylw, label = r'$m\nu$')
plt.plot (alpha_a, dlnmnu_a, lw=mylw, label = r'$\mathrm{d}\ln(m\nu)/\mathrm{d}\tau$')
plt.plot (alpha_a, dlnmnu_gw_a, lw=mylw, label = r'$\mathrm{d}\ln(m\nu)_{gw}/\mathrm{d}\tau$')
plt.plot (alpha_a, Eaa0_a, lw=mylw, label = r'$E$')
plt.axhline (1.0)

plt.grid (b=True, which='both', linestyle=':', color='0.75', linewidth=0.5)
leg = plt.legend (loc="upper left")

#plt.xscale('symlog', linthreshx=1.0e-5)
#plt.yscale('symlog', linthreshy=1.0e-20)
plt.yscale('log')

plt.show ()
plt.clf ()

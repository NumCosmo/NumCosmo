#!/usr/bin/python2

import math
import gi
gi.require_version('NumCosmo', '1.0')
gi.require_version('NumCosmoMath', '1.0')

from gi.repository import GObject
from gi.repository import NumCosmo as Nc
from gi.repository import NumCosmoMath as Ncm

import scipy.stats as ss
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from tqdm import tqdm
import time
import math
import sys

#
#  Initializing the library objects, this must be called before 
#  any other library function.
#
Ncm.cfg_init ()

#
#  New homogeneous and isotropic cosmological model NcHICosmoQGRW
#
cosmo = Nc.HICosmo.new_from_name (Nc.HICosmo, "NcHICosmoQGRW")

if len (sys.argv) != 2:
  print "%s w" % (sys.argv[0])
  sys.exit (0)

w      = float (sys.argv[1])
prec   = 1.0e-7

cosmo.props.w      = w
cosmo.props.Omegar = (1.0e-6) * 1.0
cosmo.props.Omegaw = (1.0 - 1.0e-6) * 1.0
cosmo.props.xb     = 1.e30

pert = Nc.HIPertTwoFluids.new ()

pert.props.reltol = prec
#pert.set_stiff_solver (True)

lnki  = math.log (1.0e-3)
lnkf  = math.log (1.0e3)
lnk_a = np.linspace (lnki, lnkf, 20)

ci = Ncm.Vector.new (8)

k_a       = []
Ps_zeta1  = []
Ps_S1     = []
Ps_zeta2  = []
Ps_S2     = []
Ps_Pzeta1 = []
Ps_PS1    = []
Ps_Pzeta2 = []
Ps_PS2    = []

out_file = f = open ('twofluids_spectrum_%e.dat' % (w), 'w')

start_alpha1 = 1.0e-10
start_alpha2 = 1.0e-14

for lnk in tqdm (lnk_a):
  k = math.exp (lnk)
  pert.set_mode_k (k)
  k_a.append (k)

  alphaf = cosmo.abs_alpha (1.0e20)

  #print "# Evolving mode %e from %f to %f" % (k, alphai, alphaf)

  alphai = -cosmo.abs_alpha (start_alpha1 * k**2)
  pert.get_init_cond_zetaS (cosmo, alphai, 1, 0.25 * math.pi, ci)
  pert.set_init_cond (cosmo, alphai, False, ci)
  
  print "# Mode 1 k % 21.15e, state module %f" % (k, pert.get_state_mod ())

  pert.evolve (cosmo, alphaf)
  v, alphac = pert.peek_state ()

  Delta_zeta1  = k**3 * math.hypot (v.get (Nc.HIPertITwoFluidsVars.ZETA_R),  v.get (Nc.HIPertITwoFluidsVars.ZETA_I))**2  / (2.0 * math.pi**2 * cosmo.RH_planck ()**2)
  Delta_S1     = k**3 * math.hypot (v.get (Nc.HIPertITwoFluidsVars.S_R),     v.get (Nc.HIPertITwoFluidsVars.S_I))**2     / (2.0 * math.pi**2 * cosmo.RH_planck ()**2)
  Delta_Pzeta1 = k**3 * math.hypot (v.get (Nc.HIPertITwoFluidsVars.PZETA_R), v.get (Nc.HIPertITwoFluidsVars.PZETA_I))**2 / (2.0 * math.pi**2 * cosmo.RH_planck ()**2)
  Delta_PS1    = k**3 * math.hypot (v.get (Nc.HIPertITwoFluidsVars.PS_R),    v.get (Nc.HIPertITwoFluidsVars.PS_I))**2    / (2.0 * math.pi**2 * cosmo.RH_planck ()**2)
    
  Ps_zeta1.append (Delta_zeta1)
  Ps_S1.append    (Delta_S1)
  Ps_Pzeta1.append (Delta_Pzeta1)
  Ps_PS1.append    (Delta_PS1)

  alphai = -cosmo.abs_alpha (start_alpha2 * k**2)
  pert.get_init_cond_zetaS (cosmo, alphai, 2, 0.25 * math.pi, ci)
  pert.set_init_cond (cosmo, alphai, False, ci)

  print "# Mode 2 k % 21.15e, state module %f" % (k, pert.get_state_mod ())

  pert.evolve (cosmo, alphaf)
  v, alphac = pert.peek_state ()

  Delta_zeta2  = k**3 * math.hypot (v.get (Nc.HIPertITwoFluidsVars.ZETA_R),  v.get (Nc.HIPertITwoFluidsVars.ZETA_I))**2  / (2.0 * math.pi**2 * cosmo.RH_planck ()**2)
  Delta_S2     = k**3 * math.hypot (v.get (Nc.HIPertITwoFluidsVars.S_R),     v.get (Nc.HIPertITwoFluidsVars.S_I))**2     / (2.0 * math.pi**2 * cosmo.RH_planck ()**2)
  Delta_Pzeta2 = k**3 * math.hypot (v.get (Nc.HIPertITwoFluidsVars.PZETA_R), v.get (Nc.HIPertITwoFluidsVars.PZETA_I))**2 / (2.0 * math.pi**2 * cosmo.RH_planck ()**2)
  Delta_PS2    = k**3 * math.hypot (v.get (Nc.HIPertITwoFluidsVars.PS_R),    v.get (Nc.HIPertITwoFluidsVars.PS_I))**2    / (2.0 * math.pi**2 * cosmo.RH_planck ()**2)

  Ps_zeta2.append (Delta_zeta2)
  Ps_S2.append    (Delta_S2)
  Ps_zeta2.append (Delta_Pzeta2)
  Ps_S2.append    (Delta_PS2)
  
  out_file.write ("% 20.15e % 20.15e % 20.15e % 20.15e % 20.15e % 20.15e % 20.15e % 20.15e % 20.15e\n" % (k, Delta_zeta1, Delta_zeta2, Delta_S1, Delta_S2, Delta_Pzeta1, Delta_Pzeta2, Delta_PS1, Delta_PS2))
  out_file.flush ()

out_file.close ()

plt.plot (k_a, Ps_zeta1, label = r'$P^1_\zeta$')
plt.plot (k_a, Ps_S1,    label = r'$P^1_S$')
plt.plot (k_a, Ps_zeta2, label = r'$P^2_\zeta$')
plt.plot (k_a, Ps_S2,    label = r'$P^2_S$')

plt.grid ()
plt.legend (loc="upper left")
plt.xscale('log')
plt.yscale('log')

plt.show ()
plt.clf ()

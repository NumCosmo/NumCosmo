#!/usr/bin/python2

import math
import gi
gi.require_version('NumCosmo', '1.0')
gi.require_version('NumCosmoMath', '1.0')

from gi.repository import GObject
from gi.repository import NumCosmo as Nc
from gi.repository import NumCosmoMath as Ncm

import numpy as np
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

w      = 1.0e-5
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

alpha_minf = -cosmo.abs_alpha (1.0e-40)
alpha_end  = cosmo.abs_alpha (1.0e25)
alphai     = -cosmo.abs_alpha (1.0e-40)
alphaf     = 0

print "# Inital time: %f %e" % (alphai, cosmo.x_alpha (alphai))

pert.set_stiff_solver (True)

x_a = []
gammabar11_a = []
gammabar22_a = []
gammabar12_a = []
taubar12_a   = []
nu1_a        = []
nu2_a        = []

for alpha in np.linspace (alphai, alphaf, 10000):
  eom = pert.eom (cosmo, alpha)
  x = cosmo.x_alpha (alpha)
  x_a.append (x)
  
  gammabar11_a.append (math.fabs (eom.gammabar11))
  gammabar22_a.append (math.fabs (eom.gammabar22))
  gammabar12_a.append (math.fabs (eom.gammabar12))
  taubar12_a.append (math.fabs (eom.taubar))
  nu1_a.append (eom.nu1)
  nu2_a.append (eom.nu2)
  
  #print ("% 20.15f % 20.15e % 20.15e % 20.15e % 20.15e % 20.15e % 20.15e" % (alpha, eom.gammabar[0], eom.gammabar[1], eom.gammabar[2], eom.taubar, eom.nu1, eom.nu2))

plt.plot (x_a, gammabar11_a, label = r'$\bar\gamma_{11}$')
plt.plot (x_a, gammabar22_a, label = r'$\bar\gamma_{22}$')
plt.plot (x_a, gammabar12_a, label = r'$\bar\gamma_{12}$')
plt.plot (x_a, taubar12_a, label = r'$\bar\tau_{12}$')

plt.grid ()
plt.legend (loc="upper left")
plt.xscale('log')
plt.yscale('log')
plt.yticks([1.0e-10, 1.0e-5, 1.0, 1.0e5, 1.0e10])
#plt.show ()
plt.clf ()

alphai     = -cosmo.abs_alpha (1.0e-13)

ci = pert.get_init_cond (cosmo, alphai, 1, 0.25 * math.pi)

pert.set_init_cond (cosmo, alphai, ci)

pert.evolve (cosmo, 1.0)


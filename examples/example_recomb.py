#!/usr/bin/env python

try:
  import gi
  gi.require_version('NumCosmo', '1.0')
  gi.require_version('NumCosmoMath', '1.0')
except:
  pass

from math import *
import matplotlib.pyplot as plt
from gi.repository import GObject
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
cosmo = Nc.HICosmo.new_from_name (Nc.HICosmo, "NcHICosmoDEXcdm")

#
#  New homogeneous and isotropic reionization object.
#
reion = Nc.HIReionCamb.new ()

#
# Adding submodels to the main cosmological model.
#
cosmo.add_submodel (reion)

#
#  New recombination object configured to calculate up to redshift 
#  10^9 and precision 10^-7.
#
recomb = Nc.Recomb.new_from_name ("NcRecombSeager{'prec':<1.0e-7>, 'zi':<1e9>}")

#
#  Setting values for the cosmological model, those not set stay in the
#  default values. Remeber to use the _orig_ version to set the original
#  parameters in case when a reparametrization is used.
#

#
# C-like
#
cosmo.orig_param_set (Nc.HICosmoDESParams.H0,        70.0)
cosmo.orig_param_set (Nc.HICosmoDESParams.OMEGA_C,   0.25)
cosmo.orig_param_set (Nc.HICosmoDESParams.OMEGA_X,   0.70)
cosmo.orig_param_set (Nc.HICosmoDESParams.T_GAMMA0,  2.72)
cosmo.orig_param_set (Nc.HICosmoDESParams.OMEGA_B,   0.05)
cosmo.orig_param_set (Nc.HICosmoDEXCDMSParams.W,    -1.00)

#
# OO-like
#
cosmo.props.H0      = 70.0
cosmo.props.Omegab  = 0.05
cosmo.props.Omegac  = 0.25
cosmo.props.Omegax  = 0.70
cosmo.props.Tgamma0 = 2.72
cosmo.props.w       = -1.0

#
#  Preparing recomb with cosmo.
#
recomb.prepare (cosmo)

#
#  Calculating Xe, equilibrium Xe, v_tau and its derivatives.
#
x_a      = []
x_dtau_a = []
Xe_a     = []
Xefi_a   = []
XHI_a    = []
XHII_a   = []
XHeI_a   = []
XHeII_a  = []
XHeIII_a = []

dtau_dlambda_a = []
x_E_a = []
x_Etaub_a = []
v_tau_a = []
dv_tau_dlambda_a = []
d2v_tau_dlambda2_a = []

for i in range (10000):
  alpha = -log (1.0e4) + (log (1.0e4) - log (1.0e2)) / 9999.0 * i
  Xe = recomb.Xe (cosmo, alpha)
  x = exp (-alpha)
  Xefi     = recomb.equilibrium_Xe (cosmo, x)

  XHI_i    = recomb.equilibrium_XHI (cosmo, x)
  XHII_i   = recomb.equilibrium_XHII (cosmo, x)
  XHeI_i   = recomb.equilibrium_XHeI (cosmo, x)
  XHeII_i  = recomb.equilibrium_XHeII (cosmo, x)
  XHeIII_i = recomb.equilibrium_XHeIII (cosmo, x)

  v_tau            = recomb.v_tau (cosmo, alpha)
  dv_tau_dlambda   = recomb.dv_tau_dlambda (cosmo, alpha)
  d2v_tau_dlambda2 = recomb.d2v_tau_dlambda2 (cosmo, alpha)

  x_a.append (x)
  Xe_a.append (Xe)
  Xefi_a.append (Xefi)

  XHI_a.append (XHI_i)
  XHII_a.append (XHII_i)
  XHeI_a.append (XHeI_i)
  XHeII_a.append (XHeII_i)
  XHeIII_a.append (XHeIII_i)

  v_tau_a.append (-v_tau)
  dv_tau_dlambda_a.append (-dv_tau_dlambda / 10.0)
  d2v_tau_dlambda2_a.append (-d2v_tau_dlambda2 / 200.0)

for i in range (10000):
  alpha = -log (1.0e10) + (log (1.0e10) - log (1.0e2)) / 9999.0 * i
  x = exp (-alpha)
  dtau_dlambda = abs (recomb.dtau_dlambda (cosmo, alpha))
  x_E = x / sqrt (cosmo.E2 (x))
  x_E_taub = x_E / dtau_dlambda

  x_dtau_a.append (x)
  x_E_a.append (x_E)
  x_Etaub_a.append (x_E_taub)
  dtau_dlambda_a.append (dtau_dlambda)

#
#  Ploting ionization history.
#

plt.title ("Ionization History")
plt.xscale('log')
plt.yscale('log')
plt.plot (x_a, Xe_a, 'r', label="Recombination")
plt.plot (x_a, Xefi_a, 'b--', label="Equilibrium")
plt.xlabel('$1+z$')
plt.ylabel('$X_{e^-}$')
plt.legend(loc=3)
plt.xlim(1e4, 1e2)
plt.ylim(1e-4,2)

plt.savefig ("recomb_Xe.svg")


plt.clf ()

#
#  Ploting all components history.
#

plt.title ("Fractions Equilibrium History")
plt.xscale('log')
# plt.yscale('log')
plt.plot (x_a, XHI_a, 'r', label=r'$X_\mathrm{HI}$')
plt.plot (x_a, XHII_a, 'g', label=r'$X_\mathrm{HII}$')
plt.plot (x_a, XHeI_a, 'b', label=r'$X_\mathrm{HeI}$')
plt.plot (x_a, XHeII_a, 'y', label=r'$X_\mathrm{HeII}$')
plt.plot (x_a, XHeIII_a, 'm', label=r'$X_\mathrm{HeIII}$')
plt.xlabel(r'$1+z$')
plt.ylabel(r'$X_*$')
plt.legend(loc=6)

plt.xlim(1e4,1e2)

plt.savefig ("recomb_eq_fractions.svg")

plt.clf ()

#
#  Ploting visibility function and derivatives.
#

(lambdam, lambdal, lambdau) = recomb.v_tau_lambda_features (cosmo, 2.0 * log (10.0))

plt.title ("Visibility Function and Derivatives")
plt.xscale ('log')
plt.xlabel ('$1+z$')
plt.plot (x_a, v_tau_a, 'r', label=r'$v_\tau$')
plt.plot (x_a, dv_tau_dlambda_a, 'b-', label=r'$\frac{1}{10}\frac{dv_\tau}{d\lambda}$')
plt.plot (x_a, d2v_tau_dlambda2_a, 'g--', label=r'$\frac{1}{200}\frac{d^2v_\tau}{d\lambda^2}$')
plt.legend()
plt.legend(loc=3)

#
#  Annotating max and width.
#

v_tau_max = -recomb.v_tau (cosmo, lambdam)

plt.annotate (r'$v_\tau^{max}$, $z=%5.2f$' % (exp(-lambdam)-1), xy=(exp(-lambdam), v_tau_max),  xycoords='data',
              xytext=(0.1, 0.95), textcoords='axes fraction',
              arrowprops=dict(facecolor='black', shrink=0.05))

v_tau_u = -recomb.v_tau (cosmo, lambdau)

plt.annotate (r'$v_\tau=10^{-2}v_\tau^{max}$, $z=%5.2f$' % (exp(-lambdau)-1), xy=(exp(-lambdau), v_tau_u),  xycoords='data',
              xytext=(0.02, 0.75), textcoords='axes fraction',
              arrowprops=dict(facecolor='black', shrink=0.05))

v_tau_l = -recomb.v_tau (cosmo, lambdal)

plt.annotate (r'$v_\tau=10^{-2}v_\tau^{max}$, $z=%5.2f$' % (exp(-lambdal)-1), xy=(exp(-lambdal), v_tau_l),  xycoords='data',
              xytext=(0.60, 0.25), textcoords='axes fraction',
              arrowprops=dict(facecolor='black', shrink=0.05))

#
#  Annotating value of lambda (tau == 1).
#

lambdastar = recomb.get_tau_drag_lambda (cosmo)

plt.annotate (r'$z^\star=%5.2f$' % (exp(-lambdastar)-1), 
              xy=(0.65, 0.95), xycoords='axes fraction')

plt.savefig ("recomb_v_tau.svg")

plt.clf ()

#
#  Ploting dtau_dlambda
#

plt.title (r'$d\tau/d\lambda$')
plt.xscale('log')
plt.yscale('log')
plt.plot (x_dtau_a, dtau_dlambda_a, 'r', label=r'$\vert\mathrm{d}\tau/\mathrm{d}\lambda\vert$')
plt.plot (x_dtau_a, x_E_a, 'b', label=r'$(1+z)/E(z)$')
plt.plot (x_dtau_a, (x_Etaub_a), 'g', label=r'$(1+z)/(E(z)\bar{\tau})$')

plt.xlabel ('$1 + z$')
plt.legend()
plt.legend(loc=2)

plt.savefig ("recomb_dtau_dlambda.svg")


#!/usr/bin/python2

from math import *
from gi.repository import GObject
import matplotlib.pyplot as plt
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
cosmo.orig_param_set (Nc.HICosmoDEParams.H0,        70.0)
cosmo.orig_param_set (Nc.HICosmoDEParams.OMEGA_C,   0.25)
cosmo.orig_param_set (Nc.HICosmoDEParams.OMEGA_X,   0.70)
cosmo.orig_param_set (Nc.HICosmoDEParams.T_GAMMA0,  2.72)
cosmo.orig_param_set (Nc.HICosmoDEParams.OMEGA_B,   0.05)
cosmo.orig_param_set (Nc.HICosmoDEParams.SPECINDEX, 1.00)
cosmo.orig_param_set (Nc.HICosmoDEParams.SIGMA8,    0.90)
cosmo.orig_param_set (Nc.HICosmoDEXCDMParams.W,    -1.00)

#
# OO-like
#
cosmo.props.H0      = 70.0
cosmo.props.Omegab  = 0.05
cosmo.props.Omegac  = 0.25
cosmo.props.Omegax  = 0.70
cosmo.props.Tgamma0 = 2.72
cosmo.props.ns      = 1.0
cosmo.props.sigma8  = 0.9
cosmo.props.w       = -1.0

#
#  Preparing recomb with cosmo.
#
recomb.prepare (cosmo)

#
#  Calculating Xe, equilibrium Xe, v_tau and its derivatives.
#
x_a = []
Xe_a = []
Xefi_a = []
v_tau_a = []
dv_tau_dlambda_a = []
d2v_tau_dlambda2_a = []

for i in range (10000):
  alpha = -log (10000.0) + (log (10000.0) - log (100.0)) / 9999.0 * i
  Xe = recomb.Xe (cosmo, alpha)
  x = exp (-alpha)
  Xefi = recomb.equilibrium_Xe (cosmo, x)
  v_tau = recomb.v_tau (cosmo, alpha)
  dv_tau_dlambda = recomb.dv_tau_dlambda (cosmo, alpha)
  d2v_tau_dlambda2 = recomb.d2v_tau_dlambda2 (cosmo, alpha)
  
  x_a.append (x)
  Xe_a.append (Xe)
  Xefi_a.append (Xefi)
  v_tau_a.append (-v_tau)
  dv_tau_dlambda_a.append (-dv_tau_dlambda / 10.0)
  d2v_tau_dlambda2_a.append (-d2v_tau_dlambda2 / 200.0)

#
#  Ploting ionization history.
#

plt.title ("Ionization History")
plt.xscale('log')
plt.plot (x_a, Xe_a, 'r', label="Recombination")
plt.plot (x_a, Xefi_a, 'b--', label="Equilibrium")
plt.xlabel('$x$')
plt.ylabel('$X_{e^-}$')
plt.legend(loc=2)

plt.savefig ("recomb_Xe.png")

plt.clf ()

(lambdam, lambdal, lambdau) = recomb.v_tau_lambda_features (cosmo, 2.0 * log (10.0))

#
#  Ploting visibility function and derivatives.
#

plt.title ("Visibility Function and Derivatives")
plt.xscale('log')
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
              xytext=(0.65, 0.25), textcoords='axes fraction',
              arrowprops=dict(facecolor='black', shrink=0.05))

#
#  Annotating value of zstar.
#

lambdastar = recomb.tau_zstar (cosmo)

plt.annotate (r'$z^\star=%5.2f$' % (exp(-lambdastar)-1), 
              xy=(0.65, 0.95), xycoords='axes fraction')

plt.savefig ("recomb_v_tau.png")


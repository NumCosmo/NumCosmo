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
#  New cosmological distance objects optimizied to perform calculations
#  up to redshift 2.0.
#
dist = Nc.Distance.new (2.0)

#
# New windown function 'NcWindowTophat'
#
wp =  Nc.Window.new_from_name ("NcWindowTophat")

#
# New transfer function 'NcTransferFuncEH' using the Einsenstein, Hu
# fitting formula.
#
tf = Nc.TransferFunc.new_from_name ("NcTransferFuncEH")

#
# New matter variance object using FFT method for internal calculations and
# the window and transfer functions defined above.
#
vp = Nc.MatterVar.new (Nc.MatterVarStrategy.FFT, wp, tf)

#
# New growth function
#
gf = Nc.GrowthFunc.new ()

#
# New multiplicity function 'NcMultiplicityFuncTinkerMean'
#
mulf = Nc.MultiplicityFunc.new_from_name ("NcMultiplicityFuncTinkerMean")

#
# New mass function object using the objects defined above.
#
mf = Nc.MassFunction.new (dist, vp, gf, mulf)

#
#  Setting values for the cosmological model, those not set stay in the
#  default values. Remember to use the _orig_ version to set the original
#  parameters when a reparametrization is used.
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
#  Printing the parameters used.
#
print "# Model parameters: ", 
cosmo.params_log_all ()

#
# Number of points to build the plots
#
np = 2000
divfac = 1.0 / (np - 1.0)

#
#  Calculating growth and its derivative in the [0, 2] redshift
#  range.
#
za = []
Da = []
dDa = []

gf.prepare (cosmo)

for i in range (0, np):
  z = 2.0 * divfac * i
  D = gf.eval (cosmo, z)
  dD = gf.eval_deriv (cosmo, z)
  za.append (z)
  Da.append (D)
  dDa.append (dD)
  print z, D

#
#  Ploting growth function.
#

plt.title ("Growth Function")
plt.plot (za, Da, 'r', label="D")
plt.plot (za, dDa, 'b--', label="dD/dz")
plt.xlabel('$z$')
plt.legend(loc=2)

plt.savefig ("growth_func.png")
plt.clf ()

#
# Calculating the transfer function and the matter power spectrum in the
# kh (in unities of h/Mpc) interval [1e-3, 1e3]
#

kha = []
Ta = []
Pma = []

tf.prepare (cosmo)

for i in range (0, np):
  lnkh = log (1e-4) +  log (1e7) * divfac * i
  kh = exp (lnkh)
  T = tf.eval (cosmo, kh)
  Pm = 1.0e3 / 7.0 * tf.matter_powerspectrum (cosmo, kh)
  kha.append (kh)
  Ta.append (T)
  Pma.append (Pm)

#
#  Ploting transfer and matter power spectrum
#

plt.title ("Transfer Function and Matter Power Spectrum")
plt.xscale('log')
plt.plot (kha, Ta, 'r', label="T(kh)")
plt.plot (kha, Pma, 'b--', label="P_m(kh)")
plt.xlabel('$k$')
plt.legend(loc=1)

plt.savefig ("transfer_func.png")
plt.clf ()

#
# Calculating the variance filtered with the tophat windown function using
# scales R in the interval [5, 50] at redshift 0.3.
# First calculates the growth function at z = 0.3 and then the spectrum
# amplitude from the sigma8 parameter.
#

vp.prepare (cosmo)

Dz = gf.eval (cosmo, 0.3)
A  = vp.sigma8_sqrtvar0 (cosmo)
prefact = A * A * Dz * Dz

Ra = []
sigma2a = []
dlnsigma2a = []

for i in range (0, np):
  lnR = log (5.0) +  log (10.0) * divfac * i
  R = exp (lnR)
  sigma2 = prefact * vp.var0 (cosmo, lnR)
  dlnsigma2 = vp.dlnvar0_dlnR (cosmo, lnR) 
  Ra.append (R)
  sigma2a.append (sigma2)
  dlnsigma2a.append (dlnsigma2)

#
#  Ploting filtered matter variance
#

plt.title ("Variance and Variance Derivative")
plt.plot (Ra, sigma2a, 'r', label='$\sigma^2(\ln(R))$')
plt.plot (Ra, dlnsigma2a, 'b--', label='$d\ln(\sigma^2)/d\ln(R)$')
plt.xlabel('$R$')
plt.legend(loc=1)

plt.savefig ("matter_var.png")
plt.clf ()

#
# Calculating the mass function integrated in the mass interval [1e14, 1e16]
# for the redhshifts in the interval [0, 2.0] and area 200 squared degree.
#

mf.set_area_sd (200.0)
mf.set_eval_limits (cosmo, log (1e14), log(1e16), 0.0, 2.0)

dndza = []

for i in range (0, np):
  dndz = mf.dn_dz (cosmo, log(1e14), log(1e16), za[i], True)
  dndza.append (dndz)

#
#  Ploting the mass function
#

plt.title ("Mass Function")
plt.plot (za, dndza, 'r', label='$dn/dz$')
plt.xlabel('$z$')
plt.legend(loc=1)

plt.savefig ("mass_function.png")
plt.clf ()

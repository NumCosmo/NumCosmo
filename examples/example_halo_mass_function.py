#!/usr/bin/python2

try:
  import gi
  gi.require_version('NumCosmo', '1.0')
  gi.require_version('NumCosmoMath', '1.0')
except:
  pass

import numpy
import math 
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
#  New homogeneous and isotropic primordial object.
#
prim = Nc.HIPrimPowerLaw.new () 

#
# Adding submodels to the main cosmological model.
#
cosmo.add_submodel (reion)
cosmo.add_submodel (prim)

#
#  New cosmological distance objects optimizied to perform calculations
#  up to redshift 2.0.
#
dist = Nc.Distance.new (2.0)

#
# New transfer function 'NcTransferFuncEH' using the Einsenstein, Hu
# fitting formula.
#
tf = Nc.TransferFunc.new_from_name ("NcTransferFuncEH")

#
# New linear matter power spectrum object based of the EH transfer function.
# 
psml = Nc.PowspecMLTransfer.new (tf)
psml.require_kmin (1.0e-3)
psml.require_kmax (1.0e3)

#
# Apply a tophat filter to the psml object, set best output interval.
#
psf = Ncm.PowspecFilter.new (psml, Ncm.PowspecFilterType.TOPHAT)
psf.set_best_lnr0 ()

#
# New multiplicity function 'NcMultiplicityFuncTinkerMean'
#
mulf = Nc.MultiplicityFunc.new_from_name ("NcMultiplicityFuncTinkerMean")

#
# New mass function object using the objects defined above.
#
mf = Nc.HaloMassFunction.new (dist, psf, mulf)

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

gf = psml.peek_gf ()
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
plt.plot (za, Da, 'r', label="D(z)")
plt.plot (za, dDa, 'b--', label="dD/dz")
plt.xlabel('$z$')
plt.legend(loc=1)
#plt.yscale ('log')

plt.savefig ("hmf_growth_func.svg")
plt.clf ()

#
# Calculating the transfer function and the matter power spectrum in the
# kh (in unities of h/Mpc) interval [1e-3, 1e3]
#

ka = []
Ta = []
Pma = []

psml.prepare (cosmo)

for i in range (0, np):
  lnk = math.log (1e-3) +  math.log (1e6) * divfac * i
  k   = math.exp (lnk)
  T   = tf.eval (cosmo, k)
  Pm  = psml.eval (cosmo, 0.0, k)
  ka.append (k)
  Ta.append (T)
  Pma.append (Pm)

#
#  Ploting transfer and matter power spectrum
#

plt.title ("Transfer Function and Linear Matter Power Spectrum")
plt.yscale('log')
plt.xscale('log')
plt.plot (ka, Ta, 'r', label="$T(k)$")
plt.plot (ka, Pma, 'b--', label="$P_m(k)$")
plt.xlabel('$k \; [\mathrm{Mpc}^{-1}]$')
plt.legend(loc=1)

plt.savefig ("hmf_transfer_func.svg")
plt.clf ()

#
# Calculating the variance filtered with the tophat windown function using
# scales R in the interval [5, 50] at redshift 0.3.
# First calculates the growth function at z = 0.3 and then the spectrum
# amplitude from the sigma8 parameter.
#
psf.prepare (cosmo)

Ra = []
sigma2a = []
dlnsigma2a = []

for i in range (0, np):
  lnR       = math.log (5.0) +  math.log (10.0) * divfac * i
  R         = math.exp (lnR)
  sigma2    = psf.eval_var_lnr (0.0, lnR)
  dlnsigma2 = psf.eval_dlnvar_dlnr (0.0, lnR) 
  Ra.append (R)
  sigma2a.append (sigma2)
  dlnsigma2a.append (dlnsigma2)

#
#  Ploting filtered matter variance
#
plt.title ("Variance and its derivative")
plt.plot (Ra, sigma2a, 'r', label='$\sigma^2(\ln(R))$')
plt.plot (Ra, dlnsigma2a, 'b--', label='$d\ln(\sigma^2)/d\ln(R)$')
plt.xlabel('$R \; [\mathrm{Mpc}]$')
plt.legend(loc=1)

plt.savefig ("hmf_matter_variance.svg")
plt.clf ()

#
# Calculating the halo mass function at z = 0.7, and integrated in the mass interval [1e14, 1e16]
# for the redhshifts in the interval [0, 2.0] and area 200 squared degrees.
#

mf.set_area_sd (200.0)
mf.set_eval_limits (cosmo, math.log (1e14), math.log(1e16), 0.0, 2.0)
mf.prepare (cosmo)

lnMa = []
dn_dlnMdza = []
dndza = []

for i in range (0, np):
  lnM = math.log (1e14) + math.log (1e2) * divfac * i
  dn_dlnMdz = mf.dn_dlnM (cosmo, lnM, 0.7)
  dndz = mf.dn_dz (cosmo, math.log (1.0e14), math.log (1.0e16), za[i], True)
  lnMa.append (lnM)
  dn_dlnMdza.append (dn_dlnMdz)
  dndza.append (dndz)

#
#  Ploting the mass function
#

plt.title ("Halo Mass Function")
plt.plot (lnMa, dn_dlnMdza, 'b', label='$z = 0.7$')
plt.yscale('log')
#plt.xscale('log')
plt.xlabel('$\ln M \; [\mathrm{M}_\odot]$')
plt.ylabel('${\mathrm{d}^2n}/{\mathrm{d}z \mathrm{d}\ln M}$')
plt.legend(loc=1)

plt.savefig ("hmf_mass_function.svg")
plt.clf ()

#
#  Ploting the mass function
#

plt.title ("Number of halos per redshift")
plt.plot (za, dndza, 'b', label='$M \in [10^{14}, 10^{16}], \; A = 200 \; [\mathrm{deg}^2]$')
plt.xlabel('$z$')
plt.ylabel('${\mathrm{d}N}/{\mathrm{d}z}$')
plt.legend(loc=1)

plt.savefig ("hmf_halos_redshift.svg")
plt.clf ()

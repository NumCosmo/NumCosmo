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
mulf = Nc.MultiplicityFuncTinker.new ()
mulf.set_mdef (Nc.MultiplicityFuncMassDef.MEAN)
mulf.set_Delta (200.0)

#
# New mass function object using the objects defined above.
#
mf = Nc.HaloMassFunction.new (dist, psf, mulf)

#
# New Cluster Mass object using Log normal distribution
#
lnMobs_min = log (1.0e14)
lnMobs_max = log (1.0e16)
cluster_m = Nc.ClusterMass.new_from_name ("NcClusterMassLnnormal{'lnMobs-min':<%20.15e>, 'lnMobs-max':<%20.15e>}" % (lnMobs_min, lnMobs_max))

#
# New Cluster Redshift object using a global gaussian distribution
#
z_min = 0.0
z_max = 0.7
cluster_z = Nc.ClusterRedshift.new_from_name ("NcClusterPhotozGaussGlobal{'pz-min':<%20.15e>, 'pz-max':<%20.15e>, 'z-bias':<0.0>, 'sigma0':<0.03>}" % (z_min, z_max))

#
# New Cluster abundance object that uses all objects above
#
cad = Nc.ClusterAbundance.new (mf, None)

#
# New NcmData object for number count calculations
#
ncdata = Nc.DataClusterNCount.new (cad)

#
#  Creating a new Modelset and set cosmo as the HICosmo model to be used
#  and cluster_m as the distribution of the mass-observable relation
#
mset = Ncm.MSet.new_array ([cosmo, cluster_z, cluster_m])

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
#  Setting values for the mass distribution model
#
cluster_m.props.bias       = 0.0
cluster_m.props.sigma      = 0.2

#
#  Printing the parameters used.
#
mset.pretty_log ()

#
# Creates a new random number generator from a pool named "example_ca_sampling"
# it implicitly creates this pool.
#
rng = Ncm.RNG.pool_get ("example_ca_sampling");

#
# Since ncdata is currently empty, run init_from_sampling
# using the objects above and an survey area of 300degsq^2
#
ncdata.init_from_sampling (mset, 270 * (pi / 180.0)**2, rng)

#
# Save to a fits file
#
ncdata.catalog_save ("ca_data.fits", True)

#
# Generate another sample by resampling from mset
#
ncdata.resample (mset, rng)

#
# Checking if it has the mass truth table, if so gets it
#
has_lnM_true = ncdata.has_lnM_true ()
print ("# Has mass truth table = ", has_lnM_true)
lnM_true = None
if ncdata.has_lnM_true ():
  lnM_true = ncdata.get_lnM_true ()

#
# Checking if it has the redshift truth table, if so gets it
#
has_z_true = ncdata.has_z_true ()
print ("# Has redshift truth table = ", has_z_true)
z_true = None
if ncdata.has_z_true ():
  z_true = ncdata.get_z_true ()

#
# Gets the mass observables and its parameters
#
lnM_obs = ncdata.get_lnM_obs ()
lnM_obs_params = ncdata.get_lnM_obs_params ()

#
# Gets the redshift observables
#
z_obs = ncdata.get_z_obs ()
z_obs_params = ncdata.get_z_obs_params ()

#
# Print everything in a file
#
nobjects = ncdata.get_len ()
print ("# There are ", nobjects, " objects in the catalog (%d, %d)" % (lnM_obs.col_len (), z_obs.col_len ()))

f = open ('ca_data.dat', 'w')

for i in range (0, nobjects):
  row = "%d " % (i)

  if has_lnM_true: 
    row += "%f " % (lnM_true.get (i))

  if has_z_true: 
    row += "%f " % (z_true.get (i))

  for j in range (0, lnM_obs.row_len ()):
    row += "%f " % (lnM_obs.get (i, j))

  for j in range (0, z_obs.row_len ()):
    row += "%f " % (z_obs.get (i, j))

  if lnM_obs_params:
    for j in range (0, lnM_obs_params.row_len ()):
      row += "%f " % (lnM_obs_params.get (i, j))
  
  if z_obs_params:
    for j in range (0, z_obs_params.row_len ()):
      row += "%f " % (z_obs_params.get (i, j))

  f.write (row)
  f.write ("\n")
  
f.close ()

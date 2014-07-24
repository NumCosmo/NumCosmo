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
# New Cluster Mass object using Log normal distribution
#
cluster_m = Nc.ClusterMass.new_from_name ("NcClusterMassLnnormal{'lnMobs-min':<31.5430441213567>, 'lnMobs-max':<33.6224856630365>}")

#
# New Cluster Redshift object using a global gaussian distribution
#
cluster_z = Nc.ClusterRedshift.new_from_name ("NcClusterPhotozGaussGlobal{'pz-min':<0.0>, 'pz-max':<1.0>, 'z-bias':<0.0>, 'sigma0':<0.03>}")

#
# New Cluster abundance object that uses all objects above
#
cad = Nc.ClusterAbundance.new (mf, None, cluster_z, cluster_m)

#
# New NcmData object for number count calculations
#
ncdata = Nc.DataClusterNCount.new (cad)

#
#  Creating a new Modelset and set cosmo as the HICosmo model to be used
#  and cluster_m as the distribution of the mass-observable relation
#
mset = Ncm.MSet ()
mset.set (cosmo)
mset.set (cluster_m)

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
ncdata.init_from_sampling (mset, cluster_z, cluster_m, 300 * (pi / 180.0)**2, rng)

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
print "# Has mass truth table = ", has_lnM_true
lnM_true = None
if ncdata.has_lnM_true ():
  lnM_true = ncdata.get_lnM_true ()

#
# Checking if it has the redshift truth table, if so gets it
#
has_z_true = ncdata.has_z_true ()
print "# Has redshift truth table = ", has_z_true
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
print "# There are ", nobjects, " objects in the catalog (%d, %d)" % (lnM_obs.col_len (), z_obs.col_len ())

f = open ('ca_data.txt', 'w')

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

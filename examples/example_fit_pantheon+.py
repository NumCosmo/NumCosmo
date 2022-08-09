#!/usr/bin/env python

try:
  import gi
  gi.require_version('NumCosmo', '1.0')
  gi.require_version('NumCosmoMath', '1.0')
except:
  pass

from math import *
import os.path
import sys
import matplotlib.pyplot as plt
from gi.repository import GLib
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
cosmo.omega_x2omega_k ()

#
#  Setting values for the cosmological model, those not set stay in the
#  default values. Remeber to use the _orig_ version to set the original
#  parameters in case when a reparametrization is used.
#

#
# OO-like
#
cosmo.param_set_by_name ("H0",     70.0)
cosmo.param_set_by_name ("Omegab",  0.05)
cosmo.param_set_by_name ("Omegac",  0.25)
cosmo.param_set_by_name ("Omegak",  0.0)
cosmo.param_set_by_name ("w",      -1.0)

#
#  Setting parameters Omega_c and w to be fitted.
#

cosmo.props.H0_fit = True
cosmo.props.Omegac_fit = True
cosmo.props.w_fit = False

#
#  Creating a new Distance object optimized to redshift 3.
#
dist = Nc.Distance (zf = 3.0)


#
# Getting SNIa sample
#
if len (sys.argv) <= 1:
    print (f"Usage {sys.argv[0]} sample_id")
    exit (-1)

snia_id = sys.argv[1]

try:
    snia_idval = Nc.DataSNIACov.get_catalog_id (snia_id)
except GLib.Error as err:
    if err.matches (Nc.data_snia_cov_error_quark (), Nc.DataSNIACovError.ID_NOT_FOUND):
        Ncm.cfg_enum_print_all (Nc.DataSNIAId, "Invalid sample id `{snia_id}', the valid ids are:")
        exit (-1)

#
# SNIa cov model
#
snia_model = Nc.SNIADistCov.new_by_id (dist, snia_idval)

#
#  Creating a new Modelset and set cosmo as the HICosmo model to be used.
#
mset = Ncm.MSet ()
mset.set (cosmo)
mset.set (snia_model)


#
#  Creating a new Data object from distance modulus catalogs.
#
snia = Nc.DataSNIACov.new_from_cat_id (snia_idval, False)
try:
    snia0 = snia.apply_filter_sh0es_z (0.01, True)
    if snia0:
        snia = snia0
except:
    pass

#
#  Creating a new Dataset and add snia to it.
#
dset = Ncm.Dataset ()
dset.append_data (snia)

#
#  Creating a Likelihood from the Dataset.
#
lh = Ncm.Likelihood (dataset = dset)

#
#  Creating a Fit object of type NLOPT using the fitting algorithm ln-neldermead to
#  fit the Modelset mset using the Likelihood lh and using a numerical differentiation
#  algorithm (NUMDIFF_FORWARD) to obtain the gradient (if needed).
#
fit = Ncm.Fit.new (Ncm.FitType.NLOPT, "ln-neldermead", lh, mset, Ncm.FitGradType.NUMDIFF_FORWARD)

#
#  Running the fitter printing messages.
#
fit.run_restart (Ncm.FitRunMsgs.SIMPLE, 1.0e-3, 0.0, None, None)

#
#  Printing fitting informations.
#
fit.log_info ()

#
#  Calculating the parameters covariance using numerical differentiation.
#
fit.numdiff_m2lnL_covar ()

#
#  Printing the covariance matrix.
# 
fit.log_covar ()

#
# Additional functions
#
mfunc_oa = Ncm.ObjArray.new ()

mfunc_Omegam   = Ncm.MSetFuncList.new ("NcHICosmo:Omega_m0", None)
mfunc_oa.add (mfunc_Omegam)

#
# Setting single thread calculation.
#
Ncm.func_eval_set_max_threads (4)
Ncm.func_eval_log_pool_stats ()

#
# New Gaussian prior to provide the initial points for the chain.
# It was created with size 0 (number of parameters), but once 
# initialized with mset the correct size is assigned. 
#
# The initial sampler will use a diagonal covariance with the
# diagonal terms being the parameters scale set by each model.
#
init_sampler = Ncm.MSetTransKernGauss.new (0)
init_sampler.set_mset (mset)
init_sampler.set_prior_from_mset ()
init_sampler.set_cov_from_rescale (1.0)

#
# Creates the ESMCMC walker object, this object is responsible
# for moving the walkers in each interation, the stretch move
# is affine invariant and therefore gives good results even for
# very correlated parametric space.
# 
nwalkers = 300
walker = Ncm.FitESMCMCWalkerAPES.new (nwalkers, mset.fparams_len ())

fitscat = "example_fit_{snia_id}_nwalkers%d.fits" % (nwalkers)

if os.path.exists (fitscat):
    lmcat = Ncm.MSetCatalog.new_from_file_ro (fitscat, 0)
    mcat_len = lmcat.len ()
    
    if mcat_len > nwalkers * 50:
    
        last_e = [lmcat.peek_row (mcat_len - nwalkers + i) for i in range (nwalkers)]
        nadd_vals = lmcat.nadd_vals ()
        ncols = lmcat.ncols ()
        nvar = ncols - nadd_vals
    
        k = Ncm.StatsDistKernelST.new (nvar, 1.0)
        sd = Ncm.StatsDistVKDE.new (k, Ncm.StatsDistCV.SPLIT)
        sd.reset ()
        m2lnL = []
        for row in last_e:
            m2lnL.append (row.get (0))
            sd.add_obs (row.get_subvector (nadd_vals, nvar))

        m2lnL_v = Ncm.Vector.new_array (m2lnL)
        sd.prepare_interp (m2lnL_v)
        ovs = sd.get_over_smooth ()
        walker.set_over_smooth (ovs)
        print (f"Setting over smooth to {ovs}")
        del lmcat


#
# Initialize the ESMCMC object using the objects above. It will
# use 50 walkers, i.e., each point in the MCMC chain contains
# 50 points in the parametric space. Each step uses the last point
# in the chain (the last 50 parametric points) to calculate the
# proposal points.
#
esmcmc  = Ncm.FitESMCMC.new_funcs_array (fit, nwalkers, init_sampler, walker, Ncm.FitRunMsgs.SIMPLE, mfunc_oa)

#
# These methods enable the auto-trim options on ESMCMC. This option 
# makes the sampler check the chains' health and trim any unnecessary 
# burn-in part. We set the number of divisions to 100 so we test the
# chains in blocks of n/100. The last method asserts that each 2min
# the catalog will be checked.
#
#esmcmc.set_auto_trim (True)
#esmcmc.set_auto_trim_div (100)
#esmcmc.set_max_runs_time (2.0 * 60.0)
esmcmc.set_nthreads (4)
esmcmc.set_data_file (fitscat)

#
# Running the esmcmc, it will first calculate 1000 points, after that
# it will estimate the error in the parameters mean. Using the current
# errors the algorithm tries to calculated how many extra steps are 
# necessary to obtain the required error `10^-3' in every parameters,
# and it will run such extra steps. It will repeat this procedure
# until it attains the required error in every parameter.
# 
#
esmcmc.start_run ()
esmcmc.run_lre (500, 1.0e-3)
esmcmc.end_run ()

#
# Calculates the parameter means and covariance and set it into 
# the fit object and then print.
# 
esmcmc.mean_covar ()
fit.log_covar ()

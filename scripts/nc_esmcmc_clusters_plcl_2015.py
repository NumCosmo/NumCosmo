#!/usr/bin/python2
   
from math import *
from gi.repository import GObject
import matplotlib.pyplot as plt
from gi.repository import NumCosmo as Nc
from gi.repository import NumCosmoMath as Ncm

    """
     Perform a Emsemble Sampler Markov Chain Monte Carlo (ESMCMC) analysis of the parameters of the Planck-CLASH 
     mass-observables distribution and the selection function.
     Probe: cluster pseudo counts
    """

Ncm.cfg_init ()

NT = 3           # Number of threads
NClusters = 21   # Number of clusters
NChains = 50     # Number of chains

cosmo = Nc.HICosmo.new_from_name (Nc.HICosmo, "NcHICosmoDEXcdm")
dist = Nc.Distance.new (4.0)
wp = Nc.Window.new_from_name ("NcWindowTophat")
tf = Nc.TransferFunc.new_from_name ("NcTransferFuncEH")
vp = Nc.MatterVar.new (Nc.MatterVarStrategy.FFT, wp, tf)
gf = Nc.GrowthFunc.new ()
mulf = Nc.MultiplicityFunc.new_from_name ("NcMultiplicityFuncTinkerCrit{'Delta':<500.0>}")
mf = Nc.MassFunction.new (dist, vp, gf, mulf)
clusterm = Nc.ClusterMass.new_from_name ("NcClusterMassPlCL{}")
cpc = Nc.ClusterPseudoCounts.new (mf, NClusters)

# Planck + BAO + WP + High l (page 31, table 4 1502.01589)
cosmo.props.H0      = 67.74
cosmo.props.Omegab  = 0.049
cosmo.props.Omegac  = 0.259
cosmo.props.Omegax  = 0.691
cosmo.props.Tgamma0 = 2.72
cosmo.props.ns      = 0.967
cosmo.props.sigma8  = 0.816
cosmo.props.w       = -1.0
cosmo.omega_x2omega_k ()
cosmo.param_set_by_name ('Omegak', 0.0)

vp.prepare (cosmo)

# Building Model Set
mset = Ncm.MSet.empty_new ()
mset.set (cosmo)
mset.set (clusterm)
mset.set (cpc)

# Defining which parameters to fit
clusterm.param_set_ftype (0, Ncm.ParamType.FREE)
clusterm.param_set_ftype (1, Ncm.ParamType.FREE)
clusterm.param_set_ftype (2, Ncm.ParamType.FREE)
clusterm.param_set_ftype (3, Ncm.ParamType.FIXED)
clusterm.param_set_ftype (4, Ncm.ParamType.FIXED)
clusterm.param_set_ftype (5, Ncm.ParamType.FREE)
clusterm.param_set_ftype (6, Ncm.ParamType.FREE)

clusterm.param_set_by_name ('Asz', 1.499)
clusterm.param_set_by_name ('Bsz', 0.274)
clusterm.param_set_by_name ('sigma_sz', 0.132)
clusterm.param_set_by_name ('Al', 1.0)
clusterm.param_set_by_name ('Bl', 0.0)
clusterm.param_set_by_name ('sigma_l', 0.165)
clusterm.param_set_by_name ('cor', -0.997)

cpc.param_set_ftype (0, Ncm.ParamType.FIXED)
cpc.param_set_ftype (1, Ncm.ParamType.FIXED)
cpc.param_set_ftype (2, Ncm.ParamType.FIXED)
cpc.param_set_ftype (3, Ncm.ParamType.FIXED)

cpc.param_set_by_name ('lnMcut', 33.76)
cpc.param_set_by_name ('sigma_Mcut', 0.011) 
cpc.param_set_by_name ('zmin', 0.188)
cpc.param_set_by_name ('Deltaz', 0.70214)


#cpc.param_set_upper_bound (2, 0.188)
#cpc.param_set_lower_bound (3, 0.7021)

# Cluster catalog
plclash = Nc.DataClusterPseudoCounts.new_from_file ('nc_data_cluster_planck_clash.obj')

dset = Ncm.Dataset.new ()

dset.append_data (plclash)
lh = Ncm.Likelihood (dataset = dset)

fit = Ncm.Fit.new (Ncm.FitType.NLOPT, 'ln-neldermead', lh, mset, Ncm.FitGradType.NUMDIFF_CENTRAL)

Ncm.func_eval_set_max_threads (NT)
Ncm.func_eval_log_pool_stats ()

init_sampler = Ncm.MSetTransKernGauss.new (0)
esmcmc = Ncm.FitESMCMC.new (fit, NChains, init_sampler, Ncm.FitESMCMCMoveType.STRETCH, Ncm.FitRunMsgs.FULL)
init_sampler.set_mset (mset)
init_sampler.set_prior_from_mset ()

esmcmc.set_nthreads (NT)
init_sampler.set_cov_from_scale ()

esmcmc.set_data_file ('test.fits') #("pl_cl_m2lnL_esmcmc_al_selec_fixed.fits")

esmcmc.start_run ()
esmcmc.run (2000)    # Number of points to be computed in each chain
esmcmc.end_run ()

esmcmc.mean_covar ()
fit.log_covar ()


















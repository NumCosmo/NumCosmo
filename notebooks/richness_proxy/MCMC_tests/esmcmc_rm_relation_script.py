from numcosmo_py import Ncm, Nc, GObject
Ncm.cfg_init()
Ncm.cfg_set_log_handler(lambda msg: sys.stdout.write(msg) and sys.stdout.flush())
from richness_mass_calib import create_richness_mass_calib
import numpy as np


ascaso = Nc.ClusterMassAscaso()
mset = Ncm.MSet()

lnrich_ext = Nc.ClusterMassLnrichExt()

def catalog_fit(DATA, LINEAR):
    rmdata = create_richness_mass_calib(DATA, mass_col_name = 'mass', redshift_col_name = 'redshift' )
    dset = Ncm.Dataset.new()
    dset.append_data(rmdata)

    lh = Ncm.Likelihood.new(dset)
    
    if LINEAR:
        fixed_parameters_ascaso = ['cut'] #fixing cut parameter
        ascaso.param_set_by_name("cut", 1e2) #Set cut parameter value 
        mset.set(ascaso)
        rmdata.m2lnL_val(mset)  
        mset.param_set_all_ftype(Ncm.ParamType.FREE) #All parameters free

        #All parameters free except cut parameters:
        for par in fixed_parameters_ascaso:
             mset["NcClusterMass"].param_set_desc(par, {"fit": False})
            
        mset.prepare_fparam_map()

    else:
        fixed_parameters_lnrich_ext = ['A0','cut', 'cutM1', 'cutZ1'] #fixing cut parameters

        lnrich_ext.param_set_by_name("cut", 1e2) #Set cut parameter value 
        mset.set(lnrich_ext)
        rmdata.m2lnL_val(mset)  
        mset.param_set_all_ftype(Ncm.ParamType.FREE) #All parameters free

        #All parameters free except cut parameters:
        for par in fixed_parameters_lnrich_ext:
             mset["NcClusterMass"].param_set_desc(par, {"fit": False})

        mset.prepare_fparam_map()

    fit = Ncm.Fit.factory( Ncm.FitType.NLOPT, "ln-neldermead", lh, mset, Ncm.FitGradType.NUMDIFF_CENTRAL )
    fit.log_info()
    fit.run_restart(Ncm.FitRunMsgs.SIMPLE, 1.0e-3, 0.0, None, None)
    fit.log_info()
        
    return fit




def esmcmc(DATA, N_WALKERS, N_RUN, MODEL, FILE_NAME):
    
    if MODEL == 'asc':
        fit = catalog_fit(DATA, LINEAR=True)
    
    else:
        fit = catalog_fit(DATA, LINEAR=False)

    
    Ncm.func_eval_set_max_threads(2)
    Ncm.func_eval_log_pool_stats()

    init_sampler = Ncm.MSetTransKernGauss.new(0)
    init_sampler.set_mset(mset)
    init_sampler.set_prior_from_mset()
    init_sampler.set_cov_from_rescale(1.0)

    apes = Ncm.FitESMCMCWalkerAPES.new(N_WALKERS, mset.fparams_len())

    esmcmc = Ncm.FitESMCMC.new(fit, N_WALKERS, init_sampler, apes, Ncm.FitRunMsgs.FULL)
    esmcmc.set_nthreads(2)
    esmcmc.set_data_file(FILE_NAME)

    esmcmc.start_run()
    esmcmc.run(N_RUN)  
    esmcmc.end_run()

    esmcmc.mean_covar()

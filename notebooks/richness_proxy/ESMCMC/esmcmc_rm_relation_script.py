import numpy as np

from numcosmo_py import Ncm, Nc, GObject
Ncm.cfg_init()
Ncm.cfg_set_log_handler(lambda msg: sys.stdout.write(msg) and sys.stdout.flush())


ascaso = Nc.ClusterMassAscaso()
lnrich_ext = Nc.ClusterMassLnrichExt()

mset = Ncm.MSet()

def catalog_fit(DATA, rich_cut, LINEAR):

    
    # ---------------------------------------------------------------------------- #

    # DataClusterMassRich object
    
    lnM_v = Ncm.Vector.new(len(DATA))
    z_v = Ncm.Vector.new(len(DATA))
    rich_v = Ncm.Vector.new(len(DATA))
    
    for i, mass in enumerate(DATA['mass']):
        lnM_v.set(i, np.log(mass))
    
    for i, z in enumerate(DATA['redshift']):
        z_v.set(i, z)
    
    for i, rich in enumerate(DATA['richness']):
        rich_v.set(i, np.log(rich))
    
    dmr = Nc.DataClusterMassRich.new()
    dmr.set_data(lnM_v, z_v, rich_v)
   
    
    # ---------------------------------------------------------------------------- #

    # Dataset
    
    dset = Ncm.Dataset.new()
    dset.append_data(dmr)

    
    # ---------------------------------------------------------------------------- #

    # Likelihood
    
    lh = Ncm.Likelihood.new(dset)

    
    # ---------------------------------------------------------------------------- #

    fixed_parameters = []
    
    # ---------------------------------------------------------------------------- #
    
    # Linear Model
    
    if LINEAR:
        
        # fixing cut parameter
        fixed_parameters_ascaso = ['cut'] 
        ascaso.param_set_by_name("cut", np.log(rich_cut))        

        # Mset 
        mset.set(ascaso) 
        
        dmr.m2lnL_val(mset)
        
        mset.param_set_all_ftype(Ncm.ParamType.FREE) #All parameters free

        #All parameters free except cut parameters:
        for par in fixed_parameters_ascaso:
             mset["NcClusterMass"].param_set_desc(par, {"fit": False})
            
        mset.prepare_fparam_map()

         

    
    # ---------------------------------------------------------------------------- #

    # Quadratic Model
    
    else:
        
        #fixing cut parameters
        fixed_parameters_lnrich_ext = ['A0','cut', 'cutM1', 'cutZ1'] 

        lnrich_ext.param_set_by_name("cut", np.log(rich_cut))

        # Mset 
        mset.set(lnrich_ext)
        
        dmr.m2lnL_val(mset)
        
        mset.param_set_all_ftype(Ncm.ParamType.FREE) #All parameters free

        #All parameters free except cut parameters:
        for par in fixed_parameters_lnrich_ext:
             mset["NcClusterMass"].param_set_desc(par, {"fit": False})

        mset.prepare_fparam_map()

         
    
    # ---------------------------------------------------------------------------- #

    fit = Ncm.Fit.factory( Ncm.FitType.NLOPT, "ln-neldermead", lh, mset, Ncm.FitGradType.NUMDIFF_CENTRAL )
    fit.log_info()
    fit.run_restart(Ncm.FitRunMsgs.SIMPLE, 1.0e-3, 0.0, None, None)
    fit.log_info()
        
    return fit



def esmcmc(DATA, rich_cut, N_WALKERS, N_RUN, MODEL, FILE_NAME):
    
    if MODEL == 'asc':
        fit = catalog_fit(DATA, rich_cut, LINEAR=True)
    
    else:
        fit = catalog_fit(DATA, rich_cut, LINEAR=False)

    
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


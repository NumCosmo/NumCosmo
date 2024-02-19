def Binningf(data_set):    
    
    z_ds = data_set["redshift_true"]
    lnm_ds = np.log(data_set["m200c"])
    
    d = 0.05
    d_m = 0.25
    
    z_bins = int((max(z_ds) - min (z_ds)) // d  + 1)
    m_bins = int((max(lnm_ds) - min (lnm_ds)) // d_m + 1)

    z_0 = 0.0
    z_1 = d
    halos_bin_z =[]

    for i in range(z_bins):
        cut_z = np.logical_and (data_set['redshift_true'] > z_0, data_set['redshift_true'] < z_1)
        halos_bin_z.append(data_set[cut_z])
        z_0 = z_0 + d
        z_1 = z_1 + d

    # mass bins
    label = []
    halos_bin_mz =[]
    for i in range(z_bins):

        lnM_0 = min(lnM)
        lnM_1 = min(lnM) + d_m
        for j in range(m_bins):

            cut = np.logical_and (np.log(halos_bin_z[i]["m200c"]) > lnM_0, np.log(halos_bin_z[i]["m200c"]) < lnM_1)
            halos_bin_mz.append(halos_bin_z[i][cut])
            label.append(f"{min(halos_bin_z[i]['redshift_true']):.3f} < z < {max(halos_bin_z[i]['redshift_true']):.3f}\n{lnM_0:.3f} < lnM < {lnM_1:.3f}")

            lnM_0 = lnM_0 + d_m
            lnM_1 = lnM_1 + d_m
    
    return halos_bin_mz





def bin_meanf(data_set):
    lnM_binned, z_binned, lnR_binned = [], [], [] 

    binned_halos = Binningf(data_set)

    for i in range(len(binned_halos)):

        halos = binned_halos[i]
        lnM_binned.append(np.log(halos["m200c"]))
        z_binned.append(halos["redshift_true"])
        lnR_binned.append(np.log(halos["richness"]))   

    lnR_mean, lnM_mean, z_mean = [np.mean(l) for l in lnR_binned if len(l) > 0], [np.mean(l) for l in lnM_binned if len(l) > 0], [np.mean(k) for k in z_binned if len(k) > 0]
    
    lnR_std = np.array([np.std(l) for l in lnR_binned if len(l) > 0])

    
    halos_mean = Table([np.exp(np.array(lnR_mean)), np.exp(np.array(lnM_mean)), z_mean],
               names=('richness', 'm200c', 'redshift_true'))
    
    return halos_mean, lnM_binned, z_binned, lnR_binned, lnR_std

def Model_fit(mod, data_set):
    
    #data_set
    dt_halos = Table(DC2_halos_m200c[1].data)
    rmdata = create_richness_mass_calib(data_set)
    
    fixed_parameters = [] 
    
    #Swicth
    match mod:
        case "ext_ln1pz":
            model = Nc.ClusterMassLnrichExt()
            fixed_parameters = [12, 13, 14] #fixing cut parameters
            
        case "ext_z":
            model = Nc.ClusterMassLnrichExt()
            model.set_properties('use_ln1pz', False)
            fixed_parameters = [12, 13, 14] #fixing cut parameters
            
        case "ascaso":
            model = Nc.ClusterMassAscaso()
            fixed_parameters = [6] #fixing cut parameter
    
    #Model
    model.param_set_by_name("cut", 1e15) #Set cut parameter value 
    mset = Ncm.MSet()
    mset.set(model)
    rmdata.m2lnL_val(mset)  
    mset.param_set_all_ftype(Ncm.ParamType.FREE) #All parameters free
    
    #Data
    dset = Ncm.Dataset.new()
    dset.append_data(rmdata)
    
    #Likelihood
    lh = Ncm.Likelihood.new(dset)
  
    #All parameters free except cut parameters:
    for par in fixed_parameters:
        mset.param_set_ftype(7000, par, Ncm.ParamType.FIXED)
    
    mset.prepare_fparam_map()
    
    #Fit
    fit = Ncm.Fit.factory( Ncm.FitType.NLOPT, "ln-neldermead", lh, mset, Ncm.FitGradType.NUMDIFF_CENTRAL )
    fit.log_info()
    fit.run_restart(Ncm.FitRunMsgs.SIMPLE, 1.0e-3, 0.0, None, None)
    fit.log_info()
    
    #Binning data_set
    bin_f= bin_meanf(data_set)
    halos_mean = bin_f[0]
    lnM_mean = np.log(halos_mean["m200c"])
    z_mean = halos_mean["redshift_true"]
    
    # Mean and std of data_set z mean and lnM mean
    lnR_mean_model = np.array([model.get_mean_richness(lnM_mean[i], z_mean[i]) for i in range(len(halos_mean))])
    lnR_std_model = np.array( [model.get_std_richness(lnM_mean[i], z_mean[i]) for i in range(len(halos_mean))])
    
    return lnR_mean_model, lnR_std_model, model
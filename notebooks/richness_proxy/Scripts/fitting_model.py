from numcosmo_py import Ncm, Nc, GObject
import numpy as np
import pandas as pd
from astropy.table import Table
from sklearn.model_selection import train_test_split


from richness_mass_calib import create_richness_mass_calib

import sys
sys.path.insert(0, "/global/homes/c/cinlima/NumCosmo/notebooks/richness_proxy/Scripts")
from bdata import BinningData

#-------------------------------------------------------------------------------------------------#
#FittingModel
#-------------------------------------------------------------------------------------------------#

class FittingModel:
    
    def __init__(self, data_set, b_z, b_m):
        self.data_set = data_set


        
#----- Training and test sets --------------------------------------------------------------------#

# X: input data.
# y: target data.

    def cvdata(self, X, y):
        
        X_train, X_test, y_train, y_test = train_test_split(
            X, y, test_size=0.3, random_state=0)

        #Training data
        data_train = pd.concat([X_train, y_train], axis=1)

        #Test  data
        data_test = pd.concat([X_test, y_test], axis=1)
        
        return data_train, data_test


    

#----- Fitting the models -----------------------------------------------------------------------#

# mod: model name ('ext_ln1pz' or 'ext_z' or 'ascaso');
# dt: data to calculate the fitted model;
# b_z: length of redshift bin;
# b_m: length of mass bin.


    def model_fit(self, mod, dt, b_z, b_m):
        
    #data_set
        rmdata = create_richness_mass_calib(self.data_set, mass_col_name = 'mass', redshift_col_name = 'redshift' )
        
        fixed_parameters = [] 
    
    #Swicth
        match mod:
            case "ext_ln1pz":
                model = Nc.ClusterMassLnrichExt(use_ln1pz = True)
                fixed_parameters = [13, 14, 15] #fixing cut parameters
                # model.param_set_by_name("muZ2", 0) #Set cut parameter value

            
            case "ext_z":
                model = Nc.ClusterMassLnrichExt(use_ln1pz = False)
                fixed_parameters = [13, 14, 15] #fixing cut parameters
            
            case "ascaso":
                model = Nc.ClusterMassAscaso()
                fixed_parameters = [6] #fixing cut parameter
                # model.param_set_by_name("sigmap2", 0) #Set cut parameter value 
                # model.param_set_by_name("mup2", 0) #Set cut parameter value 


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
    
    
    # Here we calculate the fitted model using mean data of bins:
        #Binning data
        bd = BinningData(dt)
        bin_f= bd.get_mean_bd(b_z, b_m)
       
        halos_mean = bin_f[0]
        lnM_mean = np.log(halos_mean["mass"])
        z_mean = halos_mean["redshift"]
    
        # Mean and std of data_set z mean and lnM mean
        lnR_mean_model = np.array([model.get_mean_richness(lnM_mean[i], z_mean[i]) for i in range(len(halos_mean))])
        
        lnR_std_model = np.array( [model.get_std_richness(lnM_mean[i], z_mean[i]) for i in range(len(halos_mean))])
    
        return lnR_mean_model, lnR_std_model, model
    

   

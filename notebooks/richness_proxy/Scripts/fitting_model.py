from numcosmo_py import Ncm, Nc, GObject
Ncm.cfg_init()
import numpy as np
import pandas as pd
from astropy.table import Table
from sklearn.model_selection import train_test_split

from richness_mass_calib import create_richness_mass_calib


#-------------------------------------------------------------------------------------------------#
#FittingModel

# data: data to calculate the fitted model;
# b_z: length of redshift bin;
# b_m: length of mass bin.

#-------------------------------------------------------------------------------------------------#

class FittingModel:
    
    def __init__(self, data_set):
        self.data_set = data_set
       

        
#----------------------------------------------------------------------------#
# train_test_data(X, y)
    
# X: input data.
# y: target data.
#----------------------------------------------------------------------------#

    def train_test_data(self, X, y):
        
        X_train, X_test, y_train, y_test = train_test_split(
            X, y, test_size=0.3, random_state=0)

        #Training data
        data_train = pd.concat([X_train, y_train], axis=1)

        #Test  data
        data_test = pd.concat([X_test, y_test], axis=1)
        
        return data_train, data_test


    

#----------------------------------------------------------------------------#  
# run_fit(mod, training)        
# mod: model name ('ext_ln1pz' or 'ext_z' or 'ascaso');
# training: If training = True the fitting is calculated using training data.

#----------------------------------------------------------------------------#

    def run_fit(self, mod):                       

    #data_set
        rmdata = create_richness_mass_calib(self.data_set, mass_col_name = 'mass', redshift_col_name = 'redshift' )    
    
    #Swicth
        fixed_parameters = [] 

        match mod:
            case "ext_ln1pz":
                model = Nc.ClusterMassLnrichExt(use_ln1pz = True)
                fixed_parameters = ['A0','cut', 'cutM1', 'cutZ1'] #fixing cut parameters
            
            case "ext_z":
                model = Nc.ClusterMassLnrichExt(use_ln1pz = False)
                fixed_parameters = ['A0' ,'cut', 'cutM1', 'cutZ1'] #fixing cut parameters
            
            case "ascaso":
                model = Nc.ClusterMassAscaso()
                fixed_parameters = ['cut'] #fixing cut parameter
              
    #Model
        model.param_set_by_name("cut", np.log(self.data_set['richness'].min())) #Set cut parameter value 
        # model.param_set_by_name("cut", np.log(10e-3))
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
             mset["NcClusterMass"].param_set_desc(par, {"fit": False})
    
        mset.prepare_fparam_map()
    
    #Fit
        fit = Ncm.Fit.factory( Ncm.FitType.NLOPT, "ln-neldermead", lh, mset, Ncm.FitGradType.NUMDIFF_CENTRAL )
        fit.log_info()
        fit.run_restart(Ncm.FitRunMsgs.SIMPLE, 1.0e-3, 0.0, None, None)
        fit.log_info()

        return model, fit

    
       
#----------------------------------------------------------------------------#
    
#get_mean_model(model, lnM,  z):

#----------------------------------------------------------------------------#        

    def get_mean_model(self, model, lnM, z):
        
        return np.array([model.get_mean_richness(lnM[i], z[i]) for i in range(len(lnM))])


#----------------------------------------------------------------------------#
    
#get_std_model(model, lnM,  z):

#----------------------------------------------------------------------------#        
    def get_std_model(self, model, lnM, z):
    
        return np.array([model.get_std_richness(lnM[i], z[i]) for i in range(len(lnM))])
            
            
    

   

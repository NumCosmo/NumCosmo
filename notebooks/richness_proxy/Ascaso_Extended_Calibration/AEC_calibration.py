from numcosmo_py import Ncm, Nc, GObject
Ncm.cfg_init()
import numpy as np
# from richness_mass_calib import create_richness_mass_calib


#-------------------------------------------------------------------------------------------------#
# Class: RichnessMassRelationCalibration
#
# data: data to calculate the fitted model;
#-------------------------------------------------------------------------------------------------#

class RichnessMassRelationCalibration:
    
    def __init__(self, data_set, rich_cut):
        self.data_set = data_set
        self.rich_cut = rich_cut



#-------------------------------------------------------------------------------------------------#
# create_richness_mass_calib()
# Description: Create a DataClusterMassRich object.
#-------------------------------------------------------------------------------------------------#

    def create_richness_mass_calib(self, mass_col_name: str="mass", redshift_col_name: str="redshift"):
           
        catalog_len = len(self.data_set)
    
        lnM_v = Ncm.Vector.new(catalog_len)
        z_v = Ncm.Vector.new(catalog_len)
        rich_v = Ncm.Vector.new(catalog_len)
    
        for i, mass in enumerate(self.data_set[mass_col_name]):
            lnM_v.set(i, np.log(mass))
    
        for i, z in enumerate(self.data_set[redshift_col_name]):
            z_v.set(i, z)
    
        for i, rich in enumerate(self.data_set['richness']):
            rich_v.set(i, np.log(rich))
    
        dmr = Nc.DataClusterMassRich.new()
        dmr.set_data(lnM_v, z_v, rich_v)
    
        return dmr
        
#----------------------------------------------------------------------------#  
# run_fit(mod, training)        
# mod: model name ('ext_ln1pz' or 'ext_z' or 'ascaso');
#
# return: model, fit, rmdata
#----------------------------------------------------------------------------#

    def run_fit(self, mod):                       

    #data_set
        rmdata = self.create_richness_mass_calib()    
    
    #Swicth
        fixed_parameters = [] 

        match mod:
            case "ext_ln1pz":
                model = Nc.ClusterMassLnrichExt(use_ln1pz = True, lnRichness_min=np.log(1.0), lnRichness_max=6)
                fixed_parameters = ['A0','cut', 'cutM1', 'cutZ1'] #fixing cut parameters
            
            case "ext_z":
                model = Nc.ClusterMassLnrichExt(use_ln1pz = False)
                fixed_parameters = ['A0' ,'cut', 'cutM1', 'cutZ1'] #fixing cut parameters
            
            case "ascaso":
                model = Nc.ClusterMassAscaso(lnRichness_min=np.log(1.0), lnRichness_max=6)
                fixed_parameters = ['cut'] #fixing cut parameter
              
    #Model
        model.param_set_by_name("cut", np.log(self.rich_cut)) #Set cut parameter value
        
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

        return model, fit, rmdata

          
#----------------------------------------------------------------------------#   
# get_mean_model(model, lnM,  z)
#
# return: Numpy array
#----------------------------------------------------------------------------#        

    def get_mean_model(self, model, lnM, z):
        
        return np.array([model.get_mean(lnM[i], z[i]) for i in range(len(lnM))])


#----------------------------------------------------------------------------#
#get_std_model(model, lnM,  z)
#
# return: Numpy array
#----------------------------------------------------------------------------#        
    def get_std_model(self, model, lnM, z):
    
        return np.array([model.get_std(lnM[i], z[i]) for i in range(len(lnM))])
            
            
    

   

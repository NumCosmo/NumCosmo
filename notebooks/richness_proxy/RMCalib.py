from numcosmo_py import Ncm, Nc, GObject
import numpy as np
import pandas as pd
from astropy.table import Table
from sklearn.model_selection import train_test_split


from richness_mass_calib import create_richness_mass_calib




#-------------------------------------------------------------------------------------------------#
#BinningData

# We make cuts in redshift (z) data with width d and put this bins in halos_bin_z. Then, we get each z bin and cut in lnM bins with width d_m. The bins are put in halos_bin_mz.

# For each bin in halos_bin_mz is calculated the avarage of lnR , z,  and lnM, and the results are put in halos_mean.
#-------------------------------------------------------------------------------------------------#
class BinningData:
    def __init__(self, data_set):
        self.data_set = data_set
        
    def binningf(self, d, d_m):
            
            z_ds = self.data_set["redshift_true"]
            lnm_ds = np.log(self.data_set["m200c"])
        
            # d = 0.05
            # d_m = 0.5
            
            # z_bins = 1
            # m_bins = 8
        
            z_bins = int( ( (max(z_ds) - min (z_ds)) // d )  + 1)
            m_bins = int( ( (max(lnm_ds) - min (lnm_ds)) // d_m ) + 1)
            
            #Redshift cut
            z_0 = 0.0
            z_1 = d
            # z_1 = max(self.data_set['redshift_true'])
            halos_bin_z =[]
            for i in range(z_bins):
                cut_z = np.logical_and(self.data_set['redshift_true'] > z_0, self.data_set['redshift_true'] < z_1)
                halos_bin_z.append(self.data_set[cut_z])
                z_0 = z_0 + d
                z_1 = z_1 + d

            # Mass cut
            #label = []
            halos_bin_mz = []
            lenbins = []
            for i in range(z_bins):

                lnM_0 = min(lnm_ds)
                lnM_1 = min(lnm_ds) + d_m
                for j in range(m_bins):

                    cut = np.logical_and (np.log(halos_bin_z[i]["m200c"]) > lnM_0, np.log(halos_bin_z[i]["m200c"]) < lnM_1)
                    if len(halos_bin_z[i][cut]) < 2: # cut off bins with 0 or 1 elements.
                        pass
                    
                    else:
                        #-------------------------------------------------
                        if len(halos_bin_z[i][cut]) < 6:                    # To calculate percentage of bins with 
                            lenbins.append(len(halos_bin_z[i][cut]))  # < 6 elements.
                        else:
                            pass
                        # #-------------------------------------------------
                        halos_bin_mz.append(halos_bin_z[i][cut])
                    #label.append(f"{min(halos_bin_z[i]['redshift_true']):.3f} < z < {max(halos_bin_z[i]['redshift_true']):.3f}\n{lnM_0:.3f} < lnM < {lnM_1:.3f}")                 
                    
                    lnM_0 = lnM_0 + d_m
                    lnM_1 = lnM_1 + d_m
            
            n_perc_bin = ( len(lenbins) / len(halos_bin_mz) ) * 100   # Percentage of bins with < 6 elements.
            print(f'{n_perc_bin}% of bins contains < 6 elements')
            
            return halos_bin_mz 

    
    
    def bin_meanf(self, d, d_m):
        lnM_binned, z_binned, lnR_binned = [], [], [] 

        binned_halos = self.binningf(d, d_m)

        for i in range(len(binned_halos)):

            halos = binned_halos[i]
            lnM_binned.append(np.log(halos["m200c"]))
            z_binned.append(halos["redshift_true"])
            lnR_binned.append(np.log(halos["richness"]))   

        lnR_mean, lnM_mean, z_mean = [np.mean(l) for l in lnR_binned if len(l) > 0], [np.mean(l) for l in lnM_binned if len(l) > 0], [np.mean(k) for k in z_binned if len(k) > 0]

        lnR_std = np.array([ np.sqrt( np.var(l)) for l in lnR_binned])
        

        halos_mean = Table([np.exp(np.array(lnR_mean)), np.exp(np.array(lnM_mean)), z_mean],
                   names=('richness', 'm200c', 'redshift_true'))

        return halos_mean, lnM_binned, z_binned, lnR_binned, lnR_std



    
#-------------------------------------------------------------------------------------------------#
#FittingModel

#-------------------------------------------------------------------------------------------------#

class FittingModel:
    
    def __init__(self, data_set, b_z, b_m):
        self.data_set = data_set

## Training and test sets 

    def cvdata(self, X, y):
        
        X_train, X_test, y_train, y_test = train_test_split(
            X, y, test_size=0.3, random_state=0)

        #Training data
        data_train = pd.concat([X_train, y_train], axis=1)

        #Test  data
        data_test = pd.concat([X_test, y_test], axis=1)
        
        return data_train, data_test

## Fitting the models

    def model_fit(self, mod, dt, b_z, b_m):
        #mod: Model
        #dt: data to calculate the fitted model
    
    #data_set
        rmdata = create_richness_mass_calib(self.data_set)
        
        fixed_parameters = [] 
    
    #Swicth
        match mod:
            case "ext_ln1pz":
                model = Nc.ClusterMassLnrichExt(use_ln1pz = True)
                fixed_parameters = [4, 10, 13, 14, 15] #fixing cut parameters
                model.param_set_by_name("muZ2", 0) #Set cut parameter value

            
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
    
    
    # Here we calculate the fitted model:
        #Binning data
        bd = BinningData(dt)
        bin_f= bd.bin_meanf(b_z, b_m)
        halos_mean = bin_f[0]
        print(len(halos_mean))
        lnM_mean = np.log(halos_mean["m200c"])
        z_mean = halos_mean["redshift_true"]
    
        # Mean and std of data_set z mean and lnM mean
        lnR_mean_model = np.array([model.get_mean_richness(lnM_mean[i], z_mean[i]) for i in range(len(halos_mean))])
        lnR_std_model = np.array( [model.get_std_richness(lnM_mean[i], z_mean[i]) for i in range(len(halos_mean))])
    
        return lnR_mean_model, lnR_std_model, model
    

   

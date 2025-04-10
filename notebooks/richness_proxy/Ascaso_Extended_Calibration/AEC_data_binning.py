from numcosmo_py import Ncm, Nc, GObject
Ncm.cfg_init()
import numpy as np
from astropy.table import Table


#-------------------------------------------------------------------------------------------------#
# Class: BinnedData 
#-------------------------------------------------------------------------------------------------#

class BinnedData:
    def __init__(self, data_set, z_bin_length, m_bin_length):
        self.data_set = data_set
        self.zl = z_bin_length
        self.ml = m_bin_length
        self.mz_bins = self.get_mz_bins()
        
#-------------------------------------------------------------------------------------------------#
# get_mz_bins()
#
# Returns: 
# halos_bin_mz : The bins of the data;
#-------------------------------------------------------------------------------------------------# 
    
    def get_mz_bins(self):

            #To extract data
            z_data, lnm_data = self.data_set["redshift"], np.log(self.data_set["mass"])

            #Number of bins
            z_bins = int((max(z_data) - min(z_data)) // self.zl) + 1
            m_bins = int((max(lnm_data) - min(lnm_data)) // self.ml) + 1
            
            halos_bin_mz = []
            
            for i in range(z_bins):
                cut_z = (z_data > i * self.zl) & (z_data < (i + 1) * self.zl)
                halos_z = self.data_set[cut_z]
                
                for j in range(m_bins):
                    lnM_0, lnM_1 = min(lnm_data) + j * self.ml, min(lnm_data) + (j + 1) * self.ml
                    cut_m = (np.log(halos_z["mass"]) > lnM_0) & (np.log(halos_z["mass"]) < lnM_1)
                    bin_subset = halos_z[cut_m]
                    
                    if len(bin_subset) > 1:
                        halos_bin_mz.append(bin_subset)
            
            return halos_bin_mz
        
        
#-------------------------------------------------------------------------------------------------#
# get_lnM_binned()
#
# Returns: 
# Average of lnM values in each bin.
#-------------------------------------------------------------------------------------------------#    
   
    def get_lnM_binned(self):

        return [np.log(binned_data["mass"]) for binned_data in self.mz_bins]
                       
#-------------------------------------------------------------------------------------------------#
# get_z_binned()
#
# Returns: 
# Average of z values in each bin.
#-------------------------------------------------------------------------------------------------#    
    
    def get_z_binned(self):

        return [binned_data["redshift"] for binned_data in self.mz_bins]

#-------------------------------------------------------------------------------------------------#
# get_lnR_binned()
#
# Returns: 
# Average of lnR values in each bin.
#-------------------------------------------------------------------------------------------------#    
    
    def get_lnR_binned(self):

        return [np.log(binned_data["richness"]) for binned_data in self.mz_bins]
    

#-------------------------------------------------------------------------------------------------#
# get_bins_mean()
#
# Returns: 
# Mean of lnR values in each bin.
#-------------------------------------------------------------------------------------------------#        
  
    def get_bins_mean(self):
                
        # lnR mean, lnM mean and z mean values in each bin:
        richness_mean = [np.average(binned_data["richness"], weights = binned_data["richness_err"] ) for binned_data in self.mz_bins if len(binned_data) > 0]
        
        mass_mean = [np.mean(binned_data["mass"]) for binned_data in self.mz_bins if len(binned_data) > 0]
        
        z_mean = [np.mean(binned_data["redshift"]) for binned_data in self.mz_bins if len(binned_data) > 0]
        
        return Table([richness_mean, mass_mean, z_mean], names=('richness', 'mass', 'redshift'))
    

#-------------------------------------------------------------------------------------------------#
# get_bins_std()
#
# Returns: 
# Std of lnR values in each bin.
#-------------------------------------------------------------------------------------------------#    

    def get_bins_std(self):

        lnr_mean = np.log(self.get_bins_mean()['richness'])

        std = []
        
        for i, lnrbin in enumerate(self.get_lnR_binned()):
        
            N = len(lnrbin)
            d2 = abs(lnrbin - lnr_mean[i])**2
            var = d2.sum() / (N - 1) 
            std.append(var**0.5)

        return std

        # return np.array( [np.std(lnrbin, ddof=1) for lnrbin in self.get_lnR_binned()] )




    

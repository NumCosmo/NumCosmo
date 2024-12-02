from numcosmo_py import Ncm, Nc, GObject
import numpy as np
from astropy.table import Table

#-------------------------------------------------------------------------------------------------#
#BinningData
#-------------------------------------------------------------------------------------------------#

class BinningData:
    def __init__(self, data_set):
        self.data_set = data_set

#-------------------------------------------------------------------------------------------------#
# data_binning(D_Z, D_M)
#
# D_Z: Length of redshift bin; 
# D_M: Length of mass bin.
#
# We make cuts in redshift (z) data with length D_Z and put this bins in halos_bin_z. Then, we get each z bin and cut in lnM bins with length D_M. The bins are put in halos_bin_mz.
#
# Returns: 
# halos_bin_mz : The bins of the data;
#-------------------------------------------------------------------------------------------------#    

    def data_binning(self, D_Z, D_M):
            
            #Data to slice
            z_data = self.data_set["redshift"]
            lnm_data = np.log(self.data_set["mass"])
            
            # Number of bins
            z_bins = int( ( (max(z_data) - min (z_data)) // D_Z )  + 1)
            m_bins = int( ( (max(lnm_data) - min (lnm_data)) // D_M ) + 1)
            
            
            #REDSHIFT CUTS
            z_0 = 0.0
            z_1 = D_Z
            halos_bin_z =[]
            
            for i in range(z_bins):
                cut_z = np.logical_and(z_data > z_0, z_data < z_1)
                halos_bin_z.append(self.data_set[cut_z])
                z_0 = z_0 + D_Z
                z_1 = z_1 + D_Z

            
            # MASS CUTS
            #label = []
            halos_bin_mz = []
            lenbins = []
            for i in range(z_bins):

                lnM_0 = min(lnm_data)
                lnM_1 = min(lnm_data) + D_M
                
                for j in range(m_bins):

                    cut = np.logical_and (np.log(halos_bin_z[i]["mass"]) > lnM_0, np.log(halos_bin_z[i]["mass"]) < lnM_1)
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
                    
                    lnM_0 = lnM_0 + D_M
                    lnM_1 = lnM_1 + D_M
            
            #n_perc_bin = ( len(lenbins) / len(halos_bin_mz) ) * 100   # Percentage of bins with < 6 elements.
            #print(f'{n_perc_bin}% of bins contains < 6 elements')
            
            return halos_bin_mz 


#-------------------------------------------------------------------------------------------------#
# get_mean_bd(D_Z, D_M)
#
# D_Z: Length of redshift bin; 
# D_M: Length of mass bin.
#
# For each bin in halos_bin_mz is calculated the avarage of lnR , z, and lnM, and the results are put in halos_mean.
#
# Returns: 
# halos_mean : average of lnR, lnM and z values in each bin;
# lnM_binned, z_binned, lnR_binned: lnM, z and lnR values in each bin;
# lnR_std: Std of lnR in each bin.
#-------------------------------------------------------------------------------------------------#    
   
    def get_mean_bd(self, D_Z,  D_M):
        
        lnM_binned, z_binned, lnR_binned = [], [], [] 
        binned_halos = self.data_binning(D_Z,  D_M)

        for i in range(len(binned_halos)):

            halos = binned_halos[i]
            
            # lnR, lnM and z values in each bin 
            lnM_binned.append(np.log(halos["mass"]))
            z_binned.append(halos["redshift"])
            lnR_binned.append(np.log(halos["richness"]))   
        
        # Average of lnR, lnM and z values in each bin:
        lnR_mean = [np.mean(l) for l in lnR_binned if len(l) > 0]
        lnM_mean = [np.mean(l) for l in lnM_binned if len(l) > 0]
        z_mean = [np.mean(k) for k in z_binned if len(k) > 0] 
        
        # Std of lnR in each bin
        lnR_std = np.array([np.sqrt( np.var(l)) for l in lnR_binned])

        halos_mean = Table([np.exp(np.array(lnR_mean)), np.exp(np.array(lnM_mean)), z_mean],
                   names=('richness', 'mass', 'redshift'))

        return halos_mean, lnM_binned, z_binned, lnR_binned, lnR_std



    

import numpy as np
from astropy.io import fits
from astropy.table import Table
from numcosmo_py import Ncm, Nc
import pandas as pd

import timeit
import tqdm

Ncm.cfg_init()

#Class used to match objects in the sky halo-halo, cluster-halo, cluster-cluster.

class NcSkyMatching:
    
    def __init__(self, cat1, cat1_coordinates, cat2, cat2_coordinates,cat1_properties=None, cat2_properties=None):
        
        self.cat1             = cat1
        self.cat2             = cat2
        
        self.cat1_coordinates = cat1_coordinates
        self.cat2_coordinates = cat2_coordinates
              
        self.cat1_properties  = cat1_properties 
        self.cat2_properties  = cat2_properties
    

    def ra_dec_to_theta_phi(self , ra, dec):
        """Convert RA and DEC to theta and phi."""
        ra_rad = np.radians(ra)
        dec_rad = np.radians(dec)

        theta = np.pi / 2 - dec_rad
        phi = ra_rad

        return theta, phi


    def print_fits_columns(self, i):
        """Print FITS columns."""
        match i:
            case 1:
                with fits.open(self.cat1) as hdul:
                    for hdu in hdul:
                        if isinstance(hdu, (fits.BinTableHDU, fits.TableHDU)):
                            print(hdu.columns)
                            return  hdu.data[self.cat1_coordinates['RA']], hdu.data[self.cat1_coordinates["DEC"]], hdu.data[self.cat1_coordinates["z"]]
                
            case 2:
                with fits.open(self.cat2) as hdul:
                    for hdu in hdul:
                        if isinstance(hdu, (fits.BinTableHDU, fits.TableHDU)):
                            print(hdu.columns)
                            return hdu.data[self.cat2_coordinates['RA']], hdu.data[self.cat2_coordinates["DEC"]], hdu.data[self.cat2_coordinates["z"]]


    def process_halos(self, cosmo, matching_distance, n_nearest_neighbours, match_file):
        """Process halos."""
        matched = {'ID':[],'RA':[], 'DEC':[] , 'z':[],  'ID_matched': [],  'RA_matched': [], 'DEC_matched': [], 'z_matched':[], "distances Mpc":[]}
        if self.cat1_properties != None:
                for prop in self.cat1_properties:
                    matched[prop] = []
                    
        if self.cat2_properties != None:
                for prop in self.cat2_properties:
                    matched[prop] = []
        
        file1 = self.cat1
        file2 = self.cat2
        
        hdul = fits.open(file1)
        hdu = hdul[1].data
        
        hdul2 = fits.open(file2)
        hdu2 = hdul2[1].data
        
        # Print columns for each file
        ra1, dec1, z1 = self.print_fits_columns(1)
        ra2, dec2, z2 = self.print_fits_columns(2)
        theta1, phi1 = self.ra_dec_to_theta_phi(ra1, dec1)
        theta2, phi2 = self.ra_dec_to_theta_phi(ra2, dec2)

        # z1 = z1[:1000000]
        # theta1 = theta1[:1000000]
        # phi1 = phi1[:1000000]

        snn = Ncm.SphereNN()
        dist = Nc.Distance.new(3.0)
        dist.prepare(cosmo)

        r_a = np.array([dist.comoving(cosmo, z) for z in z1])

        snn.insert_array(r_a, theta1, phi1)

        snn.rebuild()

        #
        # for theta, phi in zip(theta2, phi2):
        #    snn.knn_search(theta, phi, 10)
        # Do the above using tqdm

        for i, (theta, phi, z) in tqdm.tqdm(
            enumerate(zip(theta2, phi2, z2)), total=len(theta2)
        ):  
            r = dist.comoving(cosmo, z) 
            distances, indices = np.array(snn.knn_search_distances(r, theta, phi, n_nearest_neighbours))
            distances = np.sqrt(distances) * cosmo.RH_Mpc()
            
            matched['RA'].append(hdu2[self.cat2_coordinates['RA']][i])
            matched['DEC'].append(hdu2[self.cat2_coordinates['DEC']][i])
            matched['z'].append(hdu2[self.cat2_coordinates['z']][i])
            matched['ID'].append(i)
            
            if self.cat2_properties != None:
                for prop in self.cat2_properties:
                    matched[prop].append(hdu2[self.cat2_properties.get(prop)][i])
            
            distances_matched = []
            ID_matched        = []
            RA_matched        = []
            DEC_matched       = []
            z_matched         = []
            new_prop          = []
            for index in range(len(indices)):
                if distances[index] <= matching_distance:
                    
                    RA_matched.append(hdu[self.cat1_coordinates['RA']][int(indices[index])])
                    DEC_matched.append(hdu[self.cat1_coordinates['DEC']][int(indices[index])])
                    z_matched.append(hdu[self.cat1_coordinates['z']][int(indices[index])])
                    distances_matched.append(distances[index])
                    ID_matched.append(indices[index])
                    
                    
                    if self.cat1_properties != None:
                        for prop in self.cat1_properties:
                            new_prop.append(hdu[self.cat1_properties.get(prop)][int(indices[index])])
    
            
            
            if self.cat1_properties != None:
                for prop in self.cat1_properties:
                    matched[prop].append(new_prop)
            
            
            #indices = indices[indices[0] < 1e-8]
            #suborder = np.argsort(np.abs(z1[indices] - z))
            #indices = indices[suborder]
            matched['ID_matched'].append(ID_matched)
            matched['distances Mpc'].append(distances_matched)
            matched['RA_matched'].append(RA_matched)
            matched['DEC_matched'].append(DEC_matched)
            matched['z_matched'].append(z_matched)
            
            
        #matched_df = pd.DataFrame(matched)
        table = Table(matched)
        #table.write(match_file)
        return matched



if __name__ == "__main__":
    # Run precess halos and time it
    print(timeit.timeit(process_halos, number=1))

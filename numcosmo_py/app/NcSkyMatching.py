import numpy as np
from astropy.io import fits
from numcosmo_py import Ncm, Nc

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
                            return hdu.data[self.cat1_coordinates['RA']], hdu.data[self.cat1_coordinates["DEC"]], hdu.data[self.cat1_coordinates["z"]]
                
            case 2:
                with fits.open(self.cat2) as hdul:
                    for hdu in hdul:
                        if isinstance(hdu, (fits.BinTableHDU, fits.TableHDU)):
                            print(hdu.columns)
                            return hdu.data[self.cat2_coordinates['RA']], hdu.data[self.cat2_coordinates["DEC"]], hdu.data[self.cat2_coordinates["z"]]


    def process_halos(self):
        """Process halos."""
        file1 = self.cat1
        file2 = self.cat2

        # Print columns for each file
        ra1, dec1, z1 = self.print_fits_columns(1)
        ra2, dec2, z2 = self.print_fits_columns(2)
        theta1, phi1 = self.ra_dec_to_theta_phi(ra1, dec1)
        theta2, phi2 = self.ra_dec_to_theta_phi(ra2, dec2)

        # z1 = z1[:1000000]
        # theta1 = theta1[:1000000]
        # phi1 = phi1[:1000000]

        snn = Ncm.SphereNN()
        cosmo = Nc.HICosmoDEXcdm()
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
            indices = np.array(snn.knn_search(r, theta, phi, 100))
            # suborder = np.argsort(np.abs(z1[indices] - z))
            # indices = indices[suborder]

            continue
            for index in indices:
                sepa_angle = np.acos(
                    np.cos(theta) * np.cos(theta1[index])
                    + np.sin(theta) * np.sin(theta1[index]) * np.cos(phi - phi1[index])
                )
                print(
                    f"Iteration {i}: "
                    f"theta = {theta: .2f}, phi = {phi: .2f}, "
                    f"diff_theta = {theta1[index] - theta: .2e}, "
                    f"diff_phi = {phi1[index] - phi: .2e}, "
                    f"sepa_angle = {sepa_angle:.2e}, "
                    f"z = {z:.4f}, z_halo = {z1[index]:.4f}"
                )


if __name__ == "__main__":
    # Run precess halos and time it
    print(timeit.timeit(process_halos, number=1))

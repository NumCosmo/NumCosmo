import numpy as np
import pandas as pd
from astropy.table import Table, vstack
import pandas as pd

from numcosmo_py import Nc, Ncm, sky_match 
from numcosmo_py.sky_match import (
    BestCandidates,
    Coordinates,
    DistanceMethod,
    SelectionCriteria,
    SkyMatch,
    SkyMatchResult,
)

from astropy import units as u
from astropy.coordinates import SkyCoord, Angle
from scipy.special import erf
Ncm.cfg_init()


class MockGenerator:
    """Class to generate clusters, halos and galaxy members."""

    C_LIGHT = Ncm.C.c()

    def __init__(
        self,
        cosmo: Nc.HICosmo,
        halo_set_size: tuple[float, float] | None = 200,
        cluster_set_size: tuple[float, float] | None = 100,
        ra_interval: tuple[float, float] | None = (-10, 10),
        dec_interval: tuple[float, float] | None = (-10, 10),
        z_interval: tuple[float, float] | None = (0.2, 0.5),
        cluster_mass_interval: tuple[float, float] | None = (1e14, 1e15),
        halo_mass_interval: tuple[float, float] | None = (1e13, 1e14),
        seed: int | None = None
    ) -> None:
        """Create a new MockGenerator object."""
        self.cosmo = cosmo
        self.halo_set_size = halo_set_size
        self.cluster_set_size = cluster_set_size
        self.ra_interval = ra_interval
        self.dec_interval = dec_interval
        self.z_interval = z_interval
        self.cluster_mass_interval = cluster_mass_interval
        self.halo_mass_interval = halo_mass_interval
        self.seed = seed
        if self.seed is not None:
            np.random.seed(self.seed)
        
    @property
    def ra_min(self):
        """Return the minimum RA of the mock catalog."""
        return self.ra_interval[0]

    @property
    def ra_max(self):
        """Return the maximum RA of the mock catalog."""
        return self.ra_interval[1]
    
    @property
    def dec_min(self):
        """Return the minimum DEC of the mock catalog."""
        return self.dec_interval[0]       
    
    @property
    def dec_max(self):
        """Return the maximum DEC of the mock catalog."""   
        return self.dec_interval[1]
    
    @property
    def z_min(self):
        """Return the minimum redshift of the mock catalog."""
        return self.z_interval[0]
    
    @property
    def z_max(self):
        """Return the maximum redshift of the mock catalog."""
        return self.z_interval[1]
    
    @property
    def cluster_mass_min(self):
        """Return the minimum cluster mass of the mock catalog."""
        return self.cluster_mass_interval[0]
    
    @property
    def cluster_mass_max(self):
        """Return the maximum cluster mass of the mock catalog."""
        return self.cluster_mass_interval[1]
    
    @property
    def halo_mass_min(self):
        """Return the minimum halo mass of the mock catalog."""
        return self.halo_mass_interval[0]
    
    @property
    def halo_mass_max(self):
        """Return the maximum halo mass of the mock catalog."""
        return self.halo_mass_interval[1]         
    


    def generate_cluster_positions(self):
        """Generate random cluster positions within the specified RA, Dec, and redshift intervals.
        
        :return: ndarray cluster_ra: cluster right ascension (RA), ndarray cluster_dec: cluster declination (Dec), ndarray cluster_z: redshift (z)
        """
        
        cluster_ra = np.random.uniform(
            self.ra_min, self.ra_max, self.cluster_set_size
        )

        cluster_sin_dec = np.random.uniform(
            np.sin(np.radians(self.dec_min)), np.sin(np.radians(self.dec_max)), self.cluster_set_size
        )

        cluster_dec = np.degrees(
            np.arcsin(cluster_sin_dec)
        )

        cluster_z = np.random.uniform(
            self.z_min, self.z_max, self.cluster_set_size
        )
        
        return cluster_ra, cluster_dec, cluster_z


    def generate_cluster_logm(self):
        """Generate random log10(masses) within the specified mass interval."""

        cluster_logm = np.random.uniform(np.log10(self.cluster_mass_min), np.log10(self.cluster_mass_max), self.cluster_set_size)

        return cluster_logm
    
    

    def get_3D_coordinates(self, RA, DEC, R):
        """Get positions in 3D coordinates."""
        
        x1 = (
            R * np.cos(np.radians(DEC)) * np.cos(np.radians(RA))
        )
        x2 = (
            R * np.cos(np.radians(DEC)) * np.sin(np.radians(RA))
        )
        x3 = R * np.sin(np.radians(DEC))
         
        return x1, x2, x3
    


    def rho_crit(self,cosmo, z):
        return  2.8*10**(11)*self.cosmo.h2() * ((self.cosmo["Omegac"]+self.cosmo["Omegab"])/(1+z)**3 + (1 - self.cosmo["Omegac"]-self.cosmo["Omegab"])) 


    def get_R200c(self, mass, z):
        """Calculate R_{200c} for a given mass and redshift."""
        
        return ((3 * mass) / (4 * np.pi * 200 * self.rho_crit(self.cosmo, z)))**(1/3)
    
    def generate_halo_cluster(self, D_DIM = 2.0):

        """Generate halos and clusters."""
                
        # Generate cluster positions, redshifts, and masses
        cluster_ra, cluster_dec, cluster_z = self.generate_cluster_positions()
        cluster_logm = self.generate_cluster_logm()
        
        # Compute the cluster radii
        dist = Nc.Distance.new(100.0)
        dist.compute_inv_comoving(True)
        dist.prepare(self.cosmo)
        cluster_r = np.array(dist.comoving_array(self.cosmo, cluster_z)) * self.cosmo.RH_Mpc()
        
        #Compute the 3D coordinates of the clusters
        cluster_x1, cluster_x2, cluster_x3 = self.get_3D_coordinates(cluster_ra, cluster_dec, cluster_r)
        
        # Generate halo positions, redshifts, and masses around the clusters
        halo_x1 = cluster_x1 + np.random.uniform(-D_DIM, D_DIM, self.cluster_set_size)
        halo_x2 = cluster_x2 + np.random.uniform(-D_DIM, D_DIM, self.cluster_set_size)
        halo_x3 = cluster_x3 + np.random.uniform(-D_DIM, D_DIM, self.cluster_set_size)
        halo_ra = np.degrees(np.arctan2(halo_x2, halo_x1))
        halo_dec = np.degrees(np.arcsin(halo_x3 / cluster_r))
        halo_r = np.sqrt(halo_x1**2 + halo_x2**2 + halo_x3**2)
        halo_z = [dist.inv_comoving(self.cosmo, r / self.cosmo.RH_Mpc()) for r in halo_r]
        
        # Now for the halo masses we use the cluster's masses added a Gaussian noise
        halo_logm = cluster_logm + np.random.normal(0, 0.1, self.cluster_set_size)
        
        # Add more halos randomly
        DELTA_OBJECTS = self.halo_set_size - self.cluster_set_size
        halo_ra = np.append(halo_ra, np.random.uniform(self.ra_min, self.ra_max, DELTA_OBJECTS))
        halo_dec = np.append(halo_dec, np.random.uniform(self.dec_min, self.dec_max, DELTA_OBJECTS))
        halo_z = np.append(halo_z, np.random.uniform(self.z_min, self.z_max, DELTA_OBJECTS))
        halo_logm = np.append(
            halo_logm, np.random.uniform(self.logm_add_halo_min, self.logm_add_halo_max, DELTA_OBJECTS)
        )
        
        # Compute R200 for halos and clusters
        halos_R200c = self.get_R200c(10**halo_logm, halos_z)
        clusters_R200c = self.get_R200c(10**cluster_logm, cluster_z)

        #create tables for halos and clusters
        
        halo_id = np.array([int(i + 200000) for i in range(len(halo_z))])
        halos = Table([halo_id, halo_ra, halo_dec, halo_z, halo_logm], 
                    names=("halo_id", "RA", "DEC", "z", "Mass"), 
                    dtype=(int, float, float, float, float))
        
        cluster_id = np.array([int(i + 100000) for i in range(len(cluster_z))])
        clusters = Table([cluster_id, cluster_ra, cluster_dec, cluster_z, cluster_logm], 
                    names=("cluster_id", "RA", "DEC", "z", "Mass"), 
                    dtype=(int, float, float, float, float))

        halo_mass = 10**(halos['Mass'])
        cluster_mass = 10**(clusters['Mass'])
        halos['R200'] = ((3*(halo_mass))/(4*np.pi* 200*(rho_crit(self.cosmo, halos['z']))))**(1/3)
        clusters['R200'] = ((3*(cluster_mass))/(4*np.pi* 200*(rho_crit(self.cosmo, clusters['z']))))**(1/3)

        return halos, clusters


    
    def get_galaxy_coords(ra_c, dec_c, sep_angular_rad, phi_rad):
        """
        Calcula RA e DEC de uma galáxia usando trigonometria esférica.
        
        Parâmetros:
        ra_c, dec_c: Coordenadas do centro (em graus)
        sep_angular_rad: galaxy_R / dA (em radianos)
        phi_rad: Ângulo de posição sorteado (em radianos)
        """
        
        # 1. Converter coordenadas do centro para radianos
        ra_c_rad = np.radians(ra_c)
        dec_c_rad = np.radians(dec_c)
        
        # 2. Calcular a Declinação (delta)
        # sin(delta) = sin(dec_c)cos(sep) + cos(dec_c)sin(sep)cos(phi)
        sin_dec = (np.sin(dec_c_rad) * np.cos(sep_angular_rad) + 
                np.cos(dec_c_rad) * np.sin(sep_angular_rad) * np.cos(phi_rad))
        
        dec_gal_rad = np.arcsin(sin_dec)
        
        # 3. Calcular a Ascensão Reta (alpha)
        # Usamos o termo y e x para o arctan2(y, x)
        y = np.sin(phi_rad) * np.sin(sep_angular_rad) * np.cos(dec_c_rad)
        x = np.cos(sep_angular_rad) - np.sin(dec_c_rad) * np.sin(dec_gal_rad)
        
        ra_gal_rad = ra_c_rad + np.arctan2(y, x)
        
        # 4. Converter de volta para graus e normalizar o RA para [0, 360]
        ra_gal = np.degrees(ra_gal_rad) % 360
        dec_gal = np.degrees(dec_gal_rad)
        
        return ra_gal, dec_gal



    def HOD_model(m_halo, logMmin=12.72, sigma_logM=0.26, alpha=1.15, logM1=13.93, logM0=12.7):
        """
        Standard HOD model (Zheng et al. 2007)
        Returns number of centrals (0 or 1) and satellites (int)
        """
        # 1. Mean Central Occupancy (Error function)
        from scipy.special import erf
        mean_n_cen = 0.5 * (1 + erf((np.log10(m_halo) - logMmin) / sigma_logM))
        #n_cen = 1 if np.random.random() < mean_n_cen else 0
        n_cen = 1
        
        # 2. Mean Satellite Occupancy (Power Law)
        n_sat = 0
        if n_cen == 1:
            # Satellites only exist if a central exists
            diff = m_halo - 10**logM0
            mean_n_sat = (np.maximum(0, diff) / 10**logM1)**alpha
            if mean_n_sat > 0:
                n_sat = np.random.poisson(mean_n_sat)
                

        return n_cen, n_sat



    def generate_halo_galaxies(halos,cosmo,dist, RA_interval=[-10,10], DEC_interval=[-10,10], z_interval=[0.2,0.5], seed=74682):
        np.random.seed(seed)
        all_galaxies = []
        RA_MIN, RA_MAX = RA_interval[0], RA_interval[1]
        DEC_MIN, DEC_MAX = DEC_interval[0], DEC_interval[1]
        Z_MIN, Z_MAX = z_interval[0], z_interval[1]
        for i in range(len(halos)):
            n_cen, n_sat = HOD_model(10**halos['Mass'][i])
            total_gals = n_cen + n_sat
            if total_gals == 0:
                continue
                
            # Generate local coordinates
            # Central is at (0,0,0), Satellites follow NFW/Uniform profile
            galaxy_R = np.zeros(total_gals)
            if n_sat > 0:
                # Satellites distributed within R200
                galaxy_R[n_cen:] = halos['R200'][i] * np.random.uniform(0, 1, n_sat)**(1/3)
        
            h_ra = halos['RA'][i]
            h_dec = halos['DEC'][i]
            h_z = halos['z'][i]
            h_r200 = halos['R200'][i]
        
            galaxy_R = h_r200 * np.random.uniform(0, 1, total_gals)**(1/3)
            costheta = np.random.uniform(-1, 1, total_gals)
            galaxy_theta = np.arccos(costheta)
            galaxy_phi = np.random.uniform(0, 2*np.pi, total_gals)
        
            Hz = cosmo.H(h_z)
            
            distancia_radial = galaxy_R * np.cos(galaxy_theta) 
            delta_z = (Hz / C_LIGHT) * distancia_radial
            z_galaxy = h_z + delta_z
        
        
            dA = dist.angular_diameter(cosmo, h_z) * cosmo.RH_Mpc()
            sep_angular = galaxy_R / dA 
        
            ra_g, dec_g = get_galaxy_coords(h_ra, h_dec, sep_angular, galaxy_phi)
        
            ra_corrected = ra_g % 360.0
            
            ra_corrected = np.where(ra_corrected > 180, ra_corrected - 360, ra_corrected)
            
            ra_corrected = np.clip(ra_corrected, RA_MIN, RA_MAX)
            
            
            temp_table = Table()
        
            temp_table['halo_id'] = np.full(total_gals, i + 200000)
            temp_table['RA'] = ra_corrected
            temp_table['DEC'] = dec_g
            temp_table['z'] = z_galaxy
            temp_table['is_central'] = [True if j < n_cen else False for j in range(total_gals)]
            temp_table['halo_RA']  = np.full(total_gals, h_ra)
            temp_table['halo_DEC'] = np.full(total_gals, h_dec)
            temp_table['halo_z']   = np.full(total_gals, h_z)
            temp_table['halo_mass']   = np.full(total_gals, halos['Mass'][i])
            
            all_galaxies.append(temp_table)
        
        if all_galaxies:
            galaxy_halo_catalog = vstack(all_galaxies)
            galaxy_halo_catalog['galaxy_id'] = np.arange(len(galaxy_halo_catalog)) + 2000000
            
            # Reorder columns for readability
            ordered_cols = ['galaxy_id', 'halo_id', 'is_central', 'RA', 'DEC', 'z', 'halo_z', 'halo_mass']
            galaxy_halo_catalog = galaxy_halo_catalog[ordered_cols]
            
            print(f"Generated {len(galaxy_halo_catalog)} galaxies from {len(halos)} halos.")
        else:
            print("No galaxies were generated. Check your HOD mass thresholds.")

        return galaxy_halo_catalog



    def generate_cluster_galaxies(clusters,cosmo,dist, RA_interval=[-10,10], DEC_interval=[-10,10], z_interval=[0.2,0.5], seed=74682):
        np.random.seed(seed)
        all_galaxies = []

        RA_MIN, RA_MAX = RA_interval[0], RA_interval[1]
        DEC_MIN, DEC_MAX = DEC_interval[0], DEC_interval[1]
        Z_MIN, Z_MAX = z_interval[0], z_interval[1]
        
        for i in range(len(clusters)):
            n_cen, n_sat = HOD_model(10**clusters['Mass'][i])
            
            total_gals = n_cen + n_sat
            if total_gals == 0:
                continue
                
            # Generate local coordinates
            # Central is at (0,0,0), Satellites follow NFW/Uniform profile
            galaxy_R = np.zeros(total_gals)
            if n_sat > 0:
                # Satellites distributed within R200
                galaxy_R[n_cen:] = clusters['R200'][i] * np.random.uniform(0, 1, n_sat)**(1/3)
        
            c_ra = clusters['RA'][i] 
            c_dec = clusters['DEC'][i]
            c_z = clusters['z'][i]
            c_r200 = clusters['R200'][i]
        
            galaxy_R = c_r200 * np.random.uniform(0, 1, total_gals)**(1/3)
            costheta = np.random.uniform(-1, 1, total_gals)
            galaxy_theta = np.arccos(costheta)
            galaxy_phi = np.random.uniform(0, 2*np.pi, total_gals)
        
            Hz = cosmo.H(c_z)

            distancia_radial = galaxy_R * np.cos(galaxy_theta) 
            delta_z = (Hz / C_LIGHT) * distancia_radial
            z_galaxy = c_z + delta_z
        
        
            dA = dist.angular_diameter(cosmo, c_z) * cosmo.RH_Mpc()
            sep_angular = galaxy_R / dA 
        
            ra_g, dec_g = get_galaxy_coords(c_ra, c_dec, sep_angular, galaxy_phi)
        
            ra_corrected = ra_g % 360.0
            
            ra_corrected = np.where(ra_corrected > 180, ra_corrected - 360, ra_corrected)
            
            ra_corrected = np.clip(ra_corrected, RA_MIN, RA_MAX)
            
            
            temp_table = Table()
        
            temp_table['cluster_id']  = np.full(total_gals, i + 100000)
            temp_table['RA'] = ra_corrected
            temp_table['DEC'] = dec_g
            temp_table['z'] = z_galaxy
            temp_table['is_central'] = [True if j < n_cen else False for j in range(total_gals)]
            temp_table['cluster_RA']  = np.full(total_gals, c_ra)
            temp_table['cluster_DEC'] = np.full(total_gals, c_dec)
            temp_table['cluster_z']   = np.full(total_gals, c_z)
            temp_table['cluster_mass']   = np.full(total_gals, clusters['Mass'][i])
            
            all_galaxies.append(temp_table)
        
        if all_galaxies:
            galaxy_cluster_catalog = vstack(all_galaxies)
            galaxy_cluster_catalog['galaxy_id'] = np.arange(len(galaxy_cluster_catalog)) + 1000000
            
            # Reorder columns for readability
            ordered_cols = ['galaxy_id', 'cluster_id', 'is_central', 'RA', 'DEC', 'z', 'cluster_z','cluster_mass']
            galaxy_cluster_catalog = galaxy_cluster_catalog[ordered_cols]
            
            print(f"Generated {len(galaxy_cluster_catalog)} galaxies from {len(clusters)} clusters.")
        else:
            print("No galaxies were generated. Check your HOD mass thresholds.")

        return galaxy_cluster_catalog




    def generate_common_galaxies(halos,halo_galaxies,cluster_galaxies, cosmo, halo_properties, detections_properties, obj):
        halo_coordinates = {"RA":"RA" , "DEC":"DEC" , "z":"z", "ID":"ids"}   
        detections_coordinates =  {"RA":"RA" , "DEC":"DEC" , "z":"z"}
        
        
        halos_m = sky_match.SkyMatch(
                    query_data=halos, 
                    query_coordinates=halo_coordinates,
                    match_data=cluster_galaxies,
                    match_coordinates=detections_coordinates
                )
        halos_matched = halos_m.match_3d(cosmo, 20)
        halos_table = halos_matched.to_table_complete(query_properties=halo_properties , match_properties=detections_properties)
        
        # 1. Criar a lista de matches filtrados manualmente para evitar o erro de conversão
        # Isso substitui o explode e o filtro de distância em um único passo eficiente
        matched_data = []
        for h7 in halos_table:
            # Máscara booleana: True para galáxias dentro do R200
            mask = h7['distances'] <= h7['R200']
            
            # Pegamos apenas os índices que passaram no filtro
            indices = h7['Index_matched'][mask]
            
            # Criamos uma linha para cada galáxia encontrada, repetindo os dados do halo
            for idx in indices:
                matched_data.append({
                    'matched_index': idx,
                    '%s_id' % (obj): h7['%s_id' % (obj)],
                    '%s_z' % (obj): h7['z'],
                    '%s_mass' % (obj): h7['%s_mass' % (obj)]
                })
        
        # 2. Transformar os matches em DataFrame
        df_matched = pd.DataFrame(matched_data)
        
        # 3. Converter a tabela de galáxias para Pandas
        # Filtramos apenas as colunas 1D para evitar o erro multidimensional se elas existirem lá também
        names_gal = [name for name in cluster_galaxies.colnames if len(cluster_galaxies[name].shape) <= 1]
        df_galaxies = cluster_galaxies[names_gal].to_pandas()
        
        # 4. Join (Merge) para trazer os dados RA, DEC, z, etc., das galáxias
        # Usamos o 'matched_index' que guardamos acima para bater com o índice original da cluster_galaxies
        result_df = df_matched.merge(
            df_galaxies, 
            left_on='matched_index', 
            right_index=True
        )
        
        # 5. Organizar o catálogo final
        halo_galaxies_add = result_df[[
            'galaxy_id', '%s_id' % (obj), 'RA', 'DEC', 'z', '%s_z' % (obj), '%s_mass' % (obj)
        ]].copy()
        
        halo_galaxies_add['is_central'] = False

        halo_galaxies_add = Table.from_pandas(halo_galaxies_add)
        halo_galaxies = vstack([halo_galaxies, halo_galaxies_add])
        
        return halo_galaxies

#!/usr/bin/env python
# coding: utf-8

# In[8]:


#!/usr/bin/env python

import math
from scipy.stats import norm
import numpy
import sys
sys.path.insert(1, '../../../Programas_Cosmologia/PySSC/')
sys.path.insert(1, '../../../Programas_Cosmologia/CLASS/')
import PySSC
from classy import Class
try:
    import gi

    gi.require_version('NumCosmo', '1.0')
    gi.require_version('NumCosmoMath', '1.0')
except:
    pass
from gi.repository import GObject
from gi.repository import NumCosmo as Nc
from gi.repository import NumCosmoMath as Ncm

Ncm.cfg_init()
Ncm.cfg_set_log_handler (lambda msg: sys.stdout.write (msg) and sys.stdout.flush ())

# Creating a new class implementing our object Ncm.Data
#

class ncounts(Ncm.DataGaussCov):
    #
    # We need one vector property to save the independent variables x
    #
    ncdata = GObject.Property(type=Ncm.Vector, flags=GObject.PARAM_READWRITE)
    z_obs_bins = GObject.Property(type=Ncm.Vector , flags=GObject.PARAM_READWRITE)
    lnM_obs_bins = GObject.Property(type=Ncm.Vector , flags=GObject.PARAM_READWRITE)
    ca = GObject.Property(type=Nc.ClusterAbundance, flags=GObject.PARAM_READWRITE)
    S = GObject.Property(type=Ncm.Matrix, flags=GObject.PARAM_READWRITE)

    #
    # The contructor assigns some default values and calls the father's constructor.
    #
    def __init__ (self, len = 600):
        Ncm.DataGaussCov.__init__ (self, n_points = len)
        self.np = self.get_size()
        
        if self.np > 0:
            self.ncdata = Ncm.Vector.new(self.np)
        else:
            self.ncdata = None

        self.cov_init = False

        #
        # Initializing to sane values
        #
        if self.np > 0:
            self.cov = self.peek_cov()
            #self.cov.set_identity()
            self.ncdata.set_zero()

    #
    # Implements the virtual method get_length.
    #
        self.mask = numpy.ndarray(shape = (0,0))
        self.kernel = None
    def do_get_length(self):
        return self.np

    #
    # Implements the virtual method get_dof.
    #
    def do_get_dof(self):
        return self.np


    #
    # Implements the virtual method `prepare'.
    # This method should do all the necessary calculations using mset
    # to be able to calculate the likelihood afterwards.
    #
    def do_prepare(self, mset):
        return

    #
    # Implements the virtual method `mean_func'.
    # This method should compute the theoretical mean for the gaussian
    # distribution.
    #
    def do_mean_func(self, mset, vp):
        cosmo = mset.peek(Nc.HICosmo.id())
        cluster_m = mset.peek(Nc.ClusterMass.id())
        cluster_z = mset.peek(Nc.ClusterRedshift.id())
        size_z = self.z_obs_bins.len()
        size_lmM = self.lnM_obs_bins.len()
        self.ca.prepare(cosmo, cluster_z, cluster_m)
        for i in range(size_z-1):
            for j in range(size_lmM -1):
                lnM_obs_lb = [self.lnM_obs_bins.get( j + 0)]
                lnM_obs_ub = [self.lnM_obs_bins.get( j + 1)]
                z_obs_lb = [self.z_obs_bins.get( i + 0)]
                z_obs_ub = [self.z_obs_bins.get( i + 1)]
                vp.set(i+j, self.ca.intp_bin_d2n(cosmo, cluster_z, cluster_m, lnM_obs_lb, lnM_obs_ub, None, z_obs_lb, z_obs_ub, None))
    
        return 

    def do_cov_func(self, mset,rng):


       
        cosmo = mset.peek(Nc.HICosmo.id())
        cluster_m = mset.peek(Nc.ClusterMass.id())
        cluster_z = mset.peek(Nc.ClusterRedshift.id())
        prim = cosmo.peek_prim()
        '''
        #Beggin of the S matrix calculation
        nz       = self.kernel.shape[1]
        z_arr    = numpy.linspace(0, 2,num=nz+1)[1:]
        cosmo_fid = Class()
        #cosmo_fid.set({'h':cosmo.props.H0/100 ,'Omega_cdm':cosmo.props.Omegac ,'Omega_b':cosmo.props.Omegab ,'sigma8':0.82505858,'n_s':prim.props.n_SA,'output':'mPk'})
        cosmo_fid.set({'h':0.67 ,'Omega_cdm':0.2612 ,'Omega_b':0.0486 ,'sigma8':0.82505858,'n_s':0.9660,'output':'mPk'})
        cosmo_fid.compute()
        if self.mask.shape[0] == 0:
            S_lacasa = PySSC.Sij(z_arr,self.kernel,cosmo_Class=cosmo_fid)
        else:    
            S_lacasa = PySSC.Sij_psky(z_arr,self.kernel,mask=self.mask,cosmo_Class=cosmo_fid)
        
        self.S = Ncm.Matrix.new(S_lacasa.shape[0],S_lacasa.shape[1])

        for i in range(len(S_lacasa)):
            for j in range(len(S_lacasa[i])):
                self.S.set(i,j, S_lacasa[i][j])
        #End of the S matrix calculation
        '''
        
        bias = []
        poisson = []
        self.ca.prepare(cosmo, cluster_z, cluster_m)
        size_z = self.z_obs_bins.len()
        size_lmM = self.lnM_obs_bins.len()

        for z_bin in range(size_z-1):
            bias_bin_mass = []
            poisson_bin_mass = []
            for lnM_bin in range(size_lmM-1):
                P_bias_bin = self.ca.intp_bin_d2n_bias(cosmo, cluster_z, cluster_m, [self.lnM_obs_bins.get(lnM_bin)], [self.lnM_obs_bins.get(lnM_bin+1)], None, [self.z_obs_bins.get(z_bin)], [self.z_obs_bins.get(z_bin+1)], None) 
                P_poisson_bin = self.ca.intp_bin_d2n(cosmo, cluster_z, cluster_m, [self.lnM_obs_bins.get(lnM_bin)], [self.lnM_obs_bins.get(lnM_bin+1)], None, [self.z_obs_bins.get(z_bin)], [self.z_obs_bins.get(z_bin+1)], None)
                bias_bin_mass.append(P_bias_bin)
                poisson_bin_mass.append(P_poisson_bin)
            bias.append(bias_bin_mass)
            poisson.append(poisson_bin_mass)

        row = 0
        for alpha in range(size_lmM-1):
            for i in range(size_z-1):
                col = 0
                for beta in range(size_lmM-1):
                
                    for j in range(size_z-1):
                        if i == j and alpha == beta:
                            self.cov.set( row, col ,poisson[i][alpha] + bias[i][alpha] * bias[j][beta] * self.S.get(i,j))
                        else:
                            self.cov.set( row, col ,bias[i][alpha] * bias[j][beta] * self.S.get(i,j))

                        col +=1

                row +=1
        return True
    
    def set_lnM_obs_bins(self, lnM_obs_bins):
        self.lnM_obs_bins = lnM_obs_bins

    def set_z_obs_bins(self, z_obs_bins):
        self.z_obs_bins = z_obs_bins

    def set_mask(self,mask):
        self.mask = mask
    
    def set_kernel(self,kernel):
        self.kernel = kernel
    

    def set_cad(self, cad):
        self.ca = cad

    def set_S(self):
        #Beggin of the S matrix calculation
        nz       = self.kernel.shape[1]
        z_arr    = numpy.linspace(0, 2,num=nz+1)[1:]
        cosmo_fid = Class()
        #cosmo_fid.set({'h':cosmo.props.H0/100 ,'Omega_cdm':cosmo.props.Omegac ,'Omega_b':cosmo.props.Omegab ,'sigma8':0.82505858,'n_s':prim.props.n_SA,'output':'mPk'})
        cosmo_fid.set({'h':0.67 ,'Omega_cdm':0.2612 ,'Omega_b':0.0486 ,'sigma8':0.82505858,'n_s':0.9660,'output':'mPk'})
        cosmo_fid.compute()
        if self.mask.shape[0] == 0:
            S_lacasa = PySSC.Sij(z_arr,self.kernel,cosmo_Class=cosmo_fid)
        else:    
            S_lacasa = PySSC.Sij_psky(z_arr,self.kernel,mask=self.mask,cosmo_Class=cosmo_fid)
        
        self.S = Ncm.Matrix.new(S_lacasa.shape[0],S_lacasa.shape[1])

        for i in range(len(S_lacasa)):
            for j in range(len(S_lacasa[i])):
                self.S.set(i,j, S_lacasa[i][j])
        #End of the S matrix calculation


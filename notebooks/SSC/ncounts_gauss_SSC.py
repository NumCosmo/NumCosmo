#!/usr/bin/env python
# coding: utf-8

# In[8]:


#!/usr/bin/env python

import math
from scipy.stats import norm
import numpy as np
import sys
sys.path.insert(1, '../../../Programas_Cosmologia/PySSC/')
import PySSC
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

        if self.np > 0:
            self.ncdata = Ncm.Vector.new(self.np)
        else:
            self.ncdata = None

        self.cov_init = False

        #
        # Initializing to sane values
        #
        if self.np > 0:
            self.cov.set_identity()
            self.ncdata.set_zero()

    #
    # Implements the virtual method get_length.
    #
    def do_get_length(self):
        return self.np

    #
    # Implements the virtual method get_dof.
    #
    def do_get_dof(self):
        return self.np

    #
    # Implements the virtual method `begin'.
    # This method usually do some groundwork in the data
    # before the actual calculations. For example, if the likelihood
    # involves the decomposition of a constant matrix, it can be done
    # during `begin' once and then used afterwards.
    #
    def do_prepare (self, mset):
        return

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
        
        
        for i in range(size_z-1):
            for alpha in range(size_lmM-1):
                for j in range(size_z-1):
                    for beta in range(size_lmM-1):
                        if i == j and alpha == beta:
                            self.cov.set( i+alpha, j+beta ,poisson[i][alpha] + bias[i][alpha] * bias[j][beta] * self.S.get(i,j))
                        else:
                            self.cov.set( i+alpha, j+beta ,bias[i][alpha] * bias[j][beta] * self.S.get(i,j))


    def set_lnM_obs_bins(self, lnM_obs_bins):
        self.lnM_obs_bins = lnM_obs_bins

    def set_z_obs_bins(self, z_obs_bins):
        self.z_obs_bins = z_obs_bins

    def get_lnM_obs_bins(self):
        return self.lnM_obs_bins
    
    def get_z_obs_bins(self):
        return self.z_obs_bins
    

    def set_cad(self, cad):
        self.ca = cad

    def set_S(self, S):
        self.S = S


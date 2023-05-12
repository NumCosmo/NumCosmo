#!/usr/bin/env python
# coding: utf-8

# In[59]:


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

class ncounts_gauss(Ncm.DataGaussCov):
    #
    # We need one vector property to save the independent variables x
    #
    ncdata = GObject.Property(type=Ncm.Vector, flags=GObject.PARAM_READWRITE)
    z_obs_bins = GObject.Property(type=Ncm.Vector , flags=GObject.PARAM_READWRITE)
    lnM_obs_bins = GObject.Property(type=Ncm.Vector , flags=GObject.PARAM_READWRITE)
    ca = GObject.Property(type=Nc.ClusterAbundance, flags=GObject.PARAM_READWRITE)
    #
    # The contructor assigns some default values and calls the father's constructor.
    #
    def __init__(self, len=600):
        Ncm.DataGaussCov.__init__(self, n_points=len)

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
    def do_begin(self):
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
        size = self.z_obs_bins.len()
        nbins = int(size/2)
        self.ca = Nc.ClusterAbundance
        self.ca.prepare(cosmo, cluster_z, cluster_m)
        for i in range(nbins):
            j = 2 * i
            lnM_obs_lb = [self.lnM_obs_bins.get( j + 0)]
            lnM_obs_ub = [self.lnM_obs_bins.get( j + 1)]
            z_obs_lb = [self.z_obs_bins.get( j + 0)]
            z_obs_ub = [self.z_obs_bins.get( j + 1)]
            vp.set(i, self.ca.intp_bin_d2n(cosmo, cluster_z, cluster_m, lnM_obs_lb, lnM_obs_ub, None, z_obs_lb, z_obs_ub, None))
    
        return

    def do_cov_func(self, mset):
        
        cosmo = mset.peek(Nc.HICosmo.id())
        cluster_m = mset.peek(Nc.ClusterMass.id())
        cluster_z = mset.peek(Nc.ClusterRedshift.id())
        size = self.z_obs_bins.len()
        nbins = int(size/2)
        bias = []
        poisson = []
        self.ca.prepare(cosmo, cluster_z, cluster_m)
        for i in range(nbins):
            j = 2 * i
            lnM_obs_lb = [self.lnM_obs_bins.get( j + 0)]
            lnM_obs_ub = [self.lnM_obs_bins.get( j + 1)]
            z_obs_lb = [self.z_obs_bins.get( j + 0)]
            z_obs_ub = [self.z_obs_bins.get( j + 1)]
            bias.append(self.ca.intp_bin_d2n_bias(cosmo, cluster_z, cluster_m, lnM_obs_lb, lnM_obs_ub, None, z_obs_lb, z_obs_ub, None))
            variance = poisson.append(self.ca.intp_bin_d2n(cosmo, cluster_z, cluster_m, lnM_obs_lb, lnM_obs_ub, None, z_obs_lb, z_obs_ub, None))
        variance = np.array(variance)
        z_obs_py = np.array(self.z_obs_bins)
        self.cov = PySSC.Sij_fullsky(z_obs_py, bias, order  = 1) + variance
        return self.cov
    


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


# In[60]:


gauss = ncounts_gauss()
cosmo = Nc.HICosmoDEXcdm()


reion = Nc.HIReionCamb.new () 
prim = Nc.HIPrimPowerLaw.new () 

cosmo.add_submodel (reion)
cosmo.add_submodel (prim)

dist = Nc.Distance.new (2.0)

tf = Nc.TransferFunc.new_from_name ("NcTransferFuncEH")

psml = Nc.PowspecMLTransfer.new (tf)
psml.require_kmin (1.0e-6)
psml.require_kmax (1.0e3)

psf = Ncm.PowspecFilter.new (psml, Ncm.PowspecFilterType.TOPHAT)
psf.set_best_lnr0 ()

mulf = Nc.MultiplicityFuncTinker.new ()
mulf.set_mdef (Nc.MultiplicityFuncMassDef.CRITICAL)
mulf.set_Delta (200.0)

hmf = Nc.HaloMassFunction.new (dist, psf, mulf)
hbias_Tinker = Nc.HaloBiasTinker.new(hmf)
ca = Nc.ClusterAbundance.new(hmf,hbias_Tinker)

cluster_m = Nc.ClusterMass.new_from_name("NcClusterMassNodist{'lnM-min':<%20.15e>, 'lnM-max':<%20.15e>}" % (math.log(10)*np.log10(1e14),math.log(10)*np.log10(1e16)))
cluster_z = Nc.ClusterRedshift.new_from_name("NcClusterRedshiftNodist{'z-min': <%20.15e>, 'z-max':<%20.15e>}" % (0.25,2))


mset = Ncm.MSet.new_array([cosmo,cluster_m,cluster_z])

z_obs_bins  = Ncm.Vector.new_array(np.linspace(0,1,5))
lnM_obs_bins = Ncm.Vector.new_array(np.linspace(14,15,5))
gauss.set_z_obs_bins(z_obs_bins)
gauss.set_lnM_obs_bins(lnM_obs_bins)
gauss.set_cad(ca)

gauss.do_cov_func(mset)


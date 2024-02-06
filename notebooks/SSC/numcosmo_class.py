#!/usr/bin/env python

import math
import healpy as hp
from scipy.stats import norm
import numpy as np
import sys
import time
sys.path.insert(1, '../../../PySSC/')
sys.path.insert(1, '../../../CLASS/')
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


def get_numcosmo(z_arr,cosmo_numcosmo):
    """Auxiliary routine to run CLASS if needed, then compute arrays of comoving distance, volume etc necessary for Sij routines.

    Parameters
    ----------
    z_arr : array_like
        Input array of redshifts of size nz.

    cosmo_Class : classy.Class object, default None
        classy.Class object containing precomputed cosmology, if you already have it and do not want PySSC to lose time recomputing cosmology with CLASS.

    Returns
    -------
    tuple
        cosmo, h, comov_dist, dcomov_dist, growth
    """
    dist = Nc.Distance.new (2.0)
    nz = z_arr.size
    cosmo = cosmo_numcosmo
    dist.prepare(cosmo)
    comov_dist =  np.zeros(nz)
    dcomov_dist = np.zeros(nz)
    growth      = np.zeros(nz)                              #Growth factor
    h = cosmo.props.H0/100 #for  conversions Mpc/h <-> Mpc
    growth_func = Nc.GrowthFunc.new()
    growth_func.prepare(cosmo)

    tf = Nc.TransferFuncEH.new()

    psml = Nc.PowspecMLTransfer.new (tf)
    psml.require_kmin (1.0e-6)
    psml.require_kmax (1.0e3)

    # Define arrays of r(z), k, P(k)...
    for iz in range(nz):
        comov_dist[iz]  = cosmo.RH_Mpc()*dist.comoving(cosmo,z_arr[iz])  #Comoving distance r(z) in Mpc
        dcomov_dist[iz] = cosmo.RH_Mpc()*1/cosmo.E(z_arr[iz])           #Derivative dr/dz in Mpc
        growth[iz] = growth_func.eval(cosmo,z_arr[iz])
    
    return cosmo, h, comov_dist, dcomov_dist, growth, psml
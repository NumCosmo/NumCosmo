from numcosmo_py import Ncm, Nc, GObject

import numpy as np
import math
import sys

from astropy.io import fits
from astropy.table import Table

Ncm.cfg_init()

DC2_halos_m200c = fits.open('/global/cfs/projectdirs/lsst/groups/CL/cosmoDC2_v1.1.4/extragal/full/halos/halos_m200c_13.0.fits')
dt_halos = Table(DC2_halos_m200c[1].data)

class RichnessMassCalib(Ncm.DataGaussDiag):
    """A simple data class for the SLine model based on Ncm.DataGaussCov."""
    
    __gtype_name__ = "RichnessMassCalib"
    
    lnM_v = GObject.Property(type=Ncm.Vector, flags=GObject.ParamFlags.READWRITE)
    z_v = GObject.Property(type=Ncm.Vector, flags=GObject.ParamFlags.READWRITE)

    def __init__(self, dt_halos):
        catalog_len = len(dt_halos)
        Ncm.DataGaussDiag.__init__(self, n_points=catalog_len)
        self.lnM_v = Ncm.Vector.new(catalog_len)
        self.z_v = Ncm.Vector.new(catalog_len)
        
        for i, mass in enumerate(dt_halos['m200c']):
            self.lnM_v.set(i, np.log(mass)) 

        for i, z in enumerate(dt_halos['redshift_true']):
            self.z_v.set(i, z) 

        rich_v = self.peek_mean()
        for i, rich in enumerate(dt_halos['richness']):
            rich_v.set(i, np.log(rich))
        
        self.set_init(True)

    def do_prepare(self, _mset):  # pylint: disable-msg=arguments-differ
        return

    def do_mean_func(self, mset: Ncm.MSet, vp: Ncm.Vector) -> None:
        mid = mset.get_id_by_ns("NcClusterMass")
        ascaso = mset.peek(mid)
        assert isinstance(ascaso, Nc.ClusterMassAscaso)
        
        n_points = super().get_size()

        for i in range(n_points):
            lnM = self.lnM_v.get(i)
            z = self.z_v.get(i)
            lnR = ascaso.get_mean_richness(lnM, z)
            vp.set(i, lnR)
        

    def do_sigma_func(self, mset: Ncm.MSet, sigma: Ncm.Vector) -> None:
        mid = mset.get_id_by_ns("NcClusterMass")
        ascaso = mset.peek(mid)
        assert isinstance(ascaso, Nc.ClusterMassAscaso)
        n_points = super().get_size()

        for i in range(n_points):
            lnM = self.lnM_v.get(i)
            z = self.z_v.get(i)
            lnR = ascaso.get_std_richness(lnM, z)
            sigma.set(i, lnR)

#
# Register our new Python class PySLineGauss
#
GObject.type_register(RichnessMassCalib)


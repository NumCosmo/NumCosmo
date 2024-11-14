"""A simple data class for the SLine model based on Ncm.DataGaussCov."""


import numpy as np

from astropy.io import fits
from astropy.table import Table

from numcosmo_py import Ncm, Nc, GObject

DC2_halos_m200c = fits.open(
    "/global/cfs/projectdirs/lsst/groups/CL/cosmoDC2_v1.1.4/extragal/full/halos/halos_m200c_13.0.fits"
)
dt_halos = Table(DC2_halos_m200c[1].data)


def create_richness_mass_calib(dt_halos, mass_col_name: str="m200c", redshift_col_name: str="redshift_true"):
    """Create a RichnessMassCalib object."""

    catalog_len = len(dt_halos)

    lnM_v = Ncm.Vector.new(catalog_len)
    z_v = Ncm.Vector.new(catalog_len)
    rich_v = Ncm.Vector.new(catalog_len)

    for i, mass in enumerate(dt_halos[mass_col_name]):
        lnM_v.set(i, np.log(mass))

    for i, z in enumerate(dt_halos[redshift_col_name]):
        z_v.set(i, z)

    for i, rich in enumerate(dt_halos['richness']):
        rich_v.set(i, np.log(rich))

    dmr = Nc.DataClusterMassRich.new()
    dmr.set_data(lnM_v, z_v, rich_v)

    return dmr

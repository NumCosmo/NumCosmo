from numcosmo_py import Ncm
from numcosmo_py import Nc

import math


def test_cfg_init():
    Ncm.cfg_init()


def test_distances():
    Ncm.cfg_init()

    cosmo = Nc.HICosmo.new_from_name(Nc.HICosmo, "NcHICosmoDEXcdm{'massnu-length':<1>}")
    cosmo.set_reparam(Nc.HICosmoDEReparamCMB.new(cosmo.len()))

    dist = Nc.Distance.new(2.0)

    cosmo.orig_param_set(Nc.HICosmoDESParams.H0, 70.00)
    cosmo.orig_param_set(Nc.HICosmoDESParams.OMEGA_C, 0.25)
    cosmo.orig_param_set(Nc.HICosmoDESParams.OMEGA_X, 0.70)
    cosmo.orig_param_set(Nc.HICosmoDESParams.T_GAMMA0, 2.72)
    cosmo.orig_param_set(Nc.HICosmoDESParams.OMEGA_B, 0.05)
    cosmo.orig_param_set(Nc.HICosmoDEXCDMSParams.W, -1.10)

    cosmo.orig_vparam_set(Nc.HICosmoDEVParams.M, 0, 0.06)
    dist.prepare(cosmo)

    N = 100
    RH_Mpc = cosmo.RH_Mpc()

    for i in range(N):
        z = 1.0 / (N - 1.0) * i
        Dc = dist.comoving(cosmo, z)
        dc = RH_Mpc * Dc

        assert math.isfinite(dc)

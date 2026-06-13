"""Cosmology helpers for mock catalogs.

Thin wrappers over NumCosmo so the mock generator does not re-implement
cosmological quantities. These are the natural seam for a future move into C.
"""

import numpy as np

from numcosmo_py import Nc, Ncm


def critical_density(cosmo: Nc.HICosmo, z):
    """Critical mass density at redshift ``z`` in solar masses per Mpc^3.

    ``rho_crit(z) = rho_crit0 * E2(z)``, with ``rho_crit0`` from NumCosmo's
    constants and ``E2`` the (squared) normalized expansion function.

    :param cosmo: the cosmology.
    :param z: redshift (scalar or array-like).
    :return: critical density, matching the shape of ``z``.
    """
    rho_crit0 = Ncm.C.crit_mass_density_h2_solar_mass_Mpc3() * cosmo.h2()
    e2 = np.vectorize(cosmo.E2, otypes=[float])(z)
    return rho_crit0 * e2


def r200c(cosmo: Nc.HICosmo, mass, z):
    """Radius enclosing 200x the critical density, ``R_{200c}``, in Mpc.

    ``R_{200c} = (3 M / (4 pi 200 rho_crit(z)))^{1/3}``.

    :param cosmo: the cosmology.
    :param mass: halo mass in solar masses (scalar or array-like).
    :param z: redshift (scalar or array-like).
    """
    mass_a = np.asarray(mass, dtype=float)
    return ((3.0 * mass_a) / (4.0 * np.pi * 200.0 * critical_density(cosmo, z))) ** (
        1.0 / 3.0
    )

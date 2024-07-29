#!/usr/bin/env python
#
# example_wspline.py
#
# Mon May 22 16:00:00 2023
# Copyright  2023  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# example_wspline.py
# Copyright (C) 2023 Sandro Dias Pinto Vitenti <vitenti@uel.br>
#
# numcosmo is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# numcosmo is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program.  If not, see <http://www.gnu.org/licenses/>.

"""Example testing wspline cosmology."""

from numcosmo_py import Nc, Ncm


#
#  Initializing the library objects, this must be called before
#  any other library function.
#
Ncm.cfg_init()


def test_wspline() -> None:
    """Example testing wspline cosmology."""

    #
    #  New homogeneous and isotropic cosmological model NcHICosmoDEXcdm
    #  with one massive neutrino.
    #
    zf = 3.0
    cosmo = Nc.HICosmoDEWSpline.new(12, zf)
    cosmo2 = Nc.HICosmoDEXcdm.new()
    zf = 3000.0

    #
    #  New cosmological distance objects optimizied to perform calculations
    #  up to redshift 2.0.
    #
    dist = Nc.Distance.new(zf)

    #
    #  Setting values for the cosmological model, those not set stay in the
    #  default values. Remember to use the _orig_ version to set the original
    #  parameters when a reparametrization is used.
    #

    #
    # C-like
    #
    cosmo.orig_param_set(Nc.HICosmoDESParams.H0, 70.00)
    cosmo.orig_param_set(Nc.HICosmoDESParams.OMEGA_C, 0.25)
    cosmo.orig_param_set(Nc.HICosmoDESParams.OMEGA_X, 0.70)
    cosmo.orig_param_set(Nc.HICosmoDESParams.T_GAMMA0, 2.72)
    cosmo.orig_param_set(Nc.HICosmoDESParams.OMEGA_B, 0.05)

    cosmo2.orig_param_set(Nc.HICosmoDESParams.H0, 70.00)
    cosmo2.orig_param_set(Nc.HICosmoDESParams.OMEGA_C, 0.25)
    cosmo2.orig_param_set(Nc.HICosmoDESParams.OMEGA_X, 0.70)
    cosmo2.orig_param_set(Nc.HICosmoDESParams.T_GAMMA0, 2.72)
    cosmo2.orig_param_set(Nc.HICosmoDESParams.OMEGA_B, 0.05)
    cosmo2.orig_param_set(Nc.HICosmoDEXCDMSParams.W, -1.0)

    # alpha_a = np.array(cosmo.get_alpha().dup_array())
    # w_a = -1.0 + 1.0 * alpha_a / alpha_a[-1]
    w_a = [-1.0] * 12

    cosmo.orig_vparam_set_vector(
        Nc.HICosmoDEWSplineVParams.W, Ncm.Vector.new_array(w_a)
    )
    # cosmo.orig_vparam_set_vector (Nc.HICosmoDEWSplineVParams.W, Ncm.Vector.new_array (w_a))

    #
    #  Printing the parameters used.
    #
    print("# Model parameters: ")
    cosmo.params_log_all()

    dist.prepare(cosmo)

    #
    #  Printing some distances up to redshift 1.0.
    #

    N = 20
    for i in range(0, N):
        z = zf / (N - 1.0) * i

        w = cosmo.w_de(z)
        E2 = cosmo.E2(z)
        E2Omega_de = cosmo.E2Omega_de(z)
        dE2Omega_de_dz = cosmo.dE2Omega_de_dz(z)
        d2E2Omega_de_dz2 = cosmo.d2E2Omega_de_dz2(z)

        print(
            f"{z: 22.15f} w {w: 22.15g} E2 {E2: 22.15g} "
            f"E2Omega_de {E2Omega_de: 22.15g} "
            f"dE2Omega_de_dz {dE2Omega_de_dz: 22.15g} "
            f"d2E2Omega_de_dz2 {d2E2Omega_de_dz2: 22.15g}"
        )

        w = cosmo2.w_de(z)
        E2 = cosmo2.E2(z)
        E2Omega_de = cosmo2.E2Omega_de(z)
        dE2Omega_de_dz = cosmo2.dE2Omega_de_dz(z)
        d2E2Omega_de_dz2 = cosmo2.d2E2Omega_de_dz2(z)

        print(
            f"{z: 22.15f} w {w: 22.15g} E2 {E2: 22.15g} "
            f"E2Omega_de {E2Omega_de: 22.15g} "
            f"dE2Omega_de_dz {dE2Omega_de_dz: 22.15g} "
            f"d2E2Omega_de_dz2 {d2E2Omega_de_dz2: 22.15g}"
        )


if __name__ == "__main__":
    test_wspline()

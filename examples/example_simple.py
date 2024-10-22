#!/usr/bin/env python
#
# example_simple.py
#
# Mon May 22 16:00:00 2023
# Copyright  2023  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# example_simple.py
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

"""Simple example computing cosmological distances."""

from numcosmo_py import Nc, Ncm

#
#  Initializing the library objects, this must be called before
#  any other library function.
#
Ncm.cfg_init()


def test_distances() -> None:
    """Example computing cosmological distances."""

    #
    #  New homogeneous and isotropic cosmological model NcHICosmoDEXcdm
    #  with one massive neutrino.
    #
    cosmo = Nc.HICosmoDEXcdm(massnu_length=1)
    cosmo.set_reparam(Nc.HICosmoDEReparamCMB.new(cosmo.len()))

    #
    #  New cosmological distance objects optimizied to perform calculations
    #  up to redshift 2.0.
    #
    dist = Nc.Distance.new(2.0)

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
    cosmo.orig_param_set(Nc.HICosmoDEXCDMSParams.W, -1.10)

    cosmo.orig_vparam_set(Nc.HICosmoDEVParams.M, 0, 0.06)

    #
    # OO-like
    #
    cosmo.props.H0 = 70.00
    cosmo.props.Omegab = 0.04
    cosmo.props.Omegac = 0.25
    cosmo.props.Omegax = 0.70
    cosmo.props.Tgamma0 = 2.72
    cosmo.props.w = -1.10

    massnu_v = Ncm.Vector.new_array([0.06])
    cosmo.props.massnu = massnu_v

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
    RH_Mpc = cosmo.RH_Mpc()

    for i in range(0, N):
        z = 1.0 / (N - 1.0) * i
        Dc = dist.comoving(cosmo, z)
        dc = RH_Mpc * Dc / (1.0 + z)

        print(f"{z: 10.8f} {Dc: 22.15g} [c/H0] {dc: 22.15g} [Mpc]")


if __name__ == "__main__":
    test_distances()

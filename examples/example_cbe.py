#!/usr/bin/env python
#
# example_cbe.py
#
# Fri May 19 10:00:00 2023
# Copyright  2023  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# example_cbe.py
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

"""Example of using the Class backend to object oriented interface."""


from numcosmo_py import Nc, Ncm

#
#  Initializing the library objects, this must be called before
#  any other library function.
#
Ncm.cfg_init()


def test_cbe():
    """Compare the results from CLASS using CBE to NumCosmo's."""

    #
    #  New homogeneous and isotropic cosmological model NcHICosmoDEXcdm
    #
    cosmo = Nc.HICosmo.new_from_name(Nc.HICosmo, "NcHICosmoDEXcdm{'massnu-length':<1>}")
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
    cosmo.orig_param_set(Nc.HICosmoDESParams.ENNU, 2.0328)
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
    cosmo.props.ENnu = 2.0328
    cosmo.props.w = -1.10

    massnu_v = Ncm.Vector.new_array([0.06])
    cosmo.props.massnu = massnu_v

    cosmo.omega_x2omega_k()
    cosmo.param_set_by_name("Omegak", 0.0)

    #
    #  Printing the parameters used.
    #
    print("# Model parameters: ", end="", flush=True)
    cosmo.params_log_all()

    #
    #  New CLASS backend object.
    #
    cbe = Nc.CBE.new()
    cbe.set_calc_transfer(True)

    #
    # Submodels necessary for CLASS
    #
    reion = Nc.HIReionCamb.new()
    prim = Nc.HIPrimPowerLaw.new()
    cosmo.add_submodel(reion)
    cosmo.add_submodel(prim)

    #
    # Preparing CLASS backend
    #
    cbe.prepare(cosmo)
    dist.prepare(cosmo)

    print(f"# theta100CMB {dist.theta100CMB(cosmo): 22.15e}")
    print(f"# zt          {cosmo.zt(5.0): 22.15e}")
    print(f"# Omega_mnu0  {cosmo.Omega_mnu0(): 22.15e}")
    print(f"# Press_mnu0  {cosmo.Press_mnu0(): 22.15e}")
    print(f"# Omega_k0    {cosmo.Omega_k0(): 22.15e}")

    ztest = 1.0e4
    print(f"# E2Omega_mnu   ({ztest: 22.15e}) {cosmo.E2Omega_mnu(ztest): 22.15e}")
    print(f"# E2Press_mnu   ({ztest: 22.15e}) {cosmo.E2Press_mnu(ztest): 22.15e}")
    print(
        f"# E2Omega_mnu_d ({ztest: 22.15e}) "
        f"{cosmo.E2Omega_mnu(ztest) - 3.0 * cosmo.E2Press_mnu(ztest): 22.15e}"
    )
    print(f"# E2Omega_b     ({ztest: 22.15e}) {cosmo.E2Omega_b(ztest): 22.15e}")
    print(f"# E2Omega_c     ({ztest: 22.15e}) {cosmo.E2Omega_c(ztest): 22.15e}")
    print(f"# Neff           {cosmo.Neff(): 22.15e}")

    #
    # Printing comparison between CLASS and NumCosmo background
    #
    cbe.compare_bg(cosmo, True)


if __name__ == "__main__":
    test_cbe()

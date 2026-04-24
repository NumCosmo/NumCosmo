#!/usr/bin/env python
#
# example_deceleration_function.py
#
# Mon Apr 22 16:00:00 2026
# Copyright  2026  Pedro Henrique Costa Ribeiro
# <pedrohenriquecostaribeiro@gmail.com.br>
#
# example_deceleration_function.py
# Copyright (C) 2026 Pedro Henrique Costa Ribeiro
# <pedrohenriquecostaribeiro@gmail.com.br>
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

"""Example computing comoving distances for different kinetic models"""

import numpy as np
import matplotlib.pyplot as plt

from numcosmo_py import nc, ncm

#
#  Initializing the library objects, this must be called before
#  any other library function.
#
ncm.cfg_init()


def test_deceleration_function():
    """Example computing comoving distances for different kinetic models"""

    #
    #  Four homogeneous and isotropic cosmological model, one standard
    #  LCDM model and 3 with different deceleration functions.
    #
    cosmo_lcdm = nc.HICosmoLCDM()

    cosmo_qconst = nc.HICosmoQConst(q=-0.55)

    cosmo_qlinear = nc.HICosmoQLinear(q=-0.2, qp=1.0)

    spline = ncm.SplineCubicNotaknot.new()
    cosmo_qspline = nc.HICosmoQSpline.new(spline, 6, 10.0)

    cosmo_list = [cosmo_qconst, cosmo_qlinear, cosmo_lcdm, cosmo_qspline]

    #
    #  New cosmological distance objects optimizied to perform calculations
    #  up to redshift 2.0.
    #
    dist = nc.Distance.new(2.0)

    #
    #  Setting a plot of comoving distances for each cosmological model
    #
    plt.figure(figsize=(10, 6))
    ax = plt.subplot()

    ax.grid()
    ax.set_xlabel("z", fontsize=20)
    ax.set_ylabel(r"$D_c$ (z) [Mpc]", fontsize=20)
    ax.tick_params(axis='both', labelsize=16)
    ax.set_title("Comoving Distance for different kinetic model ", fontsize=20)


    #
    #  Preparing and computing the distance for each model up to redshift 2.0
    #
    z = np.linspace(0.0, 2.0, 100)
    for cosmo in cosmo_list:
        dist.prepare(cosmo)
        RH_Mpc = cosmo.RH_Mpc()    

        comov = [RH_Mpc * dist.comoving(cosmo, z=zi) for zi in z] 

        ax.plot(z, comov, label=cosmo.__class__.__name__, lw=3)

    plt.legend(
        fontsize = 16,
        title = "Cosmological Models",
        title_fontsize = 20
    )
    plt.savefig("Comoving_distance_kinetic_models.png", dpi=300)
    plt.show()


if __name__ == "__main__":
    test_deceleration_function()
#!/usr/bin/env python
#
# example_wl_surface_mass_density.py
#
# Mon May 22 16:00:00 2023
# Copyright  2023  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# example_wl_surface_mass_density.py
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

"""Example computing weak lensing surface mass density."""

import math
import numpy as np
import matplotlib.pyplot as plt

from numcosmo_py import Nc, Ncm


#
#  Initializing the library objects, this must be called before
#  any other library function.
#
Ncm.cfg_init()


def test_wl_surface_mass_density() -> None:
    """Example computing weak lensing surface mass density."""
    #
    #  New homogeneous and isotropic cosmological model NcHICosmoDEXcdm
    #
    cosmo = Nc.HICosmoDEXcdm()
    cosmo.omega_x2omega_k()
    cosmo.param_set_by_name("Omegak", 0.0)

    #
    #  New cosmological distance objects optimizied to perform calculations
    #  up to redshift 2.0.
    #
    dist = Nc.Distance.new(1.0)

    #
    # New matter density profile
    #
    hms = Nc.HaloMCParam.new(Nc.HaloMassSummaryMassDef.CRITICAL, 200)
    nfw = Nc.HaloDensityProfileNFW.new(hms)
    # Setting 4 as Douglas. In LCDM c = 5 corresponds to cluster masses.
    # (see Lokas and G. Mamon, astro-ph/0002395)
    nfw.param_set_by_name("cDelta", 4.0)
    nfw.param_set_by_name("MDelta", 1.0e15)
    # mdelta = 1.0e15
    # cdelta = 4.0
    # delta = 200.0

    zcluster = 1.0
    zsource = 1.5

    #
    # New weak lensing surface mass density
    #
    smd = Nc.WLSurfaceMassDensity.new(dist)

    #
    #  Setting values for the cosmological model, those not set keep their
    #  default values. Remember to use the _orig_ version to set the original
    #  parameters when a reparametrization is used.
    #
    cosmo.props.H0 = 70.0
    cosmo.props.Omegab = 0.045
    cosmo.props.Omegac = 0.255
    cosmo.props.Omegax = 0.7
    cosmo.param_set_by_name(
        "Omegak", 0.0
    )  # This line sets a flat universe, modifying Omega_x (small difference).
    # This is necessary since CLASS require the inclusion of the radiation density.

    dist.prepare(cosmo)

    npoints = 500
    r_a = np.logspace(math.log10(5.0e-3), 2.0, npoints, endpoint=True)

    Sigma = []
    meanSigma = []
    DeltaSigma = []
    convergence = []
    shear = []
    reduced_shear = []
    reduced_shear_inf = []

    for i in range(0, npoints):
        ri = r_a[i]
        Sig = smd.sigma(nfw, cosmo, ri, zcluster)
        meanSig = smd.sigma_mean(nfw, cosmo, ri, zcluster)
        kappa = smd.convergence(nfw, cosmo, ri, zsource, zcluster, zcluster)
        sh = smd.shear(nfw, cosmo, ri, zsource, zcluster, zcluster)
        reds = smd.reduced_shear(nfw, cosmo, ri, zsource, zcluster, zcluster)
        reds_inf = smd.reduced_shear_infinity(
            nfw, cosmo, ri, zsource, zcluster, zcluster
        )

        Sigma.append(Sig)
        meanSigma.append(meanSig)
        DeltaSigma.append(meanSig - Sig)
        convergence.append(kappa)
        shear.append(sh)
        reduced_shear.append(reds)
        reduced_shear_inf.append(reds_inf)

        print(ri, Sig, meanSig, kappa, sh, reds)

    plt.figure(figsize=(6, 5))  # in inches
    ax = plt.subplot()

    ax.plot(r_a, Sigma, label=r"$\Sigma(R)$")
    ax.plot(r_a, meanSigma, label=r"Mean $\overline{\Sigma}(<R)$")
    ax.plot(r_a, DeltaSigma, label=r"Excess/differential  $\Delta{\Sigma}(R)$")

    ax.set_xlabel(r"$R$ [Mpc]", fontsize=14)
    ax.set_ylabel(r"Surface Mass Density", fontsize=12)
    ax.set_xscale("log")
    ax.set_yscale("log")

    ax.annotate(
        r"$[\Sigma] = M_\odot/Mpc^2$",
        xy=(0.65, 0.8),
        xycoords="axes fraction",
        fontsize=12,
    )
    ax.set_title(r"NFW, $c_{200} = 4$, $M_{200} = 10^{15} \, M_\odot$")

    plt.legend(loc=3)

    plt.savefig("wl_smd_sigmas.svg")
    plt.show()
    plt.clf()

    ax = plt.subplot()

    ax.plot(r_a, convergence, label=r"Convergence $\kappa (R)$")
    ax.plot(r_a, shear, label=r"Shear $\gamma (R)$")
    ax.plot(r_a, reduced_shear, label=r"Reduced Shear $g(R)$")
    ax.plot(r_a, reduced_shear_inf, label=r"Reduced Shear $g_\infty(R)$")

    ax.set_xlabel(r"$R$ [Mpc]", fontsize=14)
    ax.set_xscale("log")
    ax.set_yscale("log")

    ax.set_title(r"NFW, $c_{200} = 4$, $M_{200} = 10^{15} \, M_\odot$")

    plt.legend(loc=1)

    plt.savefig("wl_smd_convergence_shear.svg")
    plt.show()


if __name__ == "__main__":
    test_wl_surface_mass_density()

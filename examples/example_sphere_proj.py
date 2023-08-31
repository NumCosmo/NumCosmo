#!/usr/bin/env python
#
# example_sphere_proj.py
#
# Mon May 22 16:00:00 2023
# Copyright  2023  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# example_sphere_proj.py
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

"""Simple example computing spherical projections."""

import math
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

from matplotlib.colors import SymLogNorm
from numcosmo_py import Nc, Ncm


def compute_spherical_projection():
    """Compute spherical projections of the power spectrum."""

    matplotlib.rcParams.update({"font.size": 11})

    #
    #  Initializing the library objects, this must be called before
    #  any other library function.
    #
    Ncm.cfg_init()

    #
    #  New homogeneous and isotropic cosmological model NcHICosmoDEXcdm
    #
    cosmo = Nc.HICosmoDEXcdm(massnu_length=1)
    cosmo.omega_x2omega_k()
    cosmo.param_set_by_name("Omegak", 0.0)
    cosmo.param_set_by_name("w", -1.0)
    cosmo.param_set_by_name("Omegab", 0.045)
    cosmo.param_set_by_name("Omegac", 0.255)
    cosmo.param_set_by_name("massnu_0", 0.06)
    cosmo.param_set_by_name("ENnu", 2.0328)

    reion = Nc.HIReionCamb.new()
    prim = Nc.HIPrimPowerLaw.new()

    cosmo.param_set_by_name("H0", 67.0)

    prim.param_set_by_name("n_SA", 0.96)
    prim.param_set_by_name("ln10e10ASA", 3.0904)

    reion.param_set_by_name("z_re", 9.9999)

    cosmo.add_submodel(reion)
    cosmo.add_submodel(prim)

    #
    #  Printing the parameters used.
    #
    print("# Model parameters: ", end=" ")
    cosmo.params_log_all()
    print(f"# Omega_X0: {cosmo.E2Omega_de(0.0): 22.15g}")

    ps_lin = Nc.PowspecMLCBE.new()

    # Redshift bounds
    z_min = 0.0
    z_max = 2.0

    # Mode bounds
    k_min = 1.0e-8
    k_max = 1.0e6

    ps_lin.set_kmin(k_min)
    ps_lin.set_kmax(k_max)
    ps_lin.require_zi(z_min)
    ps_lin.require_zf(z_max)

    ps_lin.set_intern_k_min(k_min)
    ps_lin.set_intern_k_max(50.0)
    psf = Ncm.PowspecFilter.new(ps_lin, Ncm.PowspecFilterType.TOPHAT)
    psf.set_best_lnr0()

    ps = Nc.PowspecMNLHaloFit.new(ps_lin, 2.0, 1.0e-8)

    ps.set_kmin(k_min)
    ps.set_kmax(k_max)
    ps.require_zi(z_min)
    ps.require_zf(z_max)

    old_amplitude = math.exp(prim.props.ln10e10ASA)
    prim.props.ln10e10ASA = math.log((0.81 / cosmo.sigma8(psf)) ** 2 * old_amplitude)

    ell_min = 0
    ell_max = 1

    sproj = Ncm.PowspecSphereProj.new(ps, ell_min, ell_max)

    sproj.props.reltol = 1.0e-8
    sproj.props.reltol_z = 1.0e-8
    xi_i = 1.0e-1
    xi_f = 1.0e7

    sproj.set_xi_i(xi_i)
    sproj.set_xi_f(xi_f)

    sproj.prepare(cosmo)

    dist = Nc.Distance.new(1.0e11)
    dist.compute_inv_comoving(True)

    scal = Nc.Scalefactor.new(1.0e11, dist)
    gf = Nc.GrowthFunc.new()

    dist.prepare(cosmo)
    scal.prepare(cosmo)
    gf.prepare(cosmo)

    RH_Mpc = cosmo.RH_Mpc()

    # for ell in range (ell_min, ell_max + 1):
    # (lnxi, Cell) = sproj.get_ell (ell)
    # xi_a = np.exp (np.array (lnxi))

    xi_a = np.geomspace(xi_i, xi_f, num=5000)
    z_a = np.array([dist.inv_comoving(cosmo, xi / RH_Mpc) for xi in xi_a])
    index = np.logical_and((z_a > 0.0), (z_a < 10.0))

    z_a = np.linspace(0.2, 2, num=2001)[1:]

    # print (xi_a)
    # print (z_a)

    for ell in range(ell_min, ell_min + 1):
        limber = np.array(
            [
                ps.eval(cosmo, z, (0.5 + ell) / xi) / (xi * xi * RH_Mpc / cosmo.E(z))
                for (xi, z) in zip(xi_a, z_a)
            ]
        )

        m = [
            [
                sproj.eval_Cell_xi1_xi2(
                    cosmo,
                    ell,
                    z1,
                    z2,
                    RH_Mpc * dist.comoving(cosmo, z1),
                    RH_Mpc * dist.comoving(cosmo, z2),
                )
                for z1 in z_a
            ]
            for z2 in z_a
        ]
        cov = np.asanyarray(m)
        std_ = np.sqrt(np.diag(cov))
        corr = cov / np.outer(std_, std_)
        print(len(cov))
        print(len(cov[0]))
        plt.matshow(
            cov,
            extent=np.array([z_a[0], z_a[-1], z_a[0], z_a[-1]]),
            origin="lower",
            norm=SymLogNorm(1.0e-6),
        )
        plt.xlabel(r"$z_2$")
        plt.ylabel(r"$z_1$")
        # plt.matshow (corr)
        # plt.matshow (corr, norm = SymLogNorm (1.0e-4))
        plt.title(r"$C_{%d}(z_1, z_2)$" % (ell))
        plt.colorbar()
        plt.show()
    # plt.xticks(xi_a)

    # plt.xlabel (r'$\xi \; [\mathrm{Mpc}]$')
    # plt.ylabel (r'$C_\ell(\xi, \xi)$')
    # plt.legend (loc = "best")
    # plt.xscale ('log')
    # plt.yscale ('symlog', linthreshy=1.0e-8)
    # plt.xlim ([200.0, 4000.0])

    # plt.show ()



if __name__ == "__main__":
    compute_spherical_projection()

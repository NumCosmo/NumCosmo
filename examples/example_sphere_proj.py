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

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

from matplotlib.colors import SymLogNorm

from numcosmo_py import Nc, Ncm


matplotlib.rcParams.update({"font.size": 11})

#
#  Initializing the library objects, this must be called before
#  any other library function.
#
Ncm.cfg_init()


def test_spherical_projection() -> None:
    """Example computing spherical projection."""

    #
    #  New homogeneous and isotropic cosmological model NcHICosmoDEXcdm
    #
    cosmo = Nc.HICosmoDEXcdm(massnu_length=1)
    cosmo.omega_x2omega_k()
    cosmo.param_set_by_name("Omegak", 0.0)
    cosmo.param_set_by_name("w", -1.0)
    cosmo.param_set_by_name("Omegab", 0.04909244421)
    cosmo.param_set_by_name("Omegac", 0.26580755578)
    cosmo.param_set_by_name("massnu_0", 0.06)
    cosmo.param_set_by_name("ENnu", 2.0328)

    reion = Nc.HIReionCamb.new()
    prim = Nc.HIPrimPowerLaw.new()

    cosmo.param_set_by_name("H0", 67.31)

    prim.param_set_by_name("n_SA", 0.9658)
    prim.param_set_by_name("ln10e10ASA", 3.0904)

    reion.param_set_by_name("z_re", 9.9999)

    cosmo.add_submodel(reion)
    cosmo.add_submodel(prim)

    #
    #  Printing the parameters used.
    #
    print("# Model parameters: ", end="")
    cosmo.params_log_all()
    print(f"# Omega_X0: {cosmo.E2Omega_de(0.0): 22.15g}")

    ps_lin = Nc.PowspecMLCBE.new()

    # Redshift bounds
    z_min = 0.0
    z_max = 2.0
    # zdiv = 0.49999999999

    # Mode bounds
    k_min = 1.0e-6
    k_max = 1.0e5

    ps_lin.set_kmin(k_min)
    ps_lin.set_kmax(k_max)
    ps_lin.require_zi(z_min)
    ps_lin.require_zf(z_max)

    ps_lin.set_intern_k_min(k_min)
    ps_lin.set_intern_k_max(50.0)

    ps = Nc.PowspecMNLHaloFit.new(ps_lin, 2.0, 1.0e-5)

    ps.set_kmin(k_min)
    ps.set_kmax(k_max)
    ps.require_zi(z_min)
    ps.require_zf(z_max)

    ell_min = 250
    ell_max = 250

    sproj = Ncm.PowspecSphereProj.new(ps, ell_min, ell_max)
    sproj.props.reltol = 1.0e-4

    xi_i = 1.0e1
    xi_f = 1.0e4

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
    # index = np.logical_and((z_a > 0.01), (z_a < 10.0))

    z_a = np.linspace(0.05, 0.5, num=2000)

    # print (xi_a)
    # print (z_a)

    for ell in range(ell_min, ell_min + 1):
        # limber = np.array(
        #     [
        #         ps.eval(cosmo, z, (0.5 + ell) / xi) / (xi * xi * RH_Mpc / cosmo.E(z))
        #         for (xi, z) in zip(xi_a, z_a)
        #     ]
        # )

        # plt.plot (xi_a[index], limber[index], label = r'$\ell$ = %d, limber' % ell)

        # for w_i in range (20):
        # w      = sproj.get_w (w_i)
        # Cell_s = sproj.peek_ell_spline (w_i, ell)
        # Cell   = [Cell_s.eval (math.log (xi)) for xi in xi_a]

        # for xi, Cell_i in zip (xi_a, Cell):
        # sCell_i = sproj.eval_Cell_xi1_xi2 (ell, xi / w, xi * w)
        # Cell_int_i = ps.sproj (cosmo, 1.0e-5, ell, 0.0, 0.0, xi / w, xi * w)

        # z_a    = np.array (z_a)
        # Cell   = np.array ([(gf.eval (cosmo, z)**2 * Cell_i) for (Cell_i, z) in zip (Cell, z_a)])

        # plt.plot (xi_a, Cell, label = r'$\ell$ = %d, $w$ = %f, full' % (ell, w))
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
        # print (np.diag (cov))
        std_ = np.sqrt(np.diag(cov))
        corr = cov / np.outer(std_, std_)
        # print (m)
        # print (corr)
        # plt.matshow (np.abs (corr), norm = LogNorm (vmin=1.0e-8, vmax = 1.0))
        # plt.matshow (corr, extent = np.array ([z_a[0], z_a[-1], z_a[0], z_a[-1]]),
        # origin = "lower")
        plt.matshow(
            corr,
            extent=np.array([z_a[0], z_a[-1], z_a[0], z_a[-1]]),
            origin="lower",
            norm=SymLogNorm(1.0e-8),
        )
        # plt.matshow (corr)
        # plt.matshow (corr, norm = SymLogNorm (1.0e-4))
        plt.title(rf"$C_{{{ell}}}(z_1, z_2)$")
        plt.colorbar()

    # plt.xticks(xi_a)

    # plt.xlabel (r'$\xi \; [\mathrm{Mpc}]$')
    # plt.ylabel (r'$C_\ell(\xi, \xi)$')
    # plt.legend (loc = "best")
    # plt.xscale ('log')
    # plt.yscale ('symlog', linthreshy=1.0e-8)
    # plt.xlim ([200.0, 4000.0])

    plt.show()


if __name__ == "__main__":
    test_spherical_projection()

#!/usr/bin/env python
#
# example_ps.py
#
# Mon May 22 16:00:00 2023
# Copyright  2023  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# example_ps.py
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

"""Example testing powerspectra."""

import math
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

from numcosmo_py import Nc, Ncm

matplotlib.rcParams.update({"font.size": 11})

#
#  Initializing the library objects, this must be called before
#  any other library function.
#
Ncm.cfg_init()


def test_ps() -> None:
    """Example testing powerspectra."""
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
    print("# Model parameters: ", end=" ")
    cosmo.params_log_all()
    print(f"# Omega_X0: {cosmo.E2Omega_de(0.0): 22.15g}")

    ps_cbe = Nc.PowspecMLCBE.new()
    ps_eh = Nc.PowspecMLTransfer.new(Nc.TransferFuncEH.new())

    # Redshift bounds
    z_min = 0.0
    z_max = 2.0
    zdiv = 0.49999999999

    # Mode bounds
    k_min = 1.0e-5
    k_max = 1.0e3

    nk = 2000
    nR = 2000
    Rh8 = 8.0 / cosmo.h()

    ps_cbe.set_kmin(k_min)
    ps_eh.set_kmin(k_min)

    ps_cbe.set_kmax(k_max)
    ps_eh.set_kmax(k_max)

    ps_cbe.require_zi(z_min)
    ps_cbe.require_zf(z_max)

    ps_eh.require_zi(z_min)
    ps_eh.require_zf(z_max)

    ps_eh.prepare(cosmo)
    ps_cbe.prepare(cosmo)

    for z in np.arange(z_min, z_max, (z_max - z_min) * zdiv, dtype=float):
        k_a = []
        Pk_eh_a = []
        Pk_cbe_a = []
        for lnk in np.arange(
            math.log(ps_cbe.props.kmin),
            math.log(ps_cbe.props.kmax),
            math.log(k_max / k_min) / nk,
        ):
            k = math.exp(lnk)
            k2 = k * k
            k3 = k2 * k
            k_a.append(k)
            Pk_cbe_a.append(k3 * ps_cbe.eval(cosmo, z, k))
            Pk_eh_a.append(k3 * ps_eh.eval(cosmo, z, k))

        plt.plot(k_a, Pk_cbe_a, label=f"CLASS $z = {z:.2f}$")
        plt.plot(k_a, Pk_eh_a, label=f"EH    $z = {z:.2f}$")

    plt.xlabel(r"$k \; [\mathrm{Mpc}^{-1}]$")
    plt.ylabel(r"$k^3P(k, z)$")
    plt.legend(loc="lower right")
    plt.xscale("log")
    plt.yscale("log")
    plt.title("Linear Matter Power Spectrum")
    plt.savefig("ps_cbe_eh.svg")
    plt.clf()

    for z in np.arange(z_min, z_max, (z_max - z_min) * zdiv * 0.5):
        k_a = []
        Pk_eh_a = []
        Pk_cbe_a = []
        for lnk in np.arange(
            math.log(ps_cbe.props.kmin),
            math.log(ps_cbe.props.kmax),
            math.log(k_max / k_min) / nk,
        ):
            k = math.exp(lnk)
            k2 = k * k
            k3 = k2 * k
            k_a.append(k)
            Pk_cbe_a.append(k3 * ps_cbe.eval(cosmo, z, k))
            Pk_eh_a.append(k3 * ps_eh.eval(cosmo, z, k))

        plt.plot(
            k_a,
            np.abs(1.0 - np.array(Pk_eh_a) / np.array(Pk_cbe_a)),
            label=f"$z = {z:.2f}$",
        )

    plt.xlabel(r"$k \; [\mathrm{Mpc}^{-1}]$")
    plt.ylabel(r"$1 - P_{EH}(k, z)/P_{CLASS}(k, z)$")
    plt.legend(loc="lower right")
    plt.xscale("log")
    plt.yscale("log")
    plt.title("Relative difference between the EH and CLASS PS")
    plt.savefig("ps_diff_cbe_eh.svg")
    plt.clf()

    #
    # Filtering
    #
    psf_cbe = Ncm.PowspecFilter.new(ps_cbe, Ncm.PowspecFilterType.TOPHAT)
    psf_eh = Ncm.PowspecFilter.new(ps_eh, Ncm.PowspecFilterType.TOPHAT)

    psf_cbe.set_best_lnr0()
    psf_eh.set_best_lnr0()

    psf_cbe.prepare(cosmo)
    psf_eh.prepare(cosmo)

    print(
        f"# CBE sigma8 = {psf_cbe.eval_sigma(0.0, Rh8): 20.15g}, "
        f"EH sigma8 = {psf_eh.eval_sigma(0.0, Rh8): 20.15g}"
    )
    print(f"# kmin {ps_cbe.props.kmin: 20.15g} kmax {ps_cbe.props.kmax: 20.15g}")
    print(f"# Rmin {psf_cbe.get_r_min(): 20.15g} Rmax {psf_cbe.get_r_max(): 20.15g}")

    lnRmin = math.log(psf_cbe.get_r_min())
    lnRmax = math.log(psf_cbe.get_r_max())

    # Variance of the matter density contrast
    for z in np.arange(z_min, z_max, (z_max - z_min) * zdiv, dtype=float):
        Rh_a = []
        sigma_eh_a = []
        sigma_cbe_a = []
        for lnR in np.arange(lnRmin, lnRmax, (lnRmax - lnRmin) / nR, dtype=float):
            R = math.exp(lnR)
            Rh = R * cosmo.h()

            Rh_a.append(Rh)
            sigma_cbe_a.append(psf_cbe.eval_sigma_lnr(z, lnR))
            sigma_eh_a.append(psf_eh.eval_sigma_lnr(z, lnR))

        plt.plot(Rh_a, sigma_cbe_a, label=f"CLASS $z = {z:.2f}$")
        plt.plot(Rh_a, sigma_eh_a, label=f"EH $z = {z:.2f}$")

    plt.xlabel(r"$R \; [h^{-1}\mathrm{Mpc}]$")
    plt.ylabel(r"$\sigma_M (R, z)$")
    plt.legend(loc="lower left")
    plt.xscale("log")
    plt.yscale("log")
    plt.title("Variance of the matter density contrast")
    plt.savefig("ps_var_cbe_eh.svg")
    plt.clf()

    # Derivative of the variance of the matter density contrast with respect to the
    # scale R
    for z in np.arange(z_min, z_max, (z_max - z_min) * zdiv):
        Rh_a = []
        dvar_eh_a = []
        dvar_cbe_a = []
        for lnR in np.arange(lnRmin, lnRmax, (lnRmax - lnRmin) / nR):
            R = math.exp(lnR)
            Rh = R * cosmo.h()

            Rh_a.append(Rh)
            dvar_cbe_a.append(psf_cbe.eval_dlnvar_dlnr(z, lnR))
            dvar_eh_a.append(psf_eh.eval_dlnvar_dlnr(z, lnR))

        plt.plot(Rh_a, dvar_cbe_a, label=f"CLASS $z = {z:.2f}$")
        plt.plot(Rh_a, dvar_eh_a, label=f"EH    $z = {z:.2f}$")

    plt.xlabel(r"$R \; [h^{-1}\mathrm{Mpc}]$")
    plt.ylabel(r"$\frac{\mathrm{d}\ln\sigma^2_M(R,z)}{\mathrm{d}\ln{}R}$")
    plt.title("Derivative of the matter density contrast")
    plt.legend(loc="lower left")
    plt.xscale("log")
    plt.savefig("ps_dvar_cbe_eh.svg")
    plt.clf()

    # Non-linear matter power spectrum: HALOFIT

    zmaxnl = 10.0
    z_max = zmaxnl
    pshf = Nc.PowspecMNLHaloFit.new(ps_cbe, zmaxnl, 1.0e-5)

    pshf.set_kmin(k_min)
    pshf.set_kmax(k_max)
    pshf.require_zi(z_min)
    pshf.require_zf(z_max)

    pshf.prepare(cosmo)

    for z in np.arange(z_min, z_max, (z_max - z_min) * zdiv):
        k_a = []
        Pk_cbe_a = []
        Pknl_cbe_a = []
        for lnk in np.arange(
            math.log(ps_cbe.props.kmin),
            math.log(ps_cbe.props.kmax),
            math.log(k_max / k_min) / nk,
        ):
            k = math.exp(lnk)
            k2 = k * k
            k3 = k2 * k
            k_a.append(k)
            Pk_cbe_a.append(k3 * ps_cbe.eval(cosmo, z, k))
            Pknl_cbe_a.append(k3 * pshf.eval(cosmo, z, k))

        plt.plot(k_a, Pk_cbe_a, label=f"CLASS $z = {z:.2f}$")
        plt.plot(k_a, Pknl_cbe_a, label=f"CLASS+HaloFit $z = {z:.2f}$")

    plt.xlabel(r"$k \; [\mathrm{Mpc}^{-1}]$")
    plt.ylabel(r"$k^3 P(k, z)$")
    plt.title("Linear and nonlinear matter power spectrum")
    plt.legend(loc="lower right")
    plt.xscale("log")
    plt.yscale("log")
    plt.savefig("ps_cbe_halofit.svg")
    plt.clf()

    psf_cbenl = Ncm.PowspecFilter.new(pshf, Ncm.PowspecFilterType.TOPHAT)
    psf_cbenl.set_best_lnr0()
    psf_cbenl.prepare(cosmo)

    print(f"# CBE+HaloFit sigma8 = {psf_cbenl.eval_sigma(0.0, Rh8): 20.15g}")


if __name__ == "__main__":
    test_ps()

#!/usr/bin/env python
#
# example_hiprim_tensor_modes.py
#
# Mon May 22 16:00:00 2023
# Copyright  2023  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# example_hiprim_tensor_modes.py
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

"""Example testing tensor modes in primordial models."""

import math
import numpy as np
import matplotlib.pyplot as plt

from numcosmo_py import Nc, Ncm

#
#  Initializing the library objects, this must be called before
#  any other library function.
#
Ncm.cfg_init()


def test_hiprim_tensor_modes() -> None:
    """Example testing tensor modes in primordial models."""

    #
    # Script parameters
    #
    # Maximum multipole
    lmax = 2500

    #
    # Creating a new instance of HIPrimPowerLaw
    #
    prim = Nc.HIPrimPowerLaw.new()
    r = 1.0
    prim.props.T_SA_ratio = r
    prim.props.n_T = -1.0 * r / 8.0

    #
    #  New CLASS backend precision object
    #  Lets also increase k_per_decade_primordial since we are
    #  dealing with a modified spectrum.
    #
    cbe_prec = Nc.CBEPrecision.new()
    cbe_prec.props.k_per_decade_primordial = 50.0

    #
    #  New CLASS backend object
    #
    cbe = Nc.CBE.prec_new(cbe_prec)

    #
    #  New CLASS backend object
    #
    Bcbe = Nc.HIPertBoltzmannCBE.full_new(cbe)
    Bcbe.set_TT_lmax(lmax)
    # Setting which CMB data to use
    Bcbe.set_target_Cls(Nc.DataCMBDataType.TT)
    # Setting if the lensed Cl's are going to be used or not.
    Bcbe.set_lensed_Cls(True)
    # Setting if the tensor contribution is going to be used or not.
    Bcbe.set_tensor(True)

    Bcbe.append_target_Cls(Nc.DataCMBDataType.TT)
    Bcbe.append_target_Cls(Nc.DataCMBDataType.TE)
    Bcbe.append_target_Cls(Nc.DataCMBDataType.EE)
    Bcbe.append_target_Cls(Nc.DataCMBDataType.BB)

    Bcbe.set_TT_lmax(lmax)
    Bcbe.set_TE_lmax(lmax)
    Bcbe.set_EE_lmax(lmax)
    Bcbe.set_BB_lmax(30)

    #
    #  New homogeneous and isotropic cosmological model NcHICosmoDEXcdm
    #
    cosmo = Nc.HICosmo.new_from_name(Nc.HICosmo, "NcHICosmoDEXcdm")
    cosmo.omega_x2omega_k()
    cosmo.param_set_by_name("Omegak", 0.0)

    #
    #  New homogeneous and isotropic reionization object
    #
    reion = Nc.HIReionCamb.new()

    #
    # Adding submodels to the main cosmological model.
    #
    cosmo.add_submodel(reion)
    cosmo.add_submodel(prim)

    #
    # Preparing the Class backend object
    #
    fact = 1.0e200
    ln10e10ASA = prim.props.ln10e10ASA
    T_SA_ratio = prim.props.T_SA_ratio

    prim.props.ln10e10ASA = prim.props.ln10e10ASA + math.log(1.0 / fact)
    prim.props.T_SA_ratio = prim.props.T_SA_ratio * fact
    Bcbe.prepare(cosmo)

    Cls1_TT = Ncm.Vector.new(lmax + 1)
    Cls2_TT = Ncm.Vector.new(lmax + 1)

    Cls1_TE = Ncm.Vector.new(lmax + 1)
    Cls2_TE = Ncm.Vector.new(lmax + 1)

    Cls1_EE = Ncm.Vector.new(lmax + 1)
    Cls2_EE = Ncm.Vector.new(lmax + 1)

    Cls1_BB = Ncm.Vector.new(31)
    Cls2_BB = Ncm.Vector.new(31)

    Bcbe.get_TT_Cls(Cls1_TT)
    Bcbe.get_TE_Cls(Cls1_TE)
    Bcbe.get_EE_Cls(Cls1_EE)
    Bcbe.get_BB_Cls(Cls1_BB)

    prim.props.ln10e10ASA = ln10e10ASA
    prim.props.T_SA_ratio = T_SA_ratio

    prim.props.T_SA_ratio = 1.0 / fact
    Bcbe.prepare(cosmo)

    Bcbe.get_TT_Cls(Cls2_TT)
    Bcbe.get_TE_Cls(Cls2_TE)
    Bcbe.get_EE_Cls(Cls2_EE)
    Bcbe.get_BB_Cls(Cls2_BB)

    Cls1_TT_a = Cls1_TT.dup_array()
    Cls1_TE_a = Cls1_TE.dup_array()
    Cls1_EE_a = Cls1_EE.dup_array()
    Cls1_BB_a = Cls1_BB.dup_array()

    Cls2_TT_a = Cls2_TT.dup_array()
    Cls2_TE_a = Cls2_TE.dup_array()
    Cls2_EE_a = Cls2_EE.dup_array()
    Cls2_BB_a = Cls2_BB.dup_array()

    ell = np.array(list(range(2, lmax + 1)))
    ell_BB = np.array(list(range(2, 31)))

    Cls1_TT_a = ell * (ell + 1.0) * np.array(Cls1_TT_a[2:])
    Cls1_TE_a = ell * (ell + 1.0) * np.array(Cls1_TE_a[2:])
    Cls1_EE_a = ell * (ell + 1.0) * np.array(Cls1_EE_a[2:])
    Cls1_BB_a = ell_BB * (ell_BB + 1.0) * np.array(Cls1_BB_a[2:])

    Cls2_TT_a = ell * (ell + 1.0) * np.array(Cls2_TT_a[2:])
    Cls2_TE_a = ell * (ell + 1.0) * np.array(Cls2_TE_a[2:])
    Cls2_EE_a = ell * (ell + 1.0) * np.array(Cls2_EE_a[2:])
    Cls2_BB_a = ell_BB * (ell_BB + 1.0) * np.array(Cls2_BB_a[2:])

    print("TT")
    print(Cls1_TT_a[:28])
    print(Cls2_TT_a[:28])

    print("TE")
    print(Cls1_TE_a[:28])
    print(Cls2_TE_a[:28])

    print("EE")
    print(Cls1_EE_a[:28])
    print(Cls2_EE_a[:28])

    print("BB")
    print(Cls1_EE_a[:28])
    print(Cls2_EE_a[:28])

    #
    #  Ploting the TT angular power spcetrum
    #

    plt.title(r"With and without tensor contribution to $C_\ell^\mathrm{TT}$")
    plt.plot(ell[:28], Cls1_TT_a[:28], "r", label="tensor")
    plt.plot(ell[:28], Cls2_TT_a[:28], "b--", label="scalar")

    plt.xlabel(r"$\ell$")
    plt.ylabel(r"$C_\ell$")
    plt.legend(loc="best")
    # plt.xscale ('log')
    # plt.yscale ('symlog')

    plt.savefig("hiprim_tensor_Cls_TT.pdf")
    plt.clf()

    plt.title(r"With and without tensor contribution to $C_\ell^\mathrm{TE}$")
    plt.plot(ell[:28], Cls1_TE_a[:28], "r", label="tensor")
    plt.plot(ell[:28], Cls2_TE_a[:28], "b--", label="scalar")

    plt.xlabel(r"$\ell$")
    plt.ylabel(r"$C_\ell$")
    plt.legend(loc="best")
    # plt.xscale ('log')
    # plt.yscale ('symlog')

    plt.savefig("hiprim_tensor_Cls_TE.pdf")
    plt.clf()

    plt.title(r"With and without tensor contribution to $C_\ell^\mathrm{EE}$")
    plt.plot(ell[:28], Cls1_EE_a[:28], "r", label="tensor")
    plt.plot(ell[:28], Cls2_EE_a[:28], "b--", label="scalar")

    plt.xlabel(r"$\ell$")
    plt.ylabel(r"$C_\ell$")
    plt.legend(loc="best")
    # plt.xscale ('log')
    # plt.yscale ('symlog')

    plt.savefig("hiprim_tensor_Cls_EE.pdf")
    plt.clf()

    plt.title(r"With and without tensor contribution to $C_\ell^\mathrm{BB}$")
    plt.plot(ell_BB[:28], Cls1_BB_a[:28], "r", label="tensor")
    plt.plot(ell_BB[:28], Cls2_BB_a[:28], "b--", label="scalar")

    plt.xlabel(r"$\ell$")
    plt.ylabel(r"$C_\ell$")
    plt.legend(loc="best")
    # plt.xscale ('log')
    # plt.yscale ('symlog')

    plt.savefig("hiprim_tensor_Cls_BB.pdf")
    plt.clf()

    theta = np.linspace(math.pi / 3.0, math.pi, 1000)
    smap = Ncm.SphereMap.new(128)
    smap.set_lmax(lmax)

    smap.set_Cls(Cls1_TT)
    s = smap.calc_Ctheta(1.0e-7)
    Ctheta1 = np.vectorize(s.eval)(theta)

    smap.set_Cls(Cls2_TT)
    s = smap.calc_Ctheta(1.0e-7)
    Ctheta2 = np.vectorize(s.eval)(theta)

    plt.title(r"$C^\mathrm{TT}(\theta)$")
    plt.plot(theta, Ctheta1, "r", lw=0.1, label=r"$C(\theta)^\mathrm{tensor}$")
    plt.plot(theta, Ctheta2, "b", lw=0.1, label=r"$C(\theta)^\mathrm{scalar}$")

    plt.xlabel(r"$\theta$")
    plt.ylabel(r"$C(\theta)$")
    plt.legend(loc="best")
    # plt.xscale ('log')
    # plt.yscale ('symlog')

    plt.savefig("hiprim_tensor_Ctheta_TT.pdf")
    plt.clf()


if __name__ == "__main__":
    test_hiprim_tensor_modes()

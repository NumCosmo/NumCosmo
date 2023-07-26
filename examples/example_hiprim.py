#!/usr/bin/env python
#
# example_hiprim.py
#
# Mon May 22 16:00:00 2023
# Copyright  2023  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# example_hiprim.py
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

"""Example testing primordial models."""

import numpy as np
import matplotlib.pyplot as plt

from py_hiprim_example import PyHIPrimExample
from numcosmo_py import Nc, Ncm


#
#  Initializing the library objects, this must be called before
#  any other library function.
#
Ncm.cfg_init()


def test_hiprim() -> None:
    """Example testing primordial models."""

    #
    # Script parameters
    #
    # maximum multipole
    lmax = 2500

    #
    # Creating a new instance of PyHIPrimExample
    #
    prim = PyHIPrimExample()

    print(f"# As        = {prim.props.As}")
    print(f"# P (k = 1) = {prim.SA_powspec_k(1.0)}")
    print(f"# (a, b, c) = ({prim.props.a}, {prim.props.b}, {prim.props.c})")

    #
    #  New CLASS backend precision object
    #  Let's also increase k_per_decade_primordial since we are
    #  dealing with a modified spectrum.
    #
    cbe_prec = Nc.CBEPrecision.new()
    cbe_prec.props.k_per_decade_primordial = 50.0
    cbe_prec.props.tight_coupling_approximation = 0

    #
    #  New CLASS backend object
    #
    cbe = Nc.CBE.prec_new(cbe_prec)

    Bcbe = Nc.HIPertBoltzmannCBE.full_new(cbe)
    Bcbe.set_TT_lmax(lmax)
    # Setting which CMB data to use
    Bcbe.set_target_Cls(Nc.DataCMBDataType.TT)
    # Setting if the lensed Cl's are going to be used or not.
    Bcbe.set_lensed_Cls(True)

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
    Bcbe.prepare(cosmo)

    Cls1 = Ncm.Vector.new(lmax + 1)
    Cls2 = Ncm.Vector.new(lmax + 1)

    Bcbe.get_TT_Cls(Cls1)

    prim.props.a = 0
    Bcbe.prepare(cosmo)
    Bcbe.get_TT_Cls(Cls2)

    Cls1_a = np.array(Cls1.dup_array())
    Cls2_a = np.array(Cls2.dup_array())

    Cls1_a = np.array(Cls1_a[2:])
    Cls2_a = np.array(Cls2_a[2:])

    ell = np.array(list(range(2, lmax + 1)))

    Cls1_a = ell * (ell + 1.0) * Cls1_a
    Cls2_a = ell * (ell + 1.0) * Cls2_a

    #
    #  Ploting the TT angular power spcetrum
    #

    plt.title(r"Modified and non-modified $C_\ell$")
    plt.xscale("log")
    plt.plot(ell, Cls1_a, "r", label="Modified")
    plt.plot(ell, Cls2_a, "b--", label="Non-modified")

    plt.xlabel(r"$\ell$")
    plt.ylabel(r"$C_\ell$")
    plt.legend(loc=2)

    plt.savefig("hiprim_Cls.svg")

    plt.clf()


if __name__ == "__main__":
    test_hiprim()

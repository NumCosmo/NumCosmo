#
# two_point.py
#
# Tue Mar 4 11:14:17 2025
# Copyright  2025  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# two_point.py
# Copyright (C) 2025 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

"""Two-point correlation functions."""

import numpy as np
import pyccl

from numcosmo_py import Ncm
from numcosmo_py.cosmology import Cosmology


def compute_kernel(
    tracer: pyccl.Tracer, cosmology: Cosmology, ell: float
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Compute the kernel for a given tracer."""
    cosmo = cosmology.cosmo
    dist = cosmology.dist
    Wchi_list, chi_list = tracer.get_kernel()
    assert chi_list is not None
    chi_a = np.array(chi_list[0])
    for chi_a_i, Wchi_a_i in zip(chi_list, Wchi_list):
        s = Ncm.Spline.new_array(
            Ncm.SplineCubicNotaknot.new(), chi_a_i.tolist(), Wchi_a_i.tolist(), True
        )
        Wchi_a_i[:] = np.array([s.eval(chi) for chi in chi_a])

    RH_Mpc = cosmo.RH_Mpc()
    nu = ell + 0.5

    bessel_factors_list = []

    for der_bessel in tracer.get_bessel_derivative():
        match der_bessel:
            case 0:
                bessel_factors_list.append(1.0)
            case -1:
                bessel_factors_list.append(1.0 / nu**2)
            case _:
                raise ValueError(f"Invalid Bessel derivative {der_bessel}")

    z_a = np.array([dist.inv_comoving(cosmo, chi / RH_Mpc) for chi in chi_a])
    H_Mpc_a = np.array([cosmo.E(z) / RH_Mpc for z in z_a])
    a_array = 1.0 / (1.0 + np.array(z_a))
    transfers_list = tracer.get_transfer(0.0, a_array)
    ell_factors_list = tracer.get_f_ell(ell)

    Wtotal = np.zeros_like(chi_a)
    for Wchi_a, transfer, ell_factor, bessel_factor in zip(
        Wchi_list, transfers_list, ell_factors_list, bessel_factors_list
    ):
        assert Wchi_a is not None
        assert transfer is not None
        Wtotal += Wchi_a * transfer * ell_factor * bessel_factor

    return z_a[1:-1], chi_a[1:-1], H_Mpc_a[1:-1], Wtotal[1:-1]

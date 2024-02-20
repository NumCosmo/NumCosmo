#
# numcosmo_class.py
#
# Tue Feb 20 11:08:00 2024
# Copyright  2024  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# numcosmo_class.py
# Copyright (C) 2024 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

"""Functions to compute cosmological quantities using numcosmo."""


import numpy as np

from ... import Nc


def get_numcosmo(
    z_arr: np.ndarray, cosmo: Nc.HICosmo
) -> tuple[float, np.ndarray, np.ndarray, np.ndarray, Nc.PowspecML]:
    """Compute using numcosmo the comoving distance, the derivative of the
    comoving distance and the growth factor."""

    dist = Nc.Distance.new(2.0)
    dist.prepare(cosmo)

    growth_func = Nc.GrowthFunc.new()
    growth_func.prepare(cosmo)

    tf = Nc.TransferFuncEH.new()
    psml = Nc.PowspecMLTransfer.new(tf)
    psml.require_kmin(1.0e-6)
    psml.require_kmax(1.0e3)
    psml.prepare(cosmo)

    h = cosmo.h()
    comov_dist = cosmo.RH_Mpc() * np.array([dist.comoving(cosmo, z) for z in z_arr])
    dcomov_dist = cosmo.RH_Mpc() * np.array([1.0 / cosmo.E(z) for z in z_arr])
    growth = np.array([growth_func.eval(cosmo, z) for z in z_arr])

    return h, comov_dist, dcomov_dist, growth, psml

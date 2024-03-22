#!/usr/bin/env python
#
# example_zeta_exp_potential.py
#
# Mon May 22 16:00:00 2023
# Copyright  2023  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# example_zeta_exp_potential.py
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

"""Example of scalar perturbations on exponential pontential cosmology."""

import math
import numpy as np
import matplotlib.pyplot as plt

from numcosmo_py import Nc, Ncm


#
#  Initializing the library objects, this must be called before
#  any other library function.
#
Ncm.cfg_init()


def test_zeta_exp_potential(new_config: bool = True) -> None:
    """Test scalar perturbations on exponential pontential cosmology."""

    adiab = Nc.HIPertAdiab.new()
    gw = Nc.HIPertGW.new()
    Vexp = Nc.HICosmoVexp.new()

    if new_config:  # Boa!
        Vexp.props.alphab = +1.0e-1
        Vexp.props.sigmaphi = +0.8
        Vexp.props.dphi = +0.5
        #  Vexp.props.xb       = 4.8e37
        Vexp.props.xb = 1.0e37
        Vexp.props.OmegaL = 1.0
        Vexp.props.Omegac = 1.0
        Vexp.props.H0 = 67.8
    else:
        Vexp.props.alphab = +8.3163e-2
        Vexp.props.sigmaphi = +9.0
        Vexp.props.dphi = -9.0e-4
        Vexp.props.xb = 2.0e36
        Vexp.props.OmegaL = 1.0
        Vexp.props.Omegac = 1.0
        Vexp.props.H0 = 67.8

    k = 0.05
    tc = +Vexp.tau_xe(1.0e15)
    reltol = 1.0e-14

    adiab.set_ti(Vexp.tau_min())
    adiab.set_tf(tc)
    adiab.set_k(k)
    adiab.set_reltol(reltol)

    gw.set_ti(Vexp.tau_min())
    gw.set_tf(tc)
    gw.set_k(k)
    gw.set_reltol(reltol)

    # If = lambda x: gw.eval_nu (Vexp, x, k) / gw.eval_mnu (Vexp, x, k)
    # print (integrate.quad (If, -2.0, 2.0))
    # for n in np.logspace (-40, -1, 10000):
    #  print (n, If (n), If (-n))
    # exit ()

    print("# Preparing ADIAB")
    found, adiab_ti = adiab.find_adiab_time_limit(Vexp, -20.0, -1.0e-2, 1.0e-6)
    assert found, "Cannot find initial time"
    adiab.set_init_cond_adiab(Vexp, adiab_ti)
    adiab.prepare(Vexp)

    # exit ()

    print("# Preparing GW")
    found, gw_ti = gw.find_adiab_time_limit(Vexp, -20.0, -1.0e-2, 1.0e-8)
    assert found, "Cannot find initial time"
    gw.set_init_cond_adiab(Vexp, gw_ti)
    gw.prepare(Vexp)

    # (t0, t1) = adiab.get_t0_t1 (Vexp)
    t0 = gw.get_ti()
    t1 = gw.get_tf()

    print(f"# BACKG (t0, t1) = ({Vexp.tau_min(): 21.15e}, {Vexp.tau_max(): 21.15e})")
    print(f"# ADIAB (t0, t1) = ({t0: 21.15e}, {t1: 21.15e})")

    (Delta_zeta_c, _Delta_Pzeta_c) = adiab.eval_powspec_at(Vexp, tc)
    (Delta_h_c, _Delta_Ph_c) = gw.eval_powspec_at(Vexp, tc)

    print(f"# Time of x_e = 10^15:     tau_c     = {tc: 21.15f}")
    print(f"# Power spectrum at tau_c: PS_ADIAB  = {Delta_zeta_c: 21.15e}")
    print(f"# Power spectrum at tau_c: PS_GW     = {Delta_h_c: 21.15e}")
    print(
        f"# Power spectrum at tau_c: r         = {2.0 * Delta_h_c / Delta_zeta_c: 21.15e}"
    )
    
    tau_a = np.linspace(-5.0, t1, 100000)
    Delta_zeta = np.array([adiab.eval_powspec_at(Vexp, t) for t in tau_a])
    Delta_h = np.array([gw.eval_powspec_at(Vexp, t) for t in tau_a])
    mylw = 1.0

    
    plt.plot(tau_a, np.abs(Delta_zeta[:,0]), lw=mylw, label=r"$\Delta_{\zeta}a$")
    plt.plot(tau_a, Delta_h[:,0], lw=mylw, label=r"$\Delta_h$")

    plt.grid(which="both", linestyle=":", color="0.75", linewidth=0.5)
    plt.legend(loc="upper left")

    # plt.xscale('symlog', linthreshx=1.0e-5)
    plt.yscale("log")
    # plt.yscale('log', linthreshy=1.0e-20)

    plt.show()
    plt.clf()


if __name__ == "__main__":
    test_zeta_exp_potential(new_config = False)

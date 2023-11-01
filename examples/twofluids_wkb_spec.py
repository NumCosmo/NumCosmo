#!/usr/bin/env python
#
# two_fluids_wkb_spec.py
#
# Mon May 22 16:00:00 2023
# Copyright  2023  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# two_fluids_wkb_spec.py
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

"""Compute WKB approximation for the two-fluids model spectrum."""

import sys
import math
import numpy as np

from tqdm import tqdm

import matplotlib.pyplot as plt

from numcosmo_py import Nc, Ncm


#
#  Initializing the library objects, this must be called before
#  any other library function.
#
Ncm.cfg_init()


def test_two_fluids_wkb_spec() -> None:
    """Compute WKB approximation for the two-fluids model spectrum."""
    
    #
    #  New homogeneous and isotropic cosmological model NcHICosmoQGRW
    #
    cosmo = Nc.HICosmoQGRW()

    if len(sys.argv) != 2:
        print(f"{sys.argv[0]} w")
        sys.exit(0)

    w = float(sys.argv[1])
    prec = 1.0e-7

    cosmo.props.w = w
    cosmo.props.Omegar = (1.0e-6) * 1.0
    cosmo.props.Omegaw = (1.0 - 1.0e-6) * 1.0
    cosmo.props.xb = 1.0e30

    pert = Nc.HIPertTwoFluids.new()

    pert.props.reltol = prec
    # pert.set_stiff_solver (True)

    lnki = math.log(1.0e-3)
    lnkf = math.log(1.0e3)
    lnk_a = np.linspace(lnki, lnkf, 20)

    ci = Ncm.Vector.new(8)

    k_a = []
    Ps_zeta1 = []
    Ps_S1 = []
    Ps_zeta2 = []
    Ps_S2 = []
    Ps_Pzeta1 = []
    Ps_PS1 = []
    Ps_Pzeta2 = []
    Ps_PS2 = []

    out_file = open("twofluids_spectrum_{w}.dat", "w", encoding="utf-8")

    start_alpha1 = 1.0e-10
    start_alpha2 = 1.0e-14

    for lnk in tqdm(lnk_a):
        k = math.exp(lnk)
        pert.set_mode_k(k)
        k_a.append(k)

        alphaf = cosmo.abs_alpha(1.0e10)


        alphai = -cosmo.abs_alpha(start_alpha1 * k**2)
        pert.get_init_cond_zetaS(cosmo, alphai, 1, 0.25 * math.pi, ci)
        pert.set_init_cond(cosmo, alphai, 1, False, ci)
        print(f"# Mode 1 k {k: 21.15e}, state module {pert.get_state_mod():f}")
        pert.evolve(cosmo, alphaf)  
        v, _alphac = pert.peek_state(cosmo)
        print ("# Evolving mode %e from %f to %f" % (k, alphai, alphaf))

        Delta_zeta1 = (
            k**3
            * math.hypot(
                v.get(Nc.HIPertITwoFluidsVars.ZETA_R),
                v.get(Nc.HIPertITwoFluidsVars.ZETA_I),
            )
            ** 2
            / (2.0 * math.pi**2 * cosmo.RH_planck() ** 2)
        )
        Delta_S1 = (
            k**3
            * math.hypot(
                v.get(Nc.HIPertITwoFluidsVars.S_R), v.get(Nc.HIPertITwoFluidsVars.S_I)
            )
            ** 2
            / (2.0 * math.pi**2 * cosmo.RH_planck() ** 2)
        )
        Delta_Pzeta1 = (
            k**3
            * math.hypot(
                v.get(Nc.HIPertITwoFluidsVars.PZETA_R),
                v.get(Nc.HIPertITwoFluidsVars.PZETA_I),
            )
            ** 2
            / (2.0 * math.pi**2 * cosmo.RH_planck() ** 2)
        )
        Delta_PS1 = (
            k**3
            * math.hypot(
                v.get(Nc.HIPertITwoFluidsVars.PS_R), v.get(Nc.HIPertITwoFluidsVars.PS_I)
            )
            ** 2
            / (2.0 * math.pi**2 * cosmo.RH_planck() ** 2)
        )

        Ps_zeta1.append(Delta_zeta1)
        Ps_S1.append(Delta_S1)
        Ps_Pzeta1.append(Delta_Pzeta1)
        Ps_PS1.append(Delta_PS1)

        alphai = -cosmo.abs_alpha(start_alpha2 * k**2)
        pert.get_init_cond_zetaS(cosmo, alphai, 2, 0.25 * math.pi, ci)
        pert.set_init_cond(cosmo, alphai, 2, False, ci)

        print(f"# Mode 2 k {k: 21.15e}, state module {pert.get_state_mod():f}")

        pert.evolve(cosmo, alphaf)
        v, _alphac = pert.peek_state(cosmo)

        Delta_zeta2 = (
            k**3
            * math.hypot(
                v.get(Nc.HIPertITwoFluidsVars.ZETA_R),
                v.get(Nc.HIPertITwoFluidsVars.ZETA_I),
            )
            ** 2
            / (2.0 * math.pi**2 * cosmo.RH_planck() ** 2)
        )
        Delta_S2 = (
            k**3
            * math.hypot(
                v.get(Nc.HIPertITwoFluidsVars.S_R), v.get(Nc.HIPertITwoFluidsVars.S_I)
            )
            ** 2
            / (2.0 * math.pi**2 * cosmo.RH_planck() ** 2)
        )
        Delta_Pzeta2 = (
            k**3
            * math.hypot(
                v.get(Nc.HIPertITwoFluidsVars.PZETA_R),
                v.get(Nc.HIPertITwoFluidsVars.PZETA_I),
            )
            ** 2
            / (2.0 * math.pi**2 * cosmo.RH_planck() ** 2)
        )
        Delta_PS2 = (
            k**3
            * math.hypot(
                v.get(Nc.HIPertITwoFluidsVars.PS_R), v.get(Nc.HIPertITwoFluidsVars.PS_I)
            )
            ** 2
            / (2.0 * math.pi**2 * cosmo.RH_planck() ** 2)
        )

        Ps_zeta2.append(Delta_zeta2)
        Ps_S2.append(Delta_S2)
        Ps_Pzeta2.append(Delta_Pzeta2)
        Ps_PS2.append(Delta_PS2)

        out_file.write(
            f"{k: 20.15e} {Delta_zeta1: 20.15e} {Delta_zeta2: 20.15e} {Delta_S1: 20.15e} "
            f"{Delta_S2: 20.15e} {Delta_Pzeta1: 20.15e} {Delta_Pzeta2: 20.15e} "
            f"{Delta_PS1: 20.15e} {Delta_PS2: 20.15e}\n"
        )
        out_file.flush()

    out_file.close()

    plt.plot(k_a, Ps_zeta1, label=r"$P^1_\zeta$")
    plt.plot(k_a, Ps_S1, label=r"$P^1_S$")
    plt.plot(k_a, Ps_zeta2, label=r"$P^2_\zeta$")
    plt.plot(k_a, Ps_S2, label=r"$P^2_S$")

    plt.grid()
    plt.legend(loc="upper left")
    plt.xscale("log")
    plt.yscale("log")

    plt.show()
    plt.clf()
    
test_two_fluids_wkb_spec()


'''
import pandas as pd

data = pd.read_csv('twofluids_spectrum_{w}.dat', header=None, delim_whitespace= True)
lnk = data[0]
Ps1zeta = data[1]
Ps1S = data[2]
Ps2zeta = data[3]
Ps2S = data[4]

ns1zeta = np.polyfit(lnk,Ps1zeta,1)
ns1S = np.polyfit(lnk,Ps1S,1)
ns2zeta = np.polyfit(lnk,Ps2zeta,1)
ns2S = np.polyfit(lnk,Ps2S,1)

print(ns1zeta)
print(ns1S)
print(ns2zeta)
print(ns2S)
'''
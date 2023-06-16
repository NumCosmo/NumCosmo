#!/usr/bin/env python
#
# two_fluids_wkb_mode.py
#
# Mon May 22 16:00:00 2023
# Copyright  2023  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# two_fluids_wkb_mode.py
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

"""Compute WKB approximation for the two-fluids model."""
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


def test_two_fluids_wkb_mode(mode: int = Nc.HIPertTwoFluidsCross.MODE1SUB) -> None:
    """Compute WKB approximation for the two-fluids model."""

    #
    #  New homogeneous and isotropic cosmological model NcHICosmoQGRW
    #
    cosmo = Nc.HICosmoQGRW()

    if len(sys.argv) != 3:
        print("twofluids_wkb_mode.py mode_k w")
        sys.exit(0)

    w = float(sys.argv[2])
    prec = 1.0e-6
    mode_k = float(sys.argv[1])

    cosmo.props.w = w
    cosmo.props.Omegar = 2.0 * (1.0e-5)
    cosmo.props.Omegaw = 2.0 * (1.0 - 1.0e-5)
    cosmo.props.xb = 1.0e30

    pert = Nc.HIPertTwoFluids.new()

    pert.props.reltol = prec
    pert.set_mode_k(mode_k)

    cross_size = 1.0e-5
    alpha_try = -cosmo.abs_alpha(1.0e-12 * mode_k**2)

    if mode == Nc.HIPertTwoFluidsCross.MODE1MAIN:
        alphai = pert.get_cross_time(
            cosmo, Nc.HIPertTwoFluidsCross.MODE1MAIN, alpha_try, cross_size
        )
    elif mode == Nc.HIPertTwoFluidsCross.MODE1SUB:
        alphai = pert.get_cross_time(
            cosmo, Nc.HIPertTwoFluidsCross.MODE1SUB, alpha_try, cross_size
        )
    elif mode == Nc.HIPertTwoFluidsCross.MODE2MAIN:
        alphai = pert.get_cross_time(
            cosmo, Nc.HIPertTwoFluidsCross.MODE2MAIN, alpha_try, cross_size
        )
    elif mode == Nc.HIPertTwoFluidsCross.MODE2SUB:
        alphai = pert.get_cross_time(
            cosmo, Nc.HIPertTwoFluidsCross.MODE2SUB, alpha_try, cross_size
        )
    else:
        raise ValueError("Invalid mode")

    alphaf = +cosmo.abs_alpha(1.0e20)

    print("# Mode k = {mode_k: 21.15g}")

    pert.set_stiff_solver(False)

    alpha_a = []
    gammabar11_a = []
    gammabar22_a = []
    gammabar12_a = []
    taubar12_a = []
    nu1_a = []
    nu2_a = []

    for alpha in np.linspace(alphai, alphaf, 10000):
        eom = pert.eom(cosmo, alpha)
        alpha_a.append(alpha)

        gammabar11_a.append(math.fabs(eom.gammabar11))
        gammabar22_a.append(math.fabs(eom.gammabar22))
        gammabar12_a.append(math.fabs(eom.gammabar12))
        taubar12_a.append(math.fabs(eom.taubar))
        nu1_a.append(eom.nu1)
        nu2_a.append(eom.nu2)

    print(
        "# Calculating mode 1, initial time {alphai: 20.15f} [{cosmo.x_alpha(alphai):8.2e}]: "
    )

    ci = Ncm.Vector.new(8)

    pert.get_init_cond_zetaS(cosmo, alphai, 1, 0.25 * math.pi, ci)
    pert.set_init_cond(cosmo, alphai, 3, False, ci)

    Ps_zeta1 = []
    Ps_S1 = []
    Ps_Pzeta1 = []
    Ps_PS1 = []

    Ps_zeta1.append(
        math.hypot(
            ci.get(Nc.HIPertITwoFluidsVars.ZETA_R),
            0.0 * ci.get(Nc.HIPertITwoFluidsVars.ZETA_I),
        )
        ** 2
    )
    Ps_S1.append(
        math.hypot(
            ci.get(Nc.HIPertITwoFluidsVars.S_R),
            0.0 * ci.get(Nc.HIPertITwoFluidsVars.S_I),
        )
        ** 2
    )
    Ps_Pzeta1.append(
        math.hypot(
            ci.get(Nc.HIPertITwoFluidsVars.PZETA_R),
            0.0 * ci.get(Nc.HIPertITwoFluidsVars.PZETA_I),
        )
        ** 2
    )
    Ps_PS1.append(
        math.hypot(
            ci.get(Nc.HIPertITwoFluidsVars.PS_R),
            0.0 * ci.get(Nc.HIPertITwoFluidsVars.PS_I),
        )
        ** 2
    )

    for alpha in tqdm(alpha_a[1:]):
        # for alpha in alpha_a[1:]:
        pert.evolve(cosmo, alpha)
        v, _alphac = pert.peek_state(cosmo)
        Ps_zeta1.append(
            math.hypot(
                v.get(Nc.HIPertITwoFluidsVars.ZETA_R),
                0.0 * v.get(Nc.HIPertITwoFluidsVars.ZETA_I),
            )
            ** 2
        )
        Ps_S1.append(
            math.hypot(
                v.get(Nc.HIPertITwoFluidsVars.S_R),
                0.0 * v.get(Nc.HIPertITwoFluidsVars.S_I),
            )
            ** 2
        )
        Ps_Pzeta1.append(
            math.hypot(
                v.get(Nc.HIPertITwoFluidsVars.PZETA_R),
                0.0 * v.get(Nc.HIPertITwoFluidsVars.PZETA_I),
            )
            ** 2
        )
        Ps_PS1.append(
            math.hypot(
                v.get(Nc.HIPertITwoFluidsVars.PS_R),
                0.0 * v.get(Nc.HIPertITwoFluidsVars.PS_I),
            )
            ** 2
        )
        # print ("norm = % 8.2e % 21.15f [%8.2e]" % (pert.get_state_mod () - 1.0,
        # alpha, cosmo.x_alpha (alpha)))

    # alphai = alpha_mode2main

    # print ("# Calculating mode 2, initial time % 20.15f [%8.2e]: " % (alphai,
    # cosmo.x_alpha (alphai)))

    # pert.get_init_cond_zetaS (cosmo, alphai, 2, 0.25 * math.pi, ci)
    # pert.set_init_cond (cosmo, alphai, 2, False, ci)

    # Ps_zeta2  = []
    # Ps_S2     = []
    # Ps_Pzeta2 = []
    # Ps_PS2    = []

    # alpha_a_pre = np.linspace (alphai, alpha_mode1main, 1000, endpoint = False)

    # for alpha in tqdm (alpha_a_pre[1:]):
    #   pert.evolve (cosmo, alpha)

    # for alpha in tqdm (alpha_a):
    # #for alpha in alpha_a:
    #   pert.evolve (cosmo, alpha)
    #   v, alphac = pert.peek_state (cosmo)
    #   Ps_zeta2.append  (math.hypot (v.get (Nc.HIPertITwoFluidsVars.ZETA_R),
    # v.get (Nc.HIPertITwoFluidsVars.ZETA_I))**2)
    #   Ps_S2.append     (math.hypot (v.get (Nc.HIPertITwoFluidsVars.S_R),
    # v.get (Nc.HIPertITwoFluidsVars.S_I))**2)
    #   Ps_Pzeta2.append (math.hypot (v.get (Nc.HIPertITwoFluidsVars.PZETA_R),
    # v.get (Nc.HIPertITwoFluidsVars.PZETA_I))**2)
    #   Ps_PS2.append    (math.hypot (v.get (Nc.HIPertITwoFluidsVars.PS_R),
    # v.get (Nc.HIPertITwoFluidsVars.PS_I))**2)

    # plt.plot (alpha_a, gammabar11_a, label = r'$\bar\gamma_{11}$')
    # plt.plot (alpha_a, gammabar22_a, label = r'$\bar\gamma_{22}$')
    # plt.plot (alpha_a, gammabar12_a, label = r'$\bar\gamma_{12}$')
    # plt.plot (alpha_a, taubar12_a,   label = r'$\bar\tau_{12}$')
    # plt.plot (alpha_a, nu1_a,        label = r'$\nu_{1}$')
    # plt.plot (alpha_a, nu2_a,        label = r'$\nu_{2}$')

    plt.plot(alpha_a, Ps_zeta1, label=r"$P^1_\zeta$")
    plt.plot(alpha_a, Ps_S1, label=r"$P^1_S$")
    # plt.plot (alpha_a, Ps_zeta2,     label = r'$P^2_\zeta$')
    # plt.plot (alpha_a, Ps_S2,        label = r'$P^2_S$')

    plt.plot(alpha_a, Ps_Pzeta1, label=r"$P^1_{P_\zeta}$")
    plt.plot(alpha_a, Ps_PS1, label=r"$P^1_{P_S}$")
    # plt.plot (alpha_a, Ps_Pzeta2,    label = r'$P^2_{P_\zeta}$')
    # plt.plot (alpha_a, Ps_PS2,       label = r'$P^2_{P_S}$')

    plt.grid()
    plt.legend(loc="upper left")
    # plt.xscale('log')
    plt.yscale("log")

    Delta_zeta1 = (
        mode_k**3 * Ps_zeta1.pop() / (2.0 * math.pi**2 * cosmo.RH_planck() ** 2)
    )
    # Delta_zeta2  = mode_k**3 * Ps_zeta2.pop () / (2.0 * math.pi**2 * cosmo.RH_planck ()**2)
    Delta_S1 = mode_k**3 * Ps_S1.pop() / (2.0 * math.pi**2 * cosmo.RH_planck() ** 2)
    # Delta_S2     = mode_k**3 * Ps_S2.pop () / (2.0 * math.pi**2 * cosmo.RH_planck ()**2)
    Delta_Pzeta1 = (
        mode_k**3 * Ps_Pzeta1.pop() / (2.0 * math.pi**2 * cosmo.RH_planck() ** 2)
    )
    # Delta_Pzeta2 = mode_k**3 * Ps_Pzeta2.pop () / (2.0 * math.pi**2 * cosmo.RH_planck ()**2)
    Delta_PS1 = (
        mode_k**3 * Ps_PS1.pop() / (2.0 * math.pi**2 * cosmo.RH_planck() ** 2)
    )
    # Delta_PS2    = mode_k**3 * Ps_PS2.pop () / (2.0 * math.pi**2 * cosmo.RH_planck ()**2)

    # print ("# Final values k= % 20.15g Ps_zeta1 = % 21.15e Ps_zeta2 = % 21.15e
    # Ps_S1 = % 21.15e Ps_S2 = % 21.15e" % (mode_k, Delta_zeta1,  Delta_zeta2,
    # Delta_S1,  Delta_S2))
    # print ("# Final values k= % 20.15g Ps_Pzeta1= % 21.15e Ps_Pzeta2= % 21.15e
    # Ps_PS1= % 21.15e Ps_PS2= % 21.15e" % (mode_k, Delta_Pzeta1, Delta_Pzeta2,
    # Delta_PS1, Delta_PS2))
    print(
        f"# Final values k = {mode_k: 20.15g} Ps_zeta1 = {Delta_zeta1: 21.15e} "
        f"Ps_Pzeta1 = {Delta_Pzeta1: 21.15e} Ps_S1 = {Delta_S1: 21.15e} "
        f"Ps_PS1 = {Delta_PS1: 21.15e}"
    )
    # print ("# Final values k= % 20.15g Ps_zeta2 = % 21.15e Ps_Pzeta2 = % 21.15e
    # Ps_S2 = % 21.15e Ps_PS2 = % 21.15e" % (mode_k, Delta_zeta2,  Delta_Pzeta2,
    # Delta_S2,  Delta_PS2))

    plt.show()
    plt.clf()

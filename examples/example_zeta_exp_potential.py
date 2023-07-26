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
        Vexp.props.alphab = +1.0e-20
        Vexp.props.sigmaphi = +8.0e-1
        Vexp.props.dphi = +5.0e-1
        Vexp.props.xb = 2.0e38
        Vexp.props.OmegaL = 1.0
        Vexp.props.Omegac = 1.0
        Vexp.props.H0 = 67.8

    k = 1.0e0
    tc = Vexp.tau_xe(1.0e15)
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
    adiab.prepare(Vexp)

    # exit ()

    print("# Preparing GW")
    gw.prepare(Vexp)

    # (t0, t1) = adiab.get_t0_t1 (Vexp)
    (t0, t1) = gw.get_t0_t1(Vexp)

    print(f"# BACKG (t0, t1) = ({Vexp.tau_min(): 21.15e}, {Vexp.tau_max(): 21.15e})")
    print(f"# ADIAB (t0, t1) = ({t0: 21.15e}, {t1: 21.15e})")

    (Delta_zeta_c, _Delta_Pzeta_c) = adiab.eval_Delta(Vexp, tc)
    (Delta_h_c, _Delta_Ph_c) = gw.eval_Delta(Vexp, tc)

    print(f"# Time of x_e = 10^15:     tau_c     = {tc: 21.15f}")
    print(f"# Power spectrum at tau_c: PS_ADIAB  = {Delta_zeta_c: 21.15e}")
    print(f"# Power spectrum at tau_c: PS_GW     = {Delta_h_c: 21.15e}")
    print(
        f"# Power spectrum at tau_c: r         = {2.0 * Delta_h_c / Delta_zeta_c: 21.15e}"
    )

    t_a = np.linspace(t0, t1, 100000)
    alpha_a = []

    Delta_zeta = []
    # Delta_Pzeta = []

    Delta_h = []
    # Delta_Ph = []

    zeta_a_a = []
    zeta_b_a = []
    Pzeta_a_a = []
    Pzeta_b_a = []

    h_a_a = []
    h_b_a = []
    Ph_a_a = []
    Ph_b_a = []

    nu_a = []
    m_a = []
    mnu_a = []
    dlnmnu_a = []

    epsilon_a = []
    gamma_a = []
    sin_thetab_a = []
    cos_thetab_a = []

    for t in t_a:
        (Delta_zeta_v, _Delta_Pzeta_v) = adiab.eval_Delta(Vexp, t)
        (Delta_h_v, _Delta_Ph_v) = gw.eval_Delta(Vexp, t)

        (zeta_a_v, zeta_b_v, Pzeta_a_v, Pzeta_b_v) = adiab.eval_QV(Vexp, t)
        (h_a_v, h_b_v, Ph_a_v, Ph_b_v) = gw.eval_QV(Vexp, t)

        (epsilon_v, gamma_v, sin_thetab_v, cos_thetab_v) = gw.eval_AA(Vexp, t)

        epsilon_a.append(epsilon_v)
        gamma_a.append(gamma_v)
        sin_thetab_a.append(sin_thetab_v)
        cos_thetab_a.append(cos_thetab_v)

        nu_adiab = adiab.eval_nu(Vexp, t, k)
        mnu_adiab = adiab.eval_mnu(Vexp, t, k)
        m_adiab = mnu_adiab / nu_adiab
        dlnmnu_adiab = math.fabs(adiab.eval_dlnmnu(Vexp, t, k))

        alpha = 0
        if t > 0.0:
            alpha = +0.5 * t**2 / math.log(10.0)
        else:
            alpha = -0.5 * t**2 / math.log(10.0)

        alpha_a.append(alpha)

        nu_a.append(nu_adiab)
        m_a.append(m_adiab)
        mnu_a.append(mnu_adiab)
        dlnmnu_a.append(dlnmnu_adiab)

        Delta_zeta.append(Delta_zeta_v)
        Delta_h.append(Delta_h_v)

        zeta_a_a.append(zeta_a_v)
        zeta_b_a.append(zeta_b_v)

        Pzeta_a_a.append(Pzeta_a_v)
        Pzeta_b_a.append(Pzeta_b_v)

        h_a_a.append(h_a_v)
        h_b_a.append(h_b_v)

        Ph_a_a.append(Ph_a_v)
        Ph_b_a.append(Ph_b_v)

    mylw = 1

    h_a_npa = np.array(h_a_a)
    h_b_npa = np.array(h_b_a)

    # Ph_a_npa = np.array(Ph_a_a)
    # Ph_b_npa = np.array(Ph_b_a)

    zeta_a_npa = np.array(zeta_a_a)
    zeta_b_npa = np.array(zeta_b_a)

    # Pzeta_a_npa = np.array(Pzeta_a_a)
    # Pzeta_b_npa = np.array(Pzeta_b_a)

    # plt.plot (alpha_a, Delta_zeta, lw=mylw, label = r'$\Delta_{\zeta}$')
    # plt.plot (alpha_a, Delta_h,    lw=mylw, label = r'$\Delta_{h}$')

    plt.plot(alpha_a, zeta_a_npa, lw=mylw, label=r"$\tilde{\zeta}^a$")
    plt.plot(alpha_a, zeta_b_npa, lw=mylw, label=r"$\tilde{\zeta}^b$")
    plt.plot(alpha_a, h_a_npa, lw=mylw, label=r"$\tilde{h}^a$")
    plt.plot(alpha_a, h_b_npa, lw=mylw, label=r"$\tilde{h}^b$")

    # plt.plot (alpha_a, epsilon_a,    lw=mylw, label = r'$\epsilon$')
    # plt.plot (alpha_a, gamma_a,      lw=mylw, label = r'$\gamma$')
    # plt.plot (alpha_a, sin_thetab_a, lw=mylw, label = r'$\sin(\theta_b)$')
    # plt.plot (alpha_a, cos_thetab_a, lw=mylw, label = r'$\cos(\theta_b)$')

    # plt.plot (alpha_a, nu_a,     lw=mylw, label = r'$\nu_h$')
    # plt.plot (alpha_a, m_a,      lw=mylw, label = r'$m_h$')
    # plt.plot (alpha_a, mnu_a,    lw=mylw, label = r'$m_h\nu_h$')
    # plt.plot (alpha_a, 1.0 / (np.array (mnu_a)),    lw=mylw, label = r'$(m_h\nu_h)^{-1}$')
    # plt.plot (alpha_a, dlnmnu_a, lw=mylw, label = r'$d\ln(m_h\nu_h)$')

    # plt.plot (alpha_a, (Ph_a_a * h_b_a), lw=mylw, label = r't_1')
    # plt.plot (alpha_a, (Ph_b_a * h_a_a), lw=mylw, label = r't_2')

    plt.grid(b=True, which="both", linestyle=":", color="0.75", linewidth=0.5)
    plt.legend(loc="upper left")

    # plt.xscale('symlog', linthreshx=1.0e-5)
    plt.yscale("symlog", linthreshy=1.0e23)
    # plt.yscale('log', linthreshy=1.0e-20)

    plt.show()
    plt.clf()


if __name__ == "__main__":
    test_zeta_exp_potential()

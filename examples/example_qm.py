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

"""Example showing QM wave function evolution."""

import sys
import math
from typing import Union

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as animation

from tqdm import tqdm

from numcosmo_py import Nc, Ncm

matplotlib.use("Agg")

# plt.rcParams['animation.ffmpeg_path'] = '/usr/bin/ffmpeg'

font = {"size": 20}
matplotlib.rc("font", **font)


def _blit_draw(_self, artists, bg_cache):
    # Handles blitted drawing, which renders only the artists given instead
    # of the entire figure.
    updated_ax = []
    for a in artists:
        # If we haven't cached the background for this axes object, do
        # so now. This might not always be reliable, but it's an attempt
        # to automate the process.
        if a.axes not in bg_cache:
            # bg_cache[a.axes] = a.figure.canvas.copy_from_bbox(a.axes.bbox)
            # change here
            bg_cache[a.axes] = a.figure.canvas.copy_from_bbox(a.axes.figure.bbox)
        a.axes.draw_artist(a)
        updated_ax.append(a.axes)

    # After rendering all the needed artists, blit each axes individually.
    for ax in set(updated_ax):
        # and here
        # ax.figure.canvas.blit(ax.bbox)
        ax.figure.canvas.blit(ax.figure.bbox)


# MONKEY PATCH!!
# pylint: disable-next=protected-access
matplotlib.animation.Animation._blit_draw = _blit_draw

#
#  Initializing the library objects, this must be called before
#  any other library function.
#
Ncm.cfg_init()


def test_qm() -> None:
    """Example showing QM wave function evolution."""

    if len(sys.argv) < 1 + 4:
        print("Usage: example_qm.py (gauss|exp) dist_param lambda H0")
        print(len(sys.argv))
        sys.exit()

    init_wf = sys.argv[1]
    dist_param = float(sys.argv[2])
    l1 = float(sys.argv[3])
    H0 = float(sys.argv[4])

    p = Nc.HIQG1D.new()
    # l1                 = 0.0
    # H0                 = -15.0e-1
    p.props.abstol = 0.0
    p.props.reltol = 1.0e-7
    p.props.nknots = 1200
    p.props.noboundary = False
    p.set_property("lambda", l1)
    offset = 5.0
    center = (10.0 + offset) * 0.0 + 0.0

    print(
        f"# lambda = {p.get_lambda(): 22.15g}, basis a = {p.get_basis_a(): 22.15g}, "
        f"acs a = {p.get_acs_a(): 22.15g} nu = {p.get_nu(): 22.15g} "
        f"mu = {p.get_mu(): 22.15g}"
    )

    psi0: Union[Nc.HIQG1DGauss, Nc.HIQG1DExp]
    if init_wf == "gauss":
        psi0 = Nc.HIQG1DGauss.new(center + dist_param, 1.0, 1.0, H0)
    elif init_wf == "exp":
        psi0 = Nc.HIQG1DExp.new(p.get_acs_a(), dist_param, H0)
    else:
        psi0 = Nc.HIQG1DExp.new(3.0, dist_param, H0)

    # psi0 = Nc.HIQG1DExp.new (3.0, 2.0, H0)

    # print (psi0.eval (1.0))
    npp = 1000
    sim = False
    tstep = 2.0 / 800.0
    tf = 2.0
    xf = center + 100.0
    xfp = center + 40.0
    xi = (center - 7.0) * 0.0 + 0.0
    dSdiv = 10.0

    if init_wf == "gauss":
        assert isinstance(psi0, Nc.HIQG1DGauss)
        p.set_init_cond_gauss(psi0, xi, xf)
    else:
        assert isinstance(psi0, Nc.HIQG1DExp)
        p.set_init_cond_exp(psi0, xi, xf)

    p.prepare()
    p.evol(0.0)
    n1 = p.nBohm()

    # q = p.peek_knots().dup_array()
    xa = np.linspace(xi, xfp, npp)
    # psi_base = [p.eval_psi(x) for x in xa]
    # max_Re_psi = max(np.abs(psi_base[:][0]))
    # max_Im_psi = max(np.abs(psi_base[:][1]))
    # yb = max(max_Re_psi, max_Im_psi)

    fig, [ax, ax2] = plt.subplots(1, 2, figsize=(16, 12))

    ax.set_xlim(xi, xfp)
    # ax.set_ylim (-1.5, 1.5)
    # ax.grid ()

    # ax2.set_xlim (-0.1, 2.0)
    # ax2.set_ylim (-8.0, 8.0)

    ax.set_xlabel(r"$V$")
    ax2.set_xlabel(r"$V$")
    ax2.set_ylabel(r"$P_V$")
    # ax2.grid ()

    # N = 4
    ttl = ax.text(0.1, 1.005, "", transform=ax.transAxes)
    lines = []
    lines.append(ax.plot([], [], label=r"$\sqrt{\psi^*\psi}$", animated=True)[0])
    lines.append(ax.plot([], [], label=r"$\mathrm{Re}(\psi)$", animated=True)[0])
    lines.append(ax.plot([], [], label=r"$\mathrm{Im}(\psi)$", animated=True)[0])
    lines.append(ax.plot([], [], animated=True)[0])
    # lines.append (ax.plot ([], [], label=r'$\partial_aS$', animated=True)[0])
    ax.legend(loc="best")

    fig.tight_layout()
    lines.append(ttl)

    lines.append(ax2.plot([], [], "bo", label=r"Bohm", animated=True)[0])
    lines.append(ax2.plot([], [], "ro", label=r"Semi-classical", animated=True)[0])
    lines.append(ax2.plot([], [], "go", label=r"Expected values", animated=True)[0])
    ax2.legend(loc="best")

    # ta = [0.0]
    traj = []
    trajP = []
    x_t = []
    y_t = []
    Q_t = []
    P_t = []

    Pini = p.Bohm_p(0)
    Vini = p.Bohm(0)
    Hini = ((Pini * Vini) ** 2 + l1) / Vini**2
    t0 = -0.5 * Pini * Vini / Hini

    def a_sc(t):
        return math.sqrt(4.0 * Hini * (t - t0) ** 2 + l1 / Hini)

    def p_sc(t):
        return 2.0 * Hini * (t - t0) / a_sc(t)

    Q = []
    P = []

    Q.append([p.Bohm(i) for i in range(n1)])
    P.append([p.Bohm_p(i) for i in range(n1)])

    Q.append([a_sc(0.0)])
    P.append([p_sc(0.0)])

    Q.append([p.int_xrho_0_inf()])
    P.append([p.expect_p()])

    for i in range(n1):
        lines.append(ax2.plot([], [])[0])
        qi = [p.Bohm(i)]
        pi = [p.Bohm_p(i)]
        traj.append(qi)
        trajP.append(pi)

    Q_t.append(Q)
    P_t.append(P)

    lines.append(ax2.plot([], [])[0])
    traj.append([a_sc(0.0)])
    trajP.append([p_sc(0.0)])

    lines.append(ax2.plot([], [])[0])
    traj.append([p.int_xrho_0_inf()])
    trajP.append([p.expect_p()])

    def init():
        lines[0].set_data([], [])
        lines[1].set_data([], [])
        lines[2].set_data([], [])
        lines[3].set_data([], [])
        lines[4].set_text(None)
        lines[5].set_data([], [])
        lines[6].set_data([], [])
        lines[7].set_data([], [])

        for i in range(n1):
            lines[i + 8].set_data([], [])

        for p in lines:
            p.set_visible(False)

        return lines

    psi_max = 0.0
    psi_min = 0.0
    Q_min = 1.0e20
    Q_max = 0.0
    P_min = 0.0
    P_max = 0.0
    kk_min = 1000000000

    if not sim:
        for i in tqdm(np.arange(0, math.ceil(tf / tstep))):
            tf = tstep * i
            p.evol(tf)

            # q = p.peek_knots().dup_array()
            xa = np.linspace(xi, xfp, npp)
            psi = np.array([p.eval_psi(x) for x in xa])
            rho = np.array([np.sum(psi_i**2) for psi_i in psi])
            dS = np.array([0.0] * len(xa))

            sqrt_rho = np.sqrt(rho)
            y = []
            y.append(sqrt_rho)
            y.append(psi[:, 0])
            y.append(psi[:, 1])
            y.append(dS / dSdiv)
            y.append(p.int_rho_0_inf())

            kk = 0
            for rho_kk in reversed(sqrt_rho):
                if rho_kk > 0.05:
                    break
                else:
                    kk = kk + 1

            kk_min = min(kk_min, kk)
            # pylint: disable=nested-min-max
            psi_max = max(max(psi[:, 0]), psi_max)
            psi_max = max(max(psi[:, 1]), psi_max)
            psi_max = max(np.sqrt(rho), psi_max)
            psi_min = min(min(psi[:, 0]), psi_min)
            psi_min = min(min(psi[:, 1]), psi_min)

            x_t.append(xa)
            y_t.append(y)

            Q = []
            P = []

            QBohm = [p.Bohm(i) for i in range(n1)]
            PBohm = [p.Bohm_p(i) for i in range(n1)]
            Q.append(QBohm)
            P.append(PBohm)

            QSC = [a_sc(tf)]
            PSC = [p_sc(tf)]
            Q.append(QSC)
            P.append(PSC)

            QMEAN = [p.int_xrho_0_inf()]
            PMEAN = [p.expect_p()]
            Q.append(QMEAN)
            P.append(PMEAN)

            Q_min = min(Q_min, min(QBohm))
            Q_min = min(Q_min, min(QSC))
            Q_min = min(Q_min, min(QMEAN))

            Q_max = max(Q_max, max(QBohm))
            Q_max = max(Q_max, max(QSC))
            Q_max = max(Q_max, max(QMEAN))

            P_min = min(P_min, min(PBohm))
            P_min = min(P_min, min(PSC))
            P_min = min(P_min, min(PMEAN))

            P_max = max(P_max, max(PBohm))
            P_max = max(P_max, max(PSC))
            P_max = max(P_max, max(PMEAN))

            Q_t.append(Q)
            P_t.append(P)

            for i in range(n1):
                traj[i].append(p.Bohm(i))
                trajP[i].append(p.Bohm_p(i))

            traj[n1].append(a_sc(tf))
            trajP[n1].append(p_sc(tf))

            traj[n1 + 1].append(p.int_xrho_0_inf())
            trajP[n1 + 1].append(p.expect_p())
            # print (tf)

    # ax.set_xlim (xi, xfp)
    ax.set_xlim(xi, xa[npp - kk_min])
    ax.set_ylim(psi_min * 1.1, psi_max * 1.1)
    ax.grid()

    ax2.set_xlim(Q_min * 0.9, Q_max * 1.1)
    ax2.set_ylim(P_min * 1.1, P_max * 1.1)

    ax.set_xlabel(r"$V$")
    ax2.set_xlabel(r"$V$")
    ax2.set_ylabel(r"$P_V$")
    ax2.grid()

    def animate(i):
        tf = tstep * i

        if i == 1:
            for pl in lines:
                pl.set_visible(True)

        if sim:
            p.evol(tf)

            # q = p.peek_knots().dup_array()
            xa = np.linspace(xi, xfp, npp)
            # print ("% 22.15g % 22.15g % 22.15g % 22.15g" % (q[0], q[1], q[2], q[3]))

            psi = np.array([p.eval_psi(x) for x in xa])
            rho = np.array([np.sum(psi_i**2) for psi_i in psi])
            dS = np.array([0.0] * len(xa))
            # dS  = np.array ([p.eval_dS (x) for x in xa])

            # print (dS / xa)

            lines[0].set_data(xa, np.sqrt(rho))
            lines[1].set_data(xa, psi[:, 0])
            lines[2].set_data(xa, psi[:, 1])
            lines[3].set_data(xa, dS / dSdiv)

            lines[5].set_data(
                [p.Bohm(i) for i in range(n1)], [p.Bohm_p(i) for i in range(n1)]
            )
            lines[6].set_data([a_sc(tf)], [p_sc(tf)])
            lines[7].set_data([p.int_xrho_0_inf()], [p.expect_p()])

            for i in range(n1):
                traj[i].append(p.Bohm(i))
                trajP[i].append(p.Bohm_p(i))
                lines[i + 8].set_data(traj[i], trajP[i])

            traj[n1].append(a_sc(tf))
            trajP[n1].append(p_sc(tf))
            lines[n1 + 8].set_data(traj[n1], trajP[n1])

            traj[n1 + 1].append(p.int_xrho_0_inf())
            trajP[n1 + 1].append(p.expect_p())
            lines[n1 + 8 + 1].set_data(traj[n1 + 1], trajP[n1 + 1])

            ttl.set_text(f"t = {tf: .15f}, norma = {p.int_rho_0_inf(): .15f}")

        else:
            lines[0].set_data(x_t[i], y_t[i][0])
            lines[1].set_data(x_t[i], y_t[i][1])
            lines[2].set_data(x_t[i], y_t[i][2])
            lines[3].set_data(x_t[i], y_t[i][3])

            lines[5].set_data(Q_t[i][0], P_t[i][0])
            lines[6].set_data(Q_t[i][1], P_t[i][1])
            lines[7].set_data(Q_t[i][2], P_t[i][2])

            ti = i + 1
            for j in range(n1):
                lines[j + 8].set_data(traj[j][0:ti], trajP[j][0:ti])

            lines[n1 + 8].set_data(traj[n1][0:ti], trajP[n1][0:ti])
            lines[n1 + 8 + 1].set_data(traj[n1 + 1][0:ti], trajP[n1 + 1][0:ti])

            ttl.set_text(f"t = {tf: .15f}, norma = {y_t[i][4]: .15f}")

        return lines

    anim = animation.FuncAnimation(
        fig,
        animate,
        np.arange(0, int(tf / tstep)),
        init_func=init,
        interval=1,
        blit=True,
        repeat=True,
    )

    # mywriter = animation.FFMpegWriter (fps = 24, codec='libx264')
    # mywriter = animation.FFMpegWriter (bitrate=500)

    Writer = animation.writers["ffmpeg"]
    writer = Writer(fps=15, metadata=dict(artist="Me"), bitrate=1800)

    anim.save(
        f"evol_{init_wf}_dp_{dist_param}_lambda_{l1}_H0_{H0}.mp4",
        writer=writer,
    )

    # plt.show ()


if __name__ == "__main__":
    test_qm()

#
# planck_commander.py
#
# Wed July 23 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# planck_commander.py
# Copyright (C) 2026 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

"""Build the native NcDataPlanckCommander from a Planck low-l commander clik file.

Parses the gibbs_gauss ``sigma.fits`` payload (per-l Cl->x cubic spline, mean,
covariance and mean-sigma) and maps it onto the construct properties of the
native #NcDataPlanckCommander. Validated to reproduce the clik likelihood.
"""

import os

import numpy as np
from astropy.io import fits

from numcosmo_py import Ncm, Nc

COMMANDER_RELPATH = os.path.join(
    "baseline", "plc_3.0", "low_l", "commander", "commander_dx12_v3_2_29.clik"
)


def build_commander(
    clik_path: str,
    pb: Nc.HIPertBoltzmann | None = None,
    clik_pi_compat: bool = False,
) -> Nc.DataPlanckCommander:
    """Build a native #NcDataPlanckCommander from a commander clik file + Cls source.

    With ``clik_pi_compat=True`` the Dl conversion uses clik's single-precision pi,
    reproducing the clik gibbs likelihood bit-for-bit (otherwise full double pi is
    used, which is more accurate but differs from clik by ~2e-6 in m2lnL).
    """
    lkl = os.path.join(clik_path, "clik", "lkl_0")
    lmin, lmax, delta_l = 2, 29, 1000
    nl = lmax - lmin + 1

    with fits.open(os.path.join(lkl, "_external", "sigma.fits")) as h:
        cl2x = np.array(h[0].data)  # (3, 249, 1000) = [comp, l=2..250, node]
        mu_all = np.array(h[1].data)  # (249,)
        cov_all = np.array(h[2].data)  # (249,249) covariance
        mu_sigma_all = np.array(h[3].data)  # (249,)

    nbin = cl2x.shape[2]
    sl = slice(0, nl)  # l=2..29 -> idx 0..27
    xa = cl2x[0, sl, :].ravel()  # [l_idx*nbin + node]
    ya = cl2x[1, sl, :].ravel()
    y2a = cl2x[2, sl, :].ravel()
    mu = mu_all[sl]
    mu_sigma = mu_sigma_all[sl]
    cov = cov_all[sl, sl]

    def vec(a):
        return Ncm.Vector.new_array([float(x) for x in np.asarray(a).ravel().tolist()])

    cmd = Nc.DataPlanckCommander.new(
        lmin,
        lmax,
        nbin,
        delta_l,
        vec(xa),
        vec(ya),
        vec(y2a),
        vec(mu),
        Ncm.Matrix.new_array(cov.ravel().tolist(), nl),
        vec(mu_sigma),
    )
    cmd.set_property("clik-pi-compat", clik_pi_compat)

    if pb is not None:
        cmd.set_hipert_boltzmann(pb)
    cmd.set_init(True)

    return cmd

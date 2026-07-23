#
# planck_simall.py
#
# Wed July 23 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# planck_simall.py
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

"""Build the native NcDataPlanckSimall from a Planck low-l simall clik file.

Parses the per-multipole tabulated log-probability payload (probEE/probBB/probTE
with their Dl steps) and maps it onto the construct properties of the native
#NcDataPlanckSimall. Validated to reproduce the clik likelihood.
"""

import os

import numpy as np
from astropy.io import fits

from numcosmo_py import Ncm, Nc

SIMALL_EE_RELPATH = os.path.join(
    "baseline",
    "plc_3.0",
    "low_l",
    "simall",
    "simall_100x143_offlike5_EE_Aplanck_B.clik",
)
SIMALL_BB_RELPATH = os.path.join(
    "baseline",
    "plc_3.0",
    "low_l",
    "simall",
    "simall_100x143_offlike5_BB_Aplanck_B.clik",
)
SIMALL_EEBB_RELPATH = os.path.join(
    "baseline",
    "plc_3.0",
    "low_l",
    "simall",
    "simall_100x143_offlike5_EEBB_Aplanck_B.clik",
)


def _read_mdb(node):
    """Parse a cldf ``_mdb`` file into a {key: (type, value)} mapping."""
    out = {}
    with open(os.path.join(node, "_mdb"), encoding="ascii") as f:
        for line in f:
            parts = line.split()
            if len(parts) >= 3:
                out[parts[0]] = (parts[1], parts[2])
    return out


def build_simall(
    clik_path: str, pb: Nc.HIPertBoltzmann | None = None
) -> Nc.DataPlanckSimall:
    """Build a native #NcDataPlanckSimall from a simall clik file + Cls source."""
    lkl = os.path.join(clik_path, "clik", "lkl_0")
    mdb = _read_mdb(lkl)
    lmin, lmax = int(mdb["lmin"][1]), int(mdb["lmax"][1])
    has_cl = fits.getdata(os.path.join(lkl, "has_cl")).astype(int).ravel()

    def vec(a):
        return Ncm.Vector.new_array([float(x) for x in np.asarray(a).ravel().tolist()])

    # has_cl order is [TT, EE, BB, TE, TB, EB]; simall uses EE(1), BB(2), TE(3)
    def spec(tag, idx):
        if not has_cl[idx]:
            return 0.0, None
        step = float(mdb[f"step{tag}"][1])
        prob = vec(fits.getdata(os.path.join(lkl, f"prob{tag}")))
        return step, prob

    step_ee, prob_ee = spec("EE", 1)
    step_bb, prob_bb = spec("BB", 2)
    step_te, prob_te = spec("TE", 3)

    simall = Nc.DataPlanckSimall.new(
        lmin, lmax, step_ee, prob_ee, step_bb, prob_bb, step_te, prob_te
    )

    if pb is not None:
        simall.set_hipert_boltzmann(pb)
    simall.set_init(True)

    return simall

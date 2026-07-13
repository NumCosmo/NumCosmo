#
# planck_smica.py
#
# Fri June 27 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# planck_smica.py
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

"""Build the native NcDataPlanckSmica from a Planck ``plik`` (SMICA) clik file.

Parses the cldf SMICA payload (bandpower data, criterion matrix, binning, CMB
mixing and the foreground/calibration component templates) and maps it onto the
construct properties of the native #NcDataPlanckSmica (an #NcmDataGaussCov).
The assembly is validated to reproduce the clik likelihood to machine precision.
"""

import os

import numpy as np
from astropy.io import fits

from numcosmo_py import Ncm, Nc, GLib

PLIK_TT_RELPATH = os.path.join(
    "baseline", "plc_3.0", "hi_l", "plik", "plik_rd12_HM_v22_TT.clik"
)
FREQS = [100.0, 143.0, 217.0]
# gcib muK -> MJ/sr conversions (100,143,217); gibXsz conv (100 zeroed by
# no_szxcib_100); tSZ colour corrections. All from the clik component defaults.
# NOTE: clik stores parameter values as 6-significant-figure strings ("%g"), so
# these conversion constants are truncated to 6 sig figs before use. We replicate
# that truncation to match clik to machine precision (the full-precision literals
# are 4096.68168783e-6, 2690.05218701e-6, 2067.43988919e-6).
GCIB_CONV = [0.00409668, 0.00269005, 0.00206744]
GIBXSZ_CONV = [0.0, 0.094, 1.0]
SZ_COLOR = [0.981, 0.975, 1.0]


def _arr(node, key):
    return fits.getdata(os.path.join(node, key)).ravel()


def build_smica_tt(clik_path: str, pb: Nc.HIPertBoltzmann | None = None) -> Nc.DataPlanckSmica:
    """Build a native #NcDataPlanckSmica (TT) from a plik clik file + Cls source."""
    lkl = os.path.join(clik_path, "clik", "lkl_0")
    m, nq = 3, 215
    lmin, lmax = 30, 2508

    rq_hat = _arr(lkl, "Rq_hat")                       # nq*m*m
    crit = _arr(lkl, "criterion_gauss_mat")            # np*np
    mask = _arr(lkl, "criterion_gauss_mask").astype(int).reshape(nq, m, m)
    ordering = _arr(lkl, "criterion_gauss_ordering").astype(int).reshape(-1, 2)
    a_cmb = _arr(lkl, "A_cmb")
    bin_lmin = _arr(lkl, "bin_lmin")
    bin_lmax = _arr(lkl, "bin_lmax")
    bin_ws = _arr(lkl, "bin_ws")

    # masked-vector index list: order = ordering pairs x bins (matches crit_cor)
    quad = []
    for iv, jv in ordering:
        for iq in range(nq):
            if mask[iq, iv, jv]:
                quad.append(iq * m * m + iv * m + jv)
    quad = np.array(quad, dtype=int)
    npt = quad.size

    crit = crit.reshape(npt, npt)
    data_mean = rq_hat[quad]
    cov = np.linalg.inv(crit)

    # templates (sz/ksz pre-normalized at ell=3000, matching the clik init)
    sz_t = _arr(lkl, os.path.join("component_3", "template")).copy()   # ell 2..
    sz_t /= sz_t[3000 - 2]
    ksz_t = _arr(lkl, os.path.join("component_5", "template")).copy()  # ell 0..
    ksz_t /= ksz_t[3000]
    gcib_t = _arr(lkl, os.path.join("component_1", "template"))
    gibxsz_t = _arr(lkl, os.path.join("component_2", "template"))
    dust_t = _arr(lkl, os.path.join("component_6", "template"))
    leak_t = _arr(lkl, os.path.join("component_7", "template"))
    sbpx_t = _arr(lkl, os.path.join("component_8", "template"))

    def vec(a):
        return Ncm.Vector.new_array([float(x) for x in np.asarray(a).ravel().tolist()])

    def uvar(a):
        return GLib.Variant("au", [int(x) for x in np.asarray(a).ravel().tolist()])

    # All structural config is construction-only: gather everything, then build
    # in a single call so NcDataPlanckSmica can validate lengths/bounds up front.
    smica = Nc.DataPlanckSmica(
        lmin=lmin,
        lmax=lmax,
        m_channels=m,
        nbins=nq,
        freqs=vec(FREQS),
        a_cmb=vec(a_cmb),
        sz_color=vec(SZ_COLOR),
        gcib_conv=vec(GCIB_CONV),
        gibxsz_conv=vec(GIBXSZ_CONV),
        bin_lmin=uvar(bin_lmin),
        bin_lmax=uvar(bin_lmax),
        bin_weight=vec(bin_ws),
        quad_idx=uvar(quad),
        tmpl_gcib=vec(gcib_t),
        tmpl_sz=vec(sz_t),
        tmpl_ksz=vec(ksz_t),
        tmpl_gibxsz=vec(gibxsz_t),
        tmpl_dust=vec(dust_t),
        tmpl_leak=vec(leak_t),
        tmpl_sbpx=vec(sbpx_t),
    )

    smica.set_size(npt)
    Ncm.DataGaussCov.replace_mean(smica, Ncm.Vector.new_array(data_mean.tolist()))
    cov_m = Ncm.Matrix.new(npt, npt)
    for i in range(npt):
        for j in range(npt):
            cov_m.set(i, j, cov[i, j])
    Ncm.DataGaussCov.set_cov(smica, cov_m)

    if pb is not None:
        smica.set_hipert_boltzmann(pb)
    smica.set_init(True)

    return smica

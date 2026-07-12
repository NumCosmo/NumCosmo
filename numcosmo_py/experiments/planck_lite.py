#
# planck_lite.py
#
# Fri June 27 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# planck_lite.py
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

"""Build the native NcDataPlanckPlikLite from a Planck ``plik_lite`` clik file.

This is the ingestion half of the native Planck port: the awkward cldf/Fortran
data layout is parsed here (in Python) and mapped onto the construct properties
of the native #NcDataPlanckPlikLite (an #NcmDataGaussCov), which the C runtime
then evaluates and serializes. No clik/cldf code runs in C.

The mapping is validated to reproduce the clik ``check_value`` to machine
precision (see tests/python/nc/data/test_planck_plik_lite.py).
"""

import os

import numpy as np

from numcosmo_py import Ncm, Nc

# plik v22 fixed layout.
PLMIN = 30
NBIN_TT = 215
NBIN_TOTAL = 613  # 215 TT + 199 TE + 199 EE
PLIK_LITE_TT_RELPATH = os.path.join(
    "baseline", "plc_3.0", "hi_l", "plik_lite", "plik_lite_v22_TT.clik"
)


def find_baseline_file(relpath: str) -> str | None:
    """Return the absolute path of a baseline data file, or None if absent.

    Uses the same base directory as the clik data downloader
    (``ncm_cfg_get_fullpath_base``), falling back to ``~/.numcosmo``.
    """
    try:
        base = Ncm.cfg_get_fullpath_base()
    except Exception:  # pylint: disable=broad-except
        base = os.path.expanduser("~/.numcosmo")
    candidate = os.path.join(base, relpath)
    if os.path.exists(candidate):
        return candidate
    fallback = os.path.join(os.path.expanduser("~/.numcosmo"), relpath)
    return fallback if os.path.exists(fallback) else None


def _read_fortran_unformatted_matrix(path: str, n: int) -> np.ndarray:
    """Read an ``n x n`` f64 matrix from a Fortran sequential-unformatted file.

    The record is wrapped by 4-byte length markers; the payload is column-major.
    Only the upper triangle is trusted (the clik writer leaves the lower part
    stale), so it is symmetrized from the upper triangle.
    """
    raw = np.fromfile(path, dtype=np.uint8)
    body = raw[4 : 4 + n * n * 8].view(np.float64)
    mat = body.reshape(n, n, order="F").copy()
    iu = np.triu_indices(n)
    mat[(iu[1], iu[0])] = mat[iu]
    return mat


def read_plik_lite_tt(clik_path: str) -> dict:
    """Parse the ``plik_lite`` TT clik ``_external`` payload (validated recipe).

    Returns a dict with the TT-selected data vector, covariance and per-bin
    binning operator (absolute multipole ranges + flattened averaging weights).
    """
    ext = os.path.join(clik_path, "clik", "lkl_0", "_external")

    x_data = np.loadtxt(os.path.join(ext, "cl_cmb_plik_v22.dat"))[:, 1]
    cov = _read_fortran_unformatted_matrix(
        os.path.join(ext, "c_matrix_plik_v22.dat"), NBIN_TOTAL
    )
    blmin = np.loadtxt(os.path.join(ext, "blmin.dat")).astype(int)[:NBIN_TT]
    blmax = np.loadtxt(os.path.join(ext, "blmax.dat")).astype(int)[:NBIN_TT]
    bin_w = np.loadtxt(os.path.join(ext, "bweight.dat"))

    lmin_abs = (blmin + PLMIN).tolist()
    lmax_abs = (blmax + PLMIN).tolist()
    weights: list[float] = []
    for lo, hi in zip(blmin, blmax):
        weights.extend(bin_w[lo : hi + 1].tolist())

    return {
        "x_data": x_data[:NBIN_TT],
        "cov": cov[:NBIN_TT, :NBIN_TT],
        "bin_lmin": lmin_abs,
        "bin_lmax": lmax_abs,
        "bin_weight": weights,
        "spectrum_id": [0] * NBIN_TT,  # 0 = TT
    }


def build_plik_lite_tt(
    clik_path: str, pb: Nc.HIPertBoltzmann | None = None, lmax: int = 2508
) -> Nc.DataPlanckPlikLite:
    """Build a native #NcDataPlanckPlikLite (TT) from a clik file and a Cls source.

    @clik_path: path to a ``plik_lite_v22_TT.clik`` directory.
    @pb: the theory Cls source (use a #NcHIPertBoltzmannCBE, the CLASS backend).
      May be %None to build a data-only object (e.g. for serialization).
    """
    d = read_plik_lite_tt(clik_path)
    np_bins = len(d["bin_lmin"])

    plik = Nc.DataPlanckPlikLite()
    plik.set_size(np_bins)
    Ncm.DataGaussCov.replace_mean(plik, Ncm.Vector.new_array(d["x_data"].tolist()))

    cov = Ncm.Matrix.new(np_bins, np_bins)
    cov_np = d["cov"]
    for i in range(np_bins):
        for j in range(np_bins):
            cov.set(i, j, cov_np[i, j])
    Ncm.DataGaussCov.set_cov(plik, cov)

    plik.set_property("bin-lmin", Ncm.Vector.new_array([float(x) for x in d["bin_lmin"]]))
    plik.set_property("bin-lmax", Ncm.Vector.new_array([float(x) for x in d["bin_lmax"]]))
    plik.set_property("bin-weight", Ncm.Vector.new_array(d["bin_weight"]))
    plik.set_property("spectrum-id", Ncm.Vector.new_array([float(x) for x in d["spectrum_id"]]))
    plik.set_property("lmax", lmax)
    plik.set_property("calib-name", "A_planck")
    if pb is not None:
        plik.set_hipert_boltzmann(pb)
    plik.set_init(True)

    return plik

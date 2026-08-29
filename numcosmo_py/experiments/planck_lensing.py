#
# planck_lensing.py
#
# Thu July 23 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# planck_lensing.py
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

"""Build the native NcDataPlanckLensing from a Planck 2018 lensing clik file.

Parses the ``clik_lensing`` cldf payload (itype=4, ``dts_lensing``) and maps it
onto the construct properties of the native #NcDataPlanckLensing. Handles both
the full (``hascl`` selects CMB spectra for the estimator renormalization) and
the CMB-marginalized (``hascl`` all zero, pure C_l^phiphi) baseline files.
Validated to reproduce the clik lensing likelihood.
"""

import os

import numpy as np

from numcosmo_py import Ncm, Nc, GLib

LENSING_FULL_RELPATH = os.path.join(
    "baseline",
    "plc_3.0",
    "lensing",
    "smicadx12_Dec5_ftl_mv2_ndclpp_p_teb_consext8.clik_lensing",
)
LENSING_MARGED_RELPATH = os.path.join(
    "baseline",
    "plc_3.0",
    "lensing",
    "smicadx12_Dec5_ftl_mv2_ndclpp_p_teb_consext8_CMBmarged.clik_lensing",
)


def _read_cldf_int(node_path: str, key: str) -> int:
    """Read a scalar int from a cldf node's ``_mdb``."""
    with open(os.path.join(node_path, "_mdb"), "r", encoding="latin-1") as fh:
        for line in fh:
            parts = line.split()
            if len(parts) >= 3 and parts[0] == key:
                return int(parts[2])
    raise KeyError(f"{key} not found in {node_path}/_mdb")


def _read_cldf_array(node_path: str, key: str) -> np.ndarray:
    """Read a cldf array node (FITS payload) as a flat float64 array."""
    from astropy.io import fits  # pylint: disable=import-outside-toplevel

    with fits.open(os.path.join(node_path, key)) as h:
        hdu0 = next(iter(h))
        assert hdu0 is not None, f"{node_path}/{key} has no HDU0"
        assert isinstance(hdu0, fits.PrimaryHDU)
        return np.array(hdu0.data, dtype=np.float64).ravel()


def build_lensing(
    clik_path: str,
    pb: Nc.HIPertBoltzmann | None = None,
) -> Nc.DataPlanckLensing:
    """Build a native #NcDataPlanckLensing from a lensing clik file + Cls source.

    The observed band-powers (``pp_hat``) and inverse covariance (``siginv``) are
    installed through the #NcmDataGauss mean/inv-cov; the theory C_l source is set
    separately. Works for both the full and CMB-marginalized baseline files.
    """
    node = os.path.join(clik_path, "clik_lensing")

    lmax = _read_cldf_int(node, "lmax")
    nbins = _read_cldf_int(node, "nbins")
    has_calib = bool(_read_cldf_int(node, "has_calib"))
    nl = lmax + 1

    hascl = _read_cldf_array(node, "hascl").astype(int)
    ncmb = int(hascl.sum())
    nlt = (1 + ncmb) * nl if ncmb else nl

    pp_hat = _read_cldf_array(node, "pp_hat")
    siginv = _read_cldf_array(node, "siginv").reshape(nbins, nbins)
    bins = _read_cldf_array(node, "bins").reshape(nbins, nl)
    cor0 = _read_cldf_array(node, "cor0")

    # Read the first-order correction whenever the file ships one (clik bugfix
    # 44e638b): the CMB-marginalized file also carries a phi (N0) renormalization
    # block, dropped by the old plc_3.0 loader.
    cors_mat = None
    if os.path.exists(os.path.join(node, "cors")):
        cors = _read_cldf_array(node, "cors").reshape(nbins, nlt)
        cors_mat = Ncm.Matrix.new_array(cors.ravel().tolist(), nlt)

    assert cors_mat is not None
    lens = Nc.DataPlanckLensing(
        lmax=lmax,
        nbins=nbins,
        has_calib=has_calib,
        hascl=GLib.Variant("au", [int(x) for x in hascl.tolist()]),
        bins=Ncm.Matrix.new_array(bins.ravel().tolist(), nl),
        cor0=Ncm.Vector.new_array(cor0.tolist()),
        cors=cors_mat,
    )

    lens.set_property("mean", Ncm.Vector.new_array(pp_hat.tolist()))
    lens.set_property("inv-cov", Ncm.Matrix.new_array(siginv.ravel().tolist(), nbins))

    if pb is not None:
        lens.set_hipert_boltzmann(pb)

    lens.set_init(True)

    return lens

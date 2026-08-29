#!/usr/bin/env python
#
# fixtures_planck.py
#
# Fri August 29 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# fixtures_planck.py
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

"""Synthetic Planck cldf trees for the native likelihood tests.

The Planck ``plc_3.0`` clik data is large and not redistributable, so the
data-backed tests are skipped wherever it is absent (CI included). The writers
here emit *structurally faithful but tiny* cldf trees -- the same node layout,
``_mdb`` syntax and FITS/text payload encodings the real files use -- with
made-up numbers. Feeding them to the production converters
(``numcosmo_py.experiments.planck_*``) exercises the whole ingestion path and
the native C likelihoods (prepare/mean_func/resample/serialization) anywhere,
and lets the tests check the assembled model against a closed-form expectation
that the real (opaque) data cannot provide.

The numbers are chosen so every likelihood is well posed: positive-definite
(inverse) covariances, monotonic Cl->x transforms, binning weights summing to
one and probability tables wide enough for the theory to land inside them.
"""

import os
import struct

import numpy as np
from astropy.io import fits

# ---------------------------------------------------------------------------
# cldf primitives
# ---------------------------------------------------------------------------


def write_mdb(node: str, entries: dict) -> None:
    """Write a cldf ``_mdb`` metadata file (``key type value`` per line)."""
    os.makedirs(node, exist_ok=True)
    with open(os.path.join(node, "_mdb"), "w", encoding="ascii") as fh:
        for key, value in entries.items():
            if isinstance(value, bool):
                kind, text = "int", str(int(value))
            elif isinstance(value, int):
                kind, text = "int", str(value)
            elif isinstance(value, float):
                kind, text = "float", repr(value)
            else:
                kind, text = "str", str(value)
            fh.write(f"{key} {kind} {text}\n")


def write_array(node: str, name: str, array) -> None:
    """Write a cldf array node: a FITS file holding the payload in its HDU0."""
    os.makedirs(node, exist_ok=True)
    fits.PrimaryHDU(np.ascontiguousarray(array)).writeto(
        os.path.join(node, name), overwrite=True
    )


def write_fortran_matrix(path: str, mat: np.ndarray) -> None:
    """Write a square matrix as a Fortran sequential-unformatted record."""
    payload = np.asarray(mat, dtype=np.float64).tobytes(order="F")
    with open(path, "wb") as fh:
        fh.write(struct.pack("<i", len(payload)))
        fh.write(payload)
        fh.write(struct.pack("<i", len(payload)))


# ---------------------------------------------------------------------------
# simall (low-l EE/BB/TE tabulated log-probabilities)
# ---------------------------------------------------------------------------

SIMALL_LMIN = 2
SIMALL_LMAX = 29
SIMALL_STEP = 1.0e-3
SIMALL_NSTEPS = 4000


def make_simall_cldf(
    root,
    spectra=("EE",),
    nsteps: int = SIMALL_NSTEPS,
    step: float = SIMALL_STEP,
    seed: int = 12,
) -> str:
    """Write a synthetic simall clik tree; returns its path.

    Each active spectrum gets a ``nell x nsteps`` log-probability table over
    $D_\\ell \\in [0, nsteps\\,\\mathrm{step})$. The default grid comfortably
    contains the low-l EE/BB theory, so the likelihood stays in range.
    """
    path = os.path.join(str(root), "simall_synthetic.clik")
    lkl = os.path.join(path, "clik", "lkl_0")
    nell = SIMALL_LMAX - SIMALL_LMIN + 1
    rng = np.random.default_rng(seed)

    has_cl = np.zeros(6, dtype=np.int64)
    mdb = {
        "lkl_type": "simall",
        "lmin": SIMALL_LMIN,
        "lmax": SIMALL_LMAX,
        "nell": nell,
        "unit": 1,
        "free_calib": "A_planck",
    }
    for tag, idx in (("EE", 1), ("BB", 2), ("TE", 3)):
        if tag not in spectra:
            continue
        has_cl[idx] = 1
        mdb[f"step{tag}"] = step
        mdb[f"nsteps{tag}"] = nsteps
        # Smooth, strictly negative log-probabilities: no zeros or NaNs to hide
        # a wrong table index behind.
        table = -1.0 - rng.random((nell, nsteps))
        write_array(lkl, f"prob{tag}", table)

    write_mdb(lkl, mdb)
    write_array(lkl, "has_cl", has_cl)
    write_mdb(os.path.join(path, "clik"), {"lmax": SIMALL_LMAX, "n_lkl_object": 1})
    write_mdb(path, {"clik_version": "synthetic"})

    return path


# ---------------------------------------------------------------------------
# commander (low-l TT gibbs_gauss)
# ---------------------------------------------------------------------------

COMMANDER_NL = 28  # the converter fixes l = 2..29
COMMANDER_NBIN = 60
COMMANDER_DL_MIN = 1.0
COMMANDER_DL_MAX = 5.0e3
COMMANDER_X_MIN = -4.0
COMMANDER_X_MAX = 4.0


def commander_tables():
    """Return the (xa, ya, y2a, mu, cov, mu_sigma) of the synthetic commander.

    The Cl->x transform is a straight line per multipole (zero second
    derivatives make the natural cubic spline exactly linear), so the
    likelihood has a closed form the test can check against.
    """
    xa = np.tile(
        np.linspace(COMMANDER_DL_MIN, COMMANDER_DL_MAX, COMMANDER_NBIN),
        (COMMANDER_NL, 1),
    )
    ya = np.tile(
        np.linspace(COMMANDER_X_MIN, COMMANDER_X_MAX, COMMANDER_NBIN),
        (COMMANDER_NL, 1),
    )
    y2a = np.zeros_like(xa)
    mu = np.linspace(-0.2, 0.2, COMMANDER_NL)
    # Diagonally dominant => positive definite, and not proportional to the
    # identity so a mistaken transpose/ordering would show up.
    cov = np.diag(np.linspace(1.0, 2.0, COMMANDER_NL))
    cov += 0.05 * np.eye(COMMANDER_NL, k=1) + 0.05 * np.eye(COMMANDER_NL, k=-1)
    mu_sigma = np.full(COMMANDER_NL, 2.5e3)
    return xa, ya, y2a, mu, cov, mu_sigma


def make_commander_cldf(root) -> str:
    """Write a synthetic commander clik tree; returns its path."""
    path = os.path.join(str(root), "commander_synthetic.clik")
    lkl = os.path.join(path, "clik", "lkl_0")
    ext = os.path.join(lkl, "_external")
    os.makedirs(ext, exist_ok=True)

    xa, ya, y2a, mu, cov, mu_sigma = commander_tables()
    cl2x = np.stack([xa, ya, y2a])  # (3, nl, nbin)

    fits.HDUList(
        [
            fits.PrimaryHDU(cl2x),
            fits.ImageHDU(mu),
            fits.ImageHDU(cov),
            fits.ImageHDU(mu_sigma),
        ]
    ).writeto(os.path.join(ext, "sigma.fits"), overwrite=True)

    write_mdb(
        lkl,
        {
            "lkl_type": "gibbs_gauss",
            "lmin": 2,
            "lmax": 29,
            "unit": 1,
            "free_calib": "A_planck",
        },
    )
    write_mdb(os.path.join(path, "clik"), {"lmax": 29, "n_lkl_object": 1})
    write_mdb(path, {"clik_version": "synthetic"})

    return path


# ---------------------------------------------------------------------------
# lensing (clik_lensing itype=4)
# ---------------------------------------------------------------------------

LENSING_LMAX = 120
LENSING_NBINS = 4


def lensing_tables(marged: bool, seed: int = 7):
    """Return the (hascl, bins, cor0, cors, pp_hat, siginv) of the synthetic file."""
    rng = np.random.default_rng(seed)
    nl = LENSING_LMAX + 1
    # CMB-marginalized files carry no CMB renormalization block.
    hascl = (
        np.zeros(6, dtype=np.int64)
        if marged
        else np.array([1, 1, 0, 1, 0, 0], dtype=np.int64)
    )
    ncmb = int(hascl.sum())
    nlt = (1 + ncmb) * nl if ncmb else nl

    bins = np.zeros((LENSING_NBINS, nl))
    edges = np.linspace(8, LENSING_LMAX + 1, LENSING_NBINS + 1).astype(int)
    for ib in range(LENSING_NBINS):
        lo, hi = edges[ib], edges[ib + 1]
        bins[ib, lo:hi] = 1.0 / (hi - lo)

    # The band-power projection of C_l^phiphi is O(1e-7); the renormalization
    # blocks are scaled so their contributions land on the same order, making
    # the model check sensitive to every term instead of just the largest.
    cor0 = 1.0e-8 * rng.random(LENSING_NBINS)
    cors = np.zeros((LENSING_NBINS, nlt))
    cors[:, :nl] = 0.1 * rng.standard_normal((LENSING_NBINS, nl))
    if ncmb:
        cors[:, nl:] = 1.0e-11 * rng.standard_normal((LENSING_NBINS, ncmb * nl))
    pp_hat = 1.0e-7 * (1.0 + rng.random(LENSING_NBINS))
    siginv = np.diag(np.linspace(1.0e14, 2.0e14, LENSING_NBINS))

    return hascl, bins, cor0, cors, pp_hat, siginv


def make_lensing_cldf(root, marged: bool = False, has_calib: bool = True) -> str:
    """Write a synthetic clik_lensing tree; returns its path."""
    name = "lensing_marged_synthetic" if marged else "lensing_synthetic"
    path = os.path.join(str(root), f"{name}.clik_lensing")
    node = os.path.join(path, "clik_lensing")

    hascl, bins, cor0, cors, pp_hat, siginv = lensing_tables(marged)

    write_mdb(
        path,
        {"clik_lensing_version": "synthetic"},
    )
    write_mdb(
        node,
        {
            "itype": 4,
            "lmax": LENSING_LMAX,
            "nbins": LENSING_NBINS,
            "has_calib": int(has_calib),
            "renorm": 1,
            "ren1": 1,
        },
    )
    write_array(node, "hascl", hascl)
    write_array(node, "bins", bins)
    write_array(node, "cor0", cor0)
    write_array(node, "cors", cors)
    write_array(node, "pp_hat", pp_hat)
    write_array(node, "siginv", siginv)

    return path


# ---------------------------------------------------------------------------
# plik_lite (CMB-marginalized high-l bandpowers)
# ---------------------------------------------------------------------------

PLIK_LITE_BIN_WIDTH = 2
PLIK_LITE_LMAX = 500  # covers 215 bins of width 2 starting at ell = 30


def plik_lite_tables(seed: int = 3):
    """Return the (x_all, cov_all, blmin, blmax, bweight) of the synthetic file."""
    # pylint: disable=import-outside-toplevel
    from numcosmo_py.experiments.planck_lite import NBIN_TT, NBIN_TOTAL

    rng = np.random.default_rng(seed)
    blmin = PLIK_LITE_BIN_WIDTH * np.arange(NBIN_TT)
    blmax = blmin + (PLIK_LITE_BIN_WIDTH - 1)
    bweight = np.full(int(blmax[-1]) + 1, 1.0 / PLIK_LITE_BIN_WIDTH)
    # Bandpowers and covariance of the same order as the theory C_l the test
    # cosmology produces, so chi2 stays O(n) instead of overflowing.
    x_all = 1.0 + rng.random(NBIN_TOTAL)
    cov_all = np.diag(np.linspace(0.5, 1.5, NBIN_TOTAL))
    return x_all, cov_all, blmin, blmax, bweight


def make_plik_lite_cldf(root, spectra=("TT",), seed: int = 3) -> str:
    """Write a synthetic plik_lite clik tree; returns its path."""
    path = os.path.join(str(root), "plik_lite_synthetic.clik")
    lkl = os.path.join(path, "clik", "lkl_0")
    ext = os.path.join(lkl, "_external")
    os.makedirs(ext, exist_ok=True)

    # pylint: disable=import-outside-toplevel
    from numcosmo_py.experiments.planck_lite import NBIN_TOTAL

    x_all, cov_all, blmin, blmax, bweight = plik_lite_tables(seed)

    has_cl = np.zeros(6, dtype=np.int64)
    for tag, idx in (("TT", 0), ("EE", 1), ("TE", 3)):
        if tag in spectra:
            has_cl[idx] = 1
    write_array(lkl, "has_cl", has_cl)

    np.savetxt(
        os.path.join(ext, "cl_cmb_plik_v22.dat"),
        np.column_stack([np.arange(NBIN_TOTAL), x_all, x_all]),
    )
    write_fortran_matrix(os.path.join(ext, "c_matrix_plik_v22.dat"), cov_all)
    np.savetxt(os.path.join(ext, "blmin.dat"), blmin, fmt="%d")
    np.savetxt(os.path.join(ext, "blmax.dat"), blmax, fmt="%d")
    np.savetxt(os.path.join(ext, "bweight.dat"), bweight)

    write_mdb(
        lkl,
        {"lkl_type": "plik_cmbonly", "lmin": 30, "lmax": PLIK_LITE_LMAX, "unit": 1},
    )
    write_mdb(os.path.join(path, "clik"), {"lmax": PLIK_LITE_LMAX, "n_lkl_object": 1})
    write_mdb(path, {"clik_version": "synthetic"})

    return path


# ---------------------------------------------------------------------------
# SMICA (plik high-l)
# ---------------------------------------------------------------------------

SMICA_NBINS = 215  # fixed by the converter
SMICA_LMIN = 30
SMICA_LMAX = 2508
SMICA_NELL = SMICA_LMAX - SMICA_LMIN + 1
SMICA_BIN_WIDTH = 11
SMICA_GCIB_NF = 4
SMICA_CN_NF = 12
# Bins active in the criterion mask: keeping only a few keeps the criterion
# matrix (and so the data vector) small without changing any code path.
SMICA_ACTIVE_BINS = (0, 100, 214)


def smica_binning():
    """Return the (bin_lmin, bin_lmax, bin_ws) of the synthetic binning.

    Contiguous bins of a fixed width, each with flat weights summing to one, so
    a bandpower is the plain average of the theory over the bin.
    """
    bin_lmin = SMICA_BIN_WIDTH * np.arange(SMICA_NBINS)
    bin_lmax = bin_lmin + (SMICA_BIN_WIDTH - 1)
    assert bin_lmax[-1] < SMICA_NELL
    bin_ws = np.full(SMICA_NBINS * SMICA_BIN_WIDTH, 1.0 / SMICA_BIN_WIDTH)
    return bin_lmin, bin_lmax, bin_ws


def smica_templates(pol: bool):
    """Return the synthetic foreground templates, keyed by component directory.

    Only the shapes and the positivity where the assembly divides by a template
    matter; the CMB-only check below zeroes every foreground amplitude, and the
    full-amplitude check only requires the result to stay finite.
    """
    nl_cn = SMICA_LMAX + 1
    # Beam leakage has no amplitude parameter (clik fixes it to one), so a zero
    # template is the only way to switch it off for the closed-form check.
    leak = np.zeros(nl_cn * SMICA_CN_NF * SMICA_CN_NF)
    # Strictly positive: the galactic dust component divides by the template at
    # its pivot multipole.
    dust = np.full(nl_cn * SMICA_CN_NF * SMICA_CN_NF, 2.0e-4)
    sbpx = np.full(nl_cn * SMICA_CN_NF * SMICA_CN_NF, 1.0e-6)
    gcib = np.full(3001 * SMICA_GCIB_NF * SMICA_GCIB_NF, 3.0e-2)
    sz = np.full(9999, 5.0)  # indexed [ell - 2], normalized at ell = 3000
    ksz = np.full(5001, 3.0)  # indexed [ell],     normalized at ell = 3000
    gibxsz = np.full(9998, -1.0)  # indexed [ell - 2]

    tmpl = {
        "component_1": gcib,
        "component_2": gibxsz,
        "component_3": sz,
        "component_5": ksz,
    }
    if pol:
        tmpl.update(
            {
                "component_6": dust,
                "component_9": leak,
                "component_10": np.full_like(dust, 4.0e-7),  # EE e2e cnoise
                "component_11": sbpx,
            }
        )
    else:
        tmpl.update({"component_6": dust, "component_7": leak, "component_8": sbpx})
    return tmpl


def smica_mask_and_ordering(m: int):
    """Return the (mask, ordering, quad) of the synthetic criterion selection."""
    ordering = np.array([(i, j) for i in range(m) for j in range(i, m)], dtype=np.int64)
    mask = np.zeros((SMICA_NBINS, m, m), dtype=np.int64)
    for iq in SMICA_ACTIVE_BINS:
        for i, j in ordering:
            mask[iq, i, j] = 1

    quad = np.array(
        [
            iq * m * m + i * m + j
            for i, j in ordering
            for iq in range(SMICA_NBINS)
            if mask[iq, i, j]
        ],
        dtype=np.int64,
    )
    return mask, ordering, quad


def make_smica_cldf(root, pol: bool = False, seed: int = 5) -> str:
    """Write a synthetic plik (SMICA) clik tree; returns its path."""
    path = os.path.join(str(root), f"plik_{'ttteee' if pol else 'tt'}_synthetic.clik")
    lkl = os.path.join(path, "clik", "lkl_0")
    m = 6 if pol else 3
    rng = np.random.default_rng(seed)

    bin_lmin, bin_lmax, bin_ws = smica_binning()
    mask, ordering, quad = smica_mask_and_ordering(m)
    npt = quad.size

    rq_hat = 1.0 + rng.random(SMICA_NBINS * m * m)
    crit = np.diag(np.linspace(0.5, 1.5, npt))

    write_mdb(
        lkl,
        {
            "lkl_type": "smica",
            "criterion": "gauss",
            "lmin": SMICA_LMIN,
            "lmax": SMICA_LMAX,
            "nbins": SMICA_NBINS,
            "m_channel_T": 3,
            "m_channel_P": 3 if pol else 0,
            "n_component": 12 if pol else 11,
            "unit": 1,
        },
    )
    write_array(lkl, "has_cl", np.array([1, 1, 0, 1, 0, 0], dtype=np.int64))
    write_array(lkl, "Rq_hat", rq_hat)
    write_array(lkl, "criterion_gauss_mat", crit.ravel())
    write_array(lkl, "criterion_gauss_mask", mask.ravel())
    write_array(lkl, "criterion_gauss_ordering", ordering.ravel())
    write_array(lkl, "A_cmb", np.ones(m))
    write_array(lkl, "wq", np.ones(SMICA_NBINS))

    if pol:
        # The TTTEEE file bins TT, TE and EE with the same relative operator; the
        # converter keeps the first block only, so ship all three.
        write_array(lkl, "bin_lmin", np.tile(bin_lmin, 3))
        write_array(lkl, "bin_lmax", np.tile(bin_lmax, 3))
        write_array(lkl, "bin_ws", np.tile(bin_ws, 3))
    else:
        write_array(lkl, "bin_lmin", bin_lmin)
        write_array(lkl, "bin_lmax", bin_lmax)
        write_array(lkl, "bin_ws", bin_ws)

    for comp, template in smica_templates(pol).items():
        write_array(os.path.join(lkl, comp), "template", template)

    if pol:
        # icalTP: 100T/217T plus the three polarization calibrations, each map
        # mixing with itself only (single-term weights).
        node = os.path.join(lkl, "component_12")
        w = np.zeros(2 * m * m)
        w[0::2] = 1.0
        other = np.zeros(2 * m * m, dtype=np.int64)
        for i in range(m):
            for j in range(m):
                other[(i * m + j) * 2] = i
                other[(i * m + j) * 2 + 1] = j
        write_array(node, "im", np.array([0, 2, 3, 4, 5], dtype=np.int64))
        write_array(node, "w", w)
        write_array(node, "other", other)
        write_mdb(node, {"component_type": "icalTP"})

    write_mdb(os.path.join(path, "clik"), {"lmax": SMICA_LMAX, "n_lkl_object": 1})
    write_mdb(path, {"clik_version": "synthetic"})

    return path


# ---------------------------------------------------------------------------
# shared model-set / theory helpers
# ---------------------------------------------------------------------------


def planck_mset(pol: bool = False, **params):
    """Build the (NcPlanckFI, NcHICosmo) set the synthetic likelihoods evaluate on.

    @pol selects the TTTEEE nuisance model; @params overrides individual
    nuisance parameters by name. Returns the (mset, planck model) pair.
    """
    # pylint: disable=import-outside-toplevel
    from numcosmo_py import Ncm, Nc
    from numcosmo_py.cosmology import create_cosmo, HIPrimModel

    cosmo = create_cosmo(prim_model=HIPrimModel.POWER_LAW)
    planck = Nc.PlanckFICorTTTEEE() if pol else Nc.PlanckFICorTT()
    planck.params_set_default_ftype()
    for name, value in params.items():
        planck[name] = value

    mset = Ncm.MSet.new_array([planck, cosmo])
    mset.prepare_fparam_map()

    return mset, planck


def theory_cls(pb, spectrum: str, lmax: int) -> np.ndarray:
    """Return the prepared theory $C_\\ell$ of @spectrum as a numpy array."""
    # pylint: disable=import-outside-toplevel
    from numcosmo_py import Ncm

    vec = Ncm.Vector.new(lmax + 1)
    getattr(pb, f"get_{spectrum}_Cls")(vec)

    return np.array([vec.get(i) for i in range(lmax + 1)])


def model_vector(data, mset) -> np.ndarray:
    """Return the model (mean) vector a Gaussian likelihood assembles."""
    # pylint: disable=import-outside-toplevel
    from numcosmo_py import Ncm

    size = data.get_length()
    vec = Ncm.Vector.new(size)
    data.mean_vector(mset, vec)

    return np.array([vec.get(i) for i in range(size)])

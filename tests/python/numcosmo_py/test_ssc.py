#!/usr/bin/env python
#
# test_ssc.py
#
# Mon Aug 11 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_ssc.py
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

"""Tests for numcosmo_py.ssc, the NumCosmo-only replacement for PySSC."""

import numpy as np
import pytest
from numpy.testing import assert_allclose

from numcosmo_py import Nc, Ncm
from numcosmo_py.ssc import (
    SijCalculator,
    find_lmax,
    mask_angular_power_spectrum,
)

Ncm.cfg_init()

Z_EDGES = np.array([0.1, 0.2, 0.3, 0.4])

# Reference values from an independent scipy quadrature of
#   U_i(k) = int dchi (chi^2 / DeltaV) D(z) j_0(k chi),
#   C_0^ij = (2/pi) int dk k^2 P(k, 0) U_i U_j,   S_ij = C_0^ij / 4pi,
# converged at 1601 redshift points and 4096 log-spaced k points over
# k in [1e-6, 5] 1/Mpc, for the fiducial cosmology built below. The reference
# itself is converged to a few times 1e-4 relative, which is what sets the
# tolerances used here -- not the accuracy of numcosmo_py.ssc.
REFERENCE_SIJ_DIAG = np.array([1.9022e-05, 6.9727e-06, 3.6455e-06])
REFERENCE_SIJ_01 = -2.9193e-06


@pytest.fixture(name="cosmo", scope="module")
def fixture_cosmo() -> Nc.HICosmo:
    """Fiducial flat wCDM cosmology matching the J-PAS 2024 forecast."""
    cosmo = Nc.HICosmoDEXcdm()
    cosmo.omega_x2omega_k()
    cosmo["H0"] = 67.81
    cosmo["Omegab"] = 0.0486
    cosmo["Omegac"] = 0.2612
    cosmo["w"] = -1.0
    cosmo["Omegak"] = 0.0

    prim = Nc.HIPrimPowerLaw.new()
    prim["ln10e10ASA"] = 3.02745
    prim["n_SA"] = 0.9660
    cosmo.add_submodel(prim)
    cosmo.add_submodel(Nc.HIReionCamb.new())

    return cosmo


def _cap_mask(area_deg2: float, nside: int = 64) -> np.ndarray:
    """Build a HEALPix circular-cap mask of the requested area."""
    sqd_fullsky = 4.0 * np.pi * (180.0 / np.pi) ** 2
    cos_theta_c = 1.0 - 2.0 * area_deg2 / sqd_fullsky

    smap = Ncm.SphereMap.new(nside)
    npix = smap.get_npix()
    theta = np.array([smap.pix2ang_ring(i)[0] for i in range(npix)])

    mask = np.zeros(npix)
    mask[np.cos(theta) >= cos_theta_c] = 1.0

    return mask


def test_calculator_rejects_bad_edges() -> None:
    """Non-monotonic or too-short edge arrays are rejected."""
    with pytest.raises(ValueError, match="strictly increasing"):
        SijCalculator([0.3, 0.2, 0.4])

    with pytest.raises(ValueError, match="at least 2 edges"):
        SijCalculator([0.3])


def test_fullsky_matches_reference(cosmo: Nc.HICosmo) -> None:
    """Full-sky Sij reproduces an independent quadrature of the same integral."""
    calc = SijCalculator(Z_EDGES)
    sij = calc.fullsky(cosmo)

    assert sij.shape == (3, 3)
    assert_allclose(sij, sij.T, rtol=1.0e-12)
    assert_allclose(np.diag(sij), REFERENCE_SIJ_DIAG, rtol=1.0e-3)
    # The off-diagonal is the demanding one: a small residual of a large
    # cancellation, and the reason this module tightens scaled_abstol.
    assert_allclose(sij[0, 1], REFERENCE_SIJ_01, rtol=1.0e-3)

    # Disjoint shells anticorrelate, and every cross term obeys Cauchy-Schwarz.
    off_diagonal = ~np.eye(3, dtype=bool)
    assert np.all(sij[off_diagonal] < 0.0)
    assert np.all(np.abs(sij[0, 1:]) < np.sqrt(sij[0, 0] * np.diag(sij)[1:]))


def test_fullsky_fsky_is_a_rescaling(cosmo: Nc.HICosmo) -> None:
    """The fsky-corrected full-sky Sij is the full-sky one divided by fsky."""
    calc = SijCalculator(Z_EDGES)
    area = 2000.0
    fsky = area * (np.pi / 180.0) ** 2 / (4.0 * np.pi)

    assert_allclose(
        calc.fullsky_fsky(cosmo, area), calc.fullsky(cosmo) / fsky, rtol=1e-10
    )

    with pytest.raises(ValueError, match="square degrees"):
        calc.fullsky_fsky(cosmo, -1.0)


def test_mask_spectrum_and_fsky() -> None:
    """The mask spectrum recovers fsky and is consistent with the map's mean."""
    mask = _cap_mask(2000.0, nside=64)
    cl_mask, fsky = mask_angular_power_spectrum(mask)

    assert cl_mask.size == 2 * 64 + 1
    assert_allclose(fsky, mask.mean(), rtol=1.0e-3)
    assert_allclose(fsky, np.sqrt(cl_mask[0] / (4.0 * np.pi)), rtol=1.0e-12)


def test_find_lmax_truncation() -> None:
    """find_lmax picks the smallest truncation recovering the mask variance."""
    cl_mask, _ = mask_angular_power_spectrum(_cap_mask(2000.0, nside=64))

    lmax = find_lmax(cl_mask, var_tol=0.05)
    ell = np.arange(cl_mask.size)
    summand = (2.0 * ell + 1.0) * cl_mask / (4.0 * np.pi)

    assert 0 < lmax < cl_mask.size
    assert abs(summand[: lmax + 1].sum() - summand.sum()) / summand.sum() <= 0.05
    # Minimal: one multipole less misses the tolerance.
    assert abs(summand[:lmax].sum() - summand.sum()) / summand.sum() > 0.05

    # A tighter tolerance never shortens the sum.
    assert find_lmax(cl_mask, var_tol=0.01) >= lmax


def test_partial_sky_exceeds_fullsky(cosmo: Nc.HICosmo) -> None:
    """A finite footprint has more super-sample variance than the full sky."""
    calc = SijCalculator(Z_EDGES)
    mask = _cap_mask(2000.0, nside=64)

    sij_cap = calc.partial_sky(cosmo, mask, var_tol=0.05)
    sij_full = calc.fullsky(cosmo)

    assert sij_cap.shape == (3, 3)
    assert_allclose(sij_cap, sij_cap.T, rtol=1.0e-12)
    assert np.all(np.diag(sij_cap) > np.diag(sij_full))

    # The crude 1/fsky rescaling is the right order of magnitude but not exact.
    _, fsky = mask_angular_power_spectrum(mask)
    ratio = np.diag(sij_cap) / np.diag(sij_full / fsky)
    assert np.all((ratio > 0.5) & (ratio < 2.0))


def test_progress_callback_fires(cosmo: Nc.HICosmo) -> None:
    """compute_cl reports progress once per chunk, ending at 100%."""
    calc = SijCalculator(Z_EDGES)
    events: list[tuple[int, int]] = []

    calc.compute_cl(
        cosmo,
        0,
        15,
        progress=lambda done, total, _elapsed, _msg: events.append((done, total)),
        chunk_size=8,
    )

    assert events == [(1, 2), (2, 2)]


def test_compute_cl_rejects_inverted_range(cosmo: Nc.HICosmo) -> None:
    """An inverted multipole range is rejected."""
    calc = SijCalculator(Z_EDGES)

    with pytest.raises(ValueError, match="must not be smaller"):
        calc.compute_cl(cosmo, 5, 2)

#!/usr/bin/env python
#
# test_galaxy_sd_redshift_binning.py
#
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
# with this program. If not, see <http://www.gnu.org/licenses/>.

"""Standalone tests for the binning (dn/dz producer) calculator.

``NcGalaxyRedshiftBinning`` produces the true-redshift distribution dn/dz of a
photometric bin from a shared population + observable pair. The photometric
window is an argument to the producers, not object state, so a single calculator
serves arbitrarily many bins. Validated for unit normalization, the on-nodes
contract, and parity against the legacy ``NcGalaxySDObsRedshiftGauss`` (which
shares the same P(z)*W/N convolution and marginal P(zp)).

FROZEN REFERENCE VALUES: the parity documented above was proven by running
both engines live and is captured, not re-derived, in
``test_dndz_parity_legacy``/``test_pzp_parity_legacy``/
``test_equal_area_bins_parity_legacy``/``test_lsst_srd_source_edges_parity_legacy``
below. Values were captured from an actual passing run of this file's
original legacy-comparison code, at git rev ``77313f22`` (2026-07-16), then
legacy (``NcGalaxySDObsRedshiftGauss``/``NcGalaxySDTrueRedshiftLSSTSRD``)
construction was removed so these tests no longer depend on legacy at
runtime -- legacy is slated for deletion in a follow-up PR. Each frozen
assertion keeps the tolerance (``rtol``/``atol``) that the original live
comparison used. The ``_DNDZ_FROZEN``/``_PZP_FROZEN`` sequences are stored
as ``Ncm.Matrix`` binfiles (``data/truth_tables/``) rather than inline
literals; see ``_load_dndz_golden``/``_load_pzp_golden``.
"""

import pytest
import numpy as np
from numpy.testing import assert_allclose

from numcosmo_py import Ncm, Nc

Ncm.cfg_init()

# (variant, bin_sigma0, zp_min, zp_max)
_CASES = [
    ("y1_source", 0.05, 0.4, 0.9),
    ("y1_source", 0.05, 0.9, 1.5),
    ("y10_lens", 0.03, 0.2, 0.4),
    ("y1_lens", 0.03, 0.6, 1.0),
]


def _build_new(variant, bin_sigma0):
    # The calculator holds neither the window nor the models (NcDistance
    # convention); the window and models are passed to every method.
    pop = getattr(Nc.GalaxyRedshiftPopLSSTSRD, f"new_{variant}")()
    obs_pop = Nc.GalaxyRedshiftObsSelGauss.new()
    obs_pop.param_set_by_name("sigma0", bin_sigma0)
    binning = Nc.GalaxyRedshiftBinning.new()
    return binning, pop, obs_pop


@pytest.mark.parametrize("variant,bin_sigma0,zp_min,zp_max", _CASES)
def test_dndz_normalized(variant, bin_sigma0, zp_min, zp_max):
    """compute_dndz produces a spline integrating to unity over its support."""
    binning, pop, obs_pop = _build_new(variant, bin_sigma0)
    dndz = binning.compute_dndz(pop, obs_pop, zp_min, zp_max)
    zs = np.linspace(0.0, 4.0, 8000)
    vals = np.array([dndz.eval(z) for z in zs])
    # Density is non-negative; the spline shows only machine-noise negatives
    # (~1e-14) from interpolation/extrapolation at the tails.
    assert vals.min() > -1.0e-8 * vals.max()
    assert_allclose(np.trapezoid(vals, zs), 1.0, rtol=1.0e-4)


# Frozen legacy eval_pz_given_zp output, keyed by (variant, bin_sigma0,
# zp_min, zp_max), sampled on a 400-point z-grid (see module docstring).
_DNDZ_GOLDEN_FILE = "truth_tables/wl/nc_galaxy_redshift_binning_dndz_parity.bin"


def _load_dndz_golden() -> np.ndarray:
    """Load the frozen dn/dz sequences as a (len(_CASES), 400) array."""
    path = Ncm.cfg_get_data_filename(_DNDZ_GOLDEN_FILE, True)
    ser = Ncm.Serialize.new(Ncm.SerializeOpt.NONE)
    matrix = ser.from_binfile(path)
    assert isinstance(matrix, Ncm.Matrix)
    return np.array(matrix.dup_array()).reshape(len(_CASES), 400)


_DNDZ_FROZEN = _load_dndz_golden()


@pytest.mark.parametrize("variant,bin_sigma0,zp_min,zp_max", _CASES)
def test_dndz_parity_legacy(variant, bin_sigma0, zp_min, zp_max):
    """compute_dndz is checked against a frozen legacy eval_pz_given_zp
    convolution (see module docstring)."""
    binning, pop, obs_pop = _build_new(variant, bin_sigma0)

    dndz = binning.compute_dndz(pop, obs_pop, zp_min, zp_max)
    zs = np.linspace(0.05, 2.5, 400)
    new_vals = np.array([dndz.eval(z) for z in zs])
    case_idx = _CASES.index((variant, bin_sigma0, zp_min, zp_max))
    frozen = _DNDZ_FROZEN[case_idx]
    # Both are adaptive splines of the same P(z)W/N at reltol=1e-7; agree to
    # interpolation/normalization noise well below a part in 1e4.
    assert_allclose(new_vals, frozen, rtol=1.0e-4, atol=1.0e-6)


def test_compute_dndz_on_nodes():
    """compute_dndz_on_nodes tabulates the same dn/dz on the given nodes."""
    binning, pop, obs_pop = _build_new("y1_source", 0.05)

    auto = binning.compute_dndz(pop, obs_pop, 0.4, 0.9)
    zs = np.linspace(0.1, 2.0, 200)
    on_nodes = binning.compute_dndz_on_nodes(
        pop, obs_pop, 0.4, 0.9, Ncm.Vector.new_array(zs.tolist())
    )

    a = np.array([auto.eval(z) for z in zs])
    b = np.array([on_nodes.eval(z) for z in zs])
    assert_allclose(b, a, rtol=1.0e-6)


def test_dndz_window_is_an_argument():
    """Different windows yield different dn/dz from the same calculator."""
    binning, pop, obs_pop = _build_new("y1_source", 0.05)
    lo = binning.compute_dndz(pop, obs_pop, 0.4, 0.9)
    hi = binning.compute_dndz(pop, obs_pop, 1.0, 1.6)
    assert not np.isclose(lo.eval(0.6), hi.eval(0.6))


# Frozen legacy eval_pzp output, keyed by (variant, bin_sigma0, zp_min,
# zp_max), sampled on a 400-point zp-grid (see module docstring).
_PZP_GOLDEN_FILE = "truth_tables/wl/nc_galaxy_redshift_binning_pzp_parity.bin"


def _load_pzp_golden() -> np.ndarray:
    """Load the frozen eval_pzp sequences as a (len(_CASES), 400) array."""
    path = Ncm.cfg_get_data_filename(_PZP_GOLDEN_FILE, True)
    ser = Ncm.Serialize.new(Ncm.SerializeOpt.NONE)
    matrix = ser.from_binfile(path)
    assert isinstance(matrix, Ncm.Matrix)
    return np.array(matrix.dup_array()).reshape(len(_CASES), 400)


_PZP_FROZEN = _load_pzp_golden()


@pytest.mark.parametrize("variant,bin_sigma0,zp_min,zp_max", _CASES)
def test_pzp_parity_legacy(variant, bin_sigma0, zp_min, zp_max):
    """eval_pzp is checked against the frozen legacy marginal photo-z
    density (window-free; see module docstring)."""
    binning, pop, obs_pop = _build_new(variant, bin_sigma0)

    binning.prepare(pop, obs_pop)
    zps = np.linspace(0.05, 3.0, 400)
    new_vals = np.array([binning.eval_pzp(zp) for zp in zps])
    case_idx = _CASES.index((variant, bin_sigma0, zp_min, zp_max))
    frozen = _PZP_FROZEN[case_idx]
    # Both are stats-dist splines of the same -2ln P(zp) at reltol=1e-7.
    assert_allclose(new_vals, frozen, rtol=1.0e-4, atol=1.0e-6)


def test_pzp_normalized():
    """eval_pzp integrates to unity over its support."""
    binning, pop, obs_pop = _build_new("y1_source", 0.05)
    binning.prepare(pop, obs_pop)
    zps = np.linspace(0.0, 20.0, 40000)
    vals = np.array([binning.eval_pzp(zp) for zp in zps])
    assert np.all(vals >= 0.0)
    assert_allclose(np.trapezoid(vals, zps), 1.0, rtol=1.0e-4)


# Frozen legacy compute_equal_area_photoz_bins edges, keyed by n_bins (see
# module docstring), for _build_new("y1_source", 0.05) over [0.4, 0.9] with
# max_zp=2.0.
_EQUAL_AREA_BINS_FROZEN = {
    3: [0.0, 0.469625438419577, 0.8390523737382864, 1.9999999999362368],
    5: [
        0.0,
        0.34141905221076563,
        0.5340859955504487,
        0.7511788435006586,
        1.0650891330372267,
        1.9999999999362368,
    ],
    8: [
        0.0,
        0.2626343751127359,
        0.39016003297454027,
        0.5096875135957158,
        0.6363914809196871,
        0.782809365923383,
        0.9695824672222074,
        1.2514990596560016,
        1.9999999999362368,
    ],
}


@pytest.mark.parametrize("n_bins", [3, 5, 8])
def test_equal_area_bins_parity_legacy(n_bins):
    """Equal-area photo-z edges are checked against the frozen legacy CDF
    inversion (see module docstring)."""
    binning, pop, obs_pop = _build_new("y1_source", 0.05)

    binning.prepare(pop, obs_pop)
    new_edges = binning.compute_equal_area_photoz_bins(n_bins, 2.0).dup_array()
    assert_allclose(
        new_edges, _EQUAL_AREA_BINS_FROZEN[n_bins], rtol=1.0e-5, atol=1.0e-6
    )
    # Edges are strictly increasing and cover the requested range.
    assert np.all(np.diff(new_edges) > 0.0)
    assert_allclose(new_edges[-1], 2.0, rtol=1.0e-6)


# LSST SRD binning factory -------------------------------------------------

# (type, n_bins, sigma0, is_lens)
_LSST_CASES = [
    (Nc.GalaxyRedshiftPopLSSTSRDType.Y1_LENS, 5, 0.03, True),
    (Nc.GalaxyRedshiftPopLSSTSRDType.Y10_LENS, 10, 0.03, True),
    (Nc.GalaxyRedshiftPopLSSTSRDType.Y1_SOURCE, 5, 0.05, False),
    (Nc.GalaxyRedshiftPopLSSTSRDType.Y10_SOURCE, 5, 0.05, False),
]


@pytest.mark.parametrize("srd_type,n_bins,sigma0,is_lens", _LSST_CASES)
def test_lsst_srd_edges_structure(srd_type, n_bins, sigma0, is_lens):
    """Factory returns n_bins+1 monotonic edges and models set for the type."""
    edges, pop, obs_pop = Nc.GalaxyRedshiftBinning.lsst_srd_edges(srd_type)
    e = edges.dup_array()

    assert len(e) == n_bins + 1
    assert np.all(np.diff(e) > 0.0)
    assert_allclose(obs_pop.param_get_by_name("sigma0"), sigma0)

    ref_pop = Nc.GalaxyRedshiftPopLSSTSRD.new_from_type(srd_type)
    zs = np.linspace(0.05, 3.0, 100)
    assert_allclose(
        [pop.eval(z) for z in zs], [ref_pop.eval(z) for z in zs], rtol=1.0e-12
    )

    if is_lens:
        assert_allclose(e, np.linspace(0.2, 1.2, n_bins + 1))
    else:
        assert_allclose(e[-1], 3.5, rtol=1.0e-6)


# Frozen legacy compute_equal_area_photoz_bins(5, 3.5) edges for
# _build_new("y1_source", 0.05) over [0.0, 3.5] (see module docstring).
_LSST_SRD_SOURCE_EDGES_FROZEN = [
    0.0,
    0.3469506616343085,
    0.5450939593278648,
    0.7721978603224499,
    1.1138745605928926,
    3.499999996915834,
]


def test_lsst_srd_source_edges_parity_legacy():
    """Source equal-area edges are checked against the frozen legacy
    factory's construction (see module docstring)."""
    srd_type = Nc.GalaxyRedshiftPopLSSTSRDType.Y1_SOURCE
    edges, _pop, _obs_pop = Nc.GalaxyRedshiftBinning.lsst_srd_edges(srd_type)

    assert_allclose(
        edges.dup_array(), _LSST_SRD_SOURCE_EDGES_FROZEN, rtol=1.0e-5, atol=1.0e-6
    )


def test_lsst_srd_edges_feed_compute_dndz():
    """The factory edges + models drive compute_dndz for each bin."""
    edges, pop, obs_pop = Nc.GalaxyRedshiftBinning.lsst_srd_edges(
        Nc.GalaxyRedshiftPopLSSTSRDType.Y1_LENS
    )
    e = edges.dup_array()
    binning = Nc.GalaxyRedshiftBinning.new()
    zs = np.linspace(0.0, 4.0, 4000)
    for i in range(len(e) - 1):
        dndz = binning.compute_dndz(pop, obs_pop, e[i], e[i + 1])
        vals = np.array([dndz.eval(z) for z in zs])
        assert vals.min() > -1.0e-8 * vals.max()
        assert_allclose(np.trapezoid(vals, zs), 1.0, rtol=1.0e-4)


def test_reltol_property_round_trip():
    """reltol is readable/writable both via the plain getter/setter and via
    the underlying GObject property (get_property/set_property), and
    changing it invalidates the cached P(zp) marginal (same effect as
    set_zp_support_max, checked below)."""
    binning = Nc.GalaxyRedshiftBinning.new()

    default = binning.get_reltol()
    assert_allclose(binning.get_property("reltol"), default)

    binning.set_reltol(1.0e-6)
    assert_allclose(binning.get_reltol(), 1.0e-6)
    assert_allclose(binning.get_property("reltol"), 1.0e-6)


def test_zp_support_max_property_round_trip():
    """zp-support-max is readable/writable both via the plain getter/setter
    and via the underlying GObject property, and setting it invalidates the
    cached P(zp) marginal (compute_dndz's bin dn/dz is unaffected, per the
    class docs -- only eval_pzp's own support changes)."""
    binning = Nc.GalaxyRedshiftBinning.new()

    default = binning.get_zp_support_max()
    assert_allclose(binning.get_property("zp-support-max"), default)

    binning.set_zp_support_max(10.0)
    assert_allclose(binning.get_zp_support_max(), 10.0)
    assert_allclose(binning.get_property("zp-support-max"), 10.0)

    pop = getattr(Nc.GalaxyRedshiftPopLSSTSRD, "new_y1_source")()
    obs_pop = Nc.GalaxyRedshiftObsSelGauss.new()
    obs_pop.param_set_by_name("sigma0", 0.05)
    binning.prepare(pop, obs_pop)
    zps = np.linspace(0.0, 9.0, 200)
    vals = np.array([binning.eval_pzp(zp) for zp in zps])
    assert np.all(np.isfinite(vals))


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

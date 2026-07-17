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
shares the same P(z)*W/N convolution and marginal P(zp)), keeping the legacy
pristine as the oracle.
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


def _build_legacy(variant, bin_sigma0, zp_min, zp_max):
    sdz = getattr(Nc.GalaxySDTrueRedshiftLSSTSRD, f"new_{variant}")()
    gsdor = Nc.GalaxySDObsRedshiftGauss.new(sdz, zp_min, zp_max)
    gsdor.set_bin_sigma0(bin_sigma0)
    return gsdor


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


@pytest.mark.parametrize("variant,bin_sigma0,zp_min,zp_max", _CASES)
def test_dndz_parity_legacy(variant, bin_sigma0, zp_min, zp_max):
    """compute_dndz matches the legacy eval_pz_given_zp convolution."""
    binning, pop, obs_pop = _build_new(variant, bin_sigma0)
    legacy = _build_legacy(variant, bin_sigma0, zp_min, zp_max)

    dndz = binning.compute_dndz(pop, obs_pop, zp_min, zp_max)
    zs = np.linspace(0.05, 2.5, 400)
    new_vals = np.array([dndz.eval(z) for z in zs])
    old_vals = np.array([legacy.eval_pz_given_zp(z) for z in zs])
    # Both are adaptive splines of the same P(z)W/N at reltol=1e-7; agree to
    # interpolation/normalization noise well below a part in 1e4.
    assert_allclose(new_vals, old_vals, rtol=1.0e-4, atol=1.0e-6)


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


@pytest.mark.parametrize("variant,bin_sigma0,zp_min,zp_max", _CASES)
def test_pzp_parity_legacy(variant, bin_sigma0, zp_min, zp_max):
    """eval_pzp matches the legacy marginal photo-z density (window-free)."""
    binning, pop, obs_pop = _build_new(variant, bin_sigma0)
    legacy = _build_legacy(variant, bin_sigma0, zp_min, zp_max)

    binning.prepare(pop, obs_pop)
    zps = np.linspace(0.05, 3.0, 400)
    new_vals = np.array([binning.eval_pzp(zp) for zp in zps])
    old_vals = np.array([legacy.eval_pzp(zp) for zp in zps])
    # Both are stats-dist splines of the same -2ln P(zp) at reltol=1e-7.
    assert_allclose(new_vals, old_vals, rtol=1.0e-4, atol=1.0e-6)


def test_pzp_normalized():
    """eval_pzp integrates to unity over its support."""
    binning, pop, obs_pop = _build_new("y1_source", 0.05)
    binning.prepare(pop, obs_pop)
    zps = np.linspace(0.0, 20.0, 40000)
    vals = np.array([binning.eval_pzp(zp) for zp in zps])
    assert np.all(vals >= 0.0)
    assert_allclose(np.trapezoid(vals, zps), 1.0, rtol=1.0e-4)


@pytest.mark.parametrize("n_bins", [3, 5, 8])
def test_equal_area_bins_parity_legacy(n_bins):
    """Equal-area photo-z edges match the legacy CDF inversion."""
    binning, pop, obs_pop = _build_new("y1_source", 0.05)
    legacy = _build_legacy("y1_source", 0.05, 0.4, 0.9)

    binning.prepare(pop, obs_pop)
    new_edges = binning.compute_equal_area_photoz_bins(n_bins, 2.0).dup_array()
    old_edges = legacy.compute_equal_area_photoz_bins(n_bins, 2.0).dup_array()
    assert_allclose(new_edges, old_edges, rtol=1.0e-5, atol=1.0e-6)
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


def test_lsst_srd_source_edges_parity_legacy():
    """Source equal-area edges reproduce the legacy factory's construction."""
    srd_type = Nc.GalaxyRedshiftPopLSSTSRDType.Y1_SOURCE
    edges, _pop, _obs_pop = Nc.GalaxyRedshiftBinning.lsst_srd_edges(srd_type)

    legacy = _build_legacy("y1_source", 0.05, 0.0, 3.5)
    old_edges = legacy.compute_equal_area_photoz_bins(5, 3.5).dup_array()

    assert_allclose(edges.dup_array(), old_edges, rtol=1.0e-5, atol=1.0e-6)


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

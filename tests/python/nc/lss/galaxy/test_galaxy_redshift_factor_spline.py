#!/usr/bin/env python
#
# test_galaxy_redshift_factor_spline.py
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

"""Standalone tests for the Spline redshift calculator scheme.

``NcGalaxyRedshiftFactorSpline`` hands out a per-galaxy pre-tabulated
$p(z)$ as-is: no ``NcmMSet`` model is involved, so there is no legacy
class sharing identical math to use as a parity oracle the way
``test_galaxy_redshift_factor_composed.py`` does. Instead, each test here
builds the spline itself from the *Composed* scheme's own joint density
(round-tripping a real, physically meaningful $p(z)$ through a spline) and
checks that Spline's ``integ``/``get_integ_lim``/``norm``/``gen`` correctly
reproduce that spline -- i.e. this validates the Spline scheme's own
correctness as a spline wrapper, not a specific physics result.
"""

import pytest
import numpy as np
from numpy.testing import assert_allclose

from numcosmo_py import Ncm, Nc

Ncm.cfg_init()

# (variant, zp, sigma0, zp_min, zp_max)
_CASES = [
    ("y1_source", 0.6, 0.05, 0.0, 20.0),  # wide-open window -> untruncated
    ("y1_source", 0.6, 0.05, 0.4, 0.9),  # narrow selection window
    ("y10_lens", 1.2, 0.08, 0.8, 1.7),
    ("y1_lens", 0.3, 0.03, 0.1, 0.5),
]

_N_GRID = 400


def _build_composed(variant, zp, sigma0, zp_min, zp_max):
    """Composed scheme + a per-galaxy data with (zp, sigma0) populated."""
    pop = getattr(Nc.GalaxyRedshiftPopLSSTSRD, f"new_{variant}")()
    obs = Nc.GalaxyRedshiftObsGauss.new()
    mset = Ncm.MSet.empty_new()
    mset.set(pop)
    mset.set(obs)
    composed = Nc.GalaxyRedshiftFactorComposed.new(zp_min, zp_max)

    data = Nc.GalaxyRedshiftFactorData.new(composed, mset)
    cols = Nc.GalaxyRedshiftFactorData.required_columns(data)
    wlobs = Nc.GalaxyWLObs.new(
        Nc.GalaxyWLObsEllipConv.TRACE_DET, Nc.WLEllipticityFrame.CELESTIAL, 1, cols
    )
    wlobs.set("z", 0, 0.0)
    wlobs.set("zp", 0, zp)
    wlobs.set("sigma0", 0, sigma0)
    data.read_row(wlobs, 0)
    return composed, mset, data


def _spline_from_composed(variant, zp, sigma0, zp_min, zp_max, n=_N_GRID):
    """Samples the Composed joint density on a grid and builds a spline of it.

    Returns the spline together with the grid and the Composed values on it,
    so callers can compare Spline's output against both the spline itself
    (exact, at the grid nodes) and the underlying physical density (close,
    within interpolation error).

    Trims a small margin off Composed's own get_integ_lim(): that support
    edge is where the *population*'s density is genuinely, exactly zero by
    construction, and NcmStatsDist1dSpline's ODE-based sampler (shared math
    with the legacy Pz class) cannot start from a zero-density edge. A real
    photo-z p(z) spline (e.g. from a photo-z code) would not sample its
    edge exactly at zero density either.
    """
    composed, mset, data = _build_composed(variant, zp, sigma0, zp_min, zp_max)
    z_lo, z_hi = composed.get_integ_lim(mset, data)
    margin = 1.0e-3 * (z_hi - z_lo)
    z_min, z_max = z_lo + margin, z_hi - margin
    integ = composed.integ(mset, False)

    zs = np.linspace(z_min, z_max, n)
    pz_vals = np.array([integ.eval(z, data) for z in zs])

    xv = Ncm.Vector.new_array(zs.tolist())
    yv = Ncm.Vector.new_array(pz_vals.tolist())
    spline = Ncm.SplineCubicNotaknot.new()
    spline.set(xv, yv, True)

    return spline, zs, pz_vals


def _build_spline_scheme(spline):
    """Spline scheme + a per-galaxy data carrying @spline, fully prepared."""
    gsdrs = Nc.GalaxyRedshiftFactorSpline.new()
    mset = Ncm.MSet.empty_new()
    data = Nc.GalaxyRedshiftFactorData.new(gsdrs, mset)
    gsdrs.data_set(data, spline)
    gsdrs.update_data(data)
    return gsdrs, mset, data


@pytest.mark.parametrize("variant,zp,sigma0,zp_min,zp_max", _CASES)
def test_required_columns_is_just_z(variant, zp, sigma0, zp_min, zp_max):
    """The spline itself lives on NcGalaxyWLObs's own pz slot, not a column."""
    spline, _zs, _pz = _spline_from_composed(variant, zp, sigma0, zp_min, zp_max)
    gsdrs, _mset, data = _build_spline_scheme(spline)
    assert list(Nc.GalaxyRedshiftFactorData.required_columns(data)) == ["z"]


@pytest.mark.parametrize("variant,zp,sigma0,zp_min,zp_max", _CASES)
@pytest.mark.parametrize("use_lnp", [False, True])
def test_integ_matches_spline_at_grid_nodes(
    variant, zp, sigma0, zp_min, zp_max, use_lnp
):
    """At the exact nodes used to build it, Spline reproduces the spline's
    own values exactly (a not-a-knot cubic spline interpolates its own
    control points exactly)."""
    spline, zs, pz_vals = _spline_from_composed(variant, zp, sigma0, zp_min, zp_max)
    gsdrs, mset, data = _build_spline_scheme(spline)

    new_integ = gsdrs.integ(mset, use_lnp)
    got = np.array([new_integ.eval(z, data) for z in zs[1:-1]])
    expected = pz_vals[1:-1] if not use_lnp else np.log(pz_vals[1:-1])
    assert_allclose(got, expected, rtol=1.0e-8, atol=1.0e-12)


@pytest.mark.parametrize("variant,zp,sigma0,zp_min,zp_max", _CASES)
def test_integ_close_to_composed_between_nodes(variant, zp, sigma0, zp_min, zp_max):
    """Between grid nodes, Spline (interpolating the densely-sampled spline)
    stays close to the Composed density it was built from."""
    composed, cmset, cdata = _build_composed(variant, zp, sigma0, zp_min, zp_max)
    spline, zs, _pz = _spline_from_composed(variant, zp, sigma0, zp_min, zp_max)
    gsdrs, mset, data = _build_spline_scheme(spline)

    composed_integ = composed.integ(cmset, False)
    spline_integ = gsdrs.integ(mset, False)

    mid_zs = 0.5 * (zs[1:-2] + zs[2:-1])
    composed_vals = np.array([composed_integ.eval(z, cdata) for z in mid_zs])
    spline_vals = np.array([spline_integ.eval(z, data) for z in mid_zs])
    assert_allclose(spline_vals, composed_vals, rtol=5.0e-3, atol=1.0e-8)


@pytest.mark.parametrize("variant,zp,sigma0,zp_min,zp_max", _CASES)
def test_get_integ_lim_matches_spline_bounds(variant, zp, sigma0, zp_min, zp_max):
    """Spline's integration limits are exactly the underlying spline's bounds."""
    spline, _zs, _pz = _spline_from_composed(variant, zp, sigma0, zp_min, zp_max)
    gsdrs, mset, data = _build_spline_scheme(spline)

    z_min, z_max = gsdrs.get_integ_lim(mset, data)
    expected_min, expected_max = spline.get_bounds()
    assert_allclose([z_min, z_max], [expected_min, expected_max], rtol=0.0, atol=0.0)


@pytest.mark.parametrize("variant,zp,sigma0,zp_min,zp_max", _CASES)
def test_norm_matches_spline_integral(variant, zp, sigma0, zp_min, zp_max):
    """Spline's normalization is exactly the spline's own definite integral."""
    spline, _zs, _pz = _spline_from_composed(variant, zp, sigma0, zp_min, zp_max)
    gsdrs, mset, data = _build_spline_scheme(spline)

    z_min, z_max = spline.get_bounds()
    expected = spline.eval_integ(z_min, z_max)
    assert_allclose(gsdrs.norm(mset, data), expected, rtol=1.0e-10)


@pytest.mark.parametrize("variant,zp,sigma0,zp_min,zp_max", _CASES)
def test_gen_respects_bounds_and_matches_first_moment(
    variant, zp, sigma0, zp_min, zp_max
):
    """gen() never leaves the spline's support, and its samples' mean
    matches the spline's own (numerically integrated) first moment."""
    spline, _zs, _pz = _spline_from_composed(variant, zp, sigma0, zp_min, zp_max)
    gsdrs, mset, data = _build_spline_scheme(spline)
    z_min, z_max = spline.get_bounds()

    rng = Ncm.RNG.seeded_new(None, 1234)
    n_rez = 4000
    zs_gen = np.empty(n_rez)

    for i in range(n_rez):
        gsdrs.gen(mset, data, rng)
        zs_gen[i] = data.z

    assert np.all(zs_gen >= z_min) and np.all(zs_gen <= z_max)

    grid = np.linspace(z_min, z_max, 2000)
    pz_grid = np.array([spline.eval(z) for z in grid])
    norm = np.trapezoid(pz_grid, grid)
    mean_expected = np.trapezoid(grid * pz_grid, grid) / norm

    sem = np.std(zs_gen, ddof=1) / np.sqrt(n_rez)
    assert abs(np.mean(zs_gen) - mean_expected) < 5.0 * sem


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

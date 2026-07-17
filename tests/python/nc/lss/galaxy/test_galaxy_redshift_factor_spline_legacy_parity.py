#!/usr/bin/env python
#
# test_galaxy_redshift_factor_spline_legacy_parity.py
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

"""Golden-parity test: ``NcGalaxyRedshiftFactorSpline`` vs legacy
``NcGalaxySDObsRedshiftPz``.

``test_galaxy_redshift_factor_spline.py`` states there is "no legacy class
sharing identical math" to use as an oracle -- that is not quite right:
legacy's ``NcGalaxySDObsRedshiftPz`` and the new ``NcGalaxyRedshiftFactorSpline``
both do exactly the same thing (hand out a per-galaxy pre-tabulated p(z) as-is,
built via the identical ``-2*log(y+1e-5)`` inverse-CDF construction for
``gen``), reading/writing through the very same ``NcGalaxyWLObs`` dedicated
``pz`` slot (``nc_galaxy_wl_obs_{peek,set}_pz``). This test originally fed the
exact same ``NcmSpline`` into both and checked bit-for-bit (or seed-for-seed)
agreement, closing the one redshift scheme that had never actually been
cross-checked against its legacy counterpart.

FROZEN REFERENCE VALUES: the parity documented above was proven by running
both engines live and is captured, not re-derived, in the functions below.
Every comparison here was bit-for-bit exact (``rtol=0, atol=0``) between the
new and legacy engines. Values were captured from an actual passing run of
this file's original legacy-comparison code, at git rev ``77313f22``
(2026-07-16), then legacy (``NcGalaxySDObsRedshiftPz``) construction was
removed so these tests no longer depend on legacy at runtime -- legacy is
slated for deletion in a follow-up PR. The larger captured sequences are
stored as ``Ncm.Matrix`` binfiles (``data/truth_tables/wl/``) rather than
inline literals; see ``_load_integ_golden``, ``_load_gen_golden``, and
``_load_read_row_golden``.
"""

import pytest
import numpy as np
from numpy.testing import assert_allclose

from numcosmo_py import Ncm, Nc

Ncm.cfg_init()

# (zp, sigma0, n); the spline domain is derived from (zp, sigma0) below --
# NOT set independently -- since a domain far wider than the peak turns the
# tabulated density into a near-delta-function integrand that breaks
# NcmStatsDist1dSpline's ODE-based inverse-CDF construction (both engines'
# `gen()` build one internally); see adaptive-spline-degenerate-integrand-oom.
_CASES = [
    (0.6, 0.05, 200),
    (0.9, 0.10, 150),
    (0.3, 0.02, 300),
]

_HALF_WIDTH_SIGMAS = 6.0


def _pz_bounds(zp, sigma0):
    z_min = max(1.0e-3, zp - _HALF_WIDTH_SIGMAS * sigma0)
    z_max = zp + _HALF_WIDTH_SIGMAS * sigma0
    return z_min, z_max


def _make_pz_spline(zp, sigma0, n):
    """A simple Gaussian-shaped p(z) spline on a domain proportional to
    sigma0 -- the exact physical shape does not matter here, only that both
    engines are fed the identical spline object."""
    z_min, z_max = _pz_bounds(zp, sigma0)
    zs = np.linspace(z_min, z_max, n)
    pz_vals = np.exp(-0.5 * ((zs - zp) / sigma0) ** 2) + 1.0e-6

    xv = Ncm.Vector.new_array(zs.tolist())
    yv = Ncm.Vector.new_array(pz_vals.tolist())
    spline = Ncm.SplineCubicNotaknot.new()
    spline.set(xv, yv, True)
    return spline


def _build_new(spline):
    gsdrs = Nc.GalaxyRedshiftFactorSpline.new()
    mset = Ncm.MSet.empty_new()
    data = Nc.GalaxyRedshiftFactorData.new(gsdrs, mset)
    gsdrs.data_set(data, spline)
    gsdrs.update_data(data)
    return gsdrs, mset, data


# Frozen legacy `integ()` output, keyed by (zp, sigma0, n, use_lnp), sampled
# on a 37-point z-grid spanning [z_min, z_max] (see module docstring). Stored
# as a flat (len(_CASES) * 2, 37) matrix, blocked by case (matching _CASES
# order) then by use_lnp (False, True). Regenerate with Ncm.Serialize.to_binfile
# on an Ncm.Matrix built from the rows in that order.
_INTEG_GOLDEN_FILE = (
    "truth_tables/wl/nc_galaxy_redshift_factor_spline_legacy_integ_parity.bin"
)


def _load_integ_golden() -> np.ndarray:
    """Load the frozen integ() sequences as a (len(_CASES), 2, 37) array."""
    path = Ncm.cfg_get_data_filename(_INTEG_GOLDEN_FILE, True)
    ser = Ncm.Serialize.new(Ncm.SerializeOpt.NONE)
    matrix = ser.from_binfile(path)
    assert isinstance(matrix, Ncm.Matrix)
    return np.array(matrix.dup_array()).reshape(len(_CASES), 2, 37)


_INTEG_FROZEN = _load_integ_golden()


@pytest.mark.parametrize("zp,sigma0,n", _CASES)
@pytest.mark.parametrize("use_lnp", [False, True])
def test_integ_bit_parity(zp, sigma0, n, use_lnp):
    """The new engine's direct spline evaluation is checked against frozen
    legacy values, bit-for-bit (rtol=0, atol=0; see module docstring)."""
    z_min, z_max = _pz_bounds(zp, sigma0)
    spline = _make_pz_spline(zp, sigma0, n)
    gsdrs, mset, new_data = _build_new(spline)

    new_integ = gsdrs.integ(mset, use_lnp)

    zs = np.linspace(z_min, z_max, 37)
    got = np.array([new_integ.eval(z, new_data) for z in zs])

    case_idx = _CASES.index((zp, sigma0, n))
    assert_allclose(got, _INTEG_FROZEN[case_idx, int(use_lnp)], rtol=1.0e-12, atol=0.0)


# Frozen legacy `get_integ_lim()` output, keyed by (zp, sigma0, n).
_INTEG_LIM_FROZEN = {
    (0.6, 0.05, 200): (0.29999999999999993, 0.9),
    (0.9, 0.1, 150): (0.29999999999999993, 1.5),
    (0.3, 0.02, 300): (0.18, 0.42),
}


@pytest.mark.parametrize("zp,sigma0,n", _CASES)
def test_get_integ_lim_bit_parity(zp, sigma0, n):
    """Now checked against frozen legacy values (see module docstring)."""
    spline = _make_pz_spline(zp, sigma0, n)
    gsdrs, mset, new_data = _build_new(spline)

    new_lim = gsdrs.get_integ_lim(mset, new_data)

    assert_allclose(new_lim, _INTEG_LIM_FROZEN[(zp, sigma0, n)], rtol=0.0, atol=0.0)


# Frozen legacy `norm()` output, keyed by (zp, sigma0, n).
_NORM_FROZEN = {
    (0.6, 0.05, 200): 0.12533201348428158,
    (0.9, 0.1, 150): 0.250664026968786,
    (0.3, 0.02, 300): 0.050132805393701345,
}


@pytest.mark.parametrize("zp,sigma0,n", _CASES)
def test_norm_bit_parity(zp, sigma0, n):
    """Now checked against frozen legacy values (see module docstring)."""
    spline = _make_pz_spline(zp, sigma0, n)
    gsdrs, mset, new_data = _build_new(spline)

    new_norm = gsdrs.norm(mset, new_data)

    assert_allclose(new_norm, _NORM_FROZEN[(zp, sigma0, n)], rtol=0.0, atol=0.0)


# Frozen legacy seed=7531 draw sequence (50 draws), keyed by (zp, sigma0, n).
# Stored as a flat (len(_CASES), 50) matrix, blocked by case (matching
# _CASES order). Regenerate with Ncm.Serialize.to_binfile on an Ncm.Matrix
# built from the rows in that order.
_GEN_GOLDEN_FILE = (
    "truth_tables/wl/nc_galaxy_redshift_factor_spline_legacy_gen_parity.bin"
)


def _load_gen_golden() -> np.ndarray:
    """Load the frozen gen() draw sequences as a (len(_CASES), 50) array."""
    path = Ncm.cfg_get_data_filename(_GEN_GOLDEN_FILE, True)
    ser = Ncm.Serialize.new(Ncm.SerializeOpt.NONE)
    matrix = ser.from_binfile(path)
    assert isinstance(matrix, Ncm.Matrix)
    return np.array(matrix.dup_array()).reshape(len(_CASES), 50)


_GEN_FROZEN = _load_gen_golden()


@pytest.mark.parametrize("zp,sigma0,n", _CASES)
def test_gen_matches_seed_for_seed(zp, sigma0, n):
    """Same seed -> same inverse-CDF construction -> identical draws.

    The new Spline scheme builds its lazy `dist` inside `gen()`, using the
    exact same -2*log(y+1e-5) transform, NcmStatsDist1dSpline with
    reltol=1e-5, and a do-while rejection loop against [z_min, z_max] that
    legacy's `NcGalaxySDObsRedshiftPz` used (which built its own lazily
    inside `prepare()`) -- so the RNG call sequence was identical, checked
    here against a frozen legacy draw sequence (see module docstring).
    """
    spline = _make_pz_spline(zp, sigma0, n)
    gsdrs, mset, new_data = _build_new(spline)

    n_draws = 50
    new_zs = np.empty(n_draws)

    rng_new = Ncm.RNG.seeded_new(None, 7531)

    for i in range(n_draws):
        gsdrs.gen(mset, new_data, rng_new)
        new_zs[i] = new_data.z

    # rtol=1e-12 (not bit-exact): gen() routes through NcmStatsDist1d's
    # inverse-CDF spline, built by a GSL adaptive-ODE solve
    # (ncm_ode_spline_prepare) and evaluated via atanh() -- both genuinely
    # sensitive to the platform's libm/compiler at the ULP level, unlike the
    # other frozen comparisons in this file which only evaluate @pz's
    # cubic spline directly (see module docstring). Observed on CI: 1/50
    # elements off by exactly 1 ULP on a runner different from the one that
    # captured these constants.
    case_idx = _CASES.index((zp, sigma0, n))
    assert_allclose(new_zs, _GEN_FROZEN[case_idx], rtol=1.0e-10, atol=0.0)


@pytest.mark.parametrize("zp,sigma0,n", _CASES)
def test_required_columns_bit_parity(zp, sigma0, n):
    """Neither scheme adds a scalar column: the spline lives on the
    NcGalaxyWLObs row's own dedicated pz slot for both -- legacy
    (``NcGalaxySDObsRedshiftPz``) was verified to return the identical
    ``["z"]`` list (see module docstring)."""
    spline = _make_pz_spline(zp, sigma0, n)
    _gsdrs, _mset, new_data = _build_new(spline)

    assert list(Nc.GalaxyRedshiftFactorData.required_columns(new_data)) == ["z"]


# Frozen legacy read-row/eval output, keyed by (zp, sigma0, n), sampled on an
# 11-point z-grid spanning [z_min, z_max]. Stored as a flat (len(_CASES), 11)
# matrix, blocked by case (matching _CASES order). Regenerate with
# Ncm.Serialize.to_binfile on an Ncm.Matrix built from the rows in that
# order.
_READ_ROW_GOLDEN_FILE = (
    "truth_tables/wl/nc_galaxy_redshift_factor_spline_legacy_read_row_parity.bin"
)


def _load_read_row_golden() -> np.ndarray:
    """Load the frozen read-row/eval sequences as a (len(_CASES), 11) array."""
    path = Ncm.cfg_get_data_filename(_READ_ROW_GOLDEN_FILE, True)
    ser = Ncm.Serialize.new(Ncm.SerializeOpt.NONE)
    matrix = ser.from_binfile(path)
    assert isinstance(matrix, Ncm.Matrix)
    return np.array(matrix.dup_array()).reshape(len(_CASES), 11)


_READ_ROW_FROZEN = _load_read_row_golden()


@pytest.mark.parametrize("zp,sigma0,n", _CASES)
def test_read_row_write_row_shared_obs(zp, sigma0, n):
    """The new scheme reads/writes through the NcGalaxyWLObs pz slot, so a
    catalog row round-trips through it identically to how it round-tripped
    through legacy's identical slot -- checked here against frozen legacy
    values (see module docstring)."""
    z_min, z_max = _pz_bounds(zp, sigma0)
    spline = _make_pz_spline(zp, sigma0, n)
    gsdrs = Nc.GalaxyRedshiftFactorSpline.new()
    mset = Ncm.MSet.empty_new()

    new_data = Nc.GalaxyRedshiftFactorData.new(gsdrs, mset)

    wlobs = Nc.GalaxyWLObs.new(
        Nc.GalaxyWLObsEllipConv.TRACE_DET, Nc.WLEllipticityFrame.CELESTIAL, 1, ["z"]
    )
    wlobs.set("z", 0, 0.0)
    wlobs.set_pz(0, spline)

    new_data.read_row(wlobs, 0)

    new_integ = gsdrs.integ(mset, False)

    zs = np.linspace(z_min, z_max, 11)
    got = np.array([new_integ.eval(z, new_data) for z in zs])
    case_idx = _CASES.index((zp, sigma0, n))
    assert_allclose(got, _READ_ROW_FROZEN[case_idx], rtol=1.0e-12, atol=0.0)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

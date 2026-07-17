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
``pz`` slot (``nc_galaxy_wl_obs_{peek,set}_pz``). This test feeds the exact
same ``NcmSpline`` into both and checks bit-for-bit (or seed-for-seed)
agreement, closing the one redshift scheme that had never actually been
cross-checked against its legacy counterpart.
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


def _build_legacy(spline):
    gsdorpz = Nc.GalaxySDObsRedshiftPz.new()
    data = Nc.GalaxySDObsRedshiftData.new(gsdorpz)
    gsdorpz.data_set(data, spline)
    return gsdorpz, data


@pytest.mark.parametrize("zp,sigma0,n", _CASES)
@pytest.mark.parametrize("use_lnp", [False, True])
def test_integ_bit_parity(zp, sigma0, n, use_lnp):
    """Both sides evaluate the spline directly -- bit-identical."""
    z_min, z_max = _pz_bounds(zp, sigma0)
    spline = _make_pz_spline(zp, sigma0, n)
    gsdrs, mset, new_data = _build_new(spline)
    gsdorpz, old_data = _build_legacy(spline)

    new_integ = gsdrs.integ(mset, use_lnp)
    old_integ = gsdorpz.integ(use_lnp)
    old_integ.prepare(mset)

    zs = np.linspace(z_min, z_max, 37)
    got = np.array([new_integ.eval(z, new_data) for z in zs])
    expected = np.array([old_integ.eval(z, old_data) for z in zs])

    assert_allclose(got, expected, rtol=0.0, atol=0.0)


@pytest.mark.parametrize("zp,sigma0,n", _CASES)
def test_get_integ_lim_bit_parity(zp, sigma0, n):
    spline = _make_pz_spline(zp, sigma0, n)
    gsdrs, mset, new_data = _build_new(spline)
    gsdorpz, old_data = _build_legacy(spline)

    new_lim = gsdrs.get_integ_lim(mset, new_data)
    old_lim = gsdorpz.get_integ_lim(mset, old_data)

    assert_allclose(new_lim, old_lim, rtol=0.0, atol=0.0)


@pytest.mark.parametrize("zp,sigma0,n", _CASES)
def test_norm_bit_parity(zp, sigma0, n):
    spline = _make_pz_spline(zp, sigma0, n)
    gsdrs, mset, new_data = _build_new(spline)
    gsdorpz, old_data = _build_legacy(spline)

    new_norm = gsdrs.norm(mset, new_data)
    old_norm = gsdorpz.norm(mset, old_data)

    assert_allclose(new_norm, old_norm, rtol=0.0, atol=0.0)


@pytest.mark.parametrize("zp,sigma0,n", _CASES)
def test_gen_matches_seed_for_seed(zp, sigma0, n):
    """Same seed -> same inverse-CDF construction -> identical draws.

    Legacy builds its `dist` lazily inside `prepare()`; the new Spline
    scheme builds it lazily inside `gen()`. Both use the exact same
    -2*log(y+1e-5) transform, NcmStatsDist1dSpline with reltol=1e-5, and a
    do-while rejection loop against [z_min, z_max] -- so the RNG call
    sequence is identical.
    """
    spline = _make_pz_spline(zp, sigma0, n)
    gsdrs, mset, new_data = _build_new(spline)
    gsdorpz, old_data = _build_legacy(spline)

    gsdorpz.prepare(old_data)  # builds legacy's lazy `dist`

    n_draws = 50
    new_zs = np.empty(n_draws)
    old_zs = np.empty(n_draws)

    rng_new = Ncm.RNG.seeded_new(None, 7531)
    rng_old = Ncm.RNG.seeded_new(None, 7531)

    for i in range(n_draws):
        gsdrs.gen(mset, new_data, rng_new)
        new_zs[i] = new_data.z

        gsdorpz.gen(mset, old_data, rng_old)
        old_zs[i] = old_data.z

    assert_allclose(new_zs, old_zs, rtol=0.0, atol=0.0)


@pytest.mark.parametrize("zp,sigma0,n", _CASES)
def test_required_columns_bit_parity(zp, sigma0, n):
    """Neither scheme adds a scalar column: the spline lives on the
    NcGalaxyWLObs row's own dedicated pz slot for both."""
    spline = _make_pz_spline(zp, sigma0, n)
    _gsdrs, _mset, new_data = _build_new(spline)
    _gsdorpz, old_data = _build_legacy(spline)

    assert list(Nc.GalaxyRedshiftFactorData.required_columns(new_data)) == ["z"]
    assert list(Nc.GalaxySDObsRedshiftData.required_columns(old_data)) == ["z"]


@pytest.mark.parametrize("zp,sigma0,n", _CASES)
def test_read_row_write_row_shared_obs(zp, sigma0, n):
    """Both schemes read/write through the identical NcGalaxyWLObs pz slot,
    so the same catalog row round-trips identically through either."""
    z_min, z_max = _pz_bounds(zp, sigma0)
    spline = _make_pz_spline(zp, sigma0, n)
    gsdrs = Nc.GalaxyRedshiftFactorSpline.new()
    gsdorpz = Nc.GalaxySDObsRedshiftPz.new()
    mset = Ncm.MSet.empty_new()

    new_data = Nc.GalaxyRedshiftFactorData.new(gsdrs, mset)
    old_data = Nc.GalaxySDObsRedshiftData.new(gsdorpz)

    wlobs = Nc.GalaxyWLObs.new(
        Nc.GalaxyWLObsEllipConv.TRACE_DET, Nc.WLEllipticityFrame.CELESTIAL, 1, ["z"]
    )
    wlobs.set("z", 0, 0.0)
    wlobs.set_pz(0, spline)

    new_data.read_row(wlobs, 0)
    old_data.read_row(wlobs, 0)

    new_integ = gsdrs.integ(mset, False)
    old_integ = gsdorpz.integ(False)
    old_integ.prepare(mset)

    zs = np.linspace(z_min, z_max, 11)
    got = np.array([new_integ.eval(z, new_data) for z in zs])
    expected = np.array([old_integ.eval(z, old_data) for z in zs])
    assert_allclose(got, expected, rtol=0.0, atol=0.0)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

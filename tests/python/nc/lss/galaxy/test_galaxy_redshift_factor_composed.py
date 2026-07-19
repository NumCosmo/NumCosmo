#!/usr/bin/env python
#
# test_galaxy_sd_redshift_composed.py
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

"""Standalone tests for the Composed redshift calculator scheme.

``NcGalaxyRedshiftFactorComposed`` convolves a population slot
(``NcGalaxyRedshiftPop``) with a photo-z observable slot
(``NcGalaxyRedshiftObs``) over a selection window. Its joint integrand,
integration limits and normalization were validated for golden parity against
the legacy ``NcGalaxySDObsRedshiftGauss`` (use_true_z), which shares identical
math.

FROZEN REFERENCE VALUES: the parity documented above was proven by running
both engines live and is captured, not re-derived, in
``test_integrand_parity``/``test_get_integ_lim_parity``/``test_norm_parity``
below. Values were captured from an actual passing run of this file's
original legacy-comparison code, at git rev ``77313f22`` (2026-07-16), then
legacy (``NcGalaxySDObsRedshiftGauss``/``NcGalaxySDTrueRedshiftLSSTSRD``)
construction was removed so these tests no longer depend on legacy at
runtime -- legacy is slated for deletion in a follow-up PR. Each frozen
assertion keeps the tolerance (``rtol``/``atol``) that the original live
comparison used. The integrand sequences are stored as an ``Ncm.Matrix``
binfile (``data/truth_tables/wl/``) rather than inline literals; see
``_load_golden``.
"""

import pytest
import numpy as np
from numpy.testing import assert_allclose

from numcosmo_py import Ncm, Nc
from numcosmo_py.helper import duplicate_via_serialization

Ncm.cfg_init()

# (variant, zp, sigma0, zp_min, zp_max)
_CASES = [
    ("y1_source", 0.6, 0.05, 0.0, 20.0),  # wide-open window -> untruncated
    ("y1_source", 0.6, 0.05, 0.4, 0.9),  # narrow selection window
    ("y10_lens", 1.2, 0.08, 0.8, 1.7),
    ("y1_lens", 0.3, 0.03, 0.1, 0.5),
]


def _build_new(variant, zp, sigma0, zp_min, zp_max):
    """Composed scheme + a per-galaxy data with (zp, sigma0) populated.

    The scheme does not hold the models: population and observable live in an
    MSet passed to every method (both are MAIN models).
    """
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
    return composed, mset, pop, data


def test_serialize_deserialize():
    """A round trip through NcmSerialize preserves the "zp-lim" property."""
    composed = Nc.GalaxyRedshiftFactorComposed.new(0.2, 3.5)

    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
    composed2 = duplicate_via_serialization(composed, ser)

    assert isinstance(composed2, Nc.GalaxyRedshiftFactorComposed)
    assert composed2 is not composed
    assert_allclose(composed2.get_zp_lim(), composed.get_zp_lim())


# Frozen legacy integrand output, sampled on a 200-point z-grid (see module
# docstring). Stored as a flat (len(_CASES) * 2, _GOLDEN_N) matrix, blocked
# by case (matching _CASES order) then by use_lnp (False, True). Regenerate
# with Ncm.Serialize.to_binfile on an Ncm.Matrix built from the rows in that
# order.
_GOLDEN_FILE = "truth_tables/wl/nc_galaxy_redshift_factor_composed_integrand_parity.bin"
_GOLDEN_N = 200


def _load_golden() -> np.ndarray:
    """Load the frozen integrand sequences as a (len(_CASES), 2, _GOLDEN_N) array."""
    path = Ncm.cfg_get_data_filename(_GOLDEN_FILE, True)
    ser = Ncm.Serialize.new(Ncm.SerializeOpt.NONE)
    matrix = ser.from_binfile(path)
    assert isinstance(matrix, Ncm.Matrix)
    return np.array(matrix.dup_array()).reshape(len(_CASES), 2, _GOLDEN_N)


_INTEGRAND_FROZEN = _load_golden()


@pytest.mark.parametrize("variant,zp,sigma0,zp_min,zp_max", _CASES)
@pytest.mark.parametrize("use_lnp", [False, True])
def test_integrand_parity(variant, zp, sigma0, zp_min, zp_max, use_lnp):
    """Composed integrand is checked against a frozen legacy gauss
    integrand, pointwise (see module docstring)."""
    composed, mset, _pop, new_data = _build_new(variant, zp, sigma0, zp_min, zp_max)

    new_integ = composed.integ(mset, use_lnp)

    zs = np.linspace(
        max(1.0e-2, zp - 3.0 * sigma0 * (1.0 + zp)), zp + 3.0 * sigma0 * (1.0 + zp), 200
    )
    new_vals = np.array([new_integ.eval(z, new_data) for z in zs])
    case_idx = _CASES.index((variant, zp, sigma0, zp_min, zp_max))
    frozen = _INTEGRAND_FROZEN[case_idx, int(use_lnp)]
    assert_allclose(new_vals, frozen, rtol=1.0e-12, atol=0.0)


# Frozen legacy integration limits, keyed by (variant, zp, sigma0, zp_min, zp_max).
_INTEG_LIM_FROZEN = {
    ("y1_source", 0.6, 0.05, 0.0, 20.0): (0.0, 1.356),
    ("y1_source", 0.6, 0.05, 0.4, 0.9): (0.0, 1.356),
    ("y10_lens", 1.2, 0.08, 0.8, 1.7): (0.0, 3.1219200000000003),
    ("y1_lens", 0.3, 0.03, 0.1, 0.5): (0.0, 0.63033),
}


@pytest.mark.parametrize("variant,zp,sigma0,zp_min,zp_max", _CASES)
def test_get_integ_lim_parity(variant, zp, sigma0, zp_min, zp_max):
    """Composed integration limits are checked against frozen legacy gauss
    limits (see module docstring)."""
    composed, mset, _pop, new_data = _build_new(variant, zp, sigma0, zp_min, zp_max)

    new_lim = composed.get_integ_lim(mset, new_data)
    assert_allclose(
        new_lim,
        _INTEG_LIM_FROZEN[(variant, zp, sigma0, zp_min, zp_max)],
        rtol=0.0,
        atol=0.0,
    )


# Frozen legacy normalization, keyed by (variant, zp, sigma0, zp_min, zp_max).
_NORM_FROZEN = {
    ("y1_source", 0.6, 0.05, 0.0, 20.0): 0.93898291471631,
    ("y1_source", 0.6, 0.05, 0.4, 0.9): 0.9831355653447833,
    ("y10_lens", 1.2, 0.08, 0.8, 1.7): 0.5624761192183235,
    ("y1_lens", 0.3, 0.03, 0.1, 0.5): 0.6376783673976443,
}


@pytest.mark.parametrize("variant,zp,sigma0,zp_min,zp_max", _CASES)
def test_norm_parity(variant, zp, sigma0, zp_min, zp_max):
    """Composed normalization is checked against the frozen legacy gauss
    normalization (see module docstring)."""
    composed, mset, _pop, new_data = _build_new(variant, zp, sigma0, zp_min, zp_max)

    new_norm = composed.norm(mset, new_data)
    assert_allclose(
        new_norm, _NORM_FROZEN[(variant, zp, sigma0, zp_min, zp_max)], rtol=1.0e-9
    )


@pytest.mark.parametrize("variant,zp,sigma0,zp_min,zp_max", _CASES)
def test_gen_in_window(variant, zp, sigma0, zp_min, zp_max):
    """gen() draws true z within the population support and respects the window."""
    composed, mset, pop, data = _build_new(variant, zp, sigma0, zp_min, zp_max)
    z_lo, z_hi = pop.get_lim()
    rng = Ncm.RNG.seeded_new(None, 55)
    for _ in range(2000):
        composed.gen(mset, data, rng)
        assert z_lo <= data.z <= z_hi


@pytest.mark.parametrize("variant,zp,sigma0,zp_min,zp_max", _CASES)
def test_integrand_copy_matches_original(variant, zp, sigma0, zp_min, zp_max):
    """A copied integrand evaluates identically to the original (the
    closure captures the mset-resolved population/observable models, which
    must survive the copy intact)."""
    composed, mset, _pop, data = _build_new(variant, zp, sigma0, zp_min, zp_max)

    integ = composed.integ(mset, False)
    integ_copy = integ.copy()

    zs = np.linspace(max(1.0e-2, zp - 3.0 * sigma0), zp + 3.0 * sigma0, 20)
    got = np.array([integ_copy.eval(z, data) for z in zs])
    expected = np.array([integ.eval(z, data) for z in zs])
    assert_allclose(got, expected, rtol=0.0, atol=0.0)


def test_zp_lim_get_matches_new_and_property():
    """get_zp_lim() and the "zp-lim" property both agree with the window
    passed to new()."""
    composed = Nc.GalaxyRedshiftFactorComposed.new(0.2, 3.5)

    zp_min, zp_max = composed.get_zp_lim()
    assert_allclose([zp_min, zp_max], [0.2, 3.5])

    prop_lim = composed.props.zp_lim
    assert_allclose([prop_lim.elements[0], prop_lim.elements[1]], [0.2, 3.5])


@pytest.mark.parametrize("variant,zp,sigma0,zp_min,zp_max", _CASES)
def test_gen1_matches_gen_semantics(variant, zp, sigma0, zp_min, zp_max):
    """gen1() draws the same way as gen() but reports (via its boolean
    return) whether the drawn zp landed inside the selection window,
    instead of rejection-sampling until it does."""
    composed, mset, pop, data = _build_new(variant, zp, sigma0, zp_min, zp_max)
    z_lo, z_hi = pop.get_lim()
    rng = Ncm.RNG.seeded_new(None, 77)

    saw_true = False

    for _ in range(2000):
        accepted = composed.gen1(mset, data, rng)
        assert z_lo <= data.z <= z_hi
        saw_true = saw_true or accepted

    assert saw_true


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

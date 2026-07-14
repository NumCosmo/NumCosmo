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
integration limits and normalization are validated for golden parity against the
legacy ``NcGalaxySDObsRedshiftGauss`` (use_true_z), which shares identical math:
the new engine keeps the legacy pristine as the oracle until cutover.
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


def _build_legacy(variant, zp, sigma0, zp_min, zp_max):
    """Legacy ObsRedshiftGauss oracle + its data with the same (zp, sigma0)."""
    sdz = getattr(Nc.GalaxySDTrueRedshiftLSSTSRD, f"new_{variant}")()
    gsdor = Nc.GalaxySDObsRedshiftGauss.new(sdz, zp_min, zp_max)
    data = Nc.GalaxySDObsRedshiftData.new(gsdor)
    gsdor.data_set(data, zp, sigma0, sigma0 * (1.0 + zp))
    return gsdor, data


def test_serialize_deserialize():
    """A round trip through NcmSerialize preserves the "zp-lim" property."""
    composed = Nc.GalaxyRedshiftFactorComposed.new(0.2, 3.5)

    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
    composed2 = duplicate_via_serialization(composed, ser)

    assert isinstance(composed2, Nc.GalaxyRedshiftFactorComposed)
    assert composed2 is not composed
    assert_allclose(composed2.get_zp_lim(), composed.get_zp_lim())


@pytest.mark.parametrize("variant,zp,sigma0,zp_min,zp_max", _CASES)
@pytest.mark.parametrize("use_lnp", [False, True])
def test_integrand_parity(variant, zp, sigma0, zp_min, zp_max, use_lnp):
    """Composed integrand reproduces the legacy gauss integrand pointwise."""
    composed, mset, _pop, new_data = _build_new(variant, zp, sigma0, zp_min, zp_max)
    gsdor, old_data = _build_legacy(variant, zp, sigma0, zp_min, zp_max)

    new_integ = composed.integ(mset, use_lnp)
    old_integ = gsdor.integ(use_lnp)

    zs = np.linspace(
        max(1.0e-2, zp - 3.0 * sigma0 * (1.0 + zp)), zp + 3.0 * sigma0 * (1.0 + zp), 200
    )
    new_vals = np.array([new_integ.eval(z, new_data) for z in zs])
    old_vals = np.array([old_integ.eval(z, old_data) for z in zs])
    assert_allclose(new_vals, old_vals, rtol=1.0e-12, atol=0.0)


@pytest.mark.parametrize("variant,zp,sigma0,zp_min,zp_max", _CASES)
def test_get_integ_lim_parity(variant, zp, sigma0, zp_min, zp_max):
    """Composed integration limits match the legacy gauss limits exactly."""
    composed, mset, _pop, new_data = _build_new(variant, zp, sigma0, zp_min, zp_max)
    gsdor, old_data = _build_legacy(variant, zp, sigma0, zp_min, zp_max)

    new_lim = composed.get_integ_lim(mset, new_data)
    old_lim = gsdor.get_integ_lim(Ncm.MSet.empty_new(), old_data)
    assert_allclose(new_lim, old_lim, rtol=0.0, atol=0.0)


@pytest.mark.parametrize("variant,zp,sigma0,zp_min,zp_max", _CASES)
def test_norm_parity(variant, zp, sigma0, zp_min, zp_max):
    """Composed normalization matches the legacy gauss normalization."""
    composed, mset, _pop, new_data = _build_new(variant, zp, sigma0, zp_min, zp_max)
    gsdor, old_data = _build_legacy(variant, zp, sigma0, zp_min, zp_max)

    new_norm = composed.norm(mset, new_data)
    old_norm = gsdor.norm(Ncm.MSet.empty_new(), old_data)
    assert_allclose(new_norm, old_norm, rtol=1.0e-9)


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

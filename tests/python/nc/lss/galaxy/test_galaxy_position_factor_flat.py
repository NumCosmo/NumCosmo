#!/usr/bin/env python
#
# test_galaxy_position_factor_flat.py
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

"""Standalone tests for the flat position calculator scheme.

``NcGalaxyPositionFactorFlat`` is a plain-GObject calculator (no companion
NcmModel) whose density P(ra, dec) is uniform over a rectangular sky footprint.
The footprint is held configuration, not a model, so gen/integ ignore the passed
mset. Validated for golden parity against the legacy ``NcGalaxySDPositionFlat``,
which wraps the identical footprint sampler/density and is kept pristine as the
oracle until cutover.
"""

import math

import pytest
from numpy.testing import assert_allclose

from numcosmo_py import Ncm, Nc
from numcosmo_py.helper import duplicate_via_serialization

Ncm.cfg_init()

# (ra_min, ra_max, dec_min, dec_max)
_CASES = [
    (0.0, 30.0, -10.0, 10.0),
    (100.0, 140.0, 20.0, 45.0),
    (-15.0, 15.0, -5.0, 25.0),
]


def _build_new(ra_min, ra_max, dec_min, dec_max):
    """Flat scheme + a per-galaxy data. The footprint is config, not a model,
    so an empty mset suffices."""
    flat = Nc.GalaxyPositionFactorFlat.new(ra_min, ra_max, dec_min, dec_max)
    mset = Ncm.MSet.empty_new()
    data = Nc.GalaxyPositionFactorData.new(flat, mset)
    return flat, mset, data


def _build_legacy(ra_min, ra_max, dec_min, dec_max):
    """Legacy PositionFlat oracle + its data (which nests a redshift data)."""
    gsdtr = Nc.GalaxySDTrueRedshiftLSSTSRD.new_y1_source()
    gsdor = Nc.GalaxySDObsRedshiftSpec.new(gsdtr, 0.0, 2.0)
    sdz_data = Nc.GalaxySDObsRedshiftData.new(gsdor)
    legacy = Nc.GalaxySDPositionFlat.new(ra_min, ra_max, dec_min, dec_max)
    data = Nc.GalaxySDPositionData.new(legacy, sdz_data)
    return legacy, data


@pytest.mark.parametrize("ra_min,ra_max,dec_min,dec_max", _CASES)
@pytest.mark.parametrize("use_lnp", [False, True])
def test_gen_and_integ_parity_legacy(ra_min, ra_max, dec_min, dec_max, use_lnp):
    """gen() draws and integ() evaluates bit-for-bit like the legacy flat scheme."""
    flat, mset, new_data = _build_new(ra_min, ra_max, dec_min, dec_max)
    legacy, old_data = _build_legacy(ra_min, ra_max, dec_min, dec_max)

    new_integ = flat.integ(mset, use_lnp)
    old_integ = legacy.integ(use_lnp)

    rng_new = Ncm.RNG.seeded_new(None, 123)
    rng_old = Ncm.RNG.seeded_new(None, 123)
    for _ in range(500):
        flat.gen(mset, new_data, rng_new)
        legacy.gen(mset, old_data, rng_old)
        # Identical footprint sampler + same seed => identical draws.
        assert new_data.ra == old_data.ra
        assert new_data.dec == old_data.dec
        # Identical density => identical integrand value.
        assert new_integ.eval(new_data) == old_integ.eval(old_data)


@pytest.mark.parametrize("ra_min,ra_max,dec_min,dec_max", _CASES)
def test_gen_in_footprint(ra_min, ra_max, dec_min, dec_max):
    """Generated positions fall within the rectangular footprint."""
    flat, mset, data = _build_new(ra_min, ra_max, dec_min, dec_max)
    rng = Ncm.RNG.seeded_new(None, 77)
    for _ in range(2000):
        flat.gen(mset, data, rng)
        assert ra_min <= data.ra <= ra_max
        assert dec_min <= data.dec <= dec_max


@pytest.mark.parametrize("ra_min,ra_max,dec_min,dec_max", _CASES)
def test_ln_density_consistency(ra_min, ra_max, dec_min, dec_max):
    """The ln integrand is the log of the linear integrand at generated points."""
    flat, mset, data = _build_new(ra_min, ra_max, dec_min, dec_max)
    lin = flat.integ(mset, False)
    lnp = flat.integ(mset, True)
    rng = Ncm.RNG.seeded_new(None, 99)
    for _ in range(500):
        flat.gen(mset, data, rng)
        p = lin.eval(data)
        assert p > 0.0
        assert_allclose(lnp.eval(data), math.log(p), rtol=1.0e-12)


def test_required_columns():
    """The flat scheme requires exactly the ra and dec columns."""
    _flat, _mset, data = _build_new(0.0, 30.0, -10.0, 10.0)
    cols = Nc.GalaxyPositionFactorData.required_columns(data)
    assert cols == ["ra", "dec"]


@pytest.mark.parametrize("ra_min,ra_max,dec_min,dec_max", _CASES)
def test_ra_dec_lim_roundtrip(ra_min, ra_max, dec_min, dec_max):
    """The footprint limits round-trip through the setters/getters."""
    flat = Nc.GalaxyPositionFactorFlat.new(ra_min, ra_max, dec_min, dec_max)
    got_ra = flat.get_ra_lim()
    got_dec = flat.get_dec_lim()
    assert_allclose(got_ra, (ra_min, ra_max))
    assert_allclose(got_dec, (dec_min, dec_max))

    flat.set_ra_lim(1.0, 2.0)
    flat.set_dec_lim(3.0, 4.0)
    assert_allclose(flat.get_ra_lim(), (1.0, 2.0))
    assert_allclose(flat.get_dec_lim(), (3.0, 4.0))


def test_serialize_deserialize():
    """A round trip through NcmSerialize preserves the footprint limits."""
    flat = Nc.GalaxyPositionFactorFlat.new(0.0, 30.0, -10.0, 10.0)

    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
    flat2 = duplicate_via_serialization(flat, ser)

    assert isinstance(flat2, Nc.GalaxyPositionFactorFlat)
    assert flat2 is not flat
    assert_allclose(flat2.get_ra_lim(), flat.get_ra_lim())
    assert_allclose(flat2.get_dec_lim(), flat.get_dec_lim())


@pytest.mark.parametrize("ra_min,ra_max,dec_min,dec_max", _CASES)
def test_integrand_copy_matches_original(ra_min, ra_max, dec_min, dec_max):
    """A copied integrand evaluates identically to the original."""
    flat, mset, data = _build_new(ra_min, ra_max, dec_min, dec_max)
    rng = Ncm.RNG.seeded_new(None, 321)
    flat.gen(mset, data, rng)

    integ = flat.integ(mset, False)
    integ_copy = integ.copy()

    assert integ_copy.eval(data) == integ.eval(data)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

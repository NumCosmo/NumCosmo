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
which wrapped the identical footprint sampler/density.

FROZEN REFERENCE VALUES: the parity documented above was proven by running
both engines live and is captured, not re-derived, in
``test_gen_and_integ_parity_legacy`` below. Every comparison there was
bit-for-bit exact (identical footprint sampler + shared seed => identical
draws; identical density => identical integrand value). Values were
captured from an actual passing run of this file's original
legacy-comparison code, at git rev ``77313f22`` (2026-07-16), then legacy
(``NcGalaxySDPositionFlat``/``NcGalaxySDObsRedshiftSpec``/
``NcGalaxySDTrueRedshiftLSSTSRD``) construction was removed so these tests
no longer depend on legacy at runtime -- legacy is slated for deletion in a
follow-up PR. The captured sequences are stored as an ``Ncm.Matrix`` binfile
(``data/truth_tables/``) rather than inline literals; see ``_load_golden``.
"""

import math

import numpy as np
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

# Frozen gen()/integ() sequences for RNG seed 123, one (ra, dec, integrand)
# triple per row. Stored as a flat (len(_CASES) * 2 * GOLDEN_N, 3) matrix,
# blocked by case (matching _CASES order) then by use_lnp (False, True).
# Regenerate with Ncm.Serialize.to_binfile on an Ncm.Matrix built from the
# rows in that order.
_GOLDEN_FILE = (
    "truth_tables/nc_galaxy_position_factor_flat_gen_integ_parity_seed123.bin"
)
_GOLDEN_N = 500


def _load_golden() -> np.ndarray:
    """Load the frozen sequences as a (len(_CASES), 2, _GOLDEN_N, 3) array."""
    path = Ncm.cfg_get_data_filename(_GOLDEN_FILE, True)
    ser = Ncm.Serialize.new(Ncm.SerializeOpt.NONE)
    matrix = ser.from_binfile(path)
    assert isinstance(matrix, Ncm.Matrix)
    return np.array(matrix.dup_array()).reshape(len(_CASES), 2, _GOLDEN_N, 3)


def _build_new(ra_min, ra_max, dec_min, dec_max):
    """Flat scheme + a per-galaxy data. The footprint is config, not a model,
    so an empty mset suffices."""
    flat = Nc.GalaxyPositionFactorFlat.new(ra_min, ra_max, dec_min, dec_max)
    mset = Ncm.MSet.empty_new()
    data = Nc.GalaxyPositionFactorData.new(flat, mset)
    return flat, mset, data


_GEN_INTEG_PARITY_FROZEN = _load_golden()


@pytest.mark.parametrize("ra_min,ra_max,dec_min,dec_max", _CASES)
@pytest.mark.parametrize("use_lnp", [False, True])
def test_gen_and_integ_parity_legacy(ra_min, ra_max, dec_min, dec_max, use_lnp):
    """gen() draws and integ() evaluations are checked against a frozen
    legacy draw/integrand sequence (rtol=1e-12: dec is drawn via asin() and
    integ's density evaluates cos(dec) -- not bit-exact across platforms;
    ra is a plain linear transform of a uniform draw, but kept at the same
    tolerance for simplicity -- see module docstring)."""
    flat, mset, new_data = _build_new(ra_min, ra_max, dec_min, dec_max)

    new_integ = flat.integ(mset, use_lnp)

    rng_new = Ncm.RNG.seeded_new(None, 123)
    case_idx = _CASES.index((ra_min, ra_max, dec_min, dec_max))
    frozen = _GEN_INTEG_PARITY_FROZEN[case_idx, int(use_lnp)]
    for expected_ra, expected_dec, expected_integ in frozen:
        flat.gen(mset, new_data, rng_new)
        # Identical footprint sampler + same seed => identical draws.
        assert_allclose(new_data.ra, expected_ra, rtol=1.0e-12, atol=0.0)
        assert_allclose(new_data.dec, expected_dec, rtol=1.0e-12, atol=0.0)
        # Identical density => identical integrand value.
        assert_allclose(
            new_integ.eval(new_data), expected_integ, rtol=1.0e-12, atol=0.0
        )


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

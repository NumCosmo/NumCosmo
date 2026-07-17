#!/usr/bin/env python
#
# test_galaxy_sd_redshift_population.py
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

"""Standalone unit tests for the population redshift distribution models.

These exercise ``NcGalaxyRedshiftPop`` (LSST-SRD variant) in isolation
from any calculator: the marginal ``P(z|I)`` is validated for normalization over
``z``, consistency of ``ln_eval`` with ``eval`` and of ``gen()`` with the pdf, and
(historically) golden parity against the legacy ``NcGalaxySDTrueRedshiftLSSTSRD``
oracle -- the two shared identical math.

FROZEN REFERENCE VALUES: the parity documented above was proven by running
both engines live and is captured, not re-derived, in ``test_parity_eval``/
``test_parity_gen`` below. Every comparison here was bit-for-bit exact
(``rtol=0, atol=0``) between the new and legacy engines. Values were
captured from an actual passing run of this file's original
legacy-comparison code, at git rev ``77313f22`` (2026-07-16), then legacy
(``NcGalaxySDTrueRedshiftLSSTSRD``) construction was removed so these tests
no longer depend on legacy at runtime -- legacy is slated for deletion in a
follow-up PR. The captured sequences are stored as ``Ncm.Matrix`` binfiles
(``data/truth_tables/wl/``) rather than inline literals; see
``_load_eval_golden``/``_load_gen_golden``.
"""

import pytest
import numpy as np
from numpy.testing import assert_allclose

from numcosmo_py import Ncm, Nc

Ncm.cfg_init()

# (type ctor, alpha, beta, z0) for each LSST-SRD parametrization.
_VARIANTS = ["y1_source", "y1_lens", "y10_source", "y10_lens"]


def _new(variant):
    """Build a new Population LSST-SRD model for the named parametrization."""
    return getattr(Nc.GalaxyRedshiftPopLSSTSRD, f"new_{variant}")()


@pytest.mark.parametrize("variant", _VARIANTS)
def test_normalization(variant):
    """eval(z) integrates to unity over the [z_min, z_max] support."""
    model = _new(variant)
    z_min, z_max = model.get_lim()
    zs = np.linspace(z_min, z_max, 200001)
    p = np.array([model.eval(z) for z in zs])
    assert_allclose(np.trapezoid(p, zs), 1.0, rtol=1.0e-5)


@pytest.mark.parametrize("variant", _VARIANTS)
def test_ln_eval_consistency(variant):
    """ln_eval(z) equals log(eval(z)) across the support."""
    model = _new(variant)
    z_min, z_max = model.get_lim()
    zs = np.linspace(max(z_min, 1.0e-3), z_max, 512)
    ln_p = np.array([model.ln_eval(z) for z in zs])
    p = np.array([model.eval(z) for z in zs])
    assert_allclose(ln_p, np.log(p), rtol=1.0e-12)


@pytest.mark.parametrize("variant", _VARIANTS)
def test_gen_consistency(variant):
    """Empirical <z> from gen() matches the pdf <z>."""
    model = _new(variant)
    z_min, z_max = model.get_lim()
    zs = np.linspace(z_min, z_max, 200001)
    p = np.array([model.eval(z) for z in zs])
    mean_z = np.trapezoid(zs * p, zs)

    rng = Ncm.RNG.seeded_new(None, 1234)
    n_samples = 200000
    s = np.array([model.gen(rng) for _ in range(n_samples)])
    # ~5 sigma on the mean; the distribution std is O(z0/alpha) ~ O(1).
    assert_allclose(s.mean(), mean_z, atol=5.0 * s.std() / np.sqrt(n_samples))


_EVAL_GOLDEN_FILE = "truth_tables/wl/nc_galaxy_redshift_pop_eval_parity.bin"


def _load_eval_golden() -> np.ndarray:
    """Load the frozen eval/ln_eval sequences as a (len(_VARIANTS), 2, 1024) array."""
    path = Ncm.cfg_get_data_filename(_EVAL_GOLDEN_FILE, True)
    ser = Ncm.Serialize.new(Ncm.SerializeOpt.NONE)
    matrix = ser.from_binfile(path)
    assert isinstance(matrix, Ncm.Matrix)
    return np.array(matrix.dup_array()).reshape(len(_VARIANTS), 2, 1024)


_PARITY_EVAL_FROZEN = _load_eval_golden()


@pytest.mark.parametrize("variant", _VARIANTS)
def test_parity_eval(variant):
    """eval / ln_eval are checked against the frozen legacy TrueRedshift
    oracle (rtol=1e-12: eval/ln_eval call pow()/exp() directly -- not
    bit-exact across platforms, see module docstring)."""
    model = _new(variant)
    z_min, z_max = model.get_lim()
    zs = np.linspace(max(z_min, 1.0e-3), z_max, 1024)
    new_p = np.array([model.eval(z) for z in zs])
    idx = _VARIANTS.index(variant)
    assert_allclose(new_p, _PARITY_EVAL_FROZEN[idx, 0], rtol=1.0e-12, atol=1.0e-12)

    new_ln = np.array([model.ln_eval(z) for z in zs])
    assert_allclose(new_ln, _PARITY_EVAL_FROZEN[idx, 1], rtol=1.0e-12, atol=1.0e-12)


_GEN_GOLDEN_FILE = "truth_tables/wl/nc_galaxy_redshift_pop_gen_parity.bin"


def _load_gen_golden() -> np.ndarray:
    """Load the frozen gen() draw sequences as a (len(_VARIANTS), 5000) array."""
    path = Ncm.cfg_get_data_filename(_GEN_GOLDEN_FILE, True)
    ser = Ncm.Serialize.new(Ncm.SerializeOpt.NONE)
    matrix = ser.from_binfile(path)
    assert isinstance(matrix, Ncm.Matrix)
    return np.array(matrix.dup_array()).reshape(len(_VARIANTS), 5000)


_PARITY_GEN_FROZEN = _load_gen_golden()


@pytest.mark.parametrize("variant", _VARIANTS)
def test_parity_gen(variant):
    """gen() is checked against the frozen legacy oracle draw sequence,
    under a shared seed (rtol=1e-12: gen() draws from a GSL gamma-variate
    sampler, ncm_rng_gamma_gen(), whose rejection loop is not bit-exact
    across platforms -- see module docstring)."""
    model = _new(variant)
    rng_new = Ncm.RNG.seeded_new(None, 99)
    new_s = np.array([model.gen(rng_new) for _ in range(5000)])
    idx = _VARIANTS.index(variant)
    assert_allclose(new_s, _PARITY_GEN_FROZEN[idx], rtol=1.0e-12, atol=1.0e-12)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

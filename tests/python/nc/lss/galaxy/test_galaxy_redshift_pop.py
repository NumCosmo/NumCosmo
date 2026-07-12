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
golden parity against the legacy ``NcGalaxySDTrueRedshiftLSSTSRD`` oracle (the two
share identical math; this pins the parallel-dev copy until the legacy is removed).
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


def _new_legacy(variant):
    """Build the legacy TrueRedshift LSST-SRD oracle for the same parametrization."""
    return getattr(Nc.GalaxySDTrueRedshiftLSSTSRD, f"new_{variant}")()


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


@pytest.mark.parametrize("variant", _VARIANTS)
def test_parity_eval(variant):
    """eval / ln_eval reproduce the legacy TrueRedshift oracle exactly."""
    model = _new(variant)
    legacy = _new_legacy(variant)
    z_min, z_max = model.get_lim()
    zs = np.linspace(max(z_min, 1.0e-3), z_max, 1024)
    new_p = np.array([model.eval(z) for z in zs])
    old_p = np.array([legacy.integ(z) for z in zs])
    assert_allclose(new_p, old_p, rtol=0.0, atol=0.0)

    new_ln = np.array([model.ln_eval(z) for z in zs])
    old_ln = np.array([legacy.ln_integ(z) for z in zs])
    assert_allclose(new_ln, old_ln, rtol=0.0, atol=0.0)


@pytest.mark.parametrize("variant", _VARIANTS)
def test_parity_gen(variant):
    """gen() reproduces the legacy oracle bit-for-bit under a shared seed."""
    model = _new(variant)
    legacy = _new_legacy(variant)
    rng_new = Ncm.RNG.seeded_new(None, 99)
    rng_old = Ncm.RNG.seeded_new(None, 99)
    new_s = np.array([model.gen(rng_new) for _ in range(5000)])
    old_s = np.array([legacy.gen(rng_old) for _ in range(5000)])
    assert_allclose(new_s, old_s, rtol=0.0, atol=0.0)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

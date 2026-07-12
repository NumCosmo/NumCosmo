#!/usr/bin/env python
#
# test_galaxy_sd_shape_intrinsic.py
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

"""Standalone unit tests for the intrinsic-ellipticity distribution models.

These exercise ``NcGalaxySDShapeIntrinsic`` (Gauss and Beta variants) in
isolation from any calculator or likelihood: the distribution of
``x = |chi_I|^2 in [0, 1]`` is validated for normalization, internal
consistency of ``e_rms()`` against the pdf, and consistency of ``gen()``
against the pdf.
"""

import pytest
import numpy as np
from numpy.testing import assert_allclose
from scipy.special import roots_jacobi  # type: ignore # pylint: disable=no-name-in-module

from numcosmo_py import Ncm, Nc

Ncm.cfg_init()


def _make(name):
    """Build a prepared (model, data) pair for the requested variant."""
    if name == "gauss":
        model = Nc.GalaxySDShapeIntrinsicGauss.new()
    elif name == "beta":
        model = Nc.GalaxySDShapeIntrinsicBeta.new()
    else:
        raise ValueError(name)
    data = Nc.GalaxySDShapeIntrinsicData.new(model)
    model.prepare(data)
    return model, data


def _gj_moment(model, data, power, n=128):
    r"""Integrate ``x**power * eval_p(x)`` over [0, 1] via Gauss-Jacobi.

    ``eval_p = (1-x)^a x^b * eval_residual(x)``, so integrating the residual
    against the Jacobi weight (which absorbs the integrable endpoint
    singularities) is exact for both the Gauss (a=b=0) and Beta variants.
    """
    a, b = data.jacobi_a, data.jacobi_b
    t, w = roots_jacobi(n, a, b)  # weight (1-t)^a (1+t)^b on [-1, 1]
    x = 0.5 * (t + 1.0)
    scale = 2.0 ** (-(a + b + 1.0))  # remap to (1-x)^a x^b on [0, 1], dx = dt/2
    res = np.array([model.eval_residual(data, xi) for xi in x])
    return scale * np.sum(w * (x**power) * res)


@pytest.mark.parametrize("name", ["gauss", "beta"])
def test_normalization(name):
    """eval_p integrates to unity over x in [0, 1]."""
    model, data = _make(name)
    assert_allclose(_gj_moment(model, data, 0), 1.0, rtol=1.0e-6)


@pytest.mark.parametrize("name", ["gauss", "beta"])
def test_e_rms_consistency(name):
    """e_rms() equals sqrt(0.5 * <x>) computed from the pdf itself."""
    model, data = _make(name)
    mean_x = _gj_moment(model, data, 1)
    assert_allclose(model.e_rms(data), np.sqrt(0.5 * mean_x), rtol=1.0e-6)


@pytest.mark.parametrize("name", ["gauss", "beta"])
def test_gen_consistency(name):
    """Empirical <e1^2 + e2^2> from gen() matches the pdf <x>."""
    model, data = _make(name)
    mean_x = _gj_moment(model, data, 1)
    rng = Ncm.RNG.seeded_new(None, 42)
    n_samples = 200000
    acc = 0.0
    for _ in range(n_samples):
        e1, e2 = model.gen(data, rng)
        acc += e1 * e1 + e2 * e2
    # Monte-Carlo std of the mean ~ sqrt(Var(x)/N); x in [0,1] => Var <= 0.25,
    # so 4 sigma ~ 4 * 0.5 / sqrt(N) ~ 4.5e-3. Use a comfortable absolute band.
    assert_allclose(acc / n_samples, mean_x, atol=5.0e-3)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

#!/usr/bin/env python
#
# test_galaxy_shape_pop.py
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

These exercise ``NcGalaxyShapePop`` (Gauss and Beta variants) in
isolation from any calculator or likelihood: the distribution of
``x = |chi_I|^2 in [0, 1]`` is validated for normalization, internal
consistency of ``e_rms()`` against the pdf, and consistency of ``gen()``
against the pdf.
"""

import pytest
import numpy as np
from numpy.testing import assert_allclose
from scipy import integrate

from numcosmo_py import Ncm, Nc

Ncm.cfg_init()


def _make(name):
    """Build a prepared (model, data) pair for the requested variant."""
    if name == "gauss":
        model = Nc.GalaxyShapePopGauss.new()
        data = Nc.GalaxyShapePopData.new(model)
    elif name == "beta":
        model = Nc.GalaxyShapePopBeta.new()
        data = Nc.GalaxyShapePopData.new(model)
    elif name == "gauss_local":
        model = Nc.GalaxyShapePopGaussLocal.new()
        data = Nc.GalaxyShapePopData.new(model)
        model.data_set(data, 0.2)
    else:
        raise ValueError(name)
    model.prepare(data)
    return model, data


def _moment(model, data, power):
    r"""Integrate ``x**power * eval_p(x)`` over [0, 1] via adaptive quadrature.

    ``eval_p`` is each model's own normalized density, evaluated directly (no
    weight/residual decomposition): QUADPACK's adaptive subdivision handles
    the mild integrable endpoint behavior Beta can have without needing any
    special-purpose quadrature.
    """
    result, _ = integrate.quad(
        lambda x: (x**power) * model.eval_p(data, x),
        0.0,
        1.0,
        epsabs=1.0e-12,
        epsrel=1.0e-10,
    )
    return result


@pytest.mark.parametrize("name", ["gauss", "beta", "gauss_local"])
def test_normalization(name):
    """eval_p integrates to unity over x in [0, 1]."""
    model, data = _make(name)
    assert_allclose(_moment(model, data, 0), 1.0, rtol=1.0e-6)


@pytest.mark.parametrize("name", ["gauss", "beta", "gauss_local"])
def test_e_rms_consistency(name):
    """e_rms() equals sqrt(0.5 * <x>) computed from the pdf itself."""
    model, data = _make(name)
    mean_x = _moment(model, data, 1)
    assert_allclose(model.e_rms(data), np.sqrt(0.5 * mean_x), rtol=1.0e-6)


@pytest.mark.parametrize("name", ["gauss", "beta", "gauss_local"])
def test_gen_consistency(name):
    """Empirical <e1^2 + e2^2> from gen() matches the pdf <x>."""
    model, data = _make(name)
    mean_x = _moment(model, data, 1)
    rng = Ncm.RNG.seeded_new(None, 42)
    n_samples = 200000
    acc = 0.0
    for _ in range(n_samples):
        e1, e2 = model.gen(data, rng)
        acc += e1 * e1 + e2 * e2
    # Monte-Carlo std of the mean ~ sqrt(Var(x)/N); x in [0,1] => Var <= 0.25,
    # so 4 sigma ~ 4 * 0.5 / sqrt(N) ~ 4.5e-3. Use a comfortable absolute band.
    assert_allclose(acc / n_samples, mean_x, atol=5.0e-3)


@pytest.mark.parametrize("name", ["gauss", "beta", "gauss_local"])
def test_eval_p_rho2_matches_eval_p(name):
    """eval_p_rho2(rho2) agrees with eval_p(x) at x = rho2/(1+rho2), whether
    the model overrides eval_p_rho2 (Beta) or inherits the generic default
    (Gauss, GaussLocal)."""
    model, data = _make(name)
    for rho2 in (1.0e-3, 0.1, 1.0, 5.0, 50.0):
        x = rho2 / (1.0 + rho2)
        assert_allclose(
            model.eval_p_rho2(data, rho2), model.eval_p(data, x), rtol=1.0e-10
        )


@pytest.mark.parametrize(
    "alpha,beta",
    [(3000.0, 7000.0), (500.0, 500.0), (6000.0, 4000.0), (990.0, 10.0)],
)
def test_beta_concentrated_no_overflow(alpha, beta):
    """eval_p/eval_p_rho2 stay finite and normalized for concentrated Beta
    populations (large alpha, beta).

    Evaluating norm=1/B(alpha,beta) and x**(alpha-1) as separate factors
    overflows/underflows for alpha/beta of a few hundred or more (already
    reachable within the class's own declared range, e.g. alpha=beta=500
    pushes -ln B(alpha,beta) to ~695, right at exp()'s overflow edge),
    silently producing NaN (0*inf). eval_p must instead accumulate in
    log-space.
    """
    model = Nc.GalaxyShapePopBeta.new()
    model["alpha"] = alpha
    model["beta"] = beta
    data = Nc.GalaxyShapePopData.new(model)
    model.prepare(data)

    mode = max((alpha - 1.0) / (alpha + beta - 2.0), 1.0e-6)
    for x in (mode * 0.999, mode, min(mode * 1.001, 1.0 - 1.0e-9)):
        p = model.eval_p(data, x)
        assert np.isfinite(p)
        assert p >= 0.0

    assert_allclose(_moment(model, data, 0), 1.0, rtol=1.0e-5)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

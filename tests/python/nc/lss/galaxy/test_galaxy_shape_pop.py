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
isolation from any calculator or likelihood: ``eval_p(r)`` is each
model's own r-marginal density, ``r = |chi_I| in [0, 1]``, normalized so
that ``integral_0^1 eval_p(r) dr = 1``. Validated here for normalization,
internal consistency of ``e_rms()`` against the pdf, and consistency of
``gen()`` against the pdf.
"""

import subprocess
import sys

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
    r"""Compute ``E[x**power]`` with ``x = r**2``, via adaptive quadrature in r.

    ``eval_p(r)`` is each model's own normalized r-marginal density, so
    ``E[x**power] = E[r**(2*power)] = integral_0^1 r**(2*power) * eval_p(r) dr``.
    Evaluated directly (no weight/residual decomposition): QUADPACK's
    adaptive subdivision handles the mild integrable endpoint behavior Beta
    can have without needing any special-purpose quadrature.
    """
    result, _ = integrate.quad(
        lambda r: (r ** (2 * power)) * model.eval_p(data, r),
        0.0,
        1.0,
        epsabs=1.0e-12,
        epsrel=1.0e-10,
    )
    return result


@pytest.mark.parametrize("name", ["gauss", "beta", "gauss_local"])
def test_normalization(name):
    """eval_p integrates to unity over r in [0, 1]."""
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


@pytest.mark.parametrize(
    "alpha,beta",
    [(3000.0, 7000.0), (500.0, 500.0), (6000.0, 4000.0), (990.0, 10.0)],
)
def test_beta_concentrated_no_overflow(alpha, beta):
    """eval_p stays finite and normalized for concentrated Beta populations
    (large alpha, beta).

    Evaluating norm=1/B(alpha,beta) and r**(alpha-1) as separate factors
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

    r_mode = (alpha - 1.0) / (alpha + beta - 2.0)
    mode = max(r_mode, 1.0e-6)
    for r in (mode * 0.999, mode, min(mode * 1.001, 1.0 - 1.0e-9)):
        p = model.eval_p(data, r)
        assert np.isfinite(p)
        assert p >= 0.0

    assert_allclose(_moment(model, data, 0), 1.0, rtol=1.0e-5)


@pytest.mark.parametrize(
    "alpha,beta",
    [(2.0, 3.0), (1.05, 4.7833333333333333), (1.4, 1.6), (50.0, 25.0)],
)
def test_beta_e_rms_mode(alpha, beta):
    """get_e_rms matches its closed-form definition (r=|chi_I| ~
    Beta(alpha,beta), x=r^2) and the first moment of eval_p itself;
    get_mode matches the natural mode(r), literally the argmax of
    eval_p(r) itself under the r-native contract."""
    model = Nc.GalaxyShapePopBeta.new()
    model["alpha"] = alpha
    model["beta"] = beta
    data = Nc.GalaxyShapePopData.new(model)
    model.prepare(data)

    e_rms = model.get_e_rms()
    mode = model.get_mode()

    mean_r2 = alpha * (alpha + 1.0) / ((alpha + beta) * (alpha + beta + 1.0))
    assert_allclose(e_rms, np.sqrt(0.5 * mean_r2))

    if alpha > 1.0 and beta > 1.0:
        assert_allclose(mode, (alpha - 1.0) / (alpha + beta - 2.0))
    else:
        assert_allclose(mode, 0.0)

    assert_allclose(2.0 * e_rms**2, _moment(model, data, 1), rtol=1.0e-5)


def test_eval_p_array_matches_eval_p():
    """eval_p_array() batches eval_p() over an array of r values; the GI
    binding always starts from an unallocated output array (Python has no
    way to pass a preallocated GArray through the (out callee-allocates)
    convention), so this also exercises the eval_p_array() default
    implementation's own allocate-on-first-use branch."""
    model, data = _make("gauss")
    rs = [0.05, 0.2, 0.5, 0.9]

    p_array = model.eval_p_array(data, rs)
    p_single = [model.eval_p(data, r) for r in rs]

    assert_allclose(p_array, p_single, rtol=1.0e-12)


def test_get_sigma_raises_for_unsupported_population():
    """get_sigma() is a per-@data capability, not every model resolves it
    (e.g. Beta has no Gaussian width): nc_galaxy_shape_pop_get_sigma()
    reports this with a fatal g_error, a real process abort, so it is
    checked in a subprocess (see test_galaxy_shape_factor_cgf.py's
    test_non_gaussian_population_gate for the same pattern)."""
    script = (
        "from numcosmo_py import Nc, Ncm\n"
        "Ncm.cfg_init()\n"
        "model = Nc.GalaxyShapePopBeta.new()\n"
        "data = Nc.GalaxyShapePopData.new(model)\n"
        "model.prepare(data)\n"
        "model.get_sigma(data)\n"
    )
    result = subprocess.run(
        [sys.executable, "-c", script],
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode != 0
    assert "does not support a Gaussian width sigma" in result.stderr


def test_beta_mset_func_list():
    """The NcGalaxyShapePopBeta:e_rms/mode NcmMSetFuncList entries match the
    direct get_e_rms()/get_mode() accessors."""
    model = Nc.GalaxyShapePopBeta.new()
    model["alpha"] = 3.0
    model["beta"] = 7.0
    mset = Ncm.MSet.new_array([model])

    e_rms_func = Ncm.MSetFuncList.new("NcGalaxyShapePopBeta:e_rms", None)
    mode_func = Ncm.MSetFuncList.new("NcGalaxyShapePopBeta:mode", None)

    assert_allclose(e_rms_func.eval0(mset), model.get_e_rms())
    assert_allclose(mode_func.eval0(mset), model.get_mode())


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

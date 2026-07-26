#!/usr/bin/env python
#
# test_catalog_stats.py
#
# Sun Jul 26 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_catalog_stats.py
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
# with this program.  If not, see <http://www.gnu.org/licenses/>.

"""Unit tests for numcosmo_py.catalog_stats.

The asymmetric-bound helpers are tested against a small hand-written fake
distribution (not a real Ncm.StatsDist1d) so the expected numbers can be
computed analytically instead of relying on KDE estimation noise.
"""

import math

import numpy as np
import pytest

from numcosmo_py import Ncm
from numcosmo_py.catalog_stats import (
    SIGMA_LEVELS,
    DerivedStat,
    asymmetric_bounds,
    median_bounds,
    parse_variable_bindings,
    resolve_param,
    safe_eval_mode,
    stat_center_and_bounds,
)

Ncm.cfg_init()

NADD_VALS = 2


@pytest.fixture(name="mset")
def fixture_mset():
    """Return an MSet for a 3-parameter MVN model, all parameters free."""
    model = Ncm.ModelMVND.new(dim=3)
    mset = Ncm.MSet.new_array([model])
    mset.param_set_all_ftype(Ncm.ParamType.FREE)
    mset.prepare_fparam_map()
    return mset


class _FakeQuantileDist:
    """A fake with only the two quantile methods `median_bounds` needs."""

    def __init__(self, finv):
        self._finv = finv

    def eval_inv_pdf(self, u: float) -> float:
        """Quantile function."""
        return self._finv(u)

    def eval_inv_pdf_tail(self, v: float) -> float:
        """Quantile function evaluated from the right tail."""
        return self._finv(1.0 - v)


class _FakeDist(_FakeQuantileDist):
    """A fake exposing the full interface `asymmetric_bounds` needs.

    Uses finv(u) = u**2 on [0, 1], whose CDF is F(x) = sqrt(x). This is
    monotonic and asymmetric enough to distinguish lower from upper bounds.
    """

    def __init__(self):
        super().__init__(lambda u: u**2)

    def get_xi(self) -> float:
        """Domain lower bound."""
        return 0.0

    def get_xf(self) -> float:
        """Domain upper bound."""
        return 1.0

    def eval_pdf(self, x: float) -> float:
        """CDF, i.e. P(X <= x)."""
        return math.sqrt(x)

    def eval_p(self, x: float) -> float:
        """Density, i.e. d/dx CDF(x). Monotonically decreasing on (0, 1]."""
        return math.inf if x <= 0.0 else 0.5 / math.sqrt(x)


def test_derived_stat_values():
    """DerivedStat exposes the three expected string values."""
    assert {s.value for s in DerivedStat} == {"median", "mode", "bestfit"}


def test_sigma_levels_are_standard_two_sided_cis():
    """SIGMA_LEVELS are the usual 1/2/3-sigma two-sided confidence levels."""
    assert len(SIGMA_LEVELS) == 3
    assert SIGMA_LEVELS[0] == pytest.approx(0.682689492137086, rel=1e-12)
    assert SIGMA_LEVELS[1] == pytest.approx(0.954499736103642, rel=1e-12)
    assert SIGMA_LEVELS[2] == pytest.approx(0.997300203936740, rel=1e-12)


def test_resolve_param_by_name(mset):
    """A free-parameter name resolves to its Ncm.MSetPIndex and column index."""
    pi, pindex = resolve_param(mset, NADD_VALS, "mu_0")
    assert pi is not None
    assert pindex == NADD_VALS + mset.fparam_get_fpi(pi.mid, pi.pid)
    assert pindex == NADD_VALS


def test_resolve_param_by_numeric_index_in_add_vals(mset):
    """A numeric index inside the add-values range resolves with pi=None."""
    pi, pindex = resolve_param(mset, NADD_VALS, "0")
    assert pi is None
    assert pindex == 0


def test_resolve_param_by_numeric_index_in_fparams(mset):
    """A numeric index inside the free-parameters range resolves a valid pi."""
    pi, pindex = resolve_param(mset, NADD_VALS, str(NADD_VALS + 1))
    assert pi is not None
    assert pindex == NADD_VALS + 1


def test_resolve_param_numeric_index_out_of_range(mset):
    """An out-of-range numeric index raises a clear ValueError."""
    total_len = mset.fparams_len() + NADD_VALS
    with pytest.raises(ValueError, match="Invalid parameter index"):
        resolve_param(mset, NADD_VALS, str(total_len))


def test_resolve_param_unknown_name(mset):
    """An unknown, non-numeric name raises a clear ValueError."""
    with pytest.raises(ValueError, match="not found"):
        resolve_param(mset, NADD_VALS, "does_not_exist")


def test_parse_variable_bindings_basic(mset):
    """Multiple "name=parameter" bindings resolve to the right column indices."""
    var_pindex = parse_variable_bindings(
        mset, NADD_VALS, ["x = mu_0", "y=mu_1"], "--variable"
    )
    _, expected_x = resolve_param(mset, NADD_VALS, "mu_0")
    _, expected_y = resolve_param(mset, NADD_VALS, "mu_1")
    assert var_pindex == {"x": expected_x, "y": expected_y}


def test_parse_variable_bindings_bare_shorthand(mset):
    """A bare "parameter" binding is shorthand for "parameter=parameter"."""
    var_pindex = parse_variable_bindings(
        mset, NADD_VALS, ["mu_0", "y=mu_1"], "--variable"
    )
    _, expected_mu0 = resolve_param(mset, NADD_VALS, "mu_0")
    _, expected_y = resolve_param(mset, NADD_VALS, "mu_1")
    assert var_pindex == {"mu_0": expected_mu0, "y": expected_y}


@pytest.mark.parametrize("binding", ["=mu_0", "x=", "   =   "])
def test_parse_variable_bindings_invalid_syntax(mset, binding):
    """A binding missing a name or a parameter raises a ValueError."""
    with pytest.raises(ValueError, match="expected name=parameter"):
        parse_variable_bindings(mset, NADD_VALS, [binding], "--my-option")


def test_parse_variable_bindings_error_mentions_option_name(mset):
    """The error message names the CLI option that produced the bad binding."""
    with pytest.raises(ValueError, match="--my-option"):
        parse_variable_bindings(mset, NADD_VALS, ["=bad"], "--my-option")


def test_parse_variable_bindings_bare_unknown_parameter(mset):
    """A bare binding to a nonexistent parameter still fails via resolve_param."""
    with pytest.raises(ValueError, match="not found"):
        parse_variable_bindings(mset, NADD_VALS, ["does_not_exist"], "--variable")


def test_median_bounds_asymmetric():
    """median_bounds() matches the analytic quantiles of a known fake CDF."""
    finv = lambda u: u**2  # noqa: E731
    fake = _FakeQuantileDist(finv)
    median = finv(0.5)

    for pa in SIGMA_LEVELS:
        v = (1.0 - pa) * 0.5
        lo, hi = median_bounds(fake, pa, median)
        assert lo == pytest.approx(median - finv(v))
        assert hi == pytest.approx(finv(1.0 - v) - median)
        assert lo > 0.0
        assert hi > 0.0


def test_asymmetric_bounds_normal_case():
    """Away from the domain edges, bounds are symmetric-in-CDF around center."""
    fake = _FakeDist()
    center = 0.25  # this is finv(0.5), i.e. the median of the fake distribution.

    for pa in SIGMA_LEVELS:
        pa_2 = pa * 0.5
        lo, hi = asymmetric_bounds(fake, pa, center)
        assert lo == pytest.approx(center - (0.5 - pa_2) ** 2)
        assert hi == pytest.approx((0.5 + pa_2) ** 2 - center)


def test_asymmetric_bounds_right_skewed_near_minimum():
    """Center close to xi clips the lower bound at the domain edge."""
    fake = _FakeDist()
    center = 0.001  # eval_pdf(center) = sqrt(0.001) ~= 0.0316

    for pa in SIGMA_LEVELS:
        lo, hi = asymmetric_bounds(fake, pa, center)
        assert lo == pytest.approx(center - 0.0)
        assert hi == pytest.approx(pa**2 - center)


def test_asymmetric_bounds_left_skewed_near_maximum():
    """Center close to xf clips the upper bound at the domain edge."""
    fake = _FakeDist()
    center = 0.999  # eval_pdf(center) = sqrt(0.999) ~= 0.9995

    for pa in SIGMA_LEVELS:
        lo, hi = asymmetric_bounds(fake, pa, center)
        assert lo == pytest.approx(center - (1.0 - pa) ** 2)
        assert hi == pytest.approx(1.0 - center)


def _epdf_from_samples(samples) -> Ncm.StatsDist1dEPDF:
    epdf = Ncm.StatsDist1dEPDF.new(1.0e-3)
    for x in samples:
        epdf.add_obs(float(x))
    epdf.prepare()
    return epdf


def test_safe_eval_mode_interior_matches_eval_mode():
    """For a well-behaved (interior-mode) posterior, matches eval_mode()."""
    rng = np.random.default_rng(1)
    epdf = _epdf_from_samples(rng.normal(loc=0.0, scale=1.0, size=3000))

    assert safe_eval_mode(epdf) == pytest.approx(epdf.eval_mode(), abs=1e-6)


def test_safe_eval_mode_boundary_case_does_not_crash():
    """A monotonic (boundary-peaked) density does not crash the process.

    Calling `epdf.eval_mode()` directly here reproduces a real bug found
    while dogfooding `catalog derived-error` on a skewed halo-mass
    posterior: NumCosmo's GSL bracketing minimizer requires the mode to be
    strictly interior to [xi, xf], and aborts the whole process (fatal
    g_error, not a catchable exception) when the true maximum sits at the
    domain edge instead. safe_eval_mode() must detect that case up front
    and avoid ever calling eval_mode() on it.
    """
    rng = np.random.default_rng(2)
    epdf = _epdf_from_samples(rng.exponential(scale=1.0, size=3000))

    mode = safe_eval_mode(epdf)
    assert mode == pytest.approx(epdf.get_xi(), abs=1e-6)


def test_safe_eval_mode_trivial_domain():
    """A degenerate (xi == xf) domain returns that point without scanning."""
    epdf = _epdf_from_samples([5.0] * 5)
    assert epdf.get_xi() == epdf.get_xf()
    assert safe_eval_mode(epdf) == pytest.approx(epdf.get_xi())


def test_stat_center_and_bounds_dispatch():
    """stat_center_and_bounds() picks the right center/bound-fn per DerivedStat."""
    fake = _FakeDist()
    bestfit_center = 0.42

    median_center, median_bnds = stat_center_and_bounds(
        DerivedStat.MEDIAN, fake, bestfit_center
    )
    assert median_center == pytest.approx(fake.eval_inv_pdf(0.5))
    assert median_bnds == [
        median_bounds(fake, pa, median_center) for pa in SIGMA_LEVELS
    ]

    mode_center, mode_bnds = stat_center_and_bounds(
        DerivedStat.MODE, fake, bestfit_center
    )
    assert mode_center == pytest.approx(safe_eval_mode(fake))
    assert mode_bnds == [
        asymmetric_bounds(fake, pa, mode_center) for pa in SIGMA_LEVELS
    ]

    bf_center, bf_bnds = stat_center_and_bounds(
        DerivedStat.BESTFIT, fake, bestfit_center
    )
    assert bf_center == bestfit_center
    assert bf_bnds == [
        asymmetric_bounds(fake, pa, bestfit_center) for pa in SIGMA_LEVELS
    ]

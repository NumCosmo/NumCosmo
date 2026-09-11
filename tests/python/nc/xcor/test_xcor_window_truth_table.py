#!/usr/bin/env python
#
# test_xcor_window_truth_table.py
#
# Wed Aug 27 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_xcor_window_truth_table.py
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

"""Check the radial integral against certified Arb values.

The radial integral

.. math:: I_\\ell(k) = \\int W(\\chi)\\, j_\\ell(k\\chi)\\, \\mathrm{d}\\chi

is where three quarters of an ``NcXcorSolver`` run is spent, and it has no
independent check inside NumCosmo -- every route to it goes through the same
Levin machinery. ``data/truth_tables/xcor/xcor_window_ilk.json.gz`` holds it
computed in Arb instead, certified to a relative ball radius of 2e-26, for
every analytic window on a grid that follows the peak across multipoles. This
test is the guard to run before and after touching the sbessel integrators.

Two companion tables, ``xcor_window_ilk_d1.json.gz`` and ``_d2.json.gz``, carry
the same integral with :math:`j_\\ell^{(d)}` in place of :math:`j_\\ell`. They are
what certifies the components weighted by a Bessel derivative -- the
redshift-space distortion term above all.

Regenerating the tables needs FLINT; running this test does not. See
``tests/tools/make_xcor_window_truth_table.py``.
"""

import gzip
import json
import pathlib

import numpy as np
import pytest
from numpy.testing import assert_allclose

from gi.repository import GLib

from numcosmo_py import Nc, Ncm

Ncm.cfg_init()

# Pinned to one worker under --dist loadgroup: this file is one of the
# xcor lane's memory peaks, and an xdist worker is its own session, so
# without this its cost is paid once per worker rather than once.
pytestmark = [pytest.mark.xcor, pytest.mark.xdist_group("window_truth")]

TRUTH_TABLE = "truth_tables/xcor/xcor_window_ilk.json.gz"


def _cases(name: str) -> list[str]:
    """The cases a committed table carries, read at collection time.

    Parametrizing from the table rather than from a list in this file is what makes
    an added case actually tested: the generator, the table and the suite then
    cannot disagree about what is certified. Collection must not fail when a table
    is absent -- a derivative order may not have been generated yet -- so a missing
    file yields no cases and the tests for it are simply not collected.
    """
    try:
        path = pathlib.Path(Ncm.cfg_get_data_filename(name, True))
    except (RuntimeError, GLib.Error):
        return []

    with gzip.open(path, "rt") as handle:
        return sorted(json.load(handle)["shapes"])


CASES = _cases(TRUTH_TABLE)

# The same integral with a Bessel-derivative weight, one table per order. These
# certify the path a redshift-space distortion term takes: its component carries
# NcXcorKernelComponent:bessel-deriv = 2, and the library reaches int W j_ell'' by
# integration by parts onto the window fit, a route nothing else in the tree checks
# against a proven value. The Arb side keeps the derivative on the Bessel factor
# instead, so the two agree only if both are right.
DERIV_TRUTH_TABLE = {
    1: "truth_tables/xcor/xcor_window_ilk_d1.json.gz",
    2: "truth_tables/xcor/xcor_window_ilk_d2.json.gz",
}

# Measured, not chosen: over all 252 entries the worst relative deviation where
# |I| is within three orders of its own peak is 6.9e-10, and the worst
# deviation measured against the peak is 5.7e-11. Both tolerances sit about a
# factor of 15 above that. The integrator's nominal reltol of 1e-13 is *not*
# the right anchor -- it does not reach it, and cannot: the documented
# conditioning floor of these integrands is around 2e-10.
RTOL = 1.0e-8

# Applied as a fraction of each multipole's own peak. Far below the peak a
# relative criterion is meaningless -- the grid reaches |I|/peak = 1e-70 -- so
# this is what those entries are actually held to.
ATOL_FRAC = 1.0e-11


def _shape_of(entry, case: str) -> str:
    """The closed form a case carries.

    Several cases are the same form at different centres and widths -- the LSST-SRD
    locations, their paired top-hats -- so the case key is not the shape name. Tables
    written before that distinction existed carry no ``shape`` field and the key is
    the form.
    """
    return entry.get("shape", case)


def _make_kernel(shape, dist, ps, sbi, ctor):
    """Build one analytic kernel from the constructor arguments in the table."""
    match shape:
        case "gauss":
            return Nc.XcorKernelAnalyticGauss.new_full(dist, ps, *ctor, sbi)
        case "tophat":
            return Nc.XcorKernelAnalyticTophat.new_full(dist, ps, *ctor, sbi)
        case "tophat_smooth":
            return Nc.XcorKernelAnalyticTophatSmooth.new_full(dist, ps, *ctor, sbi)
        case "student_t":
            return Nc.XcorKernelAnalyticStudentT.new_full(dist, ps, *ctor, sbi)
        case "power_exp":
            return Nc.XcorKernelAnalyticPowerExp.new_full(dist, ps, *ctor, sbi)
        case "lensing":
            return Nc.XcorKernelAnalyticLensing.new_full(dist, ps, *ctor, sbi)
        case "multi":
            mu, sigma, weight, n_sigma = ctor

            return Nc.XcorKernelAnalyticMulti.new_full(
                dist,
                ps,
                Ncm.Vector.new_array(mu),
                Ncm.Vector.new_array(sigma),
                Ncm.Vector.new_array(weight),
                n_sigma,
                sbi,
            )
        case _:
            raise ValueError(f"truth table names an unknown shape {shape!r}")


@pytest.fixture(name="truth_table", scope="module")
def fixture_truth_table() -> dict:
    """Load the certified table once for the module."""
    path = pathlib.Path(Ncm.cfg_get_data_filename(TRUTH_TABLE, True))

    with gzip.open(path, "rt") as f:
        return json.load(f)


@pytest.fixture(name="cosmo_bits", scope="module")
def fixture_cosmo_bits() -> tuple:
    """A prepared distance and a power spectrum, both constructor filler.

    The analytic windows are closed forms in chi alone: neither object reaches
    the value this test compares. They are built here rather than through
    ``Cosmology.default()`` because that imports the FFTW wisdom file, which
    costs more than every assertion in this file put together.
    """
    cosmo = Nc.HICosmoDEXcdm.new()
    dist = Nc.Distance.new(5.0)
    dist.prepare(cosmo)
    ps = Ncm.PowspecAnalytic.new(
        Ncm.PowspecAnalyticShape.BBKS, Ncm.PowspecAnalyticGrowth.LCDM
    )

    return cosmo, dist, ps


def test_table_is_certified_far_below_the_tolerance(truth_table: dict) -> None:
    """The reference must be a reference: its own uncertainty cannot matter.

    Guards a regeneration that silently lowered the precision target -- the
    comparison would still pass while no longer checking anything. The target is per
    multipole, because near the turning point the precision needed grows steeply with
    ell: the high multipoles are certified to 1e-18 rather than 1e-25. Both checks
    below are what makes that safe, the second one being the one that matters.
    """
    targets = truth_table["target_rel"]
    targets = (
        targets if isinstance(targets, list) else [targets] * len(truth_table["ells"])
    )

    # Every target must sit far below the comparison this file makes, or the reference
    # stops being one.
    assert max(targets) < 1.0e-3 * RTOL

    for index, (ell, target) in enumerate(zip(truth_table["ells"], targets)):
        worst = max(
            (
                radius / abs(float(value))
                for shape in truth_table["shapes"].values()
                for value, radius in zip(shape["table"][index], shape["radius"][index])
                if float(value) != 0.0
            ),
            default=0.0,
        )

        assert worst < target, f"ell = {ell} certified only to {worst:.2e}"


@pytest.mark.parametrize("shape", CASES)
def test_radial_integral_matches_arb(
    shape: str, truth_table: dict, cosmo_bits: tuple
) -> None:
    """Compare I_ell(k) against Arb for one window, over every ell and k."""
    cosmo, dist, ps = cosmo_bits
    entry = truth_table["shapes"][shape]

    kernel = _make_kernel(
        _shape_of(entry, shape), dist, ps, Ncm.SBesselIntegratorLevin.new(0, 8), entry["ctor"]
    )
    kernel.set_l_limber(-1)
    kernel.prepare(cosmo)

    supports = [kernel.get_comp_support(comp) for comp in range(kernel.get_n_comps())]

    for index, ell in enumerate(truth_table["ells"]):
        expected = np.array([float(value) for value in entry["table"][index]])
        peak = np.abs(expected).max()

        integrator = Ncm.SBesselIntegratorLevin.new(ell, ell)

        got = np.array(
            [
                sum(
                    integrator.integrate_ell(
                        # The window is a hard zero outside its support and the
                        # k*chi round trip puts nodes a hair beyond the edge, a
                        # cliff the Chebyshev fit cannot resolve. The library's
                        # own path clamps for this reason; so does this one.
                        lambda _ud, chi, _k, c=comp, lo=low, hi=high: (
                            kernel.eval_W_comp(c, min(max(chi, lo), hi))
                        ),
                        low,
                        high,
                        k,
                        ell,
                        None,
                    )
                    for comp, (low, high) in enumerate(supports)
                )
                for k in entry["kvals"][index]
            ]
        )

        assert_allclose(
            got,
            expected,
            rtol=RTOL,
            atol=ATOL_FRAC * peak,
            err_msg=f"{shape} at ell = {ell}",
        )


@pytest.fixture(name="deriv_truth_tables", scope="module")
def fixture_deriv_truth_tables() -> dict:
    """Load the certified Bessel-derivative tables once for the module."""
    tables = {}

    for deriv, name in DERIV_TRUTH_TABLE.items():
        path = pathlib.Path(Ncm.cfg_get_data_filename(name, True))

        with gzip.open(path, "rt") as f:
            tables[deriv] = json.load(f)

    return tables


# Measured, not chosen, over both derivative orders, all seven shapes, every
# multipole of the tables and both block sizes. Worst relative deviation where
# |I| is within three orders of its own peak: 1.7e-10 at d = 1, 2.1e-9 at d = 2
# (student_t, ell = 2). Worst deviation measured against the multipole's own
# peak: 1.7e-10 at d = 1, 2.8e-9 at d = 2 (multi, ell = 2, on an entry sitting at
# 9e-5 of the peak). Both allowances below sit about a factor of four above the
# d = 2 numbers.
#
# The order matters and the direction is the expected one. At d = 1 the
# derivative costs nothing -- 1.7e-10 is the d = 0 floor of this same grid. At
# d = 2 it costs about an order, because j_ell'' is reached from a combination
# of neighbouring orders whose leading terms cancel; the two worst cells are the
# two where cancellation is structural anyway, ell = 2 (fewest oscillations to
# average over) and the multi window (disjoint bump groups). That is a property
# of the integrand, not a defect: it is the same conditioning floor
# #NcXcorKernel:peak-epsilon documents, one derivative deeper.
DERIV_RTOL = 1.0e-8
DERIV_ATOL_FRAC = 1.0e-8


@pytest.mark.parametrize("deriv", sorted(DERIV_TRUTH_TABLE))
def test_deriv_table_is_certified_far_below_the_tolerance(
    deriv: int, deriv_truth_tables: dict
) -> None:
    """Each derivative table must be a reference, on the same terms as the d = 0 one."""
    table = deriv_truth_tables[deriv]

    assert table["bessel_deriv"] == deriv

    targets = table["target_rel"]
    targets = targets if isinstance(targets, list) else [targets] * len(table["ells"])

    assert max(targets) < 1.0e-3 * DERIV_RTOL

    for index, (ell, target) in enumerate(zip(table["ells"], targets)):
        worst = max(
            (
                radius / abs(float(value))
                for shape in table["shapes"].values()
                for value, radius in zip(shape["table"][index], shape["radius"][index])
                if float(value) != 0.0
            ),
            default=0.0,
        )

        assert worst < target, f"deriv {deriv}, ell = {ell} certified only to {worst:.2e}"


@pytest.mark.parametrize("n_block", [1, 8])
@pytest.mark.parametrize("deriv", sorted(DERIV_TRUTH_TABLE))
@pytest.mark.parametrize("shape", CASES)
def test_radial_integral_deriv_matches_arb(
    shape: str, deriv: int, n_block: int, deriv_truth_tables: dict, cosmo_bits: tuple
) -> None:
    """int W j_ell^(d) against Arb, for the weights a derivative component carries.

    This is the certified check the RSD path lacked. A redshift-space term reaches
    the solver as a component with #NcXcorKernelComponent:bessel-deriv set to 2, and
    the library never differentiates a Bessel function to get there: it integrates by
    parts, so the derivative lands on the window's own Chebyshev fit -- the ODE forcing
    becomes x F'(x) or x F''(x) read off in the C^(2) basis -- plus the boundary terms
    that leaves. The reference does the opposite, evaluating j_ell^(d) pointwise from
    the recurrences in the *order* and integrating W j_ell^(d) directly. The two share
    the window and nothing else, so agreement here is evidence rather than a tautology.

    Both block sizes run because production RSD is batched: n_block = 1 isolates the
    single-multipole path and 8 exercises the shared operator the ell block reuses.
    """
    cosmo, dist, ps = cosmo_bits
    table = deriv_truth_tables[deriv]
    entry = table["shapes"][shape]

    kernel = _make_kernel(
        _shape_of(entry, shape), dist, ps, Ncm.SBesselIntegratorLevin.new(0, 8), entry["ctor"]
    )
    kernel.set_l_limber(-1)
    kernel.prepare(cosmo)

    supports = [kernel.get_comp_support(comp) for comp in range(kernel.get_n_comps())]

    for index, ell in enumerate(table["ells"]):
        expected = np.array([float(value) for value in entry["table"][index]])
        peak = np.abs(expected).max()

        integrator = Ncm.SBesselIntegratorLevin.new(ell, ell + n_block - 1)
        block = Ncm.Vector.new(n_block)
        got = []

        for k in entry["kvals"][index]:
            total = 0.0

            for comp, (low, high) in enumerate(supports):
                # The same clamp the d = 0 tests make, for the same reason: the
                # window is a hard zero outside its support and the k*chi round
                # trip puts nodes a hair beyond the edge.
                integrator.integrate_deriv(
                    lambda chi, _k, c=comp, lo=low, hi=high: (
                        kernel.eval_W_comp(c, min(max(chi, lo), hi))
                    ),
                    low,
                    high,
                    k,
                    deriv,
                    block,
                )
                total += block.get(0)

            got.append(total)

        assert_allclose(
            np.array(got),
            expected,
            rtol=DERIV_RTOL,
            atol=DERIV_ATOL_FRAC * peak,
            err_msg=f"{shape} at ell = {ell}, deriv = {deriv}, block = {n_block}",
        )


# The tau constraint of the panel solve. "rule" is the production setting: a
# panel takes it when its oscillation count exceeds tau-constraint-min-osc. "forced"
# bypasses the rule and turns the tau constraint on everywhere, including panels below
# the turning point where it is invalid; the guard in the integrator must catch
# those and redo them with Dirichlet data, so this mode is the guard's own test.
TAU_CONSTRAINT_MODES = ["rule", "forced"]
TAU_CONSTRAINT_MIN_OSC = 200.0


def _set_constraint(integrator: Ncm.SBesselIntegratorLevin, mode: str) -> None:
    if mode == "rule":
        integrator.set_tau_constraint_min_osc(TAU_CONSTRAINT_MIN_OSC)
    else:
        # The rule reapplies its own decision to every operator it configures, so
        # it has to be off for the solver-wide flag to reach the shallow panels.
        integrator.set_tau_constraint_min_osc(0.0)
        integrator.set_tau_constraint(True)


@pytest.mark.parametrize("mode", TAU_CONSTRAINT_MODES)
@pytest.mark.parametrize("shape", CASES)
def test_radial_integral_matches_arb_with_tau_constraint(
    shape: str, mode: str, truth_table: dict, cosmo_bits: tuple
) -> None:
    """I_ell(k) against Arb with the tau constraint, one multipole at a time."""
    cosmo, dist, ps = cosmo_bits
    entry = truth_table["shapes"][shape]

    kernel = _make_kernel(
        _shape_of(entry, shape), dist, ps, Ncm.SBesselIntegratorLevin.new(0, 8), entry["ctor"]
    )
    kernel.set_l_limber(-1)
    kernel.prepare(cosmo)

    supports = [kernel.get_comp_support(comp) for comp in range(kernel.get_n_comps())]
    fallbacks = 0

    for index, ell in enumerate(truth_table["ells"]):
        expected = np.array([float(value) for value in entry["table"][index]])
        peak = np.abs(expected).max()

        integrator = Ncm.SBesselIntegratorLevin.new(ell, ell)
        _set_constraint(integrator, mode)

        got = np.array(
            [
                sum(
                    integrator.integrate_ell(
                        lambda _ud, chi, _k, c=comp, lo=low, hi=high: (
                            kernel.eval_W_comp(c, min(max(chi, lo), hi))
                        ),
                        low,
                        high,
                        k,
                        ell,
                        None,
                    )
                    for comp, (low, high) in enumerate(supports)
                )
                for k in entry["kvals"][index]
            ]
        )
        fallbacks += integrator.get_n_constraint_fallbacks()

        assert_allclose(
            got,
            expected,
            rtol=RTOL,
            atol=ATOL_FRAC * peak,
            err_msg=f"{shape} at ell = {ell}, tau constraint ({mode})",
        )

    if mode == "forced":
        # The table samples below the turning point for every shape; without the
        # guard those entries were wrong by 1e10 and more.
        assert fallbacks > 0, "the guard never fired on a table that requires it"


@pytest.mark.parametrize("mode", TAU_CONSTRAINT_MODES)
@pytest.mark.parametrize("shape", CASES)
def test_radial_integral_batched_matches_arb_with_tau_constraint(
    shape: str, mode: str, truth_table: dict, cosmo_bits: tuple
) -> None:
    """Same comparison through the batched path, blocks of eight multipoles."""
    cosmo, dist, ps = cosmo_bits
    entry = truth_table["shapes"][shape]
    n_block = 8

    kernel = _make_kernel(
        _shape_of(entry, shape), dist, ps, Ncm.SBesselIntegratorLevin.new(0, 8), entry["ctor"]
    )
    kernel.set_l_limber(-1)
    kernel.prepare(cosmo)

    supports = [kernel.get_comp_support(comp) for comp in range(kernel.get_n_comps())]
    fallbacks = 0

    for index, ell in enumerate(truth_table["ells"]):
        expected = np.array([float(value) for value in entry["table"][index]])
        peak = np.abs(expected).max()

        integrator = Ncm.SBesselIntegratorLevin.new(ell, ell + n_block - 1)
        _set_constraint(integrator, mode)
        block = Ncm.Vector.new(n_block)
        got = []

        for k in entry["kvals"][index]:
            total = 0.0

            for comp, (low, high) in enumerate(supports):
                integrator.integrate(
                    lambda chi, _k, c=comp, lo=low, hi=high: (
                        kernel.eval_W_comp(c, min(max(chi, lo), hi))
                    ),
                    low,
                    high,
                    k,
                    block,
                )
                total += block.get(0)

            got.append(total)

        fallbacks += integrator.get_n_constraint_fallbacks()

        assert_allclose(
            np.array(got),
            expected,
            rtol=RTOL,
            atol=ATOL_FRAC * peak,
            err_msg=f"{shape} at ell = {ell}, batched tau constraint ({mode})",
        )

    if mode == "forced":
        assert fallbacks > 0, "the guard never fired on a table that requires it"


# Worst deviation measured per shape, closure at reltol = peak-epsilon = 1e-6,
# with roughly a factor of three of headroom. Re-measure these after changing
# NC_XCOR_KERNEL_CHEB_PANEL_K_CAP: a smaller cap makes more, lower-order panels,
# which converge to the requested tolerance by a different route.
#
# Read the spread against the spline's, not on its own. Held to this same allowance
# the Chebyshev closure sits at 0.33x of it on every shape, which is the headroom
# above and not a marginal pass, while the spline needs 1.2x, 0.2x, 11x, 32x, 1.8x,
# 7.4x and 264x of it: it misses on five of the seven shapes, and beats the Chebyshev
# closure only on the plain tophat.
# Where a number here is large it is the *sampling* that binds, not the fit:
# per multipole the closure sits on the sampling floor wherever the convergence
# criterion lets it reach, which a worst-over-ell figure like this cannot show.
CLOSURE_TOL = {
    "gauss_mid": 9.0e-4,
    "tophat_mid": 7.0e-7,
    "student_t_mid": 2.0e-8,
    "power_exp_mid": 2.0e-6,
    "lensing_gal": 7.0e-6,
    "multi_mid": 4.0e-4,
    "tophat_smooth_mid": 4.0e-8,
}

# The absolute half of the same criterion, as a fraction of the block's peak, and the
# reason no point has to be excluded. A relative bound alone is unreachable wherever
# I_ell is small, which #NcXcorKernel:peak-epsilon documents: its criterion is
# absolute, so "the corresponding relative error can become large where the resulting
# C_l is small", and the useful precision is capped by cancellation in the radial
# integral that produces W_i(k). Measured, that is exactly what happens -- on
# power_exp at ell = 50 the point at 1.08e-3 of peak reads 1.7e-2 relative and does
# not move between reltol 1e-6 and 1e-9 while its neighbours reach 1e-13. An error
# that no tolerance moves belongs in an atol, not in a widened rtol or a discarded
# sample. Measured per shape with roughly a factor of three of headroom.
CLOSURE_ATOL = {
    "gauss_mid": 1.9e-6,
    "tophat_mid": 1.6e-5,
    "student_t_mid": 7.9e-6,
    "power_exp_mid": 4.9e-6,
    "lensing_gal": 1.8e-6,
    "multi_mid": 6.3e-13,
    "tophat_smooth_mid": 5.5e-8,
}


@pytest.mark.parametrize("shape", sorted(CLOSURE_TOL))
def test_chebyshev_closure_matches_arb(
    shape: str, truth_table: dict, cosmo_bits: tuple
) -> None:
    """The closure itself against Arb, not just the radial integral under it.

    A closure holds C * sqrt(P(k)) * I_ell(k) for a constant C, so dividing the
    certified I_ell and sqrt(P) out has to leave a constant. How far that ratio
    wanders across k is the closure's own fitting error, measured against values
    with proven radii rather than against another of NumCosmo's paths.

    """
    cosmo, dist, ps = cosmo_bits
    entry = truth_table["shapes"][shape]
    RH = Nc.HICosmo.RH_Mpc(cosmo)

    worst = 0.0
    worst_at = ""
    compared = 0

    for index, ell in enumerate(truth_table["ells"]):
        kernel = _make_kernel(
            _shape_of(entry, shape), dist, ps, Ncm.SBesselIntegratorLevin.new(0, 8), entry["ctor"]
        )
        kernel.set_l_limber(-1)
        kernel.set_property("reltol", 1.0e-6)
        kernel.set_property("peak-epsilon", 1.0e-6)
        kernel.prepare(cosmo)

        integrand = kernel.get_eval_vectorized_full(
            cosmo,
            ell,
            ell,
            Ncm.SBesselIntegratorLevin.new(ell, ell),
            Nc.XcorKernelClosure.CHEBYSHEV,
        )
        k_min, k_max = integrand.get_range()

        expected = np.array([float(value) for value in entry["table"][index]])

        # Outside the fitted domain the closure extrapolates, so only the fitted
        # range is comparable. Every point inside it is kept: the criterion below
        # carries an absolute term, so a sample where I_ell is a thousandth of its
        # peak no longer has to clear a relative bound that the closure's own
        # sampling floor puts out of reach.
        kept = [
            (k, value)
            for k, value in zip(entry["kvals"][index], expected)
            if k_min < k * RH < k_max
        ]

        if len(kept) < 3:
            continue

        reference = np.array(
            [value * np.sqrt(ps.eval(cosmo, 0.0, k)) for k, value in kept]
        )
        got = np.array([integrand.eval_array(k * RH)[0] for k, _ in kept])
        # The closure carries an overall constant, so the comparison is up to one.
        held = np.median(got / reference) * reference
        allowed = (
            CLOSURE_TOL[shape] * np.abs(held) + CLOSURE_ATOL[shape] * np.abs(held).max()
        )
        excess = (np.abs(got - held) / allowed).max()

        if excess > worst:
            worst, worst_at = excess, f"ell = {ell}"

        compared += len(kept)

    assert compared > 0, f"{shape} shared no k with the closure's fitted range"
    assert worst < 1.0, (
        f"{shape}: closure exceeds rtol {CLOSURE_TOL[shape]:.1e} plus atol "
        f"{CLOSURE_ATOL[shape]:.1e} of peak by {worst:.3f}x at {worst_at}, "
        f"over {compared} points"
    )

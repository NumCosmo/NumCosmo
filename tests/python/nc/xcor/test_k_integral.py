#!/usr/bin/env python
#
# test_k_integral.py
#
# Fri Aug 28 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_k_integral.py
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

"""The outer k-integral: each kernel-space method against its own truth.

Every method is handed frozen closures and compared to a reference built from
*the closures that method actually integrates* -- per ell block for
KERNEL_EXACT and KERNEL_CUBATURE, per multipole for KERNEL_GSL. Comparing all
three against one block-closure reference charges GSL for a closure it never
built, and on a strongly cancelling pair that difference is a factor of three,
not a rounding.

What is measured here is therefore the quadrature alone. The closure's own
error against certified Arb values is ``test_xcor_window_truth_table.py``'s
subject, and the two together are what ``bench_k_integral.py`` reports.

The case matrix lives in ``cases_k_integral.py``. Every tolerance below was
read off that bench driver; none was chosen.

Errors are taken against the *block's peak*, not pointwise. On a far-separated
pair the smallest multipole is five orders below the largest, and a relative
tolerance there is a much harder request than the same number on an auto
spectrum -- one no method meets and none needs to, because that entry
contributes nothing to any likelihood.
"""

import collections
import gzip
import json
import os
import pathlib

import numpy as np
import pytest

from numcosmo_py import Nc, Ncm
from xcor import cases_k_integral as cases

pytest_plugins = ["python.fixtures_xcor"]

# xdist_group pins this whole file to one worker under --dist loadgroup.
# The frozen fixture below caches up to 102 Frozen objects, each holding two
# kernels, two closures and their references -- 3.5 GB by the end of the file.
# That is per *worker*, since an xdist worker is its own session, so plain
# --dist load builds one copy per worker and exhausts memory on a many-core
# machine (measured: OOM at 12 workers, and every run at 24). Grouping costs
# this file its internal parallelism and nothing else -- every other file still
# distributes test by test.
pytestmark = [pytest.mark.xcor, pytest.mark.xdist_group("k_integral")]

CLOSURES = {
    "spline": Nc.XcorKernelClosure.SPLINE,
    "chebyshev": Nc.XcorKernelClosure.CHEBYSHEV,
}

METHODS = {
    "exact": Nc.XcorMethod.KERNEL_EXACT,
    "cubature": Nc.XcorMethod.KERNEL_CUBATURE,
    "gsl": Nc.XcorMethod.KERNEL_GSL,
    "gsl_block": Nc.XcorMethod.KERNEL_GSL_BLOCK,
}

# KERNEL_GSL calls nc_xcor_kernel_get_eval once per multipole; every other
# method calls ..._get_eval_vectorized_full once per block. KERNEL_GSL_BLOCK
# runs KERNEL_GSL's rule over KERNEL_EXACT's closure, which is what makes a
# quadrature comparison a comparison of quadratures.
PER_MULTIPOLE = frozenset({"gsl"})

# The methods reachable through Nc.Xcor.integrate_block(), i.e. the ones with
# an entry in the NcXcorKQuad table.
BLOCK_METHODS = sorted(set(METHODS) - PER_MULTIPOLE)

# Measured over the case matrix at reltol = scaled_abstol = 1e-4, as the worst
# deviation from the matching reference relative to the block's peak, with an
# order of headroom. Measured worst, on cases.ELLS_SUITE:
#
#   exact    spline     8.1e-12       cubature spline     1.7e-05
#   exact    chebyshev  2.3e-08       cubature chebyshev  1.1e-05
#   gsl      spline     4.4e-12       gsl      chebyshev  1.6e-03
#   gsl_blk  spline     3.2e-11       gsl_blk  chebyshev  3.2e-04
#
# These are an order or two above what a [2, 20, 200] ladder reported, because
# that ladder missed the l = 4-10 region where every method is worst; see
# cases.ELLS_SUITE. Read the two exact rows as the claim they are: GL(5) on the
# merged knot set is exact, and so is qagp on those same knots -- both sit at
# the reference's own floor. Cubature is not converging to the closure, it is
# stopping at the relative tolerance it was asked for.
#
# gsl/chebyshev at 1.6e-03 is the largest error anywhere in the matrix, at
# l = 10. It is the one combination where the per-multipole closure and the
# Chebyshev fit interact badly, and it is why that tolerance is so much looser
# than its spline counterpart.
#
# gsl_block runs gsl's rule on exact's closure and lands within an order of
# both: the rule is not what separates them. On splines it is the looser of the
# two by about eight times, and that is the block closure rather than the
# quadrature -- a block is fitted to an L2 norm over all its multipoles, so a
# multipole that is sub-dominant within its block is held only to the block's
# norm, while gsl fits each one on its own. See nc_xcor_compute_full().
TOLERANCE = {
    ("exact", "spline"): 1.0e-10,
    ("exact", "chebyshev"): 5.0e-7,
    ("cubature", "spline"): 2.0e-4,
    ("cubature", "chebyshev"): 2.0e-4,
    ("gsl", "spline"): 5.0e-11,
    ("gsl", "chebyshev"): 2.0e-2,
    ("gsl_block", "spline"): 5.0e-10,
    ("gsl_block", "chebyshev"): 5.0e-3,
}

# The reference's own convergence, likewise measured on cases.ELLS_SUITE: the
# worst any cell still moved at the top of the escalation, 9.0e-13 on spline
# closures and 1.3e-08 on Chebyshev ones, both at l = 6. On Chebyshev the floor
# is the conditioning of evaluating a degree-128 polynomial rather than the
# quadrature, and l = 6 is where that bites -- 130x worse than the 1.0e-10 at
# l = 2 that a [2, 20, 200] ladder saw. It bounds what any test here can claim
# about a Chebyshev closure at that multipole: exact/chebyshev's own 2.3e-08
# sits less than a factor of two above it, so that entry bounds the method
# rather than discriminating between methods.
REFERENCE_FLOOR = {"spline": 1.0e-11, "chebyshev": 1.0e-7}


class Frozen:
    """Two prepared kernels, their closures for one block, and the references."""

    def __init__(self, case: str, closure: str, lmin: int) -> None:
        self.pair = cases.PAIRS_BY_CASE[case]
        self.settings = cases.Settings(closure=CLOSURES[closure])
        self.lmin = lmin
        self.lmax = lmin + self.settings.ell_batch_size - 1
        self.cosmo, self.dist, self.ps = cases.make_cosmo_bits()
        self.RH = Nc.HICosmo.RH_Mpc(self.cosmo)

        # One integrator for the whole block, shared by both kernels and both
        # closures -- what nc_xcor_solver_solve() does, and what makes the
        # stored decomposition reusable across kernels. Giving each kernel its
        # own retains four integrators per key and never exercises that path.
        # Scoped to this Frozen rather than to the block: reuse is order
        # sensitive, and pytest-randomly means a wider scope would make cost,
        # and any non-exactness, depend on execution order.
        self.sbi = Ncm.SBesselIntegratorLevin.new(lmin, self.lmax)

        self.kernel_a = cases.build_kernel(
            self.pair.kernel_a, self.cosmo, self.dist, self.ps, self.settings, self.sbi
        )
        self.kernel_b = (
            self.kernel_a
            if self.pair.isauto
            else cases.build_kernel(
                self.pair.kernel_b,
                self.cosmo,
                self.dist,
                self.ps,
                self.settings,
                self.sbi,
            )
        )
        self.integrand_a = cases.build_integrand(
            self.kernel_a, self.cosmo, lmin, self.lmax, self.settings, self.sbi
        )
        self.integrand_b = (
            None
            if self.pair.isauto
            else cases.build_integrand(
                self.kernel_b, self.cosmo, lmin, self.lmax, self.settings, self.sbi
            )
        )
        self.reference = cases.reference_cl(self.RH, self.integrand_a, self.integrand_b)
        self.cancellation = cases.cancellation_ratio(self.integrand_a, self.integrand_b)
        self._per_multipole = None

    @property
    def per_multipole(self) -> np.ndarray:
        """The same integral from one closure per multipole, built on demand."""
        if self._per_multipole is None:
            self._per_multipole = cases.per_multipole_reference(
                self.RH,
                self.kernel_a,
                None if self.pair.isauto else self.kernel_b,
                self.cosmo,
                self.lmin,
                self.lmax,
                self.settings,
            )

        return self._per_multipole

    def compute(self, method: str) -> np.ndarray:
        """One block through the library, by the named kernel-space method."""
        xcor = Nc.Xcor.new(self.dist, self.ps, METHODS[method])
        xcor.set_closure_type(self.settings.closure)
        xcor.set_ell_batch_size(self.settings.ell_batch_size)
        xcor.set_reltol(self.settings.reltol)
        xcor.prepare(self.cosmo)

        vp = Ncm.Vector.new(self.lmax - self.lmin + 1)
        xcor.compute(self.kernel_a, self.kernel_b, self.cosmo, self.lmin, self.lmax, vp)

        return np.array(vp.dup_array())

    def integrate_block(self, method: str) -> np.ndarray:
        """The same block, driven straight off the closures built above.

        No kernel, no closure build, no batching: this is the outer quadrature
        alone, which is the only way to compare two methods on one integrand
        rather than on two separately fitted ones.
        """
        xcor = Nc.Xcor.new(self.dist, self.ps, METHODS[method])
        xcor.set_closure_type(self.settings.closure)
        xcor.set_ell_batch_size(self.settings.ell_batch_size)
        xcor.set_reltol(self.settings.reltol)
        xcor.prepare(self.cosmo)

        vp = Ncm.Vector.new(self.lmax - self.lmin + 1)
        xcor.integrate_block(
            self.integrand_a,
            self.integrand_b,
            self.lmin,
            self.lmax,
            self.pair.isauto,
            METHODS[method],
            vp,
            None,
        )

        return np.array(vp.dup_array())

    def truth_for(self, method: str) -> np.ndarray:
        """The reference built on the closure this method actually integrates."""
        return self.per_multipole if method in PER_MULTIPOLE else self.reference.cl

    def peak_error(self, got: np.ndarray, truth: np.ndarray) -> float:
        """Worst deviation over the block, against the block's own peak."""
        return float(np.abs(got - truth).max() / np.abs(truth).max())


# How many Frozen the module fixture keeps alive at once. Each holds two
# kernels, two closures and their references -- 18 MB for an easy key, 105 MB
# for the worst (ell = 200). Unbounded, the 102 keys of the matrix reach 3.5 GB,
# and since an xdist worker is its own session that is paid per worker: enough
# to exhaust a 30 GB machine at 12 workers and every time at 24. Bounded, the
# cost is a rebuild whenever an evicted key comes back, at 0.194 s each.
# 0 means unbounded.
FROZEN_CACHE_MAXSIZE = int(os.environ.get("XCOR_FROZEN_CACHE", "8"))


@pytest.fixture(name="frozen", scope="module")
def fixture_frozen():
    """Build each (case, closure, block) once and hand it to every test.

    Least-recently-used beyond FROZEN_CACHE_MAXSIZE, because holding all of
    them is what makes this file the memory ceiling of the whole xcor lane.
    """
    cache: collections.OrderedDict[tuple[str, str, int], Frozen] = (
        collections.OrderedDict()
    )

    def get(case: str, closure: str, lmin: int) -> Frozen:
        key = (case, closure, lmin)

        if key in cache:
            cache.move_to_end(key)

            return cache[key]

        cache[key] = Frozen(case, closure, lmin)

        while FROZEN_CACHE_MAXSIZE and len(cache) > FROZEN_CACHE_MAXSIZE:
            cache.popitem(last=False)

        return cache[key]

    return get


@pytest.mark.parametrize("case", [pair.case for pair in cases.PAIRS])
@pytest.mark.parametrize("closure", sorted(CLOSURES))
@pytest.mark.parametrize("lmin", cases.ELLS_SUITE)
def test_reference_is_a_reference(frozen, case: str, closure: str, lmin: int) -> None:
    """The reference's own convergence, before anything is compared to it.

    Raising the Gauss-Legendre order on a cell has to stop moving the answer,
    or the comparisons below measure the reference rather than the method. A
    spline cell is a cubic and settles at the first order; a Chebyshev cell
    carries up to 129 coefficients and settles when the order passes its
    degree.
    """
    reference = frozen(case, closure, lmin).reference

    assert reference.n_cells > 0
    assert reference.worst_cell_move < REFERENCE_FLOOR[closure]


@pytest.mark.parametrize("case", [pair.case for pair in cases.PAIRS])
@pytest.mark.parametrize("closure", sorted(CLOSURES))
@pytest.mark.parametrize("method", sorted(METHODS))
@pytest.mark.parametrize("lmin", cases.ELLS_SUITE)
def test_method_matches_reference(
    frozen, case: str, closure: str, method: str, lmin: int
) -> None:
    """Every kernel-space method against the quadrature's own truth."""
    state = frozen(case, closure, lmin)
    error = state.peak_error(state.compute(method), state.truth_for(method))
    tolerance = TOLERANCE[(method, closure)]

    assert error < tolerance, (
        f"{case} ({state.pair.regime}): {method} on a {closure} closure "
        f"deviates by {error:.3e} of the block's peak at ell "
        f"{state.lmin}-{state.lmax}, above the measured {tolerance:.1e}"
    )


@pytest.mark.parametrize("case", [pair.case for pair in cases.PAIRS])
@pytest.mark.parametrize("closure", sorted(CLOSURES))
@pytest.mark.parametrize("method", BLOCK_METHODS)
@pytest.mark.parametrize("lmin", cases.ELLS_SUITE)
def test_block_entry_point_matches_compute(
    frozen, case: str, closure: str, method: str, lmin: int
) -> None:
    """integrate_block() on shared closures is what compute() runs internally.

    Exact equality, not a tolerance: compute() builds the closures and calls
    the same NcXcorKQuad entry this does, so any difference is the two paths
    disagreeing about which quadrature or which range they are on. That makes
    this the sharpest check on the table -- it removes the closure rebuild
    that every other comparison has to allow for.
    """
    state = frozen(case, closure, lmin)
    direct = state.integrate_block(method)
    through_compute = state.compute(method)

    np.testing.assert_array_equal(
        direct,
        through_compute,
        err_msg=(
            f"{case}: {method} on a {closure} closure disagrees with itself "
            f"between Nc.Xcor.integrate_block() and nc_xcor_compute() at ell "
            f"{state.lmin}-{state.lmax}"
        ),
    )


@pytest.mark.parametrize("method", sorted(METHODS))
def test_method_introspection_is_consistent(method: str) -> None:
    """What the table says about a method matches what it does."""
    meth = METHODS[method]

    assert Nc.xcor_method_is_kernel_space(meth)
    assert Nc.xcor_method_get_name(meth).startswith("NC_XCOR_METHOD_KERNEL_")

    # Only KERNEL_EXACT reports an error, and it is the only one whose block
    # quadrature carries no tolerance of its own to report instead.
    assert Nc.xcor_method_has_error_estimate(meth) == (method == "exact")


def test_limber_z_methods_have_no_block_quadrature() -> None:
    """The redshift-space tier is a different approximation, not a rule.

    integrate_block() refuses it, and refuses KERNEL_GSL too -- that one is
    kernel-space but fits its closure one multipole at a time, so there is no
    block for a caller to hand it.
    """
    for meth in (
        Nc.XcorMethod.LIMBER_Z_GSL,
        Nc.XcorMethod.LIMBER_Z_CUBATURE,
    ):
        assert not Nc.xcor_method_is_kernel_space(meth)
        assert not Nc.xcor_method_has_error_estimate(meth)

    assert Nc.xcor_method_is_kernel_space(Nc.XcorMethod.KERNEL_GSL)
    assert not Nc.xcor_method_has_error_estimate(Nc.XcorMethod.KERNEL_GSL)


@pytest.mark.parametrize("case", [pair.case for pair in cases.PAIRS if pair.isauto])
@pytest.mark.parametrize("closure", sorted(CLOSURES))
def test_auto_spectra_cannot_cancel(frozen, case: str, closure: str) -> None:
    """W^2 >= 0, so an auto spectrum's integrand has nothing to cancel.

    This is why the matrix is not all auto spectra: a method's per-node error
    reaches the answer undivided here, and by up to seven orders more on a
    far-separated cross pair. A case set that lost its cross pairs would still
    pass every tolerance above and measure nothing.
    """
    assert frozen(case, closure, 2).cancellation == pytest.approx(1.0, abs=1.0e-12)


@pytest.mark.parametrize("closure", sorted(CLOSURES))
def test_far_separated_bins_cancel_by_orders(frozen, closure: str) -> None:
    """And the cross pair the matrix exists for still does cancel.

    Tail against tail at ell = 2 is where a quadrature's per-node error is
    amplified the most; by ell = 200 each window is confined near its own
    turning point and the overlap is genuinely small rather than cancelling.
    """
    low = frozen("X3", closure, 2).cancellation[0]

    assert low > 1.0e3


def test_narrow_shell_caps_the_spline_closure() -> None:
    """A 56 Mpc hard shell: the case the spectral closure exists for.

    The knob that matters is ``scaled-abstol``, not ``reltol``. Measured on
    this shell, four decades of ``reltol`` (1e-4 to 1e-8) change the answer by
    nothing at all for either closure and cost nothing either; every decade of
    ``scaled-abstol`` moves both. So the comparison is taken along that axis,
    at the library's own documented floor of 1e-6 -- below it the setter warns,
    because the value is measured against the peak of W(k) while the C_ell
    integrand is k^2 W_a W_b, so it enters squared.

    At that floor the spline closure sits at ~2e-5 and has nowhere left to go;
    the spectral closure is at machine zero. That gap is the justification for
    carrying two closure types, so it is asserted rather than remembered.
    """
    cosmo, dist, ps = cases.make_cosmo_bits()
    RH = Nc.HICosmo.RH_Mpc(cosmo)
    pair = cases.PAIRS_BY_CASE["N1"]

    def solve(closure: Nc.XcorKernelClosure, scaled_abstol: float) -> np.ndarray:
        settings = cases.Settings(
            reltol=1.0e-6, scaled_abstol=scaled_abstol, closure=closure
        )
        kernel = cases.build_kernel(pair.kernel_a, cosmo, dist, ps, settings)

        return cases.reference_cl(
            RH, cases.build_integrand(kernel, cosmo, 2, 9, settings), None
        ).cl

    truth = solve(Nc.XcorKernelClosure.CHEBYSHEV, 1.0e-6)
    peak = np.abs(truth).max()

    def error(values: np.ndarray) -> float:
        return float(np.abs(values - truth).max() / peak)

    # The spline closure at the floor, and one decade above it. It improves
    # with scaled-abstol, so this is not a stall -- it is a cap, and the floor
    # is where the cap bites.
    spline_floor = error(solve(Nc.XcorKernelClosure.SPLINE, 1.0e-6))
    spline_above = error(solve(Nc.XcorKernelClosure.SPLINE, 1.0e-5))

    assert spline_floor > 1.0e-6
    assert spline_above > spline_floor

    # The spectral closure, one decade *above* the floor, is already ahead.
    assert error(solve(Nc.XcorKernelClosure.CHEBYSHEV, 1.0e-5)) < spline_floor


def test_reltol_is_not_what_moves_the_narrow_shell() -> None:
    """And the knob that looks like the accuracy knob is inert here.

    Guards the reading above: if a future change makes ``reltol`` matter on
    this shell, the comparison in ``test_narrow_shell_caps_the_spline_closure``
    is taken along the wrong axis and has to be revisited.
    """
    cosmo, dist, ps = cases.make_cosmo_bits()
    RH = Nc.HICosmo.RH_Mpc(cosmo)
    pair = cases.PAIRS_BY_CASE["N1"]

    def solve(reltol: float) -> np.ndarray:
        settings = cases.Settings(
            reltol=reltol,
            scaled_abstol=1.0e-5,
            closure=Nc.XcorKernelClosure.SPLINE,
        )
        kernel = cases.build_kernel(pair.kernel_a, cosmo, dist, ps, settings)

        return cases.reference_cl(
            RH, cases.build_integrand(kernel, cosmo, 2, 9, settings), None
        ).cl

    loose = solve(1.0e-4)
    tight = solve(1.0e-8)

    assert np.abs(tight - loose).max() / np.abs(loose).max() < 1.0e-2


# ---------------------------------------------------------------- Level 2 ---
#
# Against certified Arb values rather than against a reference of our own. The
# table is data/truth_tables/xcor/xcor_kquad.json.gz, generated offline by
# tests/tools/nc_xcor_kquad_arb.c; only regeneration needs FLINT.
#
# Deviations are measured against each pair's own scale -- its largest
# certified |C_ell| over the multipoles in the table -- not pointwise. A
# far-separated pair reaches |C_ell| = 7e-20 at ell = 50, fifteen orders below
# the auto-spectrum scale, and a pointwise relative error there reads 2e+09
# while meaning nothing.

TRUTH_TABLE = "truth_tables/xcor/xcor_kquad.json.gz"

# Measured over the 43 certified entries, with headroom. These are loose
# because one regime is genuinely bad and the table says so rather than hiding
# it: on the far-separated pair at ell = 2 the spline closure is wrong by 435%
# of the pair's own scale. The chebyshev worst is that same pair at ell = 50,
# 1.04e-2 of scale, on a certified value eleven orders below the pair's peak.
# The tight statement is the comparative test below.
#
# A loose gate has a cost, paid once already: the X9 entries were certified
# with the generator's side-b leak (kdep applied to both windows), putting
# the reference 5.3% off at ell = 50, and both closures' agreement with each
# other while sitting 2.9e-2 from the table stayed under the old chebyshev
# gate of 1e-1. A tolerance-independent, closure-independent deviation is a
# wrong reference, not a bad fit -- worth a ladder check before widening
# either number here.
CERTIFIED_TOLERANCE = {"spline": 1.5e1, "chebyshev": 3.0e-2}


@pytest.fixture(name="certified", scope="module")
def fixture_certified() -> dict:
    """The certified C_ell table, loaded once."""
    path = pathlib.Path(Ncm.cfg_get_data_filename(TRUTH_TABLE, True))

    with gzip.open(path, "rt") as handle:
        return json.load(handle)


def _pair_scale(table: dict) -> dict:
    """Each pair's largest certified |C_ell| over the multipoles present."""
    scale: dict[str, float] = {}

    for entry in table["cases"].values():
        value = abs(float(entry["value"]))
        scale[entry["case"]] = max(scale.get(entry["case"], 0.0), value)

    return scale


def _library_cl(case: str, ell: int, closure: str) -> float:
    """What the library returns for one entry, by its exact method."""
    pair = cases.PAIRS_BY_CASE[case]
    settings = cases.Settings(closure=CLOSURES[closure])
    cosmo, dist, ps = cases.make_cosmo_bits()

    kernel_a = cases.build_kernel(pair.kernel_a, cosmo, dist, ps, settings)
    kernel_b = (
        kernel_a
        if pair.isauto
        else cases.build_kernel(pair.kernel_b, cosmo, dist, ps, settings)
    )

    xcor = Nc.Xcor.new(dist, ps, Nc.XcorMethod.KERNEL_EXACT)
    xcor.set_closure_type(settings.closure)
    xcor.set_ell_batch_size(1)
    xcor.set_reltol(settings.reltol)
    xcor.prepare(cosmo)

    vp = Ncm.Vector.new(1)
    xcor.compute(kernel_a, kernel_b, cosmo, ell, ell, vp)

    return vp.get(0)


def test_certified_table_is_certified(certified: dict) -> None:
    """The reference has to be a reference: its own radius cannot matter.

    Guards a regeneration that quietly lowered the precision target -- every
    assertion below would still pass while checking nothing.
    """
    worst = max(
        entry["radius"] / abs(float(entry["value"]))
        for entry in certified["cases"].values()
    )

    assert certified["cases"]
    assert worst < 1.0e-8


@pytest.mark.parametrize("closure", sorted(CLOSURES))
def test_closure_matches_certified_c_ell(certified: dict, closure: str) -> None:
    """Both closures against Arb, over every entry the table certifies."""
    scale = _pair_scale(certified)
    worst, worst_at = 0.0, ""

    for entry in certified["cases"].values():
        got = _library_cl(entry["case"], entry["ell"], closure)
        deviation = abs(got - float(entry["value"])) / scale[entry["case"]]

        if deviation > worst:
            worst, worst_at = deviation, f"{entry['case']} ell={entry['ell']}"

    assert worst < CERTIFIED_TOLERANCE[closure], (
        f"{closure} deviates from the certified C_ell by {worst:.3e} of the "
        f"pair's scale at {worst_at}"
    )


def test_spectral_closure_is_closer_to_certified_truth(certified: dict) -> None:
    """The reason two closure types exist, measured against proven values.

    Level 1 cannot state this: it compares each method to a reference built on
    the closure under test, so a closure that is wrong in the same way as its
    reference looks right. Only a certified value separates them.
    """
    scale = _pair_scale(certified)
    spline, chebyshev = [], []

    for entry in certified["cases"].values():
        truth = float(entry["value"])
        norm = scale[entry["case"]]
        spline.append(
            abs(_library_cl(entry["case"], entry["ell"], "spline") - truth) / norm
        )
        chebyshev.append(
            abs(_library_cl(entry["case"], entry["ell"], "chebyshev") - truth) / norm
        )

    closer = sum(1 for s, c in zip(spline, chebyshev) if c < s)

    # Measured: closer in 36 of 43, median 4.0e-6 against 1.8e-5. (An earlier
    # 34-of-43 was taken against X9 entries certified with the generator's
    # side-b leak; the library was right and the table was wrong there.)
    assert closer >= 0.7 * len(spline)
    assert np.median(chebyshev) < 0.5 * np.median(spline)

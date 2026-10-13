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

KERNEL_EXACT is handed frozen closures and compared to a reference built from
the closures it integrates, one per ell block.

What is measured here is therefore the quadrature alone. The closure's own
error against certified Arb values is the subject of the C test
``tests/c/nc/xcor/test_nc_xcor_window_dense.c``, and the two together are what
``bench_k_integral.py`` reports.

The case matrix lives in ``cases_k_integral.py``. Every tolerance below was
read off that bench driver; none was chosen.

Errors are taken against the *block's peak*, not pointwise. On a far-separated
pair the smallest multipole is five orders below the largest, and a relative
tolerance there is a much harder request than the same number on an auto
spectrum -- one no method meets and none needs to, because that entry
contributes nothing to any likelihood.
"""

from typing import Callable
import collections
import gzip
import json
import os
import pathlib

import numpy as np
import pytest

from xcor import cases_k_integral as cases
from numcosmo_py import Nc, Ncm

pytest_plugins = ["python.fixtures_xcor"]

# Under --dist loadgroup every test that uses the frozen fixture is pinned to one
# worker through FROZEN_GROUP below; the module-level group holds the rest. An
# xdist worker is its own session, so the module-scoped frozen cache exists once
# per worker: unbounded it reaches 3.5 GB (measured OOM at 12 workers under plain
# --dist load), bounded by FROZEN_CACHE_MAXSIZE it is under 1 GB.
pytestmark = [pytest.mark.xcor, pytest.mark.xdist_group("k_integral")]

# The worker group of every test that takes `frozen`; a function mark is closer than
# the module's, so it is the one xdist reads.
FROZEN_GROUP = pytest.mark.xdist_group("k_integral_frozen")

METHODS = {
    "exact": Nc.XcorMethod.KERNEL_EXACT,
}

# Measured over the case matrix at reltol = peak_epsilon = 1e-4, as the worst
# deviation from the matching reference relative to the block's peak, with an
# order of headroom. Measured worst, on cases.ELLS_SUITE: exact/chebyshev
# 2.3e-08. That is an order or two above what a [2, 20, 200] ladder reported,
# because that ladder missed the l = 4-10 region where every method is worst; see
# cases.ELLS_SUITE.
TOLERANCE = {
    "exact": 5.0e-7,
}

# The reference's own convergence, likewise measured on cases.ELLS_SUITE: the
# worst any cell still moved at the top of the escalation, 1.3e-08, at l = 6. The
# floor is the conditioning of evaluating a degree-128 polynomial rather than the
# quadrature, and l = 6 is where that bites -- 130x worse than the 1.0e-10 at
# l = 2 that a [2, 20, 200] ladder saw. It bounds what any test here can claim
# about a Chebyshev closure at that multipole: exact/chebyshev's own 2.3e-08
# sits less than a factor of two above it, so that entry bounds the method
# rather than discriminating between methods.
REFERENCE_FLOOR = 1.0e-7

# The pairs and multipoles the quadrature is checked on; bench_k_integral.py runs
# the whole case matrix. One pair stands for each kind of integrand the outer
# quadrature meets: Gaussian and top-hat autos, an overlapping and a far-separated
# cross, the RSD weight, a lensing kernel, a disconnected support, and X9, the only
# non-separable kernel, whose closure is built differently. The multipoles are the
# low end, l = 6, where every method and the reference are worst (the l = 4-10
# region a [2, 20, 200] ladder missed), and l = 50 above it.
CHECKED_CASES = ("A1", "A2", "X1", "X10", "R2", "X6", "A5", "X9")
CHECKED_ELLS = (2, 6, 50)


class Frozen:
    """Hold two prepared kernels, their closures for one block, and the references."""

    def __init__(self, case: str, lmin: int) -> None:
        """Build the kernels and closures for one (case, block) once."""
        self.pair = cases.PAIRS_BY_CASE[case]
        self.settings = cases.Settings()
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

    def compute(self, method: str) -> np.ndarray:
        """Compute one block through the library, by the named kernel-space method."""
        xcor = Nc.Xcor.new(self.dist, self.ps, METHODS[method])
        xcor.set_ell_batch_size(self.settings.ell_batch_size)
        xcor.set_reltol(self.settings.reltol)
        xcor.prepare(self.cosmo)

        vp = Ncm.Vector.new(self.lmax - self.lmin + 1)
        xcor.compute(self.kernel_a, self.kernel_b, self.cosmo, self.lmin, self.lmax, vp)

        return np.array(vp.dup_array())

    def integrate_block(self, method: str) -> np.ndarray:
        """Integrate one block straight off the closures built above.

        No kernel, no closure build, no batching: this is the outer quadrature
        alone, which is the only way to compare two methods on one integrand
        rather than on two separately fitted ones.
        """
        xcor = Nc.Xcor.new(self.dist, self.ps, METHODS[method])
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

    def peak_error(self, got: np.ndarray, truth: np.ndarray) -> float:
        """Return the worst deviation over the block, relative to the block's peak."""
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
def fixture_frozen() -> Callable[[str, int], Frozen]:
    """Build each (case, block) once and hand it to every test.

    Least-recently-used beyond FROZEN_CACHE_MAXSIZE, because holding all of
    them is what makes this file the memory ceiling of the whole xcor lane.
    """
    cache: collections.OrderedDict[tuple[str, int], Frozen] = collections.OrderedDict()

    def get(case: str, lmin: int) -> Frozen:
        """Return the Frozen for one key, building it on a miss."""
        key = (case, lmin)

        if key in cache:
            cache.move_to_end(key)

            return cache[key]

        cache[key] = Frozen(case, lmin)

        while FROZEN_CACHE_MAXSIZE and len(cache) > FROZEN_CACHE_MAXSIZE:
            cache.popitem(last=False)

        return cache[key]

    return get


@pytest.mark.parametrize("case", CHECKED_CASES)
@FROZEN_GROUP
@pytest.mark.parametrize("lmin", CHECKED_ELLS)
def test_reference_is_a_reference(frozen, case: str, lmin: int) -> None:
    """Check that the reference converges before any method is compared to it.

    Raising the Gauss-Legendre order on a cell has to stop moving the answer,
    or the comparisons below measure the reference rather than the method. A
    Chebyshev cell carries up to 129 coefficients and settles when the order
    passes its degree.
    """
    reference = frozen(case, lmin).reference

    assert reference.n_cells > 0
    assert reference.worst_cell_move < REFERENCE_FLOOR


# Parameter order is cache order. pytest exhausts the top decorator fastest and
# the bottom one slowest, and the frozen fixture keeps eight keys: with `method`
# innermost the methods of one (case, lmin) run back to back on one
# build. With `case` innermost the cases stream through the cache for each method
# in turn and every key is rebuilt once per method.
@pytest.mark.parametrize("method", sorted(METHODS))
@pytest.mark.parametrize("case", CHECKED_CASES)
@FROZEN_GROUP
@pytest.mark.parametrize("lmin", CHECKED_ELLS)
def test_method_matches_reference(frozen, case: str, method: str, lmin: int) -> None:
    """Check every kernel-space method against the quadrature's own reference."""
    state = frozen(case, lmin)
    error = state.peak_error(state.compute(method), state.reference.cl)
    tolerance = TOLERANCE[method]

    assert error < tolerance, (
        f"{case} ({state.pair.regime}): {method} "
        f"deviates by {error:.3e} of the block's peak at ell "
        f"{state.lmin}-{state.lmax}, above the measured {tolerance:.1e}"
    )


@pytest.mark.parametrize("method", sorted(METHODS))
@pytest.mark.parametrize("case", CHECKED_CASES)
@FROZEN_GROUP
@pytest.mark.parametrize("lmin", CHECKED_ELLS)
def test_block_entry_point_matches_compute(
    frozen, case: str, method: str, lmin: int
) -> None:
    """Check that integrate_block() and compute() agree exactly on shared closures.

    Exact equality, not a tolerance: compute() builds the closures and calls
    the same NcXcorKQuad entry this does, so any difference is the two paths
    disagreeing about which quadrature or which range they are on. That makes
    this the sharpest check on the table -- it removes the closure rebuild
    that every other comparison has to allow for.
    """
    state = frozen(case, lmin)
    direct = state.integrate_block(method)
    through_compute = state.compute(method)

    np.testing.assert_array_equal(
        direct,
        through_compute,
        err_msg=(
            f"{case}: {method} disagrees with itself "
            f"between Nc.Xcor.integrate_block() and nc_xcor_compute() at ell "
            f"{state.lmin}-{state.lmax}"
        ),
    )


@pytest.mark.parametrize("case", [pair.case for pair in cases.PAIRS if pair.isauto])
@FROZEN_GROUP
def test_auto_spectra_cannot_cancel(frozen, case: str) -> None:
    """Check that an auto spectrum does not cancel: W^2 >= 0 leaves nothing to.

    This is why the matrix is not all auto spectra: a method's per-node error
    reaches the answer undivided here, and by up to seven orders more on a
    far-separated cross pair. A case set that lost its cross pairs would still
    pass every tolerance above and measure nothing.
    """
    assert frozen(case, 2).cancellation == pytest.approx(1.0, abs=1.0e-12)


@FROZEN_GROUP
def test_far_separated_bins_cancel_by_orders(frozen) -> None:
    """Check that the far-separated cross pair cancels by orders of magnitude.

    Tail against tail at ell = 2 is where a quadrature's per-node error is
    amplified the most; by ell = 200 each window is confined near its own
    turning point and the overlap is small rather than cancelling.
    """
    low = frozen("X3", 2).cancellation[0]

    assert low > 1.0e3


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

# Measured over the 43 certified entries, with headroom. The worst is the
# far-separated pair at ell = 50, 1.04e-2 of scale, on a certified value eleven
# orders below the pair's peak.
#
# A loose gate has a cost, paid once already: the X9 entries were certified
# with the generator's side-b leak (kdep applied to both windows), putting
# the reference 5.3% off at ell = 50, and the closure sitting 2.9e-2 from the
# table stayed under the old gate of 1e-1. A tolerance-independent deviation is
# a wrong reference, not a bad fit -- worth a ladder check before widening the
# number here.
CERTIFIED_TOLERANCE = 3.0e-2

# Tail x tail with hard edges is its own regime rather than a harder version of
# the roster, and X14 is the one pair that leaves the general budget behind.
# Both top-hat transforms decay as power laws, so their product stays
# significant far out in k where a Chebyshev closure resolves the integrand
# worst: 1.34e-1 of the pair's scale at ell = 2, against 1.37e-2 for X10, the
# next pair down and the same configuration with Gaussian edges. Carving it out
# keeps 3e-2 a real statement about the other thirty-four.
CERTIFIED_TOLERANCE_BY_CASE = {"X14": 2.0e-1}


def _certified_budget(case: str) -> float:
    """Return the certified tolerance for one case, with headroom."""
    return CERTIFIED_TOLERANCE_BY_CASE.get(case, CERTIFIED_TOLERANCE)


@pytest.fixture(name="certified", scope="module")
def fixture_certified() -> dict:
    """Load the certified C_ell table once per module."""
    path = pathlib.Path(Ncm.cfg_get_data_filename(TRUTH_TABLE, True))

    with gzip.open(path, "rt") as handle:
        return json.load(handle)


def _pair_scale(table: dict) -> dict:
    """Return each pair's largest certified |C_ell| over the multipoles present."""
    scale: dict[str, float] = {}

    for entry in table["cases"].values():
        value = abs(float(entry["value"]))
        scale[entry["case"]] = max(scale.get(entry["case"], 0.0), value)

    return scale


def _library_cl(case: str, ell: int) -> float:
    """Return the library's C_ell for one (case, ell) entry."""
    pair = cases.PAIRS_BY_CASE[case]
    settings = cases.Settings()
    cosmo, dist, ps = cases.make_cosmo_bits()

    kernel_a = cases.build_kernel(pair.kernel_a, cosmo, dist, ps, settings)
    kernel_b = (
        kernel_a
        if pair.isauto
        else cases.build_kernel(pair.kernel_b, cosmo, dist, ps, settings)
    )

    xcor = Nc.Xcor.new(dist, ps, Nc.XcorMethod.KERNEL_EXACT)
    xcor.set_ell_batch_size(1)
    xcor.set_reltol(settings.reltol)
    xcor.prepare(cosmo)

    vp = Ncm.Vector.new(1)
    xcor.compute(kernel_a, kernel_b, cosmo, ell, ell, vp)

    return vp.get(0)


def test_certified_table_is_certified(certified: dict) -> None:
    """Check that every certified radius is far below the tolerance used here.

    Guards a regeneration that quietly lowered the precision target -- every
    assertion below would still pass while checking nothing.
    """
    worst = max(
        entry["radius"] / abs(float(entry["value"]))
        for entry in certified["cases"].values()
    )

    assert certified["cases"]
    assert worst < 1.0e-8


def test_closure_matches_certified_c_ell(certified: dict) -> None:
    """Check the closure against Arb over every entry the table certifies."""
    scale = _pair_scale(certified)
    worst: dict[str, tuple[float, int]] = {}

    for entry in certified["cases"].values():
        got = _library_cl(entry["case"], entry["ell"])
        deviation = abs(got - float(entry["value"])) / scale[entry["case"]]
        case = entry["case"]

        if deviation > worst.get(case, (0.0, 0))[0]:
            worst[case] = (deviation, entry["ell"])

    over = [
        (
            f"{case} ell={ell} by {deviation:.3e} "
            f"(budget {_certified_budget(case):.1e})"
        )
        for case, (deviation, ell) in sorted(worst.items())
        if deviation >= _certified_budget(case)
    ]

    assert not over, (
        f"the closure deviates from the certified C_ell, as a fraction of the "
        f"pair's scale: " + "; ".join(over)
    )

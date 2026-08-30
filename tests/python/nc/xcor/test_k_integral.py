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

import gzip
import json
import pathlib

import numpy as np
import pytest

from numcosmo_py import Nc, Ncm
from xcor import cases_k_integral as cases

pytest_plugins = ["python.fixtures_xcor"]

pytestmark = pytest.mark.xcor

CLOSURES = {
    "spline": Nc.XcorKernelClosure.SPLINE,
    "chebyshev": Nc.XcorKernelClosure.CHEBYSHEV,
}

METHODS = {
    "exact": Nc.XcorMethod.KERNEL_EXACT,
    "cubature": Nc.XcorMethod.KERNEL_CUBATURE,
    "gsl": Nc.XcorMethod.KERNEL_GSL,
}

# KERNEL_GSL calls nc_xcor_kernel_get_eval once per multipole; the other two
# call ..._get_eval_vectorized_full once per block.
PER_MULTIPOLE = frozenset({"gsl"})

# Measured over the whole matrix at reltol = scaled_abstol = 1e-4, as the worst
# deviation from the matching reference relative to the block's peak, with an
# order of headroom. Measured worst / median, for the record:
#
#   exact    spline     5.3e-13 / 9.4e-16       cubature spline     1.7e-05 / 3.2e-07
#   exact    chebyshev  4.7e-10 / 1.2e-15       cubature chebyshev  4.1e-05 / 1.6e-07
#   gsl      spline     1.1e-12 / 5.7e-16       gsl      chebyshev  1.7e-04 / 1.2e-07
#
# Read the two exact columns as the claim they are: GL(5) on the merged knot
# set is exact, and so is qagp on those same knots -- both sit at the
# reference's own floor. Cubature is not converging to the closure, it is
# stopping at the relative tolerance it was asked for, which is why it is eight
# orders looser at identical cost.
#
# The exact/chebyshev worst is N3 at ell ~ 205, where the reference is itself
# only good to 5.4e-10; that entry bounds the pair rather than discriminating
# between them, and there is nothing below it to find without an independent
# reference.
TOLERANCE = {
    ("exact", "spline"): 5.0e-12,
    ("exact", "chebyshev"): 1.0e-8,
    ("cubature", "spline"): 2.0e-4,
    ("cubature", "chebyshev"): 5.0e-4,
    ("gsl", "spline"): 1.0e-11,
    ("gsl", "chebyshev"): 2.0e-3,
}

# The reference's own convergence, likewise measured: 6.0e-14 on spline
# closures and 5.4e-10 on Chebyshev ones, where the floor is the conditioning
# of evaluating a degree-128 polynomial rather than the quadrature.
REFERENCE_FLOOR = {"spline": 1.0e-12, "chebyshev": 5.0e-9}


class Frozen:
    """Two prepared kernels, their closures for one block, and the references."""

    def __init__(self, case: str, closure: str, lmin: int) -> None:
        self.pair = cases.PAIRS_BY_CASE[case]
        self.settings = cases.Settings(closure=CLOSURES[closure])
        self.lmin = lmin
        self.lmax = lmin + self.settings.ell_batch_size - 1
        self.cosmo, self.dist, self.ps = cases.make_cosmo_bits()
        self.RH = Nc.HICosmo.RH_Mpc(self.cosmo)

        self.kernel_a = cases.build_kernel(
            self.pair.kernel_a, self.cosmo, self.dist, self.ps, self.settings
        )
        self.kernel_b = (
            self.kernel_a
            if self.pair.isauto
            else cases.build_kernel(
                self.pair.kernel_b, self.cosmo, self.dist, self.ps, self.settings
            )
        )
        self.integrand_a = cases.build_integrand(
            self.kernel_a, self.cosmo, lmin, self.lmax, self.settings
        )
        self.integrand_b = (
            None
            if self.pair.isauto
            else cases.build_integrand(
                self.kernel_b, self.cosmo, lmin, self.lmax, self.settings
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

    def truth_for(self, method: str) -> np.ndarray:
        """The reference built on the closure this method actually integrates."""
        return self.per_multipole if method in PER_MULTIPOLE else self.reference.cl

    def peak_error(self, got: np.ndarray, truth: np.ndarray) -> float:
        """Worst deviation over the block, against the block's own peak."""
        return float(np.abs(got - truth).max() / np.abs(truth).max())


@pytest.fixture(name="frozen", scope="module")
def fixture_frozen():
    """Build each (case, closure, block) once and hand it to every test."""
    cache: dict[tuple[str, str, int], Frozen] = {}

    def get(case: str, closure: str, lmin: int) -> Frozen:
        key = (case, closure, lmin)

        if key not in cache:
            cache[key] = Frozen(case, closure, lmin)

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
# of the pair's own scale. The tight statement is the comparative test below.
CERTIFIED_TOLERANCE = {"spline": 1.5e1, "chebyshev": 1.0e-1}


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

    # Measured: closer in 34 of 43, median 9.0e-6 against 5.9e-5.
    assert closer >= 0.7 * len(spline)
    assert np.median(chebyshev) < 0.5 * np.median(spline)

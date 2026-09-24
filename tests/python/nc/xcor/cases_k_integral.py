#
# cases_k_integral.py
#
# Fri Aug 28 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# cases_k_integral.py
# Copyright (C) 2026 Sandro Dias Pinto Vitenti <vitenti@uel.br>
#
# numcosmo is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# numcosmo is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program.  If not, see <http://www.gnu.org/licenses/>.

"""Case matrix and reference for the outer k-integral.

The outer integral is

    C_ell^AB = 2 / (pi RH^3) * INT dkappa kappa^2 W1_ell(kappa) W2_ell(kappa)

in the internal variable kappa = k RH, over the intersection of the two
closures' fitted ranges. The three kernel-space methods differ only in how they
evaluate that integral from the same two closures.

This module is the single source of truth for what is measured: the kernels,
the pairs, the multipoles, the tolerance settings, and the reference. It is
imported by ``test_k_integral.py`` and by ``bench_k_integral.py``; neither
carries a case of its own.

Two properties of the case set drive everything downstream:

- **Auto spectra are structurally easy.** W^2 >= 0, so the integrand never
  cancels and truncation can only bias one way. A suite of auto spectra cannot
  exercise the failure mode that matters.
- **The damage is at low ell.** Far-separated bins cancel by three orders at
  ell = 2 and not at all by ell = 50, so a sweep starting at ell ~ 10 misses
  the whole problem.

The reference (``reference_cl``) is deliberately *not* independent of the
closures: it is handed the same frozen ``NcXcorKernelIntegrand`` objects the
methods are handed, so what it measures is the outer quadrature alone. The
closures' own error against certified Arb values is measured separately, in
``test_xcor_window_truth_table.py``.
"""

from __future__ import annotations

import dataclasses
import typing

import numpy as np

from numcosmo_py import Nc, Ncm
from numcosmo_py.certified.xcor import kernels
from numcosmo_py.certified.xcor import windows as W
from numcosmo_py.certified.xcor.kernels import make_cosmo_bits

# The multipoles. Log-ish and weighted low, because that is where cancellation
# lives; ell = 2 is not decoration.
ELLS: typing.Final[list[int]] = [2, 3, 4, 6, 10, 20, 50, 100, 200]

# What the committed suite runs, as leading multipoles of a block. Chosen from
# the bench sweep over the full ELLS above, by where the methods are actually
# worst -- which is NOT the high end. Measured worst deviation of the block's
# peak, over the case matrix:
#
#   exact/chebyshev  2.3e-08 at l=6     gsl/chebyshev   1.6e-03 at l=10
#   exact/spline     8.1e-12 at l=6     gsl_block/spl   3.6e-11 at l=4
#
# against 4.7e-10, 1.7e-04, 5.3e-13 and 8.6e-12 at l=200. The error peaks
# around l = 4-10 and falls away by l = 50, so a ladder of [2, 20, 200] tested
# the easy end hard and the hard end not at all. l = 2 stays for the
# cancellation reason above; 6 and 10 are the peak; 50 keeps a high-l point
# without paying for l = 200, whose Levin operator is the most expensive of the
# ladder (its closure is the *smallest* -- 159 knots against 263 at l = 2 --
# so the cost is the operator, not the fit).
#
# Knowingly not covered: cubature/chebyshev is worst at l = 20 (4.1e-05 against
# 1.1e-05 here), and gsl_block/chebyshev at l = 200 (4.6e-04 against 3.2e-04).
# Both are within a factor of four of a point that is covered.
ELLS_SUITE: typing.Final[list[int]] = [2, 6, 10, 50]

# Gauss-Legendre orders the reference escalates through, per cell. A spline
# cell is a cubic and a k^2-weighted product of two is degree 8, so it is exact
# at the first order; a Chebyshev cell carries up to 129 coefficients, so the
# escalation is what the top of this list is for.
_GL_ORDERS: typing.Final[tuple[int, ...]] = (8, 16, 32, 64, 128, 264)

# The order of the first pass, which exists only to learn the block's peak.
_SCALE_ORDER: typing.Final[int] = 64

_GL_CACHE: dict[int, tuple[np.ndarray, np.ndarray]] = {}


def _gl(order: int) -> tuple[np.ndarray, np.ndarray]:
    """Gauss-Legendre nodes and weights on [-1, 1], cached across cases."""
    if order not in _GL_CACHE:
        _GL_CACHE[order] = np.polynomial.legendre.leggauss(order)

    return _GL_CACHE[order]


@dataclasses.dataclass(frozen=True)
class Settings:
    """Everything a case is run at, carried with the case rather than the runner.

    A ``peak_epsilon`` copied out of a test file once turned a 0.04 s top-hat
    into 19.9 s, so these are recorded in every row the bench driver emits.

    ``l_limber = -1`` is not a detail: at the library default of 0 every
    multipole takes Limber, Limber keeps the spline closure whatever
    ``closure`` says, and a Chebyshev sweep silently measures splines.
    """

    reltol: float = 1.0e-4
    peak_epsilon: float = 1.0e-4
    l_limber: int = -1
    ell_batch_size: int = 8
    closure: Nc.XcorKernelClosure = Nc.XcorKernelClosure.SPLINE

    def replace(self, **kwargs: typing.Any) -> "Settings":
        """A copy with some fields changed."""
        return dataclasses.replace(self, **kwargs)


@dataclasses.dataclass(frozen=True)
class KernelSpec:
    """One window of the pair matrix: its vocabulary name and what it stresses.

    The window itself is the vocabulary's (:mod:`numcosmo_py.certified.xcor.windows`),
    which the library kernel and the Arb generator are both driven from; this roster
    adds only the reason each window is in the matrix.
    """

    name: str
    stresses: str

    def build(
        self, dist: Nc.Distance, ps: Ncm.Powspec, sbi: Ncm.SBesselIntegrator
    ) -> Nc.XcorKernelRadial:
        """Construct this window's kernel, not yet prepared."""
        return kernels.build_kernel(self.name, dist, ps, sbi)

    def arb_args(self, side: str) -> list[str]:
        """This window as arguments to ``nc_xcor_kquad_arb``."""
        return W.arb_args(self.name, side=side)


KERNELS: typing.Final[dict[str, KernelSpec]] = {
    spec.name: spec
    for spec in (
        KernelSpec("gauss_mid", "baseline"),
        KernelSpec("gauss_thin", "thin bin"),
        KernelSpec("gauss_thin_shift", "thin bin, displaced"),
        # Overlaps gauss_mid rather than sitting nearer than it; `gauss_near` is
        # the SRD Y10 lens bin at 1095 Mpc.
        KernelSpec("gauss_over", "overlapping bin"),
        KernelSpec("gauss_low", "near bin"),
        KernelSpec("gauss_high", "far bin"),
        KernelSpec("tophat", "hard edges, 2000 Mpc wide"),
        # The regime the Chebyshev closure exists for. On a shell this narrow
        # the adaptive spline closure's accuracy is capped by peak-epsilon at
        # the library's own documented floor, where it sits at ~2e-5 while the
        # spectral closure reaches machine zero. No case built from wide bins
        # reaches it.
        KernelSpec("shell_narrow", "56 Mpc hard shell"),
        KernelSpec("shell_wide", "400 Mpc hard shell"),
        KernelSpec("tophat_smooth", "a real tomographic bin"),
        KernelSpec("student_t", "power-law tail"),
        KernelSpec("multi_disjoint", "disconnected support"),
        KernelSpec("power_exp", "skewed and broad"),
        KernelSpec("lensing", "broad, smoothed source"),
        # The equal-variance top-hat partners of the SRD Gaussians.
        KernelSpec(
            "tophat_near", "hard edges at the SRD Y10 lens 0 centre and variance"
        ),
        KernelSpec(
            "tophat_far", "hard edges at the SRD Y10 lens 9 centre and variance"
        ),
        KernelSpec("tophat_broad", "hard edges at the source-bin centre and variance"),
        # The SRD's own locations and widths, so the far-separated regime below is
        # a cross a 3x2pt analysis actually computes rather than a constructed one.
        KernelSpec("srd_lens0", "LSST-SRD Y10 lens bin 0"),
        KernelSpec("srd_lens9", "LSST-SRD Y10 lens bin 9: sigma/chi = 0.043"),
        # A source bin is not a lens bin, and the difference is the shape rather
        # than the width: measured on the SRD's own dn/dz in chi, the Y10 lens bins
        # are Gaussian to skew 0.05, while Y1 source bin 4 reaches skew 0.82. This
        # is a skewed window at that location and width in round numbers, not a fit
        # to the bin -- the roster has to span the regimes in forms Arb certifies
        # exactly, not track a survey's dn/dz through its data releases.
        KernelSpec("srd_source4", "skewed window at source-bin location and width"),
        # Bessel-derivative weights. The RSD term of NcXcorKernelGal is a component
        # of order 2, but its window is a spline over a sampled dn/dz and cannot be
        # certified to 1e-10; the same weight on closed-form windows gives certified
        # C_ell for the j_ell'' path, on every shape rather than on one.
        KernelSpec("srd_lens9_rsd", "j_ell'' on the SRD's thinnest bin"),
        KernelSpec("gauss_rsd", "j_ell'' weight: the redshift-space term"),
        KernelSpec(
            "gauss_deriv1", "j_ell' weight: first order, the control for the second"
        ),
        KernelSpec("tophat_rsd", "j_ell'' on hard edges"),
        KernelSpec("tophat_smooth_rsd", "j_ell'' on a smoothed edge"),
        KernelSpec("student_t_rsd", "j_ell'' on a power-law tail"),
        KernelSpec("power_exp_rsd", "j_ell'' on a skewed, broad window"),
        KernelSpec("lensing_rsd", "j_ell'' on a broad, smoothed source"),
        KernelSpec("multi_disjoint_rsd", "j_ell'' on disconnected support"),
        KernelSpec("shell_wide_rsd", "j_ell'' on a wide hard shell"),
        KernelSpec("shell_narrow_rsd", "j_ell'' on a narrow hard shell"),
        KernelSpec("kdep", "non-separable W: sqrt(P) no longer factors out"),
    )
}


@dataclasses.dataclass(frozen=True)
class PairSpec:
    """One spectrum to compute, and the regime it exists to probe."""

    case: str
    kernel_a: str
    kernel_b: str
    regime: str

    @property
    def isauto(self) -> bool:
        """Whether both sides are the same kernel."""
        return self.kernel_a == self.kernel_b


PAIRS: typing.Final[list[PairSpec]] = [
    PairSpec("A1", "gauss_mid", "gauss_mid", "positive integrand, baseline"),
    PairSpec("A2", "tophat", "tophat", "positive, hard-edge cost"),
    PairSpec("A3", "student_t", "student_t", "positive, power-law tail"),
    PairSpec("A4", "lensing", "lensing", "positive, broad"),
    PairSpec("A5", "multi_disjoint", "multi_disjoint", "positive, disconnected"),
    PairSpec("X1", "gauss_mid", "gauss_over", "overlapping, benign"),
    PairSpec("X2", "gauss_mid", "gauss_low", "separated peaks"),
    PairSpec("X3", "gauss_low", "gauss_high", "far separated, tail x tail"),
    PairSpec("X4", "gauss_thin", "gauss_thin_shift", "thin x thin"),
    PairSpec("X5", "gauss_mid", "tophat", "smooth x discontinuous"),
    PairSpec("X6", "gauss_mid", "lensing", "narrow x broad"),
    PairSpec("X7", "tophat", "tophat_smooth", "edge x smoothed edge"),
    PairSpec("X8", "power_exp", "gauss_high", "skewed x narrow, tail-dominated"),
    PairSpec("X9", "kdep", "gauss_mid", "non-separable W"),
    # One derivative-weighted pair per shape, each mirroring that shape's own
    # d = 0 pair above so the two layouts can be read side by side. R1 and R3 are
    # the controls on the gauss window: the auto pair squares the weight, and R3
    # carries first order, which is what makes the second-order numbers readable.
    PairSpec("R1", "gauss_rsd", "gauss_rsd", "j_ell'' auto, the RSD weight squared"),
    PairSpec("R2", "gauss_rsd", "gauss_mid", "j_ell'' x j_ell, density x RSD"),
    PairSpec("R3", "gauss_deriv1", "gauss_mid", "j_ell' x j_ell, first order"),
    PairSpec("R4", "tophat_rsd", "tophat", "j_ell'' x j_ell, hard edges (mirrors A2)"),
    PairSpec(
        "R5",
        "tophat_smooth_rsd",
        "tophat_smooth",
        "j_ell'' x j_ell, smoothed edge (mirrors X7)",
    ),
    PairSpec(
        "R6", "student_t_rsd", "student_t", "j_ell'' x j_ell, power-law tail (A3)"
    ),
    PairSpec(
        "R7", "power_exp_rsd", "gauss_high", "j_ell'' x j_ell, skewed x narrow (X8)"
    ),
    PairSpec("R8", "lensing_rsd", "gauss_mid", "j_ell'' x j_ell, broad x narrow (X6)"),
    PairSpec(
        "R9",
        "multi_disjoint_rsd",
        "multi_disjoint",
        "j_ell'' x j_ell, disconnected support (A5)",
    ),
    PairSpec(
        "R10", "shell_narrow_rsd", "shell_narrow", "j_ell'' x j_ell, hard shell (N1)"
    ),
    PairSpec(
        "R11", "shell_wide_rsd", "shell_wide", "j_ell'' x j_ell, wide hard shell (N2)"
    ),
    PairSpec(
        "R12",
        "shell_narrow_rsd",
        "shell_wide",
        "j_ell'' x j_ell, narrow x wide hard shells (N3)",
    ),
    # Tail x tail, from the survey rather than constructed: bins 0 and 9 of the
    # Y10 lens sample have disjoint supports, so their spectrum lives entirely in
    # the product of two exponential tails -- the regime where the spline closure
    # returns the wrong sign at loose tolerance (section 14). X3 probes the same
    # corner with chosen parameters; these say it is not a corner case.
    # The separation ladder: one pair per rung, the ratio of the two turning
    # points being 1.4, 2.3 and 3.5. Only the separation differs, so a method's
    # behaviour across the three is a statement about the configuration.
    PairSpec("X12", "srd_lens0", "gauss_mid", "adjacent bins, ratio 1.4"),
    PairSpec("X13", "gauss_mid", "gauss_high", "bins apart, ratio 2.3"),
    PairSpec("X10", "srd_lens0", "srd_lens9", "tail x tail, SRD lens 0 x lens 9"),
    # The same tail x tail configuration with hard edges instead of Gaussian ones:
    # these two top-hats carry the *same centres and the same variances* as X10's
    # Gaussians, so the pair isolates the edge and nothing else. It is a different
    # regime rather than a harder version of the same one -- a top-hat's transform
    # decays as a power law, so two of them multiply to something that stays
    # significant over a far wider k range than two exponential tails do, and the
    # certified k range has to follow it out there.
    #
    # Note the existing tophat/tophat disjoint pairs, N3 and R12, are not this:
    # shell_narrow and shell_wide sit at 1028 and 1300 Mpc, a separation ratio of
    # 1.26, which is the adjacent regime with a 44 Mpc gap between the supports.
    PairSpec(
        "X14", "tophat_near", "tophat_far", "tail x tail with hard edges (mirrors X10)"
    ),
    PairSpec(
        "X11", "srd_lens9", "srd_source4", "thin x skewed broad, SRD lens 9 x source 4"
    ),
    PairSpec(
        "R13",
        "srd_lens9_rsd",
        "srd_lens0",
        "j_ell'' x j_ell, tail x tail on SRD bins",
    ),
    PairSpec("N1", "shell_narrow", "shell_narrow", "narrow hard shell, auto"),
    PairSpec("N2", "shell_wide", "shell_wide", "wide hard shell, auto"),
    PairSpec("N3", "shell_narrow", "shell_wide", "narrow x wide hard shells"),
]

PAIRS_BY_CASE: typing.Final[dict[str, PairSpec]] = {p.case: p for p in PAIRS}


def build_kernel(
    name: str,
    cosmo: Nc.HICosmo,
    dist: Nc.Distance,
    ps: Ncm.Powspec,
    settings: Settings,
    sbi: Ncm.SBesselIntegrator | None = None,
) -> Nc.XcorKernelRadial:
    """Build and prepare one kernel of the matrix at the given settings.

    Pass @sbi to share one Levin integrator across the kernels of a block, the
    way nc_xcor_solver_solve() does -- it dups one integrator per block and
    hands the same one to every kernel, which is what makes the stored
    Givens-rotation decomposition reusable across kernels (plan doc
    dev-notes/xcor_ultralevin_batching_plan.md sec. 6.1). With %None each kernel
    gets its own, which is safe under concurrency -- the sbessel ODE solver is
    not reentrant, so a shared one corrupts memory across concurrent blocks --
    but disables that reuse entirely.
    """
    kernel = KERNELS[name].build(
        dist, ps, sbi if sbi is not None else Ncm.SBesselIntegratorLevin.new(0, 8)
    )
    kernel.set_l_limber(settings.l_limber)
    kernel.set_property("reltol", settings.reltol)
    kernel.set_property("peak-epsilon", settings.peak_epsilon)
    kernel.prepare(cosmo)

    return kernel


def build_integrand(
    kernel: Nc.XcorKernelRadial,
    cosmo: Nc.HICosmo,
    lmin: int,
    lmax: int,
    settings: Settings,
    sbi: Ncm.SBesselIntegrator | None = None,
) -> Nc.XcorKernelIntegrand:
    """Freeze one closure for an ell block, and check it is the one asked for.

    A Chebyshev closure carries panels and a spline closure carries knots, so
    the two are told apart by what came back rather than by what was requested
    -- which is the only way to catch the silent fallback to splines when the
    block is taken under Limber.
    """
    integrand = kernel.get_eval_vectorized_full(
        cosmo,
        lmin,
        lmax,
        sbi if sbi is not None else Ncm.SBesselIntegratorLevin.new(lmin, lmax),
        settings.closure,
    )

    wanted_panels = settings.closure == Nc.XcorKernelClosure.CHEBYSHEV
    got_panels = integrand.get_n_panels() > 0

    if wanted_panels != got_panels:
        raise AssertionError(
            f"asked for a {settings.closure.value_nick} closure and got "
            f"{'panels' if got_panels else 'knots'}; l_limber is "
            f"{kernel.get_l_limber()} and multipoles taken under Limber keep "
            f"the spline closure"
        )

    return integrand


def breakpoints(integrand: Nc.XcorKernelIntegrand) -> np.ndarray:
    """Where the closure stops being one polynomial: panel edges, or knots."""
    n_panels = integrand.get_n_panels()

    if n_panels > 0:
        edges = []

        for i in range(n_panels):
            _, a, b = integrand.peek_panel(i)
            edges.extend((a, b))

        return np.unique(np.array(edges))

    knots = integrand.peek_knots()

    if knots is None:
        a, b = integrand.get_range()

        return np.array([a, b])

    return np.unique(np.array(knots.dup_array()))


def _neumaier_sum(terms: list[np.ndarray]) -> np.ndarray:
    """Compensated summation over cells.

    The cells of a far-separated cross pair cancel by three orders, which is
    the only precision hazard in this reference and the one place plain
    accumulation would show.
    """
    total = np.zeros_like(terms[0])
    compensation = np.zeros_like(terms[0])

    for term in terms:
        candidate = total + term
        big = np.abs(total) >= np.abs(term)
        compensation += np.where(
            big, (total - candidate) + term, (term - candidate) + total
        )
        total = candidate

    return total + compensation


def _cell_values(
    integrand1: Nc.XcorKernelIntegrand,
    integrand2: typing.Optional[Nc.XcorKernelIntegrand],
    a: float,
    b: float,
    order: int,
) -> np.ndarray:
    """kappa^2 W1 W2 over one cell by a fixed Gauss-Legendre rule."""
    mid = 0.5 * (a + b)
    half = 0.5 * (b - a)
    nodes, weights = _gl(order)
    kappa = mid + half * nodes
    values = np.array([integrand1.eval_array(k) for k in kappa])

    if integrand2 is None:
        product = values * values
    else:
        product = values * np.array([integrand2.eval_array(k) for k in kappa])

    return half * np.einsum("n,n,nl->l", weights, kappa * kappa, product, optimize=True)


def _cell_integral(
    integrand1: Nc.XcorKernelIntegrand,
    integrand2: typing.Optional[Nc.XcorKernelIntegrand],
    a: float,
    b: float,
    atol: float,
) -> tuple[np.ndarray, float, int]:
    """One cell, by Gauss-Legendre raised until the answer stops moving.

    The stopping test is *absolute*, against the whole block's peak rather than
    against this cell's own size. A cell whose own value nearly cancels -- an
    oscillation caught between two nodes of the panel structure -- can never
    settle a relative criterion, and would report a converged reference as
    unconverged while contributing nothing.
    """
    previous = None
    move = np.inf
    order = _GL_ORDERS[0]

    for order in _GL_ORDERS:
        current = _cell_values(integrand1, integrand2, a, b, order)

        if previous is not None:
            move = float(np.abs(current - previous).max())

            if move < atol:
                return current, move, order

        previous = current

    return typing.cast(np.ndarray, previous), move, order


@dataclasses.dataclass
class Reference:
    """A reference C_ell with the diagnostics needed to believe it."""

    cl: np.ndarray
    worst_cell_move: float  # as a fraction of the block's peak
    max_order: int
    n_cells: int
    k_min: float
    k_max: float


def reference_cl(
    RH: float,
    integrand1: Nc.XcorKernelIntegrand,
    integrand2: typing.Optional[Nc.XcorKernelIntegrand],
    rtol: float = 1.0e-14,
) -> Reference:
    """The outer integral over frozen closures, to the accuracy doubles allow.

    Composite Gauss-Legendre on the two closures' merged breakpoints, the order
    raised per cell until refinement stops moving it, accumulated with
    compensated summation. Inside a cell each closure is one polynomial, so
    raising the order converges rather than merely refines.

    ``integrand2`` is None for an auto spectrum. ``RH`` is the Hubble radius
    in Mpc, the internal variable's scale.

    This is the k-integral's own truth: the same closures the methods are
    handed, so a disagreement is the quadrature's and nothing else.
    """
    a1, b1 = integrand1.get_range()

    if integrand2 is None:
        k_min, k_max = a1, b1
    else:
        a2, b2 = integrand2.get_range()
        k_min, k_max = max(a1, a2), min(b1, b2)

    nell = integrand1.get_len()

    if k_min >= k_max:
        return Reference(np.zeros(nell), 0.0, 0, 0, k_min, k_max)

    edges = breakpoints(integrand1)

    if integrand2 is not None:
        edges = np.unique(np.concatenate([edges, breakpoints(integrand2)]))

    edges = np.unique(np.clip(edges[(edges > k_min) & (edges < k_max)], k_min, k_max))
    edges = np.concatenate([[k_min], edges, [k_max]])

    cells = [(a, b) for a, b in zip(edges[:-1], edges[1:]) if b > a]

    # A first pass at a fixed order, only to learn the block's peak. Every
    # tolerance downstream is against that peak, so the reference's own
    # stopping test has to be too.
    estimate = _neumaier_sum(
        [_cell_values(integrand1, integrand2, a, b, _SCALE_ORDER) for a, b in cells]
    )
    scale = float(np.abs(estimate).max())

    if scale == 0.0:
        return Reference(np.zeros(nell), 0.0, _SCALE_ORDER, len(cells), k_min, k_max)

    terms = []
    worst_move = 0.0
    max_order = 0

    for a, b in cells:
        value, move, order = _cell_integral(integrand1, integrand2, a, b, rtol * scale)
        terms.append(value)
        worst_move = max(worst_move, move)
        max_order = max(max_order, order)

    const_factor = 2.0 / (np.pi * RH**3)

    return Reference(
        const_factor * _neumaier_sum(terms),
        worst_move / scale,
        max_order,
        len(cells),
        k_min,
        k_max,
    )


def per_multipole_reference(
    RH: float,
    kernel_a: Nc.XcorKernelRadial,
    kernel_b: typing.Optional[Nc.XcorKernelRadial],
    cosmo: Nc.HICosmo,
    lmin: int,
    lmax: int,
    settings: Settings,
) -> np.ndarray:
    """The same integral, from one closure per multipole.

    This is what ``NC_XCOR_METHOD_KERNEL_GSL`` integrates: it calls
    ``nc_xcor_kernel_get_eval`` per multipole rather than
    ``..._get_eval_vectorized_full`` per block, so a block-closure reference
    would charge it for a closure it never built. On a strongly cancelling pair
    the two closures do not agree at the library's default tolerances -- see
    ``test_k_integral.py`` -- so which one the reference is built on is the
    difference between measuring a quadrature and measuring that disagreement.
    """
    values = []

    for ell in range(lmin, lmax + 1):
        integrand_a = build_integrand(kernel_a, cosmo, ell, ell, settings)
        integrand_b = (
            None
            if kernel_b is None
            else build_integrand(kernel_b, cosmo, ell, ell, settings)
        )
        values.append(reference_cl(RH, integrand_a, integrand_b).cl[0])

    return np.array(values)


def cancellation_ratio(
    integrand1: Nc.XcorKernelIntegrand,
    integrand2: typing.Optional[Nc.XcorKernelIntegrand],
    n_samples: int = 20000,
) -> np.ndarray:
    """INT |kappa^2 W1 W2| / |INT kappa^2 W1 W2|, per multipole.

    The number of digits any quadrature loses before it starts: a method whose
    per-node relative error is eps delivers at best C eps. It is 1 for every
    auto spectrum by construction and reaches ~1e3 for far-separated bins at
    ell = 2, so it is what makes the two comparable. Sampled on a fixed grid --
    indicative for the oscillatory cases, not a quadrature.
    """
    a1, b1 = integrand1.get_range()

    if integrand2 is None:
        k_min, k_max = a1, b1
    else:
        a2, b2 = integrand2.get_range()
        k_min, k_max = max(a1, a2), min(b1, b2)

    nell = integrand1.get_len()

    if k_min >= k_max:
        return np.ones(nell)

    kappa = np.linspace(k_min, k_max, n_samples)
    values = np.array([integrand1.eval_array(k) for k in kappa])

    if integrand2 is None:
        product = values * values
    else:
        product = values * np.array([integrand2.eval_array(k) for k in kappa])

    weighted = (kappa * kappa)[:, None] * product
    signed = np.abs(np.trapezoid(weighted, kappa, axis=0))
    absolute = np.trapezoid(np.abs(weighted), kappa, axis=0)

    return np.where(signed > 0.0, absolute / np.maximum(signed, 1.0e-300), np.inf)

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

# The multipoles. Log-ish and weighted low, because that is where cancellation
# lives; ell = 2 is not decoration.
ELLS: typing.Final[list[int]] = [2, 3, 4, 6, 10, 20, 50, 100, 200]

# What the committed suite runs, as leading multipoles of a block: the two ends
# and one in between. The full list is the bench driver's, where cost is not
# charged to every CI run.
ELLS_SUITE: typing.Final[list[int]] = [2, 20, 200]

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

    A ``scaled_abstol`` copied out of a test file once turned a 0.04 s top-hat
    into 19.9 s, so these are recorded in every row the bench driver emits.

    ``l_limber = -1`` is not a detail: at the library default of 0 every
    multipole takes Limber, Limber keeps the spline closure whatever
    ``closure`` says, and a Chebyshev sweep silently measures splines.
    """

    reltol: float = 1.0e-4
    scaled_abstol: float = 1.0e-4
    l_limber: int = -1
    ell_batch_size: int = 8
    closure: Nc.XcorKernelClosure = Nc.XcorKernelClosure.SPLINE

    def replace(self, **kwargs: typing.Any) -> "Settings":
        """A copy with some fields changed."""
        return dataclasses.replace(self, **kwargs)


@dataclasses.dataclass(frozen=True)
class KernelSpec:
    """One analytic window, as data: what to build and what it stresses.

    ``shape`` and ``params`` are the window itself, kept as plain data because
    two independent implementations are driven from them -- the library, via
    ``builder``, and the Arb reference generator, via its command line. A
    window defined twice is a window that can drift.
    """

    name: str
    builder: typing.Callable[..., Nc.XcorKernelAnalytic]
    stresses: str
    shape: str = ""
    params: typing.Mapping[str, typing.Any] = dataclasses.field(default_factory=dict)

    def arb_args(self, side: str) -> list[str]:
        """This window as arguments to ``nc_xcor_kquad_arb``."""
        out = [f"--{side}:shape={self.shape}"]

        for key, value in self.params.items():
            if isinstance(value, (list, tuple)):
                value = ",".join(repr(float(v)) for v in value)

            out.append(f"--{side}:{key}={value}")

        return out


def _gauss(chi_mean: float, chi_sigma: float, n_sigma: float = 4.0):
    def build(dist, ps, sbi):
        return Nc.XcorKernelAnalyticGauss.new_full(
            dist, ps, chi_mean, chi_sigma, n_sigma, sbi
        )

    return build


def _tophat(chi_lower: float, chi_upper: float):
    def build(dist, ps, sbi):
        return Nc.XcorKernelAnalyticTophat.new_full(dist, ps, chi_lower, chi_upper, sbi)

    return build


def _tophat_smooth(
    chi_lower: float, chi_upper: float, chi_sigma: float, n_sigma: float
):
    def build(dist, ps, sbi):
        return Nc.XcorKernelAnalyticTophatSmooth.new_full(
            dist, ps, chi_lower, chi_upper, chi_sigma, n_sigma, sbi
        )

    return build


def _student_t(chi_mean: float, chi_scale: float, nu: float, n_scale: float):
    def build(dist, ps, sbi):
        return Nc.XcorKernelAnalyticStudentT.new_full(
            dist, ps, chi_mean, chi_scale, nu, n_scale, sbi
        )

    return build


def _power_exp(
    chi_scale: float, alpha: float, beta: float, chi_lower: float, chi_upper: float
):
    def build(dist, ps, sbi):
        return Nc.XcorKernelAnalyticPowerExp.new_full(
            dist, ps, chi_scale, alpha, beta, chi_lower, chi_upper, sbi
        )

    return build


def _lensing(chi_lower: float, chi_source_lower: float, chi_source_upper: float):
    def build(dist, ps, sbi):
        return Nc.XcorKernelAnalyticLensing.new_full(
            dist, ps, chi_lower, chi_source_lower, chi_source_upper, sbi
        )

    return build


def _multi(
    means: list[float], sigmas: list[float], weights: list[float], n_sigma: float
):
    def build(dist, ps, sbi):
        return Nc.XcorKernelAnalyticMulti.new_full(
            dist,
            ps,
            Ncm.Vector.new_array(means),
            Ncm.Vector.new_array(sigmas),
            Ncm.Vector.new_array(weights),
            n_sigma,
            sbi,
        )

    return build


def _gauss_kdep(
    chi_mean: float,
    chi_sigma: float,
    n_sigma: float,
    amplitude: float,
    k_transition: float,
    chi_ref: float,
):
    def build(dist, ps, sbi):
        return Nc.XcorKernelAnalyticGauss(
            dist=dist,
            powspec=ps,
            chi_mean=chi_mean,
            chi_sigma=chi_sigma,
            n_sigma=n_sigma,
            integrator=sbi,
            scale_dependence=Nc.XcorKernelAnalyticKDepGrowth.new(
                amplitude, k_transition, chi_ref
            ),
        )

    return build


KERNELS: typing.Final[dict[str, KernelSpec]] = {
    spec.name: spec
    for spec in (
        KernelSpec(
            "gauss_mid",
            _gauss(1500.0, 300.0),
            "baseline",
            "gauss",
            {"chi-mean": 1500.0, "chi-sigma": 300.0, "n-sigma": 4.0},
        ),
        KernelSpec(
            "gauss_thin",
            _gauss(1500.0, 50.0),
            "thin bin",
            "gauss",
            {"chi-mean": 1500.0, "chi-sigma": 50.0, "n-sigma": 4.0},
        ),
        KernelSpec(
            "gauss_thin_shift",
            _gauss(1650.0, 50.0),
            "thin bin, displaced",
            "gauss",
            {"chi-mean": 1650.0, "chi-sigma": 50.0, "n-sigma": 4.0},
        ),
        KernelSpec(
            "gauss_near",
            _gauss(1800.0, 300.0),
            "overlapping bin",
            "gauss",
            {"chi-mean": 1800.0, "chi-sigma": 300.0, "n-sigma": 4.0},
        ),
        KernelSpec(
            "gauss_low",
            _gauss(600.0, 100.0),
            "near bin",
            "gauss",
            {"chi-mean": 600.0, "chi-sigma": 100.0, "n-sigma": 4.0},
        ),
        KernelSpec(
            "gauss_high",
            _gauss(3500.0, 150.0),
            "far bin",
            "gauss",
            {"chi-mean": 3500.0, "chi-sigma": 150.0, "n-sigma": 4.0},
        ),
        KernelSpec(
            "tophat",
            _tophat(500.0, 2500.0),
            "hard edges, 2000 Mpc wide",
            "tophat",
            {"chi-lower": 500.0, "chi-upper": 2500.0},
        ),
        # The regime the Chebyshev closure exists for. On a shell this narrow
        # the adaptive spline closure's accuracy is capped by scaled-abstol at
        # the library's own documented floor, where it sits at ~2e-5 while the
        # spectral closure reaches machine zero. No case built from wide bins
        # reaches it.
        KernelSpec(
            "shell_narrow",
            _tophat(1000.0, 1056.0),
            "56 Mpc hard shell",
            "tophat",
            {"chi-lower": 1000.0, "chi-upper": 1056.0},
        ),
        KernelSpec(
            "shell_wide",
            _tophat(1100.0, 1500.0),
            "400 Mpc hard shell",
            "tophat",
            {"chi-lower": 1100.0, "chi-upper": 1500.0},
        ),
        KernelSpec(
            "tophat_smooth",
            _tophat_smooth(1000.0, 2000.0, 150.0, 6.0),
            "a real tomographic bin",
            "tophat_smooth",
            {
                "chi-lower": 1000.0,
                "chi-upper": 2000.0,
                "chi-sigma": 150.0,
                "n-sigma": 6.0,
            },
        ),
        KernelSpec(
            "student_t",
            _student_t(1500.0, 200.0, 2.0, 6.0),
            "power-law tail",
            "student_t",
            {"chi-mean": 1500.0, "chi-scale": 200.0, "nu": 2.0, "n-scale": 6.0},
        ),
        KernelSpec(
            "multi_disjoint",
            _multi([600.0, 2600.0], [100.0, 150.0], [1.0, 1.0], 4.0),
            "disconnected support",
            "multi",
            {
                "mu": [600.0, 2600.0],
                "sigma": [100.0, 150.0],
                "weight": [1.0, 1.0],
                "n-sigma": 4.0,
            },
        ),
        KernelSpec(
            "power_exp",
            _power_exp(1200.0, 2.0, 1.5, 50.0, 4000.0),
            "skewed and broad",
            "power_exp",
            {
                "chi-scale": 1200.0,
                "alpha": 2.0,
                "beta": 1.5,
                "chi-lower": 50.0,
                "chi-upper": 4000.0,
            },
        ),
        KernelSpec(
            "lensing",
            _lensing(50.0, 2000.0, 3000.0),
            "broad, smoothed source",
            "lensing",
            {"chi-lower": 50.0, "chi-source-lower": 2000.0, "chi-source-upper": 3000.0},
        ),
        KernelSpec(
            "kdep",
            _gauss_kdep(1500.0, 300.0, 4.0, 0.3, 0.05, 3000.0),
            "non-separable W: sqrt(P) no longer factors out",
            "gauss",
            {
                "chi-mean": 1500.0,
                "chi-sigma": 300.0,
                "n-sigma": 4.0,
                "kdep-amplitude": 0.3,
                "kdep-k-transition": 0.05,
                "kdep-chi-ref": 3000.0,
            },
        ),
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
    PairSpec("X1", "gauss_mid", "gauss_near", "overlapping, benign"),
    PairSpec("X2", "gauss_mid", "gauss_low", "separated peaks"),
    PairSpec("X3", "gauss_low", "gauss_high", "far separated, tail x tail"),
    PairSpec("X4", "gauss_thin", "gauss_thin_shift", "thin x thin"),
    PairSpec("X5", "gauss_mid", "tophat", "smooth x discontinuous"),
    PairSpec("X6", "gauss_mid", "lensing", "narrow x broad"),
    PairSpec("X7", "tophat", "tophat_smooth", "edge x smoothed edge"),
    PairSpec("X8", "power_exp", "gauss_high", "skewed x narrow, tail-dominated"),
    PairSpec("X9", "kdep", "gauss_mid", "non-separable W"),
    PairSpec("N1", "shell_narrow", "shell_narrow", "narrow hard shell, auto"),
    PairSpec("N2", "shell_wide", "shell_wide", "wide hard shell, auto"),
    PairSpec("N3", "shell_narrow", "shell_wide", "narrow x wide hard shells"),
]

PAIRS_BY_CASE: typing.Final[dict[str, PairSpec]] = {p.case: p for p in PAIRS}


def make_cosmo_bits() -> tuple[Nc.HICosmo, Nc.Distance, Ncm.Powspec]:
    """A cosmology, a distance, and a closed-form power spectrum.

    ``NcmPowspecAnalytic`` rather than a transfer-function spline: sqrt(P) sits
    *inside* the function the closure fits, so a splined P would put its own
    interpolation error under everything measured here, and it has no closed
    form for a reference to reproduce. BBKS is the look-alike because it is
    already nothing but a fitting formula -- elementary, hence exactly
    reproducible in extended precision, and it keeps the ln^2 k tail that a
    rational form loses.
    """
    cosmo = Nc.HICosmoDEXcdm.new()
    dist = Nc.Distance.new(5.0)
    dist.prepare(cosmo)
    ps = Ncm.PowspecAnalytic.new(
        Ncm.PowspecAnalyticShape.BBKS, Ncm.PowspecAnalyticGrowth.LCDM
    )

    return cosmo, dist, ps


def build_kernel(
    name: str,
    cosmo: Nc.HICosmo,
    dist: Nc.Distance,
    ps: Ncm.Powspec,
    settings: Settings,
) -> Nc.XcorKernelAnalytic:
    """Build and prepare one kernel of the matrix at the given settings.

    Each kernel gets its own Levin integrator: the sbessel ODE solver is not
    reentrant and a shared one corrupts memory across concurrent blocks.
    """
    kernel = KERNELS[name].builder(dist, ps, Ncm.SBesselIntegratorLevin.new(0, 8))
    kernel.set_l_limber(settings.l_limber)
    kernel.set_property("reltol", settings.reltol)
    kernel.set_property("scaled-abstol", settings.scaled_abstol)
    kernel.prepare(cosmo)

    return kernel


def build_integrand(
    kernel: Nc.XcorKernelAnalytic,
    cosmo: Nc.HICosmo,
    lmin: int,
    lmax: int,
    settings: Settings,
) -> Nc.XcorKernelIntegrand:
    """Freeze one closure for an ell block, and check it is the one asked for.

    A Chebyshev closure carries panels and a spline closure carries knots, so
    the two are told apart by what came back rather than by what was requested
    -- which is the only way to catch the silent fallback to splines when the
    block is taken under Limber.
    """
    integrand = kernel.get_eval_vectorized_full(
        cosmo, lmin, lmax, Ncm.SBesselIntegratorLevin.new(lmin, lmax), settings.closure
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
    kernel_a: Nc.XcorKernelAnalytic,
    kernel_b: typing.Optional[Nc.XcorKernelAnalytic],
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

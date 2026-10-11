"""Check a Chebyshev closure against W_ell(k) computed directly.

The direct value shares nothing with the closure but the components themselves: for
each component c,

    W_ell(k) = sum_c prefactor_c(k, ell) int kernel_c(chi, k) j_ell^(d_c)(k chi) dchi,

with its own Levin integrator over the component's chi limits, as the kernel's
non-Limber evaluator defines it, but without the evaluator, its truncation or any
fit. Points are taken inside every panel, between its Chebyshev-Lobatto nodes and
away from its edges, so each lies on a definite side of every cut, and on a
log-spaced grid over the closure's range, which a missed feature cannot avoid.

The bound at a point is what the closure was asked for: the acceptance test of its
panel, max(reltol vscale, abstol) with vscale the largest |W| on the panel's nodes
and abstol peak-epsilon times the smallest multipole peak, plus adaptive-epsilon
times the block's norm, below which a component's tail is truncated to zero.
"""

import math
from dataclasses import dataclass

import numpy as np

from numcosmo_py import Nc, Ncm


def direct_window(
    kernel: Nc.XcorKernel, cosmo: Nc.HICosmo, ell: int, k: float
) -> float:
    """W_ell(k) of a non-Limber kernel from its components, without its evaluator."""
    integrator = Ncm.SBesselIntegratorLevin.new(ell, ell)
    out = Ncm.Vector.new(1)
    total = 0.0

    for comp in kernel.get_component_list():
        chi_min, chi_max, _, _ = comp.get_limits(cosmo)
        integrator.integrate_deriv(
            lambda chi, kk, c=comp: c.eval_kernel(cosmo, chi, kk),
            chi_min,
            chi_max,
            k,
            comp.get_bessel_deriv(),
            out,
        )
        total += comp.eval_prefactor(cosmo, k, ell) * out.get(0)

    return total


@dataclass
class ClosureCheck:
    """The worst deviation of a closure from the direct values, over its bound."""

    worst: float
    at_ell: int
    at_k: float
    n_points: int


def _panel_points(a: float, b: float, n: int) -> list[float]:
    """Three points of the panel [a, b], halfway in angle between Chebyshev-Lobatto
    nodes: next to each edge and in the middle."""
    m = n - 1
    mid, half = 0.5 * (a + b), 0.5 * (b - a)

    return [mid + half * math.cos((j + 0.5) * math.pi / m) for j in (0, m // 2, m - 1)]


def check_against_direct(
    kernel: Nc.XcorKernel,
    cosmo: Nc.HICosmo,
    closure: Nc.XcorKernelIntegrand,
    lmin: int,
    n_grid: int = 200,
) -> ClosureCheck:
    """Compare a non-Limber Chebyshev closure of the multipoles lmin, lmin + 1, ...
    with the direct values, and return the worst deviation over its bound."""
    panels = closure.peek_panels()
    n_rows = panels.get_fdim()
    reltol = kernel.get_reltol()
    lo, hi = closure.get_range()

    points: list[tuple[float, int]] = []
    vscale: list[np.ndarray] = []

    for i in range(panels.get_n_panels()):
        coeffs, a, b = panels.get_panel(i)
        n = coeffs.ncols()
        nodes = [
            0.5 * (a + b) + 0.5 * (b - a) * math.cos(j * math.pi / (n - 1))
            for j in range(n)
        ]
        vscale.append(
            np.abs(np.array([closure.eval_array(x) for x in nodes])).max(axis=0)
        )
        points.extend((x, i) for x in _panel_points(a, b, n))

    for x in np.geomspace(lo, hi, n_grid + 2)[1:-1]:
        points.append((float(x), panels.find_panel(float(x))))

    got = np.array([closure.eval_array(x) for x, _ in points])
    direct = np.array(
        [
            [direct_window(kernel, cosmo, lmin + r, x) for r in range(n_rows)]
            for x, _ in points
        ]
    )

    peaks = np.abs(direct).max(axis=0)
    abstol = kernel.get_peak_epsilon() * peaks.min()
    truncation = kernel.get_adaptive_epsilon() * math.sqrt(float((peaks**2).sum()))

    worst = ClosureCheck(0.0, lmin, lo, len(points))

    for (x, i), g, d in zip(points, got, direct):
        bound = np.maximum(reltol * vscale[i], abstol) + truncation
        ratio = np.abs(g - d) / bound
        r = int(np.argmax(ratio))

        if ratio[r] > worst.worst:
            worst = ClosureCheck(float(ratio[r]), lmin + r, x, len(points))

    return worst

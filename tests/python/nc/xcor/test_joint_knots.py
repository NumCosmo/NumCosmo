#!/usr/bin/env python
#
# test_joint_knots.py
#
# Mon Aug 11 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_joint_knots.py
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

"""Shared-knot (joint) sampling of several kernels, and exact GL(5) assembly.

Each component of a kernel integrand is a cubic spline in k, so on every knot
panel the outer integrand k^2 W_i W_j is a degree-8 polynomial and a 5-node
Gauss-Legendre rule integrates it exactly. That only helps across kernel pairs
if the kernels share one knot set: independently adapted splines resolve the
same oscillation at nearby but distinct abscissas, so unioning two of them
costs as much as their sum.
"""

import numpy as np
import pytest
from numpy.testing import assert_allclose

from numcosmo_py import Nc, Ncm
from numcosmo_py.cosmology import Cosmology

pytest_plugins = [
    "python.fixtures_xcor",
]

pytestmark = pytest.mark.xcor

Z_BINS = [(0.1, 0.2), (0.3, 0.4), (0.6, 0.7)]

GL5_X, GL5_W = np.polynomial.legendre.leggauss(5)


def _kernels(cosmology: Cosmology) -> list[Nc.XcorKernel]:
    """Build non-Limber top-hat kernels sharing a cosmology."""
    out = []

    for z_lower, z_upper in Z_BINS:
        kernel = Nc.XcorKernelClusterTophat(
            dist=cosmology.dist,
            powspec=cosmology.ps_ml,
            z_lower=z_lower,
            z_upper=z_upper,
            integrator=Ncm.SBesselIntegratorLevin.new(0, 8),
        )
        kernel.set_l_limber(-1)
        kernel.set_scaled_abstol(1.0e-8)
        kernel.prepare(cosmology.cosmo)
        out.append(kernel)

    return out


def _knots(integrand: Nc.XcorKernelIntegrand) -> np.ndarray:
    """Extract an integrand's knots as a numpy array."""
    knots = integrand.peek_knots()

    return np.array([knots.get(i) for i in range(knots.len())])


def test_separate_integrand_exposes_knots(cosmology: Cosmology) -> None:
    """A spline-backed integrand reports the knots it is represented on."""
    kernel = _kernels(cosmology)[0]
    integrand = kernel.get_eval(cosmology.cosmo, 0)

    knots = _knots(integrand)
    k_min, k_max = integrand.get_range()

    assert knots.size > 2
    assert np.all(np.diff(knots) > 0.0)
    assert_allclose([knots[0], knots[-1]], [k_min, k_max], rtol=1.0e-14)


def test_joint_shares_one_knot_set(cosmology: Cosmology) -> None:
    """Joint sampling puts every kernel on one knot set, sized like the largest."""
    kernels = _kernels(cosmology)
    sbi = Ncm.SBesselIntegratorLevin.new(0, 8)

    separate = [_knots(k.get_eval(cosmology.cosmo, 0)) for k in kernels]
    joint = Nc.XcorKernel.get_eval_vectorized_joint(kernels, cosmology.cosmo, 0, 0, sbi)
    joint_knots = _knots(joint)

    assert joint.get_len() == len(kernels)  # one component per kernel at a single ell
    assert np.all(np.diff(joint_knots) > 0.0)

    # A jointly built set is a max-density set, not a sum: comfortably larger
    # than any single kernel's, but far below the sum of them all.
    per_kernel_max = max(k.size for k in separate)
    per_kernel_sum = sum(k.size for k in separate)

    assert joint_knots.size >= per_kernel_max
    assert joint_knots.size < 0.75 * per_kernel_sum

    # The shared domain is the union of the kernels' own domains.
    assert joint_knots[0] <= min(k[0] for k in separate) * (1.0 + 1.0e-12)
    assert joint_knots[-1] >= max(k[-1] for k in separate) * (1.0 - 1.0e-12)


def test_component_support_is_narrower_than_shared_domain(
    cosmology: Cosmology,
) -> None:
    """Joint components report their own support, not the shared domain."""
    kernels = _kernels(cosmology)
    sbi = Ncm.SBesselIntegratorLevin.new(0, 8)

    joint = Nc.XcorKernel.get_eval_vectorized_joint(kernels, cosmology.cosmo, 0, 0, sbi)
    shared_min, shared_max = joint.get_range()

    supports = []
    for c in range(joint.get_len()):
        ok, k_min, k_max = joint.get_component_range(c)
        assert ok
        assert shared_min <= k_min < k_max <= shared_max
        supports.append((k_min, k_max))

    # The shared domain is the union, so at least one component must stop short
    # of it -- otherwise there would be nothing to mask.
    assert min(s[1] for s in supports) < shared_max


def test_joint_gl5_assembly_matches_compute(cosmology: Cosmology) -> None:
    """Exact GL(5) on the shared panels reproduces nc_xcor_compute() for all pairs.

    This is the whole point of the shared knot set: one pass of 5 nodes per
    panel yields every pair at once as C_l = U^T W U, with no adaptive outer
    quadrature and hence no way for it to fail to converge.

    Each pair is restricted to the intersection of the two components' supports.
    Outside its own support a component's spline holds extrapolation, not a
    small number, and the k^2 weight turns that into the dominant contribution.
    """
    kernels = _kernels(cosmology)
    cosmo = cosmology.cosmo
    nbins = len(kernels)
    sbi = Ncm.SBesselIntegratorLevin.new(0, 8)

    joint = Nc.XcorKernel.get_eval_vectorized_joint(kernels, cosmo, 0, 0, sbi)
    edges = _knots(joint)
    supports = [joint.get_component_range(c)[1:] for c in range(nbins)]

    xcor = Nc.Xcor.new(cosmology.dist, cosmology.ps_ml, Nc.XcorMethod.KERNEL_CUBATURE)
    xcor.prepare(cosmo)
    const = 2.0 / (np.pi * cosmo.RH_Mpc() ** 3)

    cl = np.zeros((nbins, nbins))
    for i in range(nbins):
        for j in range(i, nbins):
            k_min = max(supports[i][0], supports[j][0])
            k_max = min(supports[i][1], supports[j][1])
            panel = edges[(edges >= k_min) & (edges <= k_max)]

            lo, hi = panel[:-1], panel[1:]
            nodes = np.concatenate(
                [0.5 * (b + a) + 0.5 * (b - a) * GL5_X for a, b in zip(lo, hi)]
            )
            weights = np.concatenate([0.5 * (b - a) * GL5_W for a, b in zip(lo, hi)])

            values = np.array([joint.eval_array(k) for k in nodes])
            cl[i, j] = cl[j, i] = np.sum(
                weights * nodes**2 * values[:, i] * values[:, j]
            )

            result = Ncm.Vector.new(1)
            xcor.compute(kernels[i], kernels[j], cosmo, 0, 0, result)

            assert_allclose(const * cl[i, j], result.get(0), rtol=5.0e-3)

    # Cauchy-Schwarz holds across the block.
    for i in range(nbins):
        for j in range(nbins):
            assert abs(cl[i, j]) <= np.sqrt(cl[i, i] * cl[j, j]) * (1.0 + 1.0e-12)


def test_kernel_fixed_matches_adaptive_methods(cosmology: Cosmology) -> None:
    """KERNEL_FIXED reproduces the adaptive kernel-space methods, with no tolerance."""
    kernels = _kernels(cosmology)
    cosmo = cosmology.cosmo
    nbins = len(kernels)

    fixed = Nc.Xcor.new(cosmology.dist, cosmology.ps_ml, Nc.XcorMethod.KERNEL_FIXED)
    fixed.prepare(cosmo)
    cubature = Nc.Xcor.new(
        cosmology.dist, cosmology.ps_ml, Nc.XcorMethod.KERNEL_CUBATURE
    )
    cubature.prepare(cosmo)

    got, expected = Ncm.Vector.new(1), Ncm.Vector.new(1)
    for i in range(nbins):
        for j in range(i, nbins):
            fixed.compute(kernels[i], kernels[j], cosmo, 0, 0, got)
            cubature.compute(kernels[i], kernels[j], cosmo, 0, 0, expected)

            assert_allclose(got.get(0), expected.get(0), rtol=5.0e-3)


def test_kernel_fixed_over_ell_block(cosmology: Cosmology) -> None:
    """KERNEL_FIXED handles a multipole block in one sweep of the shared panels.

    Compared only where both methods are above the block's absolute floor. The
    two disagree on entries far below it, and for a reason worth knowing:
    KERNEL_CUBATURE intersects the integrands' full ranges while KERNEL_FIXED
    intersects their per-component supports, which are narrower. For separated
    bins at high multipole the cross spectrum lives on that overlap edge -- a
    tail-times-tail case -- and the per-component criterion, which cannot know
    the pair needs the tail, clips it. See dev-notes/xcor_joint_knots.md sec 6.
    """
    kernels = _kernels(cosmology)
    cosmo = cosmology.cosmo
    lmin, lmax = 2, 5
    n_l = lmax - lmin + 1

    fixed = Nc.Xcor.new(cosmology.dist, cosmology.ps_ml, Nc.XcorMethod.KERNEL_FIXED)
    fixed.prepare(cosmo)
    cubature = Nc.Xcor.new(
        cosmology.dist, cosmology.ps_ml, Nc.XcorMethod.KERNEL_CUBATURE
    )
    cubature.prepare(cosmo)

    got, expected = Ncm.Vector.new(n_l), Ncm.Vector.new(n_l)
    fixed.compute(kernels[0], kernels[2], cosmo, lmin, lmax, got)
    cubature.compute(kernels[0], kernels[2], cosmo, lmin, lmax, expected)

    got_a = np.array(got.dup_array())
    expected_a = np.array(expected.dup_array())

    assert_allclose(
        got_a, expected_a, rtol=2.0e-2, atol=1.0e-3 * np.max(np.abs(expected_a))
    )


def test_kernel_fixed_through_solver(cosmology: Cosmology) -> None:
    """The solver drives KERNEL_FIXED consistently with nc_xcor_compute().

    Agreement is to sampling accuracy rather than bitwise: the solver samples
    every registered kernel onto one knot set per block, while nc_xcor_compute()
    samples only the pair, so the two build different (each individually valid)
    knot sets. The separated-bin cross is the most sensitive to that, being a
    small residual of a large cancellation.
    """
    kernels = _kernels(cosmology)
    cosmo = cosmology.cosmo

    fixed = Nc.Xcor.new(cosmology.dist, cosmology.ps_ml, Nc.XcorMethod.KERNEL_FIXED)
    fixed.prepare(cosmo)

    solver = Nc.XcorSolver.new()
    ids = [solver.register_kernel(k) for k in kernels]
    solver.request_cl(ids[0], ids[0], 0, 3)
    solver.request_cl(ids[0], ids[2], 0, 3)
    solver.plan_blocks(8)
    solver.solve(fixed, cosmo)

    for request, (i, j) in enumerate([(0, 0), (0, 2)]):
        expected = Ncm.Vector.new(4)
        fixed.compute(kernels[i], kernels[j], cosmo, 0, 3, expected)

        assert_allclose(
            np.array(solver.get_result(request).dup_array()),
            np.array(expected.dup_array()),
            rtol=2.0e-3,
            atol=1.0e-50,
        )


def test_joint_covers_ell_block(cosmology: Cosmology) -> None:
    """A joint integrand over an ell block is ordered kernel-major."""
    kernels = _kernels(cosmology)
    sbi = Ncm.SBesselIntegratorLevin.new(0, 8)
    lmin, lmax = 2, 5
    n_l = lmax - lmin + 1

    joint = Nc.XcorKernel.get_eval_vectorized_joint(
        kernels, cosmology.cosmo, lmin, lmax, sbi
    )

    assert joint.get_len() == len(kernels) * n_l

    k_mid = np.sqrt(joint.get_range()[0] * joint.get_range()[1])
    values = np.array(joint.eval_array(k_mid)).reshape(len(kernels), n_l)

    assert np.all(np.isfinite(values))

    # Component ik * n_l + il must be kernel ik evaluated at lmin + il.
    for ik, kernel in enumerate(kernels):
        single = kernel.get_eval_vectorized_full(
            cosmology.cosmo, lmin, lmax, Ncm.SBesselIntegratorLevin.new(0, 8)
        )
        assert_allclose(
            values[ik], np.array(single.eval_array(k_mid)), rtol=1.0e-3, atol=1.0e-30
        )

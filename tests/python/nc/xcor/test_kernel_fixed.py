#!/usr/bin/env python
#
# test_kernel_fixed.py
#
# Mon Aug 11 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_kernel_fixed.py
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

"""NC_XCOR_METHOD_KERNEL_FIXED: exact GL(5) on per-pair knot unions.

Each kernel is sampled independently -- the same closures KERNEL_CUBATURE
builds and NcXcorSolver caches -- so a pair's two splines live on different
abscissas. On the common refinement of those abscissas each spline is still a
single cubic piece per panel, so k^2 W_i W_j is degree 8 there and a 5-node
Gauss-Legendre rule integrates it exactly. Merging two knot sets is all the
coupling exactness needs; sampling the kernels onto one shared abscissa is not
required and costs about twice as much to produce.
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


def _kernels(cosmology: Cosmology, l_limber: int = -1) -> list[Nc.XcorKernel]:
    """Build top-hat kernels sharing a cosmology, in the requested tier."""
    out = []

    for z_lower, z_upper in Z_BINS:
        kernel = Nc.XcorKernelClusterTophat(
            dist=cosmology.dist,
            powspec=cosmology.ps_ml,
            z_lower=z_lower,
            z_upper=z_upper,
            integrator=Ncm.SBesselIntegratorLevin.new(0, 8),
        )
        kernel.set_l_limber(l_limber)
        # Library default scaled-abstol, deliberately: this used to override it
        # to 1e-8 and no assertion here needed it. scaled-abstol floors W(k)
        # against its own peak, but the C_ell integrand is k^2 W_a W_b, so the
        # floor enters squared -- the 1e-4 default is already 1e-12 there. Every
        # test in this file passes at 1e-8 (14.3 s), 1e-6 (5.5 s) and the
        # default (1.5 s), so the override bought nothing and cost 10x.
        kernel.prepare(cosmology.cosmo)
        out.append(kernel)

    return out


def _knots(integrand: Nc.XcorKernelIntegrand) -> np.ndarray:
    """Extract an integrand's knots as a numpy array."""
    knots = integrand.peek_knots()

    return np.array([knots.get(i) for i in range(knots.len())])


def test_integrand_exposes_knots(cosmology: Cosmology) -> None:
    """A spline-backed integrand reports the knots it is represented on."""
    kernel = _kernels(cosmology)[0]
    integrand = kernel.get_eval(cosmology.cosmo, 0)

    knots = _knots(integrand)
    k_min, k_max = integrand.get_range()

    assert knots.size > 2
    assert np.all(np.diff(knots) > 0.0)
    assert_allclose([knots[0], knots[-1]], [k_min, k_max], rtol=1.0e-14)


def test_kernels_keep_their_own_knot_sets(cosmology: Cosmology) -> None:
    """Kernels are sampled separately, so their knot sets differ.

    This is the premise of the union quadrature: two independently adapted
    splines resolve the same oscillation at nearby but distinct abscissas, so
    the pair has to be integrated on the merged set rather than on either one.
    """
    kernels = _kernels(cosmology)
    sets = [_knots(k.get_eval(cosmology.cosmo, 0)) for k in kernels]

    for i in range(len(sets)):
        for j in range(i + 1, len(sets)):
            union = np.union1d(sets[i], sets[j])

            # Genuinely different abscissas: the union is strictly larger than
            # either, and close to their sum.
            assert union.size > max(sets[i].size, sets[j].size)
            assert union.size > 0.75 * (sets[i].size + sets[j].size)


def test_union_gl5_matches_compute(cosmology: Cosmology) -> None:
    """Exact GL(5) on the merged knot set reproduces nc_xcor_compute().

    Reimplements the C quadrature in numpy, as an independent check that the
    merged panels really are panels on which both splines are single cubics.
    """
    kernels = _kernels(cosmology)
    cosmo = cosmology.cosmo
    nbins = len(kernels)

    igs = [k.get_eval(cosmo, 0) for k in kernels]
    edges = [_knots(ig) for ig in igs]

    fixed = Nc.Xcor.new(cosmology.dist, cosmology.ps_ml, Nc.XcorMethod.KERNEL_FIXED)
    fixed.prepare(cosmo)
    const = 2.0 / (np.pi * cosmo.RH_Mpc() ** 3)

    for i in range(nbins):
        for j in range(i, nbins):
            lo = max(igs[i].get_range()[0], igs[j].get_range()[0])
            hi = min(igs[i].get_range()[1], igs[j].get_range()[1])

            panel = np.union1d(edges[i], edges[j])
            panel = panel[(panel >= lo) & (panel <= hi)]

            a, b = panel[:-1], panel[1:]
            nodes = (0.5 * (a + b)[:, None] + 0.5 * (b - a)[:, None] * GL5_X).ravel()
            weights = (0.5 * (b - a)[:, None] * GL5_W).ravel()

            wi = np.array([igs[i].eval_array(k)[0] for k in nodes])
            wj = np.array([igs[j].eval_array(k)[0] for k in nodes])
            got = const * np.sum(weights * nodes**2 * wi * wj)

            result = Ncm.Vector.new(1)
            fixed.compute(kernels[i], kernels[j], cosmo, 0, 0, result)

            assert_allclose(got, result.get(0), rtol=1.0e-6)


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

            assert_allclose(got.get(0), expected.get(0), rtol=1.0e-4)


def test_kernel_fixed_over_ell_block(cosmology: Cosmology) -> None:
    """KERNEL_FIXED handles a multipole block, on the same range as cubature.

    Both methods now intersect the two integrands' fitted domains, so the
    separated-bin cross spectrum agrees far into the tail -- unlike the earlier
    shared-abscissa implementation, which masked each component to a narrower
    per-component support and clipped exactly this case.
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

    assert_allclose(
        np.array(got.dup_array()), np.array(expected.dup_array()), rtol=1.0e-3
    )


def test_kernel_fixed_batches_wide_ell_range(cosmology: Cosmology) -> None:
    """KERNEL_FIXED batches a range wider than NC_XCOR_KERNEL_MAX_ELL_BLOCK.

    A single k-space closure is capped at 64 multipoles, so an unbatched sweep
    aborted the process on any wider request -- while KERNEL_CUBATURE, which the
    two are meant to be interchangeable through, sliced the range into
    NcXcor:ell-batch-size sub-blocks and handled it. The batching is what makes
    the request legal, so the result must also be independent of where the
    batch boundaries fall.
    """
    kernel = _kernels(cosmology)[0]
    cosmo = cosmology.cosmo
    lmax = 70  # > NC_XCOR_KERNEL_MAX_ELL_BLOCK

    fixed = Nc.Xcor.new(cosmology.dist, cosmology.ps_ml, Nc.XcorMethod.KERNEL_FIXED)
    fixed.prepare(cosmo)

    whole = Ncm.Vector.new(lmax + 1)
    fixed.compute(kernel, kernel, cosmo, 0, lmax, whole)
    values = np.array(whole.dup_array())

    assert np.all(np.isfinite(values))
    assert np.all(values > 0.0)

    # Batches share no state, so a sub-range computed on its own reproduces the
    # corresponding slice of the wide run exactly.
    for lo, hi in [(0, 7), (8, 15), (64, 70)]:
        part = Ncm.Vector.new(hi - lo + 1)
        fixed.compute(kernel, kernel, cosmo, lo, hi, part)

        assert_allclose(
            np.array(part.dup_array()), values[lo : hi + 1], rtol=1.0e-9, atol=1.0e-50
        )


def test_kernel_fixed_through_solver(cosmology: Cosmology) -> None:
    """The solver drives KERNEL_FIXED consistently with nc_xcor_compute().

    Both go through the same per-kernel closures, so this is exact agreement
    rather than agreement to sampling accuracy.
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
            rtol=1.0e-9,
            atol=1.0e-50,
        )


def _limber_pair(cosmology: Cosmology, l_limber: list[int]) -> list[Nc.XcorKernel]:
    """Overlapping top-hat kernels with the given per-kernel l_limber."""
    out = []

    for (z_lower, z_upper), l_lim in zip([(0.1, 0.3), (0.2, 0.4)], l_limber):
        kernel = Nc.XcorKernelClusterTophat(
            dist=cosmology.dist,
            powspec=cosmology.ps_ml,
            z_lower=z_lower,
            z_upper=z_upper,
            integrator=Ncm.SBesselIntegratorLevin.new(0, 8),
        )
        kernel.set_l_limber(l_lim)
        kernel.prepare(cosmology.cosmo)
        out.append(kernel)

    return out


@pytest.mark.parametrize("l_limber", [[0, 0], [0, -1], [-1, 0]])
def test_kernel_fixed_respects_l_limber(
    cosmology: Cosmology, l_limber: list[int]
) -> None:
    """Each kernel is evaluated in its own l_limber tier."""
    k1, k2 = _limber_pair(cosmology, l_limber)
    cosmo = cosmology.cosmo

    fixed = Nc.Xcor.new(cosmology.dist, cosmology.ps_ml, Nc.XcorMethod.KERNEL_FIXED)
    fixed.prepare(cosmo)
    cubature = Nc.Xcor.new(
        cosmology.dist, cosmology.ps_ml, Nc.XcorMethod.KERNEL_CUBATURE
    )
    cubature.prepare(cosmo)

    got, expected = Ncm.Vector.new(1), Ncm.Vector.new(1)
    for ell in (10, 200):
        fixed.compute(k1, k2, cosmo, ell, ell, got)
        cubature.compute(k1, k2, cosmo, ell, ell, expected)

        assert expected.get(0) != 0.0
        assert_allclose(got.get(0), expected.get(0), rtol=5.0e-3)


def test_kernel_fixed_all_limber_needs_no_integrator(cosmology: Cosmology) -> None:
    """An all-Limber block performs no Bessel integral, so it needs no integrator."""
    kernels = _limber_pair(cosmology, [0, 0])
    cosmo = cosmology.cosmo

    fixed = Nc.Xcor.new(cosmology.dist, cosmology.ps_ml, Nc.XcorMethod.KERNEL_FIXED)
    fixed.prepare(cosmo)

    solver = Nc.XcorSolver.new()
    ids = [solver.register_kernel(k) for k in kernels]
    solver.request_cl(ids[0], ids[1], 10, 13)
    solver.plan_blocks(8)
    solver.solve(fixed, cosmo)

    expected = Ncm.Vector.new(4)
    fixed.compute(kernels[0], kernels[1], cosmo, 10, 13, expected)

    assert_allclose(
        np.array(solver.get_result(0).dup_array()),
        np.array(expected.dup_array()),
        rtol=1.0e-9,
    )

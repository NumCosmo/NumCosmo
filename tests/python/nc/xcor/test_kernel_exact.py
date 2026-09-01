#!/usr/bin/env python
#
# test_kernel_exact.py
#
# Mon Aug 11 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_kernel_exact.py
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

"""NC_XCOR_METHOD_KERNEL_EXACT: exact GL(5) on per-pair knot unions.

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
    integrand = kernel.get_eval(cosmology.cosmo, 0, Nc.XcorKernelClosure.SPLINE)

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
    sets = [
        _knots(k.get_eval(cosmology.cosmo, 0, Nc.XcorKernelClosure.SPLINE))
        for k in kernels
    ]

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

    igs = [k.get_eval(cosmo, 0, Nc.XcorKernelClosure.SPLINE) for k in kernels]
    edges = [_knots(ig) for ig in igs]

    fixed = Nc.Xcor.new(cosmology.dist, cosmology.ps_ml, Nc.XcorMethod.KERNEL_EXACT)
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


def test_kernel_exact_matches_adaptive_methods(cosmology: Cosmology) -> None:
    """KERNEL_EXACT reproduces the adaptive kernel-space methods, with no tolerance."""
    kernels = _kernels(cosmology)
    cosmo = cosmology.cosmo
    nbins = len(kernels)

    fixed = Nc.Xcor.new(cosmology.dist, cosmology.ps_ml, Nc.XcorMethod.KERNEL_EXACT)
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


def test_kernel_exact_over_ell_block(cosmology: Cosmology) -> None:
    """KERNEL_EXACT handles a multipole block, on the same range as cubature.

    Both methods now intersect the two integrands' fitted domains, so the
    separated-bin cross spectrum agrees far into the tail -- unlike the earlier
    shared-abscissa implementation, which masked each component to a narrower
    per-component support and clipped exactly this case.
    """
    kernels = _kernels(cosmology)
    cosmo = cosmology.cosmo
    lmin, lmax = 2, 5
    n_l = lmax - lmin + 1

    fixed = Nc.Xcor.new(cosmology.dist, cosmology.ps_ml, Nc.XcorMethod.KERNEL_EXACT)
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


def test_kernel_exact_batches_wide_ell_range(cosmology: Cosmology) -> None:
    """KERNEL_EXACT batches a range wider than NC_XCOR_KERNEL_MAX_ELL_BLOCK.

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

    fixed = Nc.Xcor.new(cosmology.dist, cosmology.ps_ml, Nc.XcorMethod.KERNEL_EXACT)
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


def test_kernel_exact_through_solver(cosmology: Cosmology) -> None:
    """The solver drives KERNEL_EXACT consistently with nc_xcor_compute().

    Both go through the same per-kernel closures, so this is exact agreement
    rather than agreement to sampling accuracy.
    """
    kernels = _kernels(cosmology)
    cosmo = cosmology.cosmo

    fixed = Nc.Xcor.new(cosmology.dist, cosmology.ps_ml, Nc.XcorMethod.KERNEL_EXACT)
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
def test_kernel_exact_respects_l_limber(
    cosmology: Cosmology, l_limber: list[int]
) -> None:
    """Each kernel is evaluated in its own l_limber tier."""
    k1, k2 = _limber_pair(cosmology, l_limber)
    cosmo = cosmology.cosmo

    fixed = Nc.Xcor.new(cosmology.dist, cosmology.ps_ml, Nc.XcorMethod.KERNEL_EXACT)
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


def test_kernel_exact_all_limber_needs_no_integrator(cosmology: Cosmology) -> None:
    """An all-Limber block performs no Bessel integral, so it needs no integrator."""
    kernels = _limber_pair(cosmology, [0, 0])
    cosmo = cosmology.cosmo

    fixed = Nc.Xcor.new(cosmology.dist, cosmology.ps_ml, Nc.XcorMethod.KERNEL_EXACT)
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


def test_error_estimate_is_small_on_an_auto_spectrum(cosmology: Cosmology) -> None:
    """An auto spectrum has a positive integrand, so there is nothing to amplify.

    The estimate must land far below one -- these C_ell have significant digits
    -- which is the calibration point every other case here is read against.
    """
    kernel = _kernels(cosmology)[0]
    cosmo = cosmology.cosmo
    lmin, lmax = 2, 9
    nell = lmax - lmin + 1

    exact = Nc.Xcor.new(cosmology.dist, cosmology.ps_ml, Nc.XcorMethod.KERNEL_EXACT)
    exact.prepare(cosmo)

    vp, vp_err = Ncm.Vector.new(nell), Ncm.Vector.new(nell)
    exact.compute_full(kernel, None, cosmo, lmin, lmax, vp, vp_err)

    cl = np.array(vp.dup_array())
    rel = np.abs(np.array(vp_err.dup_array()) / cl)

    assert np.all(cl > 0.0)
    # Both halves of the fit criterion are relative to something, so the
    # estimate sits a fixed few orders above the tolerances rather than at
    # them -- the peak-scaled floor is the larger of the two contributions.
    assert np.all(rel > kernel.get_reltol())
    assert np.all(rel < 1.0e-2)


def test_error_estimate_scales_with_the_fit_criterion(cosmology: Cosmology) -> None:
    """Both halves of the criterion enter linearly, so the estimate must too.

    This is what makes the estimate a lever rather than a label: tightening the
    kernels buys down the reported error in proportion.
    """
    cosmo = cosmology.cosmo
    lmin, lmax = 2, 9
    nell = lmax - lmin + 1

    exact = Nc.Xcor.new(cosmology.dist, cosmology.ps_ml, Nc.XcorMethod.KERNEL_EXACT)
    exact.prepare(cosmo)

    def worst_relative(tol):
        kernel = Nc.XcorKernelClusterTophat(
            dist=cosmology.dist,
            powspec=cosmology.ps_ml,
            z_lower=Z_BINS[0][0],
            z_upper=Z_BINS[0][1],
            reltol=tol,
            scaled_abstol=tol,
            integrator=Ncm.SBesselIntegratorLevin.new(0, 8),
        )
        kernel.set_l_limber(-1)
        kernel.prepare(cosmo)

        vp, vp_err = Ncm.Vector.new(nell), Ncm.Vector.new(nell)
        exact.compute_full(kernel, None, cosmo, lmin, lmax, vp, vp_err)

        return np.abs(np.array(vp_err.dup_array()) / np.array(vp.dup_array())).max()

    loose = worst_relative(1.0e-4)
    tight = worst_relative(1.0e-6)

    # Linear to within a factor of two; the closures are rebuilt on different
    # knots, so the two runs do not integrate quite the same function.
    assert_allclose(loose / tight, 100.0, rtol=1.0)


def test_error_estimate_grows_with_cancellation(cosmology: Cosmology) -> None:
    """The estimate must separate a well-conditioned pair from a cancelling one.

    This is the whole point of reporting it: the two C_ell vectors look equally
    respectable, and only the estimate says one of them has no digits left.
    """
    near, _, far = _kernels(cosmology)
    cosmo = cosmology.cosmo
    lmin, lmax = 2, 9
    nell = lmax - lmin + 1

    exact = Nc.Xcor.new(cosmology.dist, cosmology.ps_ml, Nc.XcorMethod.KERNEL_EXACT)
    exact.prepare(cosmo)

    def relative_estimate(k1, k2):
        vp, vp_err = Ncm.Vector.new(nell), Ncm.Vector.new(nell)
        exact.compute_full(k1, k2, cosmo, lmin, lmax, vp, vp_err)

        return np.abs(np.array(vp_err.dup_array()) / np.array(vp.dup_array()))

    auto = relative_estimate(near, near)
    cross = relative_estimate(near, far)

    # Against a distant bin the cross spectrum is built from tail against tail.
    # The separation is what is asserted: every multipole of the cross is two
    # orders worse than the worst of the auto, and the worst of them has no
    # digits left at all. A cellwise floor is deliberately not asserted -- the
    # achieved-residual estimate puts the two lowest multipoles at 0.36 and
    # 0.78, where the tolerance-only ceiling it replaced said 12 and 25.
    assert np.all(auto < 1.0e-2)
    assert np.all(cross > 100.0 * auto.max())
    assert cross.max() > 1.0


def test_error_estimate_is_nan_for_methods_that_do_not_provide_one(
    cosmology: Cosmology,
) -> None:
    """A method with no estimate must say so, not report zero error."""
    kernel = _kernels(cosmology)[0]
    cosmo = cosmology.cosmo
    nell = 3

    for method in (
        Nc.XcorMethod.KERNEL_CUBATURE,
        Nc.XcorMethod.KERNEL_GSL,
        Nc.XcorMethod.LIMBER_Z_CUBATURE,
    ):
        xcor = Nc.Xcor.new(cosmology.dist, cosmology.ps_ml, method)
        xcor.prepare(cosmo)

        vp, vp_err = Ncm.Vector.new(nell), Ncm.Vector.new(nell)
        xcor.compute_full(kernel, None, cosmo, 2, 4, vp, vp_err)

        assert np.all(np.isnan(np.array(vp_err.dup_array()))), method


def test_compute_full_without_an_error_vector_matches_compute(
    cosmology: Cosmology,
) -> None:
    """The estimate is opt-in, and asking for it must not perturb the C_ell."""
    k1, k2 = _kernels(cosmology)[:2]
    cosmo = cosmology.cosmo
    nell = 4

    exact = Nc.Xcor.new(cosmology.dist, cosmology.ps_ml, Nc.XcorMethod.KERNEL_EXACT)
    exact.prepare(cosmo)

    plain, with_err, discarded = (Ncm.Vector.new(nell) for _ in range(3))
    exact.compute(k1, k2, cosmo, 2, 5, plain)
    exact.compute_full(k1, k2, cosmo, 2, 5, with_err, None)
    exact.compute_full(k1, k2, cosmo, 2, 5, discarded, Ncm.Vector.new(nell))

    assert_allclose(np.array(with_err.dup_array()), np.array(plain.dup_array()))
    assert_allclose(np.array(discarded.dup_array()), np.array(plain.dup_array()))


def test_closure_records_the_residual_it_achieved(cosmology: Cosmology) -> None:
    """The closure reports the fit error it reached, per knot interval.

    This is what the error estimate is built on. The record is a matrix aligned
    with the knots -- row i owns the interval [k_i, k_i+1], so the last row is
    the trailing edge and carries no interval.
    """
    kernel = _kernels(cosmology)[0]
    integrand = kernel.get_eval_vectorized_full(
        cosmology.cosmo,
        2,
        9,
        Ncm.SBesselIntegratorLevin.new(0, 8),
        Nc.XcorKernelClosure.SPLINE,
    )

    residuals = integrand.peek_residuals()
    assert residuals is not None

    knots = _knots(integrand)
    assert residuals.nrows() == knots.size
    assert residuals.ncols() == 8

    values = np.array(residuals.dup_array()).reshape(residuals.nrows(), -1)
    assert np.all(np.isnan(values[-1]))
    assert np.all(np.isfinite(values[:-1]))
    assert np.all(values[:-1] >= 0.0)

    # It must be an achieved residual, not the tolerance that was requested.
    # The criterion refinement stops at is reltol * ||W(k)||_2 + a * W_max, so
    # the worst interval sits just under it -- that is where refinement quit --
    # while the typical one is orders below. A record echoing the tolerance
    # back would be flat instead.
    window = np.array([integrand.eval_array(k) for k in knots])
    criterion = kernel.get_reltol() * np.linalg.norm(window, axis=1).max()
    criterion += kernel.get_scaled_abstol() * np.abs(window).max()

    assert values[:-1].max() < criterion
    assert np.median(values[:-1]) < 0.05 * criterion


def test_tracking_off_leaves_no_record_and_a_looser_estimate(
    cosmology: Cosmology,
) -> None:
    """Turning the record off falls back to the tolerance-only ceiling."""
    cosmo = cosmology.cosmo
    lmin, lmax = 2, 9
    nell = lmax - lmin + 1

    exact = Nc.Xcor.new(cosmology.dist, cosmology.ps_ml, Nc.XcorMethod.KERNEL_EXACT)
    exact.prepare(cosmo)

    def estimate(track):
        kernel = _kernels(cosmology)[0]
        kernel.set_track_fit_residual(track)
        kernel.prepare(cosmo)

        integrand = kernel.get_eval_vectorized_full(
            cosmo,
            lmin,
            lmax,
            Ncm.SBesselIntegratorLevin.new(0, 8),
            Nc.XcorKernelClosure.SPLINE,
        )
        assert (integrand.peek_residuals() is not None) == track

        vp, vp_err = Ncm.Vector.new(nell), Ncm.Vector.new(nell)
        exact.compute_full(kernel, None, cosmo, lmin, lmax, vp, vp_err)

        return np.abs(np.array(vp_err.dup_array()) / np.array(vp.dup_array()))

    assert _kernels(cosmology)[0].get_track_fit_residual(), "on by default"

    achieved = estimate(True)
    ceiling = estimate(False)

    # Measured 12-858x against a reltol=1e-10 reference for this pair, against
    # 1.2-50x for the achieved residual. Asserted loosely: the point is that
    # the two are different objects, not the exact factor.
    assert np.all(achieved < ceiling)
    assert np.all(ceiling / achieved > 5.0)


def test_achieved_residual_estimate_still_bounds_the_true_error(
    cosmology: Cosmology,
) -> None:
    """Sharper must not mean wrong: the estimate still has to cover the truth.

    Measured against a closure built two orders tighter, which is the only
    reference available for a fit error -- the quadrature above it is exact.
    """
    cosmo = cosmology.cosmo
    lmin, lmax = 2, 9
    nell = lmax - lmin + 1

    exact = Nc.Xcor.new(cosmology.dist, cosmology.ps_ml, Nc.XcorMethod.KERNEL_EXACT)
    exact.prepare(cosmo)

    def auto_spectrum(reltol, scaled_abstol):
        kernel = Nc.XcorKernelClusterTophat(
            dist=cosmology.dist,
            powspec=cosmology.ps_ml,
            z_lower=Z_BINS[0][0],
            z_upper=Z_BINS[0][1],
            reltol=reltol,
            scaled_abstol=scaled_abstol,
            integrator=Ncm.SBesselIntegratorLevin.new(0, 8),
        )
        kernel.set_l_limber(-1)
        kernel.prepare(cosmo)

        vp, vp_err = Ncm.Vector.new(nell), Ncm.Vector.new(nell)
        exact.compute_full(kernel, None, cosmo, lmin, lmax, vp, vp_err)

        return np.array(vp.dup_array()), np.array(vp_err.dup_array())

    cl_ref, _ = auto_spectrum(1.0e-8, 1.0e-6)
    cl, err = auto_spectrum(1.0e-4, 1.0e-4)

    true_rel = np.abs(cl - cl_ref) / np.abs(cl_ref)
    est_rel = np.abs(err / cl)

    assert np.all(est_rel > true_rel)


def _closure(cosmology: Cosmology, closure_type, reltol=1.0e-4, scaled_abstol=1.0e-4):
    """A cluster top-hat closure in the requested representation."""
    kernel = Nc.XcorKernelClusterTophat(
        dist=cosmology.dist,
        powspec=cosmology.ps_ml,
        z_lower=Z_BINS[0][0],
        z_upper=Z_BINS[0][1],
        integrator=Ncm.SBesselIntegratorLevin.new(0, 8),
        reltol=reltol,
        scaled_abstol=scaled_abstol,
    )
    kernel.set_l_limber(-1)
    kernel.prepare(cosmology.cosmo)

    return kernel.get_eval_vectorized_full(
        cosmology.cosmo, 2, 9, Ncm.SBesselIntegratorLevin.new(0, 8), closure_type
    )


def test_chebyshev_closure_reports_a_spectral_representation(
    cosmology: Cosmology,
) -> None:
    """A Chebyshev closure carries coefficients; a spline one does not.

    The two are parallel implementations rather than a replacement: each reports what it has,
    and a caller that wants coefficients asks for them.
    """
    cheb = _closure(cosmology, Nc.XcorKernelClosure.CHEBYSHEV)
    spline = _closure(cosmology, Nc.XcorKernelClosure.SPLINE)

    n_panels = cheb.get_n_panels()
    assert n_panels > 0

    # Panels are contiguous and ascending, and together they tile the range.
    lo, hi = cheb.get_range()
    prev_b = lo
    for i in range(n_panels):
        coeffs, a, b = cheb.peek_panel(i)
        assert coeffs.nrows() == 8
        assert coeffs.ncols() > 1
        assert a == pytest.approx(prev_b)
        assert b > a
        prev_b = b
    assert prev_b == pytest.approx(hi)

    # The single-expansion accessor answers only when there is one panel, since
    # that is the case a caller working on coefficients can use directly.
    assert cheb.peek_spectral()[0] == (n_panels == 1)

    # A spline closure reports knots and no expansion; the Chebyshev one the
    # reverse. Neither is asked to pretend it has the other.
    assert cheb.peek_knots() is None
    assert spline.peek_knots() is not None
    assert not spline.peek_spectral()[0]
    assert spline.get_n_panels() == 0


def test_chebyshev_closure_agrees_with_the_spline_one(cosmology: Cosmology) -> None:
    """Same sampled function, same domain, so the two must describe one W.

    Both are graded against a spline closure built six orders tighter. At equal
    *requested* tolerance the two land in the same place -- panels converge to
    what was asked and no further, so the Chebyshev closure is comparable here
    rather than better. Where it wins is samples for that accuracy, and how
    cheaply it tightens: spectral convergence rather than h^4.
    """
    cheb = _closure(cosmology, Nc.XcorKernelClosure.CHEBYSHEV)
    spline = _closure(cosmology, Nc.XcorKernelClosure.SPLINE)
    ref = _closure(cosmology, Nc.XcorKernelClosure.SPLINE, 1.0e-10, 1.0e-7)

    lo = max(c.get_range()[0] for c in (cheb, spline, ref)) * 1.001
    hi = min(c.get_range()[1] for c in (cheb, spline, ref)) * 0.999
    ks = np.geomspace(lo, hi, 400)

    truth = np.array([ref.eval_array(k) for k in ks])
    got_cheb = np.array([cheb.eval_array(k) for k in ks])
    got_spline = np.array([spline.eval_array(k) for k in ks])
    peak = np.abs(truth).max(axis=0)

    err_cheb = np.abs(got_cheb - truth).max(axis=0) / peak
    err_spline = np.abs(got_spline - truth).max(axis=0) / peak

    # Both meet the tolerance they were asked for, on every multipole.
    assert np.all(err_cheb < 1.0e-3)
    assert np.all(err_spline < 1.0e-3)

    # And they describe the same function: the gap between them is no larger
    # than either one's own distance from the reference.
    assert np.all(
        np.abs(got_cheb - got_spline).max(axis=0) / peak
        < 10.0 * np.maximum(err_cheb, err_spline)
    )

    # The Chebyshev closure reaches that on fewer samples than the spline needs
    # knots -- the accuracy per evaluation is the point, evaluations being
    # radial solves.
    n_cheb = sum(cheb.peek_panel(i)[0].ncols() for i in range(cheb.get_n_panels()))
    assert n_cheb < spline.peek_knots().len()


def test_spline_closure_is_the_default(cosmology: Cosmology) -> None:
    """The representation every method is calibrated against stays the default.

    The choice belongs to NcXcor rather than to a kernel: KERNEL_EXACT
    integrates a pair on the common refinement of the two closures' panels and
    needs both to be of the same kind, so one setting governs every kernel in a
    computation and a mixed pair cannot be expressed.
    """
    xcor = Nc.Xcor.new(cosmology.dist, cosmology.ps_ml, Nc.XcorMethod.KERNEL_EXACT)

    assert xcor.get_closure_type() == Nc.XcorKernelClosure.SPLINE

    xcor.set_closure_type(Nc.XcorKernelClosure.CHEBYSHEV)
    assert xcor.get_closure_type() == Nc.XcorKernelClosure.CHEBYSHEV


def test_limber_multipoles_keep_the_spline_closure(cosmology: Cosmology) -> None:
    """Chebyshev is refused where the window is not entire in k.

    Under Limber a multipole is supported on its own band and zero outside it,
    so the block's window carries a step per multipole. A Chebyshev series
    converges on this kernel because W_l(k) is entire in k; a step is not, so
    the expansion cannot converge and a panel splitter would bisect until it
    gave up -- which is exactly what this did before the closure builder
    started refusing.
    """
    integrands = []

    for closure_type in (
        Nc.XcorKernelClosure.SPLINE,
        Nc.XcorKernelClosure.CHEBYSHEV,
    ):
        kernel = Nc.XcorKernelClusterTophat(
            dist=cosmology.dist,
            powspec=cosmology.ps_ml,
            z_lower=Z_BINS[0][0],
            z_upper=Z_BINS[0][1],
            integrator=Ncm.SBesselIntegratorLevin.new(0, 8),
            reltol=1.0e-4,
            scaled_abstol=1.0e-4,
        )
        # Every multipole in the block falls on the Limber side.
        kernel.set_l_limber(0)
        kernel.prepare(cosmology.cosmo)

        integrands.append(
            kernel.get_eval_vectorized_full(
                cosmology.cosmo,
                2,
                9,
                Ncm.SBesselIntegratorLevin.new(0, 8),
                closure_type,
            )
        )

    for integrand in integrands:
        assert integrand.peek_knots() is not None, "expected a spline closure"
        assert integrand.get_n_panels() == 0

    # The same closure either way, not merely the same kind of closure.
    spline, cheb = integrands
    assert spline.get_range() == cheb.get_range()

    lo, hi = spline.get_range()
    for k in np.geomspace(lo * 1.001, hi * 0.999, 32):
        assert_allclose(cheb.eval_array(k), spline.eval_array(k), rtol=1.0e-14)


def test_spectral_pair_is_integrated_exactly(cosmology: Cosmology) -> None:
    """KERNEL_EXACT integrates Chebyshev closures on their merged panel edges.

    Two closures on the common refinement of their panels are each a single
    polynomial per cell, so the outer integral is a bilinear form in the
    coefficients rather than a quadrature. That is a second, independent route
    to the same C_ell as feeding the same closures to an adaptive rule, which
    is what makes the agreement below worth asserting: different arithmetic,
    same answer.
    """
    cosmo = cosmology.cosmo
    lmin, lmax = 2, 9
    nell = lmax - lmin + 1

    def kernels(tol=1.0e-4):
        out = []
        for z_lower, z_upper in Z_BINS[:2]:
            kernel = Nc.XcorKernelClusterTophat(
                dist=cosmology.dist,
                powspec=cosmology.ps_ml,
                z_lower=z_lower,
                z_upper=z_upper,
                integrator=Ncm.SBesselIntegratorLevin.new(0, 8),
                reltol=tol,
                scaled_abstol=tol,
            )
            kernel.set_l_limber(-1)
            kernel.prepare(cosmo)
            out.append(kernel)
        return out

    def compute(method, ks, auto, closure_type):
        xcor = Nc.Xcor.new(cosmology.dist, cosmology.ps_ml, method)
        xcor.set_closure_type(closure_type)
        xcor.prepare(cosmo)
        vp = Ncm.Vector.new(nell)
        xcor.compute(ks[0], None if auto else ks[1], cosmo, lmin, lmax, vp)
        return np.array(vp.dup_array())

    # One set of kernels: the representation is now the computation's choice,
    # so the same kernels serve both routes.
    ks = kernels()
    reference = kernels(1.0e-6)
    cheb_t = Nc.XcorKernelClosure.CHEBYSHEV
    spline_t = Nc.XcorKernelClosure.SPLINE

    for auto in (True, False):
        exact = compute(Nc.XcorMethod.KERNEL_EXACT, ks, auto, cheb_t)
        cubature = compute(Nc.XcorMethod.KERNEL_CUBATURE, ks, auto, cheb_t)

        assert np.all(np.isfinite(exact))

        # Same closures, unrelated arithmetic: a coefficient-space bilinear
        # form against an adaptive rule in k. The agreement is limited by
        # cubature, not by the exact route -- NcXcor:reltol is 1e-6, and the
        # worst multipole here lands at 1.1e-6, on the smallest and most
        # strongly cancelling C_ell of the block.
        assert_allclose(exact, cubature, rtol=1.0e-5)

        # And it is the closer of the two to a spline closure built two orders
        # tighter. Comparing against the spline at *this* tolerance would be
        # comparing against the less accurate answer: on a cancelling cross
        # spectrum the spline at 1e-4 is 8% out by l = 9, which is the error
        # this representation exists to remove.
        truth = compute(Nc.XcorMethod.KERNEL_EXACT, reference, auto, spline_t)
        err_exact = np.abs(exact / truth - 1.0)
        err_spline = np.abs(
            compute(Nc.XcorMethod.KERNEL_EXACT, ks, auto, spline_t) / truth - 1.0
        )

        # Measured against a spline closure at 1e-8: at this tolerance the
        # Chebyshev route is 11x closer on the auto spectrum and 3.6x on the
        # cross, and at 1e-6 it is 5400x and 15x. No absolute bound is asserted
        # on the cross -- at 1e-4 the cancellation dominates for both routes,
        # putting them at 3.5e-2 and 1.3e-1 respectively, so a threshold there
        # would be testing the tolerance rather than the representation.
        assert err_exact.max() < err_spline.max()


def test_panel_order_cap_is_tunable(cosmology: Cosmology) -> None:
    """The cap is a fitted heuristic, so a caller has to be able to sweep it.

    Its default came from a sweep on two kernel families at one multipole range;
    the optimum depends on how a window's phase is spread over its domain, which
    belongs to the kernel rather than to the library. It is a property and not a
    compile-time constant precisely so that assumption can be tested.
    """

    def closure(cap):
        kernel = Nc.XcorKernelClusterTophat(
            dist=cosmology.dist,
            powspec=cosmology.ps_ml,
            z_lower=Z_BINS[0][0],
            z_upper=Z_BINS[0][1],
            integrator=Ncm.SBesselIntegratorLevin.new(0, 8),
            reltol=1.0e-4,
            scaled_abstol=1.0e-4,
            panel_order_cap=cap,
        )
        kernel.set_l_limber(-1)
        kernel.prepare(cosmology.cosmo)

        return kernel, kernel.get_eval_vectorized_full(
            cosmology.cosmo,
            2,
            9,
            Ncm.SBesselIntegratorLevin.new(0, 8),
            Nc.XcorKernelClosure.CHEBYSHEV,
        )

    default_kernel, default_closure = closure(0)
    assert default_kernel.get_panel_order_cap() == 0

    # Zero selects the default, so it must land exactly where naming it does.
    _, explicit = closure(5)
    assert default_closure.get_n_panels() == explicit.get_n_panels()

    # A lower cap gives more, smaller panels; a higher one fewer, larger.
    _, coarse = closure(7)
    _, fine = closure(4)
    assert fine.get_n_panels() > default_closure.get_n_panels()
    assert coarse.get_n_panels() < default_closure.get_n_panels()

    for i in range(coarse.get_n_panels()):
        assert coarse.peek_panel(i)[0].ncols() <= (1 << 7) + 1

    # And they describe one function: the cap changes how the domain is carved,
    # not what is being fitted. Compared against the block's own peak rather
    # than pointwise -- in the tail both are near zero and a relative
    # comparison there measures the smaller of two small numbers.
    lo, hi = default_closure.get_range()
    ks = np.geomspace(lo * 1.001, hi * 0.999, 64)
    got_fine = np.array([fine.eval_array(k) for k in ks])
    got_coarse = np.array([coarse.eval_array(k) for k in ks])
    peak = np.abs(got_fine).max()

    assert np.abs(got_fine - got_coarse).max() < 1.0e-4 * peak

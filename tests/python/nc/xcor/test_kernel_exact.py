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

"""NC_XCOR_METHOD_KERNEL_EXACT: exact integrals of pairs of Chebyshev closures.

Each kernel is sampled independently -- the closures NcXcorSolver caches --
so a pair's two closures have different panels. On the common refinement of
the two sets of edges both are polynomials, so k^2 W_i W_j is integrated
exactly as a bilinear form in their coefficients.
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
        # Library default peak-epsilon, deliberately: this used to override it
        # to 1e-8 and no assertion here needed it. peak-epsilon floors W(k)
        # against its own peak, but the C_ell integrand is k^2 W_a W_b, so the
        # floor enters squared -- the 1e-4 default is already 1e-12 there. Every
        # test in this file passes at 1e-8 (14.3 s), 1e-6 (5.5 s) and the
        # default (1.5 s), so the override bought nothing and cost 10x.
        kernel.prepare(cosmology.cosmo)
        out.append(kernel)

    return out


def test_kernel_exact_batches_wide_ell_range(cosmology: Cosmology) -> None:
    """KERNEL_EXACT batches a range wider than NC_XCOR_KERNEL_MAX_ELL_BLOCK.

    A single k-space closure is capped at 64 multipoles, so the range is sliced
    into NcXcor:ell-batch-size sub-blocks. The batching is what makes the request
    legal, so the result must also be independent of where the batch boundaries
    fall.
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


@pytest.mark.parametrize(
    "closure_type,bounds",
    [
        (Nc.XcorKernelClosure.CHEBYSHEV, (100.0, 1.0e4)),
    ],
)
def test_error_estimate_scales_with_the_fit_criterion(
    cosmology: Cosmology, closure_type, bounds
) -> None:
    """Tightening the kernels must reduce the reported error.

    This is what makes the estimate a lever rather than a label. The factor by
    which it falls is the representation's own convergence rate. The Chebyshev
    closure answers a tighter request with a higher order, and converges
    spectrally in that order, so a request 100x tighter reduces the estimate by
    more than 100x: 205x here.
    """
    cosmo = cosmology.cosmo
    lmin, lmax = 2, 9
    nell = lmax - lmin + 1

    exact = Nc.Xcor.new(cosmology.dist, cosmology.ps_ml, Nc.XcorMethod.KERNEL_EXACT)
    exact.set_closure_type(closure_type)
    exact.prepare(cosmo)

    def worst_relative(tol):
        kernel = Nc.XcorKernelClusterTophat(
            dist=cosmology.dist,
            powspec=cosmology.ps_ml,
            z_lower=Z_BINS[0][0],
            z_upper=Z_BINS[0][1],
            reltol=tol,
            peak_epsilon=tol,
            integrator=Ncm.SBesselIntegratorLevin.new(0, 8),
        )
        kernel.set_l_limber(-1)
        kernel.prepare(cosmo)

        vp, vp_err = Ncm.Vector.new(nell), Ncm.Vector.new(nell)
        exact.compute_full(kernel, None, cosmo, lmin, lmax, vp, vp_err)

        return np.abs(np.array(vp_err.dup_array()) / np.array(vp.dup_array())).max()

    loose = worst_relative(1.0e-4)
    tight = worst_relative(1.0e-6)
    lo, hi = bounds

    assert lo < loose / tight < hi


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

    def auto_spectrum(reltol, peak_epsilon):
        kernel = Nc.XcorKernelClusterTophat(
            dist=cosmology.dist,
            powspec=cosmology.ps_ml,
            z_lower=Z_BINS[0][0],
            z_upper=Z_BINS[0][1],
            reltol=reltol,
            peak_epsilon=peak_epsilon,
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


def _closure(cosmology: Cosmology, closure_type, reltol=1.0e-4, peak_epsilon=1.0e-4):
    """A cluster top-hat closure in the requested representation."""
    kernel = Nc.XcorKernelClusterTophat(
        dist=cosmology.dist,
        powspec=cosmology.ps_ml,
        z_lower=Z_BINS[0][0],
        z_upper=Z_BINS[0][1],
        integrator=Ncm.SBesselIntegratorLevin.new(0, 8),
        reltol=reltol,
        peak_epsilon=peak_epsilon,
    )
    kernel.set_l_limber(-1)
    kernel.prepare(cosmology.cosmo)

    return kernel.get_eval_vectorized_full(
        cosmology.cosmo, 2, 9, Ncm.SBesselIntegratorLevin.new(0, 8), closure_type
    )


def test_chebyshev_closure_reports_a_spectral_representation(
    cosmology: Cosmology,
) -> None:
    """A Chebyshev closure carries its panels: one row per multipole, contiguous, tiling
    the closure's range."""
    cheb = _closure(cosmology, Nc.XcorKernelClosure.CHEBYSHEV)

    panels = cheb.peek_panels()
    n_panels = panels.get_n_panels()
    assert n_panels > 0

    # Panels are contiguous and ascending, and together they tile the range.
    lo, hi = cheb.get_range()
    prev_b = lo
    for i in range(n_panels):
        coeffs, a, b = panels.get_panel(i)
        assert coeffs.nrows() == 8
        assert coeffs.ncols() > 1
        assert a == pytest.approx(prev_b)
        assert b > a
        prev_b = b
    assert prev_b == pytest.approx(hi)


def test_chebyshev_closure_is_the_default(cosmology: Cosmology) -> None:
    """The more accurate representation is the default one.

    The choice belongs to NcXcor rather than to a kernel because it is the
    computation that gets switched over: the two are alternative fits to the
    same sampled function, and comparing them only means something when every
    kernel in a run uses the same one.
    """
    xcor = Nc.Xcor.new(cosmology.dist, cosmology.ps_ml, Nc.XcorMethod.KERNEL_EXACT)

    assert xcor.get_closure_type() == Nc.XcorKernelClosure.CHEBYSHEV


def test_limber_multipoles_take_the_requested_closure(cosmology: Cosmology) -> None:
    """Under Limber the Chebyshev closure is built in u = k / nu.

    A Limber window is zero outside the multipole's band in k. In u the band is
    the same for every multipole of the block, with the band edges as panel cuts
    and membership decided per panel, so every panel is smooth and the series
    converges.
    """

    def integrand(closure_type):
        kernel = Nc.XcorKernelClusterTophat(
            dist=cosmology.dist,
            powspec=cosmology.ps_ml,
            z_lower=Z_BINS[0][0],
            z_upper=Z_BINS[0][1],
            integrator=Ncm.SBesselIntegratorLevin.new(0, 8),
            reltol=1.0e-4,
            peak_epsilon=1.0e-4,
        )
        # Every multipole in the block falls on the Limber side.
        kernel.set_l_limber(0)
        kernel.prepare(cosmology.cosmo)

        return kernel.get_eval_vectorized_full(
            cosmology.cosmo, 2, 9, Ncm.SBesselIntegratorLevin.new(0, 8), closure_type
        )

    cheb = integrand(Nc.XcorKernelClosure.CHEBYSHEV)

    panels = cheb.peek_panels()
    assert panels.get_n_panels() > 0

    # The Chebyshev closure is in u = k / nu, where every multipole has the same
    # band, so each multipole's band edges in k are the panel edges times its
    # scale.
    edges = np.array(
        [panels.get_panel(0)[1]]
        + [panels.get_panel(i)[2] for i in range(panels.get_n_panels())]
    )
    for i in range(cheb.get_len()):
        k_min, k_max = cheb.get_range_comp(i)
        scale = cheb.get_scale(i)
        assert scale == 2 + i + 0.5
        assert np.isclose(edges * scale, k_min, rtol=1.0e-12).any()
        assert np.isclose(edges * scale, k_max, rtol=1.0e-12).any()


def test_spectral_path_reports_an_error_estimate(cosmology: Cosmology) -> None:
    """The panel route carries vp_err, and the estimate still covers the truth.

    The estimate is the closures' own fit error propagated through the pair's
    conditioning, a property of the closures rather than of the quadrature
    above them, and it bounds the deviation from a tighter run.
    """
    cosmo = cosmology.cosmo
    lmin, lmax = 2, 9
    nell = lmax - lmin + 1

    def auto_spectrum(reltol, peak_epsilon):
        kernel = Nc.XcorKernelClusterTophat(
            dist=cosmology.dist,
            powspec=cosmology.ps_ml,
            z_lower=Z_BINS[0][0],
            z_upper=Z_BINS[0][1],
            reltol=reltol,
            peak_epsilon=peak_epsilon,
            integrator=Ncm.SBesselIntegratorLevin.new(0, 8),
        )
        kernel.set_l_limber(-1)
        kernel.prepare(cosmo)

        xcor = Nc.Xcor.new(cosmology.dist, cosmology.ps_ml, Nc.XcorMethod.KERNEL_EXACT)
        xcor.prepare(cosmo)
        assert xcor.get_closure_type() == Nc.XcorKernelClosure.CHEBYSHEV

        vp, vp_err = Ncm.Vector.new(nell), Ncm.Vector.new(nell)
        xcor.compute_full(kernel, None, cosmo, lmin, lmax, vp, vp_err)

        return np.array(vp.dup_array()), np.array(vp_err.dup_array())

    cl, err = auto_spectrum(1.0e-4, 1.0e-4)

    assert np.all(np.isfinite(err))
    assert np.all(err > 0.0)

    cl_ref, _ = auto_spectrum(1.0e-8, 1.0e-6)
    true_rel = np.abs(cl - cl_ref) / np.abs(cl_ref)

    assert np.all(np.abs(err / cl) > true_rel)


def test_panel_level_cap_is_tunable(cosmology: Cosmology) -> None:
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
            peak_epsilon=1.0e-4,
        )
        if cap is not None:
            kernel.peek_panels_builder().set_level_cap(cap)
        kernel.set_l_limber(-1)
        kernel.prepare(cosmology.cosmo)

        return kernel, kernel.get_eval_vectorized_full(
            cosmology.cosmo,
            2,
            9,
            Ncm.SBesselIntegratorLevin.new(0, 8),
            Nc.XcorKernelClosure.CHEBYSHEV,
        )

    default_kernel, default_closure = closure(None)
    assert default_kernel.peek_panels_builder().get_level_cap() == 5
    default_panels = default_closure.peek_panels()

    # A lower cap gives more, smaller panels; a higher one fewer, larger.
    _, coarse = closure(7)
    _, fine = closure(4)
    coarse_panels = coarse.peek_panels()
    assert fine.peek_panels().get_n_panels() > default_panels.get_n_panels()
    assert coarse_panels.get_n_panels() < default_panels.get_n_panels()

    for i in range(coarse_panels.get_n_panels()):
        assert coarse_panels.get_panel(i)[0].ncols() <= (1 << 7) + 1

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

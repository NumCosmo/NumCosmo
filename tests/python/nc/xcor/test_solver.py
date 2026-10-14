#!/usr/bin/env python
#
# test_solver.py
#
# Thu Aug 06 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_solver.py
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

"""Tests for NcXcorSolver.solve(): the per-block shared-closure kernel-space path and the
direct-delegation fallback for every other method.

Registration, requests and the ell-block planner are checked in C,
tests/c/nc/xcor/test_nc_xcor_solver.c.
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


def test_solve_kernel_exact_matches_compute(
    cosmology: Cosmology,
    kernel_tsz: Nc.XcorKernel,
    kernel_cmb_lens: Nc.XcorKernel,
    kernel_gal_bin0: Nc.XcorKernel,
) -> None:
    """solve()'s KERNEL_EXACT shared-closure path matches nc_xcor_compute().

    Three kernels, three requests all in one 8-ell block: an auto (A), a
    cross reusing A (A-B), and a cross reusing B (B-C) -- exercises the
    per-block integrand cache being shared across multiple requests for the
    same kernel, not just a single pair.
    """
    lmin, lmax = 20, 27  # matches Nc.Xcor's default ell_batch_size (8)
    kernels = [kernel_tsz, kernel_cmb_lens, kernel_gal_bin0]

    original_l_limber = [k.get_l_limber() for k in kernels]
    original_integrator = [k.peek_integrator() for k in kernels]

    try:
        for k in kernels:
            k.props.integrator = Ncm.SBesselIntegratorLevin.new(0, 8)
            k.set_l_limber(0)  # tier 2: kernel-Limber
            k.prepare(cosmology.cosmo)

        xc = Nc.Xcor.new(cosmology.dist, cosmology.ps_ml, Nc.XcorMethod.KERNEL_EXACT)
        xc.prepare(cosmology.cosmo)

        solver = Nc.XcorSolver.new()
        id_a = solver.register_kernel(kernel_tsz)
        id_b = solver.register_kernel(kernel_cmb_lens)
        id_c = solver.register_kernel(kernel_gal_bin0)
        solver.request_cl(id_a, id_a, lmin, lmax)
        solver.request_cl(id_a, id_b, lmin, lmax)
        solver.request_cl(id_b, id_c, lmin, lmax)

        solver.plan_blocks(8)
        solver.solve(xc, cosmology.cosmo)

        pairs = [
            (kernel_tsz, kernel_tsz),
            (kernel_tsz, kernel_cmb_lens),
            (kernel_cmb_lens, kernel_gal_bin0),
        ]

        for i, (k1, k2) in enumerate(pairs):
            expected = Ncm.Vector.new(lmax - lmin + 1)
            xc.compute(k1, k2, cosmology.cosmo, lmin, lmax, expected)

            got = np.array(solver.get_result(i).dup_array())
            assert_allclose(
                got, np.array(expected.dup_array()), rtol=1.0e-6, atol=1.0e-50
            )
    finally:
        for k, l_limber, integrator in zip(
            kernels, original_l_limber, original_integrator
        ):
            k.set_l_limber(l_limber)
            k.props.integrator = integrator


def test_solve_spans_multiple_blocks(
    cosmology: Cosmology, kernel_tsz: Nc.XcorKernel
) -> None:
    """A request spanning more than one block is stitched together correctly."""
    lmin, lmax = 20, 39  # 20 ells: 3 blocks at the default block size (8)

    original_l_limber = kernel_tsz.get_l_limber()
    original_integrator = kernel_tsz.peek_integrator()

    try:
        kernel_tsz.props.integrator = Ncm.SBesselIntegratorLevin.new(0, 8)
        kernel_tsz.set_l_limber(0)
        kernel_tsz.prepare(cosmology.cosmo)

        xc = Nc.Xcor.new(cosmology.dist, cosmology.ps_ml, Nc.XcorMethod.KERNEL_EXACT)
        xc.prepare(cosmology.cosmo)

        solver = Nc.XcorSolver.new()
        kernel_id = solver.register_kernel(kernel_tsz)
        solver.request_cl(kernel_id, kernel_id, lmin, lmax)
        solver.plan_blocks(8)
        assert solver.get_n_blocks() > 1  # actually exercises stitching

        solver.solve(xc, cosmology.cosmo)

        expected = Ncm.Vector.new(lmax - lmin + 1)
        xc.compute(kernel_tsz, kernel_tsz, cosmology.cosmo, lmin, lmax, expected)

        got = np.array(solver.get_result(0).dup_array())
        assert_allclose(got, np.array(expected.dup_array()), rtol=1.0e-6, atol=1.0e-50)
    finally:
        kernel_tsz.set_l_limber(original_l_limber)
        kernel_tsz.props.integrator = original_integrator


def test_solve_block_parallel_stress(
    cosmology: Cosmology,
    kernel_tsz: Nc.XcorKernel,
    kernel_cmb_lens: Nc.XcorKernel,
    kernel_gal_bin0: Nc.XcorKernel,
) -> None:
    """solve() over many blocks/kernels/requests matches nc_xcor_compute(),
    exercising the OpenMP block-parallel path (dev-notes/xcor_ultralevin_batching_plan.md
    sec. 6.3) with real concurrency, not just a single-block sanity check.

    Each thread gets its own #NcmSerialize-duplicated clone of every
    registered kernel; this is the one test that would (nondeterministically)
    catch a mistake there -- a shared #NcmSBesselIntegratorLevin instance
    racing across two concurrently-processed blocks.
    """
    lmin, lmax = 2, 121  # 120 ells: 15 blocks at the default block size (8)
    kernels = [kernel_tsz, kernel_cmb_lens, kernel_gal_bin0]

    original_l_limber = [k.get_l_limber() for k in kernels]
    original_integrator = [k.peek_integrator() for k in kernels]

    try:
        for k in kernels:
            k.props.integrator = Ncm.SBesselIntegratorLevin.new(0, 8)
            k.set_l_limber(0)  # tier 2: kernel-Limber, cheap enough for 15 blocks
            k.prepare(cosmology.cosmo)

        xc = Nc.Xcor.new(cosmology.dist, cosmology.ps_ml, Nc.XcorMethod.KERNEL_EXACT)
        xc.prepare(cosmology.cosmo)

        solver = Nc.XcorSolver.new()
        id_a = solver.register_kernel(kernel_tsz)
        id_b = solver.register_kernel(kernel_cmb_lens)
        id_c = solver.register_kernel(kernel_gal_bin0)
        solver.request_cl(id_a, id_a, lmin, lmax)
        solver.request_cl(id_a, id_b, lmin, lmax)
        solver.request_cl(id_b, id_c, lmin, lmax)
        solver.request_cl(id_a, id_c, lmin, lmax)

        solver.plan_blocks(8)
        assert solver.get_n_blocks() >= 15
        solver.solve(xc, cosmology.cosmo)

        pairs = [
            (kernel_tsz, kernel_tsz),
            (kernel_tsz, kernel_cmb_lens),
            (kernel_cmb_lens, kernel_gal_bin0),
            (kernel_tsz, kernel_gal_bin0),
        ]

        for i, (k1, k2) in enumerate(pairs):
            expected = Ncm.Vector.new(lmax - lmin + 1)
            xc.compute(k1, k2, cosmology.cosmo, lmin, lmax, expected)

            got = np.array(solver.get_result(i).dup_array())
            assert_allclose(
                got, np.array(expected.dup_array()), rtol=1.0e-6, atol=1.0e-50
            )
    finally:
        for k, l_limber, integrator in zip(
            kernels, original_l_limber, original_integrator
        ):
            k.set_l_limber(l_limber)
            k.props.integrator = integrator


def test_solve_tier3_duplicated_kernel_shrinking_last_block(
    cosmology: Cosmology,
) -> None:
    """Regression test for a heap-buffer overflow in
    _ensure_operator_size() (numcosmo/ncm/specfunc/ncm_sbessel_ode_solver.c),
    found via solve()'s block-parallel path (dev-notes/xcor_ultralevin_batching_plan.md
    sec 6.3): each OpenMP thread evaluates its own NcmSerialize-duplicated
    kernel clone across every block it's assigned -- including transitions
    where the last block in an ell range is narrower than the others (e.g.
    lmax - lmin + 1 not a multiple of the block size). That combination
    (tier 3 / true non-Limber, a duplicated kernel's shared Levin
    integrator, and a block-width transition) corrupted the integrator's
    internal operator buffer; confirmed via valgrind pointing at the exact
    memcpy in _ensure_operator_size before the fix. Not reproducible at
    tier 2 (kernel-Limber, no real Levin/ODE-solve machinery involved) or
    with an ell range that's an exact multiple of the block size (no
    narrowing transition) -- both are exercised elsewhere in this file
    without incident.

    Runs solve() directly (not a manual duplicate, per the actual code
    path) with an ell range deliberately not a multiple of 8, so
    plan_blocks(8) produces a narrower final block.
    """
    lmin, lmax = 2, 24  # 23 ells: blocks [2,9], [10,17], [18,24] (7, not 8)

    integrator = Ncm.SBesselIntegratorLevin.new(lmin, lmax)
    integrator.set_reltol(1.0e-2)

    mu, sigma, z_len = 0.5, 0.08, 1000
    z_a = np.linspace(max(mu - 8 * sigma, 0.0), mu + 8 * sigma, z_len)
    nz_a = np.exp(-(((z_a - mu) / sigma) ** 2) / 2.0) / np.sqrt(2.0 * np.pi * sigma**2)

    kernel = Nc.XcorKernelWeakLensing(
        dist=cosmology.dist,
        powspec=cosmology.ps_ml,
        dndz=Ncm.SplineCubicNotaknot.new_full(
            Ncm.Vector.new_array(z_a), Ncm.Vector.new_array(nz_a), True
        ),
        nbar=3.0,
        intr_shear=7.0,
        integrator=integrator,
        # Matched to the integrator above: a closure cannot be fitted to more
        # precision than the samples carry, and the library refuses the pairing.
        reltol=1.0e-2,
        peak_epsilon=1.0e-2,
    )
    kernel.set_l_limber(-1)  # tier 3: true non-Limber
    kernel.prepare(cosmology.cosmo)

    xc = Nc.Xcor.new(cosmology.dist, cosmology.ps_ml, Nc.XcorMethod.KERNEL_EXACT)
    # The integrator above is deliberately loose to keep this test quick; the
    # outer integral cannot ask for more precision than that closure carries.
    xc.props.reltol = 1.0e-2
    xc.prepare(cosmology.cosmo)

    solver = Nc.XcorSolver.new()
    kernel_id = solver.register_kernel(kernel)
    solver.request_cl(kernel_id, kernel_id, lmin, lmax)
    solver.plan_blocks(8)
    assert solver.get_n_blocks() == 3

    solver.solve(xc, cosmology.cosmo)  # must not crash

    result = np.array(solver.get_result(0).dup_array())
    assert np.all(np.isfinite(result))
    assert np.all(result > 0.0)  # a weak-lensing auto-spectrum is positive


def _tier3_wl_kernel(cosmology: Cosmology, lmin: int, lmax: int) -> Nc.XcorKernel:
    """A deliberately loose tier-3 weak-lensing kernel, cheap enough to solve."""
    integrator = Ncm.SBesselIntegratorLevin.new(lmin, lmax)
    integrator.set_reltol(1.0e-2)

    mu, sigma, z_len = 0.5, 0.08, 1000
    z_a = np.linspace(max(mu - 8 * sigma, 0.0), mu + 8 * sigma, z_len)
    nz_a = np.exp(-(((z_a - mu) / sigma) ** 2) / 2.0) / np.sqrt(2.0 * np.pi * sigma**2)

    kernel = Nc.XcorKernelWeakLensing(
        dist=cosmology.dist,
        powspec=cosmology.ps_ml,
        dndz=Ncm.SplineCubicNotaknot.new_full(
            Ncm.Vector.new_array(z_a), Ncm.Vector.new_array(nz_a), True
        ),
        nbar=3.0,
        intr_shear=7.0,
        integrator=integrator,
        reltol=1.0e-2,
        peak_epsilon=1.0e-2,
    )
    kernel.set_l_limber(-1)
    kernel.prepare(cosmology.cosmo)
    return kernel


def test_peek_block_integrator_before_and_after_solve(cosmology: Cosmology) -> None:
    """Blocks acquire their pinned integrator on solve(), not before."""
    lmin, lmax = 2, 17
    kernel = _tier3_wl_kernel(cosmology, lmin, lmax)

    xc = Nc.Xcor.new(cosmology.dist, cosmology.ps_ml, Nc.XcorMethod.KERNEL_EXACT)
    xc.props.reltol = 1.0e-2
    xc.prepare(cosmology.cosmo)

    solver = Nc.XcorSolver.new()
    kernel_id = solver.register_kernel(kernel)
    solver.request_cl(kernel_id, kernel_id, lmin, lmax)
    solver.plan_blocks(8)
    assert solver.get_n_blocks() == 2

    # Nothing is pinned until solve() runs, and out-of-range indices are NULL.
    assert solver.peek_block_integrator(0) is None
    assert solver.peek_block_integrator(99) is None

    solver.solve(xc, cosmology.cosmo)
    first_cl = np.array(solver.get_result(0).dup_array())

    first = [solver.peek_block_integrator(b) for b in range(solver.get_n_blocks())]
    assert all(isinstance(sbi, Ncm.SBesselIntegrator) for sbi in first)

    # Each integrator is pinned to its own block's multipole range, and the
    # blocks do not share one.
    for sbi, b in zip(first, range(solver.get_n_blocks())):
        block_lmin, block_lmax = solver.get_block(b)
        ell_min, ell_max = sbi.get_ell_range()
        assert (ell_min, ell_max) == (block_lmin, block_lmax)
    assert first[0] is not first[1]

    # A second solve() reuses them rather than rebuilding: that reuse is what
    # keeps each block's ODE operators' QR factorisation alive across calls.
    solver.solve(xc, cosmology.cosmo)
    second = [solver.peek_block_integrator(b) for b in range(solver.get_n_blocks())]
    assert all(a is b for a, b in zip(first, second))

    # And it is still the same answer.
    assert_allclose(np.array(solver.get_result(0).dup_array()), first_cl, rtol=0.0)


def test_set_integrator_overrides_and_resets_the_cache(cosmology: Cosmology) -> None:
    """An explicitly set prototype is used, and replacing it drops the cache."""
    lmin, lmax = 2, 17
    kernel = _tier3_wl_kernel(cosmology, lmin, lmax)

    xc = Nc.Xcor.new(cosmology.dist, cosmology.ps_ml, Nc.XcorMethod.KERNEL_EXACT)
    xc.props.reltol = 1.0e-2
    xc.prepare(cosmology.cosmo)

    solver = Nc.XcorSolver.new()
    kernel_id = solver.register_kernel(kernel)
    solver.request_cl(kernel_id, kernel_id, lmin, lmax)
    solver.plan_blocks(8)

    proto = Ncm.SBesselIntegratorLevin.new(lmin, lmax)
    proto.set_reltol(1.0e-2)
    solver.set_integrator(proto)

    solver.solve(xc, cosmology.cosmo)
    pinned = solver.peek_block_integrator(0)
    assert pinned is not None
    # A clone of the prototype, not the prototype itself: blocks must not share
    # one integrator's mutable state.
    assert pinned is not proto
    assert pinned.get_reltol() == pytest.approx(proto.get_reltol())

    # Setting a new prototype invalidates every pinned integrator.
    solver.set_integrator(proto)
    assert solver.peek_block_integrator(0) is None

    # Clearing it falls back to a registered kernel's own integrator.
    solver.set_integrator(None)
    solver.solve(xc, cosmology.cosmo)
    assert solver.peek_block_integrator(0) is not None


def test_solve_accepts_a_kernel_without_a_levin_integrator(
    cosmology: Cosmology,
) -> None:
    """The reltol closure check only applies to Levin integrators.

    A kernel carrying a non-Levin integrator has no Chebyshev closure whose
    precision could be compared against NcXcor:reltol, so the check must pass
    it through rather than reject it or read the wrong tolerance off it.
    """
    lmin, lmax = 2, 9
    kernel = Nc.XcorKernelWeakLensing(
        dist=cosmology.dist,
        powspec=cosmology.ps_ml,
        dndz=Ncm.SplineCubicNotaknot.new_full(
            Ncm.Vector.new_array(np.linspace(0.1, 0.9, 200)),
            Ncm.Vector.new_array(
                np.exp(-(((np.linspace(0.1, 0.9, 200) - 0.5) / 0.08) ** 2) / 2.0)
            ),
            True,
        ),
        nbar=3.0,
        intr_shear=7.0,
        integrator=Ncm.SBesselIntegratorGL.new(lmin, lmax),
    )
    kernel.set_l_limber(0)  # tier 2: kernel-Limber, no Levin closure involved
    kernel.prepare(cosmology.cosmo)

    xc = Nc.Xcor.new(cosmology.dist, cosmology.ps_ml, Nc.XcorMethod.KERNEL_EXACT)
    xc.prepare(cosmology.cosmo)

    solver = Nc.XcorSolver.new()
    kernel_id = solver.register_kernel(kernel)
    solver.request_cl(kernel_id, kernel_id, lmin, lmax)
    solver.plan_blocks(8)
    solver.solve(xc, cosmology.cosmo)

    result = np.array(solver.get_result(0).dup_array())
    assert np.all(np.isfinite(result))
    assert np.all(result > 0.0)


def test_solver_drives_spectral_closures(cosmology: Cosmology) -> None:
    """The solver reaches the block integrator by its own route.

    NcXcorSolver calls _nc_xcor_kernel_integrate_block_exact() directly rather
    than through nc_xcor_compute(), so a representation wired up only in the
    latter is invisible here -- which is exactly how the spectral path shipped
    broken for the solver at first. Same spectra either way is the assertion.
    """
    cosmo = cosmology.cosmo
    z_bins = [(0.1, 0.2), (0.3, 0.4)]

    def solve():
        solver = Nc.XcorSolver.new()
        ids = []

        for z_lower, z_upper in z_bins:
            kernel = Nc.XcorKernelClusterTophat(
                dist=cosmology.dist,
                powspec=cosmology.ps_ml,
                z_lower=z_lower,
                z_upper=z_upper,
                integrator=Ncm.SBesselIntegratorLevin.new(0, 8),
                reltol=1.0e-4,
                peak_epsilon=1.0e-4,
            )
            kernel.set_l_limber(-1)
            kernel.prepare(cosmo)
            ids.append(solver.register_kernel(kernel))

        for i, id_i in enumerate(ids):
            for id_j in ids[i:]:
                solver.request_cl(id_i, id_j, 2, 17)

        solver.plan_blocks(8)
        solver.set_integrator(Ncm.SBesselIntegratorLevin.new(0, 8))

        xcor = Nc.Xcor.new(cosmology.dist, cosmology.ps_ml, Nc.XcorMethod.KERNEL_EXACT)
        xcor.prepare(cosmo)
        solver.solve(xcor, cosmo)

        return [
            np.array(solver.get_result(r).dup_array())
            for r in range(solver.get_n_requests())
        ]

    def compute_directly():
        kernels = []

        for z_lower, z_upper in z_bins:
            kernel = Nc.XcorKernelClusterTophat(
                dist=cosmology.dist,
                powspec=cosmology.ps_ml,
                z_lower=z_lower,
                z_upper=z_upper,
                integrator=Ncm.SBesselIntegratorLevin.new(0, 8),
                reltol=1.0e-4,
                peak_epsilon=1.0e-4,
            )
            kernel.set_l_limber(-1)
            kernel.prepare(cosmo)
            kernels.append(kernel)

        xcor = Nc.Xcor.new(cosmology.dist, cosmology.ps_ml, Nc.XcorMethod.KERNEL_EXACT)
        xcor.prepare(cosmo)

        out = []
        for i in range(len(kernels)):
            for j in range(i, len(kernels)):
                vp = Ncm.Vector.new(16)
                xcor.compute(
                    kernels[i],
                    None if i == j else kernels[j],
                    cosmo,
                    2,
                    17,
                    vp,
                )
                out.append(np.array(vp.dup_array()))
        return out

    through_solver = solve()
    direct = compute_directly()

    assert len(through_solver) == 3

    # The two entry points must reach the same integrator.
    for got, expected in zip(through_solver, direct):
        assert np.all(np.isfinite(got))
        assert_allclose(got, expected, rtol=1.0e-10)


def test_replan_with_the_same_blocks_keeps_the_integrators(
    cosmology: Cosmology,
) -> None:
    """A plan that reproduces the previous blocks keeps the per-block integrators.

    They carry the factorised operators, the expensive part of a cold solve, and
    are keyed by the block's ell range only. A plan with different blocks drops
    them.
    """
    kernel = _tier3_wl_kernel(cosmology, 2, 17)
    xc = Nc.Xcor.new(cosmology.dist, cosmology.ps_ml, Nc.XcorMethod.KERNEL_EXACT)
    xc.props.reltol = 1.0e-2
    xc.prepare(cosmology.cosmo)

    solver = Nc.XcorSolver.new()
    kid = solver.register_kernel(kernel)
    solver.request_cl(kid, kid, 2, 17)
    solver.plan_blocks(8)
    solver.solve(xc, cosmology.cosmo)
    before = [solver.peek_block_integrator(b) for b in range(solver.get_n_blocks())]
    assert all(b is not None for b in before)

    solver.clear_requests()
    solver.request_cl(kid, kid, 2, 17)
    solver.plan_blocks(8)
    after = [solver.peek_block_integrator(b) for b in range(solver.get_n_blocks())]
    assert all(a is b for a, b in zip(after, before))

    solver.clear_requests()
    solver.request_cl(kid, kid, 2, 17)
    solver.plan_blocks(4)
    assert solver.get_n_blocks() == 4
    assert all(solver.peek_block_integrator(b) is None for b in range(4))

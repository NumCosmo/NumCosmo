"""Tests for NumCosmo cross-correlation integration methods.

This module tests different integration backends and methods for computing
angular power spectra using the Xcor class. Tests include:
- Different Xcor integration methods (LIMBER_Z_GSL, LIMBER_Z_CUBATURE, KERNEL_EXACT)
- Kernel-level integrator types (Limber approximation vs Levin integration)
- Xcor property management (reltol, method switching)
- Cross-correlation computation with various kernel pairs

Created: 2026-04-17
Purpose: Refactored from test_py_xcor.py for better organization
"""

from typing import cast

import numpy as np
import pytest
from numpy.testing import assert_allclose

from numcosmo_py import Nc, Ncm
from numcosmo_py.cosmology import Cosmology
from xcor import cases_k_integral as cases

pytest_plugins = [
    "python.fixtures_xcor",
]

pytestmark = pytest.mark.xcor


# Xcor Integration Method Tests


# The Limber approximation is not a usable reference for the CMB ISW
# auto-spectrum at these multipoles: its kernel is broad and peaks at low k,
# where the same-chi-peak assumption fails outright. Measured at ell 20-27,
# every kernel-space method agrees with every other to ~1e-8 while sitting a
# factor 3.3-4.2 below tier 1. Comparing the two is a statement about Limber,
# not about the quadrature, so tier 3 checks the methods against each other
# for this kernel and not against tier 1.
_TIER3_LIMBER_INVALID_KERNELS = {"kernel_cmb_isw"}


@pytest.mark.parametrize(
    "k1_name,k2_name",
    [
        ("kernel_cmb_lens", "kernel_cmb_lens"),
        ("kernel_cmb_isw", "kernel_cmb_isw"),
        ("kernel_tsz", "kernel_tsz"),
        ("kernel_gal_bin0", "kernel_gal_bin0"),
        ("kernel_wl_bin0", "kernel_wl_bin0"),
        ("kernel_cmb_lens", "kernel_cmb_isw"),  # Cross-correlation
        ("kernel_gal_bin0", "kernel_wl_bin0"),  # Galaxy-WL cross
    ],
)
@pytest.mark.parametrize(
    "l_limber", [0, -1], ids=["tier2_kernel_limber", "tier3_non_limber"]
)
def test_xcor_kernel_methods(
    cosmology: Cosmology,
    k1_name: str,
    k2_name: str,
    l_limber: int,
    request: pytest.FixtureRequest,
) -> None:
    """The kernel-space method against the z-space Limber integral.

    Checks, for both tier 2 (l_limber=0, kernel-Limber) and tier 3
    (l_limber=-1, true non-Limber), that KERNEL_EXACT, with four multipoles per
    closure, converges to LIMBER_Z_CUBATURE (tier 1) at moderate/high ell:
    exactly, for tier 2
    (mathematically the same Limber approximation, just computed via k-space
    instead of z-space); within known Limber-approximation error, for tier 3.
    The tier-3 comparison against tier 1 is dropped for the ISW
    auto-spectrum, where Limber itself is not valid at these multipoles --
    see _TIER3_LIMBER_INVALID_KERNELS above.

    Kernel fixtures are module-cached and shared with other tests, so
    l_limber is restored to its default afterward.
    """
    k1 = cast(Nc.XcorKernel, request.getfixturevalue(k1_name))
    k2 = cast(Nc.XcorKernel, request.getfixturevalue(k2_name))
    involved = {k1_name, k2_name}
    check_tier1 = not (l_limber == -1 and involved <= _TIER3_LIMBER_INVALID_KERNELS)

    lmin = 20
    lmax = 27  # small, moderate-ell block: away from the ell=0,1 edge case
    # (weak lensing's exactly-zero spin-2 prefactor there) and cheap for
    # tier 3's real Levin/ODE-solve cost.

    # Kernel fixtures share ONE session-wide NcmSBesselIntegratorLevin
    # (fixtures_xcor.py's _get_integrator(), functools.cache'd across every
    # test file). This test is the first to exercise l_limber=-1 (tier 3,
    # real Levin/ODE-solve) against it; doing so left it in a state that
    # crashed a later, unrelated test in a different file
    # (ncm_function_sample_set_expand_domain -> malloc corruption, only
    # reproducible with the full test directory or test_kernel.py run after
    # this one). Swap in dedicated, test-local integrators instead of
    # touching the shared one, restoring the originals afterward.
    original_l_limber = {k1_name: k1.get_l_limber(), k2_name: k2.get_l_limber()}
    original_integrator = {k1_name: k1.peek_integrator(), k2_name: k2.peek_integrator()}

    try:
        k1.props.integrator = Ncm.SBesselIntegratorLevin.new(0, 8)
        k2.props.integrator = Ncm.SBesselIntegratorLevin.new(0, 8)
        k1.set_l_limber(l_limber)
        k2.set_l_limber(l_limber)
        k1.prepare(cosmology.cosmo)
        k2.prepare(cosmology.cosmo)

        n = lmax - lmin + 1
        vp_limber_z = Ncm.Vector.new(n)

        xcor_lz = Nc.Xcor.new(
            cosmology.dist, cosmology.ps_ml, Nc.XcorMethod.LIMBER_Z_CUBATURE
        )
        xcor_lz.prepare(cosmology.cosmo)
        xcor_lz.compute(k1, k2, cosmology.cosmo, lmin, lmax, vp_limber_z)
        limber_z = np.array(vp_limber_z.dup_array())

        # Convergence to tier 1: near-exact for tier 2, Limber-error-sized
        # for tier 3. "Near-exact" allows a little numerical slack (observed
        # up to ~3e-4 depending on kernel shape), not literal bit-agreement.
        # Tier 3's bound is loose enough to cover cross-correlations between
        # differently-shaped kernels (observed ~4% for galaxy-counts x weak-
        # lensing at ell~20-27 -- Limber's same-chi-peak assumption is more
        # violated for cross-correlations between different tracers than for
        # auto-correlations), while still trivially catching the original
        # 15-27-orders-of-magnitude bug this test guards against.
        tier1_rtol = 1.0e-3 if l_limber == 0 else 1.0e-1

        # ISW cross-spectra need more room at tier 2, for a reason that is not
        # numerical slack. A kernel's k-space closure is fitted only out to
        # where its own amplitude falls below NcXcorKernel:adaptive-epsilon, a
        # criterion each kernel applies to itself. An auto spectrum never pays
        # for that truncation -- its integrand carries W squared, so the cut
        # sits where the contribution is doubly suppressed -- while a cross
        # spectrum is suppressed there only once, by one kernel's tail against
        # the other kernel's un-suppressed values.
        #
        # Measured for CMB lensing x ISW at ell 20-27: 1.05e-3 low against
        # tier 1, flat under refinement of the closure tolerance (1e-4 through
        # 1e-6, so it is the domain and not the closure), and collapsing to
        # 1.5e-6 once adaptive-epsilon is tightened to 1e-8. No per-kernel
        # criterion can do better in principle, since a closure is built
        # without knowing what it will be paired with. The scale that makes
        # this harmless is the pair's own: these two are 7% correlated here,
        # so 1.05e-3 of C_x is 7e-5 of sqrt(C_isw C_lens), the amplitude that
        # sets the cross spectrum's sample variance.
        if (l_limber == 0) and (involved & {"kernel_cmb_isw"}) and (k1_name != k2_name):
            tier1_rtol = 2.0e-3

        vp_kernel_exact = Ncm.Vector.new(n)
        xcor_kf = Nc.Xcor.new(
            cosmology.dist, cosmology.ps_ml, Nc.XcorMethod.KERNEL_EXACT
        )
        xcor_kf.props.ell_batch_size = 4
        xcor_kf.prepare(cosmology.cosmo)
        xcor_kf.compute(k1, k2, cosmology.cosmo, lmin, lmax, vp_kernel_exact)

        if check_tier1:
            assert_allclose(
                np.array(vp_kernel_exact.dup_array()),
                limber_z,
                rtol=tier1_rtol,
                atol=1.0e-50,
            )
        else:
            # Limber is no reference here, so the quadrature is checked on its own:
            # KERNEL_EXACT on the first block's closures against
            # cases_k_integral.reference_cl on the same closures.
            block_lmax = lmin + 3
            isauto = k1_name == k2_name
            sbi = Ncm.SBesselIntegratorLevin.new(lmin, block_lmax)
            ia = k1.get_eval_vectorized_full(cosmology.cosmo, lmin, block_lmax, sbi)
            ib = (
                None
                if isauto
                else k2.get_eval_vectorized_full(cosmology.cosmo, lmin, block_lmax, sbi)
            )
            vp_block = Ncm.Vector.new(block_lmax - lmin + 1)
            xcor_kf.integrate_block(
                ia,
                ib,
                lmin,
                block_lmax,
                isauto,
                Nc.XcorMethod.KERNEL_EXACT,
                vp_block,
                None,
            )
            reference = cases.reference_cl(cosmology.cosmo.RH_Mpc(), ia, ib)
            peak = np.abs(reference.cl).max()
            dev = np.abs(np.array(vp_block.dup_array()) - reference.cl).max() / peak
            # Measured 1.7e-15 on the ISW auto at tier 3, the reference converged to
            # the same level.
            assert reference.worst_cell_move < 1.0e-12
            assert dev < 1.0e-12
    finally:
        k1.set_l_limber(original_l_limber[k1_name])
        k2.set_l_limber(original_l_limber[k2_name])
        k1.props.integrator = original_integrator[k1_name]
        k2.props.integrator = original_integrator[k2_name]


# Measured 2.5e-10 of the block's peak on the tier-2 pair (magnification bias 0.6) and
# 1.8e-11 on the tier-3 one, at adaptive-epsilon 1e-8; an order of headroom.
SINGLE_VS_BLOCK_TOL = 1.0e-8


@pytest.mark.parametrize(
    "k1_name,k2_name,l_limber",
    [
        ("kernel_gal_bin0", "kernel_wl_bin0", 0),
        ("kernel_gal_bin0", "kernel_gal_bin0", -1),
    ],
    ids=["tier2_gal_x_wl", "tier3_gal_auto"],
)
def test_single_ell_closures_match_block(
    cosmology: Cosmology,
    k1_name: str,
    k2_name: str,
    l_limber: int,
    request: pytest.FixtureRequest,
) -> None:
    """One multipole per closure against four per closure, on the same spectra.

    The library builds its closures per multipole block, and every other test uses
    blocks; this is the one place the single-multipole closures are checked against
    them. The two build different fits of the same W_ell(k), so they agree to the
    closures' own accuracy, not exactly; the deviation is taken over the block's
    peak.

    The domain is held fixed far out (adaptive-epsilon 1e-8): at the default 1e-5 the
    magnification term's long tail in u is cut at different places for blocks of
    different widths (9.7e-6 of peak on the tier-2 pair), which is the stopping rule's
    defect, parked in NumCosmoDevNotes deferred_work.md, not the fits this test compares.
    """
    k1 = cast(Nc.XcorKernel, request.getfixturevalue(k1_name))
    k2 = cast(Nc.XcorKernel, request.getfixturevalue(k2_name))
    lmin, lmax = 20, 27
    n = lmax - lmin + 1

    original_l_limber = {k1_name: k1.get_l_limber(), k2_name: k2.get_l_limber()}
    original_integrator = {k1_name: k1.peek_integrator(), k2_name: k2.peek_integrator()}
    original_adaptive_epsilon = {
        k1_name: k1.get_adaptive_epsilon(),
        k2_name: k2.get_adaptive_epsilon(),
    }

    try:
        k1.props.integrator = Ncm.SBesselIntegratorLevin.new(0, 8)
        k2.props.integrator = Ncm.SBesselIntegratorLevin.new(0, 8)
        k1.set_l_limber(l_limber)
        k2.set_l_limber(l_limber)
        k1.set_adaptive_epsilon(1.0e-8)
        k2.set_adaptive_epsilon(1.0e-8)
        k1.prepare(cosmology.cosmo)
        k2.prepare(cosmology.cosmo)

        cl = {}
        for batch in (1, 4):
            xcor = Nc.Xcor.new(
                cosmology.dist, cosmology.ps_ml, Nc.XcorMethod.KERNEL_EXACT
            )
            xcor.props.ell_batch_size = batch
            xcor.prepare(cosmology.cosmo)
            vp = Ncm.Vector.new(n)
            xcor.compute(k1, k2, cosmology.cosmo, lmin, lmax, vp)
            cl[batch] = np.array(vp.dup_array())

        dev = np.abs(cl[1] - cl[4]).max() / np.abs(cl[4]).max()

        assert dev < SINGLE_VS_BLOCK_TOL, (
            f"{k1_name} x {k2_name} (l_limber {l_limber}): single-multipole closures "
            f"deviate from the block's by {dev:.3e} of its peak"
        )
    finally:
        k1.set_l_limber(original_l_limber[k1_name])
        k2.set_l_limber(original_l_limber[k2_name])
        k1.props.integrator = original_integrator[k1_name]
        k2.props.integrator = original_integrator[k2_name]
        k1.set_adaptive_epsilon(original_adaptive_epsilon[k1_name])
        k2.set_adaptive_epsilon(original_adaptive_epsilon[k2_name])


# Kernel-Level Integrator Tests


# Cross-Correlation Specific Tests


@pytest.mark.parametrize(
    "bin_i,bin_j",
    [(0, 0), (0, 1), (1, 1), (2, 3), (4, 4)],
)
def test_xcor_gal_bins_auto_cross(
    cosmology: Cosmology,
    bin_i: int,
    bin_j: int,
    request: pytest.FixtureRequest,
) -> None:
    """Test galaxy bin auto- and cross-correlations.

    Validates that cross-correlations between different galaxy bins
    are computed correctly and satisfy basic consistency checks.
    """
    kernel_i = cast(Nc.XcorKernelGal, request.getfixturevalue(f"kernel_gal_bin{bin_i}"))
    kernel_j = cast(Nc.XcorKernelGal, request.getfixturevalue(f"kernel_gal_bin{bin_j}"))

    xcor = Nc.Xcor.new(
        cosmology.dist,
        cosmology.ps_ml,
        Nc.XcorMethod.LIMBER_Z_CUBATURE,
    )
    xcor.prepare(cosmology.cosmo)
    kernel_i.prepare(cosmology.cosmo)
    kernel_j.prepare(cosmology.cosmo)

    lmin = 2
    lmax = 50
    cl_vec = Ncm.Vector.new(lmax - lmin + 1)

    xcor.compute(kernel_i, kernel_j, cosmology.cosmo, lmin, lmax, cl_vec)

    cl_array = np.array(cl_vec.dup_array())

    # Auto-correlations should be positive
    if bin_i == bin_j:
        assert np.all(cl_array > 0)

    # Cross-correlations should be real and finite
    assert np.all(np.isfinite(cl_array))

    # Verify symmetry: C_ij should equal C_ji
    if bin_i != bin_j:
        cl_vec_swap = Ncm.Vector.new(lmax - lmin + 1)
        xcor.compute(kernel_j, kernel_i, cosmology.cosmo, lmin, lmax, cl_vec_swap)
        cl_array_swap = np.array(cl_vec_swap.dup_array())
        assert_allclose(cl_array, cl_array_swap, rtol=1.0e-10)


@pytest.mark.parametrize(
    "bin_i,bin_j",
    [(0, 0), (0, 1), (2, 3), (4, 4)],
)
def test_xcor_wl_bins_auto_cross(
    cosmology: Cosmology,
    bin_i: int,
    bin_j: int,
    request: pytest.FixtureRequest,
) -> None:
    """Test weak lensing bin auto- and cross-correlations.

    Validates cosmic shear correlations between different source bins.
    """
    kernel_i = cast(
        Nc.XcorKernelWeakLensing, request.getfixturevalue(f"kernel_wl_bin{bin_i}")
    )
    kernel_j = cast(
        Nc.XcorKernelWeakLensing, request.getfixturevalue(f"kernel_wl_bin{bin_j}")
    )

    xcor = Nc.Xcor.new(
        cosmology.dist,
        cosmology.ps_ml,
        Nc.XcorMethod.LIMBER_Z_CUBATURE,
    )
    xcor.prepare(cosmology.cosmo)
    kernel_i.prepare(cosmology.cosmo)
    kernel_j.prepare(cosmology.cosmo)

    lmin = 2
    lmax = 50
    cl_vec = Ncm.Vector.new(lmax - lmin + 1)

    xcor.compute(kernel_i, kernel_j, cosmology.cosmo, lmin, lmax, cl_vec)

    cl_array = np.array(cl_vec.dup_array())

    # Auto-correlations should be positive
    if bin_i == bin_j:
        assert np.all(cl_array > 0)

    # All values should be finite
    assert np.all(np.isfinite(cl_array))

    # Verify symmetry
    if bin_i != bin_j:
        cl_vec_swap = Ncm.Vector.new(lmax - lmin + 1)
        xcor.compute(kernel_j, kernel_i, cosmology.cosmo, lmin, lmax, cl_vec_swap)
        cl_array_swap = np.array(cl_vec_swap.dup_array())
        assert_allclose(cl_array, cl_array_swap, rtol=1.0e-10)

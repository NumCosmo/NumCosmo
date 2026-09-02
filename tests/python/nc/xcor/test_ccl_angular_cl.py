"""Tests for the CCL-facing non-Limber angular power spectrum."""

import numpy as np
import pytest
from numpy.testing import assert_allclose
from scipy.integrate import quad
from scipy.special import spherical_jn

pytest.importorskip("pyccl")
# flake8: noqa: E402
# pylint: disable=wrong-import-position
import pyccl

from numcosmo_py import Nc, Ncm
from numcosmo_py.ccl.nc_ccl import create_nc_obj, CCLParams
from numcosmo_py.ccl.two_point import (
    angular_cl,
    bessel_transform,
    bessel_transform_block,
    tracer_component_tables,
    tracer_kernel,
    compute_kernel,
)

# Pinned to one worker under --dist loadgroup: this file is one of the
# xcor lane's memory peaks, and an xdist worker is its own session, so
# without this its cost is paid once per worker rather than once.
pytestmark = [pytest.mark.ccl, pytest.mark.xcor, pytest.mark.xdist_group("ccl_angular")]

Ncm.cfg_init()


@pytest.fixture(name="ccl_cosmo", scope="module")
def fixture_ccl_cosmo():
    CCLParams.set_high_prec_params()
    return pyccl.Cosmology(
        Omega_c=0.25,
        Omega_b=0.05,
        h=0.7,
        n_s=0.96,
        sigma8=0.8,
        transfer_function="eisenstein_hu",
        matter_power_spectrum="linear",
    )


@pytest.fixture(name="cosmology", scope="module")
def fixture_cosmology(ccl_cosmo):
    return create_nc_obj(ccl_cosmo)


@pytest.fixture(name="tracer", scope="module")
def fixture_tracer(ccl_cosmo):
    z = np.linspace(0.0, 2.0, 512)
    nz = np.exp(-0.5 * ((z - 0.8) / 0.12) ** 2)
    return pyccl.NumberCountsTracer(
        ccl_cosmo, has_rsd=False, dndz=(z, nz), bias=(z, np.ones_like(z))
    )


@pytest.mark.parametrize("k", [0.01, 0.02])
def test_bessel_transform_matches_quadrature(cosmology, tracer, k):
    """The oscillatory solve reproduces direct quadrature of the same kernel.

    The reference must carry the growth factor too: bessel_transform folds D(z) into
    the kernel because CCL's transfer does not (it evaluates the 2-D P(k, a) instead).
    """
    ell = 32
    z_a, chi_a, W_a = compute_kernel(tracer, cosmology, float(ell), reltol=1.0e-8)
    gf = Nc.GrowthFunc.new()
    gf.prepare_if_needed(cosmology.cosmo)
    W_a = W_a * np.array([gf.eval(cosmology.cosmo, float(zz)) for zz in z_a])
    sp = Ncm.SplineBSpline.new_tol(1.0e-8, 0.0)
    sp.set_array(chi_a.tolist(), W_a.tolist(), True)

    got = bessel_transform(tracer, cosmology, ell, np.array([k]), reltol=1.0e-8)[0]
    ref = quad(
        lambda c: sp.eval(c) * spherical_jn(ell, k * c),
        chi_a[0],
        chi_a[-1],
        limit=4000,
        epsabs=1.0e-16,
        epsrel=1.0e-12,
    )[0]

    assert_allclose(got, ref, rtol=1.0e-5)


@pytest.mark.parametrize("ell", [2, 8, 32])
def test_angular_cl_matches_ccl_non_limber(cosmology, ccl_cosmo, tracer, ell):
    """Agreement with CCL's own non-Limber calculation.

    Note pyccl's `l_limber` is the multipole *beyond which Limber is used*, so
    `l_limber=100000` selects non-Limber and the default `-1` selects Limber
    everywhere. Getting this backwards makes a correct result look like Limber.
    """
    got = angular_cl(cosmology, tracer, tracer, np.array([ell]))[0]
    ref = pyccl.angular_cl(ccl_cosmo, tracer, tracer, np.array([ell]), l_limber=100000)[
        0
    ]

    assert_allclose(got, ref, rtol=3.0e-3)


def test_closure_tolerance_is_converged(cosmology, tracer):
    """The library default closure tolerance is enough: 100x tighter changes little.

    Measured difference 8e-6 on a lensing tracer at ell 2; the bound below is
    for a density tracer, whose closure is easier.
    """
    a = angular_cl(cosmology, tracer, tracer, np.array([2]))[0]
    b = angular_cl(cosmology, tracer, tracer, np.array([2]), kernel_reltol=1.0e-6)[0]

    assert_allclose(a, b, rtol=2.0e-5)


def test_tracer_terms_map_to_component_kinds(cosmology, ccl_cosmo):
    """Every term of a CCL tracer becomes one component of the matching kind."""
    z = np.linspace(0.0, 2.0, 256)
    nz = np.exp(-0.5 * ((z - 0.8) / 0.12) ** 2)
    ones = np.ones_like(z)
    kinds = Nc.XcorKernelTableKind

    counts = pyccl.NumberCountsTracer(
        ccl_cosmo, has_rsd=True, dndz=(z, nz), bias=(z, ones), mag_bias=(z, 0.8 * ones)
    )
    got = [c.get_kind() for c in tracer_component_tables(counts, cosmology)]
    assert got == [kinds.DENSITY, kinds.RSD, kinds.CONVERGENCE]

    # at s = 0.4 the magnification window 5 s - 2 vanishes and CCL emits zeros
    cancelled = pyccl.NumberCountsTracer(
        ccl_cosmo, has_rsd=False, dndz=(z, nz), bias=(z, ones), mag_bias=(z, 0.4 * ones)
    )
    got = [c.get_kind() for c in tracer_component_tables(cancelled, cosmology)]
    assert got == [kinds.DENSITY]

    shear = pyccl.WeakLensingTracer(ccl_cosmo, dndz=(z, nz))
    got = [c.get_kind() for c in tracer_component_tables(shear, cosmology)]
    assert got == [kinds.SHEAR]

    kernel = tracer_kernel(counts, cosmology)
    assert kernel.get_n_components() == 3
    assert [c.get_bessel_deriv() for c in kernel.get_component_list()] == [0, 2, 0]


def test_der_bessel_unsupported_is_explicit(cosmology, ccl_cosmo):
    """Shear needs derivatives of j_l, which the kernel-folding form cannot express."""
    z = np.linspace(0.0, 2.0, 256)
    nz = np.exp(-0.5 * ((z - 0.8) / 0.12) ** 2)
    wl = pyccl.WeakLensingTracer(ccl_cosmo, dndz=(z, nz))

    if -1 in wl.get_bessel_derivative():
        pytest.skip("this CCL version folds shear into der_bessel = -1")

    with pytest.raises(NotImplementedError, match="der_bessel"):
        compute_kernel(wl, cosmology, 32.0)


def test_angular_cl_spin2_weight_matches_the_lensing_kernel(cosmology, ccl_cosmo):
    """A CCL WeakLensingTracer reproduces NumCosmo's exact weak-lensing kernel.

    The tracer becomes a SHEAR table component, so the spin-2 weight
    j_l(k chi) / (k chi)^2 and the prefactor sqrt((l+2)(l+1)l(l-1)) are the
    library's, applied exactly; its Limber value 1 / nu^2 would be 6% off at
    ell = 2. What remains is the reconstruction of CCL's table and the two
    kernels' independent closure fits: measured 3.6e-6, 6.2e-6, 8.3e-6 at
    ell 2, 5, 10 at the default closure tolerance and 1.1e-7, 1.9e-7, 3.2e-7
    at 1e-6, which is what is asserted.
    """
    ells = np.array([2, 5, 10])
    z = np.linspace(0.0, 2.5, 1000)
    nz = np.exp(-0.5 * ((z - 1.0) / 0.2) ** 2)

    wl_ccl = pyccl.WeakLensingTracer(ccl_cosmo, dndz=(z, nz), n_samples=len(z))
    dndz = Ncm.Spline.new_array(
        Ncm.SplineCubicNotaknot.new(), z.tolist(), nz.tolist(), True
    )
    wl_nc = Nc.XcorKernelWeakLensing(
        dist=cosmology.dist,
        powspec=cosmology.ps_ml,
        dndz=dndz,
        nbar=1.0,
        intr_shear=0.0,
        integrator=Ncm.SBesselIntegratorLevin.new(2, int(ells[-1])),
    )
    wl_nc.set_reltol(1.0e-6)
    wl_nc.set_scaled_abstol(1.0e-6)
    wl_nc.set_l_limber(-1)
    wl_nc.prepare(cosmology.cosmo)

    lmin, lmax = int(ells[0]), int(ells[-1])
    res = Ncm.Vector.new(lmax - lmin + 1)
    xcor = Nc.Xcor.new(cosmology.dist, cosmology.ps_ml, Nc.XcorMethod.KERNEL_EXACT)
    xcor.prepare(cosmology.cosmo)
    xcor.compute(wl_nc, wl_nc, cosmology.cosmo, lmin, lmax, res)
    expected = np.array(res.dup_array())[ells - lmin]

    got = angular_cl(cosmology, wl_ccl, wl_ccl, ells, kernel_reltol=1.0e-6)

    assert_allclose(got, expected, rtol=2.0e-6)


def test_angular_cl_convergence_kind_matches_the_cmb_lensing_kernel(
    cosmology, ccl_cosmo
):
    """A CCL CMBLensingTracer reproduces NumCosmo's exact CMB lensing kernel.

    der_bessel = -1 with der_angles = 1 becomes a CONVERGENCE component: the
    same 1 / (k chi)^2 weight as shear with the prefactor l (l + 1). The
    tracer reaches the last-scattering surface; the adapter extends the
    distance object's reach as every kernel does, so the default cosmology
    serves. Its source redshift must be the recombination value the NumCosmo
    kernel uses (z_source = 1100 instead leaves a 1e-4 difference that is the
    kernel, not the method). Measured 2.7e-9, 1e-11, 2e-11 at ell 2, 5, 10.
    """
    ells = np.array([2, 5, 10])
    z_lss = cosmology.dist.decoupling_redshift(cosmology.cosmo)

    cmb_ccl = pyccl.CMBLensingTracer(ccl_cosmo, z_source=z_lss, n_samples=1000)
    noise = Ncm.Vector.new(int(ells[-1]) + 1)
    noise.set_zero()
    cmb_nc = Nc.XcorKernelCMBLensing.new(
        cosmology.dist, cosmology.ps_ml, cosmology.recomb, noise
    )
    cmb_nc.props.integrator = Ncm.SBesselIntegratorLevin.new(2, int(ells[-1]))
    cmb_nc.set_reltol(1.0e-6)
    cmb_nc.set_scaled_abstol(1.0e-6)
    cmb_nc.set_l_limber(-1)
    cmb_nc.prepare(cosmology.cosmo)

    lmin, lmax = int(ells[0]), int(ells[-1])
    res = Ncm.Vector.new(lmax - lmin + 1)
    xcor = Nc.Xcor.new(cosmology.dist, cosmology.ps_ml, Nc.XcorMethod.KERNEL_EXACT)
    xcor.prepare(cosmology.cosmo)
    xcor.compute(cmb_nc, cmb_nc, cosmology.cosmo, lmin, lmax, res)
    expected = np.array(res.dup_array())[ells - lmin]

    got = angular_cl(cosmology, cmb_ccl, cmb_ccl, ells, kernel_reltol=1.0e-6)

    assert_allclose(got, expected, rtol=1.0e-7)


def test_batching_agrees_within_the_error_budget(cosmology, tracer):
    """A multipole's transform is the same alone or inside a block, to within the bound.

    One factorisation serves a whole block and the truncation order is a max-reduction
    across it, so a multipole batched with harder neighbours is solved at a different
    (higher) order than it needs. The two answers therefore differ, and the size of the
    difference is the method's own error bar rather than zero.

    Measured here: the cancellation ratio is C = 2.7, so the budget reltol*C is 2.7e-8,
    and the observed difference is 1.4e-7 -- about 5x the estimate. That is consistent
    with a spectral-tail criterion being an estimate rather than a rigorous bound, and
    the tolerance below is set from the observation, not from wishful thinking.

    Only the dominant wavenumber is checked: at k = 0.05 the transform is eight orders
    below its peak, where a relative comparison is meaningless.
    """
    k_a = np.array([0.01])

    alone = bessel_transform_block(tracer, cosmology, 5, 5, k_a, reltol=1.0e-8)[0]
    batched = bessel_transform_block(tracer, cosmology, 2, 9, k_a, reltol=1.0e-8)[3]

    assert_allclose(batched, alone, rtol=1.0e-6)


def test_tighter_request_narrows_the_batching_difference(cosmology, tracer):
    """The difference above is the tolerance at work, not a fixed offset.

    If batching introduced a systematic error the gap would not shrink when a tighter
    accuracy is requested; it does.
    """
    k_a = np.array([0.01])

    def gap(reltol):
        alone = bessel_transform_block(tracer, cosmology, 5, 5, k_a, reltol=reltol)[0]
        batched = bessel_transform_block(tracer, cosmology, 2, 9, k_a, reltol=reltol)[3]
        return abs(batched[0] - alone[0]) / abs(alone[0])

    assert gap(1.0e-10) < gap(1.0e-6)


def test_block_size_is_clamped(cosmology, tracer):
    """An unhonourable block size is clamped rather than silently requested."""
    from numcosmo_py import Nc

    big = angular_cl(
        cosmology,
        tracer,
        tracer,
        np.array([2, 3]),
        block_size=Nc.XCOR_KERNEL_MAX_ELL_BLOCK * 100,
    )

    assert np.all(np.isfinite(big))
    assert np.all(big > 0.0)

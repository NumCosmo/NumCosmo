"""Tests for the CCL-facing non-Limber angular power spectrum."""

import numpy as np
import pytest
from numpy.testing import assert_allclose
from scipy.integrate import quad
from scipy.special import spherical_jn

import pyccl

from numcosmo_py import Nc, Ncm
from numcosmo_py.ccl.nc_ccl import create_nc_obj, CCLParams
from numcosmo_py.ccl.two_point import (
    angular_cl,
    bessel_transform,
    bessel_transform_block,
    compute_kernel,
)

pytestmark = [pytest.mark.ccl, pytest.mark.xcor]

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
    got = angular_cl(cosmology, tracer, tracer, np.array([ell]), reltol=1.0e-8)[0]
    ref = pyccl.angular_cl(ccl_cosmo, tracer, tracer, np.array([ell]), l_limber=100000)[
        0
    ]

    assert_allclose(got, ref, rtol=3.0e-3)


def test_outer_integral_is_converged(cosmology, tracer):
    """The default samples-per-oscillation is enough: 4x more changes nothing."""
    a = angular_cl(cosmology, tracer, tracer, np.array([2]), reltol=1.0e-8, n_osc=8)[0]
    b = angular_cl(cosmology, tracer, tracer, np.array([2]), reltol=1.0e-8, n_osc=32)[0]

    assert_allclose(a, b, rtol=1.0e-5)


def test_der_bessel_unsupported_is_explicit(cosmology, ccl_cosmo):
    """Shear needs derivatives of j_l, which the kernel-folding form cannot express."""
    z = np.linspace(0.0, 2.0, 256)
    nz = np.exp(-0.5 * ((z - 0.8) / 0.12) ** 2)
    wl = pyccl.WeakLensingTracer(ccl_cosmo, dndz=(z, nz))

    if -1 in wl.get_bessel_derivative():
        pytest.skip("this CCL version folds shear into der_bessel = -1")

    with pytest.raises(NotImplementedError, match="der_bessel"):
        compute_kernel(wl, cosmology, 32.0)


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
        reltol=1.0e-6,
        block_size=Nc.XCOR_KERNEL_MAX_ELL_BLOCK * 100,
    )

    assert np.all(np.isfinite(big))
    assert np.all(big > 0.0)

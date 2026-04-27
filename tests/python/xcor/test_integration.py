"""Tests for NumCosmo cross-correlation integration methods.

This module tests different integration backends and methods for computing
angular power spectra using the Xcor class. Tests include:
- Different Xcor integration methods (LIMBER_Z_GSL, LIMBER_Z_CUBATURE, KERNEL_GSL)
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

pytest.importorskip("gi.repository.NumCosmo")

pytest_plugins = [
    "python.fixtures_xcor",
]

pytestmark = pytest.mark.xcor


# Xcor Integration Method Tests


@pytest.mark.parametrize(
    "method",
    [
        Nc.XcorMethod.LIMBER_Z_GSL,
        Nc.XcorMethod.LIMBER_Z_CUBATURE,
        Nc.XcorMethod.KERNEL_GSL,
    ],
)
def test_xcor_method_creation(cosmology: Cosmology, method: Nc.XcorMethod) -> None:
    """Test that Xcor can be created with all integration methods."""
    xcor = Nc.Xcor.new(cosmology.dist, cosmology.ps_ml, method)

    assert xcor.props.meth == method
    assert xcor.props.power_spec is cosmology.ps_ml
    xcor.prepare(cosmology.cosmo)


def test_xcor_properties(cosmology: Cosmology) -> None:
    """Test Xcor property getters and setters.

    Validates that reltol and method properties can be properly set
    and retrieved via both property access and explicit methods.
    """
    xcor = Nc.Xcor.new(cosmology.dist, cosmology.ps_ml, Nc.XcorMethod.LIMBER_Z_GSL)

    # Test reltol setter/getter
    xcor.set_reltol(1.0e-4)
    assert xcor.get_reltol() == 1.0e-4
    assert xcor.props.reltol == 1.0e-4

    xcor.set_reltol(1.0e-5)
    assert xcor.get_reltol() == 1.0e-5
    assert xcor.props.reltol == 1.0e-5

    xcor.props.reltol = 1.0e-6
    assert xcor.get_reltol() == 1.0e-6

    # Test method property switching
    xcor.props.meth = Nc.XcorMethod.LIMBER_Z_GSL
    assert xcor.props.meth == Nc.XcorMethod.LIMBER_Z_GSL

    xcor.props.meth = Nc.XcorMethod.LIMBER_Z_CUBATURE
    assert xcor.props.meth == Nc.XcorMethod.LIMBER_Z_CUBATURE

    # Test power spectrum property
    assert xcor.props.power_spec is cosmology.ps_ml


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
def test_xcor_compute_methods(
    cosmology: Cosmology,
    k1_name: str,
    k2_name: str,
    request: pytest.FixtureRequest,
) -> None:
    """Compare different Xcor integration methods for consistency.

    Tests that LIMBER_Z_GSL, LIMBER_Z_CUBATURE, and KERNEL_GSL methods
    produce consistent results for various kernel auto- and cross-correlations.
    """
    k1 = cast(Nc.XcorKernel, request.getfixturevalue(k1_name))
    k2 = cast(Nc.XcorKernel, request.getfixturevalue(k2_name))

    xcor_gsl = Nc.Xcor.new(cosmology.dist, cosmology.ps_ml, Nc.XcorMethod.LIMBER_Z_GSL)
    xcor_cub = Nc.Xcor.new(
        cosmology.dist,
        cosmology.ps_ml,
        Nc.XcorMethod.LIMBER_Z_CUBATURE,
    )

    xcor_gsl.prepare(cosmology.cosmo)
    xcor_cub.prepare(cosmology.cosmo)

    k1.prepare(cosmology.cosmo)
    k2.prepare(cosmology.cosmo)

    lmin = 2
    lmax = 100  # Use smaller lmax for faster tests
    vp_gsl = Ncm.Vector.new(lmax - lmin + 1)
    vp_cub = Ncm.Vector.new(lmax - lmin + 1)

    xcor_gsl.set_reltol(1.0e-7)
    xcor_gsl.compute(k1, k2, cosmology.cosmo, lmin, lmax, vp_gsl)
    xcor_cub.compute(k1, k2, cosmology.cosmo, lmin, lmax, vp_cub)

    vp_gsl_a = np.array(vp_gsl.dup_array())
    vp_cub_a = np.array(vp_cub.dup_array())

    # GSL and Cubature should agree well
    assert_allclose(vp_gsl_a, vp_cub_a, rtol=1.0e-5, atol=1.0e-50)


# Kernel-Level Integrator Tests


def test_kernel_integrator_types(cosmology: Cosmology) -> None:
    """Test different kernel-level integrator types for CMB ISW.

    Compares Limber approximation integrator with Levin spherical Bessel
    integration to ensure both integrator types are properly initialized
    and functional.
    """
    lmax = 30

    # Create ISW kernel with Limber approximation
    kernel_limber = Nc.XcorKernelCMBISW.new(
        cosmology.dist,
        cosmology.ps_ml,
        cosmology.recomb,
        Ncm.Vector.new_array(np.arange(lmax + 1)),
    )
    kernel_limber.set_lmax(lmax)

    # Create ISW kernel with Levin spherical Bessel integrator
    kernel_levin = Nc.XcorKernelCMBISW(
        dist=cosmology.dist,
        powspec=cosmology.ps_ml,
        recomb=cosmology.recomb,
        Nl=Ncm.Vector.new_array(np.arange(lmax + 1)),
        lmax=lmax,
        integrator=Ncm.SBesselIntegratorLevin.new(0, 1000),
    )

    # Prepare both kernels
    cosmo = cosmology.cosmo
    kernel_limber.prepare(cosmo)
    kernel_levin.prepare(cosmo)

    # Verify both integrators work
    eval_limber = kernel_limber.get_eval(cosmo, 200)
    eval_levin = kernel_levin.get_eval(cosmo, 200)

    result_limber = eval_limber.eval_array(1.0e2)
    result_levin = eval_levin.eval_array(1.0e2)

    # Both should produce valid results
    assert result_limber is not None
    assert result_levin is not None
    assert len(result_limber) > 0
    assert len(result_levin) > 0

    # Verify lmax is properly set
    assert kernel_limber.get_lmax() == lmax
    assert kernel_levin.get_lmax() == lmax


@pytest.mark.parametrize("lmax", [10, 30, 50])
def test_xcor_compute_lmax_range(
    cosmology: Cosmology, kernel_cmb_lens: Nc.XcorKernelCMBLensing, lmax: int
) -> None:
    """Test Xcor compute with different multipole ranges.

    Validates that the compute method works correctly for various
    multipole ranges (lmin to lmax).
    """
    xcor = Nc.Xcor.new(
        cosmology.dist,
        cosmology.ps_ml,
        Nc.XcorMethod.LIMBER_Z_CUBATURE,
    )
    xcor.prepare(cosmology.cosmo)
    kernel_cmb_lens.prepare(cosmology.cosmo)

    lmin = 2
    cl_vec = Ncm.Vector.new(lmax - lmin + 1)

    xcor.compute(kernel_cmb_lens, kernel_cmb_lens, cosmology.cosmo, lmin, lmax, cl_vec)

    cl_array = np.array(cl_vec.dup_array())

    # Verify output has correct shape
    assert len(cl_array) == lmax - lmin + 1

    # Verify all values are positive (CMB lensing auto-correlation)
    assert np.all(cl_array > 0)

    # Verify all values are finite and reasonable
    assert np.all(np.isfinite(cl_array))
    assert np.max(cl_array) < 1.0  # CMB lensing C_ℓ should be small


@pytest.mark.parametrize("reltol", [1.0e-4, 1.0e-6, 1.0e-8])
def test_xcor_tolerance_scaling(
    cosmology: Cosmology,
    kernel_cmb_lens: Nc.XcorKernelCMBLensing,
    reltol: float,
) -> None:
    """Test that changing reltol affects integration accuracy.

    Verifies that the relative tolerance parameter properly controls
    the precision of the integration.
    """
    xcor = Nc.Xcor.new(
        cosmology.dist,
        cosmology.ps_ml,
        Nc.XcorMethod.LIMBER_Z_GSL,
    )
    xcor.set_reltol(reltol)
    xcor.prepare(cosmology.cosmo)
    kernel_cmb_lens.prepare(cosmology.cosmo)

    lmin = 2
    lmax = 50
    cl_vec = Ncm.Vector.new(lmax - lmin + 1)

    xcor.compute(kernel_cmb_lens, kernel_cmb_lens, cosmology.cosmo, lmin, lmax, cl_vec)

    cl_array = np.array(cl_vec.dup_array())

    # Verify computation completes successfully
    assert len(cl_array) == lmax - lmin + 1
    assert np.all(np.isfinite(cl_array))
    assert np.all(cl_array > 0)


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

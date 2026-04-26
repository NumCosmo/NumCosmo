#!/usr/bin/env python
#
# test_py_xcor_kernel_observables.py
#
# Thu Apr 17 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_py_xcor_kernel_observables.py
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

"""Unit tests for XcorKernel observables.

Tests that kernels have the correct number of observables and parameters.
Modernized from test_py_xcor.py observable tests.
"""

import pytest
from numpy.testing import assert_allclose
import numpy as np

from numcosmo_py import Ncm, Nc
from numcosmo_py.cosmology import Cosmology

pytestmark = pytest.mark.xcor

Ncm.cfg_init()
pytest_plugins = ["python.fixtures_xcor"]


# Kernel observable specifications: (kernel_type, expected_obs_len, expected_params_len)
KERNEL_OBS_SPECS = {
    "kernel_cmb_lens": (Nc.XcorKernelCMBLensing, 1, 0),
    "kernel_cmb_isw": (Nc.XcorKernelCMBISW, 1, 0),
    "kernel_tsz": (Nc.XcorKerneltSZ, 1, 0),
    "kernel_gal": (Nc.XcorKernelGal, 2, 1),
    "kernel_wl": (Nc.XcorKernelWeakLensing, 2, 1),
}


@pytest.mark.parametrize(
    "kernel_fixture_name,expected_type,expected_obs,expected_params",
    [
        ("kernel_cmb_lens", Nc.XcorKernelCMBLensing, 1, 0),
        ("kernel_cmb_isw", Nc.XcorKernelCMBISW, 1, 0),
        ("kernel_tsz", Nc.XcorKerneltSZ, 1, 0),
    ],
)
def test_kernel_basic_observables(
    kernel_fixture_name: str,
    expected_type: type,
    expected_obs: int,
    expected_params: int,
    request: pytest.FixtureRequest,
) -> None:
    """Test that kernels have the correct number of observables and parameters.

    This test covers basic kernels (CMB lensing, CMB ISW, tSZ) that have simple
    observable structures without additional complexity.
    """
    kernel: Nc.XcorKernel = request.getfixturevalue(kernel_fixture_name)

    assert kernel is not None
    assert isinstance(kernel, expected_type)
    assert kernel.obs_len() == expected_obs
    assert kernel.obs_params_len() == expected_params


def test_galaxy_kernel_observables(kernel_gal: Nc.XcorKernelGal) -> None:
    """Test that galaxy kernel has the correct observables and specific features.

    Galaxy kernels have:
    - 2 observables (number counts and magnification)
    - 1 observable parameter (bias)
    - Fast update capability
    - Bias parameters that can be set and retrieved
    """
    # Basic observable counts
    assert kernel_gal is not None
    assert isinstance(kernel_gal, Nc.XcorKernelGal)
    assert kernel_gal.obs_len() == 2
    assert kernel_gal.obs_params_len() == 1

    # Test fast update mode
    kernel_gal.set_fast_update(True)
    assert kernel_gal.get_fast_update()

    kernel_gal.set_fast_update(False)
    assert not kernel_gal.get_fast_update()

    # Test bias parameter setting and getting (if available)
    if kernel_gal.vparam_len(Nc.XcorKernelGalVParams.BIAS) == 1:
        # Set bias using the vector parameter interface
        kernel_gal.orig_vparam_set(Nc.XcorKernelGalVParams.BIAS, 0, 3.21)

        # Set old bias parameters
        kernel_gal.set_bias_old(1.2345, 0.9876)

        # Retrieve and verify all bias parameters
        bias0, bias_old0, noise_bias_old0 = kernel_gal.get_bias()

        assert_allclose(bias0, 3.21, atol=0.0)
        assert_allclose(bias_old0, 1.2345, atol=0.0)
        assert_allclose(noise_bias_old0, 0.9876, atol=0.0)


def test_weak_lensing_kernel_observables(
    kernel_wl: Nc.XcorKernelWeakLensing,
) -> None:
    """Test that weak lensing kernel has the correct observables.

    Weak lensing kernels have:
    - 2 observables (cosmic shear components)
    - 1 observable parameter (multiplicative bias)
    """
    assert kernel_wl is not None
    assert isinstance(kernel_wl, Nc.XcorKernelWeakLensing)
    assert kernel_wl.obs_len() == 2
    assert kernel_wl.obs_params_len() == 1


def test_galaxy_kernel_extrapolation(
    cosmology: Cosmology, kernel_gal: Nc.XcorKernelGal
) -> None:
    """Test that galaxy kernel can extrapolate beyond the dndz range.

    This test verifies that:
    1. Kernel evaluates correctly at the edge of the dndz range
    2. Kernel can safely extrapolate beyond the dndz range
    3. All extrapolated values are finite

    This is important for numerical stability when kernels are evaluated
    at redshifts slightly beyond the defined photo-z distribution range.
    """
    cosmo = cosmology.cosmo
    dist = cosmology.dist

    # Prepare the kernel
    kernel_gal.prepare(cosmo)

    # Get the redshift array from the dndz distribution
    z_a = np.array(kernel_gal.props.dndz.peek_xv().dup_array())

    # Test evaluation at the edge of the dndz range
    edge_value = kernel_gal.eval_limber_z_full(cosmo, z_a[-1], dist, 77)
    assert np.isfinite(edge_value), "Kernel should evaluate at edge of dndz range"

    # Test extrapolation beyond the dndz range
    extra_z_a = np.linspace(z_a[-1], 2.0 * z_a[-1], 10)
    k_a = np.array(
        [kernel_gal.eval_limber_z_full(cosmo, z, dist, 77) for z in extra_z_a]
    )

    # All extrapolated values should be finite
    assert np.isfinite(
        k_a
    ).all(), "Kernel should handle extrapolation beyond dndz range"


def test_all_kernels_have_observables(kernel: Nc.XcorKernel) -> None:
    """Test that all kernel types have at least one observable.

    This is a sanity check that runs on all kernel types (parametrized via
    the kernel fixture) to ensure they properly implement the observable
    interface.
    """
    assert kernel is not None
    assert isinstance(kernel, Nc.XcorKernel)

    # All kernels must have at least one observable
    obs_len = kernel.obs_len()
    assert (
        obs_len >= 1
    ), f"Kernel {type(kernel).__name__} must have at least 1 observable"

    # Observable parameter count should be non-negative
    obs_params_len = kernel.obs_params_len()
    assert (
        obs_params_len >= 0
    ), f"Kernel {type(kernel).__name__} observable parameter count must be >= 0"


def test_kernel_type_consistency(kernel_case: tuple[str, Nc.XcorKernel]) -> None:
    """Test that kernel types match expected types based on their names.

    This ensures the fixture system is creating the correct kernel types and
    that the naming convention is consistent.
    """
    kernel_name, kernel = kernel_case

    # Map kernel name patterns to expected types
    type_map = {
        "cmb_lens": Nc.XcorKernelCMBLensing,
        "cmb_isw": Nc.XcorKernelCMBISW,
        "tsz": Nc.XcorKerneltSZ,
        "gal": Nc.XcorKernelGal,
        "wl": Nc.XcorKernelWeakLensing,
    }

    # Check if kernel name contains any of the type identifiers
    for key, expected_type in type_map.items():
        if key in kernel_name:
            assert isinstance(
                kernel, expected_type
            ), f"Kernel {kernel_name} should be of type {expected_type.__name__}"
            break

#!/usr/bin/env python
#
# test_py_xcor_kernel_noise.py
#
# Thu Apr 17 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_py_xcor_kernel_noise.py
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

"""Unit tests for XcorKernel noise behavior.

Tests that kernels correctly add noise to power spectra using the add_noise method.
Modernized from test_py_xcor.py noise tests.
"""

from typing import Callable
import pytest
from numpy.testing import assert_allclose
import numpy as np

from numcosmo_py import Ncm, Nc

pytestmark = pytest.mark.xcor

Ncm.cfg_init()
pytest_plugins = ["python.fixtures_xcor"]


def test_cmb_lens_noise(kernel_cmb_lens: Nc.XcorKernelCMBLensing) -> None:
    """Test that CMB lensing kernel adds noise correctly.

    CMB lensing noise is added starting at index 'ell_idx' with an increment
    pattern that depends on the position. For a vector of all 1.0s starting
    at index 5, the result should be np.arange(6, 16).
    """
    vp1 = Ncm.Vector.new(10)
    vp2 = Ncm.Vector.new(10)
    vp1.set_all(1.0)
    kernel_cmb_lens.add_noise(vp1, vp2, 5)
    assert_allclose(vp2.dup_array(), np.arange(6, 16), atol=0.0)


def test_cmb_isw_noise(kernel_cmb_isw: Nc.XcorKernelCMBISW) -> None:
    """Test that CMB ISW kernel adds noise correctly.

    CMB ISW noise follows the same pattern as CMB lensing - it is added
    starting at index 'ell_idx' with an increment pattern based on position.
    """
    vp1 = Ncm.Vector.new(10)
    vp2 = Ncm.Vector.new(10)
    vp1.set_all(1.0)
    kernel_cmb_isw.add_noise(vp1, vp2, 5)
    assert_allclose(vp2.dup_array(), np.arange(6, 16), atol=0.0)


def test_tsz_noise(kernel_tsz: Nc.XcorKerneltSZ) -> None:
    """Test that tSZ kernel does not add noise.

    The tSZ kernel does not modify the power spectrum when add_noise is called.
    The output should remain identical to the input.
    """
    vp1 = Ncm.Vector.new(10)
    vp2 = Ncm.Vector.new(10)
    vp1.set_all(1.0)
    kernel_tsz.add_noise(vp1, vp2, 5)
    assert_allclose(vp2.dup_array(), np.ones(10) * 1.0, atol=0.0)


def test_gal_noise(kernel_gal: Nc.XcorKernelGal) -> None:
    """Test that galaxy kernel adds constant noise.

    Galaxy kernels add a constant noise value of 2.234 to the power spectrum.
    This accounts for shot noise from the discrete galaxy distribution.
    """
    vp1 = Ncm.Vector.new(10)
    vp2 = Ncm.Vector.new(10)
    vp1.set_all(1.0)
    kernel_gal.add_noise(vp1, vp2, 5)
    assert_allclose(vp2.dup_array(), np.ones(10) * 2.234, atol=0.0)


def test_weak_lensing_noise(kernel_wl: Nc.XcorKernelWeakLensing) -> None:
    """Test that weak lensing kernel adds shape noise correctly.

    Weak lensing kernels add noise based on intrinsic ellipticity distribution
    and galaxy density. The noise formula is (intr_shear)^2 / nbar + input,
    where intr_shear=7.0 and nbar=3.0 from the fixture.
    """
    vp1 = Ncm.Vector.new(10)
    vp2 = Ncm.Vector.new(10)
    vp1.set_all(1.0)
    kernel_wl.add_noise(vp1, vp2, 5)
    expected = np.ones(10) * (7.0) ** 2 / 3.0 + 1.0
    assert_allclose(vp2.dup_array(), expected, atol=0.0)


# Parametrized tests for systematic verification


@pytest.mark.parametrize(
    "kernel_fixture_name,expected_result_factory",
    [
        (
            "kernel_cmb_lens",
            lambda: np.arange(6, 16),
        ),
        (
            "kernel_cmb_isw",
            lambda: np.arange(6, 16),
        ),
        (
            "kernel_tsz",
            lambda: np.ones(10) * 1.0,
        ),
        # Note: kernel_gal and kernel_wl are already tested via their own
        # parametrized fixtures in test_gal_noise and test_weak_lensing_noise
    ],
)
def test_kernel_noise_parametrized(
    kernel_fixture_name: str,
    expected_result_factory: Callable[[], np.ndarray],
    request: pytest.FixtureRequest,
) -> None:
    """Parametrized test for kernel noise behavior across CMB and tSZ kernels.

    This test systematically verifies that each basic kernel type (non-binned)
    adds noise correctly using a consistent test pattern. Galaxy and weak lensing
    kernels are tested separately via their parametrized fixtures.

    Args:
        kernel_fixture_name: Name of the fixture providing the kernel
        expected_result_factory: Function that returns expected result array
        request: pytest fixture request object
    """
    kernel: Nc.XcorKernel = request.getfixturevalue(kernel_fixture_name)

    vp1 = Ncm.Vector.new(10)
    vp2 = Ncm.Vector.new(10)
    vp1.set_all(1.0)

    # Add noise starting at index 5
    kernel.add_noise(vp1, vp2, 5)

    expected = expected_result_factory()
    assert_allclose(vp2.dup_array(), expected, atol=0.0)


def test_noise_with_different_ell_indices(
    kernel_cmb_lens: Nc.XcorKernelCMBLensing,
) -> None:
    """Test that noise addition respects different starting indices.

    The ell_idx parameter determines where noise addition starts in the vector.
    This test verifies that the noise pattern correctly shifts based on this index.
    """
    vp1 = Ncm.Vector.new(10)
    vp1.set_all(1.0)

    # Test different starting indices
    for start_idx in [0, 3, 5, 7]:
        vp2 = Ncm.Vector.new(10)
        kernel_cmb_lens.add_noise(vp1, vp2, start_idx)

        result = vp2.dup_array()
        expected = np.arange(start_idx + 1, start_idx + 11)

        assert_allclose(result, expected, atol=0.0)


def test_noise_with_different_input_values(kernel_wl: Nc.XcorKernelWeakLensing) -> None:
    """Test that noise is properly added to different input power spectra.

    This verifies that the noise addition is correctly implemented as an
    additive operation on top of the existing power spectrum values.
    """
    vp1 = Ncm.Vector.new(10)
    vp2 = Ncm.Vector.new(10)

    # Test with different input values
    test_values = [0.5, 1.0, 2.0, 5.0]

    for input_val in test_values:
        vp1.set_all(input_val)
        kernel_wl.add_noise(vp1, vp2, 5)

        # Weak lensing noise formula: (intr_shear)^2 / nbar + input
        expected = np.ones(10) * ((7.0) ** 2 / 3.0 + input_val)
        assert_allclose(vp2.dup_array(), expected, atol=0.0)


def test_noise_vector_size_consistency(kernel_gal: Nc.XcorKernelGal) -> None:
    """Test that noise addition works correctly with different vector sizes.

    This ensures that the add_noise method handles vectors of various lengths
    correctly, which is important for power spectra with different ell ranges.
    """
    test_sizes = [5, 10, 20, 50]

    for size in test_sizes:
        vp1 = Ncm.Vector.new(size)
        vp2 = Ncm.Vector.new(size)
        vp1.set_all(1.0)

        # Add noise starting at index 0
        kernel_gal.add_noise(vp1, vp2, 0)

        expected = np.ones(size) * 2.234
        assert_allclose(vp2.dup_array(), expected, atol=0.0)


def test_noise_does_not_modify_input(kernel_cmb_lens: Nc.XcorKernelCMBLensing) -> None:
    """Test that add_noise does not modify the input vector.

    The add_noise method should only write to the output vector (vp2)
    and leave the input vector (vp1) unchanged.
    """
    vp1 = Ncm.Vector.new(10)
    vp2 = Ncm.Vector.new(10)
    vp1.set_all(1.0)

    original_input = vp1.dup_array().copy()

    kernel_cmb_lens.add_noise(vp1, vp2, 5)

    # Verify input vector is unchanged
    assert_allclose(vp1.dup_array(), original_input, atol=0.0)

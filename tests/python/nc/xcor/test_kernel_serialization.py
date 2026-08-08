#!/usr/bin/env python
#
# test_py_xcor_kernel_serialization.py
#
# Thu Apr 17 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_py_xcor_kernel_serialization.py
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

"""Unit tests for XcorKernel serialization.

Tests that kernels can be properly serialized and deserialized, producing
functionally equivalent duplicates. Modernized from test_py_xcor.py
serialization tests.
"""

import pytest
from numpy.testing import assert_allclose
import numpy as np

from numcosmo_py import Ncm, Nc
from numcosmo_py.cosmology import Cosmology
from numcosmo_py.helper import duplicate_via_serialization

pytestmark = pytest.mark.xcor

Ncm.cfg_init()
pytest_plugins = ["python.fixtures_xcor"]


def test_kernel_serialization_basic(
    kernel: Nc.XcorKernel, cosmology: Cosmology
) -> None:
    """Test that all kernels can be serialized and deserialized.

    This test verifies that:
    1. Serialization produces a duplicate object (not the same reference)
    2. The duplicate is of the same type as the original
    3. Both objects can be prepared without errors

    This is a basic sanity check that runs on all kernel types.
    """
    # Duplicate the kernel via serialization
    kernel_dup = duplicate_via_serialization(kernel)

    # Verify duplicate was created
    assert kernel_dup is not None
    assert kernel_dup is not kernel, "Duplicate should be a different object"
    assert isinstance(kernel_dup, type(kernel)), "Duplicate should have same type"

    # Both should be preparable
    cosmo = cosmology.cosmo
    kernel.prepare(cosmo)
    kernel_dup.prepare(cosmo)


def test_kernel_serialization_evaluation(
    kernel: Nc.XcorKernel, cosmology: Cosmology
) -> None:
    """Test that serialized kernels produce identical evaluation results.

    This test verifies that after serialization and preparation, both the
    original and duplicate kernel produce identical results when evaluated
    at the same point in parameter space.
    """
    # Duplicate the kernel via serialization
    kernel_dup = duplicate_via_serialization(kernel)

    cosmo = cosmology.cosmo
    dist = cosmology.dist

    # Prepare both kernels
    kernel.prepare(cosmo)
    kernel_dup.prepare(cosmo)

    # Get valid redshift range for this kernel
    zmin, zmax, _ = kernel.get_z_range()

    # Test evaluation at a point within the valid range
    # Use z=0 if in range, otherwise use zmin
    test_z = 0.0 if zmin <= 0.0 <= zmax else zmin
    ell = 2

    result_orig = kernel.eval_limber_z_full(cosmo, test_z, dist, ell)
    result_dup = kernel_dup.eval_limber_z_full(cosmo, test_z, dist, ell)

    assert_allclose(
        result_orig,
        result_dup,
        rtol=1e-15,
        err_msg="Original and duplicated kernel should produce identical results",
    )


def test_kernel_serialization_prepare_idempotency(
    kernel: Nc.XcorKernel, cosmology: Cosmology
) -> None:
    """Test that prepare() can be called multiple times without changing results.

    This test verifies that:
    1. Calling prepare() multiple times is safe
    2. Results remain consistent after re-preparation
    3. Serialized kernels also maintain this property

    This is important for robustness in workflows where prepare() might
    be called multiple times (e.g., in optimization loops).
    """
    # Duplicate the kernel via serialization
    kernel_dup = duplicate_via_serialization(kernel)

    cosmo = cosmology.cosmo
    dist = cosmology.dist

    # Prepare both kernels
    kernel.prepare(cosmo)
    kernel_dup.prepare(cosmo)

    # Get valid redshift range
    zmin, zmax, _ = kernel.get_z_range()
    test_z = 0.0 if zmin <= 0.0 <= zmax else zmin
    ell = 2

    # Get results after first prepare
    result_orig_1 = kernel.eval_limber_z_full(cosmo, test_z, dist, ell)
    result_dup_1 = kernel_dup.eval_limber_z_full(cosmo, test_z, dist, ell)

    # Prepare for a second time
    kernel.prepare(cosmo)
    kernel_dup.prepare(cosmo)

    # Get results after second prepare
    result_orig_2 = kernel.eval_limber_z_full(cosmo, test_z, dist, ell)
    result_dup_2 = kernel_dup.eval_limber_z_full(cosmo, test_z, dist, ell)

    # Results should be identical after re-preparation
    assert_allclose(
        result_orig_1,
        result_orig_2,
        rtol=1e-15,
        err_msg="Original kernel results should not change after re-prepare",
    )
    assert_allclose(
        result_dup_1,
        result_dup_2,
        rtol=1e-15,
        err_msg="Duplicated kernel results should not change after re-prepare",
    )
    assert_allclose(
        result_orig_2,
        result_dup_2,
        rtol=1e-15,
        err_msg="Original and duplicate should still agree after re-prepare",
    )


def test_kernel_serialization_multiple_redshifts(
    kernel: Nc.XcorKernel, cosmology: Cosmology
) -> None:
    """Test serialization consistency across multiple redshifts.

    This test evaluates kernels at multiple redshift points to ensure
    serialization preserves behavior across the full valid range, not
    just at a single point.
    """
    # Duplicate the kernel via serialization
    kernel_dup = duplicate_via_serialization(kernel)

    cosmo = cosmology.cosmo
    dist = cosmology.dist

    # Prepare both kernels
    kernel.prepare(cosmo)
    kernel_dup.prepare(cosmo)

    # Get valid redshift range
    zmin, zmax, _ = kernel.get_z_range()

    # Test at multiple redshifts across the valid range
    # Use linspace for CMB kernels (which can handle z=0),
    # and avoid z=0 for others
    if zmin <= 0.0 <= zmax:
        z_test = np.linspace(max(zmin, 0.0), min(zmax, 3.0), 5)
    else:
        z_test = np.linspace(zmin, min(zmax, 3.0), 5)

    ell = 10

    for z in z_test:
        result_orig = kernel.eval_limber_z_full(cosmo, z, dist, ell)
        result_dup = kernel_dup.eval_limber_z_full(cosmo, z, dist, ell)

        assert_allclose(
            result_orig,
            result_dup,
            rtol=1e-11,
            atol=1e-20,
            err_msg=f"Results should match at z={z}",
        )


def test_kernel_serialization_multiple_ells(
    kernel: Nc.XcorKernel, cosmology: Cosmology
) -> None:
    """Test serialization consistency across multiple multipoles.

    This test evaluates kernels at multiple ell values to ensure
    serialization preserves behavior for different angular scales.
    """
    # Duplicate the kernel via serialization
    kernel_dup = duplicate_via_serialization(kernel)

    cosmo = cosmology.cosmo
    dist = cosmology.dist

    # Prepare both kernels
    kernel.prepare(cosmo)
    kernel_dup.prepare(cosmo)

    # Get valid redshift range
    zmin, zmax, _ = kernel.get_z_range()
    test_z = 0.5 if zmin <= 0.5 <= zmax else (zmin + zmax) / 2.0

    # Test at multiple ells
    ell_test = [2, 10, 50, 100, 500]

    for ell in ell_test:
        result_orig = kernel.eval_limber_z_full(cosmo, test_z, dist, ell)
        result_dup = kernel_dup.eval_limber_z_full(cosmo, test_z, dist, ell)

        assert_allclose(
            result_orig,
            result_dup,
            rtol=1e-11,
            atol=1e-20,
            err_msg=f"Results should match at ell={ell}",
        )


def test_kernel_serialization_preserves_properties(
    kernel: Nc.XcorKernel, cosmology: Cosmology
) -> None:
    """Test that serialization preserves kernel properties.

    This test verifies that getter methods return the same values for
    both the original and serialized kernel, ensuring that internal
    state is properly preserved through serialization.
    """
    # Duplicate the kernel via serialization (type-safe!)
    kernel_dup = duplicate_via_serialization(kernel)

    # Prepare both (some properties may only be accessible after prepare)
    cosmo = cosmology.cosmo
    kernel.prepare(cosmo)
    kernel_dup.prepare(cosmo)

    # Test observable counts (should be the same)
    assert kernel.obs_len() == kernel_dup.obs_len()
    assert kernel.obs_params_len() == kernel_dup.obs_params_len()

    # Test adaptive refinement properties
    assert kernel.get_reltol() == kernel_dup.get_reltol()
    assert kernel.get_scaled_abstol() == kernel_dup.get_scaled_abstol()
    assert kernel.get_max_border_expansions() == kernel_dup.get_max_border_expansions()
    assert kernel.get_max_iter() == kernel_dup.get_max_iter()
    assert kernel.get_expansion_factor() == kernel_dup.get_expansion_factor()

    # Test redshift range
    zmin_orig, zmax_orig, _ = kernel.get_z_range()
    zmin_dup, zmax_dup, _ = kernel_dup.get_z_range()
    assert_allclose(zmin_orig, zmin_dup)
    assert_allclose(zmax_orig, zmax_dup)


def test_serializer_options(kernel: Nc.XcorKernel, cosmology: Cosmology) -> None:
    """Test different serializer options produce valid results.

    This test verifies that different serialization options (CLEAN_DUP,
    NONE) both work correctly for kernel objects.
    """
    cosmo = cosmology.cosmo
    dist = cosmology.dist

    # Test with CLEAN_DUP option
    ser_clean = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
    kernel_dup_clean = duplicate_via_serialization(kernel, ser_clean)

    # Test with no special options
    ser_none = Ncm.Serialize.new(Ncm.SerializeOpt.NONE)
    kernel_dup_none = duplicate_via_serialization(kernel, ser_none)

    # Prepare all kernels
    kernel.prepare(cosmo)
    kernel_dup_clean.prepare(cosmo)
    kernel_dup_none.prepare(cosmo)

    # Get valid redshift
    zmin, zmax, _ = kernel.get_z_range()
    test_z = 0.0 if zmin <= 0.0 <= zmax else zmin
    ell = 2

    # All should produce the same result
    result_orig = kernel.eval_limber_z_full(cosmo, test_z, dist, ell)
    result_clean = kernel_dup_clean.eval_limber_z_full(cosmo, test_z, dist, ell)
    result_none = kernel_dup_none.eval_limber_z_full(cosmo, test_z, dist, ell)

    assert_allclose(result_orig, result_clean, rtol=1e-15)
    assert_allclose(result_orig, result_none, rtol=1e-15)


_REGISTRY_CHILD = """
import sys
from numcosmo_py import Ncm
from gi.repository import GObject

Ncm.cfg_init()

missing = []
for name in sys.argv[1:]:
    try:
        GObject.type_from_name(name)
    except RuntimeError:
        missing.append(name)
print(",".join(missing))
"""


def test_every_kernel_type_is_registered_in_a_fresh_process() -> None:
    """Concrete kernel types must resolve by name after cfg_init() alone.

    ncm_serialize_from_name_params() resolves a serialized name with
    g_type_from_name(), which only succeeds once the GType has been realized in
    the process. ncm_cfg_register_objects() realizes it at startup; otherwise
    merely touching the class from Python realizes it too, which is why
    dup_obj-based tests -- and even a to_string/from_string round-trip on a live
    instance -- pass for a kernel that was never registered. Reading a YAML
    experiment in a fresh process does neither, so it aborts with "object `X' is
    not registered". This ran green while NcXcorKernelCMBISW was missing from
    ncm_cfg_register_objects(), so it must not construct anything.
    """
    import subprocess
    import sys

    # Names are collected here, where realizing the classes is harmless; the
    # child only ever looks them up.
    names = sorted(
        "Nc" + n
        for n in dir(Nc)
        if n.startswith("XcorKernel")
        and not n.endswith(("SParams", "VParams", "Class", "Impl", "Integrand"))
        and n not in ("XcorKernel", "XcorKernelComponent")
    )
    assert "NcXcorKernelCMBISW" in names, "ISW kernel dropped from the sweep"

    proc = subprocess.run(
        [sys.executable, "-c", _REGISTRY_CHILD, *names],
        capture_output=True,
        text=True,
        check=False,
        timeout=600,
    )

    assert proc.returncode == 0, proc.stderr
    missing = [n for n in proc.stdout.strip().split(",") if n]
    assert not missing, (
        f"not passed to ncm_cfg_register_obj(): {missing} -- these cannot be "
        f"built from a serialized name in a fresh process"
    )


def test_kernel_reconstructible_by_type_name(
    kernel: Nc.XcorKernel, cosmology: Cosmology
) -> None:
    """Every kernel type must be reconstructible from its serialized name.

    Ncm.Serialize.dup_obj() resolves the GType from the live instance, so it
    works whether or not the type was passed to ncm_cfg_register_obj(). Going
    through a *string* does not: ncm_serialize_from_name_params() looks the name
    up in that registry and aborts if it is absent. That is the path YAML
    experiment files, the CLI and ncm_obj_array take, so a kernel missing from
    ncm_cfg_register_objects() is unusable there while every dup_obj-based test
    still passes -- which is exactly how NcXcorKernelCMBISW stayed broken.
    """
    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)

    # The GType name keeps the Nc namespace that introspection strips, and it
    # is the name the registry is keyed on.
    gtype_name = kernel.__gtype__.name
    serialized = ser.to_string(kernel, False)
    assert serialized.startswith(gtype_name)

    rebuilt = ser.from_string(serialized)

    assert rebuilt is not kernel
    assert isinstance(rebuilt, type(kernel))
    rebuilt.prepare(cosmology.cosmo)

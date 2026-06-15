#!/usr/bin/env python
#
# test_py_helper.py
#
# Thu Apr 17 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_py_helper.py
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

"""Unit tests for numcosmo_py.helper module.

Tests for utility functions in the helper module.
"""

import numpy as np
from numpy.testing import assert_allclose

from numcosmo_py import Ncm, Nc
from numcosmo_py.helper import duplicate_via_serialization, npa_to_seq
import numcosmo_py.cosmology as ncpy

Ncm.cfg_init()


def test_npa_to_seq_1d() -> None:
    """Test conversion of 1D NumPy array to sequence."""
    arr = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
    seq = npa_to_seq(arr)

    assert isinstance(seq, list)
    assert len(seq) == 5
    assert seq == [1.0, 2.0, 3.0, 4.0, 5.0]


def test_npa_to_seq_2d() -> None:
    """Test conversion of 2D NumPy array to sequence (flattened)."""
    arr = np.array([[1.0, 2.0], [3.0, 4.0]])
    seq = npa_to_seq(arr)

    assert isinstance(seq, list)
    assert len(seq) == 4
    assert seq == [1.0, 2.0, 3.0, 4.0]


def test_npa_to_seq_empty() -> None:
    """Test conversion of empty NumPy array to sequence."""
    arr = np.array([])
    seq = npa_to_seq(arr)

    assert isinstance(seq, list)
    assert len(seq) == 0
    assert seq == []


def test_duplicate_via_serialization_basic() -> None:
    """Test basic duplication of a NumCosmo object."""
    # Create a simple vector object
    vec = Ncm.Vector.new(5)
    for i in range(5):
        vec.set(i, float(i + 1))

    # Duplicate it
    vec_dup = duplicate_via_serialization(vec)

    # Verify it's a different object
    assert vec_dup is not None
    assert vec_dup is not vec
    assert isinstance(vec_dup, Ncm.Vector)

    # Verify the contents are the same
    assert vec_dup.len() == vec.len()
    for i in range(5):
        assert vec_dup.get(i) == vec.get(i)


def test_duplicate_via_serialization_type_preservation() -> None:
    """Test that duplicate_via_serialization preserves the object type."""
    # Create a matrix
    mat = Ncm.Matrix.new(3, 3)
    mat.set_identity()

    # Duplicate it
    mat_dup = duplicate_via_serialization(mat)

    # Type should be preserved
    assert isinstance(mat_dup, Ncm.Matrix)
    assert isinstance(mat_dup, type(mat))

    # Verify contents
    assert mat_dup.nrows() == mat.nrows()
    assert mat_dup.ncols() == mat.ncols()
    for i in range(3):
        for j in range(3):
            assert mat_dup.get(i, j) == mat.get(i, j)


def test_duplicate_via_serialization_with_custom_serializer() -> None:
    """Test duplication with a custom serializer instance."""
    vec = Ncm.Vector.new(3)
    vec.set_all(42.0)

    # Create custom serializer
    ser = Ncm.Serialize.new(Ncm.SerializeOpt.NONE)

    # Duplicate with custom serializer
    vec_dup = duplicate_via_serialization(vec, ser)

    assert vec_dup is not vec
    assert isinstance(vec_dup, Ncm.Vector)
    assert vec_dup.len() == 3
    assert vec_dup.get(0) == 42.0


def test_duplicate_via_serialization_complex_object() -> None:
    """Test duplication of a more complex object (cosmology object)."""
    # Create a cosmology object
    cosmo = ncpy.Cosmology.default()
    orig_h0 = cosmo.cosmo["H0"]

    # Duplicate it
    cosmo_dup = duplicate_via_serialization(cosmo.cosmo)

    # Verify it's different but equivalent
    assert cosmo_dup is not cosmo.cosmo
    assert isinstance(cosmo_dup, Nc.HICosmo)
    assert isinstance(cosmo_dup, Nc.HICosmoDE)
    assert isinstance(cosmo.cosmo, Nc.HICosmoDE)

    # Check a parameter value
    dup_h0 = cosmo_dup.props.H0
    assert_allclose(dup_h0, orig_h0, rtol=1e-15)

    # Modify the duplicate and ensure original is unchanged
    cosmo_dup.props.H0 = 75.0
    assert cosmo.cosmo.props.H0 == orig_h0
    assert cosmo_dup.props.H0 == 75.0


def test_duplicate_via_serialization_independence() -> None:
    """Test that duplicated objects are truly independent."""
    # Create a vector
    vec = Ncm.Vector.new(5)
    for i in range(5):
        vec.set(i, float(i))

    # Duplicate it
    vec_dup = duplicate_via_serialization(vec)

    # Modify the duplicate
    vec_dup.set(0, 999.0)

    # Original should be unchanged
    assert vec.get(0) == 0.0
    assert vec_dup.get(0) == 999.0


def test_duplicate_via_serialization_different_serializer_options() -> None:
    """Test that different serializer options both work correctly."""
    vec = Ncm.Vector.new(3)
    vec.set_all(7.0)

    # Test with CLEAN_DUP (default)
    ser_clean = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
    vec_dup_clean = duplicate_via_serialization(vec, ser_clean)

    # Test with NONE
    ser_none = Ncm.Serialize.new(Ncm.SerializeOpt.NONE)
    vec_dup_none = duplicate_via_serialization(vec, ser_none)

    # Both should work
    assert vec_dup_clean.get(0) == 7.0
    assert vec_dup_none.get(0) == 7.0
    assert vec_dup_clean is not vec
    assert vec_dup_none is not vec
    assert vec_dup_clean is not vec_dup_none

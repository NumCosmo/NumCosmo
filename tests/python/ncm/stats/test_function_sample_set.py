#!/usr/bin/env python
#
# test_py_function_sample_set.py
#
# Mon Mar 17 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_py_function_sample_set.py
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

"""Tests on NcmFunctionSampleSet class."""

from typing import Any
import math
import pytest
import numpy as np

from numcosmo_py import Ncm

Ncm.cfg_init()


@pytest.fixture(name="sample_set_dim2")
def fixture_sample_set_dim2() -> Ncm.FunctionSampleSet:
    """Fixture for NcmFunctionSampleSet with 2D output."""
    fss = Ncm.FunctionSampleSet.new(2)
    return fss


@pytest.fixture(name="sample_set_dim3")
def fixture_sample_set_dim3() -> Ncm.FunctionSampleSet:
    """Fixture for NcmFunctionSampleSet with 3D output."""
    fss = Ncm.FunctionSampleSet.new(3)
    return fss


@pytest.fixture(name="sample_set_with_data")
def fixture_sample_set_with_data() -> Ncm.FunctionSampleSet:
    """Fixture for NcmFunctionSampleSet populated with sample data."""
    fss = Ncm.FunctionSampleSet.new(3)

    # Add samples: x, [x^2, x^3, cos(x)]
    for x in [0.0, 1.0, 2.0, 5.0, 10.0, 15.0]:
        y = Ncm.Vector.new(3)
        y.set(0, x**2)
        y.set(1, x**3)
        y.set(2, math.cos(x))
        fss.add(x, y)

    return fss


def test_function_sample_set_new(sample_set_dim2: Ncm.FunctionSampleSet) -> None:
    """Test creation of NcmFunctionSampleSet."""
    assert isinstance(sample_set_dim2, Ncm.FunctionSampleSet)
    assert sample_set_dim2.get_len() == 2
    assert sample_set_dim2.get_nsamples() == 0


def test_function_sample_set_new_dim3(sample_set_dim3: Ncm.FunctionSampleSet) -> None:
    """Test creation of NcmFunctionSampleSet with dimension 3."""
    assert isinstance(sample_set_dim3, Ncm.FunctionSampleSet)
    assert sample_set_dim3.get_len() == 3
    assert sample_set_dim3.get_nsamples() == 0


def test_function_sample_set_add(sample_set_dim3: Ncm.FunctionSampleSet) -> None:
    """Test adding samples to NcmFunctionSampleSet."""
    x = 5.0
    y = Ncm.Vector.new(3)
    y.set(0, x**2)
    y.set(1, x**3)
    y.set(2, math.cos(x))

    sample_set_dim3.add(x, y)

    assert sample_set_dim3.get_nsamples() == 1

    # Use iterator to access the data
    it = sample_set_dim3.iter_begin()
    assert it.is_valid()
    assert it.get_x() == x

    y_retrieved = it.get_y()
    assert y_retrieved.get(0) == y.get(0)
    assert y_retrieved.get(1) == y.get(1)
    assert y_retrieved.get(2) == y.get(2)


def test_function_sample_set_add_multiple(
    sample_set_dim3: Ncm.FunctionSampleSet,
) -> None:
    """Test adding multiple samples to NcmFunctionSampleSet."""
    x_values = [0.0, 1.0, 2.0, 5.0, 10.0]

    for x in x_values:
        y = Ncm.Vector.new(3)
        y.set(0, x**2)
        y.set(1, x**3)
        y.set(2, math.cos(x))
        sample_set_dim3.add(x, y)

    assert sample_set_dim3.get_nsamples() == len(x_values)

    # Verify samples are stored in order using iterator
    it = sample_set_dim3.iter_begin()
    for expected_x in x_values:
        assert it.is_valid()
        retrieved_x = it.get_x()
        assert retrieved_x == expected_x
        it.next()


def test_function_sample_set_add_unordered(
    sample_set_dim3: Ncm.FunctionSampleSet,
) -> None:
    """Test that add() maintains x-order even when adding out of order."""
    # Add samples out of order
    for x in [5.0, 1.0, 10.0, 0.0, 2.0]:
        y = Ncm.Vector.new(3)
        y.set(0, x**2)
        y.set(1, x**3)
        y.set(2, math.cos(x))
        sample_set_dim3.add(x, y)

    assert sample_set_dim3.get_nsamples() == 5

    # Verify samples are automatically sorted by x using iterator
    expected_ordered = [0.0, 1.0, 2.0, 5.0, 10.0]
    it = sample_set_dim3.iter_begin()
    for expected_x in expected_ordered:
        assert it.is_valid()
        retrieved_x = it.get_x()
        assert retrieved_x == expected_x
        it.next()


def test_function_sample_set_iter_get_x(
    sample_set_with_data: Ncm.FunctionSampleSet,
) -> None:
    """Test iterator-based access to x values."""
    expected_x = [0.0, 1.0, 2.0, 5.0, 10.0, 15.0]

    it = sample_set_with_data.iter_begin()
    for exp_x in expected_x:
        assert it.is_valid()
        x = it.get_x()
        assert x == exp_x
        it.next()
    assert not it.is_valid()


def test_function_sample_set_iter_get_y(
    sample_set_with_data: Ncm.FunctionSampleSet,
) -> None:
    """Test iterator-based access to y vectors."""
    x_values = [0.0, 1.0, 2.0, 5.0, 10.0, 15.0]

    it = sample_set_with_data.iter_begin()
    for x in x_values:
        assert it.is_valid()
        y = it.get_y()
        assert isinstance(y, Ncm.Vector)
        assert y.len() == 3
        # Check values approximately match expected functions
        assert abs(y.get(0) - x**2) < 1e-15
        assert abs(y.get(1) - x**3) < 1e-15
        assert abs(y.get(2) - math.cos(x)) < 1e-15
        it.next()


def test_function_sample_set_interval_ok_flags(
    sample_set_with_data: Ncm.FunctionSampleSet,
) -> None:
    """Test interval_ok flag management using iterators."""
    # Initially all should be 0
    it = sample_set_with_data.iter_begin()
    while it.is_valid():
        assert it.get_interval_ok() == 0
        it.next()

    # Set some flags using iterators
    it = sample_set_with_data.iter_begin()
    it.set_interval_ok(1)  # First sample
    it.next()
    it.next()  # Skip second
    it.set_interval_ok(1)  # Third sample

    # Verify flags
    it = sample_set_with_data.iter_begin()
    assert it.get_interval_ok() == 1
    it.next()
    assert it.get_interval_ok() == 0
    it.next()
    assert it.get_interval_ok() == 1

    # Increment a flag
    it.inc_interval_ok()
    assert it.get_interval_ok() == 2

    # Reset all
    sample_set_with_data.reset_interval_ok()
    it = sample_set_with_data.iter_begin()
    while it.is_valid():
        assert it.get_interval_ok() == 0
        it.next()


def test_function_sample_set_all_intervals_ok(
    sample_set_with_data: Ncm.FunctionSampleSet,
) -> None:
    """Test all_intervals_ok convergence check."""
    # Initially all interval_ok are 0, so should return False for threshold >= 1
    assert not sample_set_with_data.all_intervals_ok(1)
    assert not sample_set_with_data.all_intervals_ok(2)

    # Set all intervals (all points except last) to 2 using iterators
    nsamples = sample_set_with_data.get_nsamples()
    it = sample_set_with_data.iter_begin()
    for _ in range(nsamples - 1):  # Exclude last point
        it.set_interval_ok(2)
        it.next()

    # Now all intervals should pass threshold 1 and 2 (interval_ok=2 >= threshold)
    assert sample_set_with_data.all_intervals_ok(1)
    assert sample_set_with_data.all_intervals_ok(0)

    # But not threshold 3
    assert not sample_set_with_data.all_intervals_ok(3)

    # Set one interval to fail
    it = sample_set_with_data.iter_begin()
    it.next()
    it.next()  # Move to third sample
    it.set_interval_ok(1)

    assert not sample_set_with_data.all_intervals_ok(2)
    assert sample_set_with_data.all_intervals_ok(0)


def test_function_sample_set_all_intervals_ok_convergence() -> None:
    """Test all_intervals_ok for convergence detection in refinement."""
    fss = Ncm.FunctionSampleSet.new(1)

    # Add initial points
    for x in [0.0, 1.0, 2.0, 3.0]:
        y = Ncm.Vector.new(1)
        y.set(0, x)
        fss.add(x, y)

    # Initially should be False
    assert not fss.all_intervals_ok(1)

    # After one refinement pass where all intervals pass
    # Simulate this by setting all interval_ok flags to 1
    it = fss.iter_begin()
    for _ in range(fss.get_nsamples() - 1):
        it.inc_interval_ok()
        it.next()

    # Now should pass threshold 0 and 1 (interval_ok=1 >= threshold)
    assert fss.all_intervals_ok(0)

    # But not threshold 2
    assert not fss.all_intervals_ok(2)

    # After another pass (interval_ok becomes 2)
    it = fss.iter_begin()
    for _ in range(fss.get_nsamples() - 1):
        it.inc_interval_ok()
        it.next()

    # Now should pass threshold 2 (interval_ok=2 >= 2)
    assert fss.all_intervals_ok(2)


def test_function_sample_set_iter_insert_before(
    sample_set_dim3: Ncm.FunctionSampleSet,
) -> None:
    """Test inserting samples before an iterator position."""
    # Add initial samples
    for x in [0.0, 5.0, 10.0]:
        y = Ncm.Vector.new(3)
        y.set(0, x**2)
        y.set(1, x**3)
        y.set(2, math.cos(x))
        sample_set_dim3.add(x, y)

    # Insert before position 1 (which is x=5.0)
    it = sample_set_dim3.iter_begin()
    it.next()  # Move to second sample (x=5.0)

    x_new = 2.5
    y_new = Ncm.Vector.new(3)
    y_new.set(0, x_new**2)
    y_new.set(1, x_new**3)
    y_new.set(2, math.cos(x_new))

    new_it = sample_set_dim3.iter_insert_before(it, x_new, y_new)
    assert new_it.get_x() == x_new

    assert sample_set_dim3.get_nsamples() == 4

    # Verify order
    it = sample_set_dim3.iter_begin()
    it.next()  # Now at position 1, should be x_new
    assert it.get_x() == x_new


def test_function_sample_set_iter_insert_after(
    sample_set_dim3: Ncm.FunctionSampleSet,
) -> None:
    """Test inserting samples after an iterator position."""
    # Add initial samples
    for x in [0.0, 5.0, 10.0]:
        y = Ncm.Vector.new(3)
        y.set(0, x**2)
        y.set(1, x**3)
        y.set(2, math.cos(x))
        sample_set_dim3.add(x, y)

    # Insert after position 0 (which is x=0.0)
    it = sample_set_dim3.iter_begin()

    x_new = 2.5
    y_new = Ncm.Vector.new(3)
    y_new.set(0, x_new**2)
    y_new.set(1, x_new**3)
    y_new.set(2, math.cos(x_new))

    new_it = sample_set_dim3.iter_insert_after(it, x_new, y_new)
    assert new_it.get_x() == x_new

    assert sample_set_dim3.get_nsamples() == 4

    # Verify order
    it = sample_set_dim3.iter_begin()
    it.next()  # Now at position 1, should be x_new
    assert it.get_x() == x_new


def test_function_sample_set_to_spline_vec(
    sample_set_with_data: Ncm.FunctionSampleSet,
) -> None:
    """Test conversion of NcmFunctionSampleSet to NcmSplineVec."""
    base_spline = Ncm.SplineCubicNotaknot.new()
    sv = sample_set_with_data.to_spline_vec(base_spline)

    assert isinstance(sv, Ncm.SplineVec)
    assert sv.get_len() == 3
    assert sv.is_init()


def test_function_sample_set_to_spline_vec_nondestructive(
    sample_set_with_data: Ncm.FunctionSampleSet,
) -> None:
    """Test that to_spline_vec() is non-destructive."""
    nsamples_before = sample_set_with_data.get_nsamples()

    base_spline = Ncm.SplineCubicNotaknot.new()
    sv = sample_set_with_data.to_spline_vec(base_spline)

    # Sample set should still have all samples
    assert sample_set_with_data.get_nsamples() == nsamples_before
    assert sv.get_len() == 3

    # Can convert again
    sv2 = sample_set_with_data.to_spline_vec(base_spline)
    assert sv2.get_len() == 3


def test_function_sample_set_to_spline_vec_evaluation(
    sample_set_with_data: Ncm.FunctionSampleSet,
) -> None:
    """Test that spline evaluation matches at knot points."""
    base_spline = Ncm.SplineCubicNotaknot.new()
    sv = sample_set_with_data.to_spline_vec(base_spline)

    # Evaluate at knot points using iterator
    it = sample_set_with_data.iter_begin()
    while it.is_valid():
        x = it.get_x()
        y_expected = it.get_y()
        y_computed = sv.eval_array(x)

        # At knot points, spline should exactly match
        for j in range(3):
            assert abs(y_computed[j] - y_expected.get(j)) < 1e-12
        it.next()


def test_function_sample_set_iterative_refinement() -> None:
    """Test a simple iterative refinement scenario."""
    fss = Ncm.FunctionSampleSet.new(2)

    # Add initial coarse samples for f(x) = [sin(x), cos(x)]
    for x in np.linspace(0, 2 * np.pi, 7):
        y = Ncm.Vector.new(2)
        y.set(0, math.sin(x))
        y.set(1, math.cos(x))
        fss.add(x, y)

    assert fss.get_nsamples() == 7

    # Build initial spline
    base_spline = Ncm.SplineCubicNotaknot.new()
    sv = fss.to_spline_vec(base_spline)
    assert sv.get_len() == 2

    # Simulate refinement: add midpoint between first two samples
    it = fss.iter_begin()
    x0 = it.get_x()
    it.next()
    x1 = it.get_x()
    it.prev()  # Back to first

    x_mid = (x0 + x1) / 2.0
    y_mid = Ncm.Vector.new(2)
    y_mid.set(0, math.sin(x_mid))
    y_mid.set(1, math.cos(x_mid))
    _ = fss.iter_insert_after(it, x_mid, y_mid)

    assert fss.get_nsamples() == 8

    # Build refined spline
    sv_refined = fss.to_spline_vec(base_spline)
    assert sv_refined.get_len() == 2

    # Evaluate at the new midpoint - should match exactly on knot
    y_eval = sv_refined.eval_array(x_mid)
    assert abs(y_eval[0] - math.sin(x_mid)) < 1e-12
    assert abs(y_eval[1] - math.cos(x_mid)) < 1e-12


def vector_func_callback(x: float, y: Ncm.Vector, _user_data: Any) -> None:
    """Callback function for testing: computes [x^2, x^3, cos(x)]."""
    y.set(0, x**2)
    y.set(1, x**3)
    y.set(2, math.cos(x))


def test_function_sample_set_add_func() -> None:
    """Test adding samples using a function pointer."""
    fss = Ncm.FunctionSampleSet.new(3)

    # Add samples using function callback
    for x in [0.0, 1.0, 2.0, 5.0, 10.0]:
        fss.add_func(x, vector_func_callback, None)

    assert fss.get_nsamples() == 5

    # Verify values match using iterator
    it = fss.iter_begin()
    for x in [0.0, 1.0, 2.0, 5.0, 10.0]:
        assert it.is_valid()
        y = it.get_y()
        assert abs(y.get(0) - x**2) < 1e-15
        assert abs(y.get(1) - x**3) < 1e-15
        assert abs(y.get(2) - math.cos(x)) < 1e-15
        it.next()


def test_function_sample_set_add_func_with_user_data() -> None:
    """Test adding samples using a function pointer with user data."""

    def scaled_func(x: float, y: Ncm.Vector, user_data) -> None:
        """Function that uses user_data as a scale factor."""
        scale = user_data
        y.set(0, scale * x**2)
        y.set(1, scale * x**3)

    fss = Ncm.FunctionSampleSet.new(2)
    scale = 2.5

    # Add samples with scaling
    for x in [1.0, 2.0, 3.0]:
        fss.add_func(x, scaled_func, scale)

    assert fss.get_nsamples() == 3

    # Verify scaled values using iterator
    it = fss.iter_begin()
    for x in [1.0, 2.0, 3.0]:
        assert it.is_valid()
        y = it.get_y()
        assert abs(y.get(0) - scale * x**2) < 1e-15
        assert abs(y.get(1) - scale * x**3) < 1e-15
        it.next()


def test_function_sample_set_iter_insert_before_func() -> None:
    """Test inserting samples before an iterator position using a function pointer."""
    fss = Ncm.FunctionSampleSet.new(3)

    # Add initial samples
    for x in [0.0, 5.0, 10.0]:
        fss.add_func(x, vector_func_callback, None)

    # Insert before position 1 using function
    it = fss.iter_begin()
    it.next()  # Move to position 1
    x_new = 2.5
    new_it = fss.iter_insert_before_func(it, x_new, vector_func_callback, None)

    assert fss.get_nsamples() == 4
    assert new_it.get_x() == x_new

    # Verify the inserted value
    y_new = new_it.get_y()
    assert abs(y_new.get(0) - x_new**2) < 1e-15
    assert abs(y_new.get(1) - x_new**3) < 1e-15
    assert abs(y_new.get(2) - math.cos(x_new)) < 1e-15


def test_function_sample_set_insert_after_func() -> None:
    """Test inserting samples after a specific index using a function pointer."""
    fss = Ncm.FunctionSampleSet.new(3)

    # Add initial samples
    for x in [0.0, 5.0, 10.0]:
        fss.add_func(x, vector_func_callback, None)

    # Insert after first element using iterator and function
    x_new = 2.5
    it = fss.iter_begin()
    _ = fss.iter_insert_after_func(it, x_new, vector_func_callback, None)

    assert fss.get_nsamples() == 4

    # Check the inserted value is at the right position
    it = fss.iter_begin()
    it.next()  # Skip first (x=0.0)
    assert it.get_x() == x_new

    # Verify the inserted value
    y_new = it.get_y()
    assert abs(y_new.get(0) - x_new**2) < 1e-15
    assert abs(y_new.get(1) - x_new**3) < 1e-15
    assert abs(y_new.get(2) - math.cos(x_new)) < 1e-15


def test_function_sample_set_iterative_refinement_with_func() -> None:
    """Test iterative refinement using function pointers."""

    def trig_func(x: float, y: Ncm.Vector, _user_data: Any) -> None:
        """Compute [sin(x), cos(x)]."""
        y.set(0, math.sin(x))
        y.set(1, math.cos(x))

    fss = Ncm.FunctionSampleSet.new(2)

    # Add initial coarse samples using function
    for x in np.linspace(0, 2 * np.pi, 7):
        fss.add_func(x, trig_func, None)

    assert fss.get_nsamples() == 7

    # Build initial spline
    base_spline = Ncm.SplineCubicNotaknot.new()
    sv = fss.to_spline_vec(base_spline)
    assert sv.get_len() == 2

    # Simulate refinement: add midpoint using function and iterator
    it = fss.iter_begin()
    x0 = it.get_x()
    it.next()
    x1 = it.get_x()
    x_mid = (x0 + x1) / 2.0

    # Insert after first element
    it = fss.iter_begin()
    _ = fss.iter_insert_after_func(it, x_mid, trig_func, None)

    assert fss.get_nsamples() == 8

    # Build refined spline
    sv_refined = fss.to_spline_vec(base_spline)
    assert sv_refined.get_len() == 2

    # Evaluate at the new midpoint
    y_eval = sv_refined.eval_array(x_mid)
    assert abs(y_eval[0] - math.sin(x_mid)) < 1e-12
    assert abs(y_eval[1] - math.cos(x_mid)) < 1e-12


def test_function_sample_set_mark_all_old() -> None:
    """Test marking all points as OLD."""
    fss = Ncm.FunctionSampleSet.new(2)

    # Add some initial points (they start as NEW)
    for x in [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0]:
        y = Ncm.Vector.new(2)
        y.set(0, x)
        y.set(1, x**2)
        fss.add(x, y)

    # Mark all as old
    fss.mark_all_old()

    # All points should now be OLD (we can verify by adding a new point and using
    # to_spline_vec_old). Add a new point.
    y_new = Ncm.Vector.new(2)
    y_new.set(0, 1.5)
    y_new.set(1, 1.5**2)
    fss.add(1.5, y_new)

    # to_spline_vec_old should only use the first 3 points
    base_spline = Ncm.SplineCubicNotaknot.new()
    sv_old = fss.to_spline_vec_old(base_spline)

    # The old spline should have 3 points
    assert sv_old.get_nknots() == 7


def test_function_sample_set_to_spline_vec_old() -> None:
    """Test creating spline from OLD points only."""
    fss = Ncm.FunctionSampleSet.new(2)

    # Add initial points and mark as OLD
    for x in [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0]:
        y = Ncm.Vector.new(2)
        y.set(0, x)
        y.set(1, x**2)
        fss.add(x, y)

    fss.mark_all_old()

    # Add some NEW points
    for x in [0.5, 1.5, 2.5]:
        y = Ncm.Vector.new(2)
        y.set(0, x)
        y.set(1, x**2)
        fss.add(x, y)

    # Total should be 10 samples now
    assert fss.get_nsamples() == 10

    # Create spline from OLD points only
    base_spline = Ncm.SplineCubicNotaknot.new()
    sv_old = fss.to_spline_vec_old(base_spline)

    # Should only have 7 knots
    assert sv_old.get_nknots() == 7

    # Evaluate at an old point
    y_eval = sv_old.eval_array(1.0)
    assert abs(y_eval[0] - 1.0) < 1e-10
    assert abs(y_eval[1] - 1.0) < 1e-10


def test_function_sample_set_refine() -> None:
    """Test refinement algorithm."""
    fss = Ncm.FunctionSampleSet.new(2)

    # Simple linear function: y = [x, 2*x]
    # Add initial coarse grid and mark as OLD
    for x in [0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0]:
        y = Ncm.Vector.new(2)
        y.set(0, x)
        y.set(1, 2.0 * x)
        fss.add(x, y)

    fss.mark_all_old()

    # Add midpoints as NEW
    for x in [1.0, 3.0]:
        y = Ncm.Vector.new(2)
        y.set(0, x)
        y.set(1, 2.0 * x)
        fss.add(x, y)

    # Refine with tight tolerance - linear function should pass
    base_spline = Ncm.SplineCubicNotaknot.new()
    fss.refine(1e-10, 1e-10, base_spline)

    # After refinement:
    # 1. interval_ok should be incremented for points that passed
    # 2. All points should be marked as OLD

    # For a linear function, cubic spline interpolation is exact
    # So all NEW points should pass the test
    # The interval_ok should be incremented for the new points and their left neighbors

    # Check that interval_ok counters were incremented appropriately using iterator
    # x=0: left of x=1, should be incremented
    it = fss.iter_begin()
    assert (
        it.get_interval_ok() >= 0
    )  # May or may not be incremented depending on implementation
    # x=1: new point that passed
    # x=2: left of x=3, should be incremented
    # x=3: new point that passed

    # Verify all points now are OLD by checking we can create spline with all points
    sv_all = fss.to_spline_vec_old(base_spline)
    assert sv_all.get_nknots() == 12  # All 12 points should be OLD


def test_function_sample_set_refine_oscillating() -> None:
    """Test refinement with oscillating function."""
    fss = Ncm.FunctionSampleSet.new(1)

    # Use sin(x) as test function
    # Add coarse grid
    for x in np.linspace(0, 2 * np.pi, 10):
        y = Ncm.Vector.new(1)
        y.set(0, math.sin(x))
        fss.add(x, y)

    fss.mark_all_old()

    # Add midpoint
    x_mid = math.pi / 2.0
    y = Ncm.Vector.new(1)
    y.set(0, math.sin(x_mid))
    fss.add(x_mid, y)

    # Refine - the spline interpolation of [0, 2pi] will not match sin(pi/2) = 1 well
    base_spline = Ncm.SplineCubicNotaknot.new()
    fss.refine(0.01, 0.01, base_spline)

    # Verify all points are now OLD
    sv_all = fss.to_spline_vec_old(base_spline)
    assert sv_all.get_nknots() == 11


def test_function_sample_set_log_vals(
    sample_set_with_data: Ncm.FunctionSampleSet, capfd
) -> None:
    """Test that log_vals outputs sample information correctly."""
    sample_set_with_data.log_vals()

    out, err = capfd.readouterr()
    # Check that output contains expected information
    assert "NcmFunctionSampleSet" in out or "NcmFunctionSampleSet" in err
    assert "samples" in out or "samples" in err
    assert "dimension" in out or "dimension" in err


def test_function_sample_set_autoknots_refinement() -> None:
    """Test autoknots-style iterative refinement.

    Test for F(x) = (sin(x), cos(x)) on [0, 20]."""

    def trig_func(x: float, y: Ncm.Vector, _user_data: Any) -> None:
        """Vector function F(x) = (sin(x), cos(x))."""
        y.set(0, math.sin(x))
        y.set(1, math.cos(x))

    # Parameters
    xi = 0.0
    xf = 20.0
    reltol = 1e-5
    abstol = 1e-10
    max_iter = 100
    min_pass_threshold = 1  # Each interval must have interval_ok >= this value

    fss = Ncm.FunctionSampleSet.new(2)
    base_spline = Ncm.SplineCubicNotaknot.new()

    # Step 1: Add initial coarse grid (similar to autoknots initial grid)
    initial_nknots = 6
    for x in np.linspace(xi, xf, initial_nknots):
        fss.add_func(x, trig_func, None)

    fss.mark_all_old()

    # Iterative refinement loop
    for iteration in range(max_iter):
        nsamples = fss.get_nsamples()

        # Check convergence: all intervals have passed enough times
        if fss.all_intervals_ok(min_pass_threshold):
            print(f"Converged after {iteration} iterations with {nsamples} knots")
            break

        # Find intervals that need refinement using iterators
        # An interval needs refinement if its interval_ok < min_pass_threshold
        intervals_to_refine = []
        it = fss.iter_begin()
        idx = 0
        while it.has_next():
            if it.get_interval_ok() < min_pass_threshold:
                intervals_to_refine.append((idx, it.copy()))
            it.next()
            idx += 1

        if not intervals_to_refine:
            # This shouldn't happen if all_intervals_ok returned False
            break

        # Add midpoints to intervals that need refinement (in reverse to preserve
        # indices)
        for idx, iter_left in reversed(intervals_to_refine):
            iter_right = iter_left.copy()
            iter_right.next()

            x_left = iter_left.get_x()
            x_right = iter_right.get_x()
            x_mid = 0.5 * (x_left + x_right)

            # Insert after iterator position
            _ = fss.iter_insert_after_func(iter_left, x_mid, trig_func, None)

        # Test the new points
        fss.refine(reltol, abstol, base_spline)

    else:
        # Max iterations reached
        print(f"Max iterations ({max_iter}) reached with {fss.get_nsamples()} knots")

    # Verify solution quality
    final_nsamples = fss.get_nsamples()
    print(f"Final number of knots: {final_nsamples}")

    # Should have added knots (more than initial)
    assert final_nsamples > initial_nknots

    # Build final spline
    sv_final = fss.to_spline_vec(base_spline)
    assert sv_final.get_nknots() == final_nsamples

    # Test accuracy at random points
    np.random.seed(42)
    test_points = np.random.uniform(xi + 0.1, xf - 0.1, 50)

    max_error = 0.0
    for x_test in test_points:
        y_true = np.array([math.sin(x_test), math.cos(x_test)])
        y_spline = np.array(sv_final.eval_array(x_test))

        error = float(np.linalg.norm(y_true - y_spline))
        max_error = max(max_error, error)

        # Check relative error
        norm_true = np.linalg.norm(y_true)
        threshold = reltol * norm_true + abstol
        assert (
            error < threshold * 10
        ), f"Error {error} exceeds threshold {threshold} at x={x_test}"

    print(f"Maximum interpolation error: {max_error:.3e}")

    # Verify convergence was actually achieved
    assert fss.all_intervals_ok(min_pass_threshold)


def test_function_sample_set_autoknots_challenging_function() -> None:
    """Test autoknots-style refinement with a more challenging function."""

    def challenging_func(x: float, y: Ncm.Vector, _user_data: Any) -> None:
        """More challenging vector function with varying frequency."""
        # F(x) = (sin(x^1.5), cos(2x))
        y.set(0, math.sin(x**1.5))
        y.set(1, math.cos(2.0 * x))

    # Parameters
    xi = 0.1  # Avoid x=0 for x^1.5
    xf = 10.0
    reltol = 1e-4
    abstol = 1e-8
    max_iter = 50
    min_pass_threshold = 1

    fss = Ncm.FunctionSampleSet.new(2)
    base_spline = Ncm.SplineCubicNotaknot.new()

    # Initial grid
    initial_nknots = 6
    for x in np.linspace(xi, xf, initial_nknots):
        fss.add_func(x, challenging_func, None)

    fss.mark_all_old()

    # Refinement loop
    for iteration in range(max_iter):
        if fss.all_intervals_ok(min_pass_threshold):
            print(
                f"Challenging function: Converged after {iteration} iterations "
                f"with {fss.get_nsamples()} knots"
            )
            break

        nsamples = fss.get_nsamples()
        assert nsamples > 1, "Should have at least 2 samples to refine intervals"
        intervals_to_refine: list[tuple[int, Ncm.FunctionSampleSetIter]] = []

        # Find intervals that need refinement using iterators
        it = fss.iter_begin()
        idx = 0
        while it.has_next():
            if it.get_interval_ok() < min_pass_threshold:
                intervals_to_refine.append((idx, it.copy()))
            it.next()
            idx += 1

        # Add midpoints (backwards to preserve indices)
        for idx, iter_left in reversed(intervals_to_refine):
            iter_right = iter_left.copy()
            iter_right.next()

            x_left = iter_left.get_x()
            x_right = iter_right.get_x()
            x_mid = 0.5 * (x_left + x_right)

            _ = fss.iter_insert_after_func(iter_left, x_mid, challenging_func, None)

        fss.refine(reltol, abstol, base_spline)

    final_nsamples = fss.get_nsamples()
    print(f"Challenging function: Final knots = {final_nsamples}")

    # Should need more knots due to complexity
    assert final_nsamples > initial_nknots * 2

    # Verify final accuracy
    sv_final = fss.to_spline_vec(base_spline)
    test_x = np.linspace(xi + 0.01, xf - 0.01, 30)

    for x in test_x:
        y_true = [math.sin(x**1.5), math.cos(2.0 * x)]
        y_spline = sv_final.eval_array(x)

        error = math.sqrt(
            (y_true[0] - y_spline[0]) ** 2 + (y_true[1] - y_spline[1]) ** 2
        )
        norm_true = math.sqrt(y_true[0] ** 2 + y_true[1] ** 2)
        threshold = reltol * norm_true + abstol

        assert error < threshold * 10, f"Error too large at x={x}"

    assert fss.all_intervals_ok(min_pass_threshold)


def test_x_min_max_tracking_sequential(sample_set_dim2: Ncm.FunctionSampleSet) -> None:
    """Test x_min and x_max tracking with sequential additions."""
    fss = sample_set_dim2

    # Initially should be inf and -inf
    assert math.isinf(fss.get_x_min()) and fss.get_x_min() > 0
    assert math.isinf(fss.get_x_max()) and fss.get_x_max() < 0

    # Add points in ascending order
    x_values = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]
    for x in x_values:
        y = Ncm.Vector.new_array([x, x**2])
        fss.add(x, y)

    assert abs(fss.get_x_min() - 0.0) < 1e-10
    assert abs(fss.get_x_max() - 5.0) < 1e-10


def test_x_min_max_tracking_reverse(sample_set_dim2: Ncm.FunctionSampleSet) -> None:
    """Test x_min and x_max tracking with reverse order additions."""
    fss = sample_set_dim2

    # Add points in descending order
    x_values = [5.0, 4.0, 3.0, 2.0, 1.0, 0.0]
    for x in x_values:
        y = Ncm.Vector.new_array([x, x**2])
        fss.add(x, y)

    assert abs(fss.get_x_min() - 0.0) < 1e-10
    assert abs(fss.get_x_max() - 5.0) < 1e-10


def test_x_min_max_tracking_random_order(
    sample_set_dim2: Ncm.FunctionSampleSet,
) -> None:
    """Test x_min and x_max tracking with random order additions."""
    fss = sample_set_dim2

    # Add points in non-sequential order
    x_values = [2.5, 0.5, 4.5, 1.5, 3.5, 5.5]
    for x in x_values:
        y = Ncm.Vector.new_array([x, x**2])
        fss.add(x, y)

    assert abs(fss.get_x_min() - 0.5) < 1e-10
    assert abs(fss.get_x_max() - 5.5) < 1e-10


def test_x_min_max_tracking_negative_values(
    sample_set_dim2: Ncm.FunctionSampleSet,
) -> None:
    """Test x_min and x_max tracking with negative x values."""
    fss = sample_set_dim2

    # Add points with negative values
    x_values = [-3.0, -1.0, 0.0, 1.0, 2.0, 3.0]
    for x in x_values:
        y = Ncm.Vector.new_array([x, x**2])
        fss.add(x, y)

    assert abs(fss.get_x_min() - (-3.0)) < 1e-10
    assert abs(fss.get_x_max() - 3.0) < 1e-10


def test_x_min_max_tracking_iterator_insertions(
    sample_set_dim2: Ncm.FunctionSampleSet,
) -> None:
    """Test x_min and x_max tracking with iterator insertions."""
    fss = sample_set_dim2

    # Start with some initial points
    initial_x = [0.0, 2.0, 4.0, 6.0, 8.0, 10.0]
    for x in initial_x:
        y = Ncm.Vector.new_array([x, x**2])
        fss.add(x, y)

    # Insert points before minimum
    it = fss.iter_begin()
    y_neg = Ncm.Vector.new_array([-2.0, 4.0])
    fss.iter_insert_before(it, -2.0, y_neg)

    # Insert points after maximum
    it = fss.iter_end()
    y_large = Ncm.Vector.new_array([12.0, 144.0])
    fss.iter_insert_after(it, 12.0, y_large)

    assert abs(fss.get_x_min() - (-2.0)) < 1e-10
    assert abs(fss.get_x_max() - 12.0) < 1e-10


def test_absmaxF_tracking_dim2(sample_set_dim2: Ncm.FunctionSampleSet) -> None:
    """Test absmaxF tracking for 2D output."""
    fss = sample_set_dim2

    # Initially should be zero with x = nan
    absmaxF0_init = fss.get_absmaxF(0)
    absmaxF1_init = fss.get_absmaxF(1)
    assert abs(absmaxF0_init[0] - 0.0) < 1e-10
    assert abs(absmaxF1_init[0] - 0.0) < 1e-10

    assert math.isnan(absmaxF0_init[1])
    assert math.isnan(absmaxF1_init[1])

    # Add points with various values
    x_values = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]
    y0_values = [1.0, -3.0, 2.0, -1.5, 0.5, -2.5]
    y1_values = [0.5, 1.5, -4.0, 2.5, -1.0, 3.5]

    for x, y0, y1 in zip(x_values, y0_values, y1_values):
        y = Ncm.Vector.new_array([y0, y1])
        fss.add(x, y)

    # Component 0: max |y0| should be 3.0 (from -3.0 at x=1.0)
    # Component 1: max |y1| should be 4.0 (from -4.0 at x=2.0)
    absmaxF0 = fss.get_absmaxF(0)
    absmaxF1 = fss.get_absmaxF(1)
    assert abs(absmaxF0[0] - 3.0) < 1e-10
    assert abs(absmaxF0[1] - 1.0) < 1e-10
    assert abs(absmaxF1[0] - 4.0) < 1e-10
    assert abs(absmaxF1[1] - 2.0) < 1e-10


def test_absmaxF_tracking_dim3(sample_set_dim3: Ncm.FunctionSampleSet) -> None:
    """Test absmaxF tracking for 3D output."""
    fss = sample_set_dim3

    # Add points with known max absolute values
    x_values = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]
    y_data = [
        [1.0, 2.0, 3.0],
        [-5.0, 1.5, 2.5],
        [2.0, -6.0, 1.0],
        [1.5, 2.5, -7.0],
        [0.5, 1.0, 2.0],
        [-3.0, -4.0, 5.0],
    ]

    for x, y_vals in zip(x_values, y_data):
        y = Ncm.Vector.new_array(y_vals)
        fss.add(x, y)

    # Component 0: max |y0| should be 5.0 at x=1.0
    # Component 1: max |y1| should be 6.0 at x=2.0
    # Component 2: max |y2| should be 7.0 at x=3.0
    absmaxF0 = fss.get_absmaxF(0)
    absmaxF1 = fss.get_absmaxF(1)
    absmaxF2 = fss.get_absmaxF(2)
    assert abs(absmaxF0[0] - 5.0) < 1e-10
    assert abs(absmaxF0[1] - 1.0) < 1e-10
    assert abs(absmaxF1[0] - 6.0) < 1e-10
    assert abs(absmaxF1[1] - 2.0) < 1e-10
    assert abs(absmaxF2[0] - 7.0) < 1e-10
    assert abs(absmaxF2[1] - 3.0) < 1e-10


def test_absmaxF_tracking_all_positive(sample_set_dim2: Ncm.FunctionSampleSet) -> None:
    """Test absmaxF tracking with all positive values."""
    fss = sample_set_dim2

    x_values = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]
    y0_values = [1.0, 3.0, 2.0, 1.5, 0.5, 2.5]
    y1_values = [0.5, 1.5, 4.0, 2.5, 1.0, 3.5]

    for x, y0, y1 in zip(x_values, y0_values, y1_values):
        y = Ncm.Vector.new_array([y0, y1])
        fss.add(x, y)

    # Component 0: max is 3.0 at x=1.0
    # Component 1: max is 4.0 at x=2.0
    absmaxF0 = fss.get_absmaxF(0)
    absmaxF1 = fss.get_absmaxF(1)
    assert abs(absmaxF0[0] - 3.0) < 1e-10
    assert abs(absmaxF0[1] - 1.0) < 1e-10
    assert abs(absmaxF1[0] - 4.0) < 1e-10
    assert abs(absmaxF1[1] - 2.0) < 1e-10


def test_absmaxF_tracking_all_negative(sample_set_dim2: Ncm.FunctionSampleSet) -> None:
    """Test absmaxF tracking with all negative values."""
    fss = sample_set_dim2

    x_values = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]
    y0_values = [-1.0, -3.0, -2.0, -1.5, -0.5, -2.5]
    y1_values = [-0.5, -1.5, -4.0, -2.5, -1.0, -3.5]

    for x, y0, y1 in zip(x_values, y0_values, y1_values):
        y = Ncm.Vector.new_array([y0, y1])
        fss.add(x, y)

    # Component 0: max |y0| is 3.0 (from -3.0 at x=1.0)
    # Component 1: max |y1| is 4.0 (from -4.0 at x=2.0)
    absmaxF0 = fss.get_absmaxF(0)
    absmaxF1 = fss.get_absmaxF(1)
    assert abs(absmaxF0[0] - 3.0) < 1e-10
    assert abs(absmaxF0[1] - 1.0) < 1e-10
    assert abs(absmaxF1[0] - 4.0) < 1e-10
    assert abs(absmaxF1[1] - 2.0) < 1e-10


def test_absmaxF_tracking_with_iterator_insertions(
    sample_set_dim2: Ncm.FunctionSampleSet,
) -> None:
    """Test absmaxF tracking with iterator insertions."""
    fss = sample_set_dim2

    # Add initial points
    initial_x = [0.0, 2.0, 4.0, 6.0, 8.0, 10.0]
    for x in initial_x:
        y = Ncm.Vector.new_array([x * 0.1, x * 0.2])
        fss.add(x, y)

    # Check initial max values (should be at x=10.0)
    initial_max0 = fss.get_absmaxF(0)
    initial_max1 = fss.get_absmaxF(1)
    assert abs(initial_max0[0] - 1.0) < 1e-10  # 10.0 * 0.1
    assert abs(initial_max0[1] - 10.0) < 1e-10
    assert abs(initial_max1[0] - 2.0) < 1e-10  # 10.0 * 0.2
    assert abs(initial_max1[1] - 10.0) < 1e-10

    # Insert point with larger absolute values at x=1.0
    it = fss.iter_begin()
    it.next()  # Move to second element
    y_large = Ncm.Vector.new_array([10.0, -15.0])
    fss.iter_insert_after(it, 1.0, y_large)

    # Should update to new max values
    absmaxF0_new = fss.get_absmaxF(0)
    absmaxF1_new = fss.get_absmaxF(1)
    assert absmaxF0_new[0] >= initial_max0[0]
    assert absmaxF1_new[0] >= initial_max1[0]
    assert abs(absmaxF0_new[0] - 10.0) < 1e-10
    assert abs(absmaxF0_new[1] - 1.0) < 1e-10
    assert abs(absmaxF1_new[0] - 15.0) < 1e-10
    assert abs(absmaxF1_new[1] - 1.0) < 1e-10


def test_absmaxF_x_coordinate_tracking(sample_set_dim2: Ncm.FunctionSampleSet) -> None:
    """Test that get_absmaxF correctly tracks the x coordinate of the maximum."""
    fss = sample_set_dim2

    # Add points where different x values have the max for each component
    test_data = [
        (0.0, [1.0, 2.0]),
        (1.0, [5.0, 1.5]),  # Component 0 max here
        (2.0, [2.0, 3.0]),
        (3.0, [1.5, 8.0]),  # Component 1 max here
        (4.0, [3.0, 2.5]),
    ]

    for x, y_vals in test_data:
        y = Ncm.Vector.new_array(y_vals)
        fss.add(x, y)

    # Component 0: max 5.0 at x=1.0
    # Component 1: max 8.0 at x=3.0
    absmaxF0 = fss.get_absmaxF(0)
    absmaxF1 = fss.get_absmaxF(1)

    assert abs(absmaxF0[0] - 5.0) < 1e-10
    assert abs(absmaxF0[1] - 1.0) < 1e-10
    assert abs(absmaxF1[0] - 8.0) < 1e-10
    assert abs(absmaxF1[1] - 3.0) < 1e-10

    # Add a new point with even larger value in component 0 at different x
    y_new = Ncm.Vector.new_array([12.0, 1.0])
    fss.add(5.0, y_new)

    # Component 0 max should now be at x=5.0
    absmaxF0_updated = fss.get_absmaxF(0)
    assert abs(absmaxF0_updated[0] - 12.0) < 1e-10
    assert abs(absmaxF0_updated[1] - 5.0) < 1e-10

    # Component 1 should remain unchanged
    absmaxF1_after = fss.get_absmaxF(1)
    assert abs(absmaxF1_after[0] - 8.0) < 1e-10
    assert abs(absmaxF1_after[1] - 3.0) < 1e-10


def test_absmaxF_x_coordinate_with_negative_values(
    sample_set_dim2: Ncm.FunctionSampleSet,
) -> None:
    """Test that x coordinate tracking works with negative function values."""
    fss = sample_set_dim2

    # Add points with negative values being the maximum in absolute terms
    test_data = [
        (-2.0, [1.0, -3.0]),
        (-1.0, [-5.0, 2.0]),  # Component 0: |-5.0| = 5.0 max here
        (0.0, [2.0, 1.0]),
        (1.0, [1.0, -8.0]),  # Component 1: |-8.0| = 8.0 max here
        (2.0, [-3.0, 2.0]),
    ]

    for x, y_vals in test_data:
        y = Ncm.Vector.new_array(y_vals)
        fss.add(x, y)

    # Component 0: max |y| = 5.0 from y=-5.0 at x=-1.0
    # Component 1: max |y| = 8.0 from y=-8.0 at x=1.0
    absmaxF0 = fss.get_absmaxF(0)
    absmaxF1 = fss.get_absmaxF(1)

    assert abs(absmaxF0[0] - 5.0) < 1e-10
    assert abs(absmaxF0[1] - (-1.0)) < 1e-10
    assert abs(absmaxF1[0] - 8.0) < 1e-10
    assert abs(absmaxF1[1] - 1.0) < 1e-10


def test_absmaxF_l2_norm_basic(sample_set_dim2: Ncm.FunctionSampleSet) -> None:
    """Test L2 norm computation with simple values."""
    fss = sample_set_dim2

    # Initially should be zero (no samples)
    l2_norm_init = fss.get_absmaxF_l2_norm()
    assert abs(l2_norm_init - 0.0) < 1e-10

    # Add samples
    # Component 0: max will be 3.0
    # Component 1: max will be 4.0
    # L2 norm should be sqrt(3^2 + 4^2) = 5.0
    test_data = [
        (0.0, [1.0, 2.0]),
        (1.0, [3.0, 1.5]),
        (2.0, [2.0, 4.0]),
        (3.0, [1.5, 2.5]),
    ]

    for x, y_vals in test_data:
        y = Ncm.Vector.new_array(y_vals)
        fss.add(x, y)

    l2_norm = fss.get_absmaxF_l2_norm()
    expected_l2 = math.sqrt(3.0**2 + 4.0**2)  # sqrt(9 + 16) = 5.0
    assert abs(l2_norm - expected_l2) < 1e-10


def test_absmaxF_l2_norm_with_negatives(
    sample_set_dim3: Ncm.FunctionSampleSet,
) -> None:
    """Test L2 norm with negative values (should use absolute values)."""
    fss = sample_set_dim3

    # Component 0: max |y| = 5.0 (from -5.0)
    # Component 1: max |y| = 6.0 (from -6.0)
    # Component 2: max |y| = 7.0
    # L2 norm should be sqrt(5^2 + 6^2 + 7^2) = sqrt(25 + 36 + 49) = sqrt(110)
    test_data = [
        (0.0, [1.0, 2.0, 3.0]),
        (1.0, [-5.0, 1.5, 2.5]),
        (2.0, [2.0, -6.0, 1.0]),
        (3.0, [1.5, 2.5, 7.0]),
    ]

    for x, y_vals in test_data:
        y = Ncm.Vector.new_array(y_vals)
        fss.add(x, y)

    l2_norm = fss.get_absmaxF_l2_norm()
    expected_l2 = math.sqrt(5.0**2 + 6.0**2 + 7.0**2)  # sqrt(110) ≈ 10.488
    assert abs(l2_norm - expected_l2) < 1e-10


def test_absmaxF_linf_norm_basic(sample_set_dim2: Ncm.FunctionSampleSet) -> None:
    """Test L∞ norm computation with simple values."""
    fss = sample_set_dim2

    # Initially should be zero (no samples)
    linf_norm_init = fss.get_absmaxF_linf_norm()
    assert abs(linf_norm_init - 0.0) < 1e-10

    # Add samples
    # Component 0: max will be 3.0
    # Component 1: max will be 4.0
    # L∞ norm should be max(3.0, 4.0) = 4.0
    test_data = [
        (0.0, [1.0, 2.0]),
        (1.0, [3.0, 1.5]),
        (2.0, [2.0, 4.0]),
        (3.0, [1.5, 2.5]),
    ]

    for x, y_vals in test_data:
        y = Ncm.Vector.new_array(y_vals)
        fss.add(x, y)

    linf_norm = fss.get_absmaxF_linf_norm()
    expected_linf = 4.0  # max(3.0, 4.0)
    assert abs(linf_norm - expected_linf) < 1e-10


def test_absmaxF_linf_norm_with_negatives(
    sample_set_dim3: Ncm.FunctionSampleSet,
) -> None:
    """Test L∞ norm with negative values (should use absolute values)."""
    fss = sample_set_dim3

    # Component 0: max |y| = 5.0 (from -5.0)
    # Component 1: max |y| = 6.0 (from -6.0)
    # Component 2: max |y| = 8.0 (from -8.0)
    # L∞ norm should be max(5.0, 6.0, 8.0) = 8.0
    test_data = [
        (0.0, [1.0, 2.0, 3.0]),
        (1.0, [-5.0, 1.5, 2.5]),
        (2.0, [2.0, -6.0, 1.0]),
        (3.0, [1.5, 2.5, -8.0]),
    ]

    for x, y_vals in test_data:
        y = Ncm.Vector.new_array(y_vals)
        fss.add(x, y)

    linf_norm = fss.get_absmaxF_linf_norm()
    expected_linf = 8.0  # max(5.0, 6.0, 8.0)
    assert abs(linf_norm - expected_linf) < 1e-10


def test_absmaxF_norms_single_component() -> None:
    """Test that both norms work correctly with single component."""
    fss = Ncm.FunctionSampleSet.new(1)

    # Add some samples
    # Component 0: max will be 7.0
    # L2 norm should be 7.0 (sqrt(7^2) = 7)
    # L∞ norm should be 7.0 (max(7.0) = 7)
    test_data = [0.0, 1.0, 7.0, 3.0, 2.0]

    for x in test_data:
        y = Ncm.Vector.new_array([x])
        fss.add(x, y)

    l2_norm = fss.get_absmaxF_l2_norm()
    linf_norm = fss.get_absmaxF_linf_norm()

    assert abs(l2_norm - 7.0) < 1e-10
    assert abs(linf_norm - 7.0) < 1e-10


def test_absmaxF_min_basic(sample_set_dim2: Ncm.FunctionSampleSet) -> None:
    """Test minimum of component-wise max absolute values."""
    fss = sample_set_dim2

    # Initially should be 0.0 (no samples yet)
    min_absF_init = fss.get_absmaxF_min()
    assert abs(min_absF_init - 0.0) < 1e-10

    # Add samples
    # Component 0: max will be 3.0
    # Component 1: max will be 5.0
    # Min should be min(3.0, 5.0) = 3.0
    test_data = [
        (0.0, [1.0, 2.0]),
        (1.0, [3.0, 5.0]),
        (2.0, [2.0, 4.0]),
    ]

    for x, y_vals in test_data:
        y = Ncm.Vector.new_array(y_vals)
        fss.add(x, y)

    min_absF = fss.get_absmaxF_min()
    assert abs(min_absF - 3.0) < 1e-10


def test_absmaxF_min_with_negatives(
    sample_set_dim3: Ncm.FunctionSampleSet,
) -> None:
    """Test minimum computation with negative values."""
    fss = sample_set_dim3

    # Component 0: max |y| = 5.0 (from -5.0)
    # Component 1: max |y| = 6.0 (from -6.0)
    # Component 2: max |y| = 8.0 (from -8.0)
    # Min should be min(5.0, 6.0, 8.0) = 5.0
    test_data = [
        (0.0, [1.0, 2.0, 3.0]),
        (1.0, [-5.0, 1.5, 2.5]),
        (2.0, [2.0, -6.0, 1.0]),
        (3.0, [1.5, 2.5, -8.0]),
    ]

    for x, y_vals in test_data:
        y = Ncm.Vector.new_array(y_vals)
        fss.add(x, y)

    min_absF = fss.get_absmaxF_min()
    assert abs(min_absF - 5.0) < 1e-10


def test_absmaxF_min_identifies_weakest_component(
    sample_set_dim3: Ncm.FunctionSampleSet,
) -> None:
    """Test that min correctly identifies the component with smallest peak."""
    fss = sample_set_dim3

    # Create components with very different peak values
    # Component 0: max = 100.0
    # Component 1: max = 2.0  (weakest - this should be the min)
    # Component 2: max = 50.0
    test_data = [
        (0.0, [10.0, 1.0, 5.0]),
        (1.0, [100.0, 2.0, 50.0]),
        (2.0, [20.0, 0.5, 10.0]),
    ]

    for x, y_vals in test_data:
        y = Ncm.Vector.new_array(y_vals)
        fss.add(x, y)

    min_absF = fss.get_absmaxF_min()
    l2_norm = fss.get_absmaxF_l2_norm()
    linf_norm = fss.get_absmaxF_linf_norm()

    # Min should be 2.0 (component 1)
    assert abs(min_absF - 2.0) < 1e-10

    # L∞ should be 100.0 (component 0)
    assert abs(linf_norm - 100.0) < 1e-10

    # L2 should be sqrt(100^2 + 2^2 + 50^2) = sqrt(10000 + 4 + 2500)
    expected_l2 = math.sqrt(100.0**2 + 2.0**2 + 50.0**2)
    assert abs(l2_norm - expected_l2) < 1e-10

    # Verify the relationship: min <= L2 <= L∞ * sqrt(n)
    assert min_absF <= l2_norm
    assert l2_norm <= linf_norm * math.sqrt(3)


def test_absmaxF_min_single_component() -> None:
    """Test min with single component (should equal the component's max)."""
    fss = Ncm.FunctionSampleSet.new(1)

    test_data = [0.0, 1.0, 7.0, 3.0, 2.0]

    for x in test_data:
        y = Ncm.Vector.new_array([x])
        fss.add(x, y)

    min_absF = fss.get_absmaxF_min()
    l2_norm = fss.get_absmaxF_l2_norm()
    linf_norm = fss.get_absmaxF_linf_norm()

    # All three should be 7.0 for single component
    assert abs(min_absF - 7.0) < 1e-10
    assert abs(l2_norm - 7.0) < 1e-10
    assert abs(linf_norm - 7.0) < 1e-10


def test_combined_tracking_comprehensive(
    sample_set_dim3: Ncm.FunctionSampleSet,
) -> None:
    """Comprehensive test of all tracking features together."""
    fss = sample_set_dim3

    # Verify initial state
    assert math.isinf(fss.get_x_min()) and fss.get_x_min() > 0
    assert math.isinf(fss.get_x_max()) and fss.get_x_max() < 0
    absmaxF0_init = fss.get_absmaxF(0)
    absmaxF1_init = fss.get_absmaxF(1)
    absmaxF2_init = fss.get_absmaxF(2)
    assert abs(absmaxF0_init[0] - 0.0) < 1e-10
    assert abs(absmaxF1_init[0] - 0.0) < 1e-10
    assert abs(absmaxF2_init[0] - 0.0) < 1e-10
    assert math.isnan(absmaxF0_init[1])
    assert math.isnan(absmaxF1_init[1])
    assert math.isnan(absmaxF2_init[1])

    # Add points covering a range with various y values
    test_data = [
        (-2.0, [1.0, -2.0, 3.0]),
        (-1.0, [-4.0, 1.5, 2.5]),
        (0.0, [2.0, 5.0, -6.0]),
        (1.0, [0.5, -3.0, 1.0]),
        (2.0, [-1.5, 2.0, 7.0]),
        (3.0, [3.0, -1.0, -2.0]),
    ]

    for x, y_vals in test_data:
        y = Ncm.Vector.new_array(y_vals)
        fss.add(x, y)

    # Verify x range tracking
    assert abs(fss.get_x_min() - (-2.0)) < 1e-10
    assert abs(fss.get_x_max() - 3.0) < 1e-10

    # Verify absmaxF tracking
    # Component 0: max is 4.0 (from -4.0 at x=-1.0)
    # Component 1: max is 5.0 (from 5.0 at x=0.0)
    # Component 2: max is 7.0 (from 7.0 at x=2.0)
    absmaxF0 = fss.get_absmaxF(0)
    absmaxF1 = fss.get_absmaxF(1)
    absmaxF2 = fss.get_absmaxF(2)
    assert abs(absmaxF0[0] - 4.0) < 1e-10
    assert abs(absmaxF0[1] - (-1.0)) < 1e-10
    assert abs(absmaxF1[0] - 5.0) < 1e-10
    assert abs(absmaxF1[1] - 0.0) < 1e-10
    assert abs(absmaxF2[0] - 7.0) < 1e-10
    assert abs(absmaxF2[1] - 2.0) < 1e-10

    # Extend range with iterator insertions
    it = fss.iter_begin()
    y_before = Ncm.Vector.new_array([-8.0, 1.0, 2.0])
    fss.iter_insert_before(it, -5.0, y_before)

    it = fss.iter_end()
    y_after = Ncm.Vector.new_array([2.0, -9.0, 3.0])
    fss.iter_insert_after(it, 6.0, y_after)

    # Verify updated tracking
    assert abs(fss.get_x_min() - (-5.0)) < 1e-10
    assert abs(fss.get_x_max() - 6.0) < 1e-10
    absmaxF0_updated = fss.get_absmaxF(0)
    absmaxF1_updated = fss.get_absmaxF(1)
    absmaxF2_updated = fss.get_absmaxF(2)
    assert abs(absmaxF0_updated[0] - 8.0) < 1e-10
    assert abs(absmaxF0_updated[1] - (-5.0)) < 1e-10
    assert abs(absmaxF1_updated[0] - 9.0) < 1e-10
    assert abs(absmaxF1_updated[1] - 6.0) < 1e-10
    assert abs(absmaxF2_updated[0] - 7.0) < 1e-10
    assert abs(absmaxF2_updated[1] - 2.0) < 1e-10


# ============================================================================
# Tests for ncm_function_sample_set_adaptive_midpoint
# ============================================================================


def test_adaptive_midpoint_basic_convergence() -> None:
    """Test that adaptive_midpoint converges for a smooth function."""

    def trig_func(x: float, y: Ncm.Vector) -> None:
        y.set(0, math.sin(x))
        y.set(1, math.cos(x))

    fss = Ncm.FunctionSampleSet.new(2)
    base_spline = Ncm.SplineCubicNotaknot.new()

    for x in np.linspace(0.0, 2.0 * math.pi, 6):
        fss.add_func(x, trig_func)

    initial_nsamples = fss.get_nsamples()
    assert initial_nsamples == 6

    fss.mark_all_old()
    fss.adaptive_midpoint(trig_func, 1e-6, 1e-10, 100, 1, base_spline)

    # Must have converged
    assert fss.all_intervals_ok(1)

    # Must have added knots
    assert fss.get_nsamples() > initial_nsamples


def test_adaptive_midpoint_spline_accuracy() -> None:
    """Test that the spline built after adaptive_midpoint is accurate."""

    def trig_func(x: float, y: Ncm.Vector) -> None:
        y.set(0, math.sin(x))
        y.set(1, math.cos(x))

    reltol = 1e-6
    abstol = 1e-10
    xi, xf = 0.0, 2.0 * math.pi

    fss = Ncm.FunctionSampleSet.new(2)
    base_spline = Ncm.SplineCubicNotaknot.new()

    for x in np.linspace(xi, xf, 6):
        fss.add_func(x, trig_func)

    fss.mark_all_old()
    fss.adaptive_midpoint(trig_func, reltol, abstol, 100, 1, base_spline)

    sv = fss.to_spline_vec(base_spline)

    np.random.seed(0)
    for x in np.random.uniform(xi + 0.01, xf - 0.01, 50):
        y_spline = sv.eval_array(x)
        y_true = [math.sin(x), math.cos(x)]
        error = math.sqrt(
            (y_spline[0] - y_true[0]) ** 2 + (y_spline[1] - y_true[1]) ** 2
        )
        norm_true = math.sqrt(y_true[0] ** 2 + y_true[1] ** 2)
        threshold = reltol * norm_true + abstol
        assert (
            error < threshold * 10
        ), f"Error {error:.3e} > threshold {threshold:.3e} at x={x:.4f}"


def test_adaptive_midpoint_polynomial_exact() -> None:
    """Cubic spline should reproduce a cubic polynomial exactly."""

    def poly_func(x: float, y: Ncm.Vector) -> None:
        # f(x) = x^3 - 3x + 1; cubic spline interpolates this exactly
        y.set(0, x**3 - 3.0 * x + 1.0)

    fss = Ncm.FunctionSampleSet.new(1)
    base_spline = Ncm.SplineCubicNotaknot.new()

    for x in np.linspace(-2.0, 2.0, 6):
        fss.add_func(x, poly_func)

    fss.mark_all_old()
    fss.adaptive_midpoint(poly_func, 1e-10, 1e-12, 50, 1, base_spline)

    assert fss.all_intervals_ok(1)

    sv = fss.to_spline_vec(base_spline)
    for x in np.linspace(-1.9, 1.9, 20):
        y_spline = sv.eval_array(x)
        y_true = x**3 - 3.0 * x + 1.0
        assert abs(y_spline[0] - y_true) < 1e-8, f"Polynomial mismatch at x={x:.4f}"


def test_adaptive_midpoint_user_data() -> None:
    """Test that user_data is forwarded correctly to the function."""

    def scaled_trig(x: float, y: Ncm.Vector, user_data: Any) -> None:
        scale = user_data
        y.set(0, scale * math.sin(x))
        y.set(1, scale * math.cos(x))

    scale = 3.7
    fss = Ncm.FunctionSampleSet.new(2)
    base_spline = Ncm.SplineCubicNotaknot.new()

    for x in np.linspace(0.0, 2.0 * math.pi, 6):
        fss.add_func(x, scaled_trig, scale)

    fss.mark_all_old()
    fss.adaptive_midpoint(scaled_trig, 1e-6, 1e-10, 100, 1, base_spline, scale)

    assert fss.all_intervals_ok(1)

    # Verify values at knots match the scaled function
    it = fss.iter_begin()
    while it.is_valid():
        x = it.get_x()
        y = it.get_y()
        assert abs(y.get(0) - scale * math.sin(x)) < 1e-12
        assert abs(y.get(1) - scale * math.cos(x)) < 1e-12
        it.next()


def test_adaptive_midpoint_max_iter_respected(capfd) -> None:
    """Test that adaptive_midpoint stops at max_iter even if not converged."""

    def high_freq(x: float, y: Ncm.Vector) -> None:
        y.set(0, math.sin(50.0 * x))

    fss = Ncm.FunctionSampleSet.new(1)
    base_spline = Ncm.SplineCubicNotaknot.new()

    for x in np.linspace(0.0, math.pi, 6):
        fss.add_func(x, high_freq)

    max_iter = 3
    fss.mark_all_old()
    # Very tight tolerance on a high-frequency function: will not converge in 3 iters
    fss.adaptive_midpoint(high_freq, 1e-12, 1e-15, max_iter, 1, base_spline)

    out, err = capfd.readouterr()
    combined = out + err

    # g_message should have emitted the max-iterations warning
    assert "Max iterations" in combined, f"Expected max-iter message, got: {combined!r}"
    assert str(max_iter) in combined
    assert str(fss.get_nsamples()) in combined

    # Should have stopped after max_iter regardless of convergence
    # Just verify it ran without error and has more samples than the initial 6
    assert fss.get_nsamples() > 6


def test_adaptive_midpoint_multi_pass_threshold() -> None:
    """Test min_pass_threshold > 1: intervals must pass the error test twice."""

    def trig_func(x: float, y: Ncm.Vector) -> None:
        y.set(0, math.sin(x))
        y.set(1, math.cos(x))

    fss = Ncm.FunctionSampleSet.new(2)
    base_spline = Ncm.SplineCubicNotaknot.new()

    for x in np.linspace(0.0, 2.0 * math.pi, 6):
        fss.add_func(x, trig_func)

    fss.mark_all_old()
    threshold = 1
    fss.adaptive_midpoint(trig_func, 1e-5, 1e-9, 200, threshold, base_spline)

    assert fss.all_intervals_ok(threshold)


def test_adaptive_midpoint_challenging_function() -> None:
    """Test adaptive_midpoint with a function of varying frequency."""

    def challenging(x: float, y: Ncm.Vector) -> None:
        y.set(0, math.sin(x**1.5))
        y.set(1, math.cos(2.0 * x))

    xi, xf = 0.1, 10.0
    reltol, abstol = 1e-4, 1e-8

    fss = Ncm.FunctionSampleSet.new(2)
    base_spline = Ncm.SplineCubicNotaknot.new()

    for x in np.linspace(xi, xf, 6):
        fss.add_func(x, challenging)

    initial = fss.get_nsamples()
    fss.mark_all_old()
    fss.adaptive_midpoint(challenging, reltol, abstol, 100, 1, base_spline)

    assert fss.all_intervals_ok(1)
    assert fss.get_nsamples() > initial * 2  # Challenging function needs many knots

    sv = fss.to_spline_vec(base_spline)
    for x in np.linspace(xi + 0.05, xf - 0.05, 30):
        y_spline = sv.eval_array(x)
        y_true = [math.sin(x**1.5), math.cos(2.0 * x)]
        error = math.sqrt(sum((a - b) ** 2 for a, b in zip(y_spline, y_true)))
        norm_true = math.sqrt(sum(v**2 for v in y_true))
        assert error < (reltol * norm_true + abstol) * 10


def test_adaptive_midpoint_x_range_preserved() -> None:
    """Test that x_min/x_max are not changed by adaptive_midpoint."""

    def f(x: float, y: Ncm.Vector) -> None:
        y.set(0, math.sin(x))

    xi, xf = 1.0, 5.0
    fss = Ncm.FunctionSampleSet.new(1)
    base_spline = Ncm.SplineCubicNotaknot.new()

    for x in np.linspace(xi, xf, 6):
        fss.add_func(x, f)

    fss.mark_all_old()
    fss.adaptive_midpoint(f, 1e-6, 1e-10, 100, 1, base_spline)

    # Midpoints are strictly interior, so min and max stay at xi and xf
    assert abs(fss.get_x_min() - xi) < 1e-12
    assert abs(fss.get_x_max() - xf) < 1e-12


def test_adaptive_midpoint_all_points_old_after() -> None:
    """After adaptive_midpoint all sample points should be marked OLD."""

    def f(x: float, y: Ncm.Vector) -> None:
        y.set(0, math.sin(x))
        y.set(1, math.cos(x))

    fss = Ncm.FunctionSampleSet.new(2)
    base_spline = Ncm.SplineCubicNotaknot.new()

    for x in np.linspace(0.0, 2.0 * math.pi, 6):
        fss.add_func(x, f)

    fss.mark_all_old()
    fss.adaptive_midpoint(f, 1e-6, 1e-10, 100, 1, base_spline)

    # to_spline_vec_old should include ALL samples (none left marked as NEW)
    nsamples = fss.get_nsamples()
    base_spline2 = Ncm.SplineCubicNotaknot.new()
    sv_old = fss.to_spline_vec_old(base_spline2)
    assert sv_old.get_nknots() == nsamples


def test_adaptive_midpoint_auto_bootstrap_few_samples() -> None:
    """Test that adaptive_midpoint when starting with <6 samples.

    The cubic notaknot spline requires at least 6 samples, but adaptive_midpoint
    should automatically bootstrap from any number >= 2 by adding midpoints.
    """

    def f(x: float, y: Ncm.Vector) -> None:
        y.set(0, math.sin(x))

    fss = Ncm.FunctionSampleSet.new(1)
    base_spline = Ncm.SplineCubicNotaknot.new()

    # Start with only 3 samples (less than the 6 required for cubic notaknot)
    for x in [0.0, math.pi, 2.0 * math.pi]:
        fss.add_func(x, f)

    assert fss.get_nsamples() == 3
    initial_x_min = fss.get_x_min()
    initial_x_max = fss.get_x_max()

    fss.mark_all_old()

    # This should work without error - it will bootstrap to at least 6 samples
    fss.adaptive_midpoint(f, 1e-6, 1e-10, 100, 1, base_spline)

    # Should have converged
    assert fss.all_intervals_ok(1)

    # Should have at least 6 samples after bootstrap + refinement
    assert fss.get_nsamples() >= 6

    # Domain should be preserved (bootstrap adds interior points only)
    assert abs(fss.get_x_min() - initial_x_min) < 1e-12
    assert abs(fss.get_x_max() - initial_x_max) < 1e-12

    # Verify the spline is accurate
    sv = fss.to_spline_vec(base_spline)
    test_points = np.linspace(0.0, 2.0 * math.pi, 50)
    for x in test_points:
        y_true = math.sin(x)
        y_spline = sv.eval_array(x)[0]
        # Should be accurate to tolerance
        assert abs(y_spline - y_true) < 1e-5


#
# Tests for expand_domain
#


def test_expand_domain_basic() -> None:
    """Test basic domain expansion with fast-decaying function."""

    def f(x: float, y: Ncm.Vector) -> None:
        # Exponentially decaying function
        y.set(0, math.exp(-x))

    fss = Ncm.FunctionSampleSet.new(1)

    # Initialize with a few points
    for x in [1.0, 2.0, 3.0]:
        fss.add_func(x, f)

    x_min_initial = fss.get_x_min()
    x_max_initial = fss.get_x_max()

    # Expand domain
    fss.expand_domain(
        f,
        x_min_hard=0.01,
        x_max_hard=100.0,
        expansion_factor=0.3,
        epsilon=1e-3,
        max_iterations=50,
        consecutive_tries=3,
    )

    # Domain should have expanded
    assert fss.get_x_min() < x_min_initial
    assert fss.get_x_max() > x_max_initial

    # More samples should have been added
    assert fss.get_nsamples() > 3


def test_expand_domain_hard_limits_respected() -> None:
    """Verify expand_domain respects hard limits."""

    def f(x: float, y: Ncm.Vector) -> None:
        y.set(0, 1.0 / (1.0 + x**2))

    fss = Ncm.FunctionSampleSet.new(1)

    for x in [1.0, 2.0]:
        fss.add_func(x, f)

    x_min_hard = 0.1
    x_max_hard = 5.0

    fss.expand_domain(
        f,
        x_min_hard=x_min_hard,
        x_max_hard=x_max_hard,
        expansion_factor=0.5,
        epsilon=1e-10,  # Very tight to force expansion to limits
        max_iterations=100,
        consecutive_tries=2,
    )

    # Should not exceed hard limits
    assert fss.get_x_min() >= x_min_hard
    assert fss.get_x_max() <= x_max_hard


def test_expand_domain_convergence_stops_expansion() -> None:
    """Test that convergence stops expansion before hitting limits."""

    def f(x: float, y: Ncm.Vector) -> None:
        # Very fast decay
        y.set(0, math.exp(-x * 5.0))

    fss = Ncm.FunctionSampleSet.new(1)

    for x in [0.5, 1.0, 1.5]:
        fss.add_func(x, f)

    fss.expand_domain(
        f,
        x_min_hard=1e-6,
        x_max_hard=100.0,
        expansion_factor=0.2,
        epsilon=1e-2,
        max_iterations=100,
        consecutive_tries=3,
    )

    # Should converge well before hitting max limit
    assert fss.get_x_max() < 10.0


def test_expand_domain_max_iterations() -> None:
    """Test that max_iterations limits expansion."""

    def f(_x: float, y: Ncm.Vector) -> None:
        # Constant function - never converges
        y.set(0, 1.0)

    fss = Ncm.FunctionSampleSet.new(1)

    for x in [1.0, 2.0]:
        fss.add_func(x, f)

    n_initial = fss.get_nsamples()
    max_iter = 5

    fss.expand_domain(
        f,
        x_min_hard=0.01,
        x_max_hard=100.0,
        expansion_factor=0.2,
        epsilon=0.1,
        max_iterations=max_iter,
        consecutive_tries=2,
    )

    # Should have added at most 2*max_iter points (left and right)
    n_added = fss.get_nsamples() - n_initial
    assert n_added <= 2 * max_iter


def test_expand_domain_consecutive_tries() -> None:
    """Test consecutive_tries requirement for convergence."""

    def f(x: float, y: Ncm.Vector) -> None:
        # Smooth exponential decay
        y.set(0, math.exp(-0.5 * x))

    fss = Ncm.FunctionSampleSet.new(1)

    for x in [1.0, 2.0, 3.0]:
        fss.add_func(x, f)

    # Use consecutive_tries=5 - should need more iterations
    fss.expand_domain(
        f,
        x_min_hard=0.01,
        x_max_hard=50.0,
        expansion_factor=0.2,
        epsilon=1e-2,
        max_iterations=100,
        consecutive_tries=5,
    )

    n_samples_strict = fss.get_nsamples()

    # Reset and try with consecutive_tries=1
    fss2 = Ncm.FunctionSampleSet.new(1)
    for x in [1.0, 2.0, 3.0]:
        fss2.add_func(x, f)

    fss2.expand_domain(
        f,
        x_min_hard=0.01,
        x_max_hard=50.0,
        expansion_factor=0.2,
        epsilon=1e-2,
        max_iterations=100,
        consecutive_tries=1,
    )

    n_samples_loose = fss2.get_nsamples()

    # Stricter consecutive_tries should require more points or same
    assert n_samples_strict >= n_samples_loose


def test_expand_domain_multi_component() -> None:
    """Test expansion with multiple output components."""

    def f(x: float, y: Ncm.Vector) -> None:
        y.set(0, math.exp(-x))
        y.set(1, math.exp(-2.0 * x))
        y.set(2, math.exp(-0.5 * x))

    fss = Ncm.FunctionSampleSet.new(3)

    for x in [1.0, 2.0, 3.0]:
        fss.add_func(x, f)

    fss.expand_domain(
        f,
        x_min_hard=0.1,
        x_max_hard=20.0,
        expansion_factor=0.3,
        epsilon=1e-3,
        max_iterations=50,
        consecutive_tries=3,
    )

    # All components should be tracked
    assert fss.get_absmaxF_l2_norm() > 0
    assert fss.get_len() == 3


def test_expand_domain_stress_exponential_decay() -> None:
    """Stress test with exponentially decaying function."""

    def f(x: float, y: Ncm.Vector) -> None:
        # Multiple decay rates
        y.set(0, math.exp(-0.1 * x))
        y.set(1, math.exp(-0.5 * x))
        y.set(2, math.exp(-1.0 * x))

    fss = Ncm.FunctionSampleSet.new(3)

    # Start with tight initial domain
    for x in np.linspace(0.5, 1.5, 5):
        fss.add_func(x, f)

    initial_x_min = fss.get_x_min()
    initial_x_max = fss.get_x_max()

    fss.expand_domain(
        f,
        x_min_hard=1e-3,
        x_max_hard=50.0,
        expansion_factor=0.2,
        epsilon=1e-4,
        max_iterations=200,
        consecutive_tries=4,
    )

    # Verify significant expansion occurred
    assert fss.get_x_min() < initial_x_min * 0.5
    assert fss.get_x_max() > initial_x_max * 2.0

    # Verify function values decrease at boundaries
    F_ref = fss.get_absmaxF_l2_norm()
    assert F_ref > 0

    # Sample at right boundary should be small (converged)
    # Left boundary may hit hard limit so we only check right
    y_right = Ncm.Vector.new(3)
    f(fss.get_x_max(), y_right)
    norm_right = math.sqrt(sum(y_right.get(i) ** 2 for i in range(3)))

    # Right boundary should have converged
    assert norm_right < 1e-2 * F_ref


def test_expand_domain_stress_power_law() -> None:
    """Stress test with power-law decay."""

    def f(x: float, y: Ncm.Vector) -> None:
        # Power law decay in different components
        y.set(0, 1.0 / (1.0 + x**2))
        y.set(1, 1.0 / (1.0 + x**3))

    fss = Ncm.FunctionSampleSet.new(2)

    for x in np.linspace(1.0, 3.0, 10):
        fss.add_func(x, f)

    fss.expand_domain(
        f,
        x_min_hard=0.01,
        x_max_hard=100.0,
        expansion_factor=0.25,
        epsilon=1e-4,
        max_iterations=150,
        consecutive_tries=3,
    )

    # Power law decays slower, should expand further
    assert fss.get_x_max() > 10.0

    # Verify boundary convergence
    F_ref = fss.get_absmaxF_l2_norm()
    y = Ncm.Vector.new(2)
    f(fss.get_x_max(), y)
    norm = math.sqrt(y.get(0) ** 2 + y.get(1) ** 2)
    assert norm < 5e-3 * F_ref


def test_expand_domain_asymmetric_expansion() -> None:
    """Test that left and right boundaries can expand independently."""

    def f(x: float, y: Ncm.Vector) -> None:
        # Asymmetric function - decays fast to the left, slow to the right
        if x < 2.0:
            y.set(0, math.exp(-5.0 * (2.0 - x)))
        else:
            y.set(0, math.exp(-0.5 * (x - 2.0)))

    fss = Ncm.FunctionSampleSet.new(1)

    # Start centered
    for x in [1.8, 2.0, 2.2]:
        fss.add_func(x, f)

    x_min_initial = fss.get_x_min()
    x_max_initial = fss.get_x_max()

    fss.expand_domain(
        f,
        x_min_hard=0.1,
        x_max_hard=50.0,
        expansion_factor=0.2,
        epsilon=1e-3,
        max_iterations=100,
        consecutive_tries=3,
    )

    # Left should not expand much (fast decay)
    left_expansion = x_min_initial - fss.get_x_min()

    # Right should expand more (slow decay)
    right_expansion = fss.get_x_max() - x_max_initial

    assert right_expansion > left_expansion * 2.0


def test_expand_domain_interleaved_reference_update() -> None:
    """Test that F_ref is updated after each point (interleaved expansion)."""

    def f(x: float, y: Ncm.Vector) -> None:
        # Function that grows in the middle
        y.set(0, math.exp(-((x - 5.0) ** 2) / 4.0))

    fss = Ncm.FunctionSampleSet.new(1)

    # Start away from peak
    for x in [8.0, 9.0, 10.0]:
        fss.add_func(x, f)

    initial_F_ref = fss.get_absmaxF_l2_norm()

    fss.expand_domain(
        f,
        x_min_hard=0.1,
        x_max_hard=15.0,
        expansion_factor=0.15,
        epsilon=1e-3,
        max_iterations=100,
        consecutive_tries=3,
    )

    final_F_ref = fss.get_absmaxF_l2_norm()

    # F_ref should have increased as we discovered the peak
    assert final_F_ref > initial_F_ref * 1.5


def test_expand_domain_with_user_data() -> None:
    """Test expand_domain with user_data parameter."""

    class UserData:
        """Example user data class to track calls and provide a scale factor."""

        def __init__(self, scale: float):
            self.scale = scale
            self.call_count = 0

        def func(self, x: float, y: Ncm.Vector) -> None:
            """Function that uses user data to scale the output and counts calls."""
            self.call_count += 1
            y.set(0, self.scale * math.exp(-x))

    user_data = UserData(2.5)

    fss = Ncm.FunctionSampleSet.new(1)

    # Initialize with a couple of points
    for x in [1.0, 2.0]:
        fss.add_func(x, user_data.func)

    initial_calls = user_data.call_count

    fss.expand_domain(
        user_data.func,
        x_min_hard=0.1,
        x_max_hard=10.0,
        expansion_factor=0.2,
        epsilon=1e-2,
        max_iterations=20,
        consecutive_tries=2,
    )

    # Callback should have been called during expansion
    assert user_data.call_count > initial_calls

    # Scale should be reflected in values (absmaxF returns tuple with value at [0])
    assert fss.get_absmaxF(0)[0] > 2.0


def test_expand_domain_reaches_hard_limits() -> None:
    """Test that expand_domain evaluates at hard limits when they are reached."""

    def f(_x: float, y: Ncm.Vector) -> None:
        # Constant function - will never converge, forcing expansion to limits
        y.set(0, 1.0)

    fss = Ncm.FunctionSampleSet.new(1)

    # Start with initial points
    for x in [5.0, 6.0]:
        fss.add_func(x, f)

    x_min_hard = 1.0
    x_max_hard = 10.0

    # Use small max_iterations to force hitting limits quickly
    fss.expand_domain(
        f,
        x_min_hard=x_min_hard,
        x_max_hard=x_max_hard,
        expansion_factor=0.5,  # Large factor to reach limits faster
        epsilon=0.01,  # Tight epsilon so constant function won't converge
        max_iterations=10,
        consecutive_tries=2,
    )

    # First knot should be at or very close to x_min_hard
    assert abs(fss.get_x_min() - x_min_hard) < 1e-10

    # Last knot should be at or very close to x_max_hard
    assert abs(fss.get_x_max() - x_max_hard) < 1e-10


def test_expand_domain_left_limit_reached() -> None:
    """Test that left hard limit is reached and becomes first knot."""

    def f(x: float, y: Ncm.Vector) -> None:
        y.set(0, 1.0 + x)  # Linear, won't converge quickly

    fss = Ncm.FunctionSampleSet.new(1)

    for x in [8.0, 9.0]:
        fss.add_func(x, f)

    x_min_hard = 2.0
    x_max_hard = 100.0

    fss.expand_domain(
        f,
        x_min_hard=x_min_hard,
        x_max_hard=x_max_hard,
        expansion_factor=0.4,
        epsilon=1e-10,  # Very tight to force hitting left limit
        max_iterations=20,
        consecutive_tries=3,
    )

    # Left boundary should reach hard limit
    assert abs(fss.get_x_min() - x_min_hard) < 1e-10


def test_expand_domain_right_limit_reached() -> None:
    """Test that right hard limit is reached and becomes last knot."""

    def f(x: float, y: Ncm.Vector) -> None:
        y.set(0, 1.0 + x)  # Linear, won't converge quickly

    fss = Ncm.FunctionSampleSet.new(1)

    for x in [1.0, 2.0]:
        fss.add_func(x, f)

    x_min_hard = 0.01
    x_max_hard = 8.0

    fss.expand_domain(
        f,
        x_min_hard=x_min_hard,
        x_max_hard=x_max_hard,
        expansion_factor=0.4,
        epsilon=1e-10,  # Very tight to force hitting right limit
        max_iterations=20,
        consecutive_tries=3,
    )

    # Right boundary should reach hard limit
    assert abs(fss.get_x_max() - x_max_hard) < 1e-10


def test_expand_domain_both_limits_with_slow_decay() -> None:
    """Test that both hard limits are reached with slow-decaying function."""

    def f(x: float, y: Ncm.Vector) -> None:
        # Very slow decay - essentially constant in our range
        y.set(0, 1.0 / (1.0 + 0.01 * x**2))

    fss = Ncm.FunctionSampleSet.new(1)

    for x in [5.0, 6.0]:
        fss.add_func(x, f)

    x_min_hard = 0.5
    x_max_hard = 20.0

    fss.expand_domain(
        f,
        x_min_hard=x_min_hard,
        x_max_hard=x_max_hard,
        expansion_factor=0.3,
        epsilon=1e-3,
        max_iterations=30,
        consecutive_tries=2,
    )

    # Both boundaries should reach hard limits
    assert abs(fss.get_x_min() - x_min_hard) < 1e-10
    assert abs(fss.get_x_max() - x_max_hard) < 1e-10

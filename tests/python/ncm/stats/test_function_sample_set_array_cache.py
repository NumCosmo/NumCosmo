#!/usr/bin/env python
#
# test_py_function_sample_set_array_cache.py
#
# Mon Mar 17 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# Tests for array caching behavior in NcmFunctionSampleSet

"""Tests for array cache reuse in NcmFunctionSampleSet.to_spline_vec()."""

import math
import numpy as np

from numcosmo_py import Ncm

Ncm.cfg_init()


def test_array_cache_multiple_points_between_knots() -> None:
    """Test adding multiple points between existing knots works with cached arrays."""
    fss = Ncm.FunctionSampleSet.new(2)
    base_spline = Ncm.SplineCubicNotaknot.new()

    # Add initial points (at least 6 for cubic notaknot)
    for x in [0.0, 5.0, 10.0, 15.0, 20.0, 25.0]:
        y = Ncm.Vector.new(2)
        y.set(0, x**2)
        y.set(1, x**3)
        fss.add(x, y)

    # Convert to spline - first use of cache
    sv1 = fss.to_spline_vec(base_spline)
    assert sv1.get_nknots() == 6

    # Add multiple points between first two knots
    for x in [1.0, 2.0, 3.0, 4.0]:
        y = Ncm.Vector.new(2)
        y.set(0, x**2)
        y.set(1, x**3)
        fss.add(x, y)

    # Convert again - should reuse and resize cache
    sv2 = fss.to_spline_vec(base_spline)
    assert sv2.get_nknots() == 10

    # Verify accuracy at test point
    x_test = 2.5
    y_expected = [x_test**2, x_test**3]
    y_interp = sv2.eval_array(x_test)

    assert abs(y_interp[0] - y_expected[0]) < 1e-10
    assert abs(y_interp[1] - y_expected[1]) < 1e-8


def test_array_cache_add_points_at_end() -> None:
    """Test adding points after the last knot works with cached arrays."""
    fss = Ncm.FunctionSampleSet.new(1)
    base_spline = Ncm.SplineCubicNotaknot.new()

    # Add initial points
    for x in np.linspace(0, 10, 6):
        y = Ncm.Vector.new(1)
        y.set(0, math.sin(x))
        fss.add(x, y)

    # First conversion
    sv1 = fss.to_spline_vec(base_spline)
    assert sv1.get_nknots() == 6

    # Add points at the end
    for x in [11.0, 12.0, 13.0, 14.0, 15.0]:
        y = Ncm.Vector.new(1)
        y.set(0, math.sin(x))
        fss.add(x, y)

    # Second conversion - cache should handle larger size
    sv2 = fss.to_spline_vec(base_spline)
    assert sv2.get_nknots() == 11

    # Verify the new points are included
    x_test = 13.5
    y_expected = math.sin(x_test)
    y_interp_array = sv2.eval_array(x_test)
    y_interp = y_interp_array[0]
    # Use slightly more lenient tolerance for interpolation
    assert abs(y_interp - y_expected) < 0.02


def test_array_cache_add_points_at_beginning() -> None:
    """Test adding points before the first knot works with cached arrays."""
    fss = Ncm.FunctionSampleSet.new(2)
    base_spline = Ncm.SplineCubicNotaknot.new()

    # Add initial points starting at x=10 (at least 6 for cubic notaknot)
    for x in [10.0, 12.0, 14.0, 16.0, 18.0, 20.0]:
        y = Ncm.Vector.new(2)
        y.set(0, x)
        y.set(1, x**2)
        fss.add(x, y)

    # First conversion
    sv1 = fss.to_spline_vec(base_spline)
    assert sv1.get_nknots() == 6

    # Add points at the beginning
    for x in [0.0, 2.5, 5.0, 7.5]:
        y = Ncm.Vector.new(2)
        y.set(0, x)
        y.set(1, x**2)
        fss.add(x, y)

    # Second conversion
    sv2 = fss.to_spline_vec(base_spline)
    assert sv2.get_nknots() == 10

    # Verify accuracy at a point in the new range (beginning)
    x_test = 3.75
    y_expected = [x_test, x_test**2]
    y_interp = sv2.eval_array(x_test)

    assert abs(y_interp[0] - y_expected[0]) < 1e-10
    assert abs(y_interp[1] - y_expected[1]) < 1e-9


def test_array_cache_extend_range_both_ends() -> None:
    """Test extending interpolation range on both ends."""
    fss = Ncm.FunctionSampleSet.new(2)
    base_spline = Ncm.SplineCubicNotaknot.new()

    def func(x: float) -> tuple[float, float]:
        return (math.sin(x), math.cos(x))

    # Initial range: [5, 15]
    for x in np.linspace(5.0, 15.0, 6):
        y = Ncm.Vector.new(2)
        y.set(0, math.sin(x))
        y.set(1, math.cos(x))
        fss.add(x, y)

    sv1 = fss.to_spline_vec(base_spline)
    assert sv1.get_nknots() == 6

    # Extend left by 10% (0.5 units)
    # Extend right by 10% (0.5 units)
    left_extension = 0.5
    right_extension = 0.5

    # Add points on the left
    for x in np.linspace(5.0 - left_extension, 5.0, 3, endpoint=False):
        y = Ncm.Vector.new(2)
        vals = func(x)
        y.set(0, vals[0])
        y.set(1, vals[1])
        fss.add(x, y)

    # Add points on the right
    for x in np.linspace(15.0, 15.0 + right_extension, 3)[1:]:
        y = Ncm.Vector.new(2)
        vals = func(x)
        y.set(0, vals[0])
        y.set(1, vals[1])
        fss.add(x, y)

    sv2 = fss.to_spline_vec(base_spline)
    assert sv2.get_nknots() == 11  # 6 original + 3 left + 2 right

    # Verify accuracy in extended ranges (just check they can be evaluated)
    x_test_left = 4.75
    x_test_right = 15.25

    y_left = sv2.eval_array(x_test_left)
    y_left_expected = func(x_test_left)
    assert abs(y_left[0] - y_left_expected[0]) < 1e-2
    assert abs(y_left[1] - y_left_expected[1]) < 1e-2

    y_right = sv2.eval_array(x_test_right)
    y_right_expected = func(x_test_right)
    assert abs(y_right[0] - y_right_expected[0]) < 1e-2
    assert abs(y_right[1] - y_right_expected[1]) < 1e-2


def test_array_cache_multiple_sequential_conversions() -> None:
    """Test multiple sequential conversions with growing sample set."""
    fss = Ncm.FunctionSampleSet.new(3)
    base_spline = Ncm.SplineCubicNotaknot.new()

    # Start with minimal points (at least 6 for cubic notaknot)
    for x in [0.0, 2.0, 4.0, 6.0, 8.0, 10.0]:
        y = Ncm.Vector.new(3)
        y.set(0, x)
        y.set(1, x**2)
        y.set(2, math.sin(x))
        fss.add(x, y)

    # Sequential conversions with progressive additions
    expected_knots = [6, 8, 11, 15, 20]

    for target_knots in expected_knots:
        # Add points to reach target
        current_knots = fss.get_nsamples()
        points_to_add = target_knots - current_knots

        for _ in range(points_to_add):
            # Add a point at a random location within range
            x = np.random.uniform(0.0, 10.0)
            y = Ncm.Vector.new(3)
            y.set(0, x)
            y.set(1, x**2)
            y.set(2, math.sin(x))
            fss.add(x, y)

        # Convert - this reuses the cache each time
        sv = fss.to_spline_vec(base_spline)
        assert sv.get_nknots() == target_knots

        # Verify the spline is valid
        x_test = 5.0
        y_test = sv.eval_array(x_test)
        assert len(y_test) == 3
        # Basic sanity checks
        assert abs(y_test[0] - x_test) < 1.0  # Linear component
        assert abs(y_test[1] - x_test**2) < 10.0  # Quadratic component


def test_array_cache_add_one_point_repeatedly() -> None:
    """Test adding one point at a time and converting repeatedly."""
    fss = Ncm.FunctionSampleSet.new(1)
    base_spline = Ncm.SplineCubicNotaknot.new()

    # Start with minimum for cubic spline
    for x in [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]:
        y = Ncm.Vector.new(1)
        y.set(0, x**2)
        fss.add(x, y)

    # Add one point at a time and convert
    for x in [0.5, 1.5, 2.5, 3.5, 4.5]:
        sv_before = fss.to_spline_vec(base_spline)
        knots_before = sv_before.get_nknots()

        # Add one point
        y = Ncm.Vector.new(1)
        y.set(0, x**2)
        fss.add(x, y)

        # Convert again
        sv_after = fss.to_spline_vec(base_spline)
        knots_after = sv_after.get_nknots()

        # Should have exactly one more knot
        assert knots_after == knots_before + 1

        # Verify accuracy at the new point
        y_interp_array = sv_after.eval_array(x)
        y_interp = y_interp_array[0]
        y_expected = x**2
        assert abs(y_interp - y_expected) < 1e-10


def test_array_cache_with_iterator_insertions() -> None:
    """Test cache works correctly with iterator-based insertions."""
    fss = Ncm.FunctionSampleSet.new(2)
    base_spline = Ncm.SplineCubicNotaknot.new()

    # Add initial sparse grid (at least 6 for cubic notaknot)
    for x in [0.0, 2.0, 4.0, 6.0, 8.0, 10.0]:
        y = Ncm.Vector.new(2)
        y.set(0, math.sin(x))
        y.set(1, math.cos(x))
        fss.add(x, y)

    sv1 = fss.to_spline_vec(base_spline)
    assert sv1.get_nknots() == 6

    # Use iterator to insert a few specific points
    # Insert after x=0.0
    it = fss.iter_begin()
    y_mid = Ncm.Vector.new(2)
    y_mid.set(0, math.sin(1.0))
    y_mid.set(1, math.cos(1.0))
    _ = fss.iter_insert_after(it, 1.0, y_mid)

    # Insert after x=4.0
    it = fss.iter_begin()
    while it.is_valid() and abs(it.get_x() - 4.0) > 1e-10:
        it.next()
    if it.is_valid():
        y_mid = Ncm.Vector.new(2)
        y_mid.set(0, math.sin(5.0))
        y_mid.set(1, math.cos(5.0))
        _ = fss.iter_insert_after(it, 5.0, y_mid)

    sv2 = fss.to_spline_vec(base_spline)
    assert sv2.get_nknots() == 8  # 6 original + 2 inserted

    # Verify accuracy
    x_test = 3.0
    y_expected = [math.sin(x_test), math.cos(x_test)]
    y_interp = sv2.eval_array(x_test)

    # Tolerance is lenient because only 8 knots over [0,10] range
    assert abs(y_interp[0] - y_expected[0]) < 0.1
    assert abs(y_interp[1] - y_expected[1]) < 0.1


def test_array_cache_range_extension_percentage() -> None:
    """Test extending range by specific percentages on both ends."""
    fss = Ncm.FunctionSampleSet.new(1)
    base_spline = Ncm.SplineCubicNotaknot.new()

    # Initial range: [10, 100]
    xi_initial = 10.0
    xf_initial = 100.0
    initial_range = xf_initial - xi_initial

    for x in np.linspace(xi_initial, xf_initial, 10):
        y = Ncm.Vector.new(1)
        y.set(0, math.log(x))
        fss.add(x, y)

    sv1 = fss.to_spline_vec(base_spline)
    assert sv1.get_nknots() == 10

    # Extend by 5% on each end
    extend_pct = 0.05
    left_extend = xi_initial - extend_pct * initial_range
    right_extend = xf_initial + extend_pct * initial_range

    # Add points in extended ranges
    for x in np.linspace(left_extend, xi_initial, 3, endpoint=False):
        y = Ncm.Vector.new(1)
        y.set(0, math.log(x))
        fss.add(x, y)

    for x in np.linspace(xf_initial, right_extend, 3)[1:]:
        y = Ncm.Vector.new(1)
        y.set(0, math.log(x))
        fss.add(x, y)

    sv2 = fss.to_spline_vec(base_spline)
    assert sv2.get_nknots() == 15  # 10 original + 3 left + 2 right

    # Verify range extension by checking evaluation works
    x_test_left = left_extend + 0.5
    x_test_right = right_extend - 0.5

    # Just verify evaluation works in extended ranges
    y_left_array = sv2.eval_array(x_test_left)
    y_left = y_left_array[0]
    y_left_expected = math.log(x_test_left)
    # Tolerance is lenient because spline extrapolation may be less accurate
    assert abs(y_left - y_left_expected) < 0.1

    y_right_array = sv2.eval_array(x_test_right)
    y_right = y_right_array[0]
    y_right_expected = math.log(x_test_right)
    assert abs(y_right - y_right_expected) < 0.1

#!/usr/bin/env python
#
# test_py_function_sample_set_auto_range.py
#
# Mon Mar 17 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# Tests for automatic range determination with localized functions

"""Tests for automatic boundary and knot determination for localized functions."""

import math
import numpy as np

from numcosmo_py import Ncm

Ncm.cfg_init()


def gaussian(x: float, center: float, std: float) -> float:
    """Evaluate Gaussian function."""
    return math.exp(-0.5 * ((x - center) / std) ** 2) / (std * math.sqrt(2.0 * math.pi))


def test_auto_range_gaussian_boundary_extension() -> None:
    """Test automatic boundary extension for localized Gaussian function.

    This test demonstrates how to use FunctionSampleSet for automatic knot
    determination when dealing with a localized function where both the
    optimal knots and the effective range need to be determined automatically.

    The algorithm:
    1. Start with initial narrow range around function center
    2. Extend boundaries by 10% if function value exceeds epsilon
       - Mark only the outermost boundary point as OLD (using add_old)
       - Mark intermediate points in extension as NEW (using regular add)
       - This ensures the spline from old points has proper support
       - Intermediate NEW points can be refined in the next iteration
    3. Refine internal knots for accuracy (operates on NEW points)
    4. Repeat until boundaries converge (function falls below epsilon)

    Key insight: Only the actual boundary point needs to be OLD. Intermediate
    points between the old and new boundary can be NEW, making them candidates
    for refinement while still maintaining valid spline support.
    """
    fss = Ncm.FunctionSampleSet.new(1)
    base_spline = Ncm.SplineCubicNotaknot.new()

    # Gaussian parameters
    center = 20.0
    std = 1.0

    # Initial range: narrow interval around center
    x_min = 18.0
    x_max = 22.0

    # Convergence parameters
    epsilon = 1e-6  # Stop extending when function < epsilon at boundaries
    reltol = 1e-4  # Relative tolerance for internal refinement
    abstol = 1e-8  # Absolute tolerance for internal refinement
    max_iterations = 20  # Safety limit

    # Initialize with at least 6 knots for cubic notaknot
    initial_knots = np.linspace(x_min, x_max, 8)
    for x in initial_knots:
        y = Ncm.Vector.new(1)
        y.set(0, gaussian(x, center, std))
        fss.add(x, y)

    print(f"\nInitial range: [{x_min:.4f}, {x_max:.4f}]")
    print(f"Initial knots: {fss.get_nsamples()}")

    # Phase 1: Extend boundaries until convergence
    iteration = 0
    boundaries_extended = True

    while boundaries_extended and iteration < max_iterations:
        iteration += 1
        boundaries_extended = False

        # Get current interval size
        current_x_min = fss.get_x_min()
        current_x_max = fss.get_x_max()
        interval_size = current_x_max - current_x_min
        extension = 0.1 * interval_size  # 10% extension

        # Check left boundary
        f_left = gaussian(current_x_min, center, std)
        if f_left > epsilon:
            # Extend left boundary
            new_x_min = current_x_min - extension
            # Add knots in extended region
            # Only the leftmost (boundary) point is marked as OLD
            # Intermediate points are NEW and can be refined
            extension_points = np.linspace(new_x_min, current_x_min, 4, endpoint=False)
            for i, x in enumerate(extension_points):
                y = Ncm.Vector.new(1)
                y.set(0, gaussian(x, center, std))
                if i == 0:  # Leftmost point is the new boundary
                    fss.add_old(x, y)  # Mark boundary as old
                else:
                    fss.add(x, y)  # Intermediate points are new
            boundaries_extended = True
            print(
                f"  Iteration {iteration}: Extended left to "
                f"{new_x_min:.4f} (f={f_left:.2e})"
            )

        # Check right boundary
        f_right = gaussian(current_x_max, center, std)
        if f_right > epsilon:
            # Extend right boundary
            new_x_max = current_x_max + extension
            # Add knots in extended region
            # Only the rightmost (boundary) point is marked as OLD
            # Intermediate points are NEW and can be refined
            extension_points = np.linspace(current_x_max, new_x_max, 4)[
                1:
            ]  # Skip current_x_max
            for i, x in enumerate(extension_points):
                y = Ncm.Vector.new(1)
                y.set(0, gaussian(x, center, std))
                if (
                    i == len(extension_points) - 1
                ):  # Rightmost point is the new boundary
                    fss.add_old(x, y)  # Mark boundary as old
                else:
                    fss.add(x, y)  # Intermediate points are new
            boundaries_extended = True
            print(
                f"  Iteration {iteration}: Extended right "
                f"to {new_x_max:.4f} (f={f_right:.2e})"
            )

        print(
            f"  After iteration {iteration}: range=[{fss.get_x_min():.4f}, "
            f"{fss.get_x_max():.4f}], knots={fss.get_nsamples()}"
        )

    # Phase 2: Refine internal knots once boundaries are established
    # Note: Boundary extensions were marked as OLD, so refine() can work directly
    print("  Refining internal knots...")
    fss.refine(reltol, abstol, base_spline)

    # Final statistics
    final_x_min = fss.get_x_min()
    final_x_max = fss.get_x_max()
    final_nknots = fss.get_nsamples()

    print(f"\nFinal range: [{final_x_min:.4f}, {final_x_max:.4f}]")
    print(f"Final knots: {final_nknots}")
    print(
        f"Range in sigmas: [{(final_x_min-center)/std:.2f}, "
        f"{(final_x_max-center)/std:.2f}]"
    )

    # Verify the final range is reasonable (should capture ~99.7% of Gaussian, i.e., ±3
    # sigma) But not too wide (shouldn't go beyond ±5 sigma where function is truly
    # negligible)
    assert (
        final_x_min < center - 2.5 * std
    ), "Left boundary should extend at least 2.5 sigma"
    assert (
        final_x_max > center + 2.5 * std
    ), "Right boundary should extend at least 2.5 sigma"
    assert (
        final_x_min > center - 6.0 * std
    ), "Left boundary shouldn't extend beyond 6 sigma"
    assert (
        final_x_max < center + 6.0 * std
    ), "Right boundary shouldn't extend beyond 6 sigma"

    # Verify boundaries have converged (function is below epsilon)
    assert gaussian(final_x_min, center, std) <= epsilon
    assert gaussian(final_x_max, center, std) <= epsilon

    # Verify the spline approximates the Gaussian well
    sv = fss.to_spline_vec(base_spline)
    test_points = np.linspace(final_x_min + 0.1, final_x_max - 0.1, 50)

    max_error = 0.0
    for x in test_points:
        y_true = gaussian(x, center, std)
        y_spline = sv.eval_array(x)[0]
        error = abs(y_true - y_spline)
        max_error = max(max_error, error)

        # Check relative error where function is significant (> 1000 * epsilon)
        if y_true > 1000 * epsilon:
            rel_error = error / y_true
            assert (
                rel_error < reltol * 100
            ), f"Relative error too large at x={x}: {rel_error}"

    print(f"Maximum absolute error: {max_error:.2e}")
    assert max_error < 1e-3, "Maximum error should be small"


def test_auto_range_asymmetric_function() -> None:
    """Test automatic range detection for an asymmetric localized function.

    Uses an exponentially modified Gaussian (EMG) which has asymmetric tails.
    """
    fss = Ncm.FunctionSampleSet.new(1)
    base_spline = Ncm.SplineCubicNotaknot.new()

    # EMG parameters: Gaussian(mu, sigma) convolved with exponential(lambda)
    mu = 20.0
    sigma = 1.0
    lambda_param = 0.5  # Controls exponential decay

    def emg(x: float) -> float:
        """Exponentially modified Gaussian."""
        # Simplified EMG for testing
        gauss = math.exp(-0.5 * ((x - mu) / sigma) ** 2) / (
            sigma * math.sqrt(2.0 * math.pi)
        )
        if x > mu:
            gauss *= math.exp(-lambda_param * (x - mu))
        return gauss

    # Initial range
    x_min = 18.0
    x_max = 22.0

    # Parameters
    epsilon = 1e-6
    reltol = 1e-3
    abstol = 1e-8
    max_iterations = 20

    # Initialize
    for x in np.linspace(x_min, x_max, 8):
        y = Ncm.Vector.new(1)
        y.set(0, emg(x))
        fss.add(x, y)

    # Phase 1: Extend boundaries
    iteration = 0
    boundaries_extended = True

    while boundaries_extended and iteration < max_iterations:
        iteration += 1
        boundaries_extended = False

        current_x_min = fss.get_x_min()
        current_x_max = fss.get_x_max()
        interval_size = current_x_max - current_x_min
        extension = 0.1 * interval_size

        # Check and extend left
        if emg(current_x_min) > epsilon:
            new_x_min = current_x_min - extension
            extension_points = np.linspace(new_x_min, current_x_min, 4, endpoint=False)
            for i, x in enumerate(extension_points):
                y = Ncm.Vector.new(1)
                y.set(0, emg(x))
                if i == 0:  # Leftmost point is the boundary
                    fss.add_old(x, y)
                else:
                    fss.add(x, y)
            boundaries_extended = True

        # Check and extend right (asymmetric tail)
        if emg(current_x_max) > epsilon:
            new_x_max = current_x_max + extension
            extension_points = np.linspace(current_x_max, new_x_max, 4)[1:]
            for i, x in enumerate(extension_points):
                y = Ncm.Vector.new(1)
                y.set(0, emg(x))
                if i == len(extension_points) - 1:  # Rightmost point is the boundary
                    fss.add_old(x, y)
                else:
                    fss.add(x, y)
            boundaries_extended = True

    # Phase 2: Refine internal knots
    print("  Refining internal knots...")
    fss.refine(reltol, abstol, base_spline)

    final_x_min = fss.get_x_min()
    final_x_max = fss.get_x_max()

    print(f"\nAsymmetric function final range: [{final_x_min:.4f}, {final_x_max:.4f}]")
    print(f"Final knots: {fss.get_nsamples()}")

    # Verify boundaries
    assert emg(final_x_min) <= epsilon
    assert emg(final_x_max) <= epsilon

    # Due to asymmetry, right extension might be different from left
    left_extent = abs(final_x_min - mu)
    right_extent = abs(final_x_max - mu)

    print(f"Left extent: {left_extent:.2f}, Right extent: {right_extent:.2f}")

    # Both extents should be reasonable but may differ
    assert left_extent > 2.0 * sigma
    assert right_extent > 2.0 * sigma
    assert left_extent < 10.0 * sigma
    assert right_extent < 10.0 * sigma


def test_auto_range_multi_component_gaussian() -> None:
    """Test automatic range detection for multi-component Gaussian functions.

    Uses two Gaussian components with different widths - range should
    accommodate the wider one.
    """
    fss = Ncm.FunctionSampleSet.new(2)
    base_spline = Ncm.SplineCubicNotaknot.new()

    # Two Gaussians centered at same point but different widths
    center = 20.0
    std_narrow = 0.5  # Narrow component
    std_wide = 2.0  # Wide component

    def multi_gaussian(x: float) -> tuple[float, float]:
        narrow = gaussian(x, center, std_narrow)
        wide = gaussian(x, center, std_wide)
        return (narrow, wide)

    # Initial range
    x_min = 18.0
    x_max = 22.0

    epsilon = 1e-6
    reltol = 1e-3
    abstol = 1e-8
    max_iterations = 20

    # Initialize
    for x in np.linspace(x_min, x_max, 8):
        y = Ncm.Vector.new(2)
        vals = multi_gaussian(x)
        y.set(0, vals[0])
        y.set(1, vals[1])
        fss.add(x, y)

    # Phase 1: Extend boundaries
    iteration = 0
    boundaries_extended = True

    while boundaries_extended and iteration < max_iterations:
        iteration += 1
        boundaries_extended = False

        current_x_min = fss.get_x_min()
        current_x_max = fss.get_x_max()
        interval_size = current_x_max - current_x_min
        extension = 0.1 * interval_size

        # Check boundary - use max of both components
        vals_left = multi_gaussian(current_x_min)
        f_left = max(vals_left)

        if f_left > epsilon:
            new_x_min = current_x_min - extension
            extension_points = np.linspace(new_x_min, current_x_min, 4, endpoint=False)
            for i, x in enumerate(extension_points):
                y = Ncm.Vector.new(2)
                vals = multi_gaussian(x)
                y.set(0, vals[0])
                y.set(1, vals[1])
                if i == 0:  # Leftmost point is the boundary
                    fss.add_old(x, y)
                else:
                    fss.add(x, y)
            boundaries_extended = True

        vals_right = multi_gaussian(current_x_max)
        f_right = max(vals_right)

        if f_right > epsilon:
            new_x_max = current_x_max + extension
            extension_points = np.linspace(current_x_max, new_x_max, 4)[1:]
            for i, x in enumerate(extension_points):
                y = Ncm.Vector.new(2)
                vals = multi_gaussian(x)
                y.set(0, vals[0])
                y.set(1, vals[1])
                if i == len(extension_points) - 1:  # Rightmost point is the boundary
                    fss.add_old(x, y)
                else:
                    fss.add(x, y)
            boundaries_extended = True

    # Phase 2: Refine internal knots
    print("  Refining internal knots...")
    fss.refine(reltol, abstol, base_spline)

    final_x_min = fss.get_x_min()
    final_x_max = fss.get_x_max()

    print(f"\nMulti-component final range: [{final_x_min:.4f}, {final_x_max:.4f}]")
    print(f"Final knots: {fss.get_nsamples()}")

    # Range should be determined by the wider Gaussian
    # Should capture at least 3 sigma of the wide component
    assert final_x_min < center - 3.0 * std_wide
    assert final_x_max > center + 3.0 * std_wide

    # But shouldn't be excessively wide
    assert final_x_min > center - 6.0 * std_wide
    assert final_x_max < center + 6.0 * std_wide

    # Verify both components are below epsilon at boundaries
    vals_left = multi_gaussian(final_x_min)
    vals_right = multi_gaussian(final_x_max)

    assert max(vals_left) <= epsilon or final_x_min <= center - 5.0 * std_wide
    assert max(vals_right) <= epsilon or final_x_max >= center + 5.0 * std_wide

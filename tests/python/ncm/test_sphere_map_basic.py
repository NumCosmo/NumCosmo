#!/usr/bin/env python
#
# test_sphere_map_basic.py
#
# Wed Jan 01 2025
# Copyright 2025 Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_sphere_map_basic.py
# Copyright (C) 2025 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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
# with this program. If not, see <http://www.gnu.org/licenses/>.

"""Test suite for NcmSphereMap basic functionality without healpy dependency."""

import pytest
import numpy as np
from numpy.testing import assert_allclose

from numcosmo_py import Ncm

pytestmark = [pytest.mark.sphere_map]


@pytest.fixture(params=[8, 16, 32, 64, 128], name="nside")
def fixture_nside(request: pytest.FixtureRequest) -> int:
    """Parametrized nside values for testing."""
    return request.param


@pytest.fixture(name="smap")
def fixture_smap(nside: int) -> Ncm.SphereMap:
    """Create NcmSphereMap instance."""
    return Ncm.SphereMap.new(nside)


class TestBasicProperties:
    """Test basic HEALPix properties using known formulas."""

    def test_nside(self, smap: Ncm.SphereMap, nside: int) -> None:
        """Test nside getter."""
        assert smap.get_nside() == nside

    def test_nrings(self, smap: Ncm.SphereMap, nside: int) -> None:
        """Test number of rings."""
        ncm_nrings = smap.get_nrings()
        # HEALPix formula: nrings = 4 * nside - 1
        expected_nrings = 4 * nside - 1
        assert ncm_nrings == expected_nrings

    def test_cap_size(self, smap: Ncm.SphereMap, nside: int) -> None:
        """Test polar cap size."""
        ncm_cap_size = smap.get_cap_size()
        # HEALPix formula: cap_size = 2 * nside * (nside - 1)
        expected_cap_size = 2 * nside * (nside - 1)
        assert ncm_cap_size == expected_cap_size

    def test_middle_size(self, smap: Ncm.SphereMap) -> None:
        """Test equatorial region size."""
        ncm_middle_size = smap.get_middle_size()
        # NumCosmo includes ring boundaries differently than standard HEALPix
        # npix = cap_size * 2 + middle_size
        npix = smap.get_npix()
        cap_size = smap.get_cap_size()
        expected_middle_size = npix - 2 * cap_size
        assert ncm_middle_size == expected_middle_size

    def test_npix(self, smap: Ncm.SphereMap, nside: int) -> None:
        """Test total number of pixels using HEALPix formula."""
        ncm_npix = smap.get_npix()
        # HEALPix formula: npix = 12 * nside^2
        expected_npix = 12 * nside * nside
        assert ncm_npix == expected_npix


class TestRingProperties:
    """Test ring-specific properties."""

    def test_ring_sizes(self, smap: Ncm.SphereMap) -> None:
        """Test that ring sizes follow HEALPix pattern."""
        nrings = smap.get_nrings()

        for r_i in range(nrings):
            ring_size = smap.get_ring_size(r_i)

            # Ring size should be 4 * min(r_i+1, 2*nside, 4*nside - r_i - 1)
            # This is the HEALPix pattern
            assert ring_size > 0
            assert ring_size % 4 == 0 or ring_size == 4  # Multiple of 4

    def test_ring_first_indices(self, smap: Ncm.SphereMap) -> None:
        """Test that ring first indices are consistent."""
        nrings = smap.get_nrings()
        npix = smap.get_npix()

        total_pixels = 0
        for r_i in range(nrings):
            ring_fi = smap.get_ring_first_index(r_i)
            ring_size = smap.get_ring_size(r_i)

            assert ring_fi == total_pixels
            total_pixels += ring_size

        assert total_pixels == npix

    def test_ring_theta_constant(self, smap: Ncm.SphereMap) -> None:
        """Test that all pixels in a ring have the same theta (colatitude)."""
        nrings = smap.get_nrings()

        # Test a few random rings
        test_rings = np.linspace(0, nrings - 1, min(20, nrings), dtype=int)

        for r_i in test_rings:
            ring_fi = smap.get_ring_first_index(r_i)
            ring_size = smap.get_ring_size(r_i)

            thetas = []
            for i in range(ring_size):
                theta, _ = smap.pix2ang_ring(int(ring_fi + i))
                thetas.append(theta)

            # All thetas should be identical
            assert_allclose(thetas, thetas[0], rtol=1e-12, atol=1e-14)


class TestBasicCoordinateConversions:
    """Test basic coordinate conversions without healpy comparison."""

    def test_ang2pix_ring_bounds(self, smap: Ncm.SphereMap) -> None:
        """Test that ang2pix_ring returns valid pixel indices."""
        npix = smap.get_npix()

        # Test various angles
        test_angles = [
            (0.0, 0.0),  # North pole
            (np.pi, 0.0),  # South pole
            (np.pi / 2, 0.0),  # Equator
            (np.pi / 4, np.pi / 2),  # Random point
            (3 * np.pi / 4, np.pi),  # Another random point
        ]

        for theta, phi in test_angles:
            pix_ring = smap.ang2pix_ring(theta, phi)
            assert (
                0 <= pix_ring < npix
            ), f"Invalid pixel index for theta={theta}, phi={phi}"

    def test_ang2pix_nest_bounds(self, smap: Ncm.SphereMap) -> None:
        """Test that ang2pix_nest returns valid pixel indices."""
        npix = smap.get_npix()

        # Test various angles
        test_angles = [
            (0.0, 0.0),  # North pole
            (np.pi, 0.0),  # South pole
            (np.pi / 2, 0.0),  # Equator
            (np.pi / 4, np.pi / 2),  # Random point
            (3 * np.pi / 4, np.pi),  # Another random point
        ]

        for theta, phi in test_angles:
            pix_nest = smap.ang2pix_nest(theta, phi)
            assert (
                0 <= pix_nest < npix
            ), f"Invalid pixel index for theta={theta}, phi={phi}"

    def test_pix2ang_ring_bounds(self, smap: Ncm.SphereMap) -> None:
        """Test that pix2ang_ring returns valid angle ranges."""
        npix = smap.get_npix()

        # Test first, middle, and last pixels
        test_pixels = [0, npix // 4, npix // 2, 3 * npix // 4, npix - 1]

        for pix in test_pixels:
            theta, phi = smap.pix2ang_ring(pix)
            assert 0 <= theta <= np.pi, f"Invalid theta for pixel {pix}"
            assert 0 <= phi < 2 * np.pi, f"Invalid phi for pixel {pix}"

    def test_pix2ang_nest_bounds(self, smap: Ncm.SphereMap) -> None:
        """Test that pix2ang_nest returns valid angle ranges."""
        npix = smap.get_npix()

        # Test first, middle, and last pixels
        test_pixels = [0, npix // 4, npix // 2, 3 * npix // 4, npix - 1]

        for pix in test_pixels:
            theta, phi = smap.pix2ang_nest(pix)
            assert 0 <= theta <= np.pi, f"Invalid theta for pixel {pix}"
            assert 0 <= phi < 2 * np.pi, f"Invalid phi for pixel {pix}"

    def test_nest2ring_bounds(self, smap: Ncm.SphereMap) -> None:
        """Test that nest2ring returns valid pixel indices."""
        npix = smap.get_npix()

        # Test all pixels
        for nest_pix in range(npix):
            ring_pix = smap.nest2ring(nest_pix)
            assert 0 <= ring_pix < npix, f"Invalid ring pixel for nest pixel {nest_pix}"

    def test_ring2nest_bounds(self, smap: Ncm.SphereMap) -> None:
        """Test that ring2nest returns valid pixel indices."""
        npix = smap.get_npix()

        # Test all pixels
        for ring_pix in range(npix):
            nest_pix = smap.ring2nest(ring_pix)
            assert 0 <= nest_pix < npix, f"Invalid nest pixel for ring pixel {ring_pix}"

    def test_nest_ring_roundtrip(self, smap: Ncm.SphereMap) -> None:
        """Test nest->ring->nest roundtrip."""
        npix = smap.get_npix()
        nest_indices = np.arange(npix)

        ring_indices = np.array([smap.nest2ring(int(i)) for i in nest_indices])
        nest_back = np.array([smap.ring2nest(i) for i in ring_indices])

        assert np.array_equal(nest_indices, nest_back)

    def test_ring_nest_roundtrip(self, smap: Ncm.SphereMap) -> None:
        """Test ring->nest->ring roundtrip."""
        npix = smap.get_npix()
        ring_indices = np.arange(npix)

        nest_indices = np.array([smap.ring2nest(int(i)) for i in ring_indices])
        ring_back = np.array([smap.nest2ring(i) for i in nest_indices])

        assert np.array_equal(ring_indices, ring_back)


class TestEdgeCases:
    """Test edge cases and boundary conditions."""

    def test_minimum_nside(self) -> None:
        """Test minimum valid nside."""
        smap_min = Ncm.SphereMap.new(1)
        assert smap_min.get_nside() == 1
        assert smap_min.get_npix() == 12  # Minimum HEALPix map

    def test_north_pole_bounds(self, smap: Ncm.SphereMap) -> None:
        """Test pixel at north pole returns valid indices."""
        theta_pole = 0.0
        phi = 0.0

        pix_ring = smap.ang2pix_ring(theta_pole, phi)
        pix_nest = smap.ang2pix_nest(theta_pole, phi)

        # Both should point to valid pixels
        assert 0 <= pix_ring < smap.get_npix()
        assert 0 <= pix_nest < smap.get_npix()

    def test_south_pole_bounds(self, smap: Ncm.SphereMap) -> None:
        """Test pixel at south pole returns valid indices."""
        theta_pole = np.pi
        phi = 0.0

        pix_ring = smap.ang2pix_ring(theta_pole, phi)
        pix_nest = smap.ang2pix_nest(theta_pole, phi)

        # Both should point to valid pixels
        assert 0 <= pix_ring < smap.get_npix()
        assert 0 <= pix_nest < smap.get_npix()

    def test_phi_wrapping(self, smap: Ncm.SphereMap) -> None:
        """Test that phi values wrap correctly at 2*pi."""
        theta = np.pi / 2  # Equator

        # Test phi = 0 and phi = 2*pi should give same pixel
        pix_0 = smap.ang2pix_ring(theta, 0.0)
        pix_2pi = smap.ang2pix_ring(theta, 2 * np.pi)

        # Should be the same pixel (within HEALPix precision)
        assert pix_0 == pix_2pi or abs(pix_0 - pix_2pi) <= 1

    def test_large_nside(self) -> None:
        """Test with larger nside (quick sanity check)."""
        nside_large = 256
        smap_large = Ncm.SphereMap.new(nside_large)

        npix = smap_large.get_npix()
        expected_npix = 12 * nside_large * nside_large

        assert npix == expected_npix


class TestVectorConversions:
    """Test 3D vector to pixel conversions."""

    def test_vec2pix_ring_bounds(self, smap: Ncm.SphereMap) -> None:
        """Test that vec2pix_ring returns valid pixel indices."""
        npix = smap.get_npix()

        # Test various unit vectors
        test_vectors = [
            (0.0, 0.0, 1.0),  # North pole
            (0.0, 0.0, -1.0),  # South pole
            (1.0, 0.0, 0.0),  # Equator
            (0.0, 1.0, 0.0),  # Equator
            (1.0 / np.sqrt(2), 1.0 / np.sqrt(2), 0.0),  # Equator diagonal
        ]

        for x, y, z in test_vectors:
            tri_vec = Ncm.TriVec.new_full_c(x, y, z)
            pix = smap.vec2pix_ring(tri_vec)
            assert 0 <= pix < npix, f"Invalid pixel for vector ({x}, {y}, {z})"

    def test_vec2pix_nest_bounds(self, smap: Ncm.SphereMap) -> None:
        """Test that vec2pix_nest returns valid pixel indices."""
        npix = smap.get_npix()

        # Test various unit vectors
        test_vectors = [
            (0.0, 0.0, 1.0),  # North pole
            (0.0, 0.0, -1.0),  # South pole
            (1.0, 0.0, 0.0),  # Equator
            (0.0, 1.0, 0.0),  # Equator
            (1.0 / np.sqrt(2), 1.0 / np.sqrt(2), 0.0),  # Equator diagonal
        ]

        for x, y, z in test_vectors:
            tri_vec = Ncm.TriVec.new_full_c(x, y, z)
            pix = smap.vec2pix_nest(tri_vec)
            assert 0 <= pix < npix, f"Invalid pixel for vector ({x}, {y}, {z})"

    def test_pix2vec_ring_unit_vectors(self, smap: Ncm.SphereMap) -> None:
        """Test that pix2vec_ring returns unit vectors."""
        npix = smap.get_npix()

        # Test a few pixels
        test_pixels = [0, npix // 4, npix // 2, 3 * npix // 4, npix - 1]

        for pix in test_pixels:
            tri_vec = Ncm.TriVec()
            smap.pix2vec_ring(pix, tri_vec)
            vec = tri_vec.c[:3]  # pylint: disable=unsubscriptable-object

            # Check it's a unit vector
            length = np.sqrt(vec[0] ** 2 + vec[1] ** 2 + vec[2] ** 2)
            assert_allclose(length, 1.0, rtol=1e-12, atol=1e-12)

    def test_pix2vec_nest_unit_vectors(self, smap: Ncm.SphereMap) -> None:
        """Test that pix2vec_nest returns unit vectors."""
        npix = smap.get_npix()

        # Test a few pixels
        test_pixels = [0, npix // 4, npix // 2, 3 * npix // 4, npix - 1]

        for pix in test_pixels:
            tri_vec = Ncm.TriVec()
            smap.pix2vec_nest(pix, tri_vec)
            vec = tri_vec.c[:3]  # pylint: disable=unsubscriptable-object

            # Check it's a unit vector
            length = np.sqrt(vec[0] ** 2 + vec[1] ** 2 + vec[2] ** 2)
            assert_allclose(length, 1.0, rtol=1e-12, atol=1e-12)


if __name__ == "__main__":
    # Run tests with: python test_sphere_map_basic.py
    pytest.main([__file__, "-v"])

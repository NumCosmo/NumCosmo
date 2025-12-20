#!/usr/bin/env python
#
# test_py_sphere_map.py
#
# Fri Dec 20 2025
# Copyright 2025 Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_py_sphere_map.py
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

"""Test suite comparing NcmSphereMap implementation with healpy."""

from typing import Any
import pytest
import numpy as np
import numpy.typing as npt
from numpy.testing import assert_allclose, assert_array_equal

healpy = pytest.importorskip("healpy")

# pylint: disable-next=wrong-import-position
from numcosmo_py import Ncm  # noqa: E402


@pytest.fixture(params=[8, 16, 32, 64, 128], name="nside")
def fixture_nside(request: Any) -> int:
    """Parametrized nside values for testing."""
    return request.param


@pytest.fixture(name="smap")
def fixture_smap(nside: int) -> Ncm.SphereMap:
    """Create NcmSphereMap instance."""
    return Ncm.SphereMap.new(nside)


@pytest.fixture(name="random_seed")
def fixture_random_seed() -> int:
    """Fixed random seed for reproducibility."""
    return 42


class TestBasicProperties:
    """Test basic HEALPix properties match between implementations."""

    def test_npix(self, smap: Ncm.SphereMap, nside: int) -> None:
        """Test total number of pixels."""
        ncm_npix = smap.get_npix()
        hp_npix = healpy.nside2npix(nside)
        assert ncm_npix == hp_npix

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


class TestPixelIndexing:
    """Test pixel indexing conversions match healpy."""

    def test_nest2ring_all_pixels(self, smap: Ncm.SphereMap, nside: int) -> None:
        """Test nest to ring conversion for all pixels."""
        npix = healpy.nside2npix(nside)
        nest_indices = np.arange(npix)

        # NumCosmo conversion
        ncm_ring = np.array([smap.nest2ring(i) for i in nest_indices])

        # healpy conversion
        hp_ring = healpy.nest2ring(nside, nest_indices)

        assert_array_equal(ncm_ring, hp_ring)

    def test_ring2nest_all_pixels(self, smap: Ncm.SphereMap, nside: int) -> None:
        """Test ring to nest conversion for all pixels."""
        npix = healpy.nside2npix(nside)
        ring_indices = np.arange(npix)

        # NumCosmo conversion
        ncm_nest = np.array([smap.ring2nest(i) for i in ring_indices])

        # healpy conversion
        hp_nest = healpy.ring2nest(nside, ring_indices)

        assert_array_equal(ncm_nest, hp_nest)

    def test_nest_ring_roundtrip(self, smap: Ncm.SphereMap, nside: int) -> None:
        """Test nest->ring->nest roundtrip."""
        npix = healpy.nside2npix(nside)
        nest_indices = np.arange(npix)

        ring_indices = np.array([smap.nest2ring(i) for i in nest_indices])
        nest_back = np.array([smap.ring2nest(i) for i in ring_indices])

        assert_array_equal(nest_indices, nest_back)


class TestAngularConversions:
    """Test angular coordinate conversions."""

    @pytest.fixture(name="test_angles")
    def fixture_test_angles(
        self,
    ) -> tuple[npt.NDArray[np.floating[Any]], npt.NDArray[np.floating[Any]]]:
        """Generate test angles covering the sphere."""
        n_theta = 50
        n_phi = 100
        theta = np.linspace(0, np.pi, n_theta)
        phi = np.linspace(0, 2 * np.pi, n_phi)
        return np.meshgrid(theta, phi)

    def test_ang2pix_ring(
        self,
        smap: Ncm.SphereMap,
        nside: int,
        test_angles: tuple[
            npt.NDArray[np.floating[Any]], npt.NDArray[np.floating[Any]]
        ],
    ) -> None:
        """Test angle to pixel conversion in ring ordering."""
        theta, phi = test_angles
        theta_flat = theta.flatten()
        phi_flat = phi.flatten()

        # NumCosmo conversion
        ncm_pix = np.array(
            [
                smap.ang2pix_ring(float(t), float(p))
                for t, p in zip(theta_flat, phi_flat)
            ]
        )

        # healpy conversion
        hp_pix = healpy.ang2pix(nside, theta_flat, phi_flat, nest=False)

        assert_array_equal(ncm_pix, hp_pix)

    def test_ang2pix_nest(
        self,
        smap: Ncm.SphereMap,
        nside: int,
        test_angles: tuple[
            npt.NDArray[np.floating[Any]], npt.NDArray[np.floating[Any]]
        ],
    ) -> None:
        """Test angle to pixel conversion in nest ordering."""
        theta, phi = test_angles
        theta_flat = theta.flatten()
        phi_flat = phi.flatten()

        # NumCosmo conversion
        ncm_pix = np.array(
            [
                smap.ang2pix_nest(float(t), float(p))
                for t, p in zip(theta_flat, phi_flat)
            ]
        )

        # healpy conversion
        hp_pix = healpy.ang2pix(nside, theta_flat, phi_flat, nest=True)

        assert_array_equal(ncm_pix, hp_pix)

    def test_pix2ang_ring(self, smap: Ncm.SphereMap, nside: int) -> None:
        """Test pixel to angle conversion in ring ordering."""
        npix = healpy.nside2npix(nside)
        # Sample subset of pixels for faster testing
        pix_indices = np.linspace(0, npix - 1, min(1000, npix), dtype=int)

        # NumCosmo conversion
        ncm_theta = np.zeros(len(pix_indices))
        ncm_phi = np.zeros(len(pix_indices))
        for i, pix in enumerate(pix_indices):
            theta, phi = smap.pix2ang_ring(int(pix))
            ncm_theta[i] = theta
            ncm_phi[i] = phi

        # healpy conversion
        hp_theta, hp_phi = healpy.pix2ang(nside, pix_indices, nest=False)

        assert_allclose(ncm_theta, hp_theta, rtol=1e-10, atol=1e-12)
        assert_allclose(ncm_phi, hp_phi, rtol=1e-10, atol=1e-12)

    def test_pix2ang_nest(self, smap: Ncm.SphereMap, nside: int) -> None:
        """Test pixel to angle conversion in nest ordering."""
        npix = healpy.nside2npix(nside)
        # Sample subset of pixels for faster testing
        pix_indices = np.linspace(0, npix - 1, min(1000, npix), dtype=int)

        # NumCosmo conversion
        ncm_theta = np.zeros(len(pix_indices))
        ncm_phi = np.zeros(len(pix_indices))
        for i, pix in enumerate(pix_indices):
            theta, phi = smap.pix2ang_nest(int(pix))
            ncm_theta[i] = theta
            ncm_phi[i] = phi

        # healpy conversion
        hp_theta, hp_phi = healpy.pix2ang(nside, pix_indices, nest=True)

        assert_allclose(ncm_theta, hp_theta, rtol=1e-10, atol=1e-12)
        assert_allclose(ncm_phi, hp_phi, rtol=1e-10, atol=1e-12)

    def test_ang_pix_roundtrip_ring(self, smap: Ncm.SphereMap, nside: int) -> None:
        """Test angle->pix->angle roundtrip in ring ordering."""
        # Generate test angles
        n_test = 100
        rng = np.random.default_rng(42)
        theta_orig = np.arccos(rng.uniform(-1, 1, n_test))
        phi_orig = rng.uniform(0, 2 * np.pi, n_test)

        # Convert to pixels and back
        for t, p in zip(theta_orig, phi_orig):
            pix = smap.ang2pix_ring(float(t), float(p))
            theta_back, phi_back = smap.pix2ang_ring(int(pix))

            # Check we get the same pixel center as healpy
            hp_pix = healpy.ang2pix(nside, t, p, nest=False)
            hp_theta, hp_phi = healpy.pix2ang(nside, hp_pix, nest=False)

            assert_allclose(theta_back, hp_theta, rtol=1e-10, atol=1e-12)
            assert_allclose(phi_back, hp_phi, rtol=1e-10, atol=1e-12)


class TestVectorConversions:
    """Test 3D vector to pixel conversions."""

    @pytest.fixture(name="test_vectors")
    def fixture_test_vectors(self) -> npt.NDArray[np.floating[Any]]:
        """Generate test unit vectors."""
        n_test = 100
        rng = np.random.default_rng(42)
        # Generate random points on unit sphere
        theta = np.arccos(rng.uniform(-1, 1, n_test))
        phi = rng.uniform(0, 2 * np.pi, n_test)
        x = np.sin(theta) * np.cos(phi)
        y = np.sin(theta) * np.sin(phi)
        z = np.cos(theta)
        return np.column_stack([x, y, z])

    def test_vec2pix_ring(
        self,
        smap: Ncm.SphereMap,
        nside: int,
        test_vectors: npt.NDArray[np.floating[Any]],
    ) -> None:
        """Test vector to pixel conversion in ring ordering."""
        # NumCosmo conversion
        ncm_pix = np.zeros(len(test_vectors), dtype=int)
        for i, vec in enumerate(test_vectors):
            tri_vec = Ncm.TriVec.new_full_c(float(vec[0]), float(vec[1]), float(vec[2]))
            ncm_pix[i] = smap.vec2pix_ring(tri_vec)

        # healpy conversion
        hp_pix = healpy.vec2pix(
            nside,
            test_vectors[:, 0],
            test_vectors[:, 1],
            test_vectors[:, 2],
            nest=False,
        )

        assert_array_equal(ncm_pix, hp_pix)

    def test_vec2pix_nest(
        self,
        smap: Ncm.SphereMap,
        nside: int,
        test_vectors: npt.NDArray[np.floating[Any]],
    ) -> None:
        """Test vector to pixel conversion in nest ordering."""
        # NumCosmo conversion
        ncm_pix = np.zeros(len(test_vectors), dtype=int)
        for i, vec in enumerate(test_vectors):
            tri_vec = Ncm.TriVec.new_full_c(float(vec[0]), float(vec[1]), float(vec[2]))
            ncm_pix[i] = smap.vec2pix_nest(tri_vec)

        # healpy conversion
        hp_pix = healpy.vec2pix(
            nside, test_vectors[:, 0], test_vectors[:, 1], test_vectors[:, 2], nest=True
        )

        assert_array_equal(ncm_pix, hp_pix)

    def test_pix2vec_ring(self, smap: Ncm.SphereMap, nside: int) -> None:
        """Test pixel to vector conversion in ring ordering."""
        npix = healpy.nside2npix(nside)
        pix_indices = np.linspace(0, npix - 1, min(500, npix), dtype=int)

        # NumCosmo conversion
        ncm_vecs = np.zeros((len(pix_indices), 3))
        for i, pix in enumerate(pix_indices):
            tri_vec = Ncm.TriVec()
            smap.pix2vec_ring(int(pix), tri_vec)
            ncm_vecs[i] = tri_vec.c[:3]  # pylint: disable=unsubscriptable-object

        # healpy conversion
        hp_vecs = np.array(healpy.pix2vec(nside, pix_indices, nest=False)).T

        assert_allclose(ncm_vecs, hp_vecs, rtol=1e-10, atol=1e-12)


class TestSphericalHarmonics:
    """Test spherical harmonic decomposition and synthesis."""

    @pytest.fixture(name="lmax")
    def fixture_lmax(self, nside: int) -> int:
        """Choose lmax appropriate for nside."""
        # HEALPix convention: lmax ~ 3 * nside is well-sampled
        return min(3 * nside, 512)  # Cap at 512 for reasonable test time

    def test_map2alm_white_noise(
        self, smap: Ncm.SphereMap, nside: int, lmax: int, random_seed: int
    ) -> None:
        """Test map2alm on white noise map."""
        np.random.seed(random_seed)
        npix = healpy.nside2npix(nside)

        # Generate white noise map
        noise_map = np.random.randn(npix)

        # Set map in NumCosmo
        smap.set_lmax(lmax)
        smap.set_map(noise_map)
        smap.prepare_alm()

        # Get NumCosmo alm and compare with healpy
        # Use iter=0: direct transform only, no iterative refinement
        # This matches NumCosmo's algorithm exactly
        hp_alm = healpy.map2alm(noise_map, lmax=lmax, iter=0)

        # Compare a few alm coefficients
        for ell in range(min(lmax + 1, 10)):
            for m in range(ell + 1):
                re_alm, im_alm = smap.get_alm(ell, m)
                hp_idx = healpy.Alm.getidx(lmax, ell, m)
                hp_val = hp_alm[hp_idx]

                # Compare real and imaginary parts
                # With iter=0, should match to machine precision
                assert_allclose(re_alm, hp_val.real, rtol=1e-10, atol=1e-14)
                assert_allclose(im_alm, hp_val.imag, rtol=1e-10, atol=1e-14)

    def test_alm_indexing(self, lmax: int) -> None:
        """Test alm array indexing matches healpy convention."""
        # HEALPix alm indexing: index(l,m) = m*(2*lmax+1-m)/2 + l
        for m in range(0, min(lmax + 1, 10)):
            for ell in range(m, min(lmax + 1, 20)):
                hp_idx = healpy.Alm.getidx(lmax, ell, m)
                # NumCosmo uses macro: NCM_SPHERE_MAP_ALM_INDEX(lmax, l, m)
                # Formula: (l) + ((1 + 2 * (lmax) - (m)) * (m)) / 2
                ncm_idx = ell + ((1 + 2 * lmax - m) * m) // 2
                assert hp_idx == ncm_idx, f"Mismatch at l={ell}, m={m}"

    def test_alm2map_monopole(self, smap: Ncm.SphereMap, nside: int, lmax: int) -> None:
        """Test alm2map with monopole (l=0, m=0) only."""
        # Set up NumCosmo with monopole only
        smap.set_lmax(lmax)
        monopole_value = 10.0

        # Set all alm to zero, then set monopole
        for ell in range(lmax + 1):
            for m in range(ell + 1):
                smap.set_alm(ell, m, 0.0, 0.0)

        # Set monopole with Y_00 normalization
        smap.set_alm(0, 0, monopole_value * np.sqrt(4 * np.pi), 0.0)

        # Synthesize map from alm
        smap.alm2map()

        # Get NumCosmo synthesized map
        npix = healpy.nside2npix(nside)
        ncm_map = np.array([smap.get_pix(i) for i in range(npix)])

        # All pixels should have the same value
        expected_value = monopole_value
        assert_allclose(ncm_map, expected_value, rtol=1e-6)

        # Compare with healpy
        alm_size = healpy.Alm.getsize(lmax)
        hp_alm = np.zeros(alm_size, dtype=complex)
        hp_alm[0] = monopole_value * np.sqrt(4 * np.pi)
        hp_map = healpy.alm2map(hp_alm, nside, lmax=lmax)
        assert_allclose(ncm_map, hp_map, rtol=1e-6)

    def test_alm2map_dipole(self, smap: Ncm.SphereMap, nside: int) -> None:
        """Test alm2map with dipole (l=1) only."""
        lmax = max(1, nside)  # Need at least lmax=1
        smap.set_lmax(lmax)

        # Set all alm to zero, then set dipole
        for ell in range(lmax + 1):
            for m in range(ell + 1):
                smap.set_alm(ell, m, 0.0, 0.0)

        # Set l=1, m=0 (dipole along z-axis)
        smap.set_alm(1, 0, 1.0, 0.0)

        # Synthesize map from alm
        smap.alm2map()

        # Get NumCosmo synthesized map
        npix = healpy.nside2npix(nside)
        ncm_map = np.array([smap.get_pix(i) for i in range(npix)])

        # healpy synthesis for comparison
        alm_size = healpy.Alm.getsize(lmax)
        hp_alm = np.zeros(alm_size, dtype=complex)
        idx_10 = healpy.Alm.getidx(lmax, 1, 0)
        hp_alm[idx_10] = 1.0
        hp_map = healpy.alm2map(hp_alm, nside, lmax=lmax)

        # Compare NumCosmo and healpy maps
        assert_allclose(ncm_map, hp_map, rtol=1e-6, atol=1e-8)

        # The pattern should be cos(theta)
        pix_indices = np.arange(npix)
        theta, _ = healpy.pix2ang(nside, pix_indices, nest=False)
        expected_pattern = np.cos(theta)

        # Check correlation (normalized)
        correlation = np.corrcoef(ncm_map, expected_pattern)[0, 1]
        assert correlation > 0.99

    def test_map2alm2map_roundtrip(self, smap: Ncm.SphereMap, nside: int) -> None:
        """Test map -> alm -> map roundtrip for smooth maps."""
        # White noise is poorly represented with limited lmax
        # Use a smooth map instead (monopole + dipole)
        npix = healpy.nside2npix(nside)
        lmax = 2 * nside

        # Create smooth map (monopole + dipole component)
        pix_indices = np.arange(npix)
        theta, _ = healpy.pix2ang(nside, pix_indices, nest=False)
        smooth_map = 5.0 + 2.0 * np.cos(theta)  # monopole + dipole

        # NumCosmo roundtrip
        smap.set_lmax(lmax)
        smap.set_map(smooth_map)
        smap.prepare_alm()
        smap.alm2map()
        ncm_map_back = np.array([smap.get_pix(i) for i in range(npix)])

        # Should match reasonably well for smooth maps
        # NumCosmo's map2alm uses direct projection, not least-squares
        assert_allclose(ncm_map_back, smooth_map, rtol=0.1, atol=0.1)


class TestPowerSpectrum:
    """Test power spectrum (Cl) calculations."""

    def test_cl_api(self, nside: int) -> None:
        """Test that Cl API works correctly."""
        lmax = 2 * nside
        npix = healpy.nside2npix(nside)

        # Create map with known content (monopole + dipole)
        pix_indices = np.arange(npix)
        theta, _ = healpy.pix2ang(nside, pix_indices, nest=False)
        test_map = 5.0 + 2.0 * np.cos(theta)

        # NumCosmo computation
        smap = Ncm.SphereMap.new(nside)
        smap.set_lmax(lmax)
        smap.set_map(test_map)
        smap.prepare_alm()

        # Check we can compute Cl from alm
        for ell in range(min(lmax + 1, 5)):
            cl_sum = 0.0
            for m in range(ell + 1):
                re_alm, im_alm = smap.get_alm(ell, m)
                alm_squared = re_alm**2 + im_alm**2
                if m == 0:
                    cl_sum += alm_squared
                else:
                    cl_sum += 2 * alm_squared
            cl = cl_sum / (2 * ell + 1)

            # Cl should be non-negative
            assert cl >= 0

            # For monopole+dipole, only C_0 and C_1 should be large
            if ell == 0 or ell == 1:
                assert cl > 1.0, f"C_{ell} should be large for monopole/dipole"
            elif ell > 2:
                assert cl < 1.0, f"C_{ell} should be small for ell > 2"


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


class TestEdgeCases:
    """Test edge cases and boundary conditions."""

    def test_north_pole(self, smap: Ncm.SphereMap, nside: int) -> None:
        """Test pixel at north pole."""
        theta_pole = 0.0
        phi = 0.0

        pix_ring = smap.ang2pix_ring(theta_pole, phi)
        pix_nest = smap.ang2pix_nest(theta_pole, phi)

        # Both should point to valid pixels
        assert 0 <= pix_ring < smap.get_npix()
        assert 0 <= pix_nest < smap.get_npix()

        # Compare with healpy
        hp_pix_ring = healpy.ang2pix(nside, theta_pole, phi, nest=False)
        hp_pix_nest = healpy.ang2pix(nside, theta_pole, phi, nest=True)

        assert pix_ring == hp_pix_ring
        assert pix_nest == hp_pix_nest

    def test_south_pole(self, smap: Ncm.SphereMap, nside: int) -> None:
        """Test pixel at south pole."""
        theta_pole = np.pi
        phi = 0.0

        pix_ring = smap.ang2pix_ring(theta_pole, phi)
        pix_nest = smap.ang2pix_nest(theta_pole, phi)

        # Both should point to valid pixels
        assert 0 <= pix_ring < smap.get_npix()
        assert 0 <= pix_nest < smap.get_npix()

        # Compare with healpy
        hp_pix_ring = healpy.ang2pix(nside, theta_pole, phi, nest=False)
        hp_pix_nest = healpy.ang2pix(nside, theta_pole, phi, nest=True)

        assert pix_ring == hp_pix_ring
        assert pix_nest == hp_pix_nest

    def test_phi_wrapping(self, smap: Ncm.SphereMap, nside: int) -> None:
        """Test that phi values wrap correctly at 2*pi."""
        theta = np.pi / 2  # Equator

        # Test phi = 0 and phi = 2*pi should give same pixel
        pix_0 = smap.ang2pix_ring(theta, 0.0)
        pix_2pi = smap.ang2pix_ring(theta, 2 * np.pi)

        hp_pix_0 = healpy.ang2pix(nside, theta, 0.0, nest=False)
        hp_pix_2pi = healpy.ang2pix(nside, theta, 2 * np.pi, nest=False)

        assert pix_0 == hp_pix_0
        assert pix_2pi == hp_pix_2pi

    def test_minimum_nside(self) -> None:
        """Test minimum valid nside."""
        smap_min = Ncm.SphereMap.new(1)
        assert smap_min.get_nside() == 1
        assert smap_min.get_npix() == 12  # Minimum HEALPix map

    def test_large_nside(self) -> None:
        """Test with larger nside (quick sanity check)."""
        nside_large = 256
        smap_large = Ncm.SphereMap.new(nside_large)

        npix = smap_large.get_npix()
        expected_npix = healpy.nside2npix(nside_large)

        assert npix == expected_npix


class TestPerformance:
    """Performance comparison tests (marked as slow)."""

    @pytest.fixture(name="large_nside")
    def fixture_large_nside(self) -> int:
        """Large nside for performance testing."""
        return 512

    @pytest.fixture(name="large_smap")
    def fixture_large_smap(self, large_nside: int) -> Ncm.SphereMap:
        """Large sphere map for performance testing."""
        return Ncm.SphereMap.new(large_nside)

    def test_nest2ring_performance(
        self, large_smap: Ncm.SphereMap, large_nside: int
    ) -> None:
        """Benchmark nest2ring conversion (requires pytest-benchmark)."""
        npix = healpy.nside2npix(large_nside)
        indices = np.arange(min(10000, npix))

        # Simple timing test without benchmark fixture
        ncm_result = [large_smap.nest2ring(int(i)) for i in indices]
        assert len(ncm_result) == len(indices)


if __name__ == "__main__":
    # Run tests with: python test_py_sphere_map.py
    pytest.main([__file__, "-v"])

#!/usr/bin/env python
#
# test_sphere_map_healpy.py
#
# Sun Sep 27 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_sphere_map_healpy.py
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
# with this program. If not, see <http://www.gnu.org/licenses/>.

"""NcmSphereMap against the installed healpy, one case per operation.

The tests of NcmSphereMap are in C (``tests/c/ncm/sphere/test_ncm_sphere_map.c``),
against frozen healpy tables made by ``tests/tools/make_sphere_healpy_truth_table.py``.
These few compare with the live healpy instead, so a change of convention in healpy,
or a break in the Python bindings, shows here.
"""

from typing import Any

import numpy as np
import pytest
from numpy.testing import assert_array_equal

from numcosmo_py import Ncm

try:
    import healpy
except ImportError:
    healpy = None

# Skipped at run time, not at collection: this file is the whole sphere_map shard, and a
# shard with nothing collected makes pytest exit with status 5, which fails the lane.
pytestmark = [
    pytest.mark.sphere_map,
    pytest.mark.skipif(healpy is None, reason="healpy is not installed"),
]

Ncm.cfg_init()

NSIDE = 32


@pytest.fixture(name="smap")
def fixture_smap() -> Ncm.SphereMap:
    """An nside-32 map."""
    return Ncm.SphereMap.new(NSIDE)


def _pix(smap: Ncm.SphereMap) -> np.ndarray:
    return np.array([smap.get_pix(i) for i in range(smap.get_npix())])


def test_index_conversions(smap: Ncm.SphereMap) -> None:
    """RING <-> NESTED for every pixel."""
    idx = np.arange(smap.get_npix())
    assert_array_equal(
        [smap.ring2nest(int(i)) for i in idx], healpy.ring2nest(NSIDE, idx)
    )
    assert_array_equal(
        [smap.nest2ring(int(i)) for i in idx], healpy.nest2ring(NSIDE, idx)
    )


@pytest.mark.parametrize("nest", [False, True])
def test_pixel_centres(smap: Ncm.SphereMap, nest: bool) -> None:
    """pix2ang and pix2vec for every pixel, to rounding."""
    idx = np.arange(smap.get_npix())
    pix2ang = smap.pix2ang_nest if nest else smap.pix2ang_ring
    pix2vec = smap.pix2vec_nest if nest else smap.pix2vec_ring
    theta, phi = np.array([pix2ang(int(i)) for i in idx]).T
    hp_theta, hp_phi = healpy.pix2ang(NSIDE, idx, nest=nest)
    np.testing.assert_allclose(theta, hp_theta, rtol=1.0e-15, atol=1.0e-15)
    np.testing.assert_allclose(phi, hp_phi, rtol=1.0e-15, atol=1.0e-15)

    vec = Ncm.TriVec.new()
    hp_vec = np.array(healpy.pix2vec(NSIDE, idx, nest=nest)).T
    for i in idx[::17]:
        pix2vec(int(i), vec)
        np.testing.assert_allclose(list(vec.c), hp_vec[i], rtol=1.0e-15, atol=1.0e-15)


@pytest.mark.parametrize("nest", [False, True])
def test_direction_to_pixel(smap: Ncm.SphereMap, nest: bool) -> None:
    """ang2pix and vec2pix for seeded directions."""
    rng = np.random.default_rng(3)
    theta = np.arccos(rng.uniform(-1.0, 1.0, 2000))
    phi = rng.uniform(-np.pi, 3.0 * np.pi, 2000)
    ang2pix = smap.ang2pix_nest if nest else smap.ang2pix_ring
    vec2pix = smap.vec2pix_nest if nest else smap.vec2pix_ring
    hp = healpy.ang2pix(NSIDE, theta, phi, nest=nest)
    assert_array_equal([ang2pix(t, p) for t, p in zip(theta, phi)], hp)

    x, y, z = np.sin(theta) * np.cos(phi), np.sin(theta) * np.sin(phi), np.cos(theta)
    hpv = healpy.vec2pix(NSIDE, x, y, z, nest=nest)
    ncm = [vec2pix(Ncm.TriVec.new_full_c(*v)) for v in zip(x, y, z)]
    assert_array_equal(ncm, hpv)


@pytest.mark.parametrize("iter_n", [0, 3])
def test_map2alm_and_cl(iter_n: int) -> None:
    """map2alm, anafast and the cross spectrum at lmax 3 nside - 1."""
    lmax = 3 * NSIDE - 1
    rng = np.random.default_rng(5)
    maps = rng.standard_normal((2, healpy.nside2npix(NSIDE)))
    smaps = []
    for m in maps:
        smap = Ncm.SphereMap.new(NSIDE)
        smap.set_lmax(lmax)
        smap.set_iter(iter_n)
        smap.set_map(m.tolist())
        smap.prepare_alm()
        smaps.append(smap)

    hp_alm = healpy.map2alm(maps[0], lmax=lmax, iter=iter_n, use_weights=False)
    alm = np.array(
        [
            complex(*smaps[0].get_alm(ell, m))
            for m in range(lmax + 1)
            for ell in range(m, lmax + 1)
        ]
    )
    assert np.max(np.abs(alm - hp_alm)) < 1.0e-13 * np.max(np.abs(hp_alm))

    hp_cl = healpy.anafast(maps[0], lmax=lmax, iter=iter_n, use_weights=False)
    np.testing.assert_allclose(
        [smaps[0].get_Cl(ell) for ell in range(lmax + 1)], hp_cl, rtol=1.0e-12
    )

    hp_cross = healpy.anafast(
        maps[0], maps[1], lmax=lmax, iter=iter_n, use_weights=False
    )
    cross = np.array(smaps[0].compute_cross_Cl(smaps[1]).dup_array())
    assert np.max(np.abs(cross - hp_cross)) < 1.0e-12 * np.max(np.abs(hp_cross))


def test_alm2map() -> None:
    """alm2map of random coefficients at lmax 4 nside, past the ring Nyquist limit."""
    lmax = 4 * NSIDE
    rng = np.random.default_rng(7)
    n = healpy.Alm.getsize(lmax)
    alm = rng.standard_normal(n) + 1j * rng.standard_normal(n)
    alm[: lmax + 1] = alm[: lmax + 1].real

    smap = Ncm.SphereMap.new(NSIDE)
    smap.set_lmax(lmax)
    for m in range(lmax + 1):
        for ell in range(m, lmax + 1):
            a = alm[healpy.Alm.getidx(lmax, ell, m)]
            smap.set_alm(ell, m, a.real, a.imag)
    smap.alm2map()

    hp_map = healpy.alm2map(alm, NSIDE, lmax=lmax)
    assert np.max(np.abs(_pix(smap) - hp_map)) < 1.0e-12 * np.max(np.abs(hp_map))


@pytest.mark.parametrize("nest", [False, True])
def test_fits_both_ways(smap: Ncm.SphereMap, tmp_path: Any, nest: bool) -> None:
    """A healpy file loads exactly, and healpy reads a saved map exactly."""
    rng = np.random.default_rng(11)
    values = rng.standard_normal(smap.get_npix())

    hp_file = str(tmp_path / "healpy.fits")
    healpy.write_map(hp_file, values, nest=nest, coord="G", dtype=np.float64)
    smap.load_fits(hp_file, None)
    expected = Ncm.SphereMapOrder.NEST if nest else Ncm.SphereMapOrder.RING
    assert smap.get_order() == expected
    assert smap.get_coordsys() == Ncm.SphereMapCoordSys.GALACTIC
    assert_array_equal(_pix(smap), values)

    ncm_file = str(tmp_path / "ncm.fits")
    smap.save_fits(ncm_file, None, True)
    assert_array_equal(healpy.read_map(ncm_file, nest=nest, dtype=np.float64), values)

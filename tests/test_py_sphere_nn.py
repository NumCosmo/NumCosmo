#!/usr/bin/env python
#
# test_py_spline_cubic_d2.py
#
# Fri May 19 10:00:00 2023
# Copyright  2023  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_py_spline_cubic_d2.py
# Copyright (C) 2023 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

"""Tests on NcmSphereNN class."""

import pytest
import numpy as np
from numpy.testing import assert_allclose

from numcosmo_py import Ncm

Ncm.cfg_init()


@pytest.fixture(name="snn")
def fixture_snn():
    """Fixture for NcmSphereNN."""
    return Ncm.SphereNN()


@pytest.fixture(name="knn", params=[1, 2, 3, 4, 5])
def fixture_knn(request):
    """Fixture for NcmSphereNN."""
    return request.param


PREC = 1.0e-9


def test_setget(snn):
    """Test set and get methods."""
    r_a = np.random.uniform(0.1, 5.0, 1000)
    cos_theta_a = np.random.uniform(-1.0, 1.0, 1000)
    phi_a = np.random.uniform(-np.pi, np.pi, 1000)
    theta_a = np.acos(cos_theta_a)

    for r, theta, phi in zip(r_a, theta_a, phi_a):
        snn.insert(r, theta, phi)

    for i, coord in enumerate(zip(r_a, theta_a, phi_a)):
        assert_allclose(snn.get(i), coord, atol=0.0, rtol=PREC)


def test_insert_array(snn):
    """Test set_array method."""
    r_a = np.random.uniform(0.1, 5.0, 1000)
    cos_theta_a = np.random.uniform(-1.0, 1.0, 1000)
    phi_a = np.random.uniform(-np.pi, np.pi, 1000)
    theta_a = np.acos(cos_theta_a)

    snn.insert_array(r_a, theta_a, phi_a)

    for i, coord in enumerate(zip(r_a, theta_a, phi_a)):
        assert_allclose(snn.get(i), coord, atol=0.0, rtol=PREC)


def test_insert_array_large(snn):
    """Test set_array method."""
    r_a = np.random.uniform(0.1, 5.0, 70_000)
    cos_theta_a = np.random.uniform(-1.0, 1.0, 70_000)
    phi_a = np.random.uniform(-np.pi, np.pi, 70_000)
    theta_a = np.acos(cos_theta_a)

    snn.insert_array(r_a, theta_a, phi_a)

    for i, coord in enumerate(zip(r_a, theta_a, phi_a)):
        assert_allclose(snn.get(i), coord, atol=0.0, rtol=PREC)


def test_rebuild(snn):
    """Test rebuild method."""
    r_a = np.random.uniform(0.1, 5.0, 1000)
    cos_theta_a = np.random.uniform(-1.0, 1.0, 1000)
    phi_a = np.random.uniform(-np.pi, np.pi, 1000)
    theta_a = np.acos(cos_theta_a)

    for r, theta, phi in zip(r_a, theta_a, phi_a):
        snn.insert(r, theta, phi)

    assert snn.get_n() == 1000

    snn.rebuild()

    for i, coord in enumerate(zip(r_a, theta_a, phi_a)):
        assert_allclose(snn.get(i), coord, atol=0.0, rtol=PREC)

    new_r_a = np.random.uniform(0.1, 5.0, 1000)
    new_cos_theta_a = np.random.uniform(-1.0, 1.0, 1000)
    new_phi_a = np.random.uniform(-np.pi, np.pi, 1000)
    new_theta_a = np.acos(new_cos_theta_a)

    for r, theta, phi in zip(new_r_a, new_theta_a, new_phi_a):
        snn.insert(r, theta, phi)

    assert snn.get_n() == 2000

    r_a = np.concatenate((r_a, new_r_a))
    cos_theta_a = np.concatenate((cos_theta_a, new_cos_theta_a))
    phi_a = np.concatenate((phi_a, new_phi_a))
    theta_a = np.concatenate((theta_a, new_theta_a))

    snn.rebuild()

    for i, coord in enumerate(zip(r_a, theta_a, phi_a)):
        assert_allclose(snn.get(i), coord, atol=0.0, rtol=PREC)


def test_knn_search_1(snn):
    """Test knn_search method."""
    r_a = np.random.uniform(0.1, 5.0, 1000)
    cos_theta_a = np.random.uniform(-1.0, 1.0, 1000)
    phi_a = np.random.uniform(-np.pi, np.pi, 1000)
    theta_a = np.acos(cos_theta_a)

    for r, theta, phi in zip(r_a, theta_a, phi_a):
        snn.insert(r, theta, phi)

    snn.rebuild()

    for i, (r, theta, phi) in enumerate(zip(r_a, theta_a, phi_a)):
        idx = snn.knn_search(r, theta, phi, 1)
        assert len(idx) == 1
        assert idx.pop() == i

    new_r_a = np.random.uniform(0.1, 5.0, 1000)
    new_cos_theta_a = np.random.uniform(-1.0, 1.0, 1000)
    new_phi_a = np.random.uniform(-np.pi, np.pi, 1000)
    new_theta_a = np.acos(new_cos_theta_a)

    for r, theta, phi in zip(new_r_a, new_theta_a, new_phi_a):
        snn.insert(r, theta, phi)

    snn.rebuild()

    for i, (r, theta, phi) in enumerate(zip(new_r_a, new_theta_a, new_phi_a)):
        idx = snn.knn_search(r, theta, phi, 1)
        assert len(idx) == 1
        assert idx.pop() == i + 1000


def test_knn_search_k(snn, knn):
    """Test knn_search method."""
    r_a = np.random.uniform(0.1, 5.0, 1000)
    cos_theta_a = np.random.uniform(-1.0, 1.0, 1000)
    phi_a = np.random.uniform(-np.pi, np.pi, 1000)
    theta_a = np.acos(cos_theta_a)

    for r, theta, phi in zip(r_a, theta_a, phi_a):
        snn.insert(r, theta, phi)

    snn.rebuild()

    for i, (r, theta, phi) in enumerate(zip(r_a, theta_a, phi_a)):
        idx = snn.knn_search(r, theta, phi, knn)
        assert len(idx) == knn
        assert i in idx

    new_r_a = np.random.uniform(0.1, 5.0, 1000)
    new_cos_theta_a = np.random.uniform(-1.0, 1.0, 1000)
    new_phi_a = np.random.uniform(-np.pi, np.pi, 1000)
    new_theta_a = np.acos(new_cos_theta_a)

    for r, theta, phi in zip(new_r_a, new_theta_a, new_phi_a):
        snn.insert(r, theta, phi)

    snn.rebuild()

    for i, (r, theta, phi) in enumerate(zip(new_r_a, new_theta_a, new_phi_a)):
        idx = snn.knn_search(r, theta, phi, knn)
        assert len(idx) == knn
        assert i + 1000 in idx


def test_knn_search_vs_search_distances(snn, knn):
    """Test knn_search method."""
    r_a = np.random.uniform(0.1, 5.0, 1000)
    cos_theta_a = np.random.uniform(-1.0, 1.0, 1000)
    phi_a = np.random.uniform(-np.pi, np.pi, 1000)
    theta_a = np.acos(cos_theta_a)

    for r, theta, phi in zip(r_a, theta_a, phi_a):
        snn.insert(r, theta, phi)

    snn.rebuild()

    for r, theta, phi in zip(r_a, theta_a, phi_a):
        idx = snn.knn_search(r, theta, phi, knn)
        dist, idx2 = snn.knn_search_distances(r, theta, phi, knn)
        assert len(idx) == len(idx2)
        assert len(dist) == len(idx2)
        assert idx == idx2


def test_knn_search_distances(snn, knn):
    """Test knn_search_distances method."""
    r_a = np.random.uniform(0.1, 5.0, 1000)
    cos_theta_a = np.random.uniform(-1.0, 1.0, 1000)
    phi_a = np.random.uniform(-np.pi, np.pi, 1000)
    theta_a = np.acos(cos_theta_a)

    for r, theta, phi in zip(r_a, theta_a, phi_a):
        snn.insert(r, theta, phi)

    snn.rebuild()

    for i, (r, theta, phi) in enumerate(zip(r_a, theta_a, phi_a)):
        dist, idx = snn.knn_search_distances(r, theta, phi, knn)
        assert len(dist) == knn
        assert len(idx) == knn
        assert i in idx
        for d, j in zip(dist, idx):
            # Euclidean distance between j and theta, phi using spherical coordinates]
            x = r * np.sin(theta) * np.cos(phi)
            y = r * np.sin(theta) * np.sin(phi)
            z = r * np.cos(theta)
            x_j = r_a[j] * np.sin(theta_a[j]) * np.cos(phi_a[j])
            y_j = r_a[j] * np.sin(theta_a[j]) * np.sin(phi_a[j])
            z_j = r_a[j] * np.cos(theta_a[j])

            assert_allclose(
                d,
                (x - x_j) ** 2 + (y - y_j) ** 2 + (z - z_j) ** 2,
                atol=0.0,
                rtol=PREC,
            )

    new_r_a = np.random.uniform(0.1, 5.0, 1000)
    new_cos_theta_a = np.random.uniform(-1.0, 1.0, 1000)
    new_phi_a = np.random.uniform(-np.pi, np.pi, 1000)
    new_theta_a = np.acos(new_cos_theta_a)

    for r, theta, phi in zip(new_r_a, new_theta_a, new_phi_a):
        snn.insert(r, theta, phi)

    snn.rebuild()

    r_a = np.concatenate((r_a, new_r_a))
    cos_theta_a = np.concatenate((cos_theta_a, new_cos_theta_a))
    phi_a = np.concatenate((phi_a, new_phi_a))
    theta_a = np.concatenate((theta_a, new_theta_a))

    for i, (r, theta, phi) in enumerate(zip(new_r_a, new_theta_a, new_phi_a)):
        dist, idx = snn.knn_search_distances(r, theta, phi, knn)
        assert len(dist) == knn
        assert len(idx) == knn
        assert i + 1000 in idx
        for d, j in zip(dist, idx):
            # Euclidean distance between j and theta, phi using spherical coordinates]
            x = r * np.sin(theta) * np.cos(phi)
            y = r * np.sin(theta) * np.sin(phi)
            z = r * np.cos(theta)
            x_j = r_a[j] * np.sin(theta_a[j]) * np.cos(phi_a[j])
            y_j = r_a[j] * np.sin(theta_a[j]) * np.sin(phi_a[j])
            z_j = r_a[j] * np.cos(theta_a[j])
            assert_allclose(
                d,
                (x - x_j) ** 2 + (y - y_j) ** 2 + (z - z_j) ** 2,
                atol=0.0,
                rtol=PREC,
            )


def test_knn_search_distances_correct(snn, knn):
    """Test knn_search_distances method."""
    r_a = np.random.uniform(0.1, 5.0, 2000)
    cos_theta_a = np.random.uniform(-1.0, 1.0, 2000)
    phi_a = np.random.uniform(-np.pi, np.pi, 2000)
    theta_a = np.acos(cos_theta_a)
    sin_theta_a = np.sin(theta_a)
    cos_phi_a = np.cos(phi_a)
    sin_phi_a = np.sin(phi_a)
    knn = 15

    x_a = r_a * sin_theta_a * cos_phi_a
    y_a = r_a * sin_theta_a * sin_phi_a
    z_a = r_a * cos_theta_a

    for r, theta, phi in zip(r_a, theta_a, phi_a):
        snn.insert(r, theta, phi)

    snn.rebuild()

    for i, (r, theta, phi) in enumerate(zip(r_a, theta_a, phi_a)):
        d2_array = (x_a[i] - x_a) ** 2 + (y_a[i] - y_a) ** 2 + (z_a[i] - z_a) ** 2
        idx = np.argsort(d2_array)[:knn]
        dist, idx2 = snn.knn_search_distances(r, theta, phi, knn)
        assert_allclose(dist, d2_array[idx], atol=0.0, rtol=1e-11)
        assert set(idx) == set(idx2)


def test_knn_search_distances_correct_repeat(snn, knn):
    """Test knn_search_distances method."""
    nsize, knn, nrepeat = 1000, 15, 500
    r_a = np.random.uniform(0.1, 5.0, nsize)
    cos_theta_a = np.random.uniform(-1.0, 1.0, nsize)
    phi_a = np.random.uniform(-np.pi, np.pi, nsize)
    theta_a = np.acos(cos_theta_a)
    sin_theta_a = np.sin(theta_a)
    cos_phi_a = np.cos(phi_a)
    sin_phi_a = np.sin(phi_a)

    r_a = np.concatenate((r_a, r_a[:nrepeat]))
    cos_theta_a = np.concatenate((cos_theta_a, cos_theta_a[:nrepeat]))
    phi_a = np.concatenate((phi_a, phi_a[:nrepeat]))
    theta_a = np.concatenate((theta_a, theta_a[:nrepeat]))
    sin_theta_a = np.sin(theta_a)
    cos_phi_a = np.cos(phi_a)
    sin_phi_a = np.sin(phi_a)

    x_a = r_a * sin_theta_a * cos_phi_a
    y_a = r_a * sin_theta_a * sin_phi_a
    z_a = r_a * cos_theta_a

    for r, theta, phi in zip(r_a, theta_a, phi_a):
        snn.insert(r, theta, phi)
    snn.rebuild()

    for i, (r, theta, phi) in enumerate(zip(r_a, theta_a, phi_a)):
        d2_array = (x_a[i] - x_a) ** 2 + (y_a[i] - y_a) ** 2 + (z_a[i] - z_a) ** 2
        idx = np.argsort(d2_array)[:knn].tolist()
        dist, idx2 = snn.knn_search_distances(r, theta, phi, knn)

        assert_allclose(dist, d2_array[idx], atol=0.0, rtol=1e-11)
        assert set(np.array(idx) % nsize) == set(np.array(idx2) % nsize)


def test_knn_search_distances_correct_rng(snn, knn):
    """Test knn_search_distances method."""
    r_a = np.random.uniform(0.1, 5.0, 2000)
    cos_theta_a = np.random.uniform(-1.0, 1.0, 2000)
    phi_a = np.random.uniform(-np.pi, np.pi, 2000)
    theta_a = np.acos(cos_theta_a)
    sin_theta_a = np.sin(theta_a)
    cos_phi_a = np.cos(phi_a)
    sin_phi_a = np.sin(phi_a)
    knn = 15

    x_a = r_a * sin_theta_a * cos_phi_a
    y_a = r_a * sin_theta_a * sin_phi_a
    z_a = r_a * cos_theta_a

    for r, theta, phi in zip(r_a, theta_a, phi_a):
        snn.insert(r, theta, phi)

    snn.rebuild()

    for _ in range(100):
        r = np.random.uniform(0.1, 5.0)
        theta = np.random.uniform(0.0, np.pi)
        phi = np.random.uniform(-np.pi, np.pi)
        x = r * np.sin(theta) * np.cos(phi)
        y = r * np.sin(theta) * np.sin(phi)
        z = r * np.cos(theta)
        d2_array = (x - x_a) ** 2 + (y - y_a) ** 2 + (z - z_a) ** 2
        idx = np.argsort(d2_array)[:knn]
        dist, idx2 = snn.knn_search_distances(r, theta, phi, knn)
        assert_allclose(dist, d2_array[idx], atol=0.0, rtol=1e-11)
        assert set(idx) == set(idx2)


def test_knn_search_distances_batch(snn, knn):
    """Test knn_search_distances_batch method."""
    r_a = np.random.uniform(0.1, 5.0, 2000)
    cos_theta_a = np.random.uniform(-1.0, 1.0, 2000)
    phi_a = np.random.uniform(-np.pi, np.pi, 2000)
    theta_a = np.acos(cos_theta_a)

    m_r_a = np.random.uniform(0.1, 5.0, 2000)
    m_cos_theta_a = np.random.uniform(-1.0, 1.0, 2000)
    m_phi_a = np.random.uniform(-np.pi, np.pi, 2000)
    m_theta_a = np.acos(m_cos_theta_a)

    knn = 15

    for r, theta, phi in zip(r_a, theta_a, phi_a):
        snn.insert(r, theta, phi)

    snn.rebuild()

    dist_list, idx_list = snn.knn_search_distances_batch(m_r_a, m_theta_a, m_phi_a, knn)
    dist = np.array(dist_list).reshape(-1, knn)
    idx = np.array(idx_list).reshape(-1, knn)

    for m_r, m_theta, m_phi, m_dist, m_idx in zip(m_r_a, m_theta_a, m_phi_a, dist, idx):
        dist2, idx2 = snn.knn_search_distances(m_r, m_theta, m_phi, knn)
        assert_allclose(m_dist, dist2, atol=0.0, rtol=1e-11)
        assert set(m_idx) == set(idx2)


def test_dump(snn, capfd):
    """Test dump method."""
    r_a = np.random.uniform(0.1, 5.0, 1000)
    cos_theta_a = np.random.uniform(-1.0, 1.0, 1000)
    phi_a = np.random.uniform(-np.pi, np.pi, 1000)
    theta_a = np.acos(cos_theta_a)

    for r, theta, phi in zip(r_a, theta_a, phi_a):
        snn.insert(r, theta, phi)

    snn.rebuild()

    snn.dump_tree()

    out, _ = capfd.readouterr()
    assert "+-------" in out

#!/usr/bin/env python
#
# test_sky_match.py
#
# Mon Jan 15 10:19:40 2024
# Copyright  2023  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_py_serialize.py
# Copyright (C) 2024 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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
#

"""Tests for the sky_match module."""

from pathlib import Path
import pytest
import numpy as np
from numpy.testing import assert_allclose
from astropy.table import Table
from astropy.io import fits

from numcosmo_py import Nc, Ncm
from numcosmo_py.sky_match import NcSkyMatching

Ncm.cfg_init()

ATOL = 1.0e-7


@pytest.fixture(name="cosmo")
def fixture_cosmo():
    """Return a cosmology object."""
    return Nc.HICosmoDEXcdm()


@pytest.fixture(name="setup_catalogs")
def fixture_setup_catalogs(tmp_path) -> tuple[Path, Path]:
    """Create temporary FITS files for testing."""
    # Create temporary FITS files for testing
    query_catalog_path = tmp_path / "query_catalog.fits"
    match_catalog_path = tmp_path / "match_catalog.fits"
    RA_min = 10.0
    RA_max = 30.0
    DEC_min = -30.0
    DEC_max = -10.0
    # Create dummy data for query catalog (1000 objects)
    query_data = Table(
        {
            "RA_query": [np.random.uniform(RA_min, RA_max) for _ in range(1000)],
            "DEC_query": [np.random.uniform(DEC_min, DEC_max) for _ in range(1000)],
            "z_query": [np.random.uniform(0.9, 1.0) for _ in range(1000)],
            "ID_query": range(1000),
        }
    )
    query_data.write(query_catalog_path, format="fits")
    # Create dummy data for match catalog (1000 objects)
    match_data = Table(
        {
            "RA_match": [np.random.uniform(RA_min, RA_max) for _ in range(1000)],
            "DEC_match": [np.random.uniform(DEC_min, DEC_max) for _ in range(1000)],
            "z_match": [np.random.uniform(0.9, 1.0) for _ in range(1000)],
            "ID_match": range(1000),
        }
    )
    match_data.write(match_catalog_path, format="fits")

    return query_catalog_path, match_catalog_path


@pytest.fixture(name="setup_catalogs_extra_columns")
def fixture_setup_catalogs_extra_columns(tmp_path) -> tuple[Path, Path]:
    """Create temporary FITS files for testing."""
    # Create temporary FITS files for testing
    query_catalog_path = tmp_path / "query_catalog.fits"
    match_catalog_path = tmp_path / "match_catalog.fits"
    RA_min = 10.0
    RA_max = 30.0
    DEC_min = -30.0
    DEC_max = -10.0
    # Create dummy data for query catalog (1000 objects)
    query_data = Table(
        {
            "RA_query": [np.random.uniform(RA_min, RA_max) for _ in range(1000)],
            "DEC_query": [np.random.uniform(DEC_min, DEC_max) for _ in range(1000)],
            "z_query": [np.random.uniform(0.9, 1.0) for _ in range(1000)],
            "ID_query": range(1000),
            "extra_query": [np.random.uniform(0.0, 1.0) for _ in range(1000)],
        }
    )
    query_data.write(query_catalog_path, format="fits")
    # Create dummy data for match catalog (1000 objects)
    match_data = Table(
        {
            "RA_match": [np.random.uniform(RA_min, RA_max) for _ in range(1000)],
            "DEC_match": [np.random.uniform(DEC_min, DEC_max) for _ in range(1000)],
            "z_match": [np.random.uniform(0.9, 1.0) for _ in range(1000)],
            "ID_match": range(1000),
            "extra_match": [np.random.uniform(0.0, 1.0) for _ in range(1000)],
        }
    )
    match_data.write(match_catalog_path, format="fits")

    return query_catalog_path, match_catalog_path


@pytest.fixture(name="setup_catalogs_nonexistent")
def fixture_setup_catalogs_nonexistent(tmp_path) -> tuple[Path, Path]:
    """Create temporary FITS files for testing."""
    # Create temporary FITS files for testing
    query_catalog_path = tmp_path / "query_catalog.fits"
    match_catalog_path = tmp_path / "match_catalog.fits"
    return query_catalog_path, match_catalog_path


@pytest.fixture(name="setup_catalogs_fits_containing_image")
def fixture_setup_catalogs_fits_containing_image(tmp_path) -> tuple[Path, Path]:
    """Create temporary FITS files for testing."""
    # Create temporary FITS files for testing
    query_catalog_path = tmp_path / "query_catalog.fits"
    match_catalog_path = tmp_path / "match_catalog.fits"
    fits_data = np.random.rand(100, 100)
    hdu = fits.PrimaryHDU(fits_data)
    hdul = fits.HDUList([hdu])
    hdul.writeto(query_catalog_path)
    hdul.writeto(match_catalog_path)

    return query_catalog_path, match_catalog_path


def test_load_fits_data(setup_catalogs):
    """Test the load_fits_data function."""
    query_catalog_path, match_catalog_path = setup_catalogs
    matching = NcSkyMatching(
        query_catalog_path=query_catalog_path,
        query_coordinates={"RA": "RA_query", "DEC": "DEC_query", "z": "z_query"},
        match_catalog_path=match_catalog_path,
        match_coordinates={"RA": "RA_match", "DEC": "DEC_match", "z": "z_match"},
    )
    assert matching is not None
    assert matching.query_data is not None
    assert matching.match_data is not None
    assert len(matching.query_data) == 1000
    assert len(matching.match_data) == 1000


def test_load_fits_data_wrong_query_map(setup_catalogs):
    """Test the load_fits_data function with wrong query map."""
    query_catalog_path, match_catalog_path = setup_catalogs
    with pytest.raises(ValueError, match="RA and DEC coordinates mapped by.*"):
        _ = NcSkyMatching(
            query_catalog_path=query_catalog_path,
            query_coordinates={
                "RA": "RA_bob",
                "DEC": "DEC_fred",
                "z": "z_wilma",
                "ID": "ID_dunga",
            },
            match_catalog_path=match_catalog_path,
            match_coordinates={"RA": "RA_match", "DEC": "DEC_match", "z": "z_match"},
        )


def test_load_fits_data_wrong_match_map(setup_catalogs):
    """Test the load_fits_data function with wrong match map."""
    query_catalog_path, match_catalog_path = setup_catalogs
    with pytest.raises(ValueError, match="RA and DEC coordinates mapped by.*"):
        _ = NcSkyMatching(
            query_catalog_path=query_catalog_path,
            query_coordinates={"RA": "RA_query", "DEC": "DEC_query", "z": "z_query"},
            match_catalog_path=match_catalog_path,
            match_coordinates={
                "RA": "RA_bob",
                "DEC": "DEC_fred",
                "z": "z_wilma",
                "ID": "ID_dunga",
            },
        )


def test_load_fits_data_nonexistent(setup_catalogs_nonexistent):
    """Test the load_fits_data function with nonexistent catalog."""
    query_catalog_path, match_catalog_path = setup_catalogs_nonexistent
    with pytest.raises(FileNotFoundError):
        _ = NcSkyMatching(
            query_catalog_path=query_catalog_path,
            query_coordinates={"RA": "RA_query", "DEC": "DEC_query", "z": "z_query"},
            match_catalog_path=match_catalog_path,
            match_coordinates={"RA": "RA_match", "DEC": "DEC_match", "z": "z_match"},
        )


def test_load_fits_data_fits_containing_image(setup_catalogs_fits_containing_image):
    """Test the load_fits_data function with FITS containing image."""
    query_catalog_path, match_catalog_path = setup_catalogs_fits_containing_image
    with pytest.raises(
        ValueError, match="No FITS table found in the provided catalog."
    ):
        _ = NcSkyMatching(
            query_catalog_path=query_catalog_path,
            query_coordinates={"RA": "RA_query", "DEC": "DEC_query", "z": "z_query"},
            match_catalog_path=match_catalog_path,
            match_coordinates={"RA": "RA_match", "DEC": "DEC_match", "z": "z_match"},
        )


def test_match_2d(setup_catalogs):
    """Test the match_2d function."""
    query_catalog_path, match_catalog_path = setup_catalogs
    matching = NcSkyMatching(
        query_catalog_path=query_catalog_path,
        query_coordinates={"RA": "RA_query", "DEC": "DEC_query", "z": "z_query"},
        match_catalog_path=match_catalog_path,
        match_coordinates={"RA": "RA_match", "DEC": "DEC_match", "z": "z_match"},
    )
    mt = matching.match_2d(
        matching_distance=np.pi / 2.0, n_nearest_neighbours=10, verbose=False
    )
    assert mt is not None

    for key in [
        "RA",
        "DEC",
        "ID",
        "RA_matched",
        "DEC_matched",
        "ID_matched",
    ]:
        assert key in mt.keys()

    for row in mt:
        row_id = row["ID"]
        qrow = matching.query_data[row_id]
        RA_query = qrow["RA_query"]
        DEC_query = qrow["DEC_query"]
        assert RA_query == row["RA"]
        assert DEC_query == row["DEC"]
        for i, (matched_id, distance) in enumerate(
            zip(row["ID_matched"], row["distances"])
        ):
            RA_match = matching.match_data[matched_id]["RA_match"]
            DEC_match = matching.match_data[matched_id]["DEC_match"]
            assert RA_match == row["RA_matched"][i]
            assert DEC_match == row["DEC_matched"][i]
            # Convert to theta, phi
            theta_match = np.pi / 2.0 - np.radians(DEC_match)
            phi_match = np.radians(RA_match)
            theta_query = np.pi / 2.0 - np.radians(DEC_query)
            phi_query = np.radians(RA_query)
            # Calculate angular distance
            # v = (sin(theta)cos(phi), sin(theta)sin(phi), cos(theta))
            # angular_distance = arccos(v1.v2)
            angular_distance = np.arccos(
                np.sin(theta_query)
                * np.sin(theta_match)
                * np.cos(phi_query - phi_match)
                + np.cos(theta_query) * np.cos(theta_match)
            )
            assert angular_distance <= np.pi
            assert_allclose(angular_distance, distance, atol=ATOL)


def test_match_3d(cosmo, setup_catalogs):
    """Test the match_3d function."""
    query_catalog_path, match_catalog_path = setup_catalogs
    matching = NcSkyMatching(
        query_catalog_path=query_catalog_path,
        query_coordinates={"RA": "RA_query", "DEC": "DEC_query", "z": "z_query"},
        match_catalog_path=match_catalog_path,
        match_coordinates={"RA": "RA_match", "DEC": "DEC_match", "z": "z_match"},
    )
    mt = matching.match_3d(
        cosmo, matching_distance=cosmo.RH_Mpc(), n_nearest_neighbours=10, verbose=False
    )
    assert mt is not None
    dist = Nc.Distance(zf=10.0)
    dist.prepare(cosmo)

    for key in [
        "RA",
        "DEC",
        "z",
        "ID",
        "RA_matched",
        "DEC_matched",
        "z_matched",
        "ID_matched",
    ]:
        assert key in mt.keys()

    for row in mt:
        row_id = row["ID"]
        qrow = matching.query_data[row_id]
        RA_query = qrow["RA_query"]
        DEC_query = qrow["DEC_query"]
        z_query = qrow["z_query"]
        assert RA_query == row["RA"]
        assert DEC_query == row["DEC"]
        assert z_query == row["z"]
        for i, (matched_id, distance) in enumerate(
            zip(row["ID_matched"], row["distances Mpc"])
        ):
            RA_match = matching.match_data[matched_id]["RA_match"]
            DEC_match = matching.match_data[matched_id]["DEC_match"]
            z_match = matching.match_data[matched_id]["z_match"]
            assert RA_match == row["RA_matched"][i]
            assert DEC_match == row["DEC_matched"][i]
            assert z_match == row["z_matched"][i]
            # Convert to theta, phi
            theta_match = np.pi / 2.0 - np.radians(DEC_match)
            phi_match = np.radians(RA_match)
            theta_query = np.pi / 2.0 - np.radians(DEC_query)
            phi_query = np.radians(RA_query)
            # Calculate angular distance
            # v_query = r(sin(theta)cos(phi), sin(theta)sin(phi), cos(theta))
            # v_match = r(sin(theta)cos(phi), sin(theta)sin(phi), cos(theta))
            # Euclidean distance = sqrt((x2-x1)^2 + (y2-y1)^2 + (z2-z1)^2)
            # Calculate euclidean distance
            distance_query = dist.comoving(cosmo, z_query) * cosmo.RH_Mpc()
            distance_match = dist.comoving(cosmo, z_match) * cosmo.RH_Mpc()
            x1_query, x2_query, x3_query = (
                distance_query * np.sin(theta_query) * np.cos(phi_query),
                distance_query * np.sin(theta_query) * np.sin(phi_query),
                distance_query * np.cos(theta_query),
            )
            x1_match, x2_match, x3_match = (
                distance_match * np.sin(theta_match) * np.cos(phi_match),
                distance_match * np.sin(theta_match) * np.sin(phi_match),
                distance_match * np.cos(theta_match),
            )
            euclidean_distance = np.sqrt(
                (x1_match - x1_query) ** 2
                + (x2_match - x2_query) ** 2
                + (x3_match - x3_query) ** 2
            )
            assert_allclose(euclidean_distance, distance, atol=ATOL)


def test_match_2d_extra_columns(setup_catalogs_extra_columns):
    """Test the match_2d function with extra columns."""
    query_catalog_path, match_catalog_path = setup_catalogs_extra_columns
    matching = NcSkyMatching(
        query_catalog_path=query_catalog_path,
        query_coordinates={"RA": "RA_query", "DEC": "DEC_query", "z": "z_query"},
        match_catalog_path=match_catalog_path,
        match_coordinates={"RA": "RA_match", "DEC": "DEC_match", "z": "z_match"},
        query_properties=["extra_query"],
        match_properties=["extra_match"],
    )
    mt = matching.match_2d(
        matching_distance=np.pi / 2.0, n_nearest_neighbours=10, verbose=False
    )
    assert mt is not None

    for key in [
        "RA",
        "DEC",
        "ID",
        "RA_matched",
        "DEC_matched",
        "ID_matched",
        "extra_query",
        "extra_match",
    ]:
        assert key in mt.keys()

    for row in mt:
        row_id = row["ID"]
        qrow = matching.query_data[row_id]
        RA_query = qrow["RA_query"]
        DEC_query = qrow["DEC_query"]
        assert RA_query == row["RA"]
        assert DEC_query == row["DEC"]
        assert qrow["extra_query"] == row["extra_query"]
        for i, (matched_id, distance) in enumerate(
            zip(row["ID_matched"], row["distances"])
        ):
            RA_match = matching.match_data[matched_id]["RA_match"]
            DEC_match = matching.match_data[matched_id]["DEC_match"]
            assert RA_match == row["RA_matched"][i]
            assert DEC_match == row["DEC_matched"][i]
            assert (
                matching.match_data[matched_id]["extra_match"] == row["extra_match"][i]
            )
            # Convert to theta, phi
            theta_match = np.pi / 2.0 - np.radians(DEC_match)
            phi_match = np.radians(RA_match)
            theta_query = np.pi / 2.0 - np.radians(DEC_query)
            phi_query = np.radians(RA_query)
            # Calculate angular distance
            # v = (sin(theta)cos(phi), sin(theta)sin(phi), cos(theta))
            # angular_distance = arccos(v1.v2)
            angular_distance = np.arccos(
                np.sin(theta_query)
                * np.sin(theta_match)
                * np.cos(phi_query - phi_match)
                + np.cos(theta_query) * np.cos(theta_match)
            )
            assert angular_distance <= np.pi
            assert_allclose(angular_distance, distance, atol=ATOL)


def test_match_3d_extra_columns(cosmo, setup_catalogs_extra_columns):
    """Test the match_3d function with extra columns."""
    query_catalog_path, match_catalog_path = setup_catalogs_extra_columns
    matching = NcSkyMatching(
        query_catalog_path=query_catalog_path,
        query_coordinates={"RA": "RA_query", "DEC": "DEC_query", "z": "z_query"},
        match_catalog_path=match_catalog_path,
        match_coordinates={"RA": "RA_match", "DEC": "DEC_match", "z": "z_match"},
        query_properties=["extra_query"],
        match_properties=["extra_match"],
    )
    mt = matching.match_3d(
        cosmo, matching_distance=cosmo.RH_Mpc(), n_nearest_neighbours=10, verbose=False
    )
    assert mt is not None
    dist = Nc.Distance(zf=10.0)
    dist.prepare(cosmo)

    for key in [
        "RA",
        "DEC",
        "z",
        "ID",
        "RA_matched",
        "DEC_matched",
        "z_matched",
        "ID_matched",
        "extra_query",
        "extra_match",
    ]:
        assert key in mt.keys()

    for row in mt:
        row_id = row["ID"]
        qrow = matching.query_data[row_id]
        RA_query = qrow["RA_query"]
        DEC_query = qrow["DEC_query"]
        z_query = qrow["z_query"]
        assert RA_query == row["RA"]
        assert DEC_query == row["DEC"]
        assert z_query == row["z"]
        assert qrow["extra_query"] == row["extra_query"]
        for i, (matched_id, distance) in enumerate(
            zip(row["ID_matched"], row["distances Mpc"])
        ):
            RA_match = matching.match_data[matched_id]["RA_match"]
            DEC_match = matching.match_data[matched_id]["DEC_match"]
            z_match = matching.match_data[matched_id]["z_match"]
            assert RA_match == row["RA_matched"][i]
            assert DEC_match == row["DEC_matched"][i]
            assert z_match == row["z_matched"][i]
            assert (
                matching.match_data[matched_id]["extra_match"] == row["extra_match"][i]
            )
            # Convert to theta, phi
            theta_match = np.pi / 2.0 - np.radians(DEC_match)
            phi_match = np.radians(RA_match)
            theta_query = np.pi / 2.0 - np.radians(DEC_query)
            phi_query = np.radians(RA_query)
            # Calculate angular distance
            # v_query = r(sin(theta)cos(phi), sin(theta)sin(phi), cos(theta))
            # v_match = r(sin(theta)cos(phi), sin(theta)sin(phi), cos(theta))
            # Euclidean distance = sqrt((x2-x1)^2 + (y2-y1)^2 + (z2-z1)^2)
            # Calculate euclidean distance
            distance_query = dist.comoving(cosmo, z_query) * cosmo.RH_Mpc()
            distance_match = dist.comoving(cosmo, z_match) * cosmo.RH_Mpc()
            x1_query, x2_query, x3_query = (
                distance_query * np.sin(theta_query) * np.cos(phi_query),
                distance_query * np.sin(theta_query) * np.sin(phi_query),
                distance_query * np.cos(theta_query),
            )
            x1_match, x2_match, x3_match = (
                distance_match * np.sin(theta_match) * np.cos(phi_match),
                distance_match * np.sin(theta_match) * np.sin(phi_match),
                distance_match * np.cos(theta_match),
            )
            euclidean_distance = np.sqrt(
                (x1_match - x1_query) ** 2
                + (x2_match - x2_query) ** 2
                + (x3_match - x3_query) ** 2
            )
            assert_allclose(euclidean_distance, distance, atol=ATOL)

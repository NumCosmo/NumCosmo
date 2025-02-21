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
from numcosmo_py.sky_match import SkyMatch, DistanceMethod, Mask

Ncm.cfg_init()

ATOL = 1.0e-7
QUERY_SIZE = 1200
MATCH_SIZE = 1000


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
    # Create dummy data for query catalog (QUERY_SIZE objects)
    query_data = Table(
        {
            "RA_query": [np.random.uniform(RA_min, RA_max) for _ in range(QUERY_SIZE)],
            "DEC_query": [
                np.random.uniform(DEC_min, DEC_max) for _ in range(QUERY_SIZE)
            ],
            "z_query": [np.random.uniform(0.9, 1.0) for _ in range(QUERY_SIZE)],
            "ID_query": range(QUERY_SIZE),
        }
    )
    query_data.write(query_catalog_path, format="fits")
    # Create dummy data for match catalog (MATCH_SIZE objects)
    match_data = Table(
        {
            "RA_match": [np.random.uniform(RA_min, RA_max) for _ in range(MATCH_SIZE)],
            "DEC_match": [
                np.random.uniform(DEC_min, DEC_max) for _ in range(MATCH_SIZE)
            ],
            "z_match": [np.random.uniform(0.9, 1.0) for _ in range(MATCH_SIZE)],
            "ID_match": range(MATCH_SIZE),
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
    # Create dummy data for query catalog (QUERY_SIZE objects)
    query_data = Table(
        {
            "RA_query": [np.random.uniform(RA_min, RA_max) for _ in range(QUERY_SIZE)],
            "DEC_query": [
                np.random.uniform(DEC_min, DEC_max) for _ in range(QUERY_SIZE)
            ],
            "z_query": [np.random.uniform(0.9, 1.0) for _ in range(QUERY_SIZE)],
            "ID_query": range(QUERY_SIZE),
            "extra_query": [np.random.uniform(0.0, 1.0) for _ in range(QUERY_SIZE)],
        }
    )
    query_data.write(query_catalog_path, format="fits")
    # Create dummy data for match catalog (MATCH_SIZE objects)
    match_data = Table(
        {
            "RA_match": [np.random.uniform(RA_min, RA_max) for _ in range(MATCH_SIZE)],
            "DEC_match": [
                np.random.uniform(DEC_min, DEC_max) for _ in range(MATCH_SIZE)
            ],
            "z_match": [np.random.uniform(0.9, 1.0) for _ in range(MATCH_SIZE)],
            "ID_match": range(MATCH_SIZE),
            "extra_match": [np.random.uniform(0.0, 1.0) for _ in range(MATCH_SIZE)],
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


@pytest.fixture(name="distance_method", params=list(DistanceMethod))
def fixture_distance_method(request):
    """Return a distance method."""
    return request.param


def test_load_fits_data(setup_catalogs):
    """Test the load_fits_data function."""
    query_catalog_path, match_catalog_path = setup_catalogs
    matching = SkyMatch.new_from_fits(
        query_catalog_path=query_catalog_path,
        query_coordinates={"RA": "RA_query", "DEC": "DEC_query", "z": "z_query"},
        match_catalog_path=match_catalog_path,
        match_coordinates={"RA": "RA_match", "DEC": "DEC_match", "z": "z_match"},
    )
    assert matching is not None
    assert matching.query_data is not None
    assert matching.match_data is not None
    assert len(matching.query_data) == QUERY_SIZE
    assert len(matching.match_data) == MATCH_SIZE


def test_load_fits_data_wrong_query_map(setup_catalogs):
    """Test the load_fits_data function with wrong query map."""
    query_catalog_path, match_catalog_path = setup_catalogs
    with pytest.raises(ValueError, match="RA and DEC coordinates mapped by.*"):
        _ = SkyMatch.new_from_fits(
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
        _ = SkyMatch.new_from_fits(
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
        _ = SkyMatch.new_from_fits(
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
        _ = SkyMatch.new_from_fits(
            query_catalog_path=query_catalog_path,
            query_coordinates={"RA": "RA_query", "DEC": "DEC_query", "z": "z_query"},
            match_catalog_path=match_catalog_path,
            match_coordinates={"RA": "RA_match", "DEC": "DEC_match", "z": "z_match"},
        )


def test_load_fits_data_missing_RA_query(setup_catalogs):
    """Test the load_fits_data function with missing RA_query column."""
    query_catalog_path, match_catalog_path = setup_catalogs
    with pytest.raises(ValueError, match="RA and DEC coordinates must be provided."):
        _ = SkyMatch.new_from_fits(
            query_catalog_path=query_catalog_path,
            query_coordinates={"DEC": "DEC_query", "z": "z_query"},
            match_catalog_path=match_catalog_path,
            match_coordinates={"RA": "RA_match", "DEC": "DEC_match", "z": "z_match"},
        )


def test_load_fits_data_missing_DEC_query(setup_catalogs):
    """Test the load_fits_data function with missing DEC_query column."""
    query_catalog_path, match_catalog_path = setup_catalogs
    with pytest.raises(ValueError, match="RA and DEC coordinates must be provided."):
        _ = SkyMatch.new_from_fits(
            query_catalog_path=query_catalog_path,
            query_coordinates={"RA": "RA_query", "z": "z_query"},
            match_catalog_path=match_catalog_path,
            match_coordinates={"RA": "RA_match", "DEC": "DEC_match", "z": "z_match"},
        )


def test_load_fits_data_missing_RA_match(setup_catalogs):
    """Test the load_fits_data function with missing RA_match column."""
    query_catalog_path, match_catalog_path = setup_catalogs
    with pytest.raises(ValueError, match="RA and DEC coordinates must be provided."):
        _ = SkyMatch.new_from_fits(
            query_catalog_path=query_catalog_path,
            query_coordinates={"RA": "RA_query", "DEC": "DEC_query", "z": "z_query"},
            match_catalog_path=match_catalog_path,
            match_coordinates={"DEC": "DEC_match", "z": "z_match"},
        )


def test_load_fits_data_missing_DEC_match(setup_catalogs):
    """Test the load_fits_data function with missing DEC_match column."""
    query_catalog_path, match_catalog_path = setup_catalogs
    with pytest.raises(ValueError, match="RA and DEC coordinates must be provided."):
        _ = SkyMatch.new_from_fits(
            query_catalog_path=query_catalog_path,
            query_coordinates={"RA": "RA_query", "DEC": "DEC_query", "z": "z_query"},
            match_catalog_path=match_catalog_path,
            match_coordinates={"RA": "RA_match", "z": "z_match"},
        )


def test_load_fits_data_wrong_RA_query(setup_catalogs):
    """Test the load_fits_data function with wrong RA_query column."""
    query_catalog_path, match_catalog_path = setup_catalogs
    with pytest.raises(
        ValueError,
        match=(
            "RA and DEC coordinates mapped by .* "
            "not found in the provided catalog .*"
        ),
    ):
        _ = SkyMatch.new_from_fits(
            query_catalog_path=query_catalog_path,
            query_coordinates={
                "RA": "RA_query_wrong",
                "DEC": "DEC_query",
                "z": "z_query",
            },
            match_catalog_path=match_catalog_path,
            match_coordinates={"RA": "RA_match", "DEC": "DEC_match", "z": "z_match"},
        )


def test_load_fits_data_wrong_DEC_query(setup_catalogs):
    """Test the load_fits_data function with wrong DEC_query column."""
    query_catalog_path, match_catalog_path = setup_catalogs
    with pytest.raises(
        ValueError,
        match=(
            "RA and DEC coordinates mapped by .* "
            "not found in the provided catalog .*"
        ),
    ):
        _ = SkyMatch.new_from_fits(
            query_catalog_path=query_catalog_path,
            query_coordinates={
                "RA": "RA_query",
                "DEC": "DEC_query_wrong",
                "z": "z_query",
            },
            match_catalog_path=match_catalog_path,
            match_coordinates={"RA": "RA_match", "DEC": "DEC_match", "z": "z_match"},
        )


def test_load_fits_data_wrong_RA_match(setup_catalogs):
    """Test the load_fits_data function with wrong RA_match column."""
    query_catalog_path, match_catalog_path = setup_catalogs
    with pytest.raises(
        ValueError,
        match=(
            "RA and DEC coordinates mapped by .* "
            "not found in the provided catalog .*"
        ),
    ):
        _ = SkyMatch.new_from_fits(
            query_catalog_path=query_catalog_path,
            query_coordinates={"RA": "RA_query", "DEC": "DEC_query", "z": "z_query"},
            match_catalog_path=match_catalog_path,
            match_coordinates={
                "RA": "RA_match_wrong",
                "DEC": "DEC_match",
                "z": "z_match",
            },
        )


def test_load_fits_data_wrong_DEC_match(setup_catalogs):
    """Test the load_fits_data function with wrong DEC_match column."""
    query_catalog_path, match_catalog_path = setup_catalogs
    with pytest.raises(
        ValueError,
        match=(
            "RA and DEC coordinates mapped by .* "
            "not found in the provided catalog .*"
        ),
    ):
        _ = SkyMatch.new_from_fits(
            query_catalog_path=query_catalog_path,
            query_coordinates={"RA": "RA_query", "DEC": "DEC_query", "z": "z_query"},
            match_catalog_path=match_catalog_path,
            match_coordinates={
                "RA": "RA_match",
                "DEC": "DEC_match_wrong",
                "z": "z_match",
            },
        )


def test_load_fits_data_wrong_z_query(setup_catalogs):
    """Test the load_fits_data function with wrong z_query column."""
    query_catalog_path, match_catalog_path = setup_catalogs
    with pytest.raises(
        ValueError,
        match="Redshift coordinate mapped by .* not found in the provided catalog .*",
    ):
        _ = SkyMatch.new_from_fits(
            query_catalog_path=query_catalog_path,
            query_coordinates={
                "RA": "RA_query",
                "DEC": "DEC_query",
                "z": "z_query_wrong",
            },
            match_catalog_path=match_catalog_path,
            match_coordinates={"RA": "RA_match", "DEC": "DEC_match", "z": "z_match"},
        )


def test_load_fits_data_wrong_z_match(setup_catalogs):
    """Test the load_fits_data function with wrong z_match column."""
    query_catalog_path, match_catalog_path = setup_catalogs
    with pytest.raises(
        ValueError,
        match="Redshift coordinate mapped by .* not found in the provided catalog .*",
    ):
        _ = SkyMatch.new_from_fits(
            query_catalog_path=query_catalog_path,
            query_coordinates={"RA": "RA_query", "DEC": "DEC_query", "z": "z_query"},
            match_catalog_path=match_catalog_path,
            match_coordinates={
                "RA": "RA_match",
                "DEC": "DEC_match",
                "z": "z_match_wrong",
            },
        )


def test_match_2d(cosmo, setup_catalogs, distance_method):
    """Test the match_2d function."""
    query_catalog_path, match_catalog_path = setup_catalogs
    matching = SkyMatch.new_from_fits(
        query_catalog_path=query_catalog_path,
        query_coordinates={"RA": "RA_query", "DEC": "DEC_query", "z": "z_query"},
        match_catalog_path=match_catalog_path,
        match_coordinates={"RA": "RA_match", "DEC": "DEC_match", "z": "z_match"},
    )
    result = matching.match_2d(
        cosmo,
        n_nearest_neighbours=10,
        distance_method=distance_method,
    )
    assert result is not None

    dist = Nc.Distance(zf=10.0)
    dist.prepare(cosmo)
    RH_Mpc = cosmo.RH_Mpc()

    match_ra = matching.match_ra
    match_dec = matching.match_dec
    match_theta = np.pi / 2.0 - np.radians(match_dec)
    match_phi = np.radians(match_ra)

    query_ra = matching.query_ra
    query_dec = matching.query_dec
    query_theta = np.pi / 2.0 - np.radians(query_dec)
    query_phi = np.radians(query_ra)

    match_theta_2d = match_theta[result.nearest_neighbours_indices]
    match_phi_2d = match_phi[result.nearest_neighbours_indices]

    query_theta_2d = np.repeat(query_theta, 10).reshape(-1, 10)
    query_phi_2d = np.repeat(query_phi, 10).reshape(-1, 10)

    angular_distance = np.arccos(
        np.sin(query_theta_2d)
        * np.sin(match_theta_2d)
        * np.cos(query_phi_2d - match_phi_2d)
        + np.cos(query_theta_2d) * np.cos(match_theta_2d)
    )
    match distance_method:
        case DistanceMethod.ANGULAR_SEPARATION:
            assert_allclose(
                angular_distance, result.nearest_neighbours_distances, atol=ATOL
            )
        case DistanceMethod.QUERY_RADIUS:
            query_r = (
                np.array(dist.angular_diameter_vector(cosmo, matching.query_z)) * RH_Mpc
            )
            distances = query_r.reshape(-1, 1) * 2.0 * np.sin(angular_distance / 2.0)
            assert_allclose(distances, result.nearest_neighbours_distances)
        case DistanceMethod.MATCH_RADIUS:
            match_r = (
                np.array(dist.angular_diameter_vector(cosmo, matching.match_z)) * RH_Mpc
            )
            distances = (
                match_r[result.nearest_neighbours_indices]
                * 2.0
                * np.sin(angular_distance / 2.0)
            )
            assert_allclose(distances, result.nearest_neighbours_distances)
        case DistanceMethod.MIN_RADIUS:
            query_r = (
                np.array(dist.angular_diameter_vector(cosmo, matching.query_z)) * RH_Mpc
            )
            match_r = (
                np.array(dist.angular_diameter_vector(cosmo, matching.match_z)) * RH_Mpc
            )
            distances = (
                np.minimum(
                    query_r.reshape(-1, 1), match_r[result.nearest_neighbours_indices]
                )
                * 2.0
                * np.sin(angular_distance / 2.0)
            )
            assert_allclose(distances, result.nearest_neighbours_distances)
        case DistanceMethod.MAX_RADIUS:
            query_r = (
                np.array(dist.angular_diameter_vector(cosmo, matching.query_z)) * RH_Mpc
            )
            match_r = (
                np.array(dist.angular_diameter_vector(cosmo, matching.match_z)) * RH_Mpc
            )
            distances = (
                np.maximum(
                    query_r.reshape(-1, 1), match_r[result.nearest_neighbours_indices]
                )
                * 2.0
                * np.sin(angular_distance / 2.0)
            )
            assert_allclose(distances, result.nearest_neighbours_distances)
        case _:
            raise ValueError(f"Invalid distance method: {distance_method}")


def test_match_2d_table(cosmo, setup_catalogs, distance_method):
    """Test the match_2d function."""
    query_catalog_path, match_catalog_path = setup_catalogs
    matching = SkyMatch.new_from_fits(
        query_catalog_path=query_catalog_path,
        query_coordinates={"RA": "RA_query", "DEC": "DEC_query", "z": "z_query"},
        match_catalog_path=match_catalog_path,
        match_coordinates={"RA": "RA_match", "DEC": "DEC_match", "z": "z_match"},
    )
    result = matching.match_2d(
        cosmo, n_nearest_neighbours=10, distance_method=distance_method
    )
    assert result is not None

    mask_array = np.random.choice(
        [True, False], size=10 * len(matching.query_ra)
    ).reshape(-1, 10)
    table = result.to_table_complete(mask=Mask(mask_array))

    for row in table:
        row_id = row["ID"]
        qrow = matching.query_data[row_id]
        RA_query = qrow["RA_query"]
        DEC_query = qrow["DEC_query"]
        assert RA_query == row["RA"]
        assert DEC_query == row["DEC"]
        assert_allclose(
            row["distances"],
            result.nearest_neighbours_distances[row_id][mask_array[row_id]],
            atol=ATOL,
        )
        matched_id = result.nearest_neighbours_indices[row_id][mask_array[row_id]]
        assert all(matched_id == row["ID_matched"])
        assert_allclose(row["RA_matched"], matching.match_data["RA_match"][matched_id])
        assert_allclose(
            row["DEC_matched"], matching.match_data["DEC_match"][matched_id]
        )


def test_match_3d(cosmo, setup_catalogs):
    """Test the match_3d function."""
    query_catalog_path, match_catalog_path = setup_catalogs
    matching = SkyMatch.new_from_fits(
        query_catalog_path=query_catalog_path,
        query_coordinates={"RA": "RA_query", "DEC": "DEC_query", "z": "z_query"},
        match_catalog_path=match_catalog_path,
        match_coordinates={"RA": "RA_match", "DEC": "DEC_match", "z": "z_match"},
    )
    result = matching.match_3d(cosmo, n_nearest_neighbours=10)
    assert result is not None
    dist = Nc.Distance(zf=10.0)
    dist.prepare(cosmo)
    RH_Mpc = cosmo.RH_Mpc()

    match_ra = matching.match_ra
    match_dec = matching.match_dec
    match_theta = np.pi / 2.0 - np.radians(match_dec)
    match_phi = np.radians(match_ra)
    match_r = np.array(dist.angular_diameter_vector(cosmo, matching.match_z)) * RH_Mpc

    query_ra = matching.query_ra
    query_dec = matching.query_dec
    query_theta = np.pi / 2.0 - np.radians(query_dec)
    query_phi = np.radians(query_ra)
    query_r = np.array(dist.angular_diameter_vector(cosmo, matching.query_z)) * RH_Mpc

    query_x1 = query_r * np.sin(query_theta) * np.cos(query_phi)
    query_x2 = query_r * np.sin(query_theta) * np.sin(query_phi)
    query_x3 = query_r * np.cos(query_theta)

    match_x1 = match_r * np.sin(match_theta) * np.cos(match_phi)
    match_x2 = match_r * np.sin(match_theta) * np.sin(match_phi)
    match_x3 = match_r * np.cos(match_theta)

    indices = result.nearest_neighbours_indices
    distances = np.sqrt(
        (query_x1.reshape(-1, 1) - match_x1[indices]) ** 2
        + (query_x2.reshape(-1, 1) - match_x2[indices]) ** 2
        + (query_x3.reshape(-1, 1) - match_x3[indices]) ** 2
    )
    assert_allclose(distances, result.nearest_neighbours_distances)


def test_match_3d_table(cosmo, setup_catalogs):
    """Test the match_3d function."""
    query_catalog_path, match_catalog_path = setup_catalogs
    matching = SkyMatch.new_from_fits(
        query_catalog_path=query_catalog_path,
        query_coordinates={"RA": "RA_query", "DEC": "DEC_query", "z": "z_query"},
        match_catalog_path=match_catalog_path,
        match_coordinates={"RA": "RA_match", "DEC": "DEC_match", "z": "z_match"},
    )
    result = matching.match_3d(cosmo, n_nearest_neighbours=10)
    assert result is not None

    mask_array = np.random.choice(
        [True, False], size=10 * len(matching.query_ra)
    ).reshape(-1, 10)
    table = result.to_table_complete(mask=Mask(mask_array))

    for row in table:
        row_id = row["ID"]
        qrow = matching.query_data[row_id]
        RA_query = qrow["RA_query"]
        DEC_query = qrow["DEC_query"]
        assert RA_query == row["RA"]
        assert DEC_query == row["DEC"]
        assert_allclose(
            row["distances"],
            result.nearest_neighbours_distances[row_id][mask_array[row_id]],
            atol=ATOL,
        )
        matched_id = result.nearest_neighbours_indices[row_id][mask_array[row_id]]
        assert all(matched_id == row["ID_matched"])
        assert_allclose(row["RA_matched"], matching.match_data["RA_match"][matched_id])
        assert_allclose(
            row["DEC_matched"], matching.match_data["DEC_match"][matched_id]
        )


def test_match_2d_extra_columns(cosmo, setup_catalogs_extra_columns):
    """Test the match_2d function with extra columns."""
    query_catalog_path, match_catalog_path = setup_catalogs_extra_columns
    matching = SkyMatch.new_from_fits(
        query_catalog_path=query_catalog_path,
        query_coordinates={"RA": "RA_query", "DEC": "DEC_query", "z": "z_query"},
        match_catalog_path=match_catalog_path,
        match_coordinates={"RA": "RA_match", "DEC": "DEC_match", "z": "z_match"},
    )
    result = matching.match_2d(cosmo, n_nearest_neighbours=10)
    assert result is not None
    mask_array = np.random.choice(
        [True, False], size=10 * len(matching.query_ra)
    ).reshape(-1, 10)
    table = result.to_table_complete(
        mask=Mask(mask_array),
        query_properties={"extra_query": "copy_extra_query"},
        match_properties={"extra_match": "copy_extra_match"},
    )

    for row in table:
        row_id = row["ID"]
        qrow = matching.query_data[row_id]
        RA_query = qrow["RA_query"]
        DEC_query = qrow["DEC_query"]
        assert RA_query == row["RA"]
        assert DEC_query == row["DEC"]
        assert_allclose(
            row["distances"],
            result.nearest_neighbours_distances[row_id][mask_array[row_id]],
            atol=ATOL,
        )
        matched_id = result.nearest_neighbours_indices[row_id][mask_array[row_id]]
        assert all(matched_id == row["ID_matched"])
        assert_allclose(row["RA_matched"], matching.match_data["RA_match"][matched_id])
        assert_allclose(
            row["DEC_matched"], matching.match_data["DEC_match"][matched_id]
        )
        assert_allclose(row["copy_extra_query"], qrow["extra_query"])
        assert_allclose(
            row["copy_extra_match"], matching.match_data["extra_match"][matched_id]
        )


def test_match_3d_missing_z_query(cosmo, setup_catalogs):
    """Test the match_3d function with missing z_query column."""
    query_catalog_path, match_catalog_path = setup_catalogs
    with pytest.raises(
        ValueError,
        match=(
            "To perform a 3D matching, the redshift "
            "must be provided for both catalogs."
        ),
    ):
        matching = SkyMatch.new_from_fits(
            query_catalog_path=query_catalog_path,
            query_coordinates={"RA": "RA_query", "DEC": "DEC_query"},
            match_catalog_path=match_catalog_path,
            match_coordinates={"RA": "RA_match", "DEC": "DEC_match", "z": "z_match"},
        )
        _ = matching.match_3d(cosmo, n_nearest_neighbours=10)


def test_match_3d_missing_z_match(cosmo, setup_catalogs):
    """Test the match_3d function with missing z_match column."""
    query_catalog_path, match_catalog_path = setup_catalogs
    with pytest.raises(
        ValueError,
        match=(
            "To perform a 3D matching, the redshift "
            "must be provided for both catalogs."
        ),
    ):
        matching = SkyMatch.new_from_fits(
            query_catalog_path=query_catalog_path,
            query_coordinates={"RA": "RA_query", "DEC": "DEC_query", "z": "z_query"},
            match_catalog_path=match_catalog_path,
            match_coordinates={"RA": "RA_match", "DEC": "DEC_match"},
        )
        _ = matching.match_3d(cosmo, n_nearest_neighbours=10)

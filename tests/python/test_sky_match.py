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
import re
import pytest
import numpy as np
from numpy.testing import assert_allclose

pytest.importorskip("astropy")
# flake8: noqa: E402
# pylint: disable=wrong-import-position

from astropy.table import Table
from astropy.io import fits

from numcosmo_py import Nc, Ncm
from numcosmo_py.sky_match import (
    SkyMatch,
    SkyMatchResult,
    DistanceMethod,
    Mask,
    SelectionCriteria,
    SharedFractionMethod,
)

Ncm.cfg_init()

ATOL = 1.0e-7
RTOL = 1.0e-5
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
            "z_err": [np.random.uniform(0.01, 0.02) for _ in range(QUERY_SIZE)],
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
            "z_err": [np.random.uniform(0.01, 0.02) for _ in range(MATCH_SIZE)],
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


def test_mask() -> None:
    """Test the mask_new function."""
    mask_all_true = Mask(mask=np.ones((10, 10), dtype=bool))
    assert mask_all_true.shape == (10, 10)
    assert mask_all_true.array.shape == (10, 10)

    mask_all_false = Mask(mask=np.zeros((10, 10), dtype=bool))
    assert mask_all_false.shape == (10, 10)
    assert mask_all_false.array.shape == (10, 10)

    mask = mask_all_false & mask_all_true
    assert mask.shape == (10, 10)
    assert mask.array.shape == (10, 10)

    assert np.all(mask.array == mask_all_false.array)

    assert np.all(mask_all_true & mask_all_true == mask_all_true)
    assert np.all(~mask_all_false == mask_all_true)


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
    with pytest.raises(ValueError, match="No table found"):
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
                angular_distance,
                result.nearest_neighbours_distances,
                atol=ATOL,
                rtol=RTOL,
            )
        case DistanceMethod.QUERY_RADIUS:
            query_r = (
                np.array(dist.angular_diameter_array(cosmo, matching.query_z)) * RH_Mpc
            )
            distances = query_r.reshape(-1, 1) * 2.0 * np.sin(angular_distance / 2.0)
            assert_allclose(
                distances, result.nearest_neighbours_distances, atol=ATOL, rtol=RTOL
            )
        case DistanceMethod.MATCH_RADIUS:
            match_r = (
                np.array(dist.angular_diameter_array(cosmo, matching.match_z)) * RH_Mpc
            )
            distances = (
                match_r[result.nearest_neighbours_indices]
                * 2.0
                * np.sin(angular_distance / 2.0)
            )
            assert_allclose(
                distances, result.nearest_neighbours_distances, atol=ATOL, rtol=RTOL
            )
        case DistanceMethod.MIN_RADIUS:
            query_r = (
                np.array(dist.angular_diameter_array(cosmo, matching.query_z)) * RH_Mpc
            )
            match_r = (
                np.array(dist.angular_diameter_array(cosmo, matching.match_z)) * RH_Mpc
            )
            distances = (
                np.minimum(
                    query_r.reshape(-1, 1), match_r[result.nearest_neighbours_indices]
                )
                * 2.0
                * np.sin(angular_distance / 2.0)
            )
            assert_allclose(
                distances, result.nearest_neighbours_distances, atol=ATOL, rtol=RTOL
            )
        case DistanceMethod.MAX_RADIUS:
            query_r = (
                np.array(dist.angular_diameter_array(cosmo, matching.query_z)) * RH_Mpc
            )
            match_r = (
                np.array(dist.angular_diameter_array(cosmo, matching.match_z)) * RH_Mpc
            )
            distances = (
                np.maximum(
                    query_r.reshape(-1, 1), match_r[result.nearest_neighbours_indices]
                )
                * 2.0
                * np.sin(angular_distance / 2.0)
            )
            assert_allclose(
                distances, result.nearest_neighbours_distances, atol=ATOL, rtol=RTOL
            )
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
            rtol=RTOL,
        )
        matched_id = result.nearest_neighbours_indices[row_id][mask_array[row_id]]
        assert all(matched_id == row["ID_matched"])
        assert_allclose(
            row["RA_matched"],
            matching.match_data["RA_match"][matched_id],
            atol=ATOL,
            rtol=RTOL,
        )
        assert_allclose(
            row["DEC_matched"],
            matching.match_data["DEC_match"][matched_id],
            atol=ATOL,
            rtol=RTOL,
        )
    # No mask
    table = result.to_table_complete()

    for row in table:
        row_id = row["ID"]
        qrow = matching.query_data[row_id]
        RA_query = qrow["RA_query"]
        DEC_query = qrow["DEC_query"]
        assert RA_query == row["RA"]
        assert DEC_query == row["DEC"]
        assert_allclose(
            row["distances"],
            result.nearest_neighbours_distances[row_id],
            atol=ATOL,
            rtol=RTOL,
        )
        matched_id = result.nearest_neighbours_indices[row_id]
        assert all(matched_id == row["ID_matched"])
        assert_allclose(
            row["RA_matched"],
            matching.match_data["RA_match"][matched_id],
            atol=ATOL,
            rtol=RTOL,
        )
        assert_allclose(
            row["DEC_matched"],
            matching.match_data["DEC_match"][matched_id],
            atol=ATOL,
            rtol=RTOL,
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
    match_r = np.array(dist.angular_diameter_array(cosmo, matching.match_z)) * RH_Mpc

    query_ra = matching.query_ra
    query_dec = matching.query_dec
    query_theta = np.pi / 2.0 - np.radians(query_dec)
    query_phi = np.radians(query_ra)
    query_r = np.array(dist.angular_diameter_array(cosmo, matching.query_z)) * RH_Mpc

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
    assert_allclose(
        distances, result.nearest_neighbours_distances, atol=ATOL, rtol=RTOL
    )


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
            rtol=RTOL,
        )
        matched_id = result.nearest_neighbours_indices[row_id][mask_array[row_id]]
        assert all(matched_id == row["ID_matched"])
        assert_allclose(
            row["RA_matched"],
            matching.match_data["RA_match"][matched_id],
            atol=ATOL,
            rtol=RTOL,
        )
        assert_allclose(
            row["DEC_matched"],
            matching.match_data["DEC_match"][matched_id],
            atol=ATOL,
            rtol=RTOL,
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
            rtol=RTOL,
        )
        matched_id = result.nearest_neighbours_indices[row_id][mask_array[row_id]]
        assert all(matched_id == row["ID_matched"])
        assert_allclose(row["RA_matched"], matching.match_data["RA_match"][matched_id])
        assert_allclose(
            row["DEC_matched"],
            matching.match_data["DEC_match"][matched_id],
            atol=ATOL,
            rtol=RTOL,
        )
        assert_allclose(
            row["copy_extra_query"], qrow["extra_query"], atol=ATOL, rtol=RTOL
        )
        assert_allclose(
            row["copy_extra_match"],
            matching.match_data["extra_match"][matched_id],
            atol=ATOL,
            rtol=RTOL,
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


def test_match_2d_missing_z_query(cosmo, setup_catalogs):
    """Test the match_3d function with missing z_query column."""
    query_catalog_path, match_catalog_path = setup_catalogs
    with pytest.raises(
        ValueError,
        match=(
            "To perform a matching, the redshift must be provided for both catalogs."
        ),
    ):
        matching = SkyMatch.new_from_fits(
            query_catalog_path=query_catalog_path,
            query_coordinates={"RA": "RA_query", "DEC": "DEC_query"},
            match_catalog_path=match_catalog_path,
            match_coordinates={"RA": "RA_match", "DEC": "DEC_match", "z": "z_match"},
        )
        _ = matching.match_2d(cosmo, n_nearest_neighbours=10)


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


def test_match_2d_missing_z_match(cosmo, setup_catalogs):
    """Test the match_3d function with missing z_match column."""
    query_catalog_path, match_catalog_path = setup_catalogs
    with pytest.raises(
        ValueError,
        match=(
            "To perform a matching, the redshift must be provided for both catalogs."
        ),
    ):
        matching = SkyMatch.new_from_fits(
            query_catalog_path=query_catalog_path,
            query_coordinates={"RA": "RA_query", "DEC": "DEC_query", "z": "z_query"},
            match_catalog_path=match_catalog_path,
            match_coordinates={"RA": "RA_match", "DEC": "DEC_match"},
        )
        _ = matching.match_2d(cosmo, n_nearest_neighbours=10)


def test_match_3d_filter_distance(cosmo, setup_catalogs):
    """Test the match_3d function with filter."""
    query_catalog_path, match_catalog_path = setup_catalogs
    matching = SkyMatch.new_from_fits(
        query_catalog_path=query_catalog_path,
        query_coordinates={"RA": "RA_query", "DEC": "DEC_query", "z": "z_query"},
        match_catalog_path=match_catalog_path,
        match_coordinates={"RA": "RA_match", "DEC": "DEC_match", "z": "z_match"},
    )
    result = matching.match_3d(cosmo, n_nearest_neighbours=10)
    assert result is not None
    min_distance = 20.0

    mask = result.filter_mask_by_distance(min_distance)
    table = result.to_table_complete(mask=mask)
    assert len(table) == QUERY_SIZE

    for row in table:
        assert all(row["distances"] < min_distance)

    table_inverted = result.to_table_complete(mask=~mask)
    assert len(table_inverted) == QUERY_SIZE

    for row in table_inverted:
        assert all(row["distances"] >= min_distance)


@pytest.mark.parametrize(
    ["sigma0", "nsigma"], [(0.01, 0.734), (0.02, 0.4), (0.007, 1.5)]
)
def test_match_3d_filter_z(cosmo, setup_catalogs, sigma0, nsigma):
    """Test the match_3d function with filter."""
    query_catalog_path, match_catalog_path = setup_catalogs
    matching = SkyMatch.new_from_fits(
        query_catalog_path=query_catalog_path,
        query_coordinates={"RA": "RA_query", "DEC": "DEC_query", "z": "z_query"},
        match_catalog_path=match_catalog_path,
        match_coordinates={"RA": "RA_match", "DEC": "DEC_match", "z": "z_match"},
    )
    result = matching.match_3d(cosmo, n_nearest_neighbours=10)
    assert result is not None

    mask = result.filter_mask_by_redshift_proximity(sigma0, nsigma)
    table = result.to_table_complete(mask=mask)
    assert len(table) == QUERY_SIZE

    for row in table:
        match_sigma_z = nsigma * sigma0 * (1.0 + row["z_matched"])
        query_sigma_z = nsigma * sigma0 * (1.0 + row["z"])
        max_z_dist = match_sigma_z + query_sigma_z
        assert all(np.abs(row["z_matched"] - row["z"]) < max_z_dist)

    table_inverted = result.to_table_complete(mask=~mask)
    assert len(table_inverted) == QUERY_SIZE

    for row in table_inverted:
        match_sigma_z = nsigma * sigma0 * (1.0 + row["z_matched"])
        query_sigma_z = nsigma * sigma0 * (1.0 + row["z"])
        max_z_dist = match_sigma_z + query_sigma_z
        assert all(np.abs(row["z_matched"] - row["z"]) >= max_z_dist)


def test_match_3d_filter_z_wrong_column_name(cosmo, setup_catalogs):
    """Test the match_3d function with filter."""
    query_catalog_path, match_catalog_path = setup_catalogs
    matching = SkyMatch.new_from_fits(
        query_catalog_path=query_catalog_path,
        query_coordinates={"RA": "RA_query", "DEC": "DEC_query", "z": "z_query"},
        match_catalog_path=match_catalog_path,
        match_coordinates={"RA": "RA_match", "DEC": "DEC_match", "z": "z_match"},
    )
    result = matching.match_3d(cosmo, n_nearest_neighbours=10)
    assert result is not None

    with pytest.raises(ValueError, match="Column kawabunga not found in match data."):
        _ = result.filter_mask_by_redshift_proximity(
            0.1, 1.0, match_sigma_z_column="kawabunga"
        )
    with pytest.raises(ValueError, match="Column kawabunga not found in query data."):
        _ = result.filter_mask_by_redshift_proximity(
            0.1, 1.0, query_sigma_z_column="kawabunga"
        )


@pytest.mark.parametrize(
    ["sigma0", "nsigma"], [(0.01, 0.734), (0.02, 0.4), (0.007, 1.5)]
)
def test_match_3d_filter_z_match_z_err(cosmo, setup_catalogs, sigma0, nsigma):
    """Test the match_3d function with filter."""
    query_catalog_path, match_catalog_path = setup_catalogs
    matching = SkyMatch.new_from_fits(
        query_catalog_path=query_catalog_path,
        query_coordinates={"RA": "RA_query", "DEC": "DEC_query", "z": "z_query"},
        match_catalog_path=match_catalog_path,
        match_coordinates={"RA": "RA_match", "DEC": "DEC_match", "z": "z_match"},
    )
    result = matching.match_3d(cosmo, n_nearest_neighbours=10)
    assert result is not None

    mask = result.filter_mask_by_redshift_proximity(
        sigma0, nsigma, match_sigma_z_column="z_err"
    )
    table = result.to_table_complete(
        mask=mask, match_properties={"z_err": "match_z_err"}
    )
    assert len(table) == QUERY_SIZE

    for row in table:
        match_sigma_z = nsigma * row["match_z_err"]
        query_sigma_z = nsigma * sigma0 * (1.0 + row["z"])
        max_z_dist = match_sigma_z + query_sigma_z
        assert all(np.abs(row["z_matched"] - row["z"]) < max_z_dist)

    table_inverted = result.to_table_complete(
        mask=~mask, match_properties={"z_err": "match_z_err"}
    )
    assert len(table_inverted) == QUERY_SIZE

    for row in table_inverted:
        match_sigma_z = nsigma * row["match_z_err"]
        query_sigma_z = nsigma * sigma0 * (1.0 + row["z"])
        max_z_dist = match_sigma_z + query_sigma_z
        assert all(np.abs(row["z_matched"] - row["z"]) >= max_z_dist)


@pytest.mark.parametrize(
    ["sigma0", "nsigma"], [(0.01, 0.734), (0.02, 0.4), (0.007, 1.5)]
)
def test_match_3d_filter_z_query_z_err(cosmo, setup_catalogs, sigma0, nsigma):
    """Test the match_3d function with filter."""
    query_catalog_path, match_catalog_path = setup_catalogs
    matching = SkyMatch.new_from_fits(
        query_catalog_path=query_catalog_path,
        query_coordinates={"RA": "RA_query", "DEC": "DEC_query", "z": "z_query"},
        match_catalog_path=match_catalog_path,
        match_coordinates={"RA": "RA_match", "DEC": "DEC_match", "z": "z_match"},
    )
    result = matching.match_3d(cosmo, n_nearest_neighbours=10)
    assert result is not None

    mask = result.filter_mask_by_redshift_proximity(
        sigma0, nsigma, query_sigma_z_column="z_err"
    )
    table = result.to_table_complete(
        mask=mask, query_properties={"z_err": "query_z_err"}
    )
    assert len(table) == QUERY_SIZE

    for row in table:
        match_sigma_z = nsigma * sigma0 * (1.0 + row["z_matched"])
        query_sigma_z = nsigma * row["query_z_err"]
        max_z_dist = match_sigma_z + query_sigma_z
        assert all(np.abs(row["z_matched"] - row["z"]) < max_z_dist)

    table_inverted = result.to_table_complete(
        mask=~mask, query_properties={"z_err": "query_z_err"}
    )
    assert len(table_inverted) == QUERY_SIZE

    for row in table_inverted:
        match_sigma_z = nsigma * sigma0 * (1.0 + row["z_matched"])
        query_sigma_z = nsigma * row["query_z_err"]
        max_z_dist = match_sigma_z + query_sigma_z
        assert all(np.abs(row["z_matched"] - row["z"]) >= max_z_dist)


@pytest.mark.parametrize(
    ["sigma0", "nsigma"], [(0.01, 0.734), (0.02, 0.4), (0.007, 1.5)]
)
def test_match_3d_filter_z_both_z_err(cosmo, setup_catalogs, sigma0, nsigma):
    """Test the match_3d function with filter."""
    query_catalog_path, match_catalog_path = setup_catalogs
    matching = SkyMatch.new_from_fits(
        query_catalog_path=query_catalog_path,
        query_coordinates={"RA": "RA_query", "DEC": "DEC_query", "z": "z_query"},
        match_catalog_path=match_catalog_path,
        match_coordinates={"RA": "RA_match", "DEC": "DEC_match", "z": "z_match"},
    )
    result = matching.match_3d(cosmo, n_nearest_neighbours=10)
    assert result is not None

    mask = result.filter_mask_by_redshift_proximity(
        sigma0, nsigma, query_sigma_z_column="z_err", match_sigma_z_column="z_err"
    )
    table = result.to_table_complete(
        mask=mask,
        query_properties={"z_err": "query_z_err"},
        match_properties={"z_err": "match_z_err"},
    )
    assert len(table) == QUERY_SIZE

    for row in table:
        match_sigma_z = nsigma * row["match_z_err"]
        query_sigma_z = nsigma * row["query_z_err"]
        max_z_dist = match_sigma_z + query_sigma_z
        assert all(np.abs(row["z_matched"] - row["z"]) < max_z_dist)

    table_inverted = result.to_table_complete(
        mask=~mask,
        query_properties={"z_err": "query_z_err"},
        match_properties={"z_err": "match_z_err"},
    )
    assert len(table_inverted) == QUERY_SIZE

    for row in table_inverted:
        match_sigma_z = nsigma * row["match_z_err"]
        query_sigma_z = nsigma * row["query_z_err"]
        max_z_dist = match_sigma_z + query_sigma_z
        assert all(np.abs(row["z_matched"] - row["z"]) >= max_z_dist)


def test_match_2d_best_distance(cosmo, setup_catalogs, distance_method):
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

    mask = result.filter_mask_by_distance(20.0)
    best = result.select_best(selection_criteria=SelectionCriteria.DISTANCES, mask=mask)
    assert best is not None

    for i, query_index in enumerate(best.query_indices):
        best_index_row = np.argmin(result.nearest_neighbours_distances[query_index])
        assert (
            result.nearest_neighbours_indices[query_index][best_index_row]
            == best.indices[i]
        )

    # Best no mask
    best = result.select_best(selection_criteria=SelectionCriteria.DISTANCES)
    assert best is not None

    for i, query_index in enumerate(best.query_indices):
        best_index_row = np.argmin(result.nearest_neighbours_distances[query_index])
        assert (
            result.nearest_neighbours_indices[query_index][best_index_row]
            == best.indices[i]
        )


def test_match_2d_best_redshift(cosmo, setup_catalogs, distance_method):
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

    mask = result.filter_mask_by_distance(20.0)
    best = result.select_best(
        selection_criteria=SelectionCriteria.REDSHIFT_PROXIMITY, mask=mask
    )
    assert best is not None

    match_z = result.sky_match.match_z
    query_z = result.sky_match.query_z

    for i, query_index in enumerate(best.query_indices):
        idx = result.nearest_neighbours_indices[query_index][mask.array[query_index]]
        best_z_index = idx[np.argmin(np.abs(match_z[idx] - query_z[query_index]))]
        assert best_z_index == best.indices[i]


def test_match_2d_best_more_massive(cosmo, setup_catalogs, distance_method):
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

    mask = result.filter_mask_by_distance(20.0)
    best = result.select_best(
        selection_criteria=SelectionCriteria.MORE_MASSIVE,
        mask=mask,
        more_massive_column="RA_match",
    )
    assert best is not None

    match_ra = result.sky_match.match_ra

    for i, query_index in enumerate(best.query_indices):
        idx = result.nearest_neighbours_indices[query_index][mask.array[query_index]]
        best_z_index = idx[np.argmax(match_ra[idx])]
        assert best_z_index == best.indices[i]


def test_match_2d_best_more_massive_wrong_column(cosmo, setup_catalogs):
    """Test the match_2d function."""
    query_catalog_path, match_catalog_path = setup_catalogs
    matching = SkyMatch.new_from_fits(
        query_catalog_path=query_catalog_path,
        query_coordinates={"RA": "RA_query", "DEC": "DEC_query", "z": "z_query"},
        match_catalog_path=match_catalog_path,
        match_coordinates={"RA": "RA_match", "DEC": "DEC_match", "z": "z_match"},
    )
    result = matching.match_2d(cosmo, n_nearest_neighbours=10)
    assert result is not None

    mask = result.filter_mask_by_distance(20.0)
    with pytest.raises(
        ValueError,
        match=re.escape(
            "A more_massive_column (Not_a_column) must be "
            "provided and present in the match data."
        ),
    ):
        _ = result.select_best(
            selection_criteria=SelectionCriteria.MORE_MASSIVE,
            mask=mask,
            more_massive_column="Not_a_column",
        )


def test_match_2d_best_distance_table(cosmo, setup_catalogs, distance_method):
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

    mask = result.filter_mask_by_distance(20.0)
    best = result.select_best(selection_criteria=SelectionCriteria.DISTANCES, mask=mask)
    assert best is not None

    table = result.to_table_best(best)
    assert len(table) == len(best.query_indices)

    for i, query_index in enumerate(best.query_indices):
        assert table["RA"][i] == result.sky_match.query_ra[query_index]
        assert table["DEC"][i] == result.sky_match.query_dec[query_index]
        assert table["z"][i] == result.sky_match.query_z[query_index]
        assert table["RA_matched"][i] == result.sky_match.match_ra[best.indices[i]]
        assert table["DEC_matched"][i] == result.sky_match.match_dec[best.indices[i]]
        assert table["z_matched"][i] == result.sky_match.match_z[best.indices[i]]

    table = result.to_table_best(
        best,
        query_properties={"ID_query": "ID_query_copy"},
        match_properties={"ID_match": "ID_match_copy"},
    )
    assert len(table) == len(best.query_indices)

    for i, query_index in enumerate(best.query_indices):
        assert table["RA"][i] == result.sky_match.query_ra[query_index]
        assert table["DEC"][i] == result.sky_match.query_dec[query_index]
        assert table["z"][i] == result.sky_match.query_z[query_index]
        assert table["RA_matched"][i] == result.sky_match.match_ra[best.indices[i]]
        assert table["DEC_matched"][i] == result.sky_match.match_dec[best.indices[i]]
        assert table["z_matched"][i] == result.sky_match.match_z[best.indices[i]]
        assert (
            table["ID_query_copy"][i]
            == result.sky_match.query_data["ID_query"][query_index]
        )
        assert (
            table["ID_match_copy"][i]
            == result.sky_match.match_data["ID_match"][best.indices[i]]
        )


def test_match_2d_best_cross(cosmo, setup_catalogs, distance_method):
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

    mask = result.filter_mask_by_distance(20.0)
    best = result.select_best(selection_criteria=SelectionCriteria.DISTANCES, mask=mask)

    inverted_matching = matching.invert_query_match()
    inverted_result = inverted_matching.match_2d(
        cosmo, n_nearest_neighbours=10, distance_method=distance_method
    )
    assert inverted_result is not None

    inverted_mask = inverted_result.filter_mask_by_distance(20.0)
    inverted_best = inverted_result.select_best(
        selection_criteria=SelectionCriteria.DISTANCES, mask=inverted_mask
    )

    cross = best.get_cross_match_indices(inverted_best)

    best_dict = best.query_match_dict
    inversed_best_dict = inverted_best.query_match_dict

    for query_i, match_i in cross.items():
        assert best_dict[query_i] == match_i
        assert inversed_best_dict[match_i] == query_i


# ---------------------------------------------------------------------------
# ID matching tests
# ---------------------------------------------------------------------------
#
# Known-answer mock. Two query objects and two match objects, each described by
# a member catalogue. Membership overlap (Jaccard) drives the matching:
#
#   q0 = {a, b, c}   t100 = {a, b, f}   shared {a, b}    -> J = 2/4 = 0.5
#   q1 = {d, e}      t101 = {d, g, h}   shared {d}       -> J = 1/4 = 0.25
#
# With all pmem = 1 the four SharedFractionMethod variants reduce to simple
# member-count ratios, so the best match per query is unambiguous.

QUERY_ID_COL = "ID_q"
MATCH_ID_COL = "ID_m"


@pytest.fixture(name="setup_id_match")
def fixture_setup_id_match():
    """Build a SkyMatch wired for ID matching with a known-answer overlap."""
    query_data = Table(
        {
            QUERY_ID_COL: [0, 1],
            "RA_query": [10.0, 20.0],
            "DEC_query": [-10.0, -20.0],
            "z_query": [0.5, 0.6],
        }
    )
    match_data = Table(
        {
            MATCH_ID_COL: [100, 101],
            "RA_match": [10.1, 20.1],
            "DEC_match": [-10.1, -20.1],
            "z_match": [0.5, 0.6],
        }
    )
    query_member = Table(
        {
            QUERY_ID_COL: [0, 0, 0, 1, 1],
            "MemID": ["a", "b", "c", "d", "e"],
            "pmem": [1.0, 1.0, 1.0, 1.0, 1.0],
        }
    )
    match_member = Table(
        {
            MATCH_ID_COL: [100, 100, 100, 101, 101, 101],
            "MemID": ["a", "b", "f", "d", "g", "h"],
            "pmem": [1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
        }
    )
    return SkyMatch(
        query_data=query_data,
        query_coordinates={"RA": "RA_query", "DEC": "DEC_query", "z": "z_query"},
        query_member_data=query_member,
        query_ids={"ID": QUERY_ID_COL, "MemberID": "MemID", "pmem": "pmem"},
        match_data=match_data,
        match_coordinates={"RA": "RA_match", "DEC": "DEC_match", "z": "z_match"},
        match_member_data=match_member,
        match_ids={"ID": MATCH_ID_COL, "MemberID": "MemID", "pmem": "pmem"},
    )


def _coeff_for_query(result, sky_match, query_row, match_external_id):
    """Return the linking coefficient stored for a given query/match pair."""
    match_ids = list(result.candidate_ids[query_row])
    coeffs = list(result.candidate_coefficients[query_row])
    return coeffs[match_ids.index(match_external_id)]


def test_match_id_coefficients(setup_id_match):
    """The Jaccard linking coefficients match the known-answer overlap."""
    result = setup_id_match.match_ID()
    assert result is not None

    assert_allclose(_coeff_for_query(result, setup_id_match, 0, 100), 0.5)
    assert_allclose(_coeff_for_query(result, setup_id_match, 1, 101), 0.25)


def test_match_id_select_best(setup_id_match):
    """select_best resolves a one-to-one matching by row index into match_data."""
    result = setup_id_match.match_ID()
    best = result.select_best()

    # Both queries have a match.
    assert list(best.query_filter) == [True, True]
    # Indices are rows into match_data; row 0 -> ID 100, row 1 -> ID 101.
    matched_ids = setup_id_match.match_id[best.indices]
    assert list(matched_ids) == [100, 101]


def test_match_id_to_table_best(setup_id_match):
    """to_table_best exposes the external IDs and linking coefficients."""
    result = setup_id_match.match_ID()
    best = result.select_best()
    table = result.to_table_best(best)

    assert len(table) == 2
    by_query = {row["query_id"]: row for row in table}
    assert by_query[0]["match_id"] == 100
    assert by_query[1]["match_id"] == 101
    assert_allclose(by_query[0]["linking_coefficient"], 0.5)
    assert_allclose(by_query[1]["linking_coefficient"], 0.25)
    # Coordinates come from the match catalogue rows.
    assert_allclose(by_query[0]["RA_matched"], 10.1)
    assert_allclose(by_query[1]["RA_matched"], 20.1)


def test_match_id_to_table_complete(setup_id_match):
    """to_table_complete lists every candidate per query with its coefficient."""
    result = setup_id_match.match_ID()
    table = result.to_table_complete()

    assert len(table) == 2
    # Each query in this mock has exactly one candidate match.
    for row in table:
        assert len(row["match_id"]) == 1
        assert len(row["linking_coefficient"]) == 1


def test_match_id_use_shared_fraction_true(setup_id_match):
    """use_shared_fraction=True must not raise and yields valid coefficients."""
    result = setup_id_match.match_ID(use_shared_fraction=True)
    assert result is not None

    for query_row in (0, 1):
        for coeff in result.candidate_coefficients[query_row]:
            assert np.isfinite(coeff)
            assert 0.0 <= coeff <= 1.0


@pytest.mark.parametrize("method", list(SharedFractionMethod))
def test_match_id_shared_fraction_methods(setup_id_match, method):
    """All shared-fraction methods run and pick the correct best match."""
    result = setup_id_match.match_ID(shared_fraction_method=method)
    best = result.select_best()

    matched_ids = setup_id_match.match_id[best.indices]
    matched_by_query = dict(zip(best.query_indices, matched_ids))
    assert matched_by_query[0] == 100
    assert matched_by_query[1] == 101


@pytest.fixture(name="setup_id_match_contested")
def fixture_setup_id_match_contested():
    """Two queries competing for a single match object.

    q0 = {a, b, c}  J(q0, t100) = 3/3 = 1.0
    q1 = {b, c, d}  J(q1, t100) = 2/4 = 0.5
    A one-to-one matching must award t100 to q0 (higher weight); q1 is left out.
    """
    query_data = Table(
        {
            QUERY_ID_COL: [0, 1],
            "RA_query": [10.0, 20.0],
            "DEC_query": [-10.0, -20.0],
            "z_query": [0.5, 0.6],
        }
    )
    match_data = Table(
        {
            MATCH_ID_COL: [100],
            "RA_match": [10.1],
            "DEC_match": [-10.1],
            "z_match": [0.5],
        }
    )
    query_member = Table(
        {
            QUERY_ID_COL: [0, 0, 0, 1, 1, 1],
            "MemID": ["a", "b", "c", "b", "c", "d"],
            "pmem": [1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
        }
    )
    match_member = Table(
        {
            MATCH_ID_COL: [100, 100, 100],
            "MemID": ["a", "b", "c"],
            "pmem": [1.0, 1.0, 1.0],
        }
    )
    return SkyMatch(
        query_data=query_data,
        query_coordinates={"RA": "RA_query", "DEC": "DEC_query", "z": "z_query"},
        query_member_data=query_member,
        query_ids={"ID": QUERY_ID_COL, "MemberID": "MemID", "pmem": "pmem"},
        match_data=match_data,
        match_coordinates={"RA": "RA_match", "DEC": "DEC_match", "z": "z_match"},
        match_member_data=match_member,
        match_ids={"ID": MATCH_ID_COL, "MemberID": "MemID", "pmem": "pmem"},
    )


def test_match_id_select_best_contested(setup_id_match_contested):
    """A shared match object is awarded to the higher-weight query only."""
    sky_match = setup_id_match_contested
    result = sky_match.match_ID()
    best = result.select_best()

    # Only one query keeps the match (the one with full overlap, q0).
    assert list(best.query_filter) == [True, False]
    matched_ids = sky_match.match_id[best.indices]
    assert list(matched_ids) == [100]


# ---------------------------------------------------------------------------
# Global distance-assignment matching (select_best_assignment)
# ---------------------------------------------------------------------------
#
# These build a SkyMatchResult directly from synthetic neighbour index/distance
# arrays so the assignment logic is exercised deterministically, independent of
# the cosmology and the kNN search.


def _make_distance_result(indices, distances, n_match=None):
    """Build a SkyMatchResult from explicit (n_query, k) index/distance arrays."""
    indices = np.array(indices, dtype=np.int64)
    distances = np.array(distances, dtype=np.float64)
    n_query = indices.shape[0]
    if n_match is None:
        n_match = int(indices.max()) + 1

    def _coords(n):
        return Table(
            {
                "RA": [float(i) for i in range(n)],
                "DEC": [-float(i) for i in range(n)],
                "z": [0.5] * n,
            }
        )

    coords = {"RA": "RA", "DEC": "DEC", "z": "z"}
    sky_match = SkyMatch(_coords(n_query), coords, _coords(n_match), coords)
    return SkyMatchResult(sky_match, indices, distances)


def test_assignment_resolves_collision():
    """Global assignment reassigns to reach a lower total than greedy collisions."""
    # q0: rows [0, 1] d [1.0, 3.0]; q1: rows [0, 2] d [1.5, 4.0]
    result = _make_distance_result([[0, 1], [0, 2]], [[1.0, 3.0], [1.5, 4.0]])

    # Greedy double-books match row 0 (both queries' nearest).
    greedy = result.select_best(selection_criteria=SelectionCriteria.DISTANCES)
    assert greedy.query_match_dict == {0: 0, 1: 0}

    # Assignment minimizes the total distance: q0->1 (3.0) + q1->0 (1.5) = 4.5,
    # which beats the greedy-consistent q0->0 (1.0) + q1->2 (4.0) = 5.0.
    best = result.select_best_assignment()
    assert best.query_match_dict == {0: 1, 1: 0}
    # One-to-one: no match row used twice.
    assert len(set(best.indices)) == len(best.indices)


def test_assignment_matches_greedy_when_disjoint():
    """With no contended match, assignment equals the greedy nearest choice."""
    result = _make_distance_result([[0, 1], [2, 3]], [[1.0, 5.0], [2.0, 6.0]])

    greedy = result.select_best(selection_criteria=SelectionCriteria.DISTANCES)
    best = result.select_best_assignment()
    assert best.query_match_dict == greedy.query_match_dict == {0: 0, 1: 2}


def test_assignment_respects_mask():
    """Masked-out candidates are not eligible for the assignment."""
    result = _make_distance_result([[0, 1], [0, 2]], [[1.0, 3.0], [1.5, 4.0]])

    # Forbid q0's second candidate (row 1); q0 can then only take row 0, forcing
    # q1 onto row 2 to keep both matched.
    mask = Mask(np.array([[True, False], [True, True]]))
    best = result.select_best_assignment(mask=mask)
    assert best.query_match_dict == {0: 0, 1: 2}


def test_assignment_one_to_one_drops_losers():
    """When several queries can only reach one match, only the closest is kept."""
    result = _make_distance_result([[0], [0], [0]], [[1.0], [2.0], [3.0]], n_match=1)

    best = result.select_best_assignment()
    assert list(best.query_filter) == [True, False, False]
    assert list(best.indices) == [0]


def test_assignment_empty_mask_yields_no_matches():
    """An all-False mask produces an empty assignment."""
    result = _make_distance_result([[0, 1], [0, 2]], [[1.0, 3.0], [1.5, 4.0]])

    mask = Mask(np.zeros((2, 2), dtype=bool))
    best = result.select_best_assignment(mask=mask)
    assert not best.query_filter.any()
    assert len(best.indices) == 0


# ---------------------------------------------------------------------------
# ID matching: validation, pmem requirements and property pass-through
# ---------------------------------------------------------------------------


def _id_sky_match(*, with_pmem=True, query_member=None, match_member=None):
    """Build an ID SkyMatch with mass columns and optional pmem member weights."""
    query_data = Table(
        {
            QUERY_ID_COL: [0, 1],
            "RA_query": [10.0, 20.0],
            "DEC_query": [-10.0, -20.0],
            "z_query": [0.5, 0.6],
            "mass_q": [2.0, 1.0],
        }
    )
    match_data = Table(
        {
            MATCH_ID_COL: [100, 101],
            "RA_match": [10.1, 20.1],
            "DEC_match": [-10.1, -20.1],
            "z_match": [0.5, 0.6],
            "mass_m": [5.0, 4.0],
        }
    )
    if query_member is None:
        query_member = {
            QUERY_ID_COL: [0, 0, 0, 1, 1],
            "MemID": ["a", "b", "c", "d", "e"],
        }
        if with_pmem:
            query_member["pmem"] = [1.0] * 5
    if match_member is None:
        match_member = {
            MATCH_ID_COL: [100, 100, 100, 101, 101, 101],
            "MemID": ["a", "b", "f", "d", "g", "h"],
        }
        if with_pmem:
            match_member["pmem"] = [1.0] * 6
    query_ids = {"ID": QUERY_ID_COL, "MemberID": "MemID"}
    match_ids = {"ID": MATCH_ID_COL, "MemberID": "MemID"}
    if with_pmem:
        query_ids["pmem"] = "pmem"
        match_ids["pmem"] = "pmem"
    return SkyMatch(
        query_data=query_data,
        query_coordinates={"RA": "RA_query", "DEC": "DEC_query", "z": "z_query"},
        query_member_data=Table(query_member),
        query_ids=query_ids,
        match_data=match_data,
        match_coordinates={"RA": "RA_match", "DEC": "DEC_match", "z": "z_match"},
        match_member_data=Table(match_member),
        match_ids=match_ids,
    )


def _minimal_id_sky_match(*, query_ids):
    """Build a tiny ID SkyMatch, parametrizing only the (possibly invalid) ID map."""
    return SkyMatch(
        query_data=Table({QUERY_ID_COL: [0, 1], "RA": [1.0, 2.0], "DEC": [-1.0, -2.0]}),
        query_coordinates={"RA": "RA", "DEC": "DEC"},
        query_member_data=Table({QUERY_ID_COL: [0, 1], "MemID": ["a", "b"]}),
        query_ids=query_ids,
        match_data=Table({MATCH_ID_COL: [9], "RA": [1.0], "DEC": [-1.0]}),
        match_coordinates={"RA": "RA", "DEC": "DEC"},
        match_member_data=Table({MATCH_ID_COL: [9], "MemID": ["a"]}),
        match_ids={"ID": MATCH_ID_COL, "MemberID": "MemID"},
    )


def test_id_check_missing_member_key():
    """Omitting MemberID from the ID map is rejected at construction."""
    with pytest.raises(ValueError, match="ID and MemberID must be provided"):
        _minimal_id_sky_match(query_ids={"ID": QUERY_ID_COL})  # no MemberID


def test_id_check_member_column_not_found():
    """A MemberID column absent from the member table is rejected."""
    with pytest.raises(ValueError, match="not found"):
        _minimal_id_sky_match(
            query_ids={"ID": QUERY_ID_COL, "MemberID": "does_not_exist"}
        )


def test_id_inconsistent_object_member_counts():
    """A member catalogue missing an object's members is rejected."""
    # query object 1 has no members, so unique member IDs (1) != object IDs (2).
    with pytest.raises(ValueError, match="number of unique 'ID'"):
        _id_sky_match(
            query_member={QUERY_ID_COL: [0, 0, 0], "MemID": ["a", "b", "c"]},
        )


@pytest.mark.parametrize(
    "method",
    [
        SharedFractionMethod.QUERY_PMEM,
        SharedFractionMethod.MATCH_PMEM,
        SharedFractionMethod.PMEM,
    ],
)
def test_match_id_pmem_method_requires_pmem(method):
    """pmem-based shared fractions raise when no pmem column is provided."""
    sky_match = _id_sky_match(with_pmem=False)
    with pytest.raises(ValueError, match="pmem column must"):
        sky_match.match_ID(use_shared_fraction=True, shared_fraction_method=method)


def test_match_id_to_table_best_properties():
    """to_table_best carries requested query and match property columns."""
    sky_match = _id_sky_match()
    result = sky_match.match_ID()
    best = result.select_best()
    table = result.to_table_best(
        best,
        query_properties={"mass_q": "MASS_Q"},
        match_properties={"mass_m": "MASS_M"},
    )
    assert "MASS_Q" in table.colnames
    assert "MASS_M" in table.colnames
    by_query = {row["query_id"]: row for row in table}
    assert by_query[0]["MASS_Q"] == 2.0  # query object 0
    assert by_query[0]["MASS_M"] == 5.0  # matched to match object 100 (row 0)


def test_match_id_to_table_complete_properties():
    """to_table_complete carries requested query and match property columns."""
    sky_match = _id_sky_match()
    result = sky_match.match_ID()
    table = result.to_table_complete(
        query_properties={"mass_q": "MASS_Q"},
        match_properties={"mass_m": "MASS_M"},
    )
    assert "MASS_Q" in table.colnames
    assert "MASS_M" in table.colnames
    by_query = {row["query_id"]: row for row in table}
    assert by_query[0]["MASS_Q"] == 2.0
    # query 0 has a single candidate (match 100, row 0) -> mass 5.0
    assert list(by_query[0]["MASS_M"]) == [5.0]

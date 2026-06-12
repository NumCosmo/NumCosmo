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

"""Distance-based matching: match_2d/3d, filters, select_best, assignment."""

import re
import pytest
import numpy as np
from numpy.testing import assert_allclose

pytest.importorskip("astropy")
# flake8: noqa: E402
# pylint: disable=wrong-import-position

from astropy.table import Table

from numcosmo_py import Nc, Ncm
from numcosmo_py.catalog import (
    SkyMatch,
    SkyMatchResult,
    DistanceMethod,
    Mask,
    SelectionCriteria,
)

Ncm.cfg_init()

ATOL = 1.0e-7
RTOL = 1.0e-5
QUERY_SIZE = 1200


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

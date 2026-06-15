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

"""Catalog loading, coordinate validation and Mask tests."""

import pytest
import numpy as np

pytest.importorskip("astropy")
# flake8: noqa: E402
# pylint: disable=wrong-import-position

from numcosmo_py import Ncm
from numcosmo_py.catalog import SkyMatch, Mask

Ncm.cfg_init()

QUERY_SIZE = 1200
MATCH_SIZE = 1000


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

    assert np.all((mask_all_true & mask_all_true).array == mask_all_true.array)
    assert np.all((~mask_all_false).array == mask_all_true.array)


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

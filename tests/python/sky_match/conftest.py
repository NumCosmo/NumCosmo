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

"""Shared fixtures for the sky_match test suite."""

from pathlib import Path
import pytest
import numpy as np

pytest.importorskip("astropy")
# flake8: noqa: E402
# pylint: disable=wrong-import-position

from astropy.table import Table
from astropy.io import fits

from numcosmo_py import Nc, Ncm
from numcosmo_py.catalog import SkyMatch, DistanceMethod

Ncm.cfg_init()

QUERY_SIZE = 1200
MATCH_SIZE = 1000
QUERY_ID_COL = "ID_q"
MATCH_ID_COL = "ID_m"


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

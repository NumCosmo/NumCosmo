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

"""ID/membership matching: match_ID, shared fractions, select_best, tables."""

import pytest
import numpy as np
from numpy.testing import assert_allclose

pytest.importorskip("astropy")
# flake8: noqa: E402
# pylint: disable=wrong-import-position

from astropy.table import Table

from numcosmo_py import Ncm
from numcosmo_py.catalog import SkyMatch, SharedFractionMethod

Ncm.cfg_init()

QUERY_ID_COL = "ID_q"
MATCH_ID_COL = "ID_m"


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

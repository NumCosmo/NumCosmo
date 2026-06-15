#!/usr/bin/env python
#
# test_halo_catalog.py
#
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
# with this program.  If not, see <http://www.gnu.org/licenses/>.
#

"""Tests for the NcHaloCatalog typed catalog with kind tag and linkage."""

import re
import pytest

from numcosmo_py import Nc, Ncm, GLib

Ncm.cfg_init()


def _cluster_catalog() -> Nc.HaloCatalog:
    """A small cluster catalog linked to halos through ``parent_id``."""
    col_names = ["cluster_id", "parent_id", "Mass"]
    col_types = [
        Ncm.CatalogColType.INT,
        Ncm.CatalogColType.INT,
        Ncm.CatalogColType.DOUBLE,
    ]
    hcat = Nc.HaloCatalog.new(
        Nc.HaloCatalogKind.CLUSTER,
        "cluster_id",
        "parent_id",
        4,
        col_names,
        col_types,
    )
    # Two clusters point at halo 10, one at halo 20, one is a fake (parent 0).
    for i, (cid, pid) in enumerate(((100, 10), (101, 10), (102, 20), (103, 0))):
        hcat.set_int("cluster_id", i, cid)
        hcat.set_int("parent_id", i, pid)
    return hcat


def test_construction_and_metadata() -> None:
    """The catalog exposes its kind and linkage column names."""
    hcat = _cluster_catalog()
    assert isinstance(hcat, Nc.HaloCatalog)
    assert isinstance(hcat, Ncm.Catalog)
    assert hcat.get_kind() == Nc.HaloCatalogKind.CLUSTER
    assert hcat.peek_id_col() == "cluster_id"
    assert hcat.peek_parent_id_col() == "parent_id"
    assert hcat.len() == 4


def test_get_id_and_parent_id() -> None:
    """Row id and parent id are read from the configured columns."""
    hcat = _cluster_catalog()
    assert hcat.get_id(0) == 100
    assert hcat.get_parent_id(0) == 10
    assert hcat.get_id(3) == 103
    assert hcat.get_parent_id(3) == 0


def test_find_children() -> None:
    """find_children returns the row indices matching a parent id."""
    hcat = _cluster_catalog()
    assert list(hcat.find_children(10)) == [0, 1]
    assert list(hcat.find_children(20)) == [2]
    assert list(hcat.find_children(999)) == []


def test_top_level_catalog_without_parent_linkage() -> None:
    """A halo catalog may omit the parent id column."""
    hcat = Nc.HaloCatalog.new(
        Nc.HaloCatalogKind.HALO,
        "halo_id",
        None,
        2,
        ["halo_id", "Mass"],
        [Ncm.CatalogColType.INT, Ncm.CatalogColType.DOUBLE],
    )
    hcat.set_int("halo_id", 0, 10)
    assert hcat.get_kind() == Nc.HaloCatalogKind.HALO
    assert hcat.peek_parent_id_col() is None
    assert hcat.get_id(0) == 10


def test_missing_linkage_raises() -> None:
    """Linkage accessors raise when their column is not configured."""
    hcat = Nc.HaloCatalog.new(
        Nc.HaloCatalogKind.HALO,
        None,
        None,
        1,
        ["Mass"],
        None,
    )

    with pytest.raises(
        GLib.Error,
        match=re.compile(
            rf"^nc-halo-catalog-error: No id column configured. "
            rf"\({int(Nc.HaloCatalogError.NO_LINKAGE)}\)$",
            re.DOTALL,
        ),
    ):
        hcat.get_id(0)

    with pytest.raises(
        GLib.Error,
        match=re.compile(
            rf"^nc-halo-catalog-error: No parent id column configured. "
            rf"\({int(Nc.HaloCatalogError.NO_LINKAGE)}\)$",
            re.DOTALL,
        ),
    ):
        hcat.get_parent_id(0)

    with pytest.raises(
        GLib.Error,
        match=re.compile(
            rf"^nc-halo-catalog-error: No parent id column configured. "
            rf"\({int(Nc.HaloCatalogError.NO_LINKAGE)}\)$",
            re.DOTALL,
        ),
    ):
        hcat.find_children(0)


def test_serialization_roundtrip() -> None:
    """Kind, linkage columns and data survive serialization."""
    hcat = _cluster_catalog()
    hcat.set("Mass", 0, 1.0e14)

    ser = Ncm.Serialize.new(Ncm.SerializeOpt.NONE)
    dup = ser.from_string(ser.to_string(hcat, True))

    assert isinstance(dup, Nc.HaloCatalog)
    assert dup.get_kind() == Nc.HaloCatalogKind.CLUSTER
    assert dup.peek_id_col() == "cluster_id"
    assert dup.peek_parent_id_col() == "parent_id"
    assert dup.get_id(0) == 100
    assert dup.get("Mass", 0) == 1.0e14
    assert list(dup.find_children(10)) == [0, 1]

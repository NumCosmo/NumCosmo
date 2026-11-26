#!/usr/bin/env python
#
# test_table.py
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

"""Tests for the NcmCatalog <-> astropy Table conversions."""

import numpy as np
from astropy.table import Table

from numcosmo_py import Ncm
from numcosmo_py.catalog import catalog_from_table, catalog_to_table
from numcosmo_py.catalog.table import catalog_column, fill_catalog

Ncm.cfg_init()


def test_to_table_preserves_dtypes() -> None:
    """catalog_to_table maps each column type to the matching numpy dtype."""
    types = [
        Ncm.CatalogColType.INT,
        Ncm.CatalogColType.DOUBLE,
        Ncm.CatalogColType.BOOL,
    ]
    catalog = Ncm.Catalog.new_full(2, ["id", "mass", "detected"], types)
    catalog.set_int("id", 0, 7)
    catalog.set_int("id", 1, 42)
    catalog.set("mass", 0, 1.5)
    catalog.set("mass", 1, -3.25)
    catalog.set_bool("detected", 0, True)
    catalog.set_bool("detected", 1, False)

    table = catalog_to_table(catalog)

    assert isinstance(table, Table)
    assert table.colnames == ["id", "mass", "detected"]
    assert np.issubdtype(table["id"].dtype, np.integer)
    assert np.issubdtype(table["mass"].dtype, np.floating)
    assert np.issubdtype(table["detected"].dtype, np.bool_)
    assert list(table["id"]) == [7, 42]
    assert list(table["mass"]) == [1.5, -3.25]
    assert list(table["detected"]) == [True, False]


def test_from_table_infers_types() -> None:
    """catalog_from_table infers column types from the table column dtypes."""
    table = Table()
    table["id"] = np.array([7, 42], dtype=np.int64)
    table["mass"] = np.array([1.5, -3.25], dtype=np.float64)
    table["detected"] = np.array([True, False], dtype=bool)

    catalog = catalog_from_table(table)

    assert catalog.len() == 2
    assert list(catalog.peek_columns()) == ["id", "mass", "detected"]
    assert catalog.get_col_type("id") == Ncm.CatalogColType.INT
    assert catalog.get_col_type("mass") == Ncm.CatalogColType.DOUBLE
    assert catalog.get_col_type("detected") == Ncm.CatalogColType.BOOL
    assert catalog.get_int("id", 1) == 42
    assert catalog.get("mass", 1) == -3.25
    assert catalog.get_bool("detected", 0) is True


def test_table_roundtrip() -> None:
    """A catalog survives a to_table/from_table round-trip unchanged."""
    types = [Ncm.CatalogColType.INT, Ncm.CatalogColType.BOOL]
    catalog = Ncm.Catalog.new_full(3, ["id", "flag"], types)
    catalog.set_int("id", 0, 1)
    catalog.set_int("id", 1, 2)
    catalog.set_int("id", 2, 3)
    catalog.set_bool("flag", 0, True)
    catalog.set_bool("flag", 2, True)

    dup = catalog_from_table(catalog_to_table(catalog))

    assert list(dup.peek_columns()) == ["id", "flag"]
    assert dup.get_col_type("id") == Ncm.CatalogColType.INT
    assert dup.get_col_type("flag") == Ncm.CatalogColType.BOOL
    for i, expected in enumerate((1, 2, 3)):
        assert dup.get_int("id", i) == expected
    assert dup.get_bool("flag", 0) is True
    assert dup.get_bool("flag", 1) is False
    assert dup.get_bool("flag", 2) is True


def test_interleaved_types_exact_ids() -> None:
    """Columns of mixed types keep their values; INT is exact above 2**53."""
    table = Table()
    table["ra"] = np.array([10.0, 20.0, 30.0])
    table["object_id"] = np.array([38558648720430529, 2**62 + 1, -3], dtype=np.int64)
    table["dec"] = np.array([-1.0, 0.0, 1.0])
    table["flag"] = np.array([False, True, False])

    catalog = catalog_from_table(table)

    assert catalog.peek_data().ncols() == 2
    assert catalog_column(catalog, "object_id").tolist() == table["object_id"].tolist()
    assert catalog_column(catalog, "dec").tolist() == [-1.0, 0.0, 1.0]
    assert catalog_column(catalog, "flag").tolist() == [False, True, False]

    back = catalog_to_table(catalog)
    for col in table.colnames:
        assert back[col].tolist() == table[col].tolist()
        assert back[col].dtype == table[col].dtype


def test_fill_catalog_sets_each_type() -> None:
    """fill_catalog writes each column through the storage of its type."""
    types = [Ncm.CatalogColType.BOOL, Ncm.CatalogColType.DOUBLE, Ncm.CatalogColType.INT]
    catalog = Ncm.Catalog.new_full(2, ["b", "x", "n"], types)

    fill_catalog(
        catalog,
        {
            "b": np.array([True, False]),
            "x": np.array([0.5, -0.5]),
            "n": np.array([2**60, 7]),
        },
    )

    assert catalog.get_bool("b", 0) is True
    assert catalog.get("x", 1) == -0.5
    assert catalog.get_int("n", 0) == 2**60

#!/usr/bin/env python
#
# test_py_catalog.py
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

"""Tests for the NcmCatalog named-column container."""

import re
import pytest

from numcosmo_py import Ncm, GLib

Ncm.cfg_init()


def test_construction_defaults_to_zero() -> None:
    """A new catalog has the requested shape and is zero-initialized."""
    catalog = Ncm.Catalog.new(4, ["ra", "dec", "z"])
    assert catalog.len() == 4
    assert catalog.ncols() == 3
    assert list(catalog.peek_columns()) == ["ra", "dec", "z"]
    for col in ("ra", "dec", "z"):
        for i in range(4):
            assert catalog.get(col, i) == 0.0


def test_column_lookup() -> None:
    """Column name to index lookup and membership."""
    catalog = Ncm.Catalog.new(1, ["a", "b", "c"])
    assert catalog.has_column("b")
    assert not catalog.has_column("missing")
    found, idx = catalog.get_index("c")
    assert found and idx == 2
    found, _ = catalog.get_index("missing")
    assert not found


def test_set_get_roundtrip() -> None:
    """Values set by (column, row) are read back at the right place."""
    catalog = Ncm.Catalog.new(3, ["x", "y"])
    catalog.set("x", 0, 1.5)
    catalog.set("y", 2, -3.25)
    assert catalog.get("x", 0) == 1.5
    assert catalog.get("y", 2) == -3.25
    assert catalog.get("y", 0) == 0.0
    # The backing matrix reflects the same values.
    data = catalog.peek_data()
    assert data.nrows() == 3 and data.ncols() == 2
    assert data.get(0, 0) == 1.5
    assert data.get(2, 1) == -3.25


def test_serialization_roundtrip() -> None:
    """Serializing and deserializing preserves columns and data."""
    catalog = Ncm.Catalog.new(2, ["ra", "dec"])
    catalog.set("ra", 0, 10.0)
    catalog.set("dec", 1, 42.0)

    ser = Ncm.Serialize.new(Ncm.SerializeOpt.NONE)
    dup = ser.from_string(ser.to_string(catalog, True))

    assert isinstance(dup, Ncm.Catalog)
    assert dup.len() == 2
    assert list(dup.peek_columns()) == ["ra", "dec"]
    assert dup.get("ra", 0) == 10.0
    assert dup.get("dec", 1) == 42.0


def test_unknown_column_get_index_returns_false() -> None:
    """get_index on an unknown column reports not-found rather than raising."""
    catalog = Ncm.Catalog.new(1, ["only"])
    found, _ = catalog.get_index("nope")
    assert not found


def test_default_column_types_are_double() -> None:
    """Columns created without explicit types default to DOUBLE."""
    catalog = Ncm.Catalog.new(2, ["a", "b"])
    assert catalog.get_col_type("a") == Ncm.CatalogColType.DOUBLE
    assert catalog.get_col_type("b") == Ncm.CatalogColType.DOUBLE


def test_new_full_assigns_column_types() -> None:
    """new_full records the per-column logical types."""
    types = [
        Ncm.CatalogColType.INT,
        Ncm.CatalogColType.DOUBLE,
        Ncm.CatalogColType.BOOL,
    ]
    catalog = Ncm.Catalog.new_full(3, ["id", "mass", "detected"], types)
    assert catalog.get_col_type("id") == Ncm.CatalogColType.INT
    assert catalog.get_col_type("mass") == Ncm.CatalogColType.DOUBLE
    assert catalog.get_col_type("detected") == Ncm.CatalogColType.BOOL


def test_typed_accessors_roundtrip() -> None:
    """Integer and boolean values round-trip through the double backing store."""
    types = [Ncm.CatalogColType.INT, Ncm.CatalogColType.BOOL]
    catalog = Ncm.Catalog.new_full(2, ["id", "flag"], types)

    catalog.set_int("id", 0, 9007199254740991)  # 2**53 - 1, largest exact int.
    catalog.set_bool("flag", 0, True)
    catalog.set_bool("flag", 1, False)

    assert catalog.get_int("id", 0) == 9007199254740991
    assert catalog.get_bool("flag", 0) is True
    assert catalog.get_bool("flag", 1) is False


def test_serialization_preserves_column_types() -> None:
    """Column types survive a serialization round-trip."""
    types = [Ncm.CatalogColType.INT, Ncm.CatalogColType.BOOL]
    catalog = Ncm.Catalog.new_full(1, ["id", "flag"], types)
    catalog.set_int("id", 0, 42)
    catalog.set_bool("flag", 0, True)

    ser = Ncm.Serialize.new(Ncm.SerializeOpt.NONE)
    dup = ser.from_string(ser.to_string(catalog, True))

    assert isinstance(dup, Ncm.Catalog)
    assert dup.get_col_type("id") == Ncm.CatalogColType.INT
    assert dup.get_col_type("flag") == Ncm.CatalogColType.BOOL
    assert dup.get_int("id", 0) == 42
    assert dup.get_bool("flag", 0) is True


def test_unknown_column_get_raises() -> None:
    """get on an unknown column raises a ValueError."""
    catalog = Ncm.Catalog.new(1, ["only"])

    with pytest.raises(
        GLib.Error,
        match=re.compile(
            rf"^ncm-catalog-error: Column 'nope' not found. "
            rf"\({int(Ncm.CatalogError.COLUMN_NOT_FOUND)}\)$",
            re.DOTALL,
        ),
    ):
        catalog.get("nope", 0)

    with pytest.raises(
        GLib.Error,
        match=re.compile(
            rf"^ncm-catalog-error: Column 'nope' not found. "
            rf"\({int(Ncm.CatalogError.COLUMN_NOT_FOUND)}\)$",
            re.DOTALL,
        ),
    ):
        catalog.get_int("nope", 0)

    with pytest.raises(
        GLib.Error,
        match=re.compile(
            rf"^ncm-catalog-error: Column 'nope' not found. "
            rf"\({int(Ncm.CatalogError.COLUMN_NOT_FOUND)}\)$",
            re.DOTALL,
        ),
    ):
        catalog.get_bool("nope", 0)

    with pytest.raises(
        GLib.Error,
        match=re.compile(
            rf"^ncm-catalog-error: Column 'nope' not found. "
            rf"\({int(Ncm.CatalogError.COLUMN_NOT_FOUND)}\)$",
            re.DOTALL,
        ),
    ):
        catalog.get_col_type("nope")


def test_unknown_column_set_raises() -> None:
    """set on an unknown column raises a ValueError."""
    catalog = Ncm.Catalog.new(1, ["only"])

    with pytest.raises(
        GLib.Error,
        match=re.compile(
            rf"^ncm-catalog-error: Column 'nope' not found. "
            rf"\({int(Ncm.CatalogError.COLUMN_NOT_FOUND)}\)$",
            re.DOTALL,
        ),
    ):
        catalog.set("nope", 0, 1.0)

    with pytest.raises(
        GLib.Error,
        match=re.compile(
            rf"^ncm-catalog-error: Column 'nope' not found. "
            rf"\({int(Ncm.CatalogError.COLUMN_NOT_FOUND)}\)$",
            re.DOTALL,
        ),
    ):
        catalog.set_int("nope", 0, 42)

    with pytest.raises(
        GLib.Error,
        match=re.compile(
            rf"^ncm-catalog-error: Column 'nope' not found. "
            rf"\({int(Ncm.CatalogError.COLUMN_NOT_FOUND)}\)$",
            re.DOTALL,
        ),
    ):
        catalog.set_bool("nope", 0, True)

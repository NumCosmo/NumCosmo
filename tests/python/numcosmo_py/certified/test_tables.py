#
# test_tables.py
#
# Thu Sep 24 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_tables.py
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

"""Conversion of the certified tables into NumCosmo binfiles keeps every certificate."""

import gzip
import json
import math
import pathlib
from fractions import Fraction

import pytest

from numcosmo_py import GLib, Ncm
from numcosmo_py.certified import tables

Ncm.cfg_init()

pytestmark = [pytest.mark.xcor]

SOURCES = [
    "truth_tables/xcor/xcor_window_ilk.json.gz",
    "truth_tables/xcor/xcor_window_ilk_d1.json.gz",
    "truth_tables/xcor/xcor_window_ilk_d2.json.gz",
    "truth_tables/xcor/xcor_pair_grid.json.gz",
    "truth_tables/xcor/xcor_pair_grid_mild.json.gz",
    "truth_tables/xcor/xcor_pair_grid_apart.json.gz",
    "truth_tables/xcor/xcor_kquad.json.gz",
]


def _source(name: str) -> pathlib.Path:
    try:
        return pathlib.Path(Ncm.cfg_get_data_filename(name, True))
    except (RuntimeError, GLib.Error):
        pytest.skip(f"{name} is not installed")


def _certified(raw: dict) -> dict[str, list[tuple[str, float]]]:
    """The certified (value, radius) pairs of each group, in the converter's order."""
    if "shapes" in raw or "pair" in raw:
        entries = raw["shapes"] if "shapes" in raw else raw["cases"]

        return {
            name: [
                (v, r)
                for vals, rads in zip(entry["table"], entry["radius"])
                for v, r in zip(vals, rads)
            ]
            for name, entry in entries.items()
        }

    out: dict[str, list[tuple[str, float]]] = {}

    for cell in sorted(raw["cases"].values(), key=lambda c: (c["case"], c["ell"])):
        out.setdefault(cell["case"], []).append((cell["value"], cell["radius"]))

    return out


def test_enclose_keeps_an_exact_double() -> None:
    """A value that is already a double keeps its radius."""
    assert tables.enclose("0.5", 1.0e-30) == (0.5, 1.0e-30)


def test_enclose_widens_by_the_rounding_error() -> None:
    """The stored interval contains the certified one for a value no double holds."""
    exact = "0.0169656654562469581726069850909"
    value, radius = tables.enclose(exact, 7.3991e-39)

    assert value == float(exact)
    assert radius > 7.3991e-39
    assert abs(Fraction(value) - Fraction(exact)) + Fraction(7.3991e-39) <= Fraction(
        radius
    )


@pytest.mark.parametrize("name", SOURCES)
def test_conversion_encloses_every_certified_cell(
    name: str, tmp_path: pathlib.Path
) -> None:
    """Every cell survives as the nearest double inside a still-certified interval."""
    source = _source(name)

    with gzip.open(source, "rt") as handle:
        raw = json.load(handle)

    path = tmp_path / "table.bin"
    tables.write(tables.convert(source), path)
    # the binary form, not the text form g_variant_print writes when binary is off
    assert not path.read_bytes().startswith(b"{")
    back = tables.read(path)
    meta = tables.metadata(back)
    certified = _certified(raw)

    assert meta["source_name"] == source.name
    assert meta["n_excluded"] == 0

    for group, (_, cells) in tables.groups(back).items():
        pairs = certified[group]

        assert len(cells) == len(pairs), group

        for stored, stored_radius, (value, radius) in zip(
            cells["value"], cells["radius"], pairs
        ):
            assert float(stored) == float(Fraction(value))
            assert abs(Fraction(float(stored)) - Fraction(value)) + Fraction(
                radius
            ) <= Fraction(float(stored_radius))
            assert math.isfinite(stored_radius)


def test_window_groups_carry_their_definition(tmp_path: pathlib.Path) -> None:
    """A window group's metadata names its form, parameters and support."""
    path = tmp_path / "window.bin"
    tables.write(tables.convert(_source(SOURCES[0])), path)
    meta, cells = tables.groups(tables.read(path))["multi_mid"]

    assert meta["form"] == "multi"
    assert list(meta["chi_mean"]) == [1000.0, 1600.0]
    assert list(meta["support"]) == [0.0, 2800.0]
    assert sorted(set(cells["ell"])) == [2, 10, 50, 200, 500, 1000]

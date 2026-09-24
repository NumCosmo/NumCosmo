#
# tables.py
#
# Thu Sep 24 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# tables.py
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

"""Certified truth tables as NumCosmo binfiles.

A table is one :class:`~gi.repository.NumCosmoMath.VarDict` serialized with
``ncm_serialize_var_dict_to_variant_file``: the table's own metadata as plain entries,
and one :class:`~gi.repository.NumCosmoMath.Catalog` per group (a window, or a pair of
windows for C_ell) under the key ``cells/<name>``, whose metadata defines the group.
C and Python read it through NumCosmo alone.

Values are stored as doubles. The generators certify each value to about thirty digits;
:func:`enclose` rounds it to the nearest double and widens the radius by the rounding
error, so ``[value - radius, value + radius]`` still provably contains the true value.
The generators' full-precision output is kept outside the repository, and every table
records the SHA-256 of the file it was converted from.

A cell whose value is not a normal double (it underflows) is not stored: no double
computation can resolve it, and against its multipole's peak it is zero. The table
records how many cells were dropped and a bound on their magnitude.
"""

import gzip
import hashlib
import json
import math
import pathlib
import sys
import typing
from fractions import Fraction

import numpy as np
from astropy.table import Table

from numcosmo_py import Ncm, dict_to_var_dict, var_dict_to_dict
from numcosmo_py.catalog.table import catalog_from_table, catalog_to_table
from numcosmo_py.certified.xcor import windows as W

CELLS_PREFIX: typing.Final[str] = "cells/"
MIN_NORMAL: typing.Final[float] = sys.float_info.min


def enclose(value: str, radius: float) -> tuple[float, float]:
    """Return the double nearest a certified value and a radius still enclosing it.

    :param value: the certified midpoint, as the generator's decimal string.
    :param radius: the certified absolute radius.
    :return: the nearest double to ``value`` and ``radius`` plus the rounding error,
        rounded up, so the stored interval contains the certified one.
    """
    # Rational arithmetic is exact for both the decimal string and every double; a
    # fixed-precision decimal would round the bound it is meant to round up.
    exact = Fraction(value)
    rounded = float(exact)
    bound = Fraction(radius) + abs(exact - Fraction(rounded))
    widened = float(bound)

    if Fraction(widened) < bound:
        widened = math.nextafter(widened, math.inf)

    return rounded, widened


def _storable(value: float) -> bool:
    return value == 0.0 or abs(value) >= MIN_NORMAL


def _window_meta(form: str, params: typing.Sequence, deriv: int) -> dict:
    """Metadata defining a window: its form, each parameter by name, the weight."""
    meta: dict[str, typing.Any] = {"form": form, "deriv": deriv}

    for name, value in zip(W.FORM_PARAMS[form], params):
        meta[name] = [float(v) for v in value] if isinstance(value, list) else value

    return meta


class _Builder:
    """Accumulates the groups of one table and the cells it drops."""

    def __init__(self, meta: dict, source: pathlib.Path) -> None:
        self.meta = dict(meta)
        self.meta["source_sha256"] = hashlib.sha256(source.read_bytes()).hexdigest()
        self.meta["source_name"] = source.name
        self.meta["byte_order"] = sys.byteorder
        self.groups: dict[str, tuple[dict, Table]] = {}
        self.n_excluded = 0
        self.excluded_below = 0.0

    def add(self, name: str, meta: dict, rows: dict[str, list]) -> None:
        """Add one group; ``rows`` maps column names to values, ``value`` included."""
        keep = []

        for i, (value, radius) in enumerate(zip(rows["value"], rows["radius"])):
            mid, rad = enclose(value, radius)

            if _storable(mid):
                keep.append((i, mid, rad))
            else:
                self.n_excluded += 1
                self.excluded_below = max(self.excluded_below, abs(mid) + rad)

        table = Table()

        for column, values in rows.items():
            if column == "value":
                table[column] = np.array([mid for _, mid, _ in keep], dtype=np.float64)
            elif column == "radius":
                table[column] = np.array([rad for _, _, rad in keep], dtype=np.float64)
            else:
                table[column] = np.array([values[i] for i, _, _ in keep])

        self.groups[name] = (meta, table)

    def finish(self) -> Ncm.VarDict:
        """Return the table as a VarDict."""
        ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
        self.meta["n_excluded"] = self.n_excluded
        self.meta["excluded_below"] = self.excluded_below
        out = dict_to_var_dict(self.meta)

        for name, (meta, table) in self.groups.items():
            catalog = catalog_from_table(table)
            catalog.set_meta(dict_to_var_dict(meta))
            out.set_object(f"{CELLS_PREFIX}{name}", ser, catalog)

        return out


def _from_window(src: dict, source: pathlib.Path) -> Ncm.VarDict:
    deriv = int(src.get("bessel_deriv", 0))
    targets = src["target_rel"]
    builder = _Builder(
        {
            "kind": "window_ilk",
            "convention": src["convention"],
            "generator": src["generator"],
            "bessel_deriv": deriv,
            "ells": [int(e) for e in src["ells"]],
            "target_rel": targets if isinstance(targets, list) else [targets],
            "x_nu_grid": [float(x) for x in src["x_nu_grid"]],
        },
        source,
    )

    for name, entry in src["shapes"].items():
        meta = _window_meta(entry.get("shape", name), entry["ctor"], deriv)
        meta["chi_ref"] = float(entry["chi_ref"])
        meta["support"] = [float(s) for s in entry["support"]]
        rows: dict[str, list] = {"ell": [], "k": [], "value": [], "radius": []}

        for ell, kvals, values, radii in zip(
            src["ells"], entry["kvals"], entry["table"], entry["radius"]
        ):
            rows["ell"] += [int(ell)] * len(kvals)
            rows["k"] += [float(k) for k in kvals]
            rows["value"] += list(values)
            rows["radius"] += [float(r) for r in radii]

        builder.add(name, meta, rows)

    return builder.finish()


def _from_pair_grid(src: dict, source: pathlib.Path) -> Ncm.VarDict:
    builder = _Builder(
        {
            "kind": "pair_grid",
            "convention": src["convention"],
            "generator": src["generator"],
            "pair_a": src["pair"][0],
            "pair_b": src["pair"][1],
            "ells": [int(e) for e in src["ells"]],
            "target_rel": [float(src["target_rel"])],
        },
        source,
    )

    for name, entry in src["cases"].items():
        meta = _window_meta(entry["form"], entry["ctor"], int(entry["deriv"]))
        rows: dict[str, list] = {"ell": [], "k": [], "value": [], "radius": []}

        for ell, kvals, values, radii in zip(
            src["ells"], src["kvals"], entry["table"], entry["radius"]
        ):
            rows["ell"] += [int(ell)] * len(kvals)
            rows["k"] += [float(k) for k in kvals]
            rows["value"] += list(values)
            rows["radius"] += [float(r) for r in radii]

        builder.add(name, meta, rows)

    return builder.finish()


def _from_kquad(src: dict, source: pathlib.Path) -> Ncm.VarDict:
    powspec = src["powspec"]
    builder = _Builder(
        {
            "kind": "kquad",
            "convention": src["convention"],
            "generator": src["generator"],
            "powspec_shape": powspec["shape"],
            "powspec_amplitude": float(powspec["amplitude"]),
            "powspec_n_s": float(powspec["n_s"]),
            "powspec_k_eq": float(powspec["k_eq"]),
            "ells": [int(e) for e in src["ells"]],
            "target_rel": [float(src["target_rel"])],
        },
        source,
    )
    by_pair: dict[str, list[dict]] = {}

    for cell in src["cases"].values():
        by_pair.setdefault(cell["case"], []).append(cell)

    for name, cells in by_pair.items():
        cells.sort(key=lambda c: c["ell"])
        meta = {
            "kernel_a": cells[0]["kernel_a"],
            "kernel_b": cells[0]["kernel_b"],
            "isauto": bool(cells[0]["isauto"]),
        }
        rows = {
            "ell": [int(c["ell"]) for c in cells],
            "value": [c["value"] for c in cells],
            "radius": [float(c["radius"]) for c in cells],
            "k_lo": [float(c["k_lo"]) for c in cells],
            "k_hi": [float(c["k_hi"]) for c in cells],
            "prec": [int(c["prec"]) for c in cells],
            "seconds": [float(c["seconds"]) for c in cells],
        }
        builder.add(name, meta, rows)

    return builder.finish()


def convert(source: pathlib.Path) -> Ncm.VarDict:
    """Convert a generator's gzipped JSON table into the binfile layout.

    :param source: a ``.json.gz`` written by ``make_xcor_window_truth_table.py``,
        ``make_xcor_pair_grid.py`` or ``make_xcor_kquad_truth_table.py``.
    :return: the table as a VarDict, ready for :func:`write`.
    :raises ValueError: if the file is none of the three.
    """
    with gzip.open(source, "rt") as handle:
        src = json.load(handle)

    if "shapes" in src:
        return _from_window(src, source)

    if "pair" in src:
        return _from_pair_grid(src, source)

    if "powspec" in src:
        return _from_kquad(src, source)

    raise ValueError(f"{source} is not a known truth-table layout")


def write(table: Ncm.VarDict, path: pathlib.Path) -> None:
    """Write a table as a NumCosmo binfile."""
    Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP).var_dict_to_variant_file(
        table, str(path), True
    )


def read(path: pathlib.Path) -> Ncm.VarDict:
    """Read a table written by :func:`write`."""
    return Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP).var_dict_from_variant_file(
        str(path), True
    )


def metadata(table: Ncm.VarDict) -> dict:
    """Return the table's own metadata, without its groups."""
    return {
        key: value
        for key, value in var_dict_to_dict(_plain(table)).items()
        if not key.startswith(CELLS_PREFIX)
    }


def _plain(table: Ncm.VarDict) -> Ncm.VarDict:
    """Return the table without its object entries, which var_dict_to_dict skips."""
    out = Ncm.VarDict.new()

    for key in table.keys():
        if not key.startswith(CELLS_PREFIX):
            found, value = table.get_variant(key)
            assert found
            out.set_variant(key, value)

    return out


def groups(table: Ncm.VarDict) -> dict[str, tuple[dict, Table]]:
    """Return every group of a table: its defining metadata and its cells.

    :param table: a table from :func:`read` or :func:`convert`.
    :return: for each group name, its metadata and an astropy table of its cells.
    """
    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
    out: dict[str, tuple[dict, Table]] = {}

    for key in sorted(table.keys()):
        if key.startswith(CELLS_PREFIX):
            found, catalog = table.get_object(key, ser)
            assert found
            meta = var_dict_to_dict(catalog.peek_meta())
            out[key[len(CELLS_PREFIX) :]] = (meta, catalog_to_table(catalog))

    return out

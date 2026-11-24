#!/usr/bin/env python3
#
# nc_galaxy_pz_update_catalog.py
#
# Tue Nov 17 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# nc_galaxy_pz_update_catalog.py
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

"""Update the p(z) objects of a serialized NcGalaxyWLObs catalog.

A catalog written before NcGalaxyPz holds, in its ``pz`` property, one NcmSpline
per galaxy with the table ``(z_j, p_j)``. This tool replaces each of them by the
NcGalaxyPzBSplineSqrt built from that table and writes a new file. Every other
property of the catalog is copied as it is.

Usage::

    python tools/nc_galaxy_pz_update_catalog.py INPUT.gvar OUTPUT.gvar

The tool edits the serialized GVariant, of type ``(sa{sv})``, and never builds
the NcGalaxyWLObs, whose ``pz`` column does not accept an NcmSpline. An entry
that already is an NcGalaxyPz is copied, so running the tool on an updated
catalog writes the same content.

This tool exists for catalogs written by NumCosmo 0.x and is dropped in 1.0.
"""

from __future__ import annotations

import argparse
import hashlib
import sys
from collections.abc import Sequence
from pathlib import Path

from gi.repository import GLib

from numcosmo_py import Nc, Ncm

# GVariant type of a serialized object: its type name and its properties.
OBJECT_TYPE = "(sa{sv})"

# GVariant type of the serialized ``pz`` property: one object per row index.
PZ_DICT_TYPE = "a{i(sa{sv})}"


def type_name(obj_variant: GLib.Variant) -> str:
    """Return the GType name of a serialized object, without its instance tag."""
    return obj_variant.get_child_value(0).get_string().split("[")[0]


def load_catalog(path: Path) -> GLib.Variant:
    """Read a serialized NcGalaxyWLObs from a binary file of NcmSerialize.

    :raises ValueError: the file holds another type of object.
    """
    data = GLib.Bytes.new(path.read_bytes())
    variant = GLib.Variant.new_from_bytes(GLib.VariantType.new(OBJECT_TYPE), data, True)
    if type_name(variant) != "NcGalaxyWLObs":
        raise ValueError(f"{path} holds a {type_name(variant)}, not a NcGalaxyWLObs")
    return variant


def properties(obj_variant: GLib.Variant) -> list[tuple[str, GLib.Variant]]:
    """Return the names and values of the properties of a serialized object."""
    props = obj_variant.get_child_value(1)
    result = []
    for i in range(props.n_children()):
        entry = props.get_child_value(i)
        name = entry.get_child_value(0).get_string()
        result.append((name, entry.get_child_value(1).get_variant()))
    return result


def update_pz(pz: GLib.Variant) -> tuple[GLib.Variant, int, float]:
    """Replace the NcmSpline entries of a serialized ``pz`` dictionary.

    :return: the new dictionary, the number of entries replaced, and the largest
        difference between a new density at a point of its table and the value
        of the table there, divided by the largest value of that table.
    :raises ValueError: ``pz`` has another GVariant type, or an entry is neither
        an NcmSpline nor an NcGalaxyPz.
    """
    if pz.get_type_string() != PZ_DICT_TYPE:
        raise ValueError(f"pz has type {pz.get_type_string()}, expected {PZ_DICT_TYPE}")

    reader = Ncm.Serialize.new(Ncm.SerializeOpt.KEEP_SHARED)
    writer = Ncm.Serialize.new(Ncm.SerializeOpt.NONE)
    builder = GLib.VariantBuilder.new(GLib.VariantType.new(PZ_DICT_TYPE))
    n_replaced = 0
    max_diff = 0.0

    for i in range(pz.n_children()):
        entry = pz.get_child_value(i)
        key = entry.get_child_value(0)
        obj_variant = entry.get_child_value(1)

        if not type_name(obj_variant).startswith("NcGalaxyPz"):
            spline = reader.from_variant(obj_variant)
            if not isinstance(spline, Ncm.Spline):
                raise ValueError(
                    f"pz entry {key.get_int32()} is a {type_name(obj_variant)}, "
                    "neither a NcmSpline nor a NcGalaxyPz"
                )
            z = spline.peek_xv().dup_array()
            p = spline.peek_yv().dup_array()
            gpz = Nc.GalaxyPzBSplineSqrt.new_from_table(
                spline.peek_xv(), spline.peek_yv()
            )
            diff = max(abs(gpz.eval_p(z_j) - p_j) for z_j, p_j in zip(z, p))
            max_diff = max(max_diff, diff / max(p))
            obj_variant = writer.to_variant(gpz)
            n_replaced += 1

        builder.add_value(GLib.Variant.new_dict_entry(key, obj_variant))

    return builder.end(), n_replaced, max_diff


def replace_pz(catalog: GLib.Variant, new_pz: GLib.Variant) -> GLib.Variant:
    """Return a serialized catalog with the value of its ``pz`` property replaced."""
    builder = GLib.VariantBuilder.new(GLib.VariantType.new("a{sv}"))

    for name, value in properties(catalog):
        builder.add_value(
            GLib.Variant.new_dict_entry(
                GLib.Variant.new_string(name),
                GLib.Variant.new_variant(new_pz if name == "pz" else value),
            )
        )

    # The stubs of pygobject declare one sequence; the override takes the elements.
    return GLib.Variant.new_tuple(  # type: ignore[call-arg]
        catalog.get_child_value(0), builder.end()
    )


def pz_rows(pz: GLib.Variant) -> list[int]:
    """Return the row indices of a serialized ``pz`` dictionary, in file order."""
    return [
        pz.get_child_value(i).get_child_value(0).get_int32()
        for i in range(pz.n_children())
    ]


def is_update_of(new: GLib.Variant, old: GLib.Variant) -> bool:
    """Tell whether ``new`` is ``old`` with its ``pz`` entries updated.

    The object name, the names and order of the properties, and the value of
    every property other than ``pz`` must be equal. The ``pz`` dictionaries must
    have the same rows, and every entry of the new one must be an NcGalaxyPz.
    """
    old_props = properties(old)
    new_props = properties(new)

    if old.get_child_value(0).get_string() != new.get_child_value(0).get_string():
        return False
    if [name for name, _ in old_props] != [name for name, _ in new_props]:
        return False

    for (name, old_value), (_, new_value) in zip(old_props, new_props):
        if name != "pz":
            if not old_value.equal(new_value):
                return False
        elif pz_rows(old_value) != pz_rows(new_value) or not all(
            type_name(new_value.get_child_value(i).get_child_value(1)).startswith(
                "NcGalaxyPz"
            )
            for i in range(new_value.n_children())
        ):
            return False

    return True


def main(argv: Sequence[str] | None = None) -> int:
    """Update one catalog file.

    :param argv: command-line arguments, those of the process when ``None``
    :return: process exit status
    :raises ValueError: the input has no ``pz`` property.
    :raises RuntimeError: the file written is not the input with ``pz`` updated.
    """
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("input", type=Path, help="catalog to read")
    parser.add_argument("output", type=Path, help="catalog to write")
    parser.add_argument(
        "--force", action="store_true", help="overwrite an existing output file"
    )
    args = parser.parse_args(argv)

    if args.output.resolve() == args.input.resolve():
        parser.error("the output must not be the input file")
    if args.output.exists() and not args.force:
        parser.error(f"{args.output} exists; use --force to overwrite it")

    Ncm.cfg_init()

    catalog = load_catalog(args.input)
    old_pz = dict(properties(catalog)).get("pz")
    if old_pz is None:
        raise ValueError("the catalog has no pz property")
    new_pz, n_replaced, max_diff = update_pz(old_pz)
    updated = replace_pz(catalog, new_pz)

    args.output.write_bytes(updated.get_data_as_bytes().get_data())
    if not is_update_of(load_catalog(args.output), catalog):
        raise RuntimeError(f"{args.output} is not {args.input} with pz updated")

    digest = hashlib.sha256(args.output.read_bytes()).hexdigest()
    n_copied = new_pz.n_children() - n_replaced
    print(f"{args.input.name}: {n_replaced} entries replaced, {n_copied} copied")
    print(f"  largest |p(z_j) - p_j| / max_j p_j: {max_diff:.3g}")
    print(f"  size: {args.input.stat().st_size} -> {args.output.stat().st_size} bytes")
    print(f"  sha256: {digest}  {args.output.name}")

    return 0


if __name__ == "__main__":
    sys.exit(main())  # pragma: no cover

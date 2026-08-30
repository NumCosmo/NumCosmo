#!/usr/bin/env python
#
# test_enum_nicks.py
#
# Fri August 29 22:00:00 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_enum_nicks.py
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


"""Tests pinning the introspected nick of every enum and flags member."""

import json

import gi

from numcosmo_py import Ncm

Ncm.cfg_init()


def collect_enum_nicks():
    """Map every enum/flags member of both namespaces to its nick.

    Keyed by C identifier, which is the one name that does not move.

    Read from the typelib rather than from the Python enum classes: whether a
    GI enum is a real enum.Enum, and so whether it exposes ``__members__``,
    depends on the PyGObject version, while the typelib is the same
    everywhere. It also reports the flags members the Python-level view drops.

    The typelib records a member's GIR name, which is its nick with dashes
    written as underscores, plus a trailing underscore when the nick collides
    with a language keyword (``lambda`` is stored as ``lambda_``). Undo both to
    recover the nick itself.
    """
    repo = gi.Repository.get_default()
    nicks = {}

    for namespace, prefix in (("NumCosmoMath", "Ncm"), ("NumCosmo", "Nc")):
        for info in repo.get_infos(namespace):
            if type(info).__name__ != "EnumInfo":
                continue

            for value in info.get_values():
                c_name = value.get_attribute("c:identifier")

                if c_name is None:
                    continue

                name = value.get_name()

                if name.endswith("_"):
                    name = name[:-1]

                nicks[c_name] = {
                    "nick": name.replace("_", "-"),
                    "type": f"{prefix}.{info.get_name()}",
                }

    return nicks


def load_truth_table():
    """The recorded nicks."""
    filename = Ncm.cfg_get_data_filename("truth_tables/enums/enum_nicks.json", True)

    with open(filename, "r", encoding="utf-8") as f:
        return json.load(f)


def test_enum_nicks_match_truth_table():
    """No enum or flags member may change its nick.

    A member's nick is its Python name, its ncm_cfg_get_enum_by_id_name_nick()
    key and what value_nick prints, so renaming one breaks user scripts
    silently. Left implicit, the stripped prefix is whatever the members of an
    enum happen to share -- which means adding a single member can rename every
    other member of that enum. Every enum in the tree therefore pins it with an
    explicit /*< prefix=... >*/ annotation on its typedef, and this test is what
    keeps that true.

    If this fails, the fix is almost never to regenerate the table: it is to
    pin the prefix of the enum whose members moved.
    """
    recorded = load_truth_table()
    current = collect_enum_nicks()

    renamed = {
        c_name: (entry["nick"], current[c_name]["nick"])
        for c_name, entry in recorded.items()
        if c_name in current and current[c_name]["nick"] != entry["nick"]
    }
    assert not renamed, f"enum members changed nick: {renamed}"

    removed = sorted(set(recorded) - set(current))
    assert not removed, f"enum members disappeared: {removed}"


def test_new_enum_members_are_recorded():
    """A new enum member must be added to the truth table deliberately.

    Split from the rename check so that adding a member fails with a message
    about the member you added, not about the table being stale.
    """
    recorded = load_truth_table()
    current = collect_enum_nicks()

    added = sorted(set(current) - set(recorded))
    assert not added, (
        f"enum members not in the truth table: {added}. Check the nick each one "
        "got -- a member that does not start with its enum's pinned prefix keeps "
        "its whole name -- then regenerate with "
        "tests/tools/make_enum_nick_truth_table.py"
    )

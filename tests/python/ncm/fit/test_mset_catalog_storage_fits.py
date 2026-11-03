#!/usr/bin/env python
#
# test_mset_catalog_storage_fits.py
#
# Tue Nov 3 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_mset_catalog_storage_fits.py
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

"""Tests on NcmMSetCatalogStorageFits: catalogs of models created by NcmModelBuilder.

The reader runs in a new process, where the builder types do not exist, so the
types the file stores are the only way to deserialize its model set.
"""

import subprocess
import sys
from pathlib import Path

from numcosmo_py import GObject, Ncm


def _write_builder_catalog(path: Path) -> None:
    """Write a catalog whose model is of a builder type with a builder parent."""
    parent = Ncm.ModelBuilder.new(Ncm.Model, "FitsBuilderParent", "Parent")
    parent.add_sparam("a", "a", 0.0, 1.0, 0.1, 0.0, 0.25, Ncm.ParamType.FREE)
    parent_type = parent.create()

    child = Ncm.ModelBuilder.new(parent_type, "FitsBuilderChild", "Child")
    child.add_sparam("b", "b", 0.0, 1.0, 0.1, 0.0, 0.75, Ncm.ParamType.FREE)
    child_type = child.create()

    GObject.new(child_type)
    mset = Ncm.MSet.new_array([child_type.pytype()])
    mset.param_set_all_ftype(Ncm.ParamType.FREE)
    mset.prepare_fparam_map()

    desc = Ncm.MSetCatalogDesc.new(
        mset, 1, False, ["m2lnL"], ["-2\\ln(L)"], 0, Ncm.MSET_CATALOG_RTYPE_ESMCMC
    )
    storage = Ncm.MSetCatalogStorageFits.new(path.as_posix(), False)
    storage.set_desc(desc)


def _read_in_new_process(path: Path, before: str = "") -> subprocess.CompletedProcess:
    """Run a reader of path in a new process, after the statements of before."""
    script = f"""
from gi.repository import GLib
from numcosmo_py import GObject, Ncm

Ncm.cfg_init()
{before}
try:
    storage = Ncm.MSetCatalogStorageFits.new({path.as_posix()!r}, True)
    model = storage.peek_desc().peek_mset().peek_array_pos(0)
except GLib.Error as error:
    print("error", error.domain, Ncm.ModelBuilderError(error.code).value_nick)
else:
    gtype = model.__gtype__
    print(gtype.name, gtype.parent.name, model.len())
    print(model.param_get(0), model.param_get(1))
"""
    return subprocess.run(
        [sys.executable, "-c", script], capture_output=True, text=True, check=False
    )


def test_builder_model_read_in_new_process(tmp_path: Path) -> None:
    """The builders stored in HDU0 recreate the model types before the model set."""
    path = tmp_path / "builders.fits"
    _write_builder_catalog(path)

    result = _read_in_new_process(path)

    assert result.returncode == 0, result.stderr
    assert result.stdout.split() == [
        "FitsBuilderChild",
        "FitsBuilderParent",
        "2",
        "0.25",
        "0.75",
    ]


def test_builder_name_clash_raises(tmp_path: Path) -> None:
    """A type of the same name but other parameters in the reader is an error."""
    path = tmp_path / "builders.fits"
    _write_builder_catalog(path)

    other = """
other = Ncm.ModelBuilder.new(Ncm.Model, "FitsBuilderParent", "Other")
other.add_sparam("c", "c", 0.0, 1.0, 0.1, 0.0, 0.5, Ncm.ParamType.FREE)
other.create()
"""
    result = _read_in_new_process(path, other)

    assert result.returncode == 0, result.stderr
    assert result.stdout.split() == ["error", "ncm-model-builder-error", "name-clash"]

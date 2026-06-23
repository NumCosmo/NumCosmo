#!/usr/bin/env python
#
# test_mset_func.py
#
# Mon Jun 22 2026
# Copyright  2026  Sandro Dias Pinto Vitenti <vitenti@uel.br>
#
# test_mset_func.py
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

"""Tests for NcmMSetFunc, focusing on the eval-x unique name/symbol identity.

A function bound to a fixed evaluation point (set_eval_x) must carry a *unique*
name/symbol that encodes that point, because the unique name identifies the
function's column in an NcmMSetCatalog. This must survive serialization:
otherwise a single function sampled over a grid (e.g. w(z) on a redshift grid)
deserializes to a set of identically named columns, which collide and make the
catalog impossible to reopen/continue.
"""

from numcosmo_py import Ncm, Nc

Ncm.cfg_init()


def _wz_func(z: float) -> Ncm.MSetFunc:
    """A w_DE(z) function (nvar=1) bound to a fixed redshift z."""
    func = Ncm.MSetFuncList.new("NcHICosmoDE:wDE_z", None)
    func.set_eval_x([z])
    return func


def test_eval_x_uname_is_unique() -> None:
    """Distinct evaluation points give distinct unique names/symbols."""
    base = Ncm.MSetFuncList.new("NcHICosmoDE:wDE_z", None)
    base_name = base.peek_name()

    f1 = _wz_func(0.3)
    f2 = _wz_func(0.7)

    # The bound functions must not fall back to the bare base name.
    assert f1.peek_uname() != base_name
    assert f2.peek_uname() != base_name
    assert f1.peek_uname() != f2.peek_uname()
    assert f1.peek_usymbol() != f2.peek_usymbol()


def test_eval_x_uname_survives_serialization() -> None:
    """The unique name/symbol is rebuilt after (de)serialization.

    Regression test: deserialization sets eval-x through the "eval-x" property,
    which previously did not rebuild the unique name, so the column fell back to
    the (non-unique) base name.
    """
    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)

    func = _wz_func(0.42)
    uname = func.peek_uname()
    usymbol = func.peek_usymbol()

    dup = ser.dup_obj(func)
    assert dup.peek_uname() == uname
    assert dup.peek_usymbol() == usymbol

    # Also via an explicit string round-trip.
    restored = ser.from_string(ser.to_string(func, True))
    assert restored.peek_uname() == uname
    assert restored.peek_usymbol() == usymbol


def test_eval_x_grid_unames_unique_after_serialization() -> None:
    """A function sampled over a grid stays column-unique through serialization."""
    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)

    z_grid = [0.0, 0.25, 0.5, 1.0, 2.33]
    oa = Ncm.ObjArray.new()
    for z in z_grid:
        oa.add(_wz_func(z))

    restored = ser.array_from_yaml(ser.array_to_yaml(oa))
    unames = [restored.get(i).peek_uname() for i in range(restored.len())]

    assert len(set(unames)) == len(z_grid)
    assert all(u != "wDE_z" for u in unames)


if __name__ == "__main__":
    test_eval_x_uname_is_unique()
    test_eval_x_uname_survives_serialization()
    test_eval_x_grid_unames_unique_after_serialization()

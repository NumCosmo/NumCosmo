#!/usr/bin/env python
#
# test_py_obj_array.py
#
# Wed Jan 17 14:34:00 2024
# Copyright  2023  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_py_obj_array.py
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

"""Test NumCosmo NcmObjArray and family."""

from numcosmo_py import Ncm, var_dict_to_dict, dict_to_var_dict


def test_dict_to_var_dict_int():
    """Test dict_to_var_dict."""
    var_dict = dict_to_var_dict({"a": 1, "b": 2})
    assert var_dict.get_int("a") == (True, 1)
    assert var_dict.get_int("b") == (True, 2)
    assert var_dict.get_int("c") == (False, 0)


def test_dict_to_var_dict_double():
    """Test dict_to_var_dict."""
    var_dict = dict_to_var_dict({"a": 1.0, "b": 2.0})
    assert var_dict.get_double("a") == (True, 1.0)
    assert var_dict.get_double("b") == (True, 2.0)
    assert var_dict.get_double("c") == (False, 0.0)


def test_dict_to_var_dict_boolean():
    """Test dict_to_var_dict."""
    var_dict = dict_to_var_dict({"a": True, "b": False})
    assert var_dict.get_boolean("a") == (True, True)
    assert var_dict.get_boolean("b") == (True, False)
    assert var_dict.get_boolean("c") == (False, False)


def test_dict_to_var_dict_str():
    """Test dict_to_var_dict."""
    var_dict = dict_to_var_dict({"a": "1", "b": "2"})
    assert var_dict.get_string("a") == (True, "1")
    assert var_dict.get_string("b") == (True, "2")
    assert var_dict.get_string("c") == (False, None)


def test_dict_to_var_dict_int_array():
    """Test dict_to_var_dict."""
    var_dict = dict_to_var_dict({"a": [1, 2], "b": [3, 4]})
    assert var_dict.get_int_array("a") == (True, [1, 2])
    assert var_dict.get_int_array("b") == (True, [3, 4])
    assert var_dict.get_int_array("c") == (False, [])


def test_dict_to_var_dict_double_array():
    """Test dict_to_var_dict."""
    var_dict = dict_to_var_dict({"a": [1.0, 2.0], "b": [3.0, 4.0]})
    assert var_dict.get_double_array("a") == (True, [1.0, 2.0])
    assert var_dict.get_double_array("b") == (True, [3.0, 4.0])
    assert var_dict.get_double_array("c") == (False, [])


def test_dict_to_var_dict_boolean_array():
    """Test dict_to_var_dict."""
    var_dict = dict_to_var_dict({"a": [True, False], "b": [False, True]})
    assert var_dict.get_boolean_array("a") == (True, [True, False])
    assert var_dict.get_boolean_array("b") == (True, [False, True])
    assert var_dict.get_boolean_array("c") == (False, [])


def test_dict_to_var_dict_mix_all():
    """Test dict_to_var_dict."""
    var_dict = dict_to_var_dict(
        {
            "a": 1,
            "b": 2.0,
            "c": True,
            "d": "1",
            "e": [1, 2],
            "f": [3.0, 4.0],
            "g": [True, False],
        }
    )
    assert var_dict.get_int("a") == (True, 1)
    assert var_dict.get_double("b") == (True, 2.0)
    assert var_dict.get_boolean("c") == (True, True)
    assert var_dict.get_string("d") == (True, "1")
    assert var_dict.get_int_array("e") == (True, [1, 2])
    assert var_dict.get_double_array("f") == (True, [3.0, 4.0])
    assert var_dict.get_boolean_array("g") == (True, [True, False])


def test_var_dict_to_dict_int():
    """Test var_dict_to_dict."""
    var_dict = Ncm.VarDict()
    var_dict.set_int("a", 1)
    var_dict.set_int("b", 2)
    dictionary = var_dict_to_dict(var_dict)
    assert dictionary["a"] == 1
    assert dictionary["b"] == 2
    assert "c" not in dictionary


def test_var_dict_to_dict_double():
    """Test var_dict_to_dict."""
    var_dict = Ncm.VarDict()
    var_dict.set_double("a", 1.0)
    var_dict.set_double("b", 2.0)
    dictionary = var_dict_to_dict(var_dict)
    assert dictionary["a"] == 1.0
    assert dictionary["b"] == 2.0
    assert "c" not in dictionary


def test_var_dict_to_dict_boolean():
    """Test var_dict_to_dict."""
    var_dict = Ncm.VarDict()
    var_dict.set_boolean("a", True)
    var_dict.set_boolean("b", False)
    dictionary = var_dict_to_dict(var_dict)
    assert dictionary["a"] is True
    assert dictionary["b"] is False
    assert "c" not in dictionary


def test_var_dict_to_dict_str():
    """Test var_dict_to_dict."""
    var_dict = Ncm.VarDict()
    var_dict.set_string("a", "1")
    var_dict.set_string("b", "2")
    dictionary = var_dict_to_dict(var_dict)
    assert dictionary["a"] == "1"
    assert dictionary["b"] == "2"
    assert "c" not in dictionary


def test_var_dict_to_dict_int_array():
    """Test var_dict_to_dict."""
    var_dict = Ncm.VarDict()
    var_dict.set_int_array("a", [1, 2])
    var_dict.set_int_array("b", [3, 4])
    dictionary = var_dict_to_dict(var_dict)
    assert dictionary["a"] == [1, 2]
    assert dictionary["b"] == [3, 4]
    assert "c" not in dictionary


def test_var_dict_to_dict_double_array():
    """Test var_dict_to_dict."""
    var_dict = Ncm.VarDict()
    var_dict.set_double_array("a", [1.0, 2.0])
    var_dict.set_double_array("b", [3.0, 4.0])
    dictionary = var_dict_to_dict(var_dict)
    assert dictionary["a"] == [1.0, 2.0]
    assert dictionary["b"] == [3.0, 4.0]
    assert "c" not in dictionary


def test_var_dict_to_dict_boolean_array():
    """Test var_dict_to_dict."""
    var_dict = Ncm.VarDict()
    var_dict.set_boolean_array("a", [True, False])
    var_dict.set_boolean_array("b", [False, True])
    dictionary = var_dict_to_dict(var_dict)
    assert dictionary["a"] == [True, False]
    assert dictionary["b"] == [False, True]
    assert "c" not in dictionary


def test_var_dict_to_dict_mix_all():
    """Test var_dict_to_dict."""
    var_dict = Ncm.VarDict()
    var_dict.set_int("a", 1)
    var_dict.set_double("b", 2.0)
    var_dict.set_boolean("c", True)
    var_dict.set_string("d", "1")
    var_dict.set_int_array("e", [1, 2])
    var_dict.set_double_array("f", [3.0, 4.0])
    var_dict.set_boolean_array("g", [True, False])
    dictionary = var_dict_to_dict(var_dict)
    assert dictionary["a"] == 1
    assert dictionary["b"] == 2.0
    assert dictionary["c"] is True
    assert dictionary["d"] == "1"
    assert dictionary["e"] == [1, 2]
    assert dictionary["f"] == [3.0, 4.0]
    assert dictionary["g"] == [True, False]

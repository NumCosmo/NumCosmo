#
# __init__.py
#
# Wed Feb 8 10:00:00 2023
# Copyright  2023  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# __init__.py
# Copyright (C) 2023 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

"""Import NumCosmo Python bindings."""

from enum import Enum
from typing import Dict, Union, Sequence, cast
import gi

gi.require_version("NumCosmo", "1.0")
gi.require_version("NumCosmoMath", "1.0")

# pyright: reportMissingModuleSource=false
# pylint:disable-next=wrong-import-position,unused-import
from gi.repository import GLib  # noqa: E402

# pylint:disable-next=wrong-import-position,unused-import
from gi.repository import GObject  # noqa: E402

# pylint:disable-next=wrong-import-position,unused-import
from . import ncm as Ncm  # noqa: E402

# pylint:disable-next=wrong-import-position,unused-import
from . import nc as Nc  # noqa: E402

__all__ = [
    "Nc",
    "Ncm",
    "GLib",
    "GObject",
    "GEnum",
    "var_dict_to_dict",
    "dict_to_var_dict",
]


class GEnum(str, Enum):
    """Enum for GObject enums."""

    def __new__(cls, value):
        """Create a new instance of the enum."""
        name = value.value_nick
        obj = str.__new__(cls, name)
        obj._value_ = name
        return obj

    def __init__(self, genum) -> None:
        """Initialize the enum."""
        super().__init__()
        self.genum = genum


def dict_to_var_dict(
    dictionary: Dict[
        str,
        Union[str, float, int, bool, Sequence[float], Sequence[int], Sequence[bool]],
    ],
) -> Ncm.VarDict:
    """Convert a dictionary to a VarDict."""
    var_dict = Ncm.VarDict.new()

    # Since there are types that are subclasses of other types, we need to
    # check the subclasses first, namely bool is a subclass of int.
    for key, value in dictionary.items():
        if isinstance(value, bool):
            var_dict.set_boolean(key, value)
        elif isinstance(value, str):
            var_dict.set_string(key, value)
        elif isinstance(value, float):
            var_dict.set_double(key, value)
        elif isinstance(value, int):
            var_dict.set_int(key, value)
        elif isinstance(value, Sequence):
            if all(isinstance(v, float) for v in value):
                var_dict.set_double_array(key, cast(Sequence[float], value))
            elif all(isinstance(v, bool) for v in value):
                var_dict.set_boolean_array(key, cast(Sequence[bool], value))
            elif all(isinstance(v, int) for v in value):
                var_dict.set_int_array(key, cast(Sequence[int], value))
            else:
                raise TypeError(f"Invalid type for sequence value: {value}")
        else:
            raise TypeError(f"Invalid type for value: {value}")

    return var_dict


def var_dict_to_dict(
    var_dict: Ncm.VarDict,
) -> Dict[
    str, Union[str, float, int, bool, Sequence[float], Sequence[int], Sequence[bool]]
]:
    """Convert a VarDict to a dictionary."""
    dictionary: Dict[
        str,
        Union[str, float, int, bool, Sequence[float], Sequence[int], Sequence[bool]],
    ] = {}

    for key in var_dict.keys():
        found: bool
        value: GLib.Variant
        found, value = var_dict.get_variant(key)

        assert found
        assert isinstance(value, GLib.Variant)

        if value.get_type_string() == "s":
            found, dictionary[key] = var_dict.get_string(key)
            assert found
        elif value.get_type_string() == "d":
            found, dictionary[key] = var_dict.get_double(key)
            assert found
        elif value.get_type_string() == "i":
            found, dictionary[key] = var_dict.get_int(key)
            assert found
        elif value.get_type_string() == "b":
            found, dictionary[key] = var_dict.get_boolean(key)
            assert found
        elif value.get_type_string() == "ad":
            found, dictionary[key] = var_dict.get_double_array(key)
            assert found
        elif value.get_type_string() == "ai":
            found, dictionary[key] = var_dict.get_int_array(key)
            assert found
        elif value.get_type_string() == "ab":
            found, dictionary[key] = var_dict.get_boolean_array(key)
            assert found
        else:
            raise TypeError(
                f"Invalid type for value: {value} variant type "
                f"{value.get_type_string()}"
            )

    return dictionary


def to_camel_case(snake_str):
    """Convert a snake case string to camel case."""
    return "".join(x.capitalize() for x in snake_str.strip().lower().split("_"))

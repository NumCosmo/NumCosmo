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
import gi

gi.require_version("NumCosmo", "1.0")
gi.require_version("NumCosmoMath", "1.0")

# pylint:disable-next=wrong-import-position,unused-import
from gi.repository import GLib  # noqa: E402

# pylint:disable-next=wrong-import-position,unused-import
from gi.repository import GObject  # noqa: E402

# pylint:disable-next=wrong-import-position,unused-import
from . import ncm as Ncm

# pylint:disable-next=wrong-import-position,unused-import
from . import nc as Nc


class GEnum(str, Enum):
    """Enum for GObject enums."""

    def __new__(cls, value):
        name = value.value_nick
        obj = str.__new__(cls, name)
        obj._value_ = name
        return obj

    def __init__(self, genum) -> None:
        super().__init__()
        self.genum = genum


__all__ = ["Nc", "Ncm", "GLib", "GObject", "GEnum"]

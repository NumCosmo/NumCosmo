#
# cosmosis_types.py
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
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program.  If not, see <http://www.gnu.org/licenses/>.

"""Choices shared by the CosmoSIS bridge and the command line.

These live apart from `numcosmo_py.external.cosmosis` because that module imports
cosmosis and firecrown, and firecrown reaches crow, CLMM and healpy. The command line
needs only these names to declare its options, so it takes them from here and leaves
the conversion machinery to be imported when a conversion is actually asked for.
"""

from enum import StrEnum, auto


class LinearMatterPowerSpectrum(StrEnum):
    """Possible linear matter power spectrum models."""

    NONE = auto()
    BBKS = auto()
    EISENSTEIN_HU = auto()
    CLASS = auto()


class NonLinearMatterPowerSpectrum(StrEnum):
    """Possible non-linear matter power spectrum models."""

    NONE = auto()
    HALOFIT = auto()

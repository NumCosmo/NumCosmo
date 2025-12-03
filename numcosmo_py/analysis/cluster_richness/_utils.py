#
# _utils.py
#
# Copyright (C) 2025 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

"""Utility functions for cluster richness analysis.

This module provides helper functions used across the package.
"""

from enum import StrEnum

from numcosmo_py import Nc


#: Format string for parameter display
PARAM_FORMAT = ".3f"
FIXED_PARAMETERS = ["cut"]


class RichnessModelType(StrEnum):
    """Enumeration of supported richness model types."""

    ASCASO = "ascaso"
    EXT = "ext"


def setup_model_fit_params(model: Nc.ClusterMassRichness) -> None:
    """Set up model parameters for fitting with default bounds.

    This function sets the "fit" flag for all parameters except the cut parameter.

    :param model: A NcClusterMassRichness subclass instance
    """
    for i in range(model.sparam_len()):
        name = model.param_name(i)
        if name in FIXED_PARAMETERS:
            continue  # Do not fit the fixed parameters
        model.param_set_desc(name, {"fit": True})


def create_richness_model(
    model_type: str | RichnessModelType = RichnessModelType.ASCASO,
    lnRichness_min: float = 0.0,
    lnRichness_max: float = 20.0,
) -> Nc.ClusterMassRichness:
    """Create a new NcClusterMassRichness model of the specified type.

    :param model_type: Type of model to create ("ascaso" or "ext")
    :param lnRichness_min: Minimum log-richness bound (default: 0.0)
    :param lnRichness_max: Maximum log-richness bound (default: 20.0)
    :return: New model instance
    :raises ValueError: If model_type is not recognized
    """
    match model_type:
        case RichnessModelType.ASCASO:
            return Nc.ClusterMassAscaso(
                lnRichness_min=lnRichness_min, lnRichness_max=lnRichness_max
            )
        case RichnessModelType.EXT:
            return Nc.ClusterMassExt(
                lnRichness_min=lnRichness_min, lnRichness_max=lnRichness_max
            )
        case _:
            raise ValueError(
                f"Unknown model type: {model_type}. Use 'ascaso' or 'ext'."
            )


def get_model_type_name(model: Nc.ClusterMassRichness) -> str:
    """Get the type name of a richness model.

    :param model: A NcClusterMassRichness subclass instance
    :return: Type name string ("ascaso", "ext", or the GObject type name)
    """
    match model:
        case Nc.ClusterMassAscaso():
            return "ascaso"
        case Nc.ClusterMassExt():
            return "ext"
        case _:
            return type(model).__name__

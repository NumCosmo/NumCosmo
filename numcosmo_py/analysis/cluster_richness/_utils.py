#
# _utils.py
#
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

"""Utility functions for cluster richness analysis.

This module provides helper functions used across the package.
"""

from numcosmo_py import Nc


#: Format string for parameter display
PARAM_FORMAT = ".3f"


def setup_model_fit_params(
    model: Nc.ClusterMassRichness,
    bounds: tuple[float, float] = (-20.0, 20.0),
    scale: float = 1.0e-1,
    abstol: float = 1.0e-50,
) -> None:
    """Set up model parameters for fitting with default bounds.

    Configures all parameters in the model with fitting bounds and tolerances.

    :param model: A NcClusterMassRichness subclass instance
    :param bounds: Lower and upper bounds for all parameters (default: (-20.0, 20.0))
    :param scale: Scale for parameter optimization (default: 1.0e-1)
    :param abstol: Absolute tolerance (default: 1.0e-50)
    """
    for i in range(model.sparam_len()):
        name = model.param_name(i)
        value = model.param_get(i)
        model.param_set_desc(
            name,
            {
                "lower-bound": bounds[0],
                "upper-bound": bounds[1],
                "scale": scale,
                "abstol": abstol,
                "fit": True,
                "value": value,
            },
        )


def create_richness_model(
    model_type: str = "ascaso",
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
    model_type_lower = model_type.lower()
    if model_type_lower == "ascaso":
        return Nc.ClusterMassAscaso(
            lnRichness_min=lnRichness_min, lnRichness_max=lnRichness_max
        )
    elif model_type_lower == "ext":
        return Nc.ClusterMassExt(
            lnRichness_min=lnRichness_min, lnRichness_max=lnRichness_max
        )
    else:
        raise ValueError(f"Unknown model type: {model_type}. Use 'ascaso' or 'ext'.")


def get_model_type_name(model: Nc.ClusterMassRichness) -> str:
    """Get the type name of a richness model.

    :param model: A NcClusterMassRichness subclass instance
    :return: Type name string ("ascaso", "ext", or the GObject type name)
    """
    if isinstance(model, Nc.ClusterMassAscaso):
        return "ascaso"
    elif isinstance(model, Nc.ClusterMassExt):
        return "ext"
    else:
        return type(model).__name__

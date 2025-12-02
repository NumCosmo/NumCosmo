#
# _parameters.py
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

"""Parameter dataclasses for the cluster mass-richness analysis.

This module defines the data structures for storing analysis results.
The model parameters are now stored directly in NcClusterMassRichness
subclass instances (NcClusterMassAscaso, NcClusterMassExt, etc.).
"""

from dataclasses import dataclass
from typing import Any, TypeVar

from numcosmo_py import Nc, Ncm


# Type variable for generic model duplication
_T = TypeVar("_T", bound=Ncm.Model)


def dup_model(model: _T) -> _T:
    """Duplicate a NumCosmo model.

    Creates a fresh serializer each time to avoid YAML anchor memory issues.
    The serializer's internal state tracks previously serialized objects to
    produce anchors, but we want independent copies without anchors.

    :param model: A NumCosmo model instance to duplicate
    :return: A duplicate of the model with the same type
    """
    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
    return ser.dup_obj(model)  # type: ignore[return-value]


def model_to_yaml(model: Ncm.Model) -> str:
    """Serialize a NumCosmo model to YAML.

    Creates a fresh serializer each time to avoid YAML anchor memory issues.

    :param model: A NumCosmo model instance to serialize
    :return: YAML string representation of the model
    """
    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
    return ser.to_yaml(model)


def model_from_yaml(yaml_str: str) -> Nc.ClusterMassRichness:
    """Deserialize a NcClusterMassRichness model from YAML.

    Creates a fresh serializer each time to avoid YAML anchor memory issues.

    :param yaml_str: YAML string representation of the model
    :return: Deserialized NcClusterMassRichness model
    """
    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
    return ser.from_yaml(yaml_str)  # type: ignore[return-value]


@dataclass
class CutAnalysisResult:
    """Result from a single cut analysis.

    Contains the cut value, number of clusters passing the cut,
    and best-fit model from different estimation methods.

    The model instances can be duplicated using Ncm.Serialize when needed:
        ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
        model_copy = ser.dup_obj(model)

    Parameters can be accessed via dictionary-like interface:
        model["mup0"] = 1.0
        value = model["mup0"]
    """

    cut: float
    n_clusters: int
    bestfit: Nc.ClusterMassRichness
    mcmc_mean: Nc.ClusterMassRichness
    mcmc_median: Nc.ClusterMassRichness
    bootstrap_mean: Nc.ClusterMassRichness
    m2lnL: float

    def get_bestfit_params_dict(self) -> dict[str, float]:
        """Get best-fit parameters as dictionary.

        :return: Dictionary with parameter names as keys
        """
        return model_params_to_dict(self.bestfit)

    def to_db_row(self, mock_seed: int) -> dict[str, Any]:
        """Convert to database row dictionary for storage.

        :param mock_seed: The mock seed identifier
        :return: Dictionary ready for database insertion
        """
        row: dict[str, Any] = {
            "mock_seed": mock_seed,
            "cut": self.cut,
            "n_clusters": self.n_clusters,
            "m2lnL": self.m2lnL,
        }
        # Add bestfit parameters
        row.update(self.get_bestfit_params_dict())
        return row


def model_params_to_dict(model: Nc.ClusterMassRichness) -> dict[str, float]:
    """Extract all fittable parameters from a model as a dictionary.

    :param model: A NcClusterMassRichness subclass instance
    :return: Dictionary mapping parameter names to values
    """
    result = {}
    sparam_len = model.sparam_len()
    for i in range(sparam_len):
        name = model.param_name(i)
        result[name] = model.param_get(i)
    return result


def model_params_from_dict(
    model: Nc.ClusterMassRichness, params: dict[str, float]
) -> None:
    """Set model parameters from a dictionary.

    :param model: A NcClusterMassRichness subclass instance to modify
    :param params: Dictionary mapping parameter names to values
    """
    for name, value in params.items():
        try:
            model[name] = value
        except Exception:  # noqa: BLE001
            pass  # Skip parameters not in this model


def copy_model_params(
    source: Nc.ClusterMassRichness, target: Nc.ClusterMassRichness
) -> None:
    """Copy parameters from source model to target model.

    Only copies parameters that exist in both models.

    :param source: Source model to copy from
    :param target: Target model to copy to
    """
    params = model_params_to_dict(source)
    model_params_from_dict(target, params)


def model_params_as_list(model: Nc.ClusterMassRichness) -> list[float]:
    """Get model parameters as an ordered list.

    :param model: A NcClusterMassRichness subclass instance
    :return: List of parameter values in model's internal order
    """
    return [model.param_get(i) for i in range(model.sparam_len())]


def model_params_from_list(model: Nc.ClusterMassRichness, values: list[float]) -> None:
    """Set model parameters from an ordered list.

    :param model: A NcClusterMassRichness subclass instance to modify
    :param values: List of parameter values in model's internal order
    """
    for i, value in enumerate(values[: model.sparam_len()]):
        model.param_set(i, value)


def get_model_param_names(model: Nc.ClusterMassRichness) -> list[str]:
    """Get list of parameter names from a model.

    :param model: A NcClusterMassRichness subclass instance
    :return: List of parameter names
    """
    return [model.param_name(i) for i in range(model.sparam_len())]

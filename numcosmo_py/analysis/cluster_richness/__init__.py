#
# __init__.py
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

"""Cluster richness analysis package.

This package provides tools for analyzing cluster mass-richness relations
using NcClusterMassRichness models (Ascaso, Ext, or custom subclasses), including:

- CutAnalyzer for analyzing datasets with progressive richness cuts
- BestfitDatabase for ACID-compliant storage of best-fit models
- MockStudy for mock data generation and analysis
- Truncated normal statistics functions for diagnostics
- Diagnostic plotting functions for model validation

Model Parameters:
    Instead of using dataclasses, this package works directly with
    NcClusterMassRichness subclass instances. Parameters can be accessed
    using dictionary-like syntax:
        model["mup0"] = 1.0
        value = model["mup0"]

    Models can be duplicated using the dup_model function:
        model_copy = dup_model(model)

    Models can be serialized/deserialized using:
        yaml_str = model_to_yaml(model)
        model = model_from_yaml(yaml_str)

Helper Functions:
    - dup_model: Duplicate a NumCosmo model (type-safe)
    - model_to_yaml: Serialize model to YAML string
    - model_from_yaml: Deserialize model from YAML string
    - model_params_to_dict: Extract parameters as dictionary
    - model_params_from_dict: Set parameters from dictionary
    - model_params_as_list: Get parameters as ordered list
    - model_params_from_list: Set parameters from ordered list
    - copy_model_params: Copy parameters between models
    - get_model_param_names: Get list of parameter names
"""

from ._parameters import (
    CutAnalysisResult,
    dup_model,
    model_to_yaml,
    model_from_yaml,
    model_params_to_dict,
    model_params_from_dict,
    copy_model_params,
    model_params_as_list,
    model_params_from_list,
    get_model_param_names,
)
from ._utils import (
    setup_model_fit_params,
    create_richness_model,
    get_model_type_name,
    RichnessModelType,
    PARAM_FORMAT,
)
from ._analyzer import CutAnalyzer
from ._database import BestfitDatabase
from ._mock_study import MockStudy
from ._truncated_stats import (
    mean_lnR_truncated,
    std_lnR_truncated,
    invert_truncated_stats,
    invert_truncated_stats_mu_from_sample,
    invert_truncated_stats_sigma_from_sample,
)
from ._diagnostics import (
    compute_binned_statistics,
    plot_mu_recovery,
    plot_sigma_recovery,
    plot_bin_counts,
    plot_mean_lnR,
    plot_empirical_vs_model_sigma,
    plot_sigma_residuals,
    plot_diagnostic_summary,
)

__all__ = [
    # Results dataclass
    "CutAnalysisResult",
    # Model utilities
    "dup_model",
    "model_to_yaml",
    "model_from_yaml",
    # Model parameter utilities
    "model_params_to_dict",
    "model_params_from_dict",
    "copy_model_params",
    "model_params_as_list",
    "model_params_from_list",
    "get_model_param_names",
    # Utils
    "setup_model_fit_params",
    "create_richness_model",
    "get_model_type_name",
    "RichnessModelType",
    "PARAM_FORMAT",
    # Analysis
    "CutAnalyzer",
    # Database
    "BestfitDatabase",
    # Mock study
    "MockStudy",
    # Truncated statistics
    "mean_lnR_truncated",
    "std_lnR_truncated",
    "invert_truncated_stats",
    "invert_truncated_stats_mu_from_sample",
    "invert_truncated_stats_sigma_from_sample",
    # Diagnostics
    "compute_binned_statistics",
    "plot_mu_recovery",
    "plot_sigma_recovery",
    "plot_bin_counts",
    "plot_mean_lnR",
    "plot_empirical_vs_model_sigma",
    "plot_sigma_residuals",
    "plot_diagnostic_summary",
]

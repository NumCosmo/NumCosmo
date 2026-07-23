#
# __init__.py
#
# Wed Jan 31 08:59:34 2024
# Copyright  2024  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# __init__.py
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

"""NumCosmo command line interface app."""

from typing import TypedDict

import typer

from .from_cosmosis import COSMOSIS
from . import from_cosmosis
from .run_fit import RunFit, RunTest
from .run_mc import RunMC
from .fisher import ComputeTheoryVector, RunFisher, RunFisherBias
from .esmcmc import RunMCMC
from .catalog import (
    AnalyzeMCMC,
    CalibrateCatalog,
    PlotCorner,
    VisualHW,
    ParameterEvolution,
    GetBestFit,
    DumpMset,
)
from .generate import (
    GeneratePlanck,
    GenerateJpasForecast,
    GenerateClusterWL,
    LoadClusterWL,
    GenerateClusterRichnessCount,
    GenerateQSpline,
    GenerateXCDM,
    GenerateDEWSpline,
)
from .cluster_richness import RunClusterRichnessAnalysis
from .inspect import InspectSummary, InspectClusterNCounts
from .xcor import ViewKernel, ListKernels

# Attempt optional import of the Firecrown-NumCosmo connector.
# If available, the import registers the connector in the GObject registry,
# enabling deserialization of models that depend on it (e.g., when loading
# from the best-fit database). Safe to ignore if Firecrown is not installed.
try:
    # pylint:disable-next=wrong-import-position,unused-import
    import firecrown.connector.numcosmo.numcosmo  # noqa: E402, F401
except ImportError:
    pass


app = typer.Typer(no_args_is_help=True, help="NumCosmo command line interface.")
app_run = typer.Typer(no_args_is_help=True, help="Run different statistical analyses.")
app_run_mcmc = typer.Typer(no_args_is_help=True, help="Run MCMC analyses.")
app_cat = typer.Typer(
    no_args_is_help=True, help="MCMC catalog analysis and calibration."
)
app_generate = typer.Typer(no_args_is_help=True, help="Generate experiment files.")
app_analysis = typer.Typer(no_args_is_help=True, help="Data analysis tools.")
app_inspect = typer.Typer(no_args_is_help=True, help="Inspect generated experiments.")
app_xcor = typer.Typer(no_args_is_help=True, help="Cross-correlation kernel analysis.")
app_xcor_kernel = typer.Typer(
    no_args_is_help=True, help="Kernel visualization and analysis."
)

app.add_typer(app_run, name="run")
app_run.add_typer(app_run_mcmc, name="mcmc")
app.add_typer(app_cat, name="catalog")
app.add_typer(app_generate, name="generate")
app.add_typer(app_analysis, name="analysis")
app.add_typer(app_inspect, name="inspect")
app.add_typer(app_xcor, name="xcor")
app_xcor.add_typer(app_xcor_kernel, name="kernel")

CMDArg = TypedDict("CMDArg", {"no_args_is_help": bool, "help": str, "name": str})

FROM_COSMOSIS_CMD: CMDArg = {
    "name": "from-cosmosis",
    "no_args_is_help": True,
    "help": (
        "Converts CosmoSIS .ini files to NumCosmo .yaml experiment files. "
        "This command is only available if CosmoSIS and Firecrown are installed."
    ),
}

RUN_FIT_CMD: CMDArg = {
    "name": "fit",
    "no_args_is_help": True,
    "help": "Computes the best fit of the model to the data.",
}

RUN_MC_CMD: CMDArg = {
    "name": "mc",
    "no_args_is_help": True,
    "help": "Computes the Monte Carlo of the model to the data.",
}

RUN_TEST_CMD: CMDArg = {
    "name": "test",
    "no_args_is_help": True,
    "help": "Loads the experiment file and computes the likelihood once.",
}

RUN_THEORY_VECTOR_CMD: CMDArg = {
    "name": "theory-vector",
    "no_args_is_help": True,
    "help": "Computes theory vector for a given experiment.",
}

RUN_FISHER_CMD: CMDArg = {
    "name": "fisher",
    "no_args_is_help": True,
    "help": "Computes the Fisher matrix of the model to the data.",
}

RUN_FISHER_BIAS_CMD: CMDArg = {
    "name": "fisher-bias",
    "no_args_is_help": True,
    "help": "Computes the Fisher matrix of the model to the data and the bias.",
}

RUN_MCMC_APES_CMD: CMDArg = {
    "name": "apes",
    "no_args_is_help": True,
    "help": "Computes the MCMC using APES.",
}

CAT_ANALYZE_CMD: CMDArg = {
    "name": "analyze",
    "no_args_is_help": True,
    "help": "Analyzes the results of a MCMC run.",
}

CAT_CALIBRATE_CMD: CMDArg = {
    "name": "calibrate",
    "no_args_is_help": True,
    "help": "Calibrate the APES sampler using a given catalog.",
}

CAT_PLOT_CORNER_CMD: CMDArg = {
    "name": "plot-corner",
    "no_args_is_help": True,
    "help": "Plots the corner plot for a given catalog.",
}

CAT_VISUAL_HW_CMD: CMDArg = {
    "name": "visual-hw",
    "no_args_is_help": True,
    "help": "Visualizes the Heidelberger and Welch convergence test.",
}

CAT_PARAM_EVOLUTION_CMD: CMDArg = {
    "name": "param-evolution",
    "no_args_is_help": True,
    "help": "Plots the parameter evolution for a given catalog.",
}

CAT_GET_BEST_FIT_CMD: CMDArg = {
    "name": "get-best-fit",
    "no_args_is_help": True,
    "help": "Get the best fit from a given catalog.",
}

CAT_DUMP_MSET_CMD: CMDArg = {
    "name": "dump-mset",
    "no_args_is_help": True,
    "help": "Dump the model-set stored in a catalog file as YAML.",
}

GEN_PLANCK_CMD: CMDArg = {
    "name": "planck18",
    "no_args_is_help": True,
    "help": "Generate Planck 2018 baseline experiments.",
}

GEN_JPAS_FORECAST_CMD: CMDArg = {
    "name": "jpas-forecast",
    "no_args_is_help": True,
    "help": "Generate JPAS 2024 forecast experiments.",
}

GEN_CLUSTER_WL_CMD: CMDArg = {
    "name": "cluster-wl",
    "no_args_is_help": True,
    "help": "Generate mock cluster weak lensing experiments.",
}

LOAD_CLUSTER_WL_CMD: CMDArg = {
    "name": "cluster-wl-load",
    "no_args_is_help": True,
    "help": "Load a cluster weak lensing experiment from a real NcGalaxyWLObs catalog.",
}

GEN_CLUSTER_RICHNESS_COUNT_CMD: CMDArg = {
    "name": "cluster-richness-count",
    "no_args_is_help": True,
    "help": "Generate mock cluster mass-richness-count (Poisson-Lognormal) experiments.",
}

GEN_QSPLINE_CMD: CMDArg = {
    "name": "qspline",
    "no_args_is_help": True,
    "help": "Generate qspline experiments.",
}

GEN_XCDM_CMD: CMDArg = {
    "name": "xcdm",
    "no_args_is_help": True,
    "help": "Generate xcdm experiments.",
}

GEN_DEWSPLINE_CMD: CMDArg = {
    "name": "de-wspline",
    "no_args_is_help": True,
    "help": "Generate DE w(z) spline experiments.",
}

ANALYSIS_CLUSTER_RICHNESS_CMD: CMDArg = {
    "name": "cluster-richness",
    "no_args_is_help": True,
    "help": "Analyze cluster mass-richness scaling relations.",
}
XCOR_KERNEL_VIEW_CMD: CMDArg = {
    "name": "view",
    "no_args_is_help": True,
    "help": "View and analyze cross-correlation kernels.",
}

XCOR_KERNEL_LIST_CMD: CMDArg = {
    "name": "list",
    "no_args_is_help": False,
    "help": "List all available kernel types and their parameters.",
}

INSPECT_SUMMARY_CMD: CMDArg = {
    "name": "summary",
    "no_args_is_help": True,
    "help": "Inspect experiment summary for likelihood and covariance diagnostics.",
}

INSPECT_CLUSTER_NCOUNTS_CMD: CMDArg = {
    "name": "cluster-ncounts",
    "no_args_is_help": True,
    "help": "Plot data-vector, covariance/correlation, and optional S_ij diagnostics.",
}

# ------------------------------------------------------------------------------
# Installing from-cosmosis command if COSMOSIS is installed and
# all prerequisites are met.
if COSMOSIS:
    app.command(**FROM_COSMOSIS_CMD)(from_cosmosis.numcosmo_from_cosmosis)
# ------------------------------------------------------------------------------
# Installing fitting related subcommands
app_run.command(**RUN_FIT_CMD)(RunFit)
app_run.command(**RUN_MC_CMD)(RunMC)
app_run.command(**RUN_TEST_CMD)(RunTest)
# ------------------------------------------------------------------------------
# Installing Fisher related subcommands
app_run.command(**RUN_THEORY_VECTOR_CMD)(ComputeTheoryVector)
app_run.command(**RUN_FISHER_CMD)(RunFisher)
app_run.command(**RUN_FISHER_BIAS_CMD)(RunFisherBias)
# ------------------------------------------------------------------------------
# Installing MCMC related subcommands
app_run_mcmc.command(**RUN_MCMC_APES_CMD)(RunMCMC)
# ------------------------------------------------------------------------------
# Installing catalog related subcommands
app_cat.command(**CAT_ANALYZE_CMD)(AnalyzeMCMC)
app_cat.command(**CAT_CALIBRATE_CMD)(CalibrateCatalog)
app_cat.command(**CAT_PLOT_CORNER_CMD)(PlotCorner)
app_cat.command(**CAT_VISUAL_HW_CMD)(VisualHW)
app_cat.command(**CAT_PARAM_EVOLUTION_CMD)(ParameterEvolution)
app_cat.command(**CAT_GET_BEST_FIT_CMD)(GetBestFit)
app_cat.command(**CAT_DUMP_MSET_CMD)(DumpMset)
# ------------------------------------------------------------------------------
# Installing experiment generation subcommands
app_generate.command(**GEN_PLANCK_CMD)(GeneratePlanck)
app_generate.command(**GEN_JPAS_FORECAST_CMD)(GenerateJpasForecast)
app_generate.command(**GEN_CLUSTER_WL_CMD)(GenerateClusterWL)
app_generate.command(**LOAD_CLUSTER_WL_CMD)(LoadClusterWL)
app_generate.command(**GEN_CLUSTER_RICHNESS_COUNT_CMD)(GenerateClusterRichnessCount)
app_generate.command(**GEN_QSPLINE_CMD)(GenerateQSpline)
app_generate.command(**GEN_XCDM_CMD)(GenerateXCDM)
app_generate.command(**GEN_DEWSPLINE_CMD)(GenerateDEWSpline)
# ------------------------------------------------------------------------------
# Installing analysis subcommands
app_analysis.command(**ANALYSIS_CLUSTER_RICHNESS_CMD)(RunClusterRichnessAnalysis)
# Installing inspect subcommands
app_inspect.command(**INSPECT_SUMMARY_CMD)(InspectSummary)
app_inspect.command(**INSPECT_CLUSTER_NCOUNTS_CMD)(InspectClusterNCounts)
# Installing xcor kernel subcommands
app_xcor_kernel.command(**XCOR_KERNEL_VIEW_CMD)(ViewKernel)
app_xcor_kernel.command(**XCOR_KERNEL_LIST_CMD)(ListKernels)

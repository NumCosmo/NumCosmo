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
)
from .generate import (
    GeneratePlanck,
    GenerateJpasForecast,
    GenerateClusterWL,
    GenerateQSpline,
    GenerateXCDM,
)

app = typer.Typer(no_args_is_help=True, help="NumCosmo command line interface.")
app_run = typer.Typer(no_args_is_help=True, help="Run different statistical analyses.")
app_run_mcmc = typer.Typer(no_args_is_help=True, help="Run MCMC analyses.")
app_cat = typer.Typer(
    no_args_is_help=True, help="MCMC catalog analysis and calibration."
)
app_generate = typer.Typer(no_args_is_help=True, help="Generate experiment files.")

app.add_typer(app_run, name="run")
app_run.add_typer(app_run_mcmc, name="mcmc")
app.add_typer(app_cat, name="catalog")
app.add_typer(app_generate, name="generate")

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
    "help": "Generate cluster weak lensing experiments.",
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
# ------------------------------------------------------------------------------
# Installing experiment generation subcommands
app_generate.command(**GEN_PLANCK_CMD)(GeneratePlanck)
app_generate.command(**GEN_JPAS_FORECAST_CMD)(GenerateJpasForecast)
app_generate.command(**GEN_CLUSTER_WL_CMD)(GenerateClusterWL)
app_generate.command(**GEN_QSPLINE_CMD)(GenerateQSpline)
app_generate.command(**GEN_XCDM_CMD)(GenerateXCDM)

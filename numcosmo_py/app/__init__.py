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

import dataclasses
from enum import Enum
from typing import Optional, Annotated, Tuple, cast, Union
from pathlib import Path
import typer
from rich.table import Table
from rich.text import Text
import numpy as np

from numcosmo_py import Ncm
from numcosmo_py.sampling import (
    check_runner_algorithm,
    FisherType,
    FitGradType,
    FitRunMessages,
    FitRunner,
    NcmFitLogger,
    set_ncm_console,
)
from numcosmo_py.interpolation.stats_dist import (
    InterpolationMethod,
    InterpolationKernel,
)

try:
    from numcosmo_py.external.cosmosis import (
        convert_likelihoods,
        create_numcosmo_mapping,
        LinearMatterPowerSpectrum,
        NonLinearMatterPowerSpectrum,
    )
except ImportError:
    COSMOSIS = False
else:
    COSMOSIS = True

app = typer.Typer(no_args_is_help=True, help="NumCosmo command line interface.")
app_run = typer.Typer(no_args_is_help=True, help="Run different statistical analyses.")
app_run_mcmc = typer.Typer(no_args_is_help=True, help="Run MCMC analyses.")
app.add_typer(app_run, name="run")
app_run.add_typer(app_run_mcmc, name="mcmc")

if COSMOSIS:
    # Only add the cosmosis commands if cosmosis is available

    def numcosmo_from_cosmosis(
        inifile: Annotated[Path, typer.Argument(help="Path to the Cosmosis ini file.")],
        *,
        outfile: Annotated[
            Optional[Path],
            typer.Option(
                help="Path to the output file, if not given,"
                " the input file name is used with the extension .yaml."
            ),
        ] = None,
        matter_ps: Annotated[
            LinearMatterPowerSpectrum,
            typer.Option(help="Matter power spectrum to use."),
        ] = LinearMatterPowerSpectrum.NONE,
        nonlin_matter_ps: Annotated[
            NonLinearMatterPowerSpectrum,
            typer.Option(help="Non-linear matter power spectrum to use."),
        ] = NonLinearMatterPowerSpectrum.NONE,
        distance_max_z: Annotated[
            float,
            typer.Option(
                help="Max distance to optimize distance computations", min=0.0
            ),
        ] = 10.0,
    ):
        """Converts a Cosmosis ini file to a NumCosmo yaml file, containing
        the same information. The NumCosmo yaml file can be used to run the
        same likelihoods in NumCosmo.

        :param inifile: Path to the Cosmosis ini file.
        :param outfile: Path to the output file, if not given, the input file name
            is used with the extension .yaml.
        :param matter_ps: Matter power spectrum to use.
        :param nonlin_matter_ps: Non-linear matter power spectrum to use.
        :param distance_max_z: Max distance to optimize distance computations.
        """

        Ncm.cfg_init()

        if outfile is None:
            outfile = Path(inifile.stem + ".yaml")

        mapping = create_numcosmo_mapping(
            matter_ps=matter_ps,
            nonlin_matter_ps=nonlin_matter_ps,
            distance_max_z=distance_max_z,
        )

        model_builders, mset, likelihood = convert_likelihoods(inifile, mapping=mapping)

        builders_file = outfile.with_suffix(".builders.yaml")

        experiment = Ncm.ObjDictStr.new()
        experiment.add("likelihood", likelihood)
        experiment.add("model-set", mset)

        ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
        ser.dict_str_to_yaml_file(model_builders, builders_file.absolute().as_posix())
        ser.dict_str_to_yaml_file(experiment, outfile.absolute().as_posix())

    app.command(name="from-cosmosis", no_args_is_help=True)(numcosmo_from_cosmosis)


@dataclasses.dataclass
class LoadExperiment:
    """Common block for commands that load an experiment."""

    experiment: Annotated[
        Path, typer.Argument(help="Path to the experiment file to fit.")
    ]
    starting_point: Annotated[
        Optional[Path],
        typer.Option(
            help=(
                "Path to the file containing the starting point for the fit. "
                "The output of a previous fit can be used."
            ),
        ),
    ] = None
    output: Annotated[
        Optional[Path],
        typer.Option(
            help="Path to the output file, if given, the computed results are written "
            "to this file, otherwise they are not saved."
        ),
    ] = None

    def __post_init__(self) -> None:
        ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)

        builders_file = self.experiment.with_suffix(".builders.yaml")
        # Builders file is optional
        if builders_file.exists():
            model_builders = ser.dict_str_from_yaml_file(
                builders_file.absolute().as_posix()
            )

            for model_builder_name in model_builders.keys():
                model_builder: Ncm.ModelBuilder = cast(
                    Ncm.ModelBuilder, model_builders.get(model_builder_name)
                )
                assert isinstance(model_builder, Ncm.ModelBuilder)
                model_builder.create()

        # We need to initialize NumCosmo after creating the model builders
        # this is necessary because when using MPI, the model builders
        # should be created in all processes before initializing NumCosmo.
        Ncm.cfg_init()
        console = set_ncm_console()

        experiment_objects = ser.dict_str_from_yaml_file(
            self.experiment.absolute().as_posix()
        )

        if experiment_objects.peek("likelihood") is None:
            raise RuntimeError("No likelihood found in experiment file")

        likelihood: Ncm.Likelihood = cast(
            Ncm.Likelihood, experiment_objects.get("likelihood")
        )
        assert isinstance(likelihood, Ncm.Likelihood)

        if experiment_objects.peek("model-set") is None:
            raise RuntimeError("No model-set found in experiment file")

        mset: Ncm.MSet = cast(Ncm.MSet, experiment_objects.get("model-set"))
        assert isinstance(mset, Ncm.MSet)
        mset.prepare_fparam_map()

        if self.starting_point is not None:
            if not self.starting_point.exists():
                raise RuntimeError(
                    f"Starting point file {self.starting_point} not found."
                )

            ser.reset(False)
            starting_dict = ser.dict_str_from_yaml_file(
                self.starting_point.absolute().as_posix()
            )
            if starting_dict.peek("model-set") is None:
                raise RuntimeError(
                    f"Starting point file {self.starting_point} does not contain "
                    f"a model-set."
                )
            saved_mset: Ncm.MSet = cast(Ncm.MSet, starting_dict.get("model-set"))
            assert isinstance(saved_mset, Ncm.MSet)
            assert isinstance(saved_mset, Ncm.MSet)
            if not mset.cmp(saved_mset, True):
                raise RuntimeError(
                    f"Starting point file {self.starting_point} "
                    f"does not match experiment."
                )
            mset.param_set_mset(saved_mset)

        self.console = console
        self.likelihood = likelihood
        self.mset = mset

        if self.output is not None:
            if self.output.exists():
                ser.reset(False)
                self.output_dict = ser.dict_str_from_yaml_file(
                    self.output.absolute().as_posix()
                )
            else:
                self.output_dict = Ncm.ObjDictStr.new()

    def end_experiment(self):
        """Ends the experiment and writes the output file."""
        if self.output is not None:
            ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
            ser.dict_str_to_yaml_file(
                self.output_dict, self.output.absolute().as_posix()
            )


@app_run.command(
    name="theory-vector",
    help="Computes theory vectory for a given experiment.",
    no_args_is_help=True,
)
@dataclasses.dataclass
class ComputeTheoryVector(LoadExperiment):
    """Computes theory vectory for a given experiment."""

    def __post_init__(self) -> None:
        super().__post_init__()

        dset: Ncm.Dataset = self.likelihood.peek_dataset()
        if not dset.has_mean_vector():
            raise RuntimeError("mean vector computation not supported by this dataset.")

        theory_vector = Ncm.Vector.new(dset.get_n())
        dset.mean_vector(self.mset, theory_vector)
        if self.output is not None:
            self.output_dict.add("theory-vector", theory_vector)
        else:
            self.console.print(theory_vector.dup_array())

        self.end_experiment()


@dataclasses.dataclass
class RunCommonOptions(LoadExperiment):
    """Common options for the run command."""

    runner: Annotated[
        FitRunner,
        typer.Option(
            help="Algorithm to use for the fit.",
        ),
    ] = FitRunner.NLOPT.value
    algorithm: Annotated[
        Optional[str],
        typer.Option(
            help="Algorithm to use for the fit.",
        ),
    ] = None
    grad_type: Annotated[
        FitGradType,
        typer.Option(
            help="Gradient type to use for the fit.",
        ),
    ] = FitGradType.NUMDIFF_FORWARD.value
    run_messages: Annotated[
        FitRunMessages,
        typer.Option(
            help="Verbosity level for the fit.",
        ),
    ] = FitRunMessages.SIMPLE.value

    def __post_init__(self) -> None:
        super().__post_init__()

        check_runner_algorithm(self.runner, self.algorithm)

        fit = Ncm.Fit.factory(
            self.runner.genum,
            self.algorithm,
            self.likelihood,
            self.mset,
            self.grad_type.genum,
        )
        fit.set_messages(self.run_messages.genum)

        fit_logger = NcmFitLogger(self.console)
        fit.set_logger(
            fit_logger.write_progress,
            fit_logger.update_progress,
            fit_logger.start_update,
            fit_logger.end_update,
        )
        self.fit = fit


@app_run.command(
    name="fit",
    help="Computes the best fit of the model to the data.",
    no_args_is_help=True,
)
@dataclasses.dataclass
class RunFit(RunCommonOptions):
    """Computes the best fit of the model to the data."""

    restart: Annotated[
        Tuple[float, float],
        typer.Option(
            help=(
                "Restart the fit until the given the value of m2lnL varies less"
                " than the given tolerance (abstol, reltol)."
            ),
        ),
    ] = (
        None,
        None,
    )  # type: ignore

    def __post_init__(self) -> None:
        super().__post_init__()
        self.fit.log_info()

        abstol, reltol = self.restart

        if abstol is None or reltol is None:
            self.fit.run(self.run_messages.genum)
        else:
            if abstol <= 0.0 and reltol <= 0.0:
                raise RuntimeError(f"Invalid tolerance for restart {self.restart}.")
            output_filename = (
                None
                if self.output is None
                else self.output.with_suffix(".tmp").absolute().as_posix()
            )
            self.fit.run_restart(
                self.run_messages.genum,
                abstol,
                reltol,
                None,
                output_filename,
            )

        if self.output is not None:
            self.output_dict.add("model-set", self.fit.peek_mset())

        self.end_experiment()


@app_run.command(
    name="test",
    help="Loads the experiment file and computes the likelihood once.",
    no_args_is_help=True,
)
@dataclasses.dataclass
class RunTest(RunCommonOptions):
    """Loads the experiment file and computes the likelihood once."""

    def __post_init__(self) -> None:
        super().__post_init__()
        self.mset.param_set_all_ftype(Ncm.ParamType.FIXED)
        self.mset.prepare_fparam_map()
        self.fit.log_info()
        self.fit.run(FitRunMessages.SIMPLE.genum)


@app_run.command(
    name="fisher",
    help="Computes the Fisher matrix of the model to the data.",
    no_args_is_help=True,
)
@dataclasses.dataclass
class RunFisher(RunCommonOptions):
    """Computes the Fisher matrix of the model to the data."""

    fisher_type: Annotated[
        FisherType,
        typer.Option(
            help="Type of Fisher matrix to compute.",
        ),
    ] = FisherType.OBSERVED

    def __post_init__(self) -> None:
        super().__post_init__()
        if self.fisher_type == FisherType.OBSERVED:
            self.fit.obs_fisher()
            self.fit.log_covar()
        elif self.fisher_type == FisherType.EXPECTED:
            self.fit.fisher()
            self.fit.log_covar()
        else:
            raise RuntimeError(f"Invalid Fisher type {self.fisher_type}.")

        if self.output is not None:
            self.output_dict.add("model-set", self.fit.peek_mset())
            self.output_dict.add("covariance", self.fit.get_covar())

        self.end_experiment()


@app_run.command(
    name="fisher-bias",
    help="Computes the Fisher matrix of the model to the data and the bias.",
    no_args_is_help=True,
)
@dataclasses.dataclass
class RunFisherBias(RunCommonOptions):
    """Computes the Fisher matrix of the model to the data and the bias."""

    theory_vector: Annotated[
        Optional[Path],
        typer.Option(
            help="Path to the theory vector file to compute the bias relative to."
        ),
    ] = None

    def __post_init__(self) -> None:
        super().__post_init__()

        if self.theory_vector is None:
            raise RuntimeError("No theory vector file given.")

        ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
        theory_vector_dict = ser.dict_str_from_yaml_file(
            self.theory_vector.absolute().as_posix()
        )

        dset = self.likelihood.peek_dataset()
        if not dset.has_mean_vector():
            raise RuntimeError(
                "mean vector computation not supported by this dataset, "
                "cannot compute bias."
            )

        theory_vector = theory_vector_dict.get("theory-vector")
        if not isinstance(theory_vector, Ncm.Vector):
            raise RuntimeError("Invalid theory vector file.")

        if theory_vector.len() != dset.get_n():
            raise RuntimeError(
                "Theory vector and dataset have different sizes, "
                "cannot compute bias."
            )

        delta_theta = self.fit.fisher_bias(theory_vector)

        if self.output is not None:
            self.output_dict.add("covariance", self.fit.get_covar())
            self.output_dict.add("delta-theta", delta_theta)
        else:
            self.console.print(delta_theta.dup_array())

        self.end_experiment()


class IniSampler(str, Enum):
    """Sampler to use for the MCMC."""

    GAUSS_MSET = "gauss-mset"
    GAUSS_COV = "gauss-cov"
    FROM_CATALOG = "from-catalog"


class Parallezation(str, Enum):
    """Parallel sampler to use for the MCMC."""

    NONE = "none"
    MPI = "mpi"
    THREADS = "threads"


@app_run_mcmc.command(
    name="apes",
    help="Computes the MCMC using APES.",
    no_args_is_help=True,
)
@dataclasses.dataclass
class RunMCMC(RunCommonOptions):
    """Computes the MCMC using APES."""

    nwalkers: Annotated[
        int,
        typer.Option(
            help="Number of walkers to use, 0 means automatic.",
            min=0,
        ),
    ] = 0

    nsamples: Annotated[
        int,
        typer.Option(
            help="Number of samples to compute.",
            min=1,
        ),
    ] = 10

    robust: Annotated[
        bool,
        typer.Option(
            help="Use robust covariance estimation. (Experimental)",
        ),
    ] = False

    interpolation_method: Annotated[
        InterpolationMethod,
        typer.Option(
            help="Interpolation method to use.",
        ),
    ] = InterpolationMethod.VKDE.value

    interpolation_kernel: Annotated[
        InterpolationKernel,
        typer.Option(
            help="Interpolation kernel to use.",
        ),
    ] = InterpolationKernel.CAUCHY.value

    over_smooth: Annotated[
        float,
        typer.Option(
            help="Over-smoothing parameter to use.",
            min=1.0e-2,
        ),
    ] = 0.2

    local_fraction: Annotated[
        Optional[float],
        typer.Option(
            help="Local fraction to use.",
            min=0.02,
        ),
    ] = None

    use_interpolation: Annotated[
        bool,
        typer.Option(
            help="Use interpolation to compute the weights of the APES approximation.",
        ),
    ] = True

    parallel: Annotated[
        Parallezation,
        typer.Option(
            help=(
                "Parallelization to use. Python likelihoods are not compatible with "
                "multi-threading."
            ),
        ),
    ] = Parallezation.NONE

    nthreads: Annotated[
        int,
        typer.Option(
            help="Number of threads to use when using multi-threading.",
            min=2,
        ),
    ] = 4

    initial_points_sampler: Annotated[
        Optional[IniSampler],
        typer.Option(
            help=(
                "Sampler to use for the initial points. "
                "Gaussian based samplers use a gaussian approximation of the "
                "posterior. The covariance matrix can be computed using the "
                "scales in the model-set or using the covariance matrix of the last "
                "fit. The catalog sampler uses a catalog of points to sample the "
                "initial points."
            ),
        ),
    ] = IniSampler.GAUSS_MSET

    initial_sampler_rescale: Annotated[
        float,
        typer.Option(
            help=(
                "Rescale factor for the Gaussian based initial samplers. "
                "The rescale factor is applied to the covariance matrix or to the "
                "scales in the model-set."
            ),
            min=0.01,
        ),
    ] = 1.0

    initial_sampler_covar: Annotated[
        Optional[Path],
        typer.Option(
            help=(
                "Path to the covariance matrix file to use for the initial points "
                "sampler."
            ),
        ),
    ] = None

    initial_catalog: Annotated[
        Optional[Path],
        typer.Option(
            help="Path to the catalog file to use for the initial points sampler.",
        ),
    ] = None

    initial_catalog_burnin: Annotated[
        int,
        typer.Option(
            help="Number of points to discard from the catalog.",
            min=0,
        ),
    ] = 0

    def __post_init__(self) -> None:
        super().__post_init__()
        self.fit.log_info()

        fparams_len = self.mset.fparams_len()

        if self.nwalkers == 0:
            self.nwalkers = 100 * (fparams_len + 1)

        init_sampler: Union[Ncm.MSetTransKernGauss, Ncm.MSetTransKernCat]
        if self.initial_points_sampler == IniSampler.GAUSS_MSET:
            init_sampler = Ncm.MSetTransKernGauss.new(0)
            init_sampler.set_mset(self.mset)
            init_sampler.set_prior_from_mset()
            init_sampler.set_cov_from_rescale(self.initial_sampler_rescale)
        elif self.initial_points_sampler == IniSampler.GAUSS_COV:
            init_sampler = Ncm.MSetTransKernGauss.new(0)
            init_sampler.set_mset(self.mset)
            init_sampler.set_prior_from_mset()
            if self.initial_sampler_covar is None:
                raise RuntimeError("No covariance file given.")

            ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
            cov_dict = ser.dict_str_from_yaml_file(
                self.initial_sampler_covar.absolute().as_posix()
            )
            cov = cov_dict.get("covariance")
            if not isinstance(cov, Ncm.Matrix):
                raise RuntimeError("Invalid covariance file.")
            cov.scale(self.initial_sampler_rescale)
            init_sampler.set_cov(cov)
        elif self.initial_points_sampler == IniSampler.FROM_CATALOG:
            if self.initial_catalog is None:
                raise RuntimeError("No catalog file given.")
            if not self.initial_catalog.exists():
                raise RuntimeError(f"Catalog file {self.initial_catalog} not found.")

            mcat = Ncm.MSetCatalog.new_from_file_ro(
                self.initial_catalog.absolute().as_posix(),
                self.initial_catalog_burnin,
            )
            init_sampler = Ncm.MSetTransKernCat.new(mcat, None)

            init_sampler.set_sampling(Ncm.MSetTransKernCatSampling.CHOOSE)
            init_sampler.set_mset(self.mset)
            init_sampler.set_prior_from_mset()
        else:
            raise RuntimeError(
                f"Invalid initial points sampler {self.initial_points_sampler}."
            )

        apes_walker = Ncm.FitESMCMCWalkerAPES.new(self.nwalkers, fparams_len)
        apes_walker.set_over_smooth(self.over_smooth)
        if self.local_fraction is not None:
            apes_walker.set_local_frac(self.local_fraction)

        apes_walker.use_interp(self.use_interpolation)
        apes_walker.set_method(self.interpolation_method.genum)
        apes_walker.set_k_type(self.interpolation_kernel.genum)

        if self.parallel == Parallezation.THREADS.value:
            apes_walker.set_use_threads(True)
        elif self.parallel == Parallezation.MPI.value:
            if Ncm.cfg_mpi_nslaves() == 0:
                raise RuntimeError(
                    "MPI parallelization requested but MPI is not initialized."
                )
        else:
            apes_walker.set_use_threads(False)

        esmcmc = Ncm.FitESMCMC.new(
            self.fit,
            self.nwalkers,
            init_sampler,
            apes_walker,
            self.run_messages.genum,
        )
        if self.parallel == Parallezation.THREADS.value:
            esmcmc.set_nthreads(self.nthreads)

        if self.output is not None:
            esmcmc.set_data_file(
                self.output.absolute().with_suffix(".mcmc.fits").as_posix()
            )

        esmcmc.start_run()
        esmcmc.run(self.nsamples)
        esmcmc.end_run()

        self.end_experiment()


@app.command(
    name="analyze",
    help="Analyzes the results of a MCMC run.",
    no_args_is_help=True,
)
@dataclasses.dataclass
class AnalyzeMCMC(LoadExperiment):
    """Analyzes the results of a MCMC run."""

    mcmc_file: Annotated[
        Optional[Path],
        typer.Argument(
            help="Path to the MCMC file.",
        ),
    ] = None

    burnin: Annotated[
        int,
        typer.Option(
            help="Number of samples to discard as burnin.",
            min=0,
        ),
    ] = 0

    info: Annotated[
        bool,
        typer.Option(
            help="Prints information about the MCMC file.",
        ),
    ] = True

    def __post_init__(self) -> None:
        super().__post_init__()

        if self.mcmc_file is None:
            raise RuntimeError("No MCMC file given.")

        if not self.mcmc_file.exists():
            raise RuntimeError(f"MCMC file {self.mcmc_file} not found.")

        mcat = Ncm.MSetCatalog.new_from_file_ro(
            self.mcmc_file.absolute().as_posix(), self.burnin
        )
        mset = mcat.peek_mset()
        mset.prepare_fparam_map()
        fparams_len = mset.fparams_len()
        nadd_vals = mcat.nadd_vals()
        total_columns = fparams_len + nadd_vals
        nchains = mcat.nchains()
        full_stats = mcat.peek_pstats()
        if nchains > 1:
            stats = mcat.peek_e_mean_stats()
        else:
            stats = mcat.peek_pstats()
        nitems = stats.nitens()
        if nitems >= 10:
            mcat.estimate_autocorrelation_tau(False)

        desc_color = "bold bright_cyan"
        values_color = "bold bright_green"
        main_table = Table(title="Catalog information")
        main_table.show_header = False

        main_table.add_column(justify="left")

        details = Table(title="Run details", expand=False)
        details.show_header = False
        details.add_column(justify="left", style=desc_color)
        details.add_column(justify="right", style=values_color)

        details.add_row("Run type", mcat.get_run_type())
        details.add_row("Size", f"{mcat.len()}")
        details.add_row("Number of Iterations", f"{mcat.max_time()}")
        details.add_row("Number of chains", f"{nchains}")
        details.add_row("Number of parameters", f"{fparams_len}")
        details.add_row("Number of extra columns", f"{nadd_vals}")
        details.add_row("Weighted", f"{mcat.weighted()}")
        main_table.add_row(details)

        # Global diagnostics

        global_diag = Table(
            title="Global Convergence Diagnostics",
            expand=True,
        )
        global_diag.add_column("Diagnostic Statistic", justify="left", style=desc_color)
        global_diag.add_column("Suggested cut-off", justify="left", style=values_color)
        global_diag.add_column("Worst parameter", justify="left", style=values_color)
        global_diag.add_column("AR model order", justify="left", style=values_color)
        global_diag.add_column("Value", justify="left", style=values_color)

        param_diag = Table(title="Parameters", expand=False, show_lines=True)
        param_diag_matrix = []

        # Parameter names
        param_diag.add_column(
            "Parameter", justify="left", style=desc_color, vertical="middle"
        )
        param_diag_matrix.append([mcat.col_full_name(i) for i in range(total_columns)])

        # Values color
        val_color = values_color
        # Parameter best fit
        best_fit_vec = mcat.get_bestfit_row()
        param_diag.add_column(
            "Best-fit", justify="left", style=val_color, vertical="middle"
        )
        param_diag_matrix.append(
            [f"{best_fit_vec.get(i): .6g}" for i in range(total_columns)]
        )

        # Parameter mean
        param_diag.add_column(
            "Mean", justify="left", style=val_color, vertical="middle"
        )
        param_diag_matrix.append(
            [f"{full_stats.get_mean(i): .6g}" for i in range(total_columns)]
        )

        # Standard Deviation

        param_diag.add_column(
            "Standard Deviation", justify="left", style=val_color, vertical="middle"
        )
        param_diag_matrix.append(
            [f"{full_stats.get_sd(i): .6g}" for i in range(total_columns)]
        )

        if nitems >= 10:
            # Mean Standard Deviation
            param_diag.add_column(
                "Mean Standard Deviation",
                justify="left",
                style=val_color,
                vertical="middle",
            )
            tau_vec = mcat.peek_autocorrelation_tau()
            mean_sd_array = [
                np.sqrt(full_stats.get_var(i) * tau_vec.get(i) / full_stats.nitens())
                for i in range(total_columns)
            ]
            param_diag_matrix.append([f"{mean_sd: .6g}" for mean_sd in mean_sd_array])

            # Autocorrelation Time
            tau_row = []
            tau_row.append("Autocorrelation time (tau)")
            tau_row.append("NA")
            tau_row.append(
                f"{tau_vec.get_max():.0f} ({mcat.col_full_name(tau_vec.get_max_index())})"
            )
            tau_row.append("NA")
            tau_row.append(f"{tau_vec.get_max():.3f}")
            global_diag.add_row(*tau_row)

            param_diag.add_column(
                "tau", justify="left", style=val_color, vertical="middle"
            )
            param_diag_matrix.append(
                [f"{tau_vec.get(i): .6g}" for i in range(total_columns)]
            )

        if nchains > 1:
            # Gelman Rubin
            gelman_rubin_row = []
            gelman_rubin_row.append("Gelman-Rubin (G&B) Shrink Factor (R-1)")
            skf = [mcat.get_param_shrink_factor(i) - 1 for i in range(total_columns)]
            gelman_rubin_row.append("NA")
            gr_worst = int(np.argmin(skf))
            gelman_rubin_row.append(
                f"{skf[gr_worst]:.3f} ({mcat.col_full_name(gr_worst)})"
            )
            gelman_rubin_row.append("NA")
            gelman_rubin_row.append(f"{mcat.get_shrink_factor() - 1:.3f}")
            global_diag.add_row(*gelman_rubin_row)

            param_diag.add_column(
                "G&R", justify="left", style=val_color, vertical="middle"
            )
            param_diag_matrix.append([f"{skf_i:.3f}" for skf_i in skf])

        # Constant Break

        cb = [stats.estimate_const_break(i) for i in range(total_columns)]
        cb_worst = int(np.argmax(cb))
        const_break_row = []
        const_break_row.append("Constant Break (CB) (iterations, points)")
        const_break_row.append(f"{cb[cb_worst]:.0f}")
        const_break_row.append(f"{cb[cb_worst]:.0f} ({mcat.col_full_name(cb_worst)})")
        const_break_row.append("NA")
        const_break_row.append(f"{cb[cb_worst]:.0f}")
        global_diag.add_row(*const_break_row)

        param_diag.add_column("CB", justify="left", style=val_color)
        param_diag_matrix.append([f"{cb_i:.0f} {cb_i * nchains:.0f}" for cb_i in cb])

        if nitems >= 10:
            # Effective sample size
            (
                ess_vec,
                ess_best_cutoff,
                ess_worst_index,
                ess_worst_order,
                ess_worst_ess,
            ) = stats.max_ess_time(100)
            ess_row = []
            ess_row.append("Effective Sample Size (ESS) (ensembles, points)")
            ess_row.append(f"{ess_best_cutoff}")
            ess_row.append(
                f"{ess_vec.get(ess_worst_index):.0f} ({mcat.col_full_name(ess_worst_index)})"
            )
            ess_row.append(f"{ess_worst_order}")
            ess_row.append(f"{ess_worst_ess:.0f}")
            global_diag.add_row(*ess_row)

            param_diag.add_column("ESS", justify="left", style=val_color)
            param_diag_matrix.append(
                [
                    f"{ess_vec.get(i):.0f} {ess_vec.get(i) * nchains:.0f}"
                    for i in range(total_columns)
                ]
            )

            # Heidelberger and Welch

            hw_pvalue = 1.0 - 0.95 ** (1.0 / fparams_len)
            (
                hw_vec,
                hw_best_cutoff,
                hw_worst_index,
                hw_worst_order,
                hw_worst_pvalue,
            ) = stats.heidel_diag(100, hw_pvalue)

            hw_row = []
            hw_row.append(f"Heidelberger and Welch p-value (>{hw_pvalue * 100.0:.1f}%)")

            if hw_best_cutoff >= 0:
                hw_row.append(f"{hw_best_cutoff}")
            else:
                hw_row.append("All tests fail")
            hw_row.append(
                f"{(1.0 - hw_worst_pvalue) * 100.0:.1f}% ({mcat.col_full_name(hw_worst_index)})"
            )
            hw_row.append(f"{hw_worst_order}")
            hw_row.append(f"{(1.0 - hw_worst_pvalue) * 100.0:.1f}%")
            global_diag.add_row(*hw_row)

            param_diag.add_column(
                "H&W",
                justify="left",
                style=val_color,
            )
            param_diag_matrix.append(
                [f"{(1.0 - hw_vec.get(i)) * 100.0:.1f}" for i in range(total_columns)]
            )

        for row in np.array(param_diag_matrix).T:
            param_diag.add_row(*row)

        # Add the global diagnostics to the main table
        main_table.add_row(global_diag)
        main_table.add_row(param_diag)

        covariance_matrix = Table(title="Covariance Matrix", expand=False)
        covariance_matrix.add_column("Parameter", justify="right", style="bold")
        for i in range(total_columns):
            covariance_matrix.add_column(
                mcat.col_name(i).split(":")[-1], justify="right"
            )

        for i in range(total_columns):
            row = [mcat.col_name(i).split(":")[-1]]
            for j in range(total_columns):
                cov_ij = full_stats.get_cor(i, j)
                cor_ij_string = f"{cov_ij * 100.0: 3.0f}%"
                styles_array = [
                    "bold bright_red",
                    "bright_red",
                    "dim bright_red",
                    "dim bright_green",
                    "bright_green",
                    "bold bright_green",
                ]
                cov_color_index = int(np.round((cov_ij + 1.0) * 2.5))

                row.append(Text(cor_ij_string, style=styles_array[cov_color_index]))
            covariance_matrix.add_row(*row)

        main_table.add_row(covariance_matrix)
        self.console.print(main_table)

        self.end_experiment()

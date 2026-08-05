#
# loading.py
#
# Wed Feb 14 18:59:34 2024
# Copyright  2024  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# loading.py
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


"""NumCosmo APP dataclasses and subcommands to load data.

This module contains dataclasses and subcommands to load data from files.
"""

import dataclasses
from typing import Optional, Annotated, cast

from pathlib import Path
import typer

from numcosmo_py import Ncm
from .logging import AppLogging


@dataclasses.dataclass(kw_only=True)
class LoadExperiment(AppLogging):
    """Load an experiment file.

    Common block for commands that load an experiment. All commands that load an
    experiment should inherit from this class.
    """

    experiment: Annotated[
        Path, typer.Argument(help="Path to the experiment file to fit.")
    ]
    product_file: Annotated[
        bool,
        typer.Option(
            "--product-file",
            "-p",
            help=(
                "If given, the product file is written, the file name is the same as "
                "the experiment file with the extension .product.yaml. "
                "This option is incompatible with the output and starting-point "
                "options since the product file contains the output and starting "
                "point."
            ),
        ),
    ] = False
    starting_point: Annotated[
        Optional[Path],
        typer.Option(
            "--starting-point",
            "-s",
            help=(
                "Path to the file containing the starting point for the fit. "
                "The output of a previous fit can be used."
            ),
        ),
    ] = None
    output: Annotated[
        Optional[Path],
        typer.Option(
            "--output",
            "-o",
            help="Path to the output file, if given, the computed results are written "
            "to this file, otherwise they are not saved.",
        ),
    ] = None

    def __post_init__(self) -> None:
        """Load the experiment file and prepare the experiment."""
        ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)

        builders_file = self.experiment.with_suffix(".builders.yaml")
        # Builders file is optional.
        #
        # The models need to be created before initializing NumCosmo, this is necessary
        # since in a MPI environment the slave processes need to have the models
        # registered before waiting for commands from the master process.
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

        # Initialize logging from parent class (also initializes NumCosmo)
        super().__post_init__()

        dataset_file = self.experiment.with_suffix(".dataset.gvar")
        if dataset_file.exists():
            dataset = ser.from_binfile(
                self.experiment.with_suffix(".dataset.gvar").absolute().as_posix()
            )
            assert isinstance(dataset, Ncm.Dataset)

        experiment_objects = ser.dict_str_from_yaml_file(
            self.experiment.absolute().as_posix()
        )

        functions_file = self.experiment.with_suffix(".functions.yaml")
        self.functions: Optional[Ncm.ObjArray] = None
        if functions_file.exists():
            functions: Ncm.ObjArray = ser.array_from_yaml_file(
                functions_file.absolute().as_posix()
            )
            assert isinstance(functions, Ncm.ObjArray)
            self.functions = functions
            for i in range(functions.len()):
                function: Ncm.MSetFunc = cast(Ncm.MSetFunc, functions.get(i))
                if not isinstance(function, Ncm.MSetFunc):
                    raise RuntimeError(f"Invalid function file {functions_file}.")

        if self.product_file:
            if self.output is not None:
                raise RuntimeError(
                    f"The product file option is incompatible with the output "
                    f"option {self.output}."
                )
            if self.starting_point is not None:
                raise RuntimeError(
                    "The product file option is incompatible with the starting-point "
                    "option."
                )
            self.output = self.experiment.with_suffix(".product.yaml")

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

        if self.output is not None:
            if self.output.exists():
                ser.reset(False)
                self.output_dict = ser.dict_str_from_yaml_file(
                    self.output.absolute().as_posix()
                )
            else:
                self.output_dict = Ncm.ObjDictStr.new()

        saved_mset = self._load_saved_mset()
        if saved_mset is not None:
            if not mset.cmp(saved_mset, True):
                raise RuntimeError(
                    f"Starting point file {self.starting_point} "
                    f"does not match experiment."
                )
            mset.param_set_mset(saved_mset)

        self.likelihood = likelihood
        self.mset = mset

    def _load_saved_mset(self) -> Optional[Ncm.MSet]:
        """Load the saved model.

        Load the saved model-set from the starting point file or the product file.
        """
        if self.starting_point is not None:
            if not self.starting_point.exists():
                raise RuntimeError(
                    f"Starting point file {self.starting_point} not found."
                )

            ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
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

            return saved_mset

        if self.product_file:
            product_mset: Ncm.MSet = cast(Ncm.MSet, self.output_dict.get("model-set"))
            if product_mset is not None:
                assert isinstance(product_mset, Ncm.MSet)
                return product_mset

        return None

    def end_experiment(self):
        """End the experiment and writes the output file."""
        if self.output is not None:
            ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
            ser.dict_str_to_yaml_file(
                self.output_dict, self.output.absolute().as_posix()
            )
        self.close_logging()


def _catalog_indices(
    mcat: Ncm.MSetCatalog,
    total_columns: int,
    include: Optional[list[str]],
    exclude: Optional[list[str]],
) -> list[int]:
    """Resolve the --include/--exclude column selection to a list of indices."""
    include = include or []
    exclude = exclude or []

    if not include and not exclude:
        return list(range(total_columns))

    indices = []
    if include and exclude:
        for i in range(total_columns):
            name = mcat.col_full_name(i)
            if any(s in name for s in include) and not any(s in name for s in exclude):
                indices.append(i)
    elif include:
        for i in range(total_columns):
            name = mcat.col_full_name(i)
            if any(s in name for s in include):
                indices.append(i)
    else:
        for i in range(total_columns):
            name = mcat.col_full_name(i)
            if not any(s in name for s in exclude):
                indices.append(i)

    return indices


@dataclasses.dataclass
class LoadedCatalog:
    """The contents of an MCMC catalog file, loaded and ready for analysis.

    Catalog files are self-sufficient: their model-set and, if present, the
    functions that computed their extra columns are embedded in the file
    itself (FITS HDU0), so no companion experiment file is needed to load
    one.
    """

    mcat: Ncm.MSetCatalog
    mset: Ncm.MSet
    functions: Optional[Ncm.ObjArray]
    fparams_len: int
    nadd_vals: int
    total_columns: int
    nchains: int
    indices: list[int]
    full_stats: Ncm.StatsVec
    stats: Ncm.StatsVec
    nitems: int


def _resolve_burnin_rows(mcmc_file: Path, burnin: int, tail: Optional[int]) -> int:
    """Resolve a --burnin/--tail request (in iterations) to a row count.

    `burnin` discards the first N iterations (ensemble steps); `tail` keeps
    only the last N instead. Peeks the catalog's row/chain counts first
    (cheap: a few FITS header keys, no model-set deserialization) to convert
    iterations to rows and validate the request before the catalog is
    actually opened.
    """
    if tail is not None and burnin != 0:
        raise typer.BadParameter("Give at most one of --burnin and --tail.")

    if burnin == 0 and tail is None:
        return 0

    nrows, nchains, _first_id = Ncm.MSetCatalog.peek_info_from_file(
        mcmc_file.absolute().as_posix()
    )
    n_iterations = nrows // nchains

    if tail is not None:
        if tail < 0:
            raise typer.BadParameter(f"--tail must be non-negative, got {tail}.")
        burnin_iterations = max(0, n_iterations - tail)
    else:
        burnin_iterations = burnin

    if burnin_iterations > n_iterations:
        raise typer.BadParameter(
            f"--burnin of {burnin_iterations} iteration(s) exceeds catalog "
            f"{mcmc_file}: it only has {n_iterations} iteration(s) "
            f"({nrows} rows, {nchains} chains)."
        )

    return burnin_iterations * nchains


def load_catalog(
    mcmc_file: Path,
    burnin: int = 0,
    tail: Optional[int] = None,
    include: Optional[list[str]] = None,
    exclude: Optional[list[str]] = None,
) -> LoadedCatalog:
    """Load an MCMC catalog file and prepare it for analysis.

    `burnin` and `tail` are in iterations (ensemble steps), not raw rows --
    see _resolve_burnin_rows().
    """
    if not mcmc_file.exists():
        raise typer.BadParameter(f"MCMC file {mcmc_file} not found.")

    burnin_rows = _resolve_burnin_rows(mcmc_file, burnin, tail)

    mcat: Ncm.MSetCatalog = Ncm.MSetCatalog.new_from_file_ro(
        mcmc_file.absolute().as_posix(), burnin_rows
    )
    assert isinstance(mcat, Ncm.MSetCatalog)

    mset: Ncm.MSet = mcat.peek_mset()
    assert isinstance(mset, Ncm.MSet)
    mset.prepare_fparam_map()

    functions: Optional[Ncm.ObjArray] = mcat.peek_functions_array()

    fparams_len = mset.fparams_len()
    nadd_vals: int = mcat.nadd_vals()
    total_columns: int = fparams_len + nadd_vals
    nchains: int = mcat.nchains()

    indices = _catalog_indices(mcat, total_columns, include, exclude)

    full_stats: Ncm.StatsVec = mcat.peek_pstats()
    assert isinstance(full_stats, Ncm.StatsVec)

    if nchains > 1:
        stats: Ncm.StatsVec = mcat.peek_e_mean_stats()
    else:
        stats = mcat.peek_pstats()
    assert isinstance(stats, Ncm.StatsVec)

    nitems: int = stats.nitens()

    return LoadedCatalog(
        mcat=mcat,
        mset=mset,
        functions=functions,
        fparams_len=fparams_len,
        nadd_vals=nadd_vals,
        total_columns=total_columns,
        nchains=nchains,
        indices=indices,
        full_stats=full_stats,
        stats=stats,
        nitems=nitems,
    )


@dataclasses.dataclass(kw_only=True)
class LoadCatalog(AppLogging):
    """Loads a single MCMC catalog file for analysis.

    No experiment file is needed: the catalog is self-sufficient (see
    load_catalog()). Commands that need to combine several catalogs (e.g.
    plot-corner) use load_catalog() directly instead of this class.
    """

    mcmc_file: Annotated[
        Path,
        typer.Argument(
            help="Path to the MCMC catalog file.",
        ),
    ]

    burnin: Annotated[
        int,
        typer.Option(
            help="Number of iterations (ensemble steps) to discard as burnin.",
            min=0,
        ),
    ] = 0

    tail: Annotated[
        Optional[int],
        typer.Option(
            help=(
                "Keep only the last N iterations (ensemble steps) instead of "
                "burning in from the start. Give at most one of --burnin/--tail."
            ),
            min=0,
        ),
    ] = None

    include: Annotated[
        Optional[list[str]],
        typer.Option(
            help="List of parameters and or model names to include in the analysis.",
        ),
    ] = None

    exclude: Annotated[
        Optional[list[str]],
        typer.Option(
            help="List of parameters and or model names to exclude from the analysis.",
        ),
    ] = None

    output: Annotated[
        Optional[Path],
        typer.Option(
            "--output",
            "-o",
            help="Path to the output file, if given.",
        ),
    ] = None

    # These are set in __post_init__ from load_catalog(), not from the CLI.
    mcat: Ncm.MSetCatalog = dataclasses.field(init=False)
    mset: Ncm.MSet = dataclasses.field(init=False)
    functions: Optional[Ncm.ObjArray] = dataclasses.field(init=False)
    fparams_len: int = dataclasses.field(init=False)
    nadd_vals: int = dataclasses.field(init=False)
    total_columns: int = dataclasses.field(init=False)
    nchains: int = dataclasses.field(init=False)
    indices: list[int] = dataclasses.field(init=False)
    full_stats: Ncm.StatsVec = dataclasses.field(init=False)
    stats: Ncm.StatsVec = dataclasses.field(init=False)
    nitems: int = dataclasses.field(init=False)

    def __post_init__(self) -> None:
        """Load the MCMC file and prepare the catalog."""
        super().__post_init__()

        loaded = load_catalog(
            self.mcmc_file, self.burnin, self.tail, self.include, self.exclude
        )
        self.mcat = loaded.mcat
        self.mset = loaded.mset
        self.functions = loaded.functions
        self.fparams_len = loaded.fparams_len
        self.nadd_vals = loaded.nadd_vals
        self.total_columns = loaded.total_columns
        self.nchains = loaded.nchains
        self.indices = loaded.indices
        self.full_stats = loaded.full_stats
        self.stats = loaded.stats
        self.nitems = loaded.nitems

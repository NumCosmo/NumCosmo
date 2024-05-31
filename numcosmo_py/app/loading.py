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
from numcosmo_py.sampling import set_ncm_console


@dataclasses.dataclass(kw_only=True)
class LoadExperiment:
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

    log_file: Annotated[
        Optional[Path],
        typer.Option(
            "--log-file",
            "-l",
            help="Path to the file where the log should be written.",
        ),
    ] = None

    def __post_init__(self) -> None:
        """Load the experiment file and prepare the experiment."""
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
        self.console_io = None
        if self.log_file:
            self.console_io = open(self.log_file, "w", encoding="utf-8")

        console = set_ncm_console(self.console_io)

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
                    "The product file option is incompatible with the output option."
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

        self.console = console
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
        if self.console_io is not None:
            self.console_io.close()


@dataclasses.dataclass(kw_only=True)
class LoadCatalog(LoadExperiment):
    """Analyzes the results of a MCMC run."""

    mcmc_file: Annotated[
        Path,
        typer.Argument(
            help="Path to the MCMC file.",
        ),
    ]

    burnin: Annotated[
        int,
        typer.Option(
            help="Number of samples to discard as burnin.",
            min=0,
        ),
    ] = 0

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

    def __post_init__(self) -> None:
        """Load the MCMC file and prepare the catalog."""
        super().__post_init__()

        if not self.mcmc_file.exists():
            raise RuntimeError(f"MCMC file {self.mcmc_file} not found.")

        self.mcat: Ncm.MSetCatalog = Ncm.MSetCatalog.new_from_file_ro(
            self.mcmc_file.absolute().as_posix(), self.burnin
        )
        assert isinstance(self.mcat, Ncm.MSetCatalog)

        self.catalog_mset: Ncm.MSet = self.mcat.peek_mset()
        assert isinstance(self.catalog_mset, Ncm.MSet)

        self.catalog_mset.prepare_fparam_map()
        self.fparams_len = self.catalog_mset.fparams_len()
        self.nadd_vals: int = self.mcat.nadd_vals()
        self.total_columns: int = self.fparams_len + self.nadd_vals
        self.nchains: int = self.mcat.nchains()

        self._extract_indices()

        self.full_stats: Ncm.StatsVec = self.mcat.peek_pstats()
        assert isinstance(self.full_stats, Ncm.StatsVec)

        if self.nchains > 1:
            self.stats: Ncm.StatsVec = self.mcat.peek_e_mean_stats()
            assert isinstance(self.stats, Ncm.StatsVec)
        else:
            self.stats = self.mcat.peek_pstats()
            assert isinstance(self.stats, Ncm.StatsVec)

        self.nitems: int = self.stats.nitens()

    def _extract_indices(self):
        """Extract the indices to include in the analysis."""
        if self.include is None:
            self.include = []
        if self.exclude is None:
            self.exclude = []
        assert self.include is not None
        assert self.exclude is not None
        if not self.include and not self.exclude:
            self.indices = list(range(self.total_columns))
        else:
            self.indices = []
            if self.include and self.exclude:
                for i in range(self.total_columns):
                    name = self.mcat.col_full_name(i)
                    if any(s in name for s in self.include) and not any(
                        s in name for s in self.exclude
                    ):
                        self.indices.append(i)
            elif self.include:
                for i in range(self.total_columns):
                    name = self.mcat.col_full_name(i)
                    if any(s in name for s in self.include):
                        self.indices.append(i)
            else:
                for i in range(self.total_columns):
                    name = self.mcat.col_full_name(i)
                    if not any(s in name for s in self.exclude):
                        self.indices.append(i)

#
# run_fit.py
#
# Wed Feb 14 18:59:34 2024
# Copyright  2024  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# run_fit.py
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

"""NumCosmo APP subcommand to fit data."""

import dataclasses
from typing import Optional, Annotated, Tuple

import typer

from .. import Ncm
from .loading import LoadExperiment
from ..sampling import (
    check_runner_algorithm,
    FitGradType,
    FitRunMessages,
    FitRunner,
    NcmFitLogger,
)


@dataclasses.dataclass
class RunCommonOptions(LoadExperiment):
    """Common options for the run command."""

    runner: Annotated[
        FitRunner,
        typer.Option(
            help="Algorithm to use for the fit.",
        ),
    ] = FitRunner.NLOPT
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
    ] = FitGradType.NUMDIFF_FORWARD
    run_messages: Annotated[
        FitRunMessages,
        typer.Option(
            help="Verbosity level for the fit.",
        ),
    ] = FitRunMessages.SIMPLE

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


@dataclasses.dataclass
class RunTest(RunCommonOptions):
    """Loads the experiment file and computes the likelihood once."""

    def __post_init__(self) -> None:
        super().__post_init__()
        self.mset.param_set_all_ftype(Ncm.ParamType.FIXED)
        self.mset.prepare_fparam_map()
        self.fit.log_info()
        self.fit.run(FitRunMessages.SIMPLE.genum)

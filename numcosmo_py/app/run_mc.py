#
# run_mc.py
#
# Tue Nov 19 09:36:55 2024
# Copyright  2024  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# run_mc.py
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

"""NumCosmo APP subcommand to run Monte Carlo Analysis."""

import dataclasses
from typing import Annotated

import typer

from .. import Ncm
from ..sampling import (
    FitMCResampleType,
    FitRunMessages,
)
from .run_fit import RunCommonOptions


@dataclasses.dataclass(kw_only=True)
class RunMC(RunCommonOptions):
    """Computes the Monte Carlo Analysis."""

    run_type: Annotated[
        FitMCResampleType,
        typer.Option(
            help="Resampling type for the fit.",
        ),
    ] = FitMCResampleType.FROM_MODEL

    run_messages: Annotated[
        FitRunMessages,
        typer.Option(
            help="Verbosity level for the fit.",
        ),
    ] = FitRunMessages.SIMPLE

    nthreads: Annotated[
        int,
        typer.Option(
            help="Number of threads to use for the fit.",
        ),
    ] = 1

    nmc: Annotated[
        int,
        typer.Option(
            help="Number of Monte Carlo samples to generate.",
        ),
    ] = 100

    seed: Annotated[
        int | None,
        typer.Option(
            help="Seed for the random number generator. "
            "If None (default), a random seed is used.",
        ),
    ] = None

    def __post_init__(self) -> None:
        """Compute Monte Carlo Analysis."""
        super().__post_init__()

        mc = Ncm.FitMC.new(
            fit=self.fit, rtype=self.run_type.genum, mtype=self.run_messages.genum
        )

        if self.output is None:
            raise typer.BadParameter(
                "Either --output or --product-file must be provided."
            )

        mc.set_nthreads(self.nthreads)
        mc.set_data_file(self.output.with_suffix(".mc.fits").absolute().as_posix())

        if self.seed is not None:
            rng = Ncm.RNG.seeded_new("mt19937", self.seed)
            mc.set_rng(rng)

        mc.start_run()
        mc.run(self.nmc)
        mc.end_run()

        self.end_experiment()

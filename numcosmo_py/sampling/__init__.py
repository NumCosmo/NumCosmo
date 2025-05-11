#
# __init__.py
#
# Wed Feb 8 10:00:00 2023
# Copyright  2023  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# __init__.py
# Copyright (C) 2023 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

"""Sampling module for numcosmo."""

from typing import Optional, Union, Type, IO
from enum import StrEnum, auto

from rich.console import Console
from rich.highlighter import RegexHighlighter
from rich.progress import (
    BarColumn,
    MofNCompleteColumn,
    Progress,
    ProgressColumn,
    SpinnerColumn,
    Task,
    TaskID,
    TextColumn,
    TimeElapsedColumn,
    TimeRemainingColumn,
)
from rich.text import Text
from rich.theme import Theme

from .. import Ncm, GEnum


class FitRunner(GEnum):
    """Fit algorithm for Ncm.Fit."""

    GSL_LS = Ncm.FitType.GSL_LS
    GSL_MM = Ncm.FitType.GSL_MM
    GSL_MMS = Ncm.FitType.GSL_MMS
    LEVMAR = Ncm.FitType.LEVMAR
    NLOPT = Ncm.FitType.NLOPT


class FitRunMessages(GEnum):
    """Fit messages for Ncm.Fit."""

    NONE = Ncm.FitRunMsgs.NONE
    SIMPLE = Ncm.FitRunMsgs.SIMPLE
    FULL = Ncm.FitRunMsgs.FULL


class FitGradType(GEnum):
    """Fit gradient type for Ncm.Fit."""

    NUMDIFF_FORWARD = Ncm.FitGradType.NUMDIFF_FORWARD
    NUMDIFF_CENTRAL = Ncm.FitGradType.NUMDIFF_CENTRAL
    NUMDIFF_ACCURATE = Ncm.FitGradType.NUMDIFF_ACCURATE


class FitMCResampleType(GEnum):
    """Fit MC resample type for Ncm.Fit."""

    FROM_MODEL = Ncm.FitMCResampleType.FROM_MODEL
    BOOTSTRAP_NOMIX = Ncm.FitMCResampleType.BOOTSTRAP_NOMIX
    BOOTSTRAP_MIX = Ncm.FitMCResampleType.BOOTSTRAP_MIX


def get_algorithms(
    runner: FitRunner,
) -> Optional[
    Union[
        Type[Ncm.FitNloptAlgorithm],
        Type[Ncm.FitLevmarAlgos],
        Type[Ncm.FitGSLMMSAlgos],
        Type[Ncm.FitGSLMMAlgos],
    ]
]:
    """Get algorithms for a given runner."""
    if runner == FitRunner.NLOPT:
        return Ncm.FitNloptAlgorithm
    if runner == FitRunner.LEVMAR:
        return Ncm.FitLevmarAlgos
    if runner == FitRunner.GSL_MMS:
        return Ncm.FitGSLMMSAlgos
    if runner == FitRunner.GSL_MM:
        return Ncm.FitGSLMMAlgos
    if runner == FitRunner.GSL_LS:
        return None
    raise RuntimeError(f"Runner {runner} not found.")


def check_runner_algorithm(runner: FitRunner, algorithm: Optional[str]):
    """Check if algorithm is valid."""
    if algorithm is not None:
        algorithms = get_algorithms(runner)
        if algorithms is None:
            raise RuntimeError(f"Runner {runner} do not support algorithms.")
        if Ncm.cfg_get_enum_by_id_name_nick(algorithms, algorithm) is None:
            Ncm.cfg_enum_print_all(algorithms, "Allowed algorithms")
            raise RuntimeError(f"Algorithm {algorithm} not found for runner {runner}.")


class NcmHighlighter(RegexHighlighter):
    """Apply style to anything that looks like an email."""

    base_style = "Ncm."
    highlights = [
        r"(?P<FitTypeFIXED>FIXED)",
        r"(?P<FitTypeFREE>FREE)",
        r"\b(?P<float>\d+(\.?\d+)?([eE][-+]?\d+)?)\b",
        r"(?P<float_signed>[-+]\d+(\.?\d+)?([eE][-+]?\d+)?)\b",
        r"(?P<datetime>\d{2}:\d{2}:\d{2}(\.\d{3,})?)",
        r"\b(?P<percentage>\d+(\.?\d+)?%)",
    ]


def set_ncm_console(file: Optional[IO[str]]) -> Console:
    """Set console for Ncm.Fit."""
    theme = Theme(
        {
            "Ncm.FitTypeFIXED": "bold red",
            "Ncm.FitTypeFREE": "bold green",
            "Ncm.datetime": "bold yellow",
            "Ncm.float": "bold cyan",
            "Ncm.float_signed": "bold cyan",
            "Ncm.percentage": "bold magenta",
        }
    )
    console = Console(
        highlighter=NcmHighlighter(), theme=theme, soft_wrap=True, file=file
    )

    Ncm.cfg_set_log_handler(lambda msg: console.print(msg, end=""))

    return console


class FitSpeedColumn(ProgressColumn):
    """Renders human readable transfer speed."""

    def render(self, task: Task) -> Text:
        """Show data transfer speed."""
        speed = task.finished_speed or task.speed
        if speed is None:
            return Text("?", style="progress.data.speed")

        return Text(f"{speed:.2f}it/s", style="progress.data.speed")


class NcmFitLogger:
    """Class implementing logging functions for Ncm.Fit."""

    def __init__(self, console: Optional[Console]) -> None:
        """Initialize NcmFitLogger."""
        self.console = console if console is not None else Console()
        self.progress = Progress(
            TextColumn("# [progress.description]{task.description}: "),
            SpinnerColumn(),
            MofNCompleteColumn(),
            BarColumn(bar_width=None),
            TimeElapsedColumn(),
            FitSpeedColumn(),
            TimeRemainingColumn(),
            transient=False,
            console=self.console,
            expand=True,
        )
        self.task: Optional[TaskID] = None

    def write_progress(self, _fit: Ncm.Fit, message: str):
        """Write progress to Rich."""
        self.console.print(message, end="")

    def update_progress(self, _fit: Ncm.Fit, n: Union[int, float]):
        """Update progress bar."""
        total = self.progress.tasks[0].total
        assert self.task is not None
        assert total is not None
        if n > total:
            new_total = 2 * total
            self.progress.update(self.task, completed=n, total=new_total)
        else:
            self.progress.update(self.task, completed=n)

    def start_update(self, _fit: Ncm.Fit, start_message: str):
        """Start updates."""
        self.task = self.progress.add_task(start_message)
        self.progress.start()
        self.progress.start_task(self.task)
        self.progress.update(self.task, completed=0, total=10)

    def end_update(self, _fit: Ncm.Fit, _end_message: str):
        """End updates."""
        assert self.task is not None

        self.progress.stop_task(self.task)
        self.progress.stop()

        self.progress.remove_task(self.task)
        self.task = None


class FisherType(StrEnum):
    """Possible linear matter power spectrum models."""

    OBSERVED = auto()
    EXPECTED = auto()

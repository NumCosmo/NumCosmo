#
# logging.py
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

"""NumCosmo APP base dataclass for logging configuration.

This module provides a base dataclass that handles console and logging
configuration for CLI commands. Commands that need logging should inherit
from this class.
"""

import dataclasses
from pathlib import Path
from typing import Annotated, Optional, IO

import typer
from rich.console import Console

from numcosmo_py import Ncm
from numcosmo_py.sampling import set_ncm_console


@dataclasses.dataclass(kw_only=True)
class AppLogging:
    """Base dataclass for CLI commands with logging support.

    Provides centralized console and log file handling. All commands that
    produce output should inherit from this class to ensure consistent
    logging behavior.

    Attributes:
        log_file: Optional path to write log output
        console: Rich Console instance for output (created in __post_init__)
        console_io: File handle for log file (if log_file is specified)
    """

    log_file: Annotated[
        Optional[Path],
        typer.Option(
            "--log-file",
            "-l",
            help="Path to the file where the log should be written.",
        ),
    ] = None

    quite: Annotated[
        bool,
        typer.Option(
            "--quite",
            "-q",
            help="Suppress console output.",
            is_flag=True,
        ),
    ] = False

    # These are set in __post_init__, not from CLI
    console: Console = dataclasses.field(init=False)
    console_io: Optional[IO[str]] = dataclasses.field(init=False, default=None)

    def __post_init__(self) -> None:
        """Initialize logging configuration.

        Sets up the console and optional log file. Also initializes NumCosmo
        configuration and sets up the NCM log handler.
        """
        Ncm.cfg_init()

        self.console_io = None
        if self.log_file:
            self.console_io = open(self.log_file, "w", encoding="utf-8")

        self.console = set_ncm_console(self.console_io, self.quite)

    def close_logging(self) -> None:
        """Close logging resources.

        Should be called when the command finishes to properly close
        any open log files.
        """
        if self.console_io is not None:
            self.console_io.close()
            self.console_io = None

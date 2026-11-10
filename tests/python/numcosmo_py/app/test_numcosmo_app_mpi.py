#!/usr/bin/env python
#
# test_numcosmo_app_mpi.py
#
# Mon Nov 2 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_numcosmo_app_mpi.py
# Copyright (C) 2026 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

"""Test the numcosmo command line under mpiexec.

The command line runs in a subprocess launched by mpiexec, whose path meson
test sets in NUMCOSMO_TEST_MPIEXEC when the library is built with MPI.
"""

import os
import subprocess
import sys
from pathlib import Path

import pytest

from numcosmo_py import Ncm

pytestmark = [
    pytest.mark.app,
    pytest.mark.skipif(
        "NUMCOSMO_TEST_MPIEXEC" not in os.environ,
        reason="needs NUMCOSMO_TEST_MPIEXEC, set by meson test for an MPI build",
    ),
]

Ncm.cfg_init()


def test_run_mcmc_apes_mpi_builder_model(builder_experiment: Path) -> None:
    """Run APES with --parallel mpi on an experiment with a builders file.

    The worker rank deserializes the model set sent by rank 0, which holds the
    model of the builders file, so the worker must create that model before it
    enters the worker loop.
    """
    output = builder_experiment.with_suffix(".out.yaml")
    command = [
        os.environ["NUMCOSMO_TEST_MPIEXEC"],
        "-n",
        "2",
        sys.executable,
        "-c",
        "from numcosmo_py.app import app; app()",
        "run",
        "mcmc",
        "apes",
        builder_experiment.as_posix(),
        "--parallel",
        "mpi",
        "--nwalkers",
        "60",
        "--nsamples",
        "120",
        "--output",
        output.as_posix(),
    ]
    # python -c puts the working directory first on sys.path; tmp_path holds no
    # numcosmo_py, so the subprocess imports the one the test imported.
    result = subprocess.run(
        command,
        cwd=builder_experiment.parent,
        capture_output=True,
        text=True,
        timeout=600,
        check=False,
    )

    assert result.returncode == 0, result.stdout + result.stderr
    assert output.with_suffix(".mcmc.fits").exists()

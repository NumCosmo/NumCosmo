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

from numcosmo_py import GObject, Ncm

pytestmark = [
    pytest.mark.app,
    pytest.mark.skipif(
        "NUMCOSMO_TEST_MPIEXEC" not in os.environ,
        reason="needs NUMCOSMO_TEST_MPIEXEC, set by meson test for an MPI build",
    ),
]

Ncm.cfg_init()


@pytest.fixture(name="builder_experiment")
def fixture_builder_experiment(tmp_path: Path) -> Path:
    """Write an experiment whose model set holds a model of its builders file."""
    model_builder = Ncm.ModelBuilder.new(
        Ncm.Model, "NcmTestAppMPIBuilderModel", "Model created by a builder"
    )
    model_builder.add_sparam("p_1", "p1", 0.0, 1.0, 0.1, 0.0, 0.5, Ncm.ParamType.FIXED)
    BuilderModel = model_builder.create()
    GObject.new(BuilderModel)

    mset = Ncm.MSet.new_array([Ncm.ModelMVND.new(2), BuilderModel.pytype()])
    for i in range(2):
        mset.param_set_ftype(Ncm.ModelMVND.id(), i, Ncm.ParamType.FREE)

    rng = Ncm.RNG.seeded_new(None, 1234)
    data_mvnd = Ncm.DataGaussCovMVND.new_full(2, 1.0, 2.0, 30.0, 0.0, 0.0, rng)
    likelihood = Ncm.Likelihood.new(dset=Ncm.Dataset.new_array([data_mvnd]))

    experiment_file = tmp_path / "builder_experiment.yaml"
    ser = Ncm.Serialize.new(Ncm.SerializeOpt.KEEP_SHARED)

    model_builders = Ncm.ObjDictStr.new()
    model_builders.add("NcmTestAppMPIBuilderModel", model_builder)
    ser.dict_str_to_yaml_file(
        model_builders, experiment_file.with_suffix(".builders.yaml").as_posix()
    )

    experiment = Ncm.ObjDictStr.new()
    experiment.add("model-set", mset)
    experiment.add("likelihood", likelihood)
    ser.dict_str_to_yaml_file(experiment, experiment_file.as_posix())

    return experiment_file


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

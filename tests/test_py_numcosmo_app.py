#!/usr/bin/env python
#
# test_py_numcosmo_app.py
#
# Wed Jan 31 10:04:00 2024
# Copyright  2023  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_py_numcosmo_app.py
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

"""Test NumCosmo app."""

from typing import Tuple
from pathlib import Path
import pytest
from typer.testing import CliRunner

from numcosmo_py import Ncm
from numcosmo_py.app import app
from numcosmo_py.sampling import FitRunner, FitGradType, FitRunMessages, FisherType

runner = CliRunner()

Ncm.cfg_init()


@pytest.fixture(name="simple_experiment")
def fixture_simple_experiment(tmp_path) -> Tuple[Path, Ncm.ObjDictStr]:
    """Create a simple experiment."""

    rng = Ncm.RNG.seeded_new(None, 1234)
    model_mvnd = Ncm.ModelMVND.new(5)
    mset = Ncm.MSet.new_array([model_mvnd])
    mset.param_set_all_ftype(Ncm.ParamType.FREE)

    data_mvnd = Ncm.DataGaussCovMVND.new_full(5, 1.0, 2.0, 30.0, 0.0, 0.0, rng)
    likelihood = Ncm.Likelihood.new(dset=Ncm.Dataset.new_array([data_mvnd]))

    tmp_file = tmp_path / "simple_experiment.yaml"

    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)

    experiment = Ncm.ObjDictStr.new()
    experiment.add("model-set", mset)
    experiment.add("likelihood", likelihood)

    ser.dict_str_to_yaml_file(experiment, tmp_file.as_posix())

    return tmp_file, experiment


@pytest.fixture(name="fit_runner", params=[e.value for e in FitRunner])
def fixture_fit_runner(request) -> str:
    """Returns a fit runner"""

    return request.param


@pytest.fixture(name="fit_grad_type", params=[e.value for e in FitGradType])
def fixture_fit_grad_type(request) -> str:
    """Returns a fit grad type"""

    return request.param


@pytest.fixture(name="fit_run_messages", params=[e.value for e in FitRunMessages])
def fixture_fit_run_messages(request) -> str:
    """Returns a fit run messages"""

    return request.param


@pytest.fixture(name="fisher_type", params=[e.value for e in FisherType])
def fixture_fisher_type(request) -> str:
    """Returns fisher type"""

    return request.param


def test_run_test(simple_experiment):
    """Test run test."""

    filename, _ = simple_experiment
    result = runner.invoke(app, ["run", "test", filename.as_posix()])
    assert result.exit_code == 0


def test_run_fit(simple_experiment):
    """Test run fit."""
    filename, _ = simple_experiment
    result = runner.invoke(app, ["run", "fit", filename.as_posix()])
    assert result.exit_code == 0


def test_run_fit_runner(simple_experiment, fit_runner):
    """Test run fit."""
    filename, _ = simple_experiment
    result = runner.invoke(
        app, ["run", "fit", filename.as_posix(), "--runner", fit_runner]
    )
    assert result.exit_code == 0


def test_run_fit_grad_type(simple_experiment, fit_grad_type):
    """Test run fit."""
    filename, _ = simple_experiment
    result = runner.invoke(
        app, ["run", "fit", filename.as_posix(), "--grad-type", fit_grad_type]
    )
    assert result.exit_code == 0


def test_run_fit_run_messages(simple_experiment, fit_run_messages):
    """Test run fit."""
    filename, _ = simple_experiment
    result = runner.invoke(
        app, ["run", "fit", filename.as_posix(), "--run-messages", fit_run_messages]
    )
    assert result.exit_code == 0


def test_run_fit_runner_grad_type(simple_experiment, fit_runner, fit_grad_type):
    """Test run fit."""
    filename, _ = simple_experiment
    result = runner.invoke(
        app,
        [
            "run",
            "fit",
            filename.as_posix(),
            "--runner",
            fit_runner,
            "--grad-type",
            fit_grad_type,
        ],
    )
    assert result.exit_code == 0


def test_run_fit_runner_grad_type_run_messages(
    simple_experiment, fit_runner, fit_grad_type, fit_run_messages
):
    """Test run fit."""
    filename, _ = simple_experiment
    result = runner.invoke(
        app,
        [
            "run",
            "fit",
            filename.as_posix(),
            "--runner",
            fit_runner,
            "--grad-type",
            fit_grad_type,
            "--run-messages",
            fit_run_messages,
        ],
    )
    assert result.exit_code == 0


def test_run_fit_starting_point(
    simple_experiment, fit_runner, fit_grad_type, fit_run_messages
):
    """Test run fit with starting point."""

    filename, _ = simple_experiment
    output = filename.with_suffix(".out.yaml")
    result = runner.invoke(
        app,
        ["run", "fit", filename.as_posix(), "--output", output.as_posix()],
    )
    assert result.exit_code == 0

    result = runner.invoke(
        app,
        [
            "run",
            "fit",
            filename.as_posix(),
            "--starting-point",
            output.as_posix(),
            "--runner",
            fit_runner,
            "--grad-type",
            fit_grad_type,
            "--run-messages",
            fit_run_messages,
        ],
    )

    assert result.exit_code == 0


def test_run_fit_restart(
    simple_experiment, fit_runner, fit_grad_type, fit_run_messages
):
    """Test run fit with starting point."""

    filename, _ = simple_experiment
    result = runner.invoke(
        app,
        [
            "run",
            "fit",
            filename.as_posix(),
            "--runner",
            fit_runner,
            "--grad-type",
            fit_grad_type,
            "--run-messages",
            fit_run_messages,
            "--restart",
            "1.0e-3",
            "0.0",
        ],
    )
    assert result.exit_code == 0


def test_run_theory_vector(simple_experiment):
    """Test run theory vector."""
    filename, _ = simple_experiment
    output = filename.with_suffix(".out.yaml")
    result = runner.invoke(
        app,
        ["run", "theory-vector", filename.as_posix(), "--output", output.as_posix()],
    )

    assert result.exit_code == 0

    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
    output_dict = ser.dict_str_from_yaml_file(output.as_posix())
    assert isinstance(output_dict, Ncm.ObjDictStr)
    v1 = output_dict.get("theory-vector")
    assert isinstance(v1, Ncm.Vector)


def test_run_theory_vector_starting_point(simple_experiment):
    """Test run theory vector with starting point."""

    filename, _ = simple_experiment
    output = filename.with_suffix(".out.yaml")
    result = runner.invoke(
        app,
        ["run", "fit", filename.as_posix(), "--output", output.as_posix()],
    )

    assert result.exit_code == 0

    result = runner.invoke(
        app,
        [
            "run",
            "theory-vector",
            filename.as_posix(),
            "--starting-point",
            output.as_posix(),
            "--output",
            output.as_posix(),
        ],
    )

    assert result.exit_code == 0

    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
    output_dict = ser.dict_str_from_yaml_file(output.as_posix())
    assert isinstance(output_dict, Ncm.ObjDictStr)
    v1 = output_dict.get("theory-vector")
    assert isinstance(v1, Ncm.Vector)


def test_run_fisher(simple_experiment, fisher_type):
    """Test run fisher"""

    filename, _ = simple_experiment
    result = runner.invoke(
        app,
        ["run", "fisher", filename.as_posix(), "--fisher-type", fisher_type],
    )

    assert result.exit_code == 0


def test_run_fisher_output(simple_experiment, fisher_type):
    """Test run fisher"""

    filename, _ = simple_experiment
    output = filename.with_suffix(".out.yaml")
    result = runner.invoke(
        app,
        [
            "run",
            "fisher",
            filename.as_posix(),
            "--fisher-type",
            fisher_type,
            "--output",
            output.as_posix(),
        ],
    )

    assert result.exit_code == 0

    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
    output_dict = ser.dict_str_from_yaml_file(output.as_posix())
    assert isinstance(output_dict, Ncm.ObjDictStr)
    v1 = output_dict.get("covariance")
    assert isinstance(v1, Ncm.Matrix)


def test_run_fisher_bias(simple_experiment):
    """Computes fisher bias based on a theory vector."""

    filename, _ = simple_experiment
    output = filename.with_suffix(".out.yaml")
    result = runner.invoke(
        app,
        ["run", "theory-vector", filename.as_posix(), "--output", output.as_posix()],
    )
    assert result.exit_code == 0

    result = runner.invoke(
        app,
        ["run", "fit", filename.as_posix(), "--output", output.as_posix()],
    )
    assert result.exit_code == 0

    result = runner.invoke(
        app,
        [
            "run",
            "fisher-bias",
            filename.as_posix(),
            "--theory-vector",
            output.as_posix(),
            "--output",
            output.as_posix(),
        ],
    )
    assert result.exit_code == 0

    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
    output_dict = ser.dict_str_from_yaml_file(output.as_posix())
    assert isinstance(output_dict, Ncm.ObjDictStr)
    v1 = output_dict.get("delta-theta")
    assert isinstance(v1, Ncm.Vector)

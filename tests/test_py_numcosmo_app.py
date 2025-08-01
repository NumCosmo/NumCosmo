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
from numcosmo_py.app.esmcmc import IniSampler, Parallelization
from numcosmo_py.app.generate import Planck18Types
from numcosmo_py.interpolation.stats_dist import CrossValidationMethod
from numcosmo_py.sampling import FitRunner, FitGradType, FitRunMessages, FisherType
from numcosmo_py.interpolation.stats_dist import (
    InterpolationMethod,
    InterpolationKernel,
)

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
    """Return a fit runner."""
    return request.param


@pytest.fixture(name="fit_grad_type", params=[e.value for e in FitGradType])
def fixture_fit_grad_type(request) -> str:
    """Return a fit grad type."""
    return request.param


@pytest.fixture(name="fit_run_messages", params=[e.value for e in FitRunMessages])
def fixture_fit_run_messages(request) -> str:
    """Return a fit run messages."""
    return request.param


@pytest.fixture(name="fisher_type", params=[e.value for e in FisherType])
def fixture_fisher_type(request) -> str:
    """Return fisher type."""
    return request.param


@pytest.fixture(
    name="interpolation_method", params=[e.value for e in InterpolationMethod]
)
def fixture_interpolation_method(request) -> str:
    """Return interpolation method."""
    return request.param


@pytest.fixture(
    name="interpolation_kernel", params=[e.value for e in InterpolationKernel]
)
def fixture_interpolation_kernel(request) -> str:
    """Return interpolation kernel."""
    return request.param


@pytest.fixture(
    name="calibration_method", params=[e.value for e in CrossValidationMethod]
)
def fixture_calibration_method(request) -> str:
    """Return calibration method."""
    return request.param


@pytest.fixture(name="planck18_type", params=[e.value for e in Planck18Types])
def fixture_planck18_type(request) -> str:
    """Return planck18 type."""
    return request.param


def test_run_test(simple_experiment):
    """Test run test."""
    filename, _ = simple_experiment
    result = runner.invoke(app, ["run", "test", filename.as_posix()])
    if result.exit_code != 0:
        raise result.exception


def test_run_fit(simple_experiment):
    """Test run fit."""
    filename, _ = simple_experiment
    result = runner.invoke(app, ["run", "fit", filename.as_posix()])
    if result.exit_code != 0:
        raise result.exception


def test_run_fit_runner(simple_experiment, fit_runner):
    """Test run fit."""
    filename, _ = simple_experiment
    result = runner.invoke(
        app, ["run", "fit", filename.as_posix(), "--runner", fit_runner]
    )
    if result.exit_code != 0:
        raise result.exception


def test_run_fit_grad_type(simple_experiment, fit_grad_type):
    """Test run fit."""
    filename, _ = simple_experiment
    result = runner.invoke(
        app, ["run", "fit", filename.as_posix(), "--grad-type", fit_grad_type]
    )
    if result.exit_code != 0:
        raise result.exception


def test_run_fit_run_messages(simple_experiment, fit_run_messages):
    """Test run fit."""
    filename, _ = simple_experiment
    result = runner.invoke(
        app, ["run", "fit", filename.as_posix(), "--run-messages", fit_run_messages]
    )
    if result.exit_code != 0:
        raise result.exception


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
    if result.exit_code != 0:
        raise result.exception


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
    if result.exit_code != 0:
        raise result.exception


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
    if result.exit_code != 0:
        raise result.exception

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

    if result.exit_code != 0:
        raise result.exception


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
    if result.exit_code != 0:
        raise result.exception


def test_run_fit_params_reltol(simple_experiment):
    """Test run fit with relative tolerance."""
    filename, _ = simple_experiment
    result = runner.invoke(
        app, ["run", "fit", filename.as_posix(), "--params-reltol", "1.0e-7"]
    )
    if result.exit_code != 0:
        raise result.exception


def test_run_fit_m2lnL_abstol(simple_experiment):
    """Test run fit with absolute tolerance."""
    filename, _ = simple_experiment
    result = runner.invoke(
        app, ["run", "fit", filename.as_posix(), "--m2lnl-abstol", "1.0e-7"]
    )
    if result.exit_code != 0:
        raise result.exception


def test_run_fit_m2lnL_reltol(simple_experiment):
    """Test run fit with relative tolerance."""
    filename, _ = simple_experiment
    result = runner.invoke(
        app, ["run", "fit", filename.as_posix(), "--m2lnl-reltol", "1.0e-7"]
    )
    if result.exit_code != 0:
        raise result.exception


def test_run_mc(simple_experiment):
    """Test run mc."""
    filename, _ = simple_experiment
    result = runner.invoke(app, ["run", "mc", "-p", filename.as_posix()])
    if result.exit_code != 0:
        raise result.exception


def test_run_mc_seed(simple_experiment):
    """Test run mc with seed."""
    filename, _ = simple_experiment
    result = runner.invoke(
        app, ["run", "mc", "-p", "--seed", "123", filename.as_posix()]
    )
    if result.exit_code != 0:
        raise result.exception


def test_run_theory_vector(simple_experiment):
    """Test run theory vector."""
    filename, _ = simple_experiment
    output = filename.with_suffix(".out.yaml")
    result = runner.invoke(
        app,
        ["run", "theory-vector", filename.as_posix(), "--output", output.as_posix()],
    )

    if result.exit_code != 0:
        raise result.exception

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

    if result.exit_code != 0:
        raise result.exception

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

    if result.exit_code != 0:
        raise result.exception

    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
    output_dict = ser.dict_str_from_yaml_file(output.as_posix())
    assert isinstance(output_dict, Ncm.ObjDictStr)
    v1 = output_dict.get("theory-vector")
    assert isinstance(v1, Ncm.Vector)


def test_run_fisher(simple_experiment, fisher_type):
    """Test run fisher."""
    filename, _ = simple_experiment
    result = runner.invoke(
        app,
        ["run", "fisher", filename.as_posix(), "--fisher-type", fisher_type],
    )

    if result.exit_code != 0:
        raise result.exception


def test_run_fisher_output(simple_experiment, fisher_type):
    """Test run fisher."""
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

    if result.exit_code != 0:
        raise result.exception

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

    if result.exit_code != 0:
        raise result.exception

    result = runner.invoke(
        app,
        ["run", "fit", filename.as_posix(), "--output", output.as_posix()],
    )
    if result.exit_code != 0:
        raise result.exception

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

    if result.exit_code != 0:
        raise result.exception

    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
    output_dict = ser.dict_str_from_yaml_file(output.as_posix())
    assert isinstance(output_dict, Ncm.ObjDictStr)
    v1 = output_dict.get("delta-theta")
    assert isinstance(v1, Ncm.Vector)


def test_run_mcmc_apes(simple_experiment):
    """Run a MCMC analysis using APES."""
    filename, _ = simple_experiment
    result = runner.invoke(app, ["run", "mcmc", "apes", filename.as_posix()])

    if result.exit_code != 0:
        raise result.exception


def test_run_mcmc_apes_threads(simple_experiment):
    """Run a MCMC analysis using APES."""
    filename, _ = simple_experiment
    result = runner.invoke(
        app,
        [
            "run",
            "mcmc",
            "apes",
            filename.as_posix(),
            "--parallel",
            Parallelization.THREADS.value,
        ],
    )

    if result.exit_code != 0:
        raise result.exception


def test_run_mcmc_apes_output(simple_experiment):
    """Run a MCMC analysis using APES."""
    filename, _ = simple_experiment
    output = filename.with_suffix(".out.yaml")
    result = runner.invoke(
        app,
        ["run", "mcmc", "apes", filename.as_posix(), "--output", output.as_posix()],
    )

    assert output.absolute().with_suffix(".mcmc.fits").exists()
    if result.exit_code != 0:
        raise result.exception


def test_run_mcmc_apes_init_gauss_cov(simple_experiment):
    """Run a MCMC analysis using APES starting at a best-fit and fisher matrix."""
    filename, _ = simple_experiment
    output = filename.with_suffix(".out.yaml")

    # Computes best-fit
    result = runner.invoke(
        app, ["run", "fit", filename.as_posix(), "--output", output.as_posix()]
    )
    if result.exit_code != 0:
        raise result.exception

    # Computes covariance
    result = runner.invoke(
        app,
        [
            "run",
            "fisher",
            filename.as_posix(),
            "--output",
            output.as_posix(),
        ],
    )
    if result.exit_code != 0:
        raise result.exception

    # Runs MCMC
    result = runner.invoke(
        app,
        [
            "run",
            "mcmc",
            "apes",
            filename.as_posix(),
            "--starting-point",
            output.as_posix(),
            "--output",
            output.as_posix(),
            "--initial-points-sampler",
            IniSampler.GAUSS_COV.value,
            "--initial-sampler-covar",
            output.as_posix(),
        ],
    )
    if result.exit_code != 0:
        raise result.exception


def test_run_mcmc_apes_init_catalog(simple_experiment):
    """Run a MCMC analysis using APES starting at previously computed catalog."""
    filename, _ = simple_experiment
    output = filename.with_suffix(".out.yaml")

    # Computes best-fit
    result = runner.invoke(
        app, ["run", "mcmc", "apes", filename.as_posix(), "--output", output.as_posix()]
    )
    if result.exit_code != 0:
        raise result.exception

    # Runs MCMC
    result = runner.invoke(
        app,
        [
            "run",
            "mcmc",
            "apes",
            filename.as_posix(),
            "--initial-points-sampler",
            IniSampler.FROM_CATALOG.value,
            "--initial-catalog",
            output.absolute().with_suffix(".mcmc.fits").as_posix(),
        ],
    )
    if result.exit_code != 0:
        raise result.exception


def test_run_mcmc_apes_method_kernel(
    simple_experiment, interpolation_method, interpolation_kernel
):
    """Run a MCMC analysis using APES."""
    filename, _ = simple_experiment
    output = filename.with_suffix(".out.yaml")
    result = runner.invoke(
        app,
        [
            "run",
            "mcmc",
            "apes",
            filename.as_posix(),
            "--output",
            output.as_posix(),
            "--interpolation-method",
            interpolation_method,
            "--interpolation-kernel",
            interpolation_kernel,
        ],
    )
    if result.exit_code != 0:
        raise result.exception


def test_run_mcmc_apes_method_kernel_no_interp(
    simple_experiment, interpolation_method, interpolation_kernel
):
    """Run a MCMC analysis using APES."""
    filename, _ = simple_experiment
    output = filename.with_suffix(".out.yaml")
    result = runner.invoke(
        app,
        [
            "run",
            "mcmc",
            "apes",
            filename.as_posix(),
            "--output",
            output.as_posix(),
            "--interpolation-method",
            interpolation_method,
            "--interpolation-kernel",
            interpolation_kernel,
            "--no-use-interpolation",
        ],
    )
    if result.exit_code != 0:
        raise result.exception


def test_run_mcmc_apes_analyze(simple_experiment):
    """Run a MCMC analysis using APES."""
    filename, _ = simple_experiment
    output = filename.with_suffix(".out.yaml")
    result = runner.invoke(
        app,
        ["run", "mcmc", "apes", filename.as_posix(), "--output", output.as_posix()],
    )

    assert output.absolute().with_suffix(".mcmc.fits").exists()
    if result.exit_code != 0:
        raise result.exception

    result = runner.invoke(
        app,
        [
            "catalog",
            "analyze",
            filename.as_posix(),
            output.absolute().with_suffix(".mcmc.fits").as_posix(),
        ],
    )

    if result.exit_code != 0:
        raise result.exception


def test_run_mcmc_apes_plot_corner(simple_experiment):
    """Run a MCMC analysis using APES."""
    filename, _ = simple_experiment
    output = filename.with_suffix(".out.yaml")
    result = runner.invoke(
        app,
        ["run", "mcmc", "apes", filename.as_posix(), "--output", output.as_posix()],
    )

    assert output.absolute().with_suffix(".mcmc.fits").exists()
    if result.exit_code != 0:
        raise result.exception

    result = runner.invoke(
        app,
        [
            "catalog",
            "plot-corner",
            filename.as_posix(),
            output.absolute().with_suffix(".mcmc.fits").as_posix(),
            "--no-show",
            "--output",
            output.absolute(),
        ],
    )

    assert output.absolute().with_suffix(".corner.pdf").exists()

    if result.exit_code != 0:
        raise result.exception


def test_run_mcmc_apes_analyze_evidence(simple_experiment):
    """Run a MCMC analysis using APES."""
    filename, _ = simple_experiment
    output = filename.with_suffix(".out.yaml")
    result = runner.invoke(
        app,
        ["run", "mcmc", "apes", filename.as_posix(), "--output", output.as_posix()],
    )

    assert output.absolute().with_suffix(".mcmc.fits").exists()
    if result.exit_code != 0:
        raise result.exception

    result = runner.invoke(
        app,
        [
            "catalog",
            "analyze",
            filename.as_posix(),
            output.absolute().with_suffix(".mcmc.fits").as_posix(),
            "--evidence",
        ],
    )

    if result.exit_code != 0:
        raise result.exception


def test_run_mcmc_apes_calibrate(simple_experiment, calibration_method):
    """Run a MCMC analysis using APES."""
    filename, _ = simple_experiment
    output = filename.with_suffix(".out.yaml")
    result = runner.invoke(
        app,
        ["run", "mcmc", "apes", filename.as_posix(), "--output", output.as_posix()],
    )

    assert output.absolute().with_suffix(".mcmc.fits").exists()
    if result.exit_code != 0:
        raise result.exception

    result = runner.invoke(
        app,
        [
            "catalog",
            "calibrate",
            filename.as_posix(),
            output.absolute().with_suffix(".mcmc.fits").as_posix(),
            "--cv-method",
            calibration_method,
        ],
    )

    if result.exit_code != 0:
        raise result.exception


def test_generate_planck(tmp_path, planck18_type):
    """Test run theory vector."""
    tmp_file = tmp_path / "planck_generated.yaml"

    result = runner.invoke(
        app,
        ["generate", "planck18", tmp_file.as_posix(), "--data-type", planck18_type],
    )

    if result.exit_code != 0:
        raise result.exception


def test_generate_planck_test(tmp_path, planck18_type):
    """Test run theory vector."""
    tmp_file = tmp_path / "planck_generated.yaml"

    result = runner.invoke(
        app,
        ["generate", "planck18", tmp_file.as_posix(), "--data-type", planck18_type],
    )

    if result.exit_code != 0:
        raise result.exception

    result = runner.invoke(app, ["run", "test", tmp_file.as_posix()])

    if result.exit_code != 0:
        raise result.exception

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
import typer
from typer.testing import CliRunner

pytest.importorskip("astropy")
pytest.importorskip("getdist")
# flake8: noqa: E402
# pylint: disable=wrong-import-position


from numcosmo_py import Ncm, Nc
from numcosmo_py.app.loading import load_catalog
from numcosmo_py.app import app
from numcosmo_py.app.esmcmc import IniSampler, Parallelization
from numcosmo_py.app.generate import Planck18Types
from numcosmo_py.interpolation.stats_dist import CrossValidationMethod
from numcosmo_py.sampling import FitRunner, FitGradType, FitRunMessages, FisherType
from numcosmo_py.interpolation.stats_dist import (
    InterpolationMethod,
    InterpolationKernel,
)

pytestmark = pytest.mark.app
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


@pytest.mark.parametrize("subcommand", [["run", "test"], ["run", "fit"]])
def test_run_bench(simple_experiment, subcommand):
    """--bench reports wall-clock time and peak RSS at the end of the run."""
    filename, _ = simple_experiment
    result = runner.invoke(app, [*subcommand, filename.as_posix(), "--bench"])
    if result.exit_code != 0:
        raise result.exception
    assert "# Benchmark:" in result.output
    assert "wall time" in result.output
    assert "peak RSS" in result.output


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


def test_run_mc_fiducial(simple_experiment, tmp_path):
    """run mc resamples from an external fiducial mset (--fiducial)."""
    filename, experiment = simple_experiment

    # A fiducial mset with shifted parameters, serialized to YAML.
    fiducial = experiment.get("model-set").dup(Ncm.Serialize.new(0))
    fiducial.param_set_all_ftype(Ncm.ParamType.FREE)
    fiducial.fparams_set_array([0.5] * fiducial.fparam_len())
    fiducial_file = tmp_path / "fiducial.mset.yaml"
    Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP).to_yaml_file(
        fiducial, fiducial_file.as_posix()
    )

    result = runner.invoke(
        app,
        [
            "run",
            "mc",
            "-p",
            "--seed",
            "123",
            "--nmc",
            "5",
            "--fiducial",
            fiducial_file.as_posix(),
            filename.as_posix(),
        ],
    )
    if result.exit_code != 0:
        raise result.exception


def test_run_mc_fiducial_not_an_mset(simple_experiment, tmp_path):
    """run mc rejects a --fiducial file that is not an NcmMSet."""
    filename, _ = simple_experiment

    # A serialized model (not wrapped in an MSet) must be rejected.
    not_mset = tmp_path / "not_mset.yaml"
    Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP).to_yaml_file(
        Ncm.ModelMVND.new(2), not_mset.as_posix()
    )

    result = runner.invoke(
        app,
        ["run", "mc", "-p", "--fiducial", not_mset.as_posix(), filename.as_posix()],
    )
    assert result.exit_code != 0


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
            output.absolute().with_suffix(".mcmc.fits").as_posix(),
        ],
    )

    if result.exit_code != 0:
        raise result.exception


def test_run_mcmc_apes_analyze_burnin_iterations(simple_experiment):
    """--burnin is interpreted in iterations (ensemble steps), not raw rows."""
    filename, _ = simple_experiment
    output = filename.with_suffix(".out.yaml")
    result = runner.invoke(
        app,
        ["run", "mcmc", "apes", filename.as_posix(), "--output", output.as_posix()],
    )
    if result.exit_code != 0:
        raise result.exception

    catalog = output.absolute().with_suffix(".mcmc.fits")

    result = runner.invoke(
        app, ["catalog", "analyze", catalog.as_posix(), "--burnin", "3"]
    )
    if result.exit_code != 0:
        raise result.exception


def test_run_mcmc_apes_analyze_burnin_beyond_catalog_size(simple_experiment):
    """--burnin larger than the catalog's iteration count fails with a clear,
    catchable error instead of aborting the process."""
    filename, _ = simple_experiment
    output = filename.with_suffix(".out.yaml")
    result = runner.invoke(
        app,
        ["run", "mcmc", "apes", filename.as_posix(), "--output", output.as_posix()],
    )
    if result.exit_code != 0:
        raise result.exception

    catalog = output.absolute().with_suffix(".mcmc.fits")

    result = runner.invoke(
        app, ["catalog", "analyze", catalog.as_posix(), "--burnin", "1000000"]
    )
    assert result.exit_code != 0
    assert "exceeds catalog" in result.output


def test_run_mcmc_apes_analyze_tail(simple_experiment):
    """--tail keeps only the last N iterations."""
    filename, _ = simple_experiment
    output = filename.with_suffix(".out.yaml")
    result = runner.invoke(
        app,
        ["run", "mcmc", "apes", filename.as_posix(), "--output", output.as_posix()],
    )
    if result.exit_code != 0:
        raise result.exception

    catalog = output.absolute().with_suffix(".mcmc.fits")

    result = runner.invoke(
        app, ["catalog", "analyze", catalog.as_posix(), "--tail", "2"]
    )
    if result.exit_code != 0:
        raise result.exception


def test_run_mcmc_apes_analyze_burnin_and_tail_rejected(simple_experiment):
    """--burnin and --tail are mutually exclusive."""
    filename, _ = simple_experiment
    output = filename.with_suffix(".out.yaml")
    result = runner.invoke(
        app,
        ["run", "mcmc", "apes", filename.as_posix(), "--output", output.as_posix()],
    )
    if result.exit_code != 0:
        raise result.exception

    catalog = output.absolute().with_suffix(".mcmc.fits")

    result = runner.invoke(
        app,
        [
            "catalog",
            "analyze",
            catalog.as_posix(),
            "--burnin",
            "1",
            "--tail",
            "2",
        ],
    )
    assert result.exit_code != 0


def test_run_mcmc_apes_analyze_exclude_only(simple_experiment):
    """--exclude alone drops only matching columns."""
    filename, _ = simple_experiment
    output = filename.with_suffix(".out.yaml")
    result = runner.invoke(
        app,
        ["run", "mcmc", "apes", filename.as_posix(), "--output", output.as_posix()],
    )
    if result.exit_code != 0:
        raise result.exception

    catalog = output.absolute().with_suffix(".mcmc.fits")
    result = runner.invoke(
        app,
        ["catalog", "analyze", catalog.as_posix(), "--exclude", "mu_4"],
    )
    if result.exit_code != 0:
        raise result.exception
    assert "mu_0" in result.output
    assert "mu_4" not in result.output


def test_run_mcmc_apes_analyze_include_only(simple_experiment):
    """--include alone selects only matching columns."""
    filename, _ = simple_experiment
    output = filename.with_suffix(".out.yaml")
    result = runner.invoke(
        app,
        ["run", "mcmc", "apes", filename.as_posix(), "--output", output.as_posix()],
    )
    if result.exit_code != 0:
        raise result.exception

    catalog = output.absolute().with_suffix(".mcmc.fits")
    result = runner.invoke(
        app,
        [
            "catalog",
            "analyze",
            catalog.as_posix(),
            "--include",
            "mu_0",
            "--include",
            "mu_1",
        ],
    )
    if result.exit_code != 0:
        raise result.exception
    assert "mu_0" in result.output
    assert "mu_1" in result.output
    assert "mu_2" not in result.output


def test_run_mcmc_apes_analyze_include_and_exclude(simple_experiment):
    """--include and --exclude combined narrow the selection further."""
    filename, _ = simple_experiment
    output = filename.with_suffix(".out.yaml")
    result = runner.invoke(
        app,
        ["run", "mcmc", "apes", filename.as_posix(), "--output", output.as_posix()],
    )
    if result.exit_code != 0:
        raise result.exception

    catalog = output.absolute().with_suffix(".mcmc.fits")
    result = runner.invoke(
        app,
        [
            "catalog",
            "analyze",
            catalog.as_posix(),
            "--include",
            "mu",
            "--exclude",
            "mu_4",
        ],
    )
    if result.exit_code != 0:
        raise result.exception
    assert "mu_0" in result.output
    assert "mu_4" not in result.output


def test_run_mcmc_apes_analyze_empty_catalog(simple_experiment):
    """--burnin consuming the whole catalog leaves an empty catalog, handled
    cleanly instead of crashing."""
    filename, _ = simple_experiment
    output = filename.with_suffix(".out.yaml")
    result = runner.invoke(
        app,
        ["run", "mcmc", "apes", filename.as_posix(), "--output", output.as_posix()],
    )
    if result.exit_code != 0:
        raise result.exception

    catalog = output.absolute().with_suffix(".mcmc.fits")
    nrows, nchains, _first_id = Ncm.MSetCatalog.peek_info_from_file(catalog.as_posix())
    n_iterations = nrows // nchains

    result = runner.invoke(
        app, ["catalog", "analyze", catalog.as_posix(), "--burnin", str(n_iterations)]
    )
    if result.exit_code != 0:
        raise result.exception
    assert "Empty catalog" in result.output


def test_run_mcmc_apes_plot_corner_too_many_plot_names(simple_experiment):
    """More --plot-name values than catalogs is a clean, catchable error."""
    filename, _ = simple_experiment
    output = filename.with_suffix(".out.yaml")
    result = runner.invoke(
        app,
        ["run", "mcmc", "apes", filename.as_posix(), "--output", output.as_posix()],
    )
    if result.exit_code != 0:
        raise result.exception

    catalog = output.absolute().with_suffix(".mcmc.fits")
    result = runner.invoke(
        app,
        [
            "catalog",
            "plot-corner",
            catalog.as_posix(),
            "--plot-name",
            "A",
            "--plot-name",
            "B",
            "--no-show",
        ],
    )
    assert result.exit_code != 0
    # Avoid asserting on a substring straddling "--plot-name": Rich highlights
    # option-looking tokens and injects ANSI codes between their characters
    # when color is on (as it is in CI, unlike a plain local terminal), which
    # would otherwise break a naive substring match.
    assert "values than catalog files" in result.output


def test_run_mcmc_apes_plot_corner_mark_bestfit(simple_experiment):
    """--mark-bestfit marks the first catalog's best-fit point."""
    filename, _ = simple_experiment
    output = filename.with_suffix(".out.yaml")
    result = runner.invoke(
        app,
        ["run", "mcmc", "apes", filename.as_posix(), "--output", output.as_posix()],
    )
    if result.exit_code != 0:
        raise result.exception

    catalog = output.absolute().with_suffix(".mcmc.fits")
    result = runner.invoke(
        app,
        [
            "catalog",
            "plot-corner",
            catalog.as_posix(),
            "--mark-bestfit",
            "--no-show",
            "--output",
            output.as_posix(),
        ],
    )
    if result.exit_code != 0:
        raise result.exception
    assert output.with_suffix(".corner.pdf").exists()


def test_catalog_visual_hw(simple_experiment, monkeypatch):
    """catalog visual-hw runs end-to-end.

    VisualHW has no --no-show option (unlike plot-corner) and calls
    plt.show() unconditionally, so it's stubbed out here to keep the test
    from blocking on a display.
    """
    monkeypatch.setattr("matplotlib.pyplot.show", lambda *args, **kwargs: None)

    filename, _ = simple_experiment
    output = filename.with_suffix(".out.yaml")
    result = runner.invoke(
        app,
        ["run", "mcmc", "apes", filename.as_posix(), "--output", output.as_posix()],
    )
    if result.exit_code != 0:
        raise result.exception

    catalog = output.absolute().with_suffix(".mcmc.fits")
    result = runner.invoke(
        app, ["catalog", "visual-hw", catalog.as_posix(), "--param-name", "mu_0"]
    )
    if result.exit_code != 0:
        raise result.exception


def test_catalog_param_evolution(simple_experiment, monkeypatch):
    """catalog param-evolution runs end-to-end. See test_catalog_visual_hw
    for why plt.show() is stubbed out."""
    monkeypatch.setattr("matplotlib.pyplot.show", lambda *args, **kwargs: None)

    filename, _ = simple_experiment
    output = filename.with_suffix(".out.yaml")
    result = runner.invoke(
        app,
        ["run", "mcmc", "apes", filename.as_posix(), "--output", output.as_posix()],
    )
    if result.exit_code != 0:
        raise result.exception

    catalog = output.absolute().with_suffix(".mcmc.fits")
    result = runner.invoke(
        app,
        [
            "catalog",
            "param-evolution",
            catalog.as_posix(),
            "--param-name",
            "mu_0",
            "--grid-size",
            "20",
        ],
    )
    if result.exit_code != 0:
        raise result.exception


def test_run_mc_analyze_single_chain(simple_experiment):
    """catalog analyze on a single-chain (run mc) catalog exercises the
    nchains == 1 branch of load_catalog()'s stats selection."""
    filename, _ = simple_experiment
    output = filename.with_suffix(".mc_out.yaml")
    result = runner.invoke(
        app,
        [
            "run",
            "mc",
            filename.as_posix(),
            "--output",
            output.as_posix(),
            "--nmc",
            "20",
        ],
    )
    if result.exit_code != 0:
        raise result.exception

    catalog = output.absolute().with_suffix(".mc.fits")
    assert catalog.exists()

    result = runner.invoke(app, ["catalog", "analyze", catalog.as_posix()])
    if result.exit_code != 0:
        raise result.exception


def test_catalog_analyze_missing_file(tmp_path):
    """catalog analyze on a nonexistent file fails with a clean, catchable
    error instead of a traceback."""
    result = runner.invoke(
        app, ["catalog", "analyze", (tmp_path / "does_not_exist.fits").as_posix()]
    )
    assert result.exit_code != 0
    assert "not found" in result.output


def test_load_catalog_negative_tail_rejected(simple_experiment):
    """load_catalog()'s --tail validation rejects a negative value. Only
    reachable via a direct call: the CLI's own min=0 constraint on --tail
    blocks this before load_catalog() is ever invoked."""
    filename, _ = simple_experiment
    output = filename.with_suffix(".out.yaml")
    result = runner.invoke(
        app,
        ["run", "mcmc", "apes", filename.as_posix(), "--output", output.as_posix()],
    )
    if result.exit_code != 0:
        raise result.exception

    catalog = output.absolute().with_suffix(".mcmc.fits")

    with pytest.raises(typer.BadParameter):
        load_catalog(catalog, tail=-1)


def test_run_mcmc_apes_embeds_functions_array(tmp_path):
    """A run's --functions (extra columns) end up embedded in the catalog
    file's HDU0, readable back with no companion .functions.yaml needed."""
    cosmo = Nc.HICosmoDEXcdm.new()
    cosmo.param_set_by_name("H0", 70.0)
    cosmo.param_set_ftype(cosmo.param_index_from_name("H0")[1], Ncm.ParamType.FREE)
    mset = Ncm.MSet.new_array([cosmo])
    mset.prepare_fparam_map()

    data = Nc.DataHubble.new_from_id(Nc.DataHubbleId.SIMON2005)
    likelihood = Ncm.Likelihood.new(dset=Ncm.Dataset.new_array([data]))

    experiment = tmp_path / "exp.yaml"
    ser = Ncm.Serialize.new(Ncm.SerializeOpt.NONE)
    experiment_dict = Ncm.ObjDictStr.new()
    experiment_dict.add("model-set", mset)
    experiment_dict.add("likelihood", likelihood)
    ser.dict_str_to_yaml_file(experiment_dict, experiment.as_posix())

    # A constant (nvar=0) derived quantity, the kind --functions is for.
    func = Ncm.MSetFuncList.new("NcHICosmoDE:wDE", None)
    funcs_oa = Ncm.ObjArray.new()
    funcs_oa.add(func)
    ser.reset(True)
    ser.array_to_yaml_file(
        funcs_oa, experiment.with_suffix(".functions.yaml").as_posix()
    )

    output = experiment.with_suffix(".out.yaml")
    result = runner.invoke(
        app,
        ["run", "mcmc", "apes", experiment.as_posix(), "--output", output.as_posix()],
    )
    if result.exit_code != 0:
        raise result.exception

    catalog = output.absolute().with_suffix(".mcmc.fits")
    mcat = Ncm.MSetCatalog.new_from_file_ro(catalog.as_posix(), 0)
    functions = mcat.peek_functions_array()

    assert functions is not None
    assert functions.len() == 1
    assert "wDE" in [mcat.col_name(i) for i in range(mcat.ncols())]


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
            output.absolute().with_suffix(".mcmc.fits").as_posix(),
            "--no-show",
            "--output",
            output.absolute(),
        ],
    )

    assert output.absolute().with_suffix(".corner.pdf").exists()

    if result.exit_code != 0:
        raise result.exception


def test_run_mcmc_apes_plot_corner_multi_catalog(simple_experiment):
    """plot-corner overlays multiple catalogs given positionally, with
    --plot-name controlling each one's legend entry."""
    filename, _ = simple_experiment

    catalogs = []
    for name in ("a", "b"):
        output = filename.with_suffix(f".out_{name}.yaml")
        result = runner.invoke(
            app,
            [
                "run",
                "mcmc",
                "apes",
                filename.as_posix(),
                "--output",
                output.as_posix(),
            ],
        )
        if result.exit_code != 0:
            raise result.exception
        catalogs.append(output.absolute().with_suffix(".mcmc.fits"))
        assert catalogs[-1].exists()

    plot_output = filename.with_suffix(".corner_multi")
    result = runner.invoke(
        app,
        [
            "catalog",
            "plot-corner",
            catalogs[0].as_posix(),
            catalogs[1].as_posix(),
            "--plot-name",
            "Run A",
            "--plot-name",
            "Run B",
            "--no-show",
            "--output",
            plot_output.as_posix(),
        ],
    )
    if result.exit_code != 0:
        raise result.exception

    assert plot_output.with_suffix(".corner.pdf").exists()


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
            output.absolute().with_suffix(".mcmc.fits").as_posix(),
            "--cv-method",
            calibration_method,
        ],
    )

    if result.exit_code != 0:
        raise result.exception


def test_run_mcmc_apes_get_best_fit(simple_experiment):
    """get-best-fit writes the catalog's best-fit parameters as a
    model-set yaml, dict-wrapped so it can be used as --starting-point."""
    filename, _ = simple_experiment
    output = filename.with_suffix(".out.yaml")
    result = runner.invoke(
        app,
        ["run", "mcmc", "apes", filename.as_posix(), "--output", output.as_posix()],
    )
    if result.exit_code != 0:
        raise result.exception

    catalog = output.absolute().with_suffix(".mcmc.fits")
    best_fit_output = filename.with_suffix(".best_fit.yaml")

    result = runner.invoke(
        app,
        [
            "catalog",
            "get-best-fit",
            catalog.as_posix(),
            "--output",
            best_fit_output.as_posix(),
        ],
    )
    if result.exit_code != 0:
        raise result.exception

    assert best_fit_output.exists()

    ser = Ncm.Serialize.new(Ncm.SerializeOpt.NONE)
    saved = ser.dict_str_from_yaml_file(best_fit_output.absolute().as_posix())
    mset = saved.get("model-set")
    assert isinstance(mset, Ncm.MSet)


def test_run_mcmc_apes_get_best_fit_requires_output(simple_experiment):
    """get-best-fit without --output fails with a clean, catchable error."""
    filename, _ = simple_experiment
    output = filename.with_suffix(".out.yaml")
    result = runner.invoke(
        app,
        ["run", "mcmc", "apes", filename.as_posix(), "--output", output.as_posix()],
    )
    if result.exit_code != 0:
        raise result.exception

    catalog = output.absolute().with_suffix(".mcmc.fits")

    result = runner.invoke(app, ["catalog", "get-best-fit", catalog.as_posix()])
    assert result.exit_code != 0
    assert "Output file not defined" in result.output


def test_generate_planck(tmp_path, planck18_type):
    """Test run theory vector."""
    tmp_file = tmp_path / "planck_generated1.yaml"

    result = runner.invoke(
        app,
        ["generate", "planck18", tmp_file.as_posix(), "--data-type", planck18_type],
    )

    if result.exit_code != 0:
        raise result.exception


def test_generate_planck_test(tmp_path, planck18_type):
    """Test run theory vector."""
    tmp_file = tmp_path / "planck_generated2.yaml"

    result = runner.invoke(
        app,
        ["generate", "planck18", tmp_file.as_posix(), "--data-type", planck18_type],
    )

    if result.exit_code != 0:
        raise result.exception

    result = runner.invoke(app, ["run", "test", tmp_file.as_posix()])

    if result.exit_code != 0:
        raise result.exception

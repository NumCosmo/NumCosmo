#!/usr/bin/env python
#
# test_inspect_galaxy_shape_integrand.py
#
# Sat Aug 1 2026
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
# with this program. If not, see <http://www.gnu.org/licenses/>.

"""CLI-level tests for the `inspect galaxy-shape-integrand` command.

Builds small mock cluster-WL experiments (via `generate cluster-wl`) and
exercises the integrand-heatmap CLI through the actual Typer app: the
different `--component` panels, the Spectral-only `--show-panels` overlay,
both ellipticity conventions, and the error paths (`_find_wl_factor()` /
`_find_galaxy_shape_pop()` lookups, an out-of-range `--galaxy-index`).
"""

from pathlib import Path

import pytest
from typer.testing import CliRunner

pytest.importorskip("astropy")
pytest.importorskip("getdist")
# flake8: noqa: E402

from numcosmo_py import Ncm
from numcosmo_py.app import app

pytestmark = pytest.mark.app
runner = CliRunner()

Ncm.cfg_init()


def _generate_wl_experiment(tmp_path: Path, name: str, shape_factor: str) -> Path:
    """Generate a small mock cluster-WL experiment with the given
    --shape-factor spec and return its experiment file path."""
    experiment_file = tmp_path / name
    result = runner.invoke(
        app,
        [
            "generate",
            "cluster-wl",
            experiment_file.as_posix(),
            "--shape-factor",
            shape_factor,
        ],
    )
    assert result.exit_code == 0, result.output
    return experiment_file


@pytest.fixture(name="spectral_experiment")
def fixture_spectral_experiment(tmp_path: Path) -> Path:
    """A mock cluster-WL experiment using the Spectral shape factor (the engine
    eval_panels()/--show-panels works with)."""
    return _generate_wl_experiment(tmp_path, "spectral_experiment.yaml", "spectral")


@pytest.fixture(name="simple_experiment")
def fixture_simple_experiment(tmp_path: Path) -> Path:
    """A generic experiment with no WL data at all, for the
    _find_wl_factor() failure path."""
    rng = Ncm.RNG.seeded_new(None, 1234)
    model_mvnd = Ncm.ModelMVND.new(5)
    mset = Ncm.MSet.new_array([model_mvnd])
    mset.param_set_all_ftype(Ncm.ParamType.FREE)

    data_mvnd = Ncm.DataGaussCovMVND.new_full(5, 1.0, 2.0, 30.0, 0.0, 0.0, rng)
    likelihood = Ncm.Likelihood.new(dset=Ncm.Dataset.new_array([data_mvnd]))

    exp_file = tmp_path / "simple_experiment.yaml"
    ser = Ncm.Serialize.new(Ncm.SerializeOpt.KEEP_SHARED)

    experiment = Ncm.ObjDictStr.new()
    experiment.add("model-set", mset)
    experiment.add("likelihood", likelihood)

    ser.dict_str_to_yaml_file(experiment, exp_file.as_posix())

    return exp_file


def test_galaxy_shape_integrand_default_saves_total(
    spectral_experiment: Path, tmp_path: Path
):
    """Default --component (total) saves a single figure."""
    out_prefix = tmp_path / "default"

    result = runner.invoke(
        app,
        [
            "inspect",
            "galaxy-shape-integrand",
            spectral_experiment.as_posix(),
            "--no-show-plot",
            "--output-prefix",
            out_prefix.as_posix(),
            "--n-grid",
            "25",
        ],
    )

    assert result.exit_code == 0, result.output
    assert (tmp_path / "default_shape_integrand.png").exists()


@pytest.mark.parametrize("component", ["noise", "population", "total"])
def test_galaxy_shape_integrand_single_component_panels(
    spectral_experiment: Path, tmp_path: Path, component: str
):
    """Each single-panel --component value renders and saves its own
    figure, with its own output filename suffix."""
    out_prefix = tmp_path / f"c_{component}"

    result = runner.invoke(
        app,
        [
            "inspect",
            "galaxy-shape-integrand",
            spectral_experiment.as_posix(),
            "--no-show-plot",
            "--component",
            component,
            "--output-prefix",
            out_prefix.as_posix(),
            "--n-grid",
            "20",
        ],
    )

    assert result.exit_code == 0, result.output
    suffix = (
        "shape_integrand" if component == "total" else f"shape_integrand_{component}"
    )
    assert (tmp_path / f"c_{component}_{suffix}.png").exists()


def test_galaxy_shape_integrand_all_component_renders_three_panels(
    spectral_experiment: Path, tmp_path: Path
):
    """--component all draws all three panels side by side under one
    shared color scale and saves a single "_all" figure."""
    out_prefix = tmp_path / "all"

    result = runner.invoke(
        app,
        [
            "inspect",
            "galaxy-shape-integrand",
            spectral_experiment.as_posix(),
            "--no-show-plot",
            "--component",
            "all",
            "--output-prefix",
            out_prefix.as_posix(),
            "--n-grid",
            "20",
        ],
    )

    assert result.exit_code == 0, result.output
    assert (tmp_path / "all_shape_integrand_all.png").exists()


def test_galaxy_shape_integrand_spectral_overlays_panels(
    spectral_experiment: Path, tmp_path: Path
):
    """With the Spectral engine the accepted panels are fetched and drawn and their
    count is printed."""
    out_prefix = tmp_path / "spectral"

    result = runner.invoke(
        app,
        [
            "inspect",
            "galaxy-shape-integrand",
            spectral_experiment.as_posix(),
            "--no-show-plot",
            "--output-prefix",
            out_prefix.as_posix(),
            "--n-grid",
            "25",
            "--g1",
            "0.2",
            "--g2",
            "0.1",
        ],
    )

    assert result.exit_code == 0, result.output
    assert (tmp_path / "spectral_shape_integrand.png").exists()
    assert "Spectral panels:" in result.output


def test_galaxy_shape_integrand_spectral_no_show_panels(spectral_experiment: Path):
    """--no-show-panels does not fetch the panels."""
    result = runner.invoke(
        app,
        [
            "inspect",
            "galaxy-shape-integrand",
            spectral_experiment.as_posix(),
            "--no-show-plot",
            "--no-show-panels",
            "--n-grid",
            "20",
        ],
    )

    assert result.exit_code == 0, result.output
    assert "Spectral panels:" not in result.output


@pytest.mark.parametrize(
    "g1, expected",
    [("1.0", "closed form and no panels"), ("1.5", "Spectral panels:")],
)
def test_galaxy_shape_integrand_spectral_large_shear(
    spectral_experiment: Path, g1: str, expected: str
):
    """|g| = 1 has no panels and says so; |g| > 1 draws the panels of the
    integral at g / |g|^2."""
    result = runner.invoke(
        app,
        [
            "inspect",
            "galaxy-shape-integrand",
            spectral_experiment.as_posix(),
            "--no-show-plot",
            "--n-grid",
            "20",
            "--g1",
            g1,
        ],
    )

    assert result.exit_code == 0, result.output
    assert expected in result.output


def test_galaxy_shape_integrand_trace_ellip_conv(tmp_path: Path):
    """The TRACE ellipticity convention selects the *_trace_ptr shear
    helpers instead of the *_trace_det_ptr ones."""
    experiment = _generate_wl_experiment(
        tmp_path, "trace_experiment.yaml", "spectral ellip_conv=trace"
    )

    result = runner.invoke(
        app,
        [
            "inspect",
            "galaxy-shape-integrand",
            experiment.as_posix(),
            "--no-show-plot",
            "--n-grid",
            "20",
        ],
    )

    assert result.exit_code == 0, result.output


def test_galaxy_shape_integrand_galaxy_index_out_of_range(spectral_experiment: Path):
    """An out-of-range --galaxy-index raises a clear error."""
    result = runner.invoke(
        app,
        [
            "inspect",
            "galaxy-shape-integrand",
            spectral_experiment.as_posix(),
            "--no-show-plot",
            "--galaxy-index",
            "10000000",
        ],
    )

    assert result.exit_code != 0
    assert result.exception is not None
    assert "out of range" in str(result.exception)


def test_galaxy_shape_integrand_requires_wl_factor(simple_experiment: Path):
    """_find_wl_factor() raises for an experiment with no
    Nc.DataClusterWLFactor in its dataset."""
    result = runner.invoke(
        app,
        [
            "inspect",
            "galaxy-shape-integrand",
            simple_experiment.as_posix(),
            "--no-show-plot",
        ],
    )

    assert result.exit_code != 0
    assert result.exception is not None
    assert "No Nc.DataClusterWLFactor data object found in dataset." in str(
        result.exception
    )


def test_galaxy_shape_integrand_requires_galaxy_shape_pop(
    spectral_experiment: Path, monkeypatch: pytest.MonkeyPatch
):
    """_find_galaxy_shape_pop() raises when no Nc.GalaxyShapePop model is
    present in the model set, even though a WL factor is."""
    monkeypatch.setattr(Ncm.MSet, "peek_array_pos", lambda self, i: None)

    result = runner.invoke(
        app,
        [
            "inspect",
            "galaxy-shape-integrand",
            spectral_experiment.as_posix(),
            "--no-show-plot",
        ],
    )

    assert result.exit_code != 0
    assert result.exception is not None
    assert "No Nc.GalaxyShapePop model found in the model set." in str(result.exception)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

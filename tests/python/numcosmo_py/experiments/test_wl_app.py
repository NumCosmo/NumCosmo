#!/usr/bin/env python
#
# test_py_cluster_wl_app.py
#
# Fri Jun 13 14:02:22 2025
# Copyright  2025  Caio Lima de Oliveira
# <caiolimadeoliveira@pm.me>
#
# test_py_cluster_wl_app.py
# Copyright (C) 2025 Caio Lima de Oliveira <caiolimadeoliveira@pm.me>
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

"""Test NumCosmo app's Cluster WL module.

Targets the ``NcGalaxy*Factor``/``NcDataClusterWLFactor`` pipeline.
"""

from typing import cast
from filecmp import cmp
from pathlib import Path
import pytest
from typer.testing import CliRunner

pytest.importorskip("astropy")
pytest.importorskip("getdist")
# flake8: noqa: E402
# pylint: disable=wrong-import-position

import numpy as np
from numpy import sin, cos, pi, log10, deg2rad
from numpy.random import uniform, choice

from numcosmo_py import (
    Ncm,
    Nc,
)
from numcosmo_py.app import app
from numcosmo_py.experiments.cluster_wl import (
    GalaxyZGen,
    ShapeFactorGen,
    GalaxyPopGen,
    HaloProfileType,
    ClusterModel,
    HaloPositionData,
    GalaxyPopGenBeta,
)

pytestmark = pytest.mark.app

runner = CliRunner()

Ncm.cfg_init()

# Columns written by NcDataClusterWLFactor.resample(): NOT the legacy
# sigma_0/sigma_z spelling -- just sigma0 (a fixed per-galaxy calibration
# input), no separately-stored sigma_z (it's sigma0*(1+z) internally).
_VARIANTS = ["y1-source", "y1-lens", "y10-source", "y10-lens"]


@pytest.fixture(name="experiment_file")
def fixture_experiment_file(tmp_path: Path) -> Path:
    """Fixture for the experiment file."""
    tmp_file = tmp_path / "experiment.yaml"
    return tmp_file


@pytest.fixture(name="experiment_file_copy")
def fixture_experiment_file_copy(tmp_path: Path) -> Path:
    """Fixture for the experiment file copy."""
    tmp_file = tmp_path / "experiment_copy.yaml"
    return tmp_file


def _build_real_wl_obs(tmp_path: Path, filename: str, with_meta: bool) -> Path:
    """A small synthetic NcGalaxyWLObs, standing in for a real catalog.

    Mirrors the column schema (and per-galaxy pz spline) of the curated
    Subaru HSC-SSP PDR1 catalogs shipped in the NumCosmo data file release,
    without requiring network access in tests. If with_meta, also carries the
    cluster_ra/cluster_dec/cluster_z/cluster_c catalog metadata those real
    catalogs are expected to ship with.
    """
    rng = np.random.default_rng(20260720)
    n = 20
    cols = [
        "epsilon_int_1",
        "epsilon_int_2",
        "epsilon_obs_1",
        "epsilon_obs_2",
        "std_noise",
        "c1",
        "c2",
        "m",
        "ra",
        "dec",
        "z",
    ]
    obs = Nc.GalaxyWLObs.new(
        Nc.GalaxyWLObsEllipConv.TRACE, Nc.WLEllipticityFrame.CELESTIAL, n, cols
    )
    ra = rng.uniform(-0.05, 0.05, n)
    dec = rng.uniform(-0.05, 0.05, n)
    z = rng.uniform(0.5, 1.5, n)
    epsilon_obs_1 = rng.normal(0.0, 0.2, n)
    epsilon_obs_2 = rng.normal(0.0, 0.2, n)
    std_noise = rng.uniform(0.1, 0.2, n)

    for i in range(n):
        obs.set("epsilon_int_1", i, 0.0)
        obs.set("epsilon_int_2", i, 0.0)
        obs.set("epsilon_obs_1", i, float(epsilon_obs_1[i]))
        obs.set("epsilon_obs_2", i, float(epsilon_obs_2[i]))
        obs.set("std_noise", i, float(std_noise[i]))
        obs.set("c1", i, 0.0)
        obs.set("c2", i, 0.0)
        obs.set("m", i, 0.0)
        obs.set("ra", i, float(ra[i]))
        obs.set("dec", i, float(dec[i]))
        obs.set("z", i, float(z[i]))

        zs = np.linspace(max(z[i] - 0.3, 0.01), z[i] + 0.3, 16)
        pz_vals = np.exp(-0.5 * ((zs - z[i]) / 0.05) ** 2)
        spline = Ncm.SplineCubicNotaknot.new()
        spline.set(
            Ncm.Vector.new_array(zs.tolist()),
            Ncm.Vector.new_array(pz_vals.tolist()),
            True,
        )
        obs.set_pz(i, spline)

    if with_meta:
        meta = Ncm.VarDict.new()
        meta.set_double("cluster_ra", 0.0)
        meta.set_double("cluster_dec", 0.0)
        meta.set_double("cluster_z", 0.8)
        meta.set_double("cluster_c", 4.0)
        obs.set_meta(meta)

    obs_file = tmp_path / filename
    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
    ser.to_binfile(obs, obs_file.absolute().as_posix())

    return obs_file


@pytest.fixture(name="real_wl_obs_file")
def fixture_real_wl_obs_file(tmp_path: Path) -> Path:
    """A synthetic NcGalaxyWLObs with no catalog metadata."""
    return _build_real_wl_obs(tmp_path, "real_wl_obs.gvar", with_meta=False)


@pytest.fixture(name="real_wl_obs_file_with_meta")
def fixture_real_wl_obs_file_with_meta(tmp_path: Path) -> Path:
    """A synthetic NcGalaxyWLObs carrying cluster_ra/dec/z/c metadata."""
    return _build_real_wl_obs(tmp_path, "real_wl_obs_with_meta.gvar", with_meta=True)


@pytest.fixture(name="redshift_dist", params=[e.value for e in GalaxyZGen])
def fixture_redshift_dist(request) -> tuple[str, dict]:
    """Fixture for the redshift distribution."""
    config = {
        "variant": choice(_VARIANTS),
        "zp_min": uniform(0.1, 0.5),
        "zp_max": uniform(0.6, 10.0),
        "sigma0": uniform(0.01, 0.1),
    }
    return request.param, config


@pytest.fixture(
    name="redshift_dist_bad", params=[e.value for e in GalaxyZGen] + ["bogus"]
)
def fixture_redshift_dist_bad(request) -> tuple[str, list[dict]]:
    """Fixture for the redshift distribution bad configuration."""
    config: list[dict] = [
        {"zp_min": 0.5, "zp_max": 0.2, "sigma0": 0.01},
        {"zp_min": "a"},
        {"zp_max": "b"},
        {"sigma0": "c"},
        {"a": None},
        {"zp_min": -0.1},
        {"zp_max": -0.01},
        {"sigma0": -0.01},
        {"variant": "bogus"},
    ]
    return request.param, config


@pytest.fixture(name="shape_factor", params=[s.value for s in ShapeFactorGen])
def fixture_shape_factor(request) -> tuple[str, dict]:
    """Fixture for the galaxy shape-marginalization factor, sweeping the scheme."""
    config = {
        "ellip_conv": choice(["trace", "trace-det"]),
        "ellip_coord": choice(["celestial", "cartesian"]),
        "std_noise": uniform(0.01, 0.1),
        "std_sigma": uniform(0.01, 0.1),
        "c1_sigma": uniform(0.01, 0.1),
        "c2_sigma": uniform(0.01, 0.1),
        "m_sigma": uniform(0.01, 0.1),
    }
    return request.param, config


@pytest.fixture(name="pop_dist")
def fixture_pop_dist() -> tuple[str, dict]:
    """Fixture for the galaxy intrinsic-ellipticity population distribution.

    Every ShapeFactorGen scheme in ``shape_factor`` above is compatible with
    a Gaussian population (see check_shape_pop_compat()), so this fixture
    doesn't need to vary with the scheme the way shape_factor does.
    """
    config = {"sigma": uniform(0.2, 0.4)}
    return GalaxyPopGen.GAUSS.value, config


@pytest.fixture(name="shape_factor_bad", params=[s.value for s in ShapeFactorGen])
def fixture_shape_factor_bad(request) -> tuple[str, list[dict]]:
    """Fixture for the galaxy shape-marginalization factor bad configuration."""
    config: list[dict] = [
        {"ellip_conv": "a"},
        {"ellip_coord": "b"},
        {"std_noise": "d"},
        {"std_sigma": "e"},
        {"c1_sigma": "f"},
        {"c2_sigma": "g"},
        {"m_sigma": "h"},
        {"a": None},
        {"std_noise": -0.01},
        {"std_sigma": -0.01},
        {"c1_sigma": -0.01},
        {"c2_sigma": -0.01},
        {"m_sigma": -0.01},
    ]
    return request.param, config


@pytest.fixture(
    name="halo_profile",
    params=[e.value for e in HaloProfileType],
)
def fixture_halo_profile(request) -> str:
    """Fixture for the halo profile configuration."""
    return request.param


@pytest.fixture(
    name="halo_profile_bad",
    params=[
        "--profile-type=bogus",
        "--profile-type=0",
    ],
)
def fixture_halo_profile_bad(request) -> str:
    """Fixture for the halo profile bad configuration."""
    return request.param


@pytest.fixture(
    name="halo_mass_summary", params=zip(uniform(1e13, 1e16, 3), uniform(1.0, 10.0, 3))
)
def fixture_halo_mass_summary(request) -> tuple[float, float]:
    """Fixture for the halo mass summary."""
    return request.param


@pytest.fixture(
    name="halo_mass_summary_bad",
    params=[
        "--cluster-mass=bogus",
        "--cluster-c=bogus",
        "--cluster-mass=0",
        "--cluster-mass=1e8",
        "--cluster-mass=1e20",
        "--cluster-mass=-5",
        "--cluster-c=0",
        "--cluster-c=-5",
        "--cluster-c=100",
    ],
)
def fixture_halo_mass_summary_bad(request) -> str:
    """Fixture for the halo mass summary bad configuration."""
    return request.param


@pytest.fixture(
    name="halo_position",
    params=zip(uniform(-180, 180, 3), uniform(-90, 90, 3), uniform(0.1, 1.0, 3)),
)
def fixture_halo_position(request) -> tuple[float, float, float]:
    """Fixture for the halo position."""
    return request.param


@pytest.fixture(
    name="halo_position_bad",
    params=[
        ["--cluster-ra=bogus"],
        ["--cluster-dec=bogus"],
        ["--cluster-z=bogus"],
        ["--cluster-ra=192"],
        ["--cluster-dec=92"],
        ["--cluster-ra=-358"],
        ["--cluster-dec=-100"],
        ["--cluster-ra=10", "--ra-min=11", "--ra-max=12"],
        ["--cluster-ra=10", "--ra-min=9", "--ra-max=8"],
        ["--cluster-dec=10", "--dec-min=11", "--dec-max=12"],
        ["--cluster-dec=10", "--dec-min=9", "--dec-max=8"],
        ["--cluster-z=-0.1"],
        ["--cluster-z=5.1"],
        ["--ra-min=bogus"],
        ["--ra-max=bogus"],
        ["--dec-min=bogus"],
        ["--dec-max=bogus"],
    ],
)
def fixture_halo_position_bad(request) -> list[str]:
    """Fixture for the halo position bad configuration."""
    return request.param


@pytest.fixture(name="radius", params=zip(uniform(0.1, 1.0, 3), uniform(1.0, 10.0, 3)))
def fixture_radius(request) -> tuple[float, float]:
    """Fixture for the radius."""
    return request.param


@pytest.fixture(
    name="radius_bad",
    params=[
        "--r-min=bogus",
        "--r-max=bogus",
        "--r-min=0",
        "--r-max=0",
        "--r-min=-0.1",
        "--r-max=-0.1",
    ],
)
def fixture_radius_bad(request) -> list[str]:
    """Fixture for the radius bad configuration."""
    return request.param


@pytest.fixture(
    name="density_bad",
    params=[
        "bogus",
        "0",
        "-5",
    ],
)
def fixture_density_bad(request) -> str:
    """Fixture for the galaxy density bad configuration."""
    return request.param


@pytest.fixture(
    name="fit_parameters",
    params=[
        "NcHaloMassSummary:log10MDelta",
        "NcHaloMassSummary:cDelta",
        "NcHaloPosition:ra",
        "NcHaloPosition:dec",
    ],
)
def fixture_fit_parameters(request) -> str:
    """Fixture for the fit parameters."""
    return request.param


@pytest.fixture(
    name="fit_parameters_bad", params=["NcHaloMassSummary:bogus", "bougs:log10MDelta"]
)
def fixture_fit_parameters_bad(request) -> str:
    """Fixture for the fit parameters bad configuration."""
    return request.param


def test_cluster_wl_app_generate_default(experiment_file):
    """Test the default generation of the cluster WL app."""
    result = runner.invoke(app, ["generate", "cluster-wl", experiment_file.as_posix()])

    assert result.exit_code == 0, result.output

    dataset_file = experiment_file.with_suffix(".dataset.gvar")

    assert (
        experiment_file.exists()
    ), f"Experiment file {experiment_file} does not exist."
    assert dataset_file.exists(), f"Dataset file {dataset_file} does not exist."


def test_cluster_wl_app_summary(experiment_file):
    """Test the summary of the cluster WL app."""
    result = runner.invoke(
        app,
        [
            "generate",
            "cluster-wl",
            experiment_file.as_posix(),
            "--summary",
        ],
    )
    assert result.exit_code == 0, result.output

    assert "Cluster Parameters" in result.output
    assert "Galaxy Sample Parameters" in result.output


def test_cluster_wl_app_generate_redshift(experiment_file, redshift_dist):
    """Test the generation of the cluster WL app with specific distributions."""
    z_dist = redshift_dist[0]
    for key, value in redshift_dist[1].items():
        z_dist += f" {key}={str(value)}"
    result = runner.invoke(
        app,
        [
            "generate",
            "cluster-wl",
            experiment_file.as_posix(),
            f"--z-dist={z_dist}",
            "--summary",
        ],
    )
    assert result.exit_code == 0, result.output

    assert "GalaxyZGenComposed" in result.output

    dataset_file = experiment_file.with_suffix(".dataset.gvar")

    assert (
        experiment_file.exists()
    ), f"Experiment file {experiment_file} does not exist."
    assert (
        dataset_file.exists()
    ), f"Dataset file {dataset_file} does not exist:\n {result.output}"

    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
    dataset = ser.from_binfile(dataset_file.as_posix())
    cluster_data = dataset.get_data(0)
    obs = cluster_data.peek_obs()

    for i in range(obs.len()):
        assert obs.get("zp", i) >= redshift_dist[1]["zp_min"]
        assert obs.get("zp", i) <= redshift_dist[1]["zp_max"]
        assert obs.get("sigma0", i) == redshift_dist[1]["sigma0"]


def test_cluster_wl_app_generate_redshift_bad(experiment_file, redshift_dist_bad):
    """Test the generation of the cluster WL app with bad distribution configuration."""
    for bad_config in redshift_dist_bad[1]:
        z_dist = redshift_dist_bad[0]
        for key, value in bad_config.items():
            z_dist += f" {key}={str(value)}"
        result = runner.invoke(
            app,
            [
                "generate",
                "cluster-wl",
                experiment_file.as_posix(),
                f"--z-dist={z_dist}",
            ],
        )
        assert result.exit_code != 0, result.output

        dataset_file = experiment_file.with_suffix(".dataset.gvar")

        assert (
            not experiment_file.exists()
        ), f"Experiment file {experiment_file} should not exist."
        assert (
            not dataset_file.exists()
        ), f"Dataset file {dataset_file} should not exist."


def test_cluster_wl_app_generate_shape(experiment_file, shape_factor, pop_dist):
    """Test the generation of the cluster WL app with specific shape factors."""
    s_factor = shape_factor[0]
    for key, value in shape_factor[1].items():
        s_factor += f" {key}={str(value)}"
    p_dist = pop_dist[0]
    for key, value in pop_dist[1].items():
        p_dist += f" {key}={str(value)}"
    result = runner.invoke(
        app,
        [
            "generate",
            "cluster-wl",
            experiment_file.as_posix(),
            f"--shape-factor={s_factor}",
            f"--pop-dist={p_dist}",
            "--summary",
        ],
    )
    assert result.exit_code == 0, result.output

    scheme = shape_factor[0]
    assert ShapeFactorGen(scheme).model_cls.__name__ in result.output

    dataset_file = experiment_file.with_suffix(".dataset.gvar")

    assert (
        experiment_file.exists()
    ), f"Experiment file {experiment_file} does not exist."
    assert dataset_file.exists(), f"Dataset file {dataset_file} does not exist."

    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
    dataset = cast(Ncm.Dataset, ser.from_binfile(dataset_file.as_posix()))
    experiment = ser.dict_str_from_yaml_file(experiment_file.as_posix())
    mset = cast(Ncm.MSet, experiment.get("model-set"))
    cluster_data = cast(Nc.DataClusterWLFactor, dataset.get_data(0))
    obs = cluster_data.peek_obs()
    s_dist_model = cast(Nc.GalaxyShapePopGauss, mset.peek_by_name("NcGalaxyShapePop"))

    assert s_dist_model["sigma"] == pop_dist[1]["sigma"]

    bound = 1.0 + 6.0 * (
        pop_dist[1]["sigma"] ** 2 + shape_factor[1]["std_noise"] ** 2
    ) ** (1 / 2)

    for i in range(obs.len()):
        assert abs(obs.get("epsilon_obs_1", i)) <= max(bound, 1.5)
        assert abs(obs.get("epsilon_obs_2", i)) <= max(bound, 1.5)
        assert obs.get("std_noise", i) <= 0.5
        assert abs(obs.get("c1", i)) <= 6.0 * shape_factor[1]["c1_sigma"]
        assert abs(obs.get("c2", i)) <= 6.0 * shape_factor[1]["c2_sigma"]
        assert abs(obs.get("m", i)) <= 1.0 + 6.0 * shape_factor[1]["m_sigma"]

        if scheme == ShapeFactorGen.VAR_ADD.value:
            # A single global sigma, no per-galaxy intrinsic-width draw, so
            # the bound above is exact rather than a generous stand-in.
            assert abs(obs.get("epsilon_obs_1", i)) <= bound
            assert abs(obs.get("epsilon_obs_2", i)) <= bound


def test_cluster_wl_app_generate_shape_bad(experiment_file, shape_factor_bad):
    """Test the generation of the cluster WL.

    Test the generation of the cluster WL app with bad shape factor
    configuration.
    """
    for bad_config in shape_factor_bad[1]:
        s_factor = shape_factor_bad[0]
        for key, value in bad_config.items():
            s_factor += f" {key}={str(value)}"
        result = runner.invoke(
            app,
            [
                "generate",
                "cluster-wl",
                experiment_file.as_posix(),
                f"--shape-factor={s_factor}",
            ],
        )
        assert result.exit_code != 0, result.output

        dataset_file = experiment_file.with_suffix(".dataset.gvar")

        assert (
            not experiment_file.exists()
        ), f"Experiment file {experiment_file} should not exist."
        assert (
            not dataset_file.exists()
        ), f"Dataset file {dataset_file} should not exist."


def test_cluster_wl_app_generate_shape_factor_bogus_type(experiment_file):
    """Test the generation of the cluster WL app with an unknown --shape-factor type."""
    result = runner.invoke(
        app,
        [
            "generate",
            "cluster-wl",
            experiment_file.as_posix(),
            "--shape-factor=bogus",
        ],
    )
    assert result.exit_code != 0, result.output


@pytest.fixture(name="pop_dist_bad")
def fixture_pop_dist_bad() -> tuple[str, list[dict]]:
    """Fixture for the galaxy population distribution bad configuration."""
    config: list[dict] = [
        {"sigma": "c"},
        {"sigma": -0.1},
        {"a": None},
    ]
    return GalaxyPopGen.GAUSS.value, config


def test_cluster_wl_app_generate_pop_bad(experiment_file, pop_dist_bad):
    """Test the generation of the cluster WL app with bad population configuration."""
    for bad_config in pop_dist_bad[1]:
        p_dist = pop_dist_bad[0]
        for key, value in bad_config.items():
            p_dist += f" {key}={str(value)}"
        result = runner.invoke(
            app,
            [
                "generate",
                "cluster-wl",
                experiment_file.as_posix(),
                f"--pop-dist={p_dist}",
            ],
        )
        assert result.exit_code != 0, result.output

        dataset_file = experiment_file.with_suffix(".dataset.gvar")

        assert (
            not experiment_file.exists()
        ), f"Experiment file {experiment_file} should not exist."
        assert (
            not dataset_file.exists()
        ), f"Dataset file {dataset_file} should not exist."


def test_cluster_wl_app_generate_pop_dist_bogus_type(experiment_file):
    """Test the generation of the cluster WL app with an unknown --pop-dist type."""
    result = runner.invoke(
        app,
        [
            "generate",
            "cluster-wl",
            experiment_file.as_posix(),
            "--pop-dist=bogus",
        ],
    )
    assert result.exit_code != 0, result.output


def test_cluster_wl_app_generate_shape_pop_incompatible(experiment_file):
    """Test that an incompatible scheme/pop_dist combination fails cleanly."""
    result = runner.invoke(
        app,
        [
            "generate",
            "cluster-wl",
            experiment_file.as_posix(),
            "--shape-factor=var_add",
            "--pop-dist=beta",
        ],
    )
    assert result.exit_code != 0, result.output
    assert "nc_galaxy_shape_pop_get_sigma" in result.output

    dataset_file = experiment_file.with_suffix(".dataset.gvar")

    assert (
        not experiment_file.exists()
    ), f"Experiment file {experiment_file} should not exist."
    assert not dataset_file.exists(), f"Dataset file {dataset_file} should not exist."


def test_cluster_model_position_data():
    """ClusterModel.position_data reports back the position it was built with.

    A plain accessor (reads the internal NcHaloPosition back out as a
    HaloPositionData), not otherwise exercised by any CLI-level test.
    """
    position = HaloPositionData(ra=12.34, dec=-55.123, z=0.2)
    model = ClusterModel(position=position)

    reported = model.position_data

    assert reported.ra == pytest.approx(position.ra)
    assert reported.dec == pytest.approx(position.dec)
    assert reported.z == pytest.approx(position.z)


def test_galaxy_pop_gen_beta_get_shape_pop():
    """GalaxyPopGenBeta.get_shape_pop() returns the live NcGalaxyShapePopBeta
    it built, with alpha/beta already set -- e.g. for get_mean()/
    get_concentration() reporting, not otherwise exercised by any CLI-level
    test (those only inspect the model after a yaml round-trip)."""
    pop_gen = GalaxyPopGenBeta(alpha=2.0, beta=5.0)
    pop = pop_gen.get_shape_pop()

    assert isinstance(pop, Nc.GalaxyShapePopBeta)
    assert pop["alpha"] == 2.0
    assert pop["beta"] == 5.0


def test_cluster_wl_load_app_missing_ra_dec_columns(experiment_file, tmp_path: Path):
    """load_cluster_wl()'s own guard: a catalog missing 'ra'/'dec' columns
    must fail cleanly, not proceed into a footprint computation that has
    nothing to compute from."""
    obs = Nc.GalaxyWLObs.new(
        Nc.GalaxyWLObsEllipConv.TRACE,
        Nc.WLEllipticityFrame.CELESTIAL,
        1,
        ["epsilon_obs_1", "epsilon_obs_2", "std_noise", "z"],
    )
    obs.set("epsilon_obs_1", 0, 0.0)
    obs.set("epsilon_obs_2", 0, 0.0)
    obs.set("std_noise", 0, 0.1)
    obs.set("z", 0, 0.8)

    no_ra_dec_file = tmp_path / "no_ra_dec.gvar"
    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
    ser.to_binfile(obs, no_ra_dec_file.as_posix())

    result = runner.invoke(
        app,
        [
            "generate",
            "cluster-wl-load",
            experiment_file.as_posix(),
            f"--data-file={no_ra_dec_file.as_posix()}",
            "--cluster-ra=0.0",
            "--cluster-dec=0.0",
            "--cluster-z=0.8",
            "--cluster-c=4.0",
            "--shape-factor=var_add ellip_conv=trace ellip_coord=celestial",
        ],
    )
    assert result.exit_code != 0, result.output
    assert "'ra' and 'dec'" in str(result.exception)


def test_cluster_wl_app_generate_pop_gauss_local(experiment_file):
    """Test the generation of the cluster WL app with --pop-dist=gauss_local.

    GaussLocal reads a per-galaxy e_rms catalog column instead of a global
    model parameter -- the only ``GalaxyPopGen`` branch not otherwise
    exercised by ``test_cluster_wl_app_generate_shape`` (which only sweeps
    Gauss). Covers ``GalaxyPopGenGaussLocal``'s constructor/model_post_init/
    register_models/get_mfuncs/write_calib.
    """
    result = runner.invoke(
        app,
        [
            "generate",
            "cluster-wl",
            experiment_file.as_posix(),
            "--pop-dist=gauss_local sigma=0.3 std_sigma=0.05",
        ],
    )
    assert result.exit_code == 0, result.output

    dataset_file = experiment_file.with_suffix(".dataset.gvar")
    assert dataset_file.exists(), f"Dataset file {dataset_file} does not exist."

    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
    dataset = cast(Ncm.Dataset, ser.from_binfile(dataset_file.as_posix()))
    cluster_data = cast(Nc.DataClusterWLFactor, dataset.get_data(0))
    obs = cluster_data.peek_obs()

    for i in range(obs.len()):
        e_rms = obs.get(Nc.GALAXY_SHAPE_POP_GAUSS_LOCAL_COL_E_RMS, i)
        assert 0.0 < e_rms <= 1.0


def test_cluster_wl_app_generate_pop_beta(experiment_file):
    """Test the generation of the cluster WL app with --pop-dist=beta.

    Beta requires a scheme that doesn't linearize around a Gaussian (see
    ``check_shape_pop_compat()``) -- fixed_quad qualifies. Covers
    ``GalaxyPopGenBeta``'s get_shape_pop/register_models/get_mfuncs, and
    ``mfunc_oa.add(func)`` in generate.py's own ``_build_experiment`` (the
    only ``GalaxyPopGen`` variant with non-empty ``get_mfuncs()``), plus
    the ``NcGalaxyShapePopBeta:mean``/``:std`` ``NcmMSetFuncList`` entries
    read back from the written ``.functions.yaml``.
    """
    result = runner.invoke(
        app,
        [
            "generate",
            "cluster-wl",
            experiment_file.as_posix(),
            "--shape-factor=fixed_quad",
            "--pop-dist=beta alpha=2.0 beta=5.0",
        ],
    )
    assert result.exit_code == 0, result.output

    dataset_file = experiment_file.with_suffix(".dataset.gvar")
    functions_file = experiment_file.with_suffix(".functions.yaml")
    assert dataset_file.exists(), f"Dataset file {dataset_file} does not exist."
    assert functions_file.exists(), f"Functions file {functions_file} does not exist."

    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
    # Dataset must be read first with the same NcmSerialize instance --
    # shared object aliases span both files.
    ser.from_binfile(dataset_file.as_posix())
    experiment = ser.dict_str_from_yaml_file(experiment_file.as_posix())
    mset = cast(Ncm.MSet, experiment.get("model-set"))
    pop_model = cast(Nc.GalaxyShapePopBeta, mset.peek_by_name("NcGalaxyShapePop"))

    assert pop_model["alpha"] == 2.0
    assert pop_model["beta"] == 5.0

    mfunc_oa = cast(Ncm.ObjArray, ser.array_from_yaml_file(functions_file.as_posix()))
    funcs = [cast(Ncm.MSetFuncList, mfunc_oa.get(i)) for i in range(mfunc_oa.len())]
    names = {(f.peek_ns(), f.peek_name()) for f in funcs}
    assert names == {
        ("NcGalaxyShapePopBeta", "mean"),
        ("NcGalaxyShapePopBeta", "std"),
    }
    for func in funcs:
        func.eval0(mset)


def test_cluster_wl_app_generate_redshift_dist_bogus_type(experiment_file):
    """Test the generation of the cluster WL app with an unknown --z-dist type."""
    result = runner.invoke(
        app,
        [
            "generate",
            "cluster-wl",
            experiment_file.as_posix(),
            "--z-dist=bogus",
        ],
    )
    assert result.exit_code != 0, result.output


def test_cluster_wl_app_generate_halo_profile(experiment_file: Path, halo_profile: str):
    """Test the generation of the cluster WL app with diferent cluster profiles."""
    # pylint: disable=unused-variable
    result = runner.invoke(
        app,
        [
            "generate",
            "cluster-wl",
            experiment_file.as_posix(),
            f"--profile-type={halo_profile}",
        ],
    )
    assert result.exit_code == 0, result.output

    dataset_file = experiment_file.with_suffix(".dataset.gvar")

    assert (
        experiment_file.exists()
    ), f"Experiment file {experiment_file} does not exist."
    assert dataset_file.exists(), f"Dataset file {dataset_file} does not exist."

    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
    _ = ser.from_binfile(dataset_file.as_posix())
    experiment = ser.dict_str_from_yaml_file(experiment_file.as_posix())
    mset: Ncm.MSet = cast(Ncm.MSet, experiment.get("model-set"))

    hp = mset.peek_by_name("NcHaloDensityProfile")

    match halo_profile:
        case HaloProfileType.NFW.value:
            assert isinstance(hp, Nc.HaloDensityProfileNFW)
        case HaloProfileType.EINASTO.value:
            assert isinstance(hp, Nc.HaloDensityProfileEinasto)
        case HaloProfileType.HERNQUIST.value:
            assert isinstance(hp, Nc.HaloDensityProfileHernquist)


def test_cluster_wl_app_generate_halo_profile_bad(experiment_file, halo_profile_bad):
    """Test the generation of the cluster WL app.

    Test the generation of the cluster WL app with bad cluster profile configuration.
    """
    result = runner.invoke(
        app,
        [
            "generate",
            "cluster-wl",
            experiment_file.as_posix(),
            halo_profile_bad,
        ],
    )
    assert result.exit_code != 0, result.output

    dataset_file = experiment_file.with_suffix(".dataset.gvar")

    assert (
        not experiment_file.exists()
    ), f"Experiment file {experiment_file} should not exist."
    assert not dataset_file.exists(), f"Dataset file {dataset_file} should not exist."


def test_cluster_wl_app_halo_mass_summary(experiment_file, halo_mass_summary):
    """Test the generation of the cluster WL app with specific halo mass summary."""
    # pylint: disable=unused-variable
    mass, c = halo_mass_summary
    result = runner.invoke(
        app,
        [
            "generate",
            "cluster-wl",
            experiment_file.as_posix(),
            f"--cluster-mass={mass}",
            f"--cluster-c={c}",
        ],
    )
    assert result.exit_code == 0, result.output

    dataset_file = experiment_file.with_suffix(".dataset.gvar")

    assert (
        experiment_file.exists()
    ), f"Experiment file {experiment_file} does not exist."
    assert dataset_file.exists(), f"Dataset file {dataset_file} does not exist."

    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
    _ = ser.from_binfile(dataset_file.as_posix())
    experiment = ser.dict_str_from_yaml_file(experiment_file.as_posix())
    mset = cast(Ncm.MSet, experiment.get("model-set"))
    hms = cast(Nc.HaloMassSummary, mset.peek_by_name("NcHaloMassSummary"))

    assert hms["log10MDelta"] == log10(mass)
    assert hms["cDelta"] == c


def test_cluster_wl_app_halo_mass_summary_bad(experiment_file, halo_mass_summary_bad):
    """Test the generation of the cluster WL app.

    Test the generation of the cluster WL app with bad halo mass summary
    configuration.
    """
    result = runner.invoke(
        app,
        [
            "generate",
            "cluster-wl",
            experiment_file.as_posix(),
            halo_mass_summary_bad,
        ],
    )
    assert result.exit_code != 0, result.output

    dataset_file = experiment_file.with_suffix(".dataset.gvar")

    assert (
        not experiment_file.exists()
    ), f"Experiment file {experiment_file} should not exist."
    assert not dataset_file.exists(), f"Dataset file {dataset_file} should not exist."


def test_cluster_wl_app_halo_position(experiment_file, halo_position):
    """Test the generation of the cluster WL app with specific RA and Dec."""
    # pylint: disable=unused-variable
    cluster_ra, cluster_dec, cluster_z = halo_position
    dec_min = cluster_dec - 0.2
    dec_max = cluster_dec + 0.2
    ra_min = cluster_ra - 0.2 / cos(deg2rad(cluster_dec))
    ra_max = cluster_ra + 0.2 / cos(deg2rad(cluster_dec))
    result = runner.invoke(
        app,
        [
            "generate",
            "cluster-wl",
            experiment_file.as_posix(),
            f"--cluster-ra={cluster_ra}",
            f"--cluster-dec={cluster_dec}",
            f"--cluster-z={cluster_z}",
            f"--ra-min={ra_min}",
            f"--ra-max={ra_max}",
            f"--dec-min={dec_min}",
            f"--dec-max={dec_max}",
        ],
    )
    assert result.exit_code == 0, result.output

    dataset_file = experiment_file.with_suffix(".dataset.gvar")

    assert (
        experiment_file.exists()
    ), f"Experiment file {experiment_file} does not exist."
    assert dataset_file.exists(), f"Dataset file {dataset_file} does not exist."

    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
    _ = ser.from_binfile(dataset_file.as_posix())
    experiment = ser.dict_str_from_yaml_file(experiment_file.as_posix())
    mset = cast(Ncm.MSet, experiment.get("model-set"))
    hp = cast(Nc.HaloPosition, mset.peek_by_name("NcHaloPosition"))

    assert hp["ra"] == cluster_ra
    assert hp["dec"] == cluster_dec
    assert hp["z"] == cluster_z


def test_cluster_wl_app_halo_position_bad(experiment_file, halo_position_bad):
    """Test the generation of the cluster WL app with specific RA and Dec."""
    result = runner.invoke(
        app,
        [
            "generate",
            "cluster-wl",
            experiment_file.as_posix(),
            *halo_position_bad,
        ],
    )
    assert result.exit_code != 0, result.output

    dataset_file = experiment_file.with_suffix(".dataset.gvar")

    assert not experiment_file.exists(), f"Experiment file {experiment_file} exists."
    assert not dataset_file.exists(), f"Dataset file {dataset_file} exists."


def test_cluster_wl_app_radius(experiment_file, radius):
    """Test the generation of the cluster WL app with specific radius."""
    r_min, r_max = radius
    result = runner.invoke(
        app,
        [
            "generate",
            "cluster-wl",
            experiment_file.as_posix(),
            f"--r-min={r_min}",
            f"--r-max={r_max}",
        ],
    )
    assert result.exit_code == 0, result.output

    dataset_file = experiment_file.with_suffix(".dataset.gvar")

    assert (
        experiment_file.exists()
    ), f"Experiment file {experiment_file} does not exist."
    assert dataset_file.exists(), f"Dataset file {dataset_file} does not exist."

    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
    dataset = cast(Ncm.Dataset, ser.from_binfile(dataset_file.as_posix()))
    experiment = ser.dict_str_from_yaml_file(experiment_file.as_posix())
    mset = cast(Ncm.MSet, experiment.get("model-set"))
    cluster_data = cast(Nc.DataClusterWLFactor, dataset.get_data(0))
    obs = cluster_data.peek_obs()
    hp = cast(Nc.HaloPosition, mset.peek_by_name("NcHaloPosition"))
    cosmo = cast(Nc.HICosmo, mset.peek_by_name("NcHICosmo"))

    hp.prepare(cosmo)

    for i in range(obs.len()):
        ra = obs.get("ra", i)
        dec = obs.get("dec", i)
        r = hp.projected_radius_from_ra_dec(cosmo, ra, dec)

        assert r >= r_min
        assert r <= r_max


def test_cluster_wl_app_radius_bad(experiment_file, radius_bad):
    """Test the generation of the cluster WL app with bad radius configuration."""
    result = runner.invoke(
        app,
        [
            "generate",
            "cluster-wl",
            experiment_file.as_posix(),
            radius_bad,
        ],
    )
    assert result.exit_code != 0, result.output

    dataset_file = experiment_file.with_suffix(".dataset.gvar")

    assert (
        not experiment_file.exists()
    ), f"Experiment file {experiment_file} should not exist."
    assert not dataset_file.exists(), f"Dataset file {dataset_file} should not exist."


def test_cluster_wl_app_galaxy_density(experiment_file, halo_position):
    """Test the generation of the cluster WL app with specific galaxy density."""
    galaxy_density = uniform(10, 50)
    cluster_ra, cluster_dec, cluster_z = halo_position
    dec_min = cluster_dec - 0.2
    dec_max = cluster_dec + 0.2
    ra_min = cluster_ra - 0.2 / cos(deg2rad(cluster_dec))
    ra_max = cluster_ra + 0.2 / cos(deg2rad(cluster_dec))
    n_galaxies = int(
        (deg2rad(ra_max) - deg2rad(ra_min))
        * (sin(deg2rad(dec_max)) - sin(deg2rad(dec_min)))
        * (180.0 / pi) ** 2.0
        * galaxy_density
        * 60
        * 60
    )
    result = runner.invoke(
        app,
        [
            "generate",
            "cluster-wl",
            experiment_file.as_posix(),
            f"--galaxy-density={galaxy_density}",
            f"--cluster-ra={cluster_ra}",
            f"--cluster-dec={cluster_dec}",
            f"--cluster-z={cluster_z}",
            f"--ra-min={ra_min}",
            f"--ra-max={ra_max}",
            f"--dec-min={dec_min}",
            f"--dec-max={dec_max}",
        ],
    )
    assert result.exit_code == 0, result.output

    dataset_file = experiment_file.with_suffix(".dataset.gvar")

    assert (
        experiment_file.exists()
    ), f"Experiment file {experiment_file} does not exist."
    assert dataset_file.exists(), f"Dataset file {dataset_file} does not exist."

    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
    dataset = ser.from_binfile(dataset_file.as_posix())
    cluster_data = dataset.get_data(0)
    obs = cluster_data.peek_obs()

    # resample() retries each row's position draw until it lands in the
    # cut -- every requested galaxy is kept, so this is exact.
    assert obs.len() == n_galaxies


def test_cluster_wl_app_galaxy_density_bad(experiment_file, density_bad):
    """Test the generation of the cluster WL app.

    Test the generation of the cluster WL app with bad galaxy density configuration.
    """
    result = runner.invoke(
        app,
        [
            "generate",
            "cluster-wl",
            experiment_file.as_posix(),
            f"--galaxy-density={density_bad}",
        ],
    )
    assert result.exit_code != 0, result.output

    dataset_file = experiment_file.with_suffix(".dataset.gvar")

    assert (
        not experiment_file.exists()
    ), f"Experiment file {experiment_file} should not exist."
    assert not dataset_file.exists(), f"Dataset file {dataset_file} should not exist."


def test_cluster_wl_app_seed(experiment_file, experiment_file_copy):
    """Test the generation of the cluster WL app with a specific seed."""
    result = runner.invoke(
        app,
        [
            "generate",
            "cluster-wl",
            experiment_file.as_posix(),
            "--seed=42",
        ],
    )
    assert result.exit_code == 0, result.output

    dataset_file = experiment_file.with_suffix(".dataset.gvar")

    assert (
        experiment_file.exists()
    ), f"Experiment file {experiment_file} does not exist."
    assert dataset_file.exists(), f"Dataset file {dataset_file} does not exist."

    result = runner.invoke(
        app,
        [
            "generate",
            "cluster-wl",
            experiment_file_copy.as_posix(),
            "--seed=42",
        ],
    )
    assert result.exit_code == 0, result.output

    dataset_file_copy = experiment_file_copy.with_suffix(".dataset.gvar")

    assert (
        experiment_file_copy.exists()
    ), f"Experiment file {experiment_file_copy} does not exist."
    assert (
        dataset_file_copy.exists()
    ), f"Dataset file {dataset_file_copy} does not exist."
    assert cmp(
        dataset_file.as_posix(), dataset_file_copy.as_posix(), shallow=False
    ), f"Dataset files {dataset_file} and {dataset_file_copy} are not equal."


def test_cluster_wl_app_file_extension_bad(experiment_file):
    """Test the generation of the cluster WL app with an invalid file extension."""
    experiment_file = experiment_file.with_suffix(".bogus")
    result = runner.invoke(
        app,
        [
            "generate",
            "cluster-wl",
            experiment_file.as_posix(),
        ],
    )
    assert result.exit_code != 0, result.output

    dataset_file = experiment_file.with_suffix(".dataset.gvar")

    assert (
        not experiment_file.exists()
    ), f"Experiment file {experiment_file} does not exist."
    assert not dataset_file.exists(), f"Dataset file {dataset_file} does not exist."


def test_cluster_wl_app_fit_parameters(experiment_file, fit_parameters):
    """Test the generation of the cluster WL app with specific fit parameters."""
    # pylint: disable=unused-variable
    result = runner.invoke(
        app,
        [
            "generate",
            "cluster-wl",
            experiment_file.as_posix(),
            f"--parameter-list={fit_parameters}",
        ],
    )
    assert result.exit_code == 0, result.output

    dataset_file = experiment_file.with_suffix(".dataset.gvar")

    assert (
        experiment_file.exists()
    ), f"Experiment file {experiment_file} does not exist."
    assert dataset_file.exists(), f"Dataset file {dataset_file} does not exist."

    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
    _ = ser.from_binfile(dataset_file.as_posix())
    experiment = ser.dict_str_from_yaml_file(experiment_file.as_posix())
    mset = cast(Ncm.MSet, experiment.get("model-set"))
    fparam = mset.fparam_full_name(0)

    assert (
        fparam == fit_parameters
    ), f"Fit parameter {fparam} does not match {fit_parameters}."


def test_cluster_wl_app_fit_parameters_bad(experiment_file, fit_parameters_bad):
    """Test the generation of the cluster WL app with bad fit parameters."""
    result = runner.invoke(
        app,
        [
            "generate",
            "cluster-wl",
            experiment_file.as_posix(),
            f"--parameter-list={fit_parameters_bad}",
        ],
    )
    assert result.exit_code != 0, result.output

    dataset_file = experiment_file.with_suffix(".dataset.gvar")

    assert (
        not experiment_file.exists()
    ), f"Experiment file {experiment_file} should not exist."
    assert not dataset_file.exists(), f"Dataset file {dataset_file} should not exist."


def test_cluster_wl_load_app_data_file(experiment_file, real_wl_obs_file_with_meta):
    """Test loading a real NcGalaxyWLObs catalog via --data-file.

    cluster_ra/dec/z/c are all taken from the catalog's own metadata --
    none are given on the command line.
    """
    result = runner.invoke(
        app,
        [
            "generate",
            "cluster-wl-load",
            experiment_file.as_posix(),
            f"--data-file={real_wl_obs_file_with_meta.as_posix()}",
            "--shape-factor=var_add ellip_conv=trace ellip_coord=celestial",
            "--summary",
        ],
    )
    assert result.exit_code == 0, result.output
    assert "Number of galaxies" in result.output

    dataset_file = experiment_file.with_suffix(".dataset.gvar")

    assert (
        experiment_file.exists()
    ), f"Experiment file {experiment_file} does not exist."
    assert dataset_file.exists(), f"Dataset file {dataset_file} does not exist."

    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
    dataset = cast(Ncm.Dataset, ser.from_binfile(dataset_file.as_posix()))
    cluster_data = cast(Nc.DataClusterWLFactor, dataset.get_data(0))

    assert cluster_data.peek_obs().len() == 20


def test_cluster_wl_load_app_no_metadata_requires_explicit_values(
    experiment_file, real_wl_obs_file
):
    """Without --cluster-ra/dec/z/c and without catalog metadata, this must
    fail with a clear error rather than silently using a meaningless value.
    """
    result = runner.invoke(
        app,
        [
            "generate",
            "cluster-wl-load",
            experiment_file.as_posix(),
            f"--data-file={real_wl_obs_file.as_posix()}",
            "--shape-factor=var_add ellip_conv=trace ellip_coord=celestial",
        ],
    )
    assert result.exit_code != 0, result.output


def test_cluster_wl_load_app_cli_override_beats_metadata(
    experiment_file, real_wl_obs_file_with_meta
):
    """An explicit --cluster-ra takes precedence over catalog metadata, and
    is still validated against the catalog's RA/Dec window (there is no
    lensing signal outside it to fit) even though metadata is present.
    """
    result = runner.invoke(
        app,
        [
            "generate",
            "cluster-wl-load",
            experiment_file.as_posix(),
            f"--data-file={real_wl_obs_file_with_meta.as_posix()}",
            "--cluster-ra=999",
            "--shape-factor=var_add ellip_conv=trace ellip_coord=celestial",
        ],
    )
    assert result.exit_code != 0, result.output


def test_cluster_wl_load_app_cluster_position_inside_window(
    experiment_file, real_wl_obs_file
):
    """An explicit cluster position inside the catalog's RA/Dec window is
    accepted (catalog has no metadata, so all four values must be given).
    """
    result = runner.invoke(
        app,
        [
            "generate",
            "cluster-wl-load",
            experiment_file.as_posix(),
            f"--data-file={real_wl_obs_file.as_posix()}",
            "--cluster-ra=0.0",
            "--cluster-dec=0.0",
            "--cluster-z=0.8",
            "--cluster-c=4.0",
            "--shape-factor=var_add ellip_conv=trace ellip_coord=celestial",
        ],
    )
    assert result.exit_code == 0, result.output


def test_cluster_wl_load_app_catalog(experiment_file):
    """Test loading a real catalog via --catalog (downloaded and cached from
    the NumCosmo datafile-release-v1.0.0 GitHub release on first use).

    Uses the smallest curated Subaru HSC-SSP PDR1 field (HWL16a-094, ~3MB,
    2200 galaxies) to keep the download light; requires network access, same
    as nc.data.test_snia_cov.test_constructor_catalog_id's own real-download
    convention. Covers LoadClusterWL._load_obs's --catalog branch and, in the
    C layer, nc_galaxy_wl_obs_new_from_catalog_id()/
    nc_galaxy_wl_obs_catalog_id_get_filename() -- otherwise untested (0%
    coverage), since every other cluster-wl-load test in this file uses
    --data-file with a small synthetic catalog instead.
    """
    result = runner.invoke(
        app,
        [
            "generate",
            "cluster-wl-load",
            experiment_file.as_posix(),
            "--catalog=094",
            "--shape-factor=var_add ellip_conv=trace ellip_coord=celestial",
            "--summary",
        ],
    )
    assert result.exit_code == 0, result.output
    assert "Number of galaxies" in result.output

    dataset_file = experiment_file.with_suffix(".dataset.gvar")
    assert dataset_file.exists(), f"Dataset file {dataset_file} does not exist."

    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
    dataset = cast(Ncm.Dataset, ser.from_binfile(dataset_file.as_posix()))
    cluster_data = cast(Nc.DataClusterWLFactor, dataset.get_data(0))

    assert cluster_data.peek_obs().len() == 2200


def test_cluster_wl_load_app_catalog_data_file_mutually_exclusive(
    experiment_file, real_wl_obs_file
):
    """--catalog and --data-file cannot be given together."""
    result = runner.invoke(
        app,
        [
            "generate",
            "cluster-wl-load",
            experiment_file.as_posix(),
            f"--data-file={real_wl_obs_file.as_posix()}",
            "--catalog=002",
        ],
    )
    assert result.exit_code != 0, result.output


def test_cluster_wl_load_app_catalog_data_file_neither_given(experiment_file):
    """Exactly one of --catalog or --data-file is required."""
    result = runner.invoke(
        app,
        [
            "generate",
            "cluster-wl-load",
            experiment_file.as_posix(),
        ],
    )
    assert result.exit_code != 0, result.output


def test_cluster_wl_load_app_data_file_wrong_type(experiment_file, tmp_path: Path):
    """--data-file pointing at a validly-serialized object of the wrong type
    must fail with a clear error (LoadClusterWL._load_obs's isinstance
    check), not silently misinterpret it as a NcGalaxyWLObs."""
    wrong_type_file = tmp_path / "not_a_wl_obs.gvar"
    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
    ser.to_binfile(Nc.HICosmoDEXcdm.new(), wrong_type_file.as_posix())

    result = runner.invoke(
        app,
        [
            "generate",
            "cluster-wl-load",
            experiment_file.as_posix(),
            f"--data-file={wrong_type_file.as_posix()}",
            "--shape-factor=var_add ellip_conv=trace ellip_coord=celestial",
        ],
    )
    assert result.exit_code != 0, result.output
    assert "NcGalaxyWLObs" in str(result.exception)


def test_cluster_wl_load_app_ellip_conv_mismatch(experiment_file, real_wl_obs_file):
    """--shape-factor ellip_conv/ellip_coord must match the catalog's own."""
    result = runner.invoke(
        app,
        [
            "generate",
            "cluster-wl-load",
            experiment_file.as_posix(),
            f"--data-file={real_wl_obs_file.as_posix()}",
            "--shape-factor=var_add ellip_conv=trace-det ellip_coord=celestial",
        ],
    )
    assert result.exit_code != 0, result.output

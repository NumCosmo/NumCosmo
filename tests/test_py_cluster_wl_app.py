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

"""Test NumCosmo app's Cluster WL module."""

from filecmp import cmp
import pytest
from typer.testing import CliRunner
from numpy import sin, cos, pi, log10, deg2rad
from numpy.random import uniform, choice

from numcosmo_py import (
    Ncm,
    Nc,
)
from numcosmo_py.app import app
from numcosmo_py.experiments.cluster_wl import (
    GalaxyZGen,
    GalaxyShapeGen,
    HaloProfileType,
)


runner = CliRunner()

Ncm.cfg_init()


@pytest.fixture(name="experiment_file")
def fixture_experiment_file(tmp_path) -> str:
    """Fixture for the experiment file."""
    tmp_file = tmp_path / "experiment.yaml"
    return tmp_file


@pytest.fixture(name="experiment_file_copy")
def fixture_experiment_file_copy(tmp_path) -> str:
    """Fixture for the experiment file copy."""
    tmp_file = tmp_path / "experiment_copy.yaml"
    return tmp_file


@pytest.fixture(name="redshift_dist", params=[e.value for e in GalaxyZGen])
def fixture_redshift_dist(request) -> str:
    """Fixture for the redshift distribution."""
    if request.param == GalaxyZGen.SPEC:
        config = {"z_min": uniform(0.1, 0.5), "z_max": uniform(0.6, 10.0)}
        bad_config = [
            {"zp_min": 0.1, "zp_max": 0.5, "sigma0": 0.01},
            {"z_min": 0.1, "z_max": 0.5, "sigma0": 0.01},
            {"zp_min": 0.1, "zp_max": 0.5},
            {"z_min": 0.5, "z_max": 0.2},
            {"z_min": "a"},
            {"z_max": "b"},
            {"a": None},
            {"z_min": -0.1},
            {"z_max": -0.01},
        ]
    else:
        config = {
            "zp_min": uniform(0.1, 0.5),
            "zp_max": uniform(0.6, 10.0),
            "sigma0": uniform(0.01, 0.1),
        }
        bad_config = [
            {"z_min": 0.1, "z_max": 0.5, "sigma0": 0.01},
            {"zp_min": 0.5, "zp_max": 0.2},
            {"zp_min": "a"},
            {"zp_max": "b"},
            {"sigma0": "c"},
            {"a": None},
            {"zp_min": -0.1},
            {"zp_max": -0.01},
            {"sigma0": -0.01},
        ]
    return request.param, config, bad_config


@pytest.fixture(name="shape_dist", params=[e.value for e in GalaxyShapeGen])
def fixture_shape_dist(request) -> str:
    """Fixture for the galaxy shape distribution."""
    if request.param == GalaxyShapeGen.GAUSS:
        config = {
            "ellip_conv": choice(["trace", "trace-det"]),
            "ellip_coord": choice(["celestial", "euclidean"]),
            "sigma": uniform(0.2, 0.5),
            "std_noise": uniform(0.01, 0.1),
        }
        bad_config = [
            {
                "ellip_conv": "trace",
                "ellip_coord": "celestial",
                "sigma": 0.15,
                "std_noise": 0.01,
                "std_shape": 0.2,
            },
            {"std_shape": 0.2},
            {"ellip_conv": "a"},
            {"ellip_coord": "b"},
            {"sigma": "c"},
            {"std_noise": "d"},
            {"a": None},
            {"sigma": -0.1},
            {"std_noise": -0.01},
        ]
    else:
        config = {
            "ellip_conv": choice(["trace", "trace-det"]),
            "ellip_coord": choice(["celestial", "euclidean"]),
            "std_shape": uniform(0.2, 0.4),
            "std_noise": uniform(0.01, 0.1),
            "std_sigma": uniform(0.01, 0.1),
            "c1_sigma": uniform(0.01, 0.1),
            "c2_sigma": uniform(0.01, 0.1),
            "m_sigma": uniform(0.01, 0.1),
        }
        bad_config = [
            {
                "ellip_conv": "trace",
                "ellip_coord": "celestial",
                "std_shape": 0.15,
                "std_noise": 0.01,
                "std_sigma": 0.01,
                "c1_sigma": 0.01,
                "c2_sigma": 0.01,
                "m_sigma": 0.01,
                "sigma": 0.15,
            },
            {"sigma": 0.15},
            {"ellip_conv": "a"},
            {"ellip_coord": "b"},
            {"std_shape": "c"},
            {"std_noise": "d"},
            {"std_sigma": "e"},
            {"c1_sigma": "f"},
            {"c2_sigma": "g"},
            {"m_sigma": "h"},
            {"a": None},
            {"std_shape": -0.1},
            {"std_noise": -0.01},
            {"std_sigma": -0.01},
            {"c1_sigma": -0.01},
            {"c2_sigma": -0.01},
            {"m_sigma": -0.01},
        ]
    return request.param, config, bad_config


@pytest.fixture(
    name="halo_profile",
    params=[e.value for e in HaloProfileType],
)
def fixture_halo_profile(request) -> str:
    """Fixture for the halo profile configuration"""
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
    params=zip(uniform(0, 360, 3), uniform(-90, 90, 3), uniform(0.1, 1.0, 3)),
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
        ["--cluster-ra=372"],
        ["--cluster-dec=92"],
        ["--cluster-ra=-10"],
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
        "--galaxy-density=bogus",
        "--galaxy-density=0",
        "--galaxy-density=-5",
    ],
)
def fixture_density_bad(request) -> str:
    """Fixture for the galaxy density bad configuration."""
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
            "--z-dist",
            z_dist,
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

    for i in range(obs.len()):
        if redshift_dist[0] == GalaxyZGen.SPEC:
            assert obs.get("z", i) >= redshift_dist[1]["z_min"]
            assert obs.get("z", i) <= redshift_dist[1]["z_max"]
        else:
            assert obs.get("zp", i) >= redshift_dist[1]["zp_min"]
            assert obs.get("zp", i) <= redshift_dist[1]["zp_max"]
            assert obs.get("sigma_0", i) == redshift_dist[1]["sigma0"]
            assert obs.get("sigma_z", i) == redshift_dist[1]["sigma0"] * (
                1 + obs.get("z", i)
            )


def test_cluster_wl_app_generate_redshift_bad(experiment_file, redshift_dist):
    """Test the generation of the cluster WL app with bad distribution configuration."""
    for bad_config in redshift_dist[2]:
        z_dist = redshift_dist[0]
        for key, value in bad_config.items():
            z_dist += f" {key}={str(value)}"
        result = runner.invoke(
            app,
            [
                "generate",
                "cluster-wl",
                experiment_file.as_posix(),
                "--z-dist",
                z_dist,
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


def test_cluster_wl_app_generate_shape(experiment_file, shape_dist):
    """Test the generation of the cluster WL app with specific shape distributions."""
    s_dist = shape_dist[0]
    for key, value in shape_dist[1].items():
        s_dist += f" {key}={str(value)}"
    result = runner.invoke(
        app,
        [
            "generate",
            "cluster-wl",
            experiment_file.as_posix(),
            "--shape-dist",
            s_dist,
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

    for i in range(obs.len()):
        if shape_dist[0] == GalaxyShapeGen.GAUSS:
            assert (
                abs(obs.get("epsilon_obs_1", i))
                <= 1.0 + 6.0 * shape_dist[1]["std_noise"]
            )
            assert (
                abs(obs.get("epsilon_obs_2", i))
                <= 1.0 + 6.0 * shape_dist[1]["std_noise"]
            )
            assert obs.get("std_noise", i) == shape_dist[1]["std_noise"]
        else:
            assert abs(obs.get("epsilon_obs_1", i)) <= 1.0 + 6.0 * (
                shape_dist[1]["std_shape"] ** 2 + shape_dist[1]["std_noise"] ** 2
            ) ** (1 / 2)
            assert abs(obs.get("epsilon_obs_2", i)) <= 1.0 + 6.0 * (
                shape_dist[1]["std_shape"] ** 2 + shape_dist[1]["std_noise"] ** 2
            ) ** (1 / 2)
            assert obs.get("std_shape", i) <= 0.45
            assert abs(obs.get("c1", i)) <= 6.0 * shape_dist[1]["c1_sigma"]
            assert abs(obs.get("c2", i)) <= 6.0 * shape_dist[1]["c2_sigma"]
            assert abs(obs.get("m", i)) <= 1.0 + 6.0 * shape_dist[1]["m_sigma"]


def test_cluster_wl_app_generate_shape_bad(experiment_file, shape_dist):
    """Test the generation of the cluster WL app with
    bad shape distribution configuration."""
    for bad_config in shape_dist[2]:
        s_dist = shape_dist[0]
        for key, value in bad_config.items():
            s_dist += f" {key}={str(value)}"
        result = runner.invoke(
            app,
            [
                "generate",
                "cluster-wl",
                experiment_file.as_posix(),
                "--shape-dist",
                s_dist,
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


def test_cluster_wl_app_generate_halo_profile(experiment_file, halo_profile):
    """Test the generation of the cluster WL app with diferent cluster profiles."""
    # pylint: disable=unused-variable
    result = runner.invoke(
        app,
        [
            "generate",
            "cluster-wl",
            experiment_file.as_posix(),
            "--profile-type",
            halo_profile,
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
    experiment = ser.dict_str_from_yaml_file(experiment_file.as_posix())
    mset = experiment.get("model-set")

    hp = mset.peek_by_name("NcHaloDensityProfile")

    match halo_profile:
        case HaloProfileType.NFW.value:
            assert isinstance(hp, Nc.HaloDensityProfileNFW)
        case HaloProfileType.EINASTO.value:
            assert isinstance(hp, Nc.HaloDensityProfileEinasto)
        case HaloProfileType.HERNQUIST.value:
            assert isinstance(hp, Nc.HaloDensityProfileHernquist)


def test_cluster_wl_app_generate_halo_profile_bad(experiment_file, halo_profile_bad):
    """Test the generation of the cluster WL app with
    bad cluster profile configuration."""
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
            "--cluster-mass",
            str(mass),
            "--cluster-c",
            str(c),
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
    experiment = ser.dict_str_from_yaml_file(experiment_file.as_posix())
    mset = experiment.get("model-set")

    hms = mset.peek_by_name("NcHaloMassSummary")

    assert hms["log10MDelta"] == log10(mass)
    assert hms["cDelta"] == c


def test_cluster_wl_app_halo_mass_summary_bad(experiment_file, halo_mass_summary_bad):
    """Test the generation of the cluster WL app with
    bad halo mass summary configuration."""
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
            "--cluster-ra",
            str(cluster_ra),
            "--cluster-dec",
            str(cluster_dec),
            "--cluster-z",
            str(cluster_z),
            "--ra-min",
            str(ra_min),
            "--ra-max",
            str(ra_max),
            "--dec-min",
            str(dec_min),
            "--dec-max",
            str(dec_max),
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
    experiment = ser.dict_str_from_yaml_file(experiment_file.as_posix())
    mset = experiment.get("model-set")

    halo_position = mset.peek_by_name("NcHaloPosition")

    assert halo_position["ra"] == cluster_ra
    assert halo_position["dec"] == cluster_dec
    assert halo_position["z"] == cluster_z


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
            "--r-min",
            str(r_min),
            "--r-max",
            str(r_max),
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
    experiment = ser.dict_str_from_yaml_file(experiment_file.as_posix())
    mset = experiment.get("model-set")
    cluster_data = dataset.get_data(0)
    obs = cluster_data.peek_obs()

    hp = mset.peek_by_name("NcHaloPosition")
    cosmo = mset.peek_by_name("NcHICosmo")

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
            "--galaxy-density",
            str(galaxy_density),
            "--cluster-ra",
            str(cluster_ra),
            "--cluster-dec",
            str(cluster_dec),
            "--cluster-z",
            str(cluster_z),
            "--ra-min",
            str(ra_min),
            "--ra-max",
            str(ra_max),
            "--dec-min",
            str(dec_min),
            "--dec-max",
            str(dec_max),
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

    assert obs.len() <= n_galaxies


def test_cluster_wl_app_galaxy_density_bad(experiment_file, density_bad):
    """Test the generation of the cluster WL app with
    bad galaxy density configuration."""
    result = runner.invoke(
        app,
        [
            "generate",
            "cluster-wl",
            experiment_file.as_posix(),
            density_bad,
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

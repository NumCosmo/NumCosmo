#!/usr/bin/env python
#
# test_py_cfg.py
#
# Fri Nov 08 15:58:22 2024
# Copyright  2024  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_py_cfg.py
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

"""Tests on ncm_cfg namespace."""

import os
import pytest
from numpy.testing import assert_allclose

from numcosmo_py import Ncm, GLib

Ncm.cfg_init()


@pytest.fixture(
    name="flag_string", params=["estimate", "measure", "patient", "exhaustive"]
)
def fixture_flag_string(request) -> str:
    """Fixture for FFTW flag string."""
    return request.param


@pytest.fixture(name="timelimit", params=[1.23, 4.56])
def fixture_timelimit(request) -> float:
    """Fixture for FFTW time limit."""
    return request.param


def test_set_fftw_default_flag(timelimit: float) -> None:
    """Test setting and getting the default FFTW flag."""
    Ncm.cfg_set_fftw_default_flag(0, timelimit)
    assert Ncm.cfg_get_fftw_default_flag() == 0
    assert_allclose(Ncm.cfg_get_fftw_timelimit(), timelimit)


def test_set_fftw_default_flag_str(flag_string: str, timelimit: float) -> None:
    """Test setting and getting the default FFTW flag."""
    Ncm.cfg_set_fftw_default_flag_str(flag_string, timelimit)
    assert Ncm.cfg_get_fftw_default_flag_str() == flag_string
    assert_allclose(Ncm.cfg_get_fftw_timelimit(), timelimit)

    flag = Ncm.cfg_get_fftw_default_flag()
    Ncm.cfg_set_fftw_default_flag(flag, timelimit)
    assert Ncm.cfg_get_fftw_default_flag_str() == flag_string
    assert Ncm.cfg_get_fftw_default_flag() == flag


def test_set_fftw_invalid_flag(timelimit: float) -> None:
    """Test setting an invalid FFTW flag."""
    with pytest.raises(GLib.Error, match="Invalid FFTW flag '123'"):
        Ncm.cfg_set_fftw_default_flag(123, timelimit)


def test_set_fftw_invalid_flag_str(timelimit: float) -> None:
    """Test setting an invalid FFTW flag."""
    with pytest.raises(GLib.Error, match="Invalid FFTW flag string 'invalid'"):
        Ncm.cfg_set_fftw_default_flag_str("invalid", timelimit)


def test_set_fftw_from_env_no_env(flag_string: str, timelimit: float) -> None:
    """Test setting FFTW flag from environment variable, but using the fallback."""
    os.environ.pop("NCM_FFTW_PLANNER", None)
    os.environ.pop("NCM_FFTW_PLANNER_TIMELIMIT", None)

    Ncm.cfg_set_fftw_default_flag_str(flag_string, timelimit + 10.0)

    Ncm.cfg_set_fftw_default_from_env(0, timelimit)
    assert Ncm.cfg_get_fftw_default_flag() == 0
    assert_allclose(Ncm.cfg_get_fftw_timelimit(), timelimit)


def test_set_fftw_from_env(flag_string: str, timelimit: float) -> None:
    """Test setting FFTW flag from environment variable."""
    os.environ["NCM_FFTW_PLANNER"] = flag_string
    os.environ["NCM_FFTW_PLANNER_TIMELIMIT"] = str(timelimit)

    Ncm.cfg_set_fftw_default_from_env(0, timelimit + 10.0)
    assert Ncm.cfg_get_fftw_default_flag_str() == flag_string
    assert_allclose(Ncm.cfg_get_fftw_timelimit(), timelimit)


def test_set_fftw_from_env_str(flag_string: str, timelimit: float) -> None:
    """Test setting FFTW flag from environment variable."""
    os.environ["NCM_FFTW_PLANNER"] = flag_string
    os.environ["NCM_FFTW_PLANNER_TIMELIMIT"] = str(timelimit)

    Ncm.cfg_set_fftw_default_from_env_str("estimate", timelimit + 10.0)
    assert Ncm.cfg_get_fftw_default_flag_str() == flag_string
    assert_allclose(Ncm.cfg_get_fftw_timelimit(), timelimit)


def test_set_fftw_from_env_invalid_timelimit_str(
    flag_string: str, timelimit: float
) -> None:
    """Test setting FFTW flag from environment variable with invalid timelimit."""
    os.environ["NCM_FFTW_PLANNER"] = flag_string
    os.environ["NCM_FFTW_PLANNER_TIMELIMIT"] = "invalid"

    with pytest.raises(GLib.Error, match="Invalid FFTW planner timelimit 'invalid'"):
        Ncm.cfg_set_fftw_default_from_env(0, timelimit)


def test_set_fftw_from_env_str_invalid_timelimit_str(
    flag_string: str, timelimit: float
) -> None:
    """Test setting FFTW flag from environment variable with invalid timelimit."""
    os.environ["NCM_FFTW_PLANNER"] = flag_string
    os.environ["NCM_FFTW_PLANNER_TIMELIMIT"] = "invalid"

    with pytest.raises(GLib.Error, match="Invalid FFTW planner timelimit 'invalid'"):
        Ncm.cfg_set_fftw_default_from_env_str("estimate", timelimit)


def test_get_version() -> None:
    """Test getting the version string."""
    version = Ncm.cfg_get_version()
    # This functions returns (10000 * major + 100 * minor + micro, major, minor, micro)
    assert version == (2200, 0, 22, 0)


def test_get_git_hash() -> None:
    """Test getting the git hash."""
    git_hash = Ncm.cfg_get_commit_hash()
    assert len(git_hash) > 0


def test_get_version_string() -> None:
    """Test getting the version string."""
    version = Ncm.cfg_get_version_string()
    assert version == "0.22.0"


def test_check_version() -> None:
    """Test checking the version."""
    version = Ncm.cfg_get_version()

    assert Ncm.cfg_version_check(version[1], version[2], version[3])
    # Check for a version that is not the current one, taking care of the
    # possibility of major, macro or micro version being zero.
    assert not Ncm.cfg_version_check(version[1] + 1, version[2], version[3])
    assert not Ncm.cfg_version_check(version[1], version[2] + 1, version[3])
    assert not Ncm.cfg_version_check(version[1], version[2], version[3] + 1)
    assert not Ncm.cfg_version_check(version[1] + 1, version[2] + 1, version[3] + 1)
    assert not Ncm.cfg_version_check(version[1] + 1, version[2], version[3] + 1)
    assert not Ncm.cfg_version_check(version[1], version[2] + 1, version[3] + 1)
    assert not Ncm.cfg_version_check(version[1] + 1, version[2] + 1, version[3])

    # Now prior versions
    if version[1] > 0:
        assert Ncm.cfg_version_check(version[1] - 1, version[2], version[3])
        assert Ncm.cfg_version_check(version[1] - 1, version[2] + 1, version[3])
        assert Ncm.cfg_version_check(version[1] - 1, version[2], version[3] + 1)
    if version[2] > 0:
        assert Ncm.cfg_version_check(version[1], version[2] - 1, version[3])
        assert Ncm.cfg_version_check(version[1], version[2] - 1, version[3] + 1)
    if version[3] > 0:
        assert Ncm.cfg_version_check(version[1], version[2], version[3] - 1)

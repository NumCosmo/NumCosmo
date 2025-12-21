#!/usr/bin/env python
#
# test_py_generate.py
#
# Wed Jan 31 10:04:00 2024
# Copyright  2023  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_py_generate.py
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

"""Test Generate classes."""

from typing import cast
from pathlib import Path
import pytest

from numcosmo_py import Ncm
import numcosmo_py.app.generate as gen
import numcosmo_py.datasets.hicosmo as hicosmo

Ncm.cfg_init()


# JPAS Forecast
def test_generate_jpas(tmp_path: Path):
    """Test JPAS forecast generation with YAML file."""
    exp_file = tmp_path / "jpas_forecast.yaml"
    _ = gen.GenerateJpasForecast(experiment=exp_file.absolute())

    assert exp_file.exists()


def test_generate_jpas_invalid_suffix(tmp_path: Path):
    """Test JPAS forecast rejects invalid file suffix."""
    exp_file = tmp_path / "jpas_forecast.txt"
    with pytest.raises(ValueError, match="Invalid experiment file suffix: .txt"):
        _ = gen.GenerateJpasForecast(experiment=exp_file.absolute())


# QSpline
def test_generate_qspline(tmp_path: Path):
    """Test QSpline generation with all data types."""
    exp_file = tmp_path / "qspline.yaml"
    bao_id = hicosmo.BAOID.ALL_COMBINED_JUN_2025
    sne_id = hicosmo.SNIaID.COV_DES_Y5_STAT_SYS
    h_id = hicosmo.HID.ALL_COMBINED_APR_2025
    _ = gen.GenerateQSpline(
        experiment=exp_file.absolute(),
        include_snia=sne_id,
        include_hubble=h_id,
        include_bao=bao_id,
    )

    assert exp_file.exists()


def test_generate_qspline_no_data(tmp_path: Path):
    """Test QSpline rejects experiment with no data."""
    exp_file = tmp_path / "qspline.yaml"
    with pytest.raises(ValueError, match="No data included in the experiment."):
        _ = gen.GenerateQSpline(experiment=exp_file.absolute())


def test_generate_qspline_invalid_suffix(tmp_path: Path):
    """Test QSpline rejects invalid file suffix."""
    exp_file = tmp_path / "qspline.txt"
    bao_id = hicosmo.BAOID.ALL_COMBINED_JUN_2025
    with pytest.raises(ValueError, match="Invalid experiment file suffix: .txt"):
        _ = gen.GenerateQSpline(experiment=exp_file.absolute(), include_bao=bao_id)


# XCDM
def test_generate_xcdm(tmp_path: Path):
    """Test XCDM generation with multiple data sets."""
    exp_file = tmp_path / "xcdm.yaml"
    bao_ids = [hicosmo.BAOID.ALL_COMBINED_JUN_2025, hicosmo.BAOID.ALL_COMBINED_JAN_2023]
    sne_ids = [hicosmo.SNIaID.COV_DES_Y5_STAT_SYS, hicosmo.SNIaID.COV_DES_Y5_STATONLY]
    h_ids = [hicosmo.HID.ALL_COMBINED_APR_2025, hicosmo.HID.ALL_COMBINED_JAN_2023]

    for bao_id, sne_id, h_id in zip(bao_ids, sne_ids, h_ids):
        _ = gen.GenerateXCDM(
            experiment=exp_file.absolute(),
            include_snia=sne_id,
            include_hubble=h_id,
            include_bao=bao_id,
        )

    assert exp_file.exists()


def test_generate_xcdm_no_data(tmp_path: Path):
    """Test XCDM rejects experiment with no data."""
    exp_file = tmp_path / "xcdm.yaml"
    with pytest.raises(ValueError, match="No data included in the experiment."):
        _ = gen.GenerateXCDM(experiment=exp_file.absolute())


def test_generate_xcdm_no_valid_bao_id(tmp_path: Path):
    """Test XCDM rejects invalid BAO ID."""
    exp_file = tmp_path / "xcdm.yaml"
    bao_id = cast(hicosmo.BAOID, str(-1))
    with pytest.raises(ValueError, match="Unknown BAO data set id: -1"):
        _ = gen.GenerateXCDM(experiment=exp_file.absolute(), include_bao=bao_id)


def test_generate_xcdm_no_valid_hubble_id(tmp_path: Path):
    """Test XCDM rejects invalid Hubble ID."""
    exp_file = tmp_path / "xcdm.yaml"
    h_id = cast(hicosmo.HID, str(-1))
    with pytest.raises(ValueError, match="Unknown Hubble data set id: -1"):
        _ = gen.GenerateXCDM(experiment=exp_file.absolute(), include_hubble=h_id)


def test_generate_xcdm_invalid_suffix(tmp_path: Path):
    """Test XCDM rejects invalid file suffix."""
    exp_file = tmp_path / "xcdm.txt"
    bao_id = hicosmo.BAOID.ALL_COMBINED_JUN_2025
    with pytest.raises(ValueError, match="Invalid experiment file suffix: .txt"):
        _ = gen.GenerateXCDM(experiment=exp_file.absolute(), include_bao=bao_id)


# DEWSpline
# QSpline
def test_generate_dewspline(tmp_path: Path):
    """Test DEWSpline generation with all data types."""
    exp_file = tmp_path / "dewspline.yaml"
    bao_id = hicosmo.BAOID.ALL_COMBINED_JUN_2025
    sne_id = hicosmo.SNIaID.COV_DES_Y5_STAT_SYS
    h_id = hicosmo.HID.ALL_COMBINED_APR_2025
    _ = gen.GenerateDEWSpline(
        experiment=exp_file.absolute(),
        include_snia=sne_id,
        include_hubble=h_id,
        include_bao=bao_id,
    )

    assert exp_file.exists()


def test_generate_dewspline_no_data(tmp_path: Path):
    """Test DEWSpline rejects experiment with no data."""
    exp_file = tmp_path / "qspline.yaml"
    with pytest.raises(ValueError, match="No data included in the experiment."):
        _ = gen.GenerateDEWSpline(experiment=exp_file.absolute())


def test_generate_dewspline_invalid_suffix(tmp_path: Path):
    """Test DEWSpline rejects invalid file suffix."""
    exp_file = tmp_path / "dewspline.txt"
    bao_id = hicosmo.BAOID.ALL_COMBINED_JUN_2025
    with pytest.raises(ValueError, match="Invalid experiment file suffix: .txt"):
        _ = gen.GenerateDEWSpline(experiment=exp_file.absolute(), include_bao=bao_id)

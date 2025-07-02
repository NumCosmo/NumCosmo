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

from pathlib import Path
import pytest

from numcosmo_py import Ncm, Nc
import numcosmo_py.app.generate as gen
import numcosmo_py.datasets.hicosmo as hicosmo

Ncm.cfg_init()


# JPAS Forecast
def test_generate_jpas(tmp_path: Path):
    exp_file = tmp_path / "jpas_forecast.yaml"
    exp = gen.GenerateJpasForecast(experiment=exp_file.absolute())

    assert exp_file.exists()


def test_generate_jpas_invalid_suffix(tmp_path: Path):
    exp_file = tmp_path / "jpas_forecast.txt"
    with pytest.raises(ValueError, match="Invalid experiment file suffix: .txt"):
        exp = gen.GenerateJpasForecast(experiment=exp_file.absolute())


# QSpline
def test_generate_qspline(tmp_path: Path):
    exp_file = tmp_path / "qspline.yaml"
    bao_id = hicosmo.BAOID.ALL_COMBINED_JUN_2025
    sne_id = hicosmo.SNIaID.COV_DES_Y5_STAT_SYS
    h_id = hicosmo.HID.ALL_COMBINED_APR_2025
    exp = gen.GenerateQSpline(
        experiment=exp_file.absolute(),
        include_snia=sne_id,
        include_hubble=h_id,
        include_bao=bao_id,
    )

    assert exp_file.exists()


def test_generate_qspline_no_data(tmp_path: Path):
    exp_file = tmp_path / "qspline.yaml"
    with pytest.raises(ValueError, match="No data included in the experiment."):
        exp = gen.GenerateQSpline(experiment=exp_file.absolute())


def test_generate_qspline_invalid_suffix(tmp_path: Path):
    exp_file = tmp_path / "qspline.txt"
    bao_id = hicosmo.BAOID.ALL_COMBINED_JUN_2025
    with pytest.raises(ValueError, match="Invalid experiment file suffix: .txt"):
        exp = gen.GenerateQSpline(experiment=exp_file.absolute(), include_bao=bao_id)


# XCDM
def test_generate_xcdm(tmp_path: Path):
    exp_file = tmp_path / "xcdm.yaml"
    bao_id = [hicosmo.BAOID.ALL_COMBINED_JUN_2025, hicosmo.BAOID.ALL_COMBINED_JAN_2023]
    sne_id = [hicosmo.SNIaID.COV_DES_Y5_STAT_SYS, hicosmo.SNIaID.COV_DES_Y5_STATONLY]
    h_id = [hicosmo.HID.ALL_COMBINED_APR_2025, hicosmo.HID.ALL_COMBINED_JAN_2023]

    for i in range(len(bao_id)):
        exp = gen.GenerateXCDM(
            experiment=exp_file.absolute(),
            include_snia=sne_id[i],
            include_hubble=h_id[i],
            include_bao=bao_id[i],
        )

    assert exp_file.exists()


def test_generate_xcdm_no_data(tmp_path: Path):
    exp_file = tmp_path / "xcdm.yaml"
    with pytest.raises(ValueError, match="No data included in the experiment."):
        exp = gen.GenerateXCDM(experiment=exp_file.absolute())


def test_generate_xcdm_no_valid_bao_id(tmp_path: Path):
    exp_file = tmp_path / "xcdm.yaml"
    bao_id = str(-1)
    with pytest.raises(ValueError, match="Unknown BAO data set id: -1"):
        exp = gen.GenerateXCDM(experiment=exp_file.absolute(), include_bao=bao_id)


def test_generate_xcdm_no_valid_hubble_id(tmp_path: Path):
    exp_file = tmp_path / "xcdm.yaml"
    h_id = str(-1)
    with pytest.raises(ValueError, match="Unknown Hubble data set id: -1"):
        exp = gen.GenerateXCDM(experiment=exp_file.absolute(), include_hubble=h_id)


def test_generate_xcdm_invalid_suffix(tmp_path: Path):
    exp_file = tmp_path / "xcdm.txt"
    bao_id = hicosmo.BAOID.ALL_COMBINED_JUN_2025
    with pytest.raises(ValueError, match="Invalid experiment file suffix: .txt"):
        exp = gen.GenerateXCDM(experiment=exp_file.absolute(), include_bao=bao_id)

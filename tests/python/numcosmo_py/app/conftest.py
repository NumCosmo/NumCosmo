#!/usr/bin/env python
#
# conftest.py
#
# Tue Nov 3 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# conftest.py
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

"""Fixtures shared by the numcosmo app tests."""

from pathlib import Path

import pytest

from numcosmo_py import GObject, Ncm


@pytest.fixture(name="builder_experiment")
def fixture_builder_experiment(tmp_path: Path) -> Path:
    """Write an experiment whose model set holds a model of its builders file.

    The model set holds a two-parameter NcmModelMVND, both parameters free, and a
    model of the type NcmTestAppBuilderModel created by the builder written to the
    .builders.yaml file next to the experiment.
    """
    model_builder = Ncm.ModelBuilder.new(
        Ncm.Model, "NcmTestAppBuilderModel", "Model created by a builder"
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
    model_builders.add("NcmTestAppBuilderModel", model_builder)
    ser.dict_str_to_yaml_file(
        model_builders, experiment_file.with_suffix(".builders.yaml").as_posix()
    )

    experiment = Ncm.ObjDictStr.new()
    experiment.add("model-set", mset)
    experiment.add("likelihood", likelihood)
    ser.dict_str_to_yaml_file(experiment, experiment_file.as_posix())

    return experiment_file

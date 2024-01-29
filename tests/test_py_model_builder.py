#!/usr/bin/env python
#
# test_py_model_builder.py
#
# Mon Jan 29 17:31:00 2024
# Copyright  2024  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_py_model_builder.py
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

"""Tests on NcmModelBuilder class."""


from numpy.testing import assert_allclose
from numcosmo_py import Ncm, GObject


def test_create_with_scalar():
    """Test create simple model."""
    mb = Ncm.ModelBuilder.new(Ncm.Model, "MyNewModel", "This is my new model")
    assert mb is not None
    assert isinstance(mb, Ncm.ModelBuilder)

    mb.add_sparam("p_1", "param1", 0.0, 1.0, 0.01, 0.0, 0.5, Ncm.ParamType.FREE)

    MyNewModel = mb.create()
    assert MyNewModel is not None
    GObject.new(MyNewModel)
    model = MyNewModel.pytype()

    assert model is not None
    assert isinstance(model, Ncm.Model)
    assert model.name() == "MyNewModel"
    assert isinstance(model, MyNewModel.pytype)


def test_create_with_vector():
    """Test create model with vector."""
    mb = Ncm.ModelBuilder.new(Ncm.Model, "MyNewModel2", "This is my new model")
    assert mb is not None
    assert isinstance(mb, Ncm.ModelBuilder)

    mb.add_sparam("p_1", "param1", 0.0, 1.0, 0.01, 0.0, 0.45, Ncm.ParamType.FREE)
    mb.add_vparam(10, "p_2", "param2", 0.0, 1.0, 0.01, 0.0, 0.5, Ncm.ParamType.FREE)

    MyNewModel2 = mb.create()
    assert MyNewModel2 is not None
    GObject.new(MyNewModel2)
    model = MyNewModel2.pytype()

    assert model is not None
    assert isinstance(model, Ncm.Model)
    assert model.name() == "MyNewModel2"
    assert isinstance(model, MyNewModel2.pytype)

    assert_allclose(model.orig_param_get(0), 0.45)

    assert_allclose(model.vparam_len(0), 10)
    assert_allclose(model.orig_vparam_get_vector(0).dup_array(), 0.5)


def test_serialization():
    """Test model builder serialization."""

    mb = Ncm.ModelBuilder.new(Ncm.Model, "MyNewModel3", "This is my new model")
    assert mb is not None
    assert isinstance(mb, Ncm.ModelBuilder)

    mb.add_sparam("p_1", "param1", 0.0, 1.0, 0.01, 0.0, 0.45, Ncm.ParamType.FREE)
    mb.add_vparam(10, "p_2", "param2", 0.0, 1.0, 0.01, 0.0, 0.5, Ncm.ParamType.FREE)

    ser = Ncm.Serialize.new(Ncm.SerializeOpt.NONE)

    mb_dup = ser.dup_obj(mb)

    assert mb_dup is not None
    assert isinstance(mb_dup, Ncm.ModelBuilder)

    MyNewModel3 = mb_dup.create()
    assert MyNewModel3 is not None
    GObject.new(MyNewModel3)
    model = MyNewModel3.pytype()

    assert model is not None
    assert isinstance(model, Ncm.Model)
    assert model.name() == "MyNewModel3"
    assert isinstance(model, MyNewModel3.pytype)

    assert_allclose(model.orig_param_get(0), 0.45)

    assert_allclose(model.vparam_len(0), 10)
    assert_allclose(model.orig_vparam_get_vector(0).dup_array(), 0.5)

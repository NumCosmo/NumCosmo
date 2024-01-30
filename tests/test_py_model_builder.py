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


from typing import cast

from numpy.testing import assert_allclose
from numcosmo_py import Ncm, GObject


def test_create_with_scalar() -> None:
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


def test_create_with_vector() -> None:
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


def test_serialization() -> None:
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


def test_get_many_sparams() -> None:
    """Test many sparams."""
    mb = Ncm.ModelBuilder.new(Ncm.Model, "MyNewModel4", "This is my new model")
    assert mb is not None
    assert isinstance(mb, Ncm.ModelBuilder)

    for i in range(100):
        mb.add_sparam(
            f"p_{i}", f"param{i}", 0.0, 1.0, 0.01, 0.0, 0.45, Ncm.ParamType.FREE
        )

    obj_array = mb.get_sparams()

    assert obj_array is not None
    assert isinstance(obj_array, Ncm.ObjArray)

    assert obj_array.len() == 100

    for i in range(100):
        sparam: Ncm.SParam = cast(Ncm.SParam, obj_array.peek(i))
        assert sparam is not None
        assert isinstance(sparam, Ncm.SParam)

        assert sparam.symbol() == f"p_{i}"
        assert sparam.name() == f"param{i}"
        assert sparam.get_fit_type() == Ncm.ParamType.FREE
        assert sparam.get_default_value() == 0.45
        assert sparam.get_lower_bound() == 0.0
        assert sparam.get_upper_bound() == 1.0
        assert sparam.get_scale() == 0.01
        assert sparam.get_absolute_tolerance() == 0.0


def test_set_many_sparams() -> None:
    """Test many sparams."""

    obj_array = Ncm.ObjArray.new()
    assert obj_array is not None

    for i in range(100):
        sparam = Ncm.SParam.new(
            f"param{i}", f"p_{i}", 0.0, 1.0, 0.01, 0.0, 0.45, Ncm.ParamType.FREE
        )
        assert sparam is not None
        assert isinstance(sparam, Ncm.SParam)
        obj_array.add(sparam)

    mb = Ncm.ModelBuilder.new(Ncm.Model, "MyNewModel5", "This is my new model")
    assert mb is not None

    mb.add_sparams(obj_array)

    MyNewModel5 = mb.create()
    assert MyNewModel5 is not None
    GObject.new(MyNewModel5)
    model = MyNewModel5.pytype()

    assert model is not None
    assert isinstance(model, Ncm.Model)
    assert model.name() == "MyNewModel5"
    assert isinstance(model, MyNewModel5.pytype)

    model.params_set_default_ftype()

    for i in range(100):
        assert model.orig_param_get(i) == 0.45
        assert model.orig_param_get_abstol(i) == 0.0
        assert model.orig_param_get_lower_bound(i) == 0.0
        assert model.orig_param_get_upper_bound(i) == 1.0
        assert model.orig_param_get_scale(i) == 0.01
        assert model.param_index_from_name(f"param{i}") == (True, i)
        assert model.param_get_ftype(i) == Ncm.ParamType.FREE

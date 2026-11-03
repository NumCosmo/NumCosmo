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

import pytest
from gi.repository import GLib
from numpy.testing import assert_allclose
from numcosmo_py import Ncm
from numcosmo_py.helper import duplicate_via_serialization, GObject
from numcosmo_py.helper import register_model_class


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


def test_register_model_class() -> None:
    """Test create simple model."""
    mb = Ncm.ModelBuilder.new(Ncm.Model, "MyNewModelRMC", "This is my new model")
    assert mb is not None
    assert isinstance(mb, Ncm.ModelBuilder)

    mb.add_sparam("p_1", "param1", 0.0, 1.0, 0.01, 0.0, 0.5, Ncm.ParamType.FREE)

    NewModel = register_model_class(mb)
    model = NewModel()

    assert model is not None
    assert isinstance(model, Ncm.Model)
    assert model.name() == "MyNewModelRMC"
    assert isinstance(model, NewModel)


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

    mb_dup = duplicate_via_serialization(mb, ser)

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


def _builder(
    name: str, default: float = 0.5, parent: str = "NcmModel"
) -> Ncm.ModelBuilder:
    """Build a model builder with one scalar parameter of default value default."""
    mb = Ncm.ModelBuilder(parent_type_string=parent, name=name, description="Test")
    mb.add_sparam("p_1", "param1", 0.0, 1.0, 0.01, 0.0, default, Ncm.ParamType.FREE)
    return mb


def test_create_same_definition_returns_existing_type() -> None:
    """A second builder with the same definition gets the existing type."""
    mb1 = _builder("MyNewModelSame")
    mb2 = _builder("MyNewModelSame")

    assert mb1.create() == mb2.create()


def test_create_name_clash_raises() -> None:
    """A builder whose parameters differ from the registered type raises."""
    mb1 = _builder("MyNewModelClash", 0.5)
    mb2 = _builder("MyNewModelClash", 0.7)

    mb1.create()
    with pytest.raises(GLib.Error, match="scalar parameter 0 .* differs"):
        mb2.create()


def test_create_parent_clash_raises() -> None:
    """A builder whose parent differs from the registered type raises."""
    mb1 = _builder("MyNewModelClashParent")
    mb2 = _builder("MyNewModelClashParent", parent="NcmModelRosenbrock")

    mb1.create()
    with pytest.raises(GLib.Error, match="already exists with parent `NcmModel'"):
        mb2.create()


def test_create_stackable_clash_raises() -> None:
    """A builder whose stackable flag differs from the registered type raises."""
    mb1 = _builder("MyNewModelClashStack")
    mb2 = Ncm.ModelBuilder(
        parent_type_string="NcmModel",
        name="MyNewModelClashStack",
        description="Test",
        stackable=True,
    )
    mb2.add_sparam("p_1", "param1", 0.0, 1.0, 0.01, 0.0, 0.5, Ncm.ParamType.FREE)

    mb1.create()
    with pytest.raises(
        GLib.Error, match="stackable FALSE, but this builder has stackable TRUE"
    ):
        mb2.create()


def test_create_unknown_parent_raises() -> None:
    """A parent name that is not a registered type raises at create."""
    mb = Ncm.ModelBuilder(parent_type_string="NoSuchModel", name="MyNewModelNoParent")

    with pytest.raises(
        GLib.Error, match="parent type `NoSuchModel' .* is not registered"
    ):
        mb.create()


def test_create_parent_not_model_raises() -> None:
    """A parent that is not a NcmModel raises at create."""
    mb = Ncm.ModelBuilder(parent_type_string="NcmVector", name="MyNewModelVectorParent")

    with pytest.raises(
        GLib.Error, match="parent type `NcmVector' .* is not a NcmModel"
    ):
        mb.create()


def test_peek_from_type() -> None:
    """A type created by a builder gives that builder back; other types give None."""
    mb1 = _builder("MyNewModelPeek")
    mb2 = _builder("MyNewModelPeek")
    gtype = mb1.create()

    assert Ncm.ModelBuilder.peek_from_type(gtype) is mb1
    assert mb2.create() == gtype
    assert Ncm.ModelBuilder.peek_from_type(gtype) is mb1
    assert Ncm.ModelBuilder.peek_from_type(Ncm.Model.__gtype__) is None
    assert Ncm.ModelBuilder.peek_from_type(Ncm.ModelRosenbrock.__gtype__) is None

    GObject.new(gtype)
    model = gtype.pytype()
    assert Ncm.ModelBuilder.peek_from_type(model.__gtype__) is mb1


def test_deserialize_before_parent_is_created() -> None:
    """Builders deserialize before their parent type exists; created in order."""
    parent = _builder("MyNewModelChainParent")
    child = Ncm.ModelBuilder(
        parent_type_string="MyNewModelChainParent",
        name="MyNewModelChainChild",
        description="Test",
    )
    child.add_sparam("p_2", "param2", 0.0, 1.0, 0.01, 0.0, 0.2, Ncm.ParamType.FREE)

    builders = Ncm.ObjArray.new()
    builders.add(parent)
    builders.add(child)

    ser = Ncm.Serialize.new(Ncm.SerializeOpt.NONE)
    yaml = ser.array_to_yaml(builders)
    ser.reset(True)

    with pytest.raises(RuntimeError, match="unknown type name"):
        GObject.type_from_name("MyNewModelChainParent")
    builders_dup = ser.array_from_yaml(yaml)

    parent_type = cast(Ncm.ModelBuilder, builders_dup.peek(0)).create()
    child_type = cast(Ncm.ModelBuilder, builders_dup.peek(1)).create()

    assert child_type.parent == parent_type
    assert Ncm.ModelBuilder.peek_from_type(parent_type) is builders_dup.peek(0)
    assert Ncm.ModelBuilder.peek_from_type(child_type) is builders_dup.peek(1)

    GObject.new(child_type)
    model = child_type.pytype()
    assert model.len() == 2
    assert model.param_name(0) == "param1"
    assert model.param_name(1) == "param2"

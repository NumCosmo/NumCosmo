#!/usr/bin/env python
#
# test_py_serialize.py
#
# Mon Jan 15 10:19:40 2024
# Copyright  2023  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_py_serialize.py
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

"""Tests for serialization."""

import pytest
import numpy as np
from numcosmo_py import Ncm, GObject, dict_to_var_dict

Ncm.cfg_init()

# pylint: disable=no-member


class GTestA(GObject.Object):
    """Test class for serialization."""

    __gtype_name__ = "GTestA"

    A_string: str = GObject.Property(
        type=GObject.TYPE_STRING,  # type: ignore
        flags=GObject.ParamFlags.READWRITE | GObject.ParamFlags.CONSTRUCT,
    )

    A_int: int = GObject.Property(
        type=GObject.TYPE_INT,  # type: ignore
        flags=GObject.ParamFlags.READWRITE | GObject.ParamFlags.CONSTRUCT,
    )

    A_double: float = GObject.Property(
        type=GObject.TYPE_DOUBLE,  # type: ignore
        flags=GObject.ParamFlags.READWRITE | GObject.ParamFlags.CONSTRUCT,
    )

    A_vector: Ncm.Vector = GObject.Property(
        type=Ncm.Vector,  # type: ignore
        flags=GObject.ParamFlags.READWRITE | GObject.ParamFlags.CONSTRUCT,
    )

    A_matrix: Ncm.Matrix = GObject.Property(
        type=Ncm.Matrix,  # type: ignore
        flags=GObject.ParamFlags.READWRITE | GObject.ParamFlags.CONSTRUCT,
    )

    A_var_dict: Ncm.VarDict = GObject.Property(
        type=Ncm.VarDict,  # type: ignore
        flags=GObject.ParamFlags.READWRITE | GObject.ParamFlags.CONSTRUCT,
    )

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.extra = 4244

    def __eq__(self, other):
        if not isinstance(other, GTestA):
            return False

        t1 = self.A_string == other.A_string
        t2 = self.A_int == other.A_int
        t3 = self.A_double == other.A_double
        t4 = np.allclose(self.A_vector.dup_array(), other.A_vector.dup_array())
        t5 = np.allclose(self.A_matrix.dup_array(), other.A_matrix.dup_array())
        t6 = self.extra == other.extra

        t7 = self.A_var_dict.len() == other.A_var_dict.len()
        t8 = all(
            self.A_var_dict.get_variant(key) == other.A_var_dict.get_variant(key)
            for key in self.A_var_dict.keys()
        )

        return t1 and t2 and t3 and t4 and t5 and t6 and t7 and t8


GObject.type_register(GTestA)


class GTestB(GObject.Object):
    """Test class for serialization with container properties."""

    __gtype_name__ = "GTestB"

    B_string = GObject.Property(
        type=GObject.TYPE_STRING,  # type: ignore
        flags=GObject.ParamFlags.READWRITE | GObject.ParamFlags.CONSTRUCT,
    )

    As = GObject.Property(
        type=GTestA,  # type: ignore
        flags=GObject.ParamFlags.READWRITE | GObject.ParamFlags.CONSTRUCT,
    )

    Aa = GObject.Property(
        type=Ncm.ObjArray,  # type: ignore
        flags=GObject.ParamFlags.READWRITE | GObject.ParamFlags.CONSTRUCT,
    )

    Aods = GObject.Property(
        type=Ncm.ObjDictStr,  # type: ignore
        flags=GObject.ParamFlags.READWRITE | GObject.ParamFlags.CONSTRUCT,
    )

    Aodi = GObject.Property(
        type=Ncm.ObjDictInt,  # type: ignore
        flags=GObject.ParamFlags.READWRITE | GObject.ParamFlags.CONSTRUCT,
    )

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, GTestB):
            return False

        t1 = self.B_string == other.B_string
        t2 = self.As == other.As
        if self.Aa.len() != other.Aa.len():
            return False
        t3 = all(self.Aa.get(i) == other.Aa.get(i) for i in range(self.Aa.len()))

        if self.Aods.len() != other.Aods.len():
            return False
        t4 = all(
            self.Aods.peek(key) == other.Aods.peek(key) for key in self.Aods.keys()
        )

        if self.Aodi.len() != other.Aodi.len():
            return False
        t5 = all(
            self.Aodi.peek(key) == other.Aodi.peek(key) for key in self.Aodi.keys()
        )
        return t1 and t2 and t3 and t4 and t5


GObject.type_register(GTestB)


class GTestC(GObject.Object):
    """Test class for serialization with container properties."""

    __gtype_name__ = "GTestC"

    Ca = GObject.Property(
        type=Ncm.ObjArray,  # type: ignore
        flags=GObject.ParamFlags.READWRITE | GObject.ParamFlags.CONSTRUCT,
    )

    Cods = GObject.Property(
        type=Ncm.ObjDictStr,  # type: ignore
        flags=GObject.ParamFlags.READWRITE | GObject.ParamFlags.CONSTRUCT,
    )

    Codi = GObject.Property(
        type=Ncm.ObjDictInt,  # type: ignore
        flags=GObject.ParamFlags.READWRITE | GObject.ParamFlags.CONSTRUCT,
    )

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, GTestC):
            return False

        if self.Ca.len() != other.Ca.len():
            return False
        t1 = all(self.Ca.get(i) == other.Ca.get(i) for i in range(self.Ca.len()))

        if self.Cods.len() != other.Cods.len():
            return False
        t2 = all(
            self.Cods.peek(key) == other.Cods.peek(key) for key in self.Cods.keys()
        )

        if self.Codi.len() != other.Codi.len():
            return False
        t3 = all(
            self.Codi.peek(key) == other.Codi.peek(key) for key in self.Codi.keys()
        )
        return t1 and t2 and t3


GObject.type_register(GTestC)


@pytest.fixture(
    name="ser_opt",
    params=[Ncm.SerializeOpt.NONE, Ncm.SerializeOpt.CLEAN_DUP],
    ids=["NONE", "CLEAN_DUP"],
)
def fixture_ser_opt(request):
    """Fixture for serialization options."""
    return request.param


@pytest.fixture(
    name="dup_mode",
    params=["variant", "yaml"],
)
def fixture_dup_mode(request):
    """Fixture for serialization options."""
    return request.param


@pytest.fixture(
    name="a1",
)
def fixture_a1():
    """Fixture for serialization options."""
    return GTestA(  # type: ignore
        A_string="test",
        A_int=42,
        A_double=3.1415,
        A_vector=Ncm.Vector.new_array([1.0, 2.0, 3.0]),
        A_matrix=Ncm.Matrix.new_array([1.0, 2.0, 3.0, 4.0], 2),
        A_var_dict=dict_to_var_dict(
            {
                "a": 1,
                "b": 2,
                "c": 3,
                "d": "string",
                "e": True,
                "f": 3.14,
                "g": 3.14e-10,
                "h": 3.14e10,
                "i": [1, 2, 3],
                "j": [1.2, 2.3, 3.4],
                "k": [True, False, True],
            }
        ),
    )


@pytest.fixture(
    name="a2",
)
def fixture_a2():
    """Fixture for serialization options."""
    return GTestA(  # type: ignore
        A_string="test2",
        A_int=43,
        A_double=3.1416,
        A_vector=Ncm.Vector.new_array([1.1, 2.1, 3.1]),
        A_matrix=Ncm.Matrix.new_array([1.1, 2.1, 3.1, 4.1], 2),
        A_var_dict=dict_to_var_dict(
            {
                "jaba_a": 10,
                "jaba_b": 22,
                "jaba_c": 3566,
                "jaba_d": "another string",
                "jaba_e": False,
                "jaba_f": 3.1433,
                "jaba_g": 3.14e-30,
                "jaba_h": 3.14e100,
                "jaba_i": [-1, 20, 30000],
                "jaba_j": [1.223, 234.3563, 64643.4],
                "jaba_k": [True, True, True],
            }
        ),
    )


@pytest.fixture(
    name="a3",
)
def fixture_a3():
    """Fixture for serialization options."""
    return GTestA(  # type: ignore
        A_string="test3",
        A_int=44,
        A_double=3.1417,
        A_vector=Ncm.Vector.new_array([1.2, 2.2, 3.2]),
        A_matrix=Ncm.Matrix.new_array([1.2, 2.2, 3.2, 4.2], 2),
        A_var_dict=dict_to_var_dict(
            {
                "jaba_1": 100,
                "jaba_2": 220,
                "jaba_3": 35660,
                "jaba_4": "another string again",
                "jaba_5": False,
                "jaba_6": 34.1433,
                "jaba_7": 36.14e-30,
                "jaba_8": 37.14e100,
                "jaba_9": [-10, 2220, 3000],
                "jaba_0": [1.23, 234.363, 643.4],
                "jaba_11": [True, True, False],
            }
        ),
    )


def test_serialize_A(ser_opt, dup_mode, a1) -> None:
    """Test function for serialization."""

    ser = Ncm.Serialize.new(ser_opt)

    if dup_mode == "variant":
        a1_dup = ser.dup_obj(a1)
    elif dup_mode == "yaml":
        a1_yaml = ser.to_yaml(a1)
        a1_dup = ser.from_yaml(a1_yaml)
    else:
        raise ValueError(f"Invalid dup_mode: {dup_mode}")

    assert isinstance(a1_dup, GTestA)

    assert a1_dup == a1


def test_serialize_B_scalar(ser_opt, dup_mode, a1) -> None:
    """Test function for serialization."""

    ods = Ncm.ObjDictStr.new()
    odi = Ncm.ObjDictInt.new()
    oa = Ncm.ObjArray.new()

    b = GTestB(B_string="test", As=a1, Aa=oa, Aods=ods, Aodi=odi)  # type: ignore

    ser = Ncm.Serialize.new(ser_opt)

    if dup_mode == "variant":
        b_dup = ser.dup_obj(b)
    elif dup_mode == "yaml":
        b_yaml = ser.to_yaml(b)
        b_dup = ser.from_yaml(b_yaml)
    else:
        raise ValueError(f"Invalid dup_mode: {dup_mode}")

    assert isinstance(b_dup, GTestB)

    assert b_dup == b


def test_serialize_B_array(ser_opt, dup_mode, a1, a2, a3) -> None:
    """Test function for serialization."""

    ods = Ncm.ObjDictStr.new()
    odi = Ncm.ObjDictInt.new()
    oa = Ncm.ObjArray.new()
    oa.add(a1)
    oa.add(a2)
    oa.add(a3)

    b = GTestB(B_string="test", As=a1, Aa=oa, Aods=ods, Aodi=odi)  # type: ignore

    ser = Ncm.Serialize.new(ser_opt)

    if dup_mode == "variant":
        b_dup = ser.dup_obj(b)
    elif dup_mode == "yaml":
        b_yaml = ser.to_yaml(b)
        b_dup = ser.from_yaml(b_yaml)
    else:
        raise ValueError(f"Invalid dup_mode: {dup_mode}")

    assert isinstance(b_dup, GTestB)

    assert b_dup == b


def test_serialize_B_dict_str(ser_opt, dup_mode, a1, a2, a3) -> None:
    """Test function for serialization."""

    ods = Ncm.ObjDictStr.new()
    ods.set("a1", a1)
    ods.set("a2", a2)
    ods.set("a3", a3)
    odi = Ncm.ObjDictInt.new()
    oa = Ncm.ObjArray.new()

    b = GTestB(B_string="test", As=a1, Aa=oa, Aods=ods, Aodi=odi)  # type: ignore

    ser = Ncm.Serialize.new(ser_opt)

    if dup_mode == "variant":
        b_dup = ser.dup_obj(b)
    elif dup_mode == "yaml":
        b_yaml = ser.to_yaml(b)
        b_dup = ser.from_yaml(b_yaml)
    else:
        raise ValueError(f"Invalid dup_mode: {dup_mode}")

    assert isinstance(b_dup, GTestB)

    assert b_dup.B_string == "test"

    assert b_dup == b


def test_serialize_B_dict_int(ser_opt, dup_mode, a1, a2, a3) -> None:
    """Test function for serialization."""

    ods = Ncm.ObjDictStr.new()
    odi = Ncm.ObjDictInt.new()
    odi.set(1001, a1)
    odi.set(1002, a2)
    odi.set(1003, a3)
    oa = Ncm.ObjArray.new()

    b = GTestB(B_string="test", As=a1, Aa=oa, Aods=ods, Aodi=odi)  # type: ignore

    ser = Ncm.Serialize.new(ser_opt)

    if dup_mode == "variant":
        b_dup = ser.dup_obj(b)
    elif dup_mode == "yaml":
        b_yaml = ser.to_yaml(b)
        b_dup = ser.from_yaml(b_yaml)
    else:
        raise ValueError(f"Invalid dup_mode: {dup_mode}")

    assert isinstance(b_dup, GTestB)

    assert b_dup == b


def test_serialize_B_cmp_objs_clean_dup(dup_mode, a1, a2, a3) -> None:
    """Test function for serialization."""

    ods = Ncm.ObjDictStr.new()
    ods.set("a1", a1)
    ods.set("a2", a2)
    ods.set("a3", a3)

    odi = Ncm.ObjDictInt.new()
    odi.set(1001, a1)
    odi.set(1002, a2)
    odi.set(1003, a3)

    oa = Ncm.ObjArray.new()
    oa.add(a1)
    oa.add(a2)
    oa.add(a3)

    b = GTestB(B_string="test", As=a1, Aa=oa, Aods=ods, Aodi=odi)  # type: ignore

    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)

    if dup_mode == "variant":
        b_dup = ser.dup_obj(b)
    elif dup_mode == "yaml":
        b_yaml = ser.to_yaml(b)
        b_dup = ser.from_yaml(b_yaml)
    else:
        raise ValueError(f"Invalid dup_mode: {dup_mode}")

    assert isinstance(b_dup, GTestB)

    assert b_dup.As is b_dup.Aods.get("a1")
    assert b_dup.As is b_dup.Aodi.get(1001)
    assert b_dup.As is b_dup.Aa.get(0)

    assert b_dup.Aods.get("a2") is b_dup.Aodi.get(1002)
    assert b_dup.Aods.get("a2") is b_dup.Aa.get(1)

    assert b_dup.Aods.get("a3") is b_dup.Aodi.get(1003)
    assert b_dup.Aods.get("a3") is b_dup.Aa.get(2)


def test_serialize_B_cmp_objs_none(dup_mode, a1, a2, a3) -> None:
    """Test function for serialization."""

    ods = Ncm.ObjDictStr.new()
    ods.set("a1", a1)
    ods.set("a2", a2)
    ods.set("a3", a3)

    odi = Ncm.ObjDictInt.new()
    odi.set(1001, a1)
    odi.set(1002, a2)
    odi.set(1003, a3)

    oa = Ncm.ObjArray.new()
    oa.add(a1)
    oa.add(a2)
    oa.add(a3)

    b = GTestB(B_string="test", As=a1, Aa=oa, Aods=ods, Aodi=odi)  # type: ignore

    ser = Ncm.Serialize.new(Ncm.SerializeOpt.NONE)

    if dup_mode == "variant":
        b_dup = ser.dup_obj(b)
    elif dup_mode == "yaml":
        b_yaml = ser.to_yaml(b)
        b_dup = ser.from_yaml(b_yaml)
    else:
        raise ValueError(f"Invalid dup_mode: {dup_mode}")

    assert isinstance(b_dup, GTestB)

    assert b_dup.As is not b_dup.Aods.get("a1")
    assert b_dup.As is not b_dup.Aodi.get(1001)
    assert b_dup.As is not b_dup.Aa.get(0)

    assert b_dup.Aods.get("a2") is not b_dup.Aodi.get(1002)
    assert b_dup.Aods.get("a2") is not b_dup.Aa.get(1)

    assert b_dup.Aods.get("a3") is not b_dup.Aodi.get(1003)
    assert b_dup.Aods.get("a3") is not b_dup.Aa.get(2)


def test_serialize_C(ser_opt, dup_mode, a1, a2, a3) -> None:
    """Test function for serialization."""

    ods1 = Ncm.ObjDictStr.new()
    ods1.set("a1", a1)
    ods1.set("a2", a2)

    ods2 = Ncm.ObjDictStr.new()
    ods2.set("a3", a3)

    odi1 = Ncm.ObjDictInt.new()
    odi1.set(1001, a1)

    odi2 = Ncm.ObjDictInt.new()
    odi2.set(1002, a2)
    odi2.set(1003, a3)

    oa1 = Ncm.ObjArray.new()
    oa1.add(a1)
    oa1.add(a2)

    oa2 = Ncm.ObjArray.new()
    oa2.add(a2)
    oa2.add(a3)

    b1 = GTestB(B_string="test", As=a1, Aa=oa1, Aods=ods1, Aodi=odi1)  # type: ignore
    b2 = GTestB(B_string="test", As=a2, Aa=oa2, Aods=ods2, Aodi=odi2)  # type: ignore

    b_ods = Ncm.ObjDictStr.new()
    b_ods.set("b1", b1)
    b_ods.set("b2", b2)

    b_odi = Ncm.ObjDictInt.new()
    b_odi.set(2001, b1)

    b_oa = Ncm.ObjArray.new()
    b_oa.add(b2)

    c = GTestC(Ca=b_oa, Cods=b_ods, Codi=b_odi)  # type: ignore

    ser = Ncm.Serialize.new(ser_opt)

    if dup_mode == "variant":
        c_dup = ser.dup_obj(c)
    elif dup_mode == "yaml":
        c_yaml = ser.to_yaml(c)
        # print(c_yaml)
        c_dup = ser.from_yaml(c_yaml)
    else:
        raise ValueError(f"Invalid dup_mode: {dup_mode}")

    assert isinstance(c_dup, GTestC)

    assert c_dup == c


def test_serialize_C_cmp_objs_none(dup_mode, a1, a2, a3) -> None:
    """Test function for serialization."""

    ods1 = Ncm.ObjDictStr.new()
    ods1.set("a1", a1)
    ods1.set("a2", a2)

    ods2 = Ncm.ObjDictStr.new()
    ods2.set("a3", a3)

    odi1 = Ncm.ObjDictInt.new()
    odi1.set(1001, a1)

    odi2 = Ncm.ObjDictInt.new()
    odi2.set(1002, a2)
    odi2.set(1003, a3)

    oa1 = Ncm.ObjArray.new()
    oa1.add(a1)
    oa1.add(a2)

    oa2 = Ncm.ObjArray.new()
    oa2.add(a2)
    oa2.add(a3)

    b1 = GTestB(B_string="test", As=a1, Aa=oa1, Aods=ods1, Aodi=odi1)  # type: ignore
    b2 = GTestB(B_string="test", As=a2, Aa=oa2, Aods=ods2, Aodi=odi2)  # type: ignore

    b_ods = Ncm.ObjDictStr.new()
    b_ods.set("b1", b1)
    b_ods.set("b2", b2)

    b_odi = Ncm.ObjDictInt.new()
    b_odi.set(2001, b1)

    b_oa = Ncm.ObjArray.new()
    b_oa.add(b2)

    c = GTestC(Ca=b_oa, Cods=b_ods, Codi=b_odi)  # type: ignore

    ser = Ncm.Serialize.new(Ncm.SerializeOpt.NONE)

    if dup_mode == "variant":
        c_dup = ser.dup_obj(c)
    elif dup_mode == "yaml":
        c_yaml = ser.to_yaml(c)
        # print(c_yaml)
        c_dup = ser.from_yaml(c_yaml)
    else:
        raise ValueError(f"Invalid dup_mode: {dup_mode}")

    assert isinstance(c_dup, GTestC)

    assert c_dup.Cods.get("b1") is not c_dup.Codi.get(2001)
    assert c_dup.Cods.get("b2") is not c_dup.Ca.get(0)


def test_serialize_C_cmp_objs_clean_dup(dup_mode, a1, a2, a3) -> None:
    """Test function for serialization."""

    ods1 = Ncm.ObjDictStr.new()
    ods1.set("a1", a1)
    ods1.set("a2", a2)

    ods2 = Ncm.ObjDictStr.new()
    ods2.set("a3", a3)

    odi1 = Ncm.ObjDictInt.new()
    odi1.set(1001, a1)

    odi2 = Ncm.ObjDictInt.new()
    odi2.set(1002, a2)
    odi2.set(1003, a3)

    oa1 = Ncm.ObjArray.new()
    oa1.add(a1)
    oa1.add(a2)

    oa2 = Ncm.ObjArray.new()
    oa2.add(a2)
    oa2.add(a3)

    b1 = GTestB(B_string="test", As=a1, Aa=oa1, Aods=ods1, Aodi=odi1)  # type: ignore
    b2 = GTestB(B_string="test", As=a2, Aa=oa2, Aods=ods2, Aodi=odi2)  # type: ignore

    b_ods = Ncm.ObjDictStr.new()
    b_ods.set("b1", b1)
    b_ods.set("b2", b2)

    b_odi = Ncm.ObjDictInt.new()
    b_odi.set(2001, b1)

    b_oa = Ncm.ObjArray.new()
    b_oa.add(b2)

    c = GTestC(Ca=b_oa, Cods=b_ods, Codi=b_odi)  # type: ignore

    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)

    if dup_mode == "variant":
        c_dup = ser.dup_obj(c)
    elif dup_mode == "yaml":
        c_yaml = ser.to_yaml(c)
        # print(c_yaml)
        c_dup = ser.from_yaml(c_yaml)
    else:
        raise ValueError(f"Invalid dup_mode: {dup_mode}")

    assert isinstance(c_dup, GTestC)

    assert c_dup.Cods.get("b1") is c_dup.Codi.get(2001)
    assert c_dup.Cods.get("b2") is c_dup.Ca.get(0)

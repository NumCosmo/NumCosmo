#!/usr/bin/env python
#
# test_py_model.py
#
# Sun Sep 08 11:27:59 2024
# Copyright  2023  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_py_model.py
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

"""Unit tests for the NcmMSet class."""

import re
import pytest
import numpy as np
from numpy.testing import assert_allclose
from numcosmo_py import Ncm, Nc, GLib, GObject

Ncm.cfg_init()


def test_ncm_model_reparam_lin_shift() -> None:
    """Test the Ncm.ModelReparamLin class."""
    rosenbrock = Ncm.ModelRosenbrock.new()

    m = Ncm.Matrix.new(2, 2)
    m.set_identity()

    v = Ncm.Vector.new(2)
    v.set(0, 1.234)
    v.set(1, 5.678)

    reparam = Ncm.ReparamLinear.new(2, m, v)
    reparam.set_param_desc_full(
        0, "mu_0", "mu_0", -1.0e2, 1.0e2, 1.0, 0.0, 0.0, Ncm.ParamType.FREE
    )
    reparam.set_param_desc_full(
        1, "mu_1", "mu_1", -1.0e2, 1.0e2, 1.0, 0.0, 0.0, Ncm.ParamType.FREE
    )

    rosenbrock.set_reparam(reparam)
    assert rosenbrock.peek_reparam() == reparam

    rosenbrock["mu_0"] = 2.2
    rosenbrock["mu_1"] = 3.3

    assert rosenbrock["mu_0"] == 2.2
    assert rosenbrock["mu_1"] == 3.3

    assert_allclose(rosenbrock.orig_param_get_by_name("x1"), 2.2 - 1.234)
    assert_allclose(rosenbrock.orig_param_get_by_name("x2"), 3.3 - 5.678)


def test_ncm_model_reparam_lin_rot() -> None:
    """Test the Ncm.ModelReparamLin class."""
    rosenbrock = Ncm.ModelRosenbrock.new()

    theta = np.pi / 5.0
    m = Ncm.Matrix.new(2, 2)
    m.set(0, 0, np.cos(theta))
    m.set(0, 1, -np.sin(theta))
    m.set(1, 0, np.sin(theta))
    m.set(1, 1, np.cos(theta))

    v = Ncm.Vector.new(2)
    v.set(0, 0.0)
    v.set(1, 0.0)

    reparam = Ncm.ReparamLinear.new(2, m, v)
    reparam.set_param_desc_full(
        0, "mu_0", "mu_0", -1.0e2, 1.0e2, 1.0, 0.0, 0.0, Ncm.ParamType.FREE
    )
    reparam.set_param_desc_full(
        1, "mu_1", "mu_1", -1.0e2, 1.0e2, 1.0, 0.0, 0.0, Ncm.ParamType.FREE
    )

    rosenbrock.set_reparam(reparam)
    assert rosenbrock.peek_reparam() == reparam

    rosenbrock["mu_0"] = 2.2
    rosenbrock["mu_1"] = 3.3

    assert rosenbrock["mu_0"] == 2.2
    assert rosenbrock["mu_1"] == 3.3

    assert_allclose(
        rosenbrock.orig_param_get_by_name("x1"),
        2.2 * np.cos(theta) + 3.3 * np.sin(theta),
    )
    assert_allclose(
        rosenbrock.orig_param_get_by_name("x2"),
        -2.2 * np.sin(theta) + 3.3 * np.cos(theta),
    )


def test_ncm_model_reparam_invalid() -> None:
    """Test the Ncm.ModelReparamLin class."""
    rosenbrock = Ncm.ModelRosenbrock.new()

    reparam = Nc.HICosmoDEReparamOk.new(6)

    with pytest.raises(
        GLib.Error,
        match=re.compile(
            rf"^ncm-model-error: ncm_model_set_reparam: model `NcmModelRosenbrock' is "
            rf"not compatible with the reparametrization `NcHICosmoDE' "
            rf"\({int(Ncm.ModelError.REPARAM_INCOMPATIBLE)}\)$",
            re.DOTALL,
        ),
    ):
        rosenbrock.set_reparam(reparam)


def test_ncm_model_reparam_param_changed() -> None:
    """Test the Ncm.ModelReparamLin class."""
    rosenbrock = Ncm.ModelRosenbrock.new()

    m = Ncm.Matrix.new(2, 2)
    m.set_identity()

    v = Ncm.Vector.new(2)
    v.set(0, 1.234)
    v.set(1, 5.678)

    reparam = Ncm.ReparamLinear.new(2, m, v)
    reparam.set_param_desc_full(
        0, "mu_0", "mu_0", -1.0e2, 1.0e2, 1.0, 0.0, 0.0, Ncm.ParamType.FREE
    )
    reparam.set_param_desc_full(
        1, "mu_1", "mu_1", -1.0e2, 1.0e2, 1.0, 0.0, 0.0, Ncm.ParamType.FREE
    )

    rosenbrock.set_reparam(reparam)
    assert rosenbrock.peek_reparam() == reparam

    with pytest.raises(
        GLib.Error,
        match=re.compile(
            rf"^ncm-model-error: ncm_model_param_index_from_name: parameter \(x1\) "
            rf"was changed by a NcmReparam, it is now named \(mu_0\). "
            rf"\({int(Ncm.ModelError.PARAM_CHANGED)}\)$",
            re.DOTALL,
        ),
    ):
        _ = rosenbrock["x1"]

    with pytest.raises(
        GLib.Error,
        match=re.compile(
            rf"^ncm-model-error: ncm_model_param_index_from_name: parameter \(x2\) "
            rf"was changed by a NcmReparam, it is now named \(mu_1\). "
            rf"\({int(Ncm.ModelError.PARAM_CHANGED)}\)$",
            re.DOTALL,
        ),
    ):
        rosenbrock["x2"] = 1.0


def test_model_id_by_type() -> None:
    """Test the Ncm.Model.id_by_type() method."""
    assert Ncm.Model.id_by_type(Nc.HICosmo) == Nc.HICosmo.id()
    assert Ncm.Model.id_by_type(Nc.HICosmoDE) == Nc.HICosmoDE.id()
    assert Ncm.Model.id_by_type(Nc.HICosmoDEXcdm) == Nc.HICosmoDEXcdm.id()

    assert Ncm.Model.id_by_type(Nc.HIReion) == Nc.HIReion.id()
    assert Ncm.Model.id_by_type(Nc.HIReionCamb) == Nc.HIReionCamb.id()


def test_model_id_by_type_invalid() -> None:
    """Test the Ncm.Model.id_by_name() method."""
    with pytest.raises(
        GLib.Error,
        match=re.compile(
            rf"^ncm-model-error: ncm_model_id_by_type: type \(GObject\) "
            rf"is not a NcmModel "
            rf"\({int(Ncm.ModelError.INVALID_TYPE)}\)$",
            re.DOTALL,
        ),
    ):
        _ = Ncm.Model.id_by_type(GObject.Object)


def test_model_setget_by_name() -> None:
    """Test the Ncm.Model.set_by_name() method."""
    rosenbrock = Ncm.ModelRosenbrock.new()

    rosenbrock.param_set_by_name("x1", 1.0)
    rosenbrock.param_set_by_name("x2", 2.0)

    assert_allclose(rosenbrock.param_get_by_name("x1"), 1.0)
    assert_allclose(rosenbrock.param_get_by_name("x2"), 2.0)


def test_model_setget_by_name_not_found() -> None:
    """Test the Ncm.Model.set_by_name() method."""
    rosenbrock = Ncm.ModelRosenbrock.new()

    with pytest.raises(
        GLib.Error,
        match=re.compile(
            rf"^ncm-model-error: ncm_model_param_get_by_name: model "
            rf"`NcmModelRosenbrock' does not have a parameter called `x3'. Use the "
            rf"method ncm_model_param_index_from_name\(\) to check if the parameter "
            rf"exists. "
            rf"\({int(Ncm.ModelError.PARAM_NAME_NOT_FOUND)}\)$",
            re.DOTALL,
        ),
    ):
        _ = rosenbrock.param_get_by_name("x3")

    with pytest.raises(
        GLib.Error,
        match=re.compile(
            rf"^ncm-model-error: ncm_model_param_set_by_name: model "
            rf"`NcmModelRosenbrock' does not have a parameter called `x3'. Use the "
            rf"method ncm_model_param_index_from_name\(\) to check if the parameter "
            rf"exists. "
            rf"\({int(Ncm.ModelError.PARAM_NAME_NOT_FOUND)}\)$",
            re.DOTALL,
        ),
    ):
        rosenbrock.param_set_by_name("x3", 1.0)


def test_model_orig_setget_by_name() -> None:
    """Test the Ncm.Model.orig_param_set_by_name() method."""
    rosenbrock = Ncm.ModelRosenbrock.new()

    rosenbrock.orig_param_set_by_name("x1", 1.0)
    rosenbrock.orig_param_set_by_name("x2", 2.0)

    assert_allclose(rosenbrock.orig_param_get_by_name("x1"), 1.0)
    assert_allclose(rosenbrock.orig_param_get_by_name("x2"), 2.0)


def test_model_orig_setget_by_name_not_found() -> None:
    """Test the Ncm.Model.orig_param_set_by_name() method."""
    rosenbrock = Ncm.ModelRosenbrock.new()

    with pytest.raises(
        GLib.Error,
        match=re.compile(
            rf"^ncm-model-error: ncm_model_orig_param_get_by_name: model "
            rf"`NcmModelRosenbrock' does not have a parameter called `x3'. "
            rf"Use the method ncm_model_orig_param_index_from_name\(\) to "
            rf"check if the parameter exists. "
            rf"\({int(Ncm.ModelError.ORIG_PARAM_NAME_NOT_FOUND)}\)$",
            re.DOTALL,
        ),
    ):
        _ = rosenbrock.orig_param_get_by_name("x3")

    with pytest.raises(
        GLib.Error,
        match=re.compile(
            rf"^ncm-model-error: ncm_model_orig_param_set_by_name: model "
            rf"`NcmModelRosenbrock' does not have a parameter called `x3'. "
            rf"Use the method ncm_model_orig_param_index_from_name\(\) to "
            rf"check if the parameter exists. "
            rf"\({int(Ncm.ModelError.ORIG_PARAM_NAME_NOT_FOUND)}\)$",
            re.DOTALL,
        ),
    ):
        rosenbrock.orig_param_set_by_name("x3", 1.0)


def test_model_setitem_getitem() -> None:
    """Test the Ncm.Model.__getitem__() and Ncm.Model.__setitem__() methods."""
    rosenbrock = Ncm.ModelRosenbrock.new()

    rosenbrock["x1"] = 1.0
    rosenbrock["x2"] = 2.0

    assert_allclose(rosenbrock["x1"], 1.0)
    assert_allclose(rosenbrock["x2"], 2.0)


def test_model_setitem_getitem_not_found() -> None:
    """Test the Ncm.Model.__getitem__() and Ncm.Model.__setitem__() methods."""
    rosenbrock = Ncm.ModelRosenbrock.new()

    with pytest.raises(
        GLib.Error,
        match=re.compile(
            rf"^ncm-model-error: Parameter named: x3 does not exist in "
            rf"NcmModelRosenbrock "
            rf"\({int(Ncm.ModelError.PARAM_NAME_NOT_FOUND)}\)$",
            re.DOTALL,
        ),
    ):
        _ = rosenbrock["x3"]

    with pytest.raises(
        GLib.Error,
        match=re.compile(
            rf"^ncm-model-error: Parameter named: x3 does not exist in "
            rf"NcmModelRosenbrock "
            rf"\({int(Ncm.ModelError.PARAM_NAME_NOT_FOUND)}\)$",
            re.DOTALL,
        ),
    ):
        rosenbrock["x3"] = 1.0


def test_model_get_desc() -> None:
    """Test the Ncm.Model.param_get_desc() method."""
    rosenbrock = Ncm.ModelRosenbrock.new()
    x1_desc = rosenbrock.param_get_desc("x1")

    assert x1_desc["name"] == rosenbrock.param_name(Ncm.ModelRosenbrockSParams.X1)
    assert x1_desc["symbol"] == rosenbrock.param_symbol(Ncm.ModelRosenbrockSParams.X1)
    assert x1_desc["lower-bound"] == rosenbrock.param_get_lower_bound(
        Ncm.ModelRosenbrockSParams.X1
    )
    assert x1_desc["upper-bound"] == rosenbrock.param_get_upper_bound(
        Ncm.ModelRosenbrockSParams.X1
    )
    assert x1_desc["scale"] == rosenbrock.param_get_scale(Ncm.ModelRosenbrockSParams.X1)
    assert x1_desc["abstol"] == rosenbrock.param_get_abstol(
        Ncm.ModelRosenbrockSParams.X1
    )
    assert (
        Ncm.ParamType.FREE if x1_desc["fit"] else Ncm.ParamType.FIXED
    ) == rosenbrock.param_get_ftype(Ncm.ModelRosenbrockSParams.X1)


def test_model_get_desc_not_found() -> None:
    """Test the Ncm.Model.param_get_desc() method."""
    rosenbrock = Ncm.ModelRosenbrock.new()

    with pytest.raises(
        GLib.Error,
        match=re.compile(
            rf"^ncm-model-error: ncm_model_param_get_desc: model `NcmModelRosenbrock' "
            rf"does not have a parameter called `x3'. "
            rf"\({int(Ncm.ModelError.PARAM_NAME_NOT_FOUND)}\)$",
            re.DOTALL,
        ),
    ):
        _ = rosenbrock.param_get_desc("x3")


def test_model_set_desc() -> None:
    """Test the Ncm.Model.param_set_desc() method."""
    rosenbrock = Ncm.ModelRosenbrock.new()

    rosenbrock.param_set_desc(
        "x1",
        {
            "lower-bound": -1.023,
            "upper-bound": 21.0,
            "scale": 155.0,
            "abstol": 1.0e-50,
            "fit": True,
            "value": 45.22,
        },
    )

    assert_allclose(rosenbrock.param_get_by_name("x1"), 45.22)
    assert_allclose(
        rosenbrock.param_get_lower_bound(Ncm.ModelRosenbrockSParams.X1), -1.023
    )
    assert_allclose(
        rosenbrock.param_get_upper_bound(Ncm.ModelRosenbrockSParams.X1), 21.0
    )
    assert_allclose(rosenbrock.param_get_scale(Ncm.ModelRosenbrockSParams.X1), 155.0)
    assert_allclose(rosenbrock.param_get_abstol(Ncm.ModelRosenbrockSParams.X1), 1.0e-50)
    assert (
        rosenbrock.param_get_ftype(Ncm.ModelRosenbrockSParams.X1) == Ncm.ParamType.FREE
    )


def test_model_set_desc_not_found() -> None:
    """Test the Ncm.Model.param_set_desc() method."""
    rosenbrock = Ncm.ModelRosenbrock.new()

    with pytest.raises(
        GLib.Error,
        match=re.compile(
            rf"^ncm-model-error: ncm_model_param_set_desc: model `NcmModelRosenbrock' "
            rf"does not have a parameter called `x3'. "
            rf"\({int(Ncm.ModelError.PARAM_NAME_NOT_FOUND)}\)$",
            re.DOTALL,
        ),
    ):
        rosenbrock.param_set_desc(
            "x3",
            {
                "lower-bound": -1.023,
                "upper-bound": 21.0,
                "scale": 155.0,
                "abstol": 1.0e-50,
                "fit": True,
                "value": 45.22,
            },
        )


def test_model_set_desc_invalid_scale() -> None:
    """Test the Ncm.Model.param_set_desc() method."""
    rosenbrock = Ncm.ModelRosenbrock.new()

    with pytest.raises(
        GLib.Error,
        match=re.compile(
            rf"^ncm-model-error: ncm_model_param_set_desc: scale must be a double. "
            rf"\({int(Ncm.ModelError.PARAM_INVALID_TYPE)}\)$",
            re.DOTALL,
        ),
    ):
        rosenbrock.param_set_desc("x1", {"scale": ("as", "234")})


def test_model_set_desc_invalid_abstol() -> None:
    """Test the Ncm.Model.param_set_desc() method."""
    rosenbrock = Ncm.ModelRosenbrock.new()

    with pytest.raises(
        GLib.Error,
        match=re.compile(
            rf"^ncm-model-error: ncm_model_param_set_desc: abstol must be a double. "
            rf"\({int(Ncm.ModelError.PARAM_INVALID_TYPE)}\)$",
            re.DOTALL,
        ),
    ):
        rosenbrock.param_set_desc("x1", {"abstol": "Not a float"})


def test_model_set_desc_invalid_fit() -> None:
    """Test the Ncm.Model.param_set_desc() method."""
    rosenbrock = Ncm.ModelRosenbrock.new()

    with pytest.raises(
        GLib.Error,
        match=re.compile(
            rf"^ncm-model-error: ncm_model_param_set_desc: fit must be a boolean. "
            rf"\({int(Ncm.ModelError.PARAM_INVALID_TYPE)}\)$",
            re.DOTALL,
        ),
    ):
        rosenbrock.param_set_desc("x1", {"fit": 1.0})


def test_model_set_desc_invalid_value() -> None:
    """Test the Ncm.Model.param_set_desc() method."""
    rosenbrock = Ncm.ModelRosenbrock.new()

    with pytest.raises(
        GLib.Error,
        match=re.compile(
            rf"^ncm-model-error: ncm_model_param_set_desc: value must be a double. "
            rf"\({int(Ncm.ModelError.PARAM_INVALID_TYPE)}\)$",
            re.DOTALL,
        ),
    ):
        rosenbrock.param_set_desc("x1", {"value": "Not a float"})


def test_model_set_desc_invalid_lower_bound() -> None:
    """Test the Ncm.Model.param_set_desc() method."""
    rosenbrock = Ncm.ModelRosenbrock.new()

    with pytest.raises(
        GLib.Error,
        match=re.compile(
            rf"^ncm-model-error: ncm_model_param_set_desc: lower-bound must be a double. "
            rf"\({int(Ncm.ModelError.PARAM_INVALID_TYPE)}\)$",
            re.DOTALL,
        ),
    ):
        rosenbrock.param_set_desc("x1", {"lower-bound": "Not a float"})


def test_model_set_desc_invalid_upper_bound() -> None:
    """Test the Ncm.Model.param_set_desc() method."""
    rosenbrock = Ncm.ModelRosenbrock.new()

    with pytest.raises(
        GLib.Error,
        match=re.compile(
            rf"^ncm-model-error: ncm_model_param_set_desc: upper-bound must be a double. "
            rf"\({int(Ncm.ModelError.PARAM_INVALID_TYPE)}\)$",
            re.DOTALL,
        ),
    ):
        rosenbrock.param_set_desc("x1", {"upper-bound": "Not a float"})


def test_model_set_desc_invalid_key() -> None:
    """Test the Ncm.Model.param_set_desc() method."""
    rosenbrock = Ncm.ModelRosenbrock.new()

    with pytest.raises(
        GLib.Error,
        match=re.compile(
            rf"^ncm-model-error: ncm_model_param_set_desc: invalid key `cow'. "
            rf"\({int(Ncm.ModelError.PARAM_INVALID_KEY)}\)$",
            re.DOTALL,
        ),
    ):
        rosenbrock.param_set_desc("x1", {"cow": 1.0})


def test_ncm_model_reparam_longer_than_model_is_resized() -> None:
    """A reparametrization written for more parameters than the model now has.

    NcmReparam carries no parameter values -- only its length, its descriptors
    and its compatible type -- so the model rebuilds it at its own length. This
    is how catalogs written before NcHICosmoDE dropped Yp still load.
    """
    cosmo = Nc.HICosmoDEXcdm.new()
    assert cosmo.len() == 7

    cosmo.set_reparam(Nc.HICosmoDEReparamOk.new(8))

    assert cosmo.len() == 7
    assert cosmo.param_name(Nc.HICosmoDESParams.OMEGA_X) == "Omegak"
    assert cosmo.peek_reparam() is not None


def test_ncm_model_reparam_shorter_than_model_is_resized() -> None:
    """Resizing is not a one-way shrink: a short reparametrization grows too."""
    cosmo = Nc.HICosmoDEXcdm.new()

    cosmo.set_reparam(Nc.HICosmoDEReparamOk.new(3))

    assert cosmo.len() == 7
    assert cosmo.param_name(Nc.HICosmoDESParams.OMEGA_X) == "Omegak"


def test_ncm_model_reparam_resized_reparametrizes_the_same_slot() -> None:
    """The rebuilt reparametrization must act, not merely attach.

    Omegak is Omegax expressed as curvature, so a resized reparametrization
    that ran old2new() reports the same value as one built at the right length.
    """
    expected = Nc.HICosmoDEXcdm.new()
    expected.set_reparam(Nc.HICosmoDEReparamOk.new(expected.len()))

    cosmo = Nc.HICosmoDEXcdm.new()
    cosmo.set_reparam(Nc.HICosmoDEReparamOk.new(8))

    assert_allclose(
        cosmo.param_get(Nc.HICosmoDESParams.OMEGA_X),
        expected.param_get(Nc.HICosmoDESParams.OMEGA_X),
    )


def test_ncm_model_reparam_without_descriptors_is_resized() -> None:
    """A reparametrization whose descriptor dictionary was dropped still fits."""
    cosmo = Nc.HICosmoDEXcdm.new()
    reparam = Nc.HICosmoDEReparamOk.new(8)
    reparam.set_property("params-desc", None)

    cosmo.set_reparam(reparam)

    assert cosmo.len() == 7


@pytest.mark.parametrize("index", [7, -1])
def test_ncm_model_reparam_descriptor_outside_model_raises(index: int) -> None:
    """A descriptor naming a slot the model does not have cannot be rebuilt.

    Descriptors are keyed by index into the model's original parameters, so a
    key outside [0, len) means the parameter it reparametrizes is gone. Index 7
    is past the end; a negative one is not a slot at all.
    """
    cosmo = Nc.HICosmoDEXcdm.new()
    reparam = Nc.HICosmoDEReparamOk.new(8)

    desc = Ncm.ObjDictInt.new()
    desc.add(
        index,
        Ncm.SParam.new("Bogus", "B", 0.0, 1.0, 0.1, 0.0, 0.5, Ncm.ParamType.FIXED),
    )
    reparam.set_property("params-desc", desc)

    with pytest.raises(
        GLib.Error,
        match=re.compile(
            rf"^ncm-model-error: ncm_model_set_reparam: reparametrization "
            rf"`NcHICosmoDEReparamOk' is for 8 parameters but model "
            rf"`NcHICosmoDEXcdm' has 7, and it reparametrizes parameter {index}, "
            rf"which no longer exists\. "
            rf"\({int(Ncm.ModelError.REPARAM_INCOMPATIBLE)}\)$",
            re.DOTALL,
        ),
    ):
        cosmo.set_reparam(reparam)


def test_ncm_model_reparam_with_extra_state_raises() -> None:
    """A subclass carrying its own state cannot be rebuilt at another length.

    Rebuilding carries over length, params-desc and compat-type -- the whole of
    a NcmReparam, but not of NcmReparamLinear, whose matrix and vector are sized
    to the old length. Refusing is the only honest answer; constructing one
    without them used to dereference NULL.
    """
    rosenbrock = Ncm.ModelRosenbrock.new()
    assert rosenbrock.len() == 2

    matrix = Ncm.Matrix.new(3, 3)
    matrix.set_identity()
    vector = Ncm.Vector.new(3)
    vector.set_all(0.0)

    reparam = Ncm.ReparamLinear.new(3, matrix, vector)
    reparam.set_compat_type("NcmModelRosenbrock")

    with pytest.raises(
        GLib.Error,
        match=re.compile(
            rf"^ncm-model-error: ncm_model_set_reparam: reparametrization "
            rf"`NcmReparamLinear' is for 3 parameters but model "
            rf"`NcmModelRosenbrock' has 2, and it carries state of its own "
            rf"\(property `(vector|matrix)'\) that no rebuild can size\. "
            rf"\({int(Ncm.ModelError.REPARAM_INCOMPATIBLE)}\)$",
            re.DOTALL,
        ),
    ):
        rosenbrock.set_reparam(reparam)


def test_ncm_model_reparam_with_extra_state_at_the_right_length_attaches() -> None:
    """The refusal above is about the length, not about NcmReparamLinear."""
    rosenbrock = Ncm.ModelRosenbrock.new()

    matrix = Ncm.Matrix.new(2, 2)
    matrix.set_identity()
    vector = Ncm.Vector.new(2)
    vector.set_all(0.0)

    reparam = Ncm.ReparamLinear.new(2, matrix, vector)
    reparam.set_compat_type("NcmModelRosenbrock")

    rosenbrock.set_reparam(reparam)

    assert rosenbrock.peek_reparam() == reparam

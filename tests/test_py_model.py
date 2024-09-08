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

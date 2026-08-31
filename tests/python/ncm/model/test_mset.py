#!/usr/bin/env python
#
# test_py_mset.py
#
# Sat Jan 20 21:16:03 2024
# Copyright  2023  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_py_mset.py
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
from numcosmo_py import Ncm, Nc, GLib

Ncm.cfg_init()


def test_mset_new_empty() -> None:
    """Test the NcmMSet constructor."""
    mset = Ncm.MSet.empty_new()
    assert mset.total_len() == 0
    assert mset.fparam_len() == 0
    assert mset.nmodels() == 0

    mset.prepare_fparam_map()

    assert mset.fparams_len() == 0


def test_mset_new_array() -> None:
    """Test the NcmMSet constructor."""
    mset = Ncm.MSet.new_array(
        [Ncm.ModelFunnel.new(5), Ncm.ModelRosenbrock.new(), Ncm.ModelMVND.new(8)]
    )
    assert mset.total_len() == 16
    assert mset.fparam_len() == 0
    assert mset.nmodels() == 3

    mset.prepare_fparam_map()

    assert mset.fparams_len() == 0


def test_mset_peek_model() -> None:
    """Test the NcmMSet peek_model function."""
    mset = Ncm.MSet.new_array(
        [Ncm.ModelFunnel.new(5), Ncm.ModelRosenbrock.new(), Ncm.ModelMVND.new(8)]
    )
    assert isinstance(mset.peek(Ncm.ModelMVND.id()), Ncm.ModelMVND)
    assert isinstance(mset.peek(Ncm.ModelFunnel.id()), Ncm.ModelFunnel)
    assert isinstance(mset.peek(Ncm.ModelRosenbrock.id()), Ncm.ModelRosenbrock)


def test_mset_fetch_model() -> None:
    """Test the NcmMSet fetch_model function."""
    mset = Ncm.MSet.new_array(
        [Ncm.ModelFunnel.new(5), Ncm.ModelRosenbrock.new(), Ncm.ModelMVND.new(8)]
    )
    assert isinstance(mset.fetch(Ncm.ModelMVND.id()), Ncm.ModelMVND)
    assert isinstance(mset.fetch(Ncm.ModelFunnel.id()), Ncm.ModelFunnel)
    assert isinstance(mset.fetch(Ncm.ModelRosenbrock.id()), Ncm.ModelRosenbrock)


def test_mset_model_peek_by_name() -> None:
    """Test the NcmMSet peek_model_by_name function."""
    mset = Ncm.MSet.new_array(
        [Ncm.ModelFunnel.new(5), Ncm.ModelRosenbrock.new(), Ncm.ModelMVND.new(8)]
    )
    assert isinstance(mset.peek_by_name("NcmModelMVND"), Ncm.ModelMVND)
    assert isinstance(mset.peek_by_name("NcmModelFunnel"), Ncm.ModelFunnel)
    assert isinstance(mset.peek_by_name("NcmModelRosenbrock"), Ncm.ModelRosenbrock)

    assert mset.peek_by_name("NcmModelMVND") == mset.peek(Ncm.ModelMVND.id())
    assert mset.peek_by_name("NcmModelFunnel") == mset.peek(Ncm.ModelFunnel.id())
    assert mset.peek_by_name("NcmModelRosenbrock") == mset.peek(
        Ncm.ModelRosenbrock.id()
    )


def test_mset_model_fetch_by_name() -> None:
    """Test the NcmMSet fetch_model_by_name function."""
    mset = Ncm.MSet.new_array(
        [Ncm.ModelFunnel.new(5), Ncm.ModelRosenbrock.new(), Ncm.ModelMVND.new(8)]
    )
    assert isinstance(mset.fetch_by_name("NcmModelMVND"), Ncm.ModelMVND)
    assert isinstance(mset.fetch_by_name("NcmModelFunnel"), Ncm.ModelFunnel)
    assert isinstance(mset.fetch_by_name("NcmModelRosenbrock"), Ncm.ModelRosenbrock)

    assert mset.fetch_by_name("NcmModelMVND") == mset.fetch(Ncm.ModelMVND.id())
    assert mset.fetch_by_name("NcmModelFunnel") == mset.fetch(Ncm.ModelFunnel.id())
    assert mset.fetch_by_name("NcmModelRosenbrock") == mset.fetch(
        Ncm.ModelRosenbrock.id()
    )


def test_mset_peek_by_name_with_stackpos() -> None:
    """Test the NcmMSet peek_by_name function."""
    mset = Ncm.MSet.new_array(
        [Ncm.ModelFunnel.new(5), Ncm.ModelRosenbrock.new(), Ncm.ModelMVND.new(8)]
    )

    assert isinstance(mset.peek_by_name("NcmModelMVND:0"), Ncm.ModelMVND)
    assert isinstance(mset.peek_by_name("NcmModelFunnel:0"), Ncm.ModelFunnel)
    assert isinstance(mset.peek_by_name("NcmModelRosenbrock:0"), Ncm.ModelRosenbrock)

    assert mset.peek_by_name("NcmModelMVND:0") == mset.peek(Ncm.ModelMVND.id())
    assert mset.peek_by_name("NcmModelFunnel:0") == mset.peek(Ncm.ModelFunnel.id())
    assert mset.peek_by_name("NcmModelRosenbrock:0") == mset.peek(
        Ncm.ModelRosenbrock.id()
    )


def test_mset_fetch_model_not_set() -> None:
    """Test the NcmMSet fetch_model function."""
    mset = Ncm.MSet.new_array([Ncm.ModelFunnel.new(5), Ncm.ModelRosenbrock.new()])
    with pytest.raises(
        GLib.Error,
        match=re.compile(
            rf"^ncm-mset-error: ncm_mset_fetch: model with id -1 not found. "
            rf"\({int(Ncm.MSetError.MODEL_NOT_SET)}\)$",
            re.DOTALL,
        ),
    ):
        _ = mset.fetch(-1)


def test_mset_fetch_by_name_not_set() -> None:
    """Test the NcmMSet fetch_by_name function."""
    mset = Ncm.MSet.new_array([Ncm.ModelFunnel.new(5), Ncm.ModelRosenbrock.new()])
    with pytest.raises(
        GLib.Error,
        match=re.compile(
            rf"^ncm-mset-error: ncm_mset_fetch_by_name: model with name `NcmModelMVND' "
            rf"not found. \({int(Ncm.MSetError.MODEL_NOT_SET)}\)$",
            re.DOTALL,
        ),
    ):
        _ = mset.fetch_by_name("NcmModelMVND")


def test_mset_fparam_map() -> None:
    """Test the NcmMSet free parameter map."""
    mset = Ncm.MSet.new_array(
        [Ncm.ModelFunnel.new(5), Ncm.ModelRosenbrock.new(), Ncm.ModelMVND.new(8)]
    )

    mset.peek(Ncm.ModelMVND.id()).param_set_ftype(0, Ncm.ParamType.FREE)
    mset.peek(Ncm.ModelMVND.id()).param_set_ftype(1, Ncm.ParamType.FREE)
    mset.peek(Ncm.ModelMVND.id()).param_set_ftype(2, Ncm.ParamType.FREE)

    mset.peek(Ncm.ModelRosenbrock.id()).param_set_ftype(0, Ncm.ParamType.FREE)

    mset.peek(Ncm.ModelFunnel.id()).param_set_ftype(3, Ncm.ParamType.FREE)
    mset.peek(Ncm.ModelFunnel.id()).param_set_ftype(4, Ncm.ParamType.FREE)

    mset.prepare_fparam_map()

    assert mset.total_len() == 16
    assert mset.fparams_len() == 6

    assert mset.fparam_get_pi(0).mid == Ncm.ModelMVND.id()
    assert mset.fparam_get_pi(0).pid == 0

    assert mset.fparam_get_pi(1).mid == Ncm.ModelMVND.id()
    assert mset.fparam_get_pi(1).pid == 1

    assert mset.fparam_get_pi(2).mid == Ncm.ModelMVND.id()
    assert mset.fparam_get_pi(2).pid == 2

    assert mset.fparam_get_pi(3).mid == Ncm.ModelRosenbrock.id()
    assert mset.fparam_get_pi(3).pid == 0

    assert mset.fparam_get_pi(4).mid == Ncm.ModelFunnel.id()
    assert mset.fparam_get_pi(4).pid == 3

    assert mset.fparam_get_pi(5).mid == Ncm.ModelFunnel.id()
    assert mset.fparam_get_pi(5).pid == 4

    assert mset.fparam_get_fpi(Ncm.ModelMVND.id(), 0) == 0
    assert mset.fparam_get_fpi(Ncm.ModelMVND.id(), 1) == 1
    assert mset.fparam_get_fpi(Ncm.ModelMVND.id(), 2) == 2
    assert mset.fparam_get_fpi(Ncm.ModelMVND.id(), 3) == -1

    assert mset.fparam_get_fpi(Ncm.ModelRosenbrock.id(), 0) == 3
    assert mset.fparam_get_fpi(Ncm.ModelRosenbrock.id(), 1) == -1

    assert mset.fparam_get_fpi(Ncm.ModelFunnel.id(), 0) == -1
    assert mset.fparam_get_fpi(Ncm.ModelFunnel.id(), 1) == -1
    assert mset.fparam_get_fpi(Ncm.ModelFunnel.id(), 2) == -1
    assert mset.fparam_get_fpi(Ncm.ModelFunnel.id(), 3) == 4
    assert mset.fparam_get_fpi(Ncm.ModelFunnel.id(), 4) == 5


def test_mset_fparam_map_change() -> None:
    """Test the NcmMSet free parameter map."""
    mset = Ncm.MSet.new_array(
        [Ncm.ModelFunnel.new(5), Ncm.ModelRosenbrock.new(), Ncm.ModelMVND.new(8)]
    )

    mset.peek(Ncm.ModelMVND.id()).param_set_ftype(0, Ncm.ParamType.FREE)
    mset.peek(Ncm.ModelMVND.id()).param_set_ftype(1, Ncm.ParamType.FREE)
    mset.peek(Ncm.ModelMVND.id()).param_set_ftype(2, Ncm.ParamType.FREE)

    mset.peek(Ncm.ModelRosenbrock.id()).param_set_ftype(0, Ncm.ParamType.FREE)

    mset.peek(Ncm.ModelFunnel.id()).param_set_ftype(3, Ncm.ParamType.FREE)
    mset.peek(Ncm.ModelFunnel.id()).param_set_ftype(4, Ncm.ParamType.FREE)

    mset.prepare_fparam_map()

    assert mset.total_len() == 16
    assert mset.fparams_len() == 6

    assert mset.fparam_get_pi(0).mid == Ncm.ModelMVND.id()
    assert mset.fparam_get_pi(0).pid == 0

    assert mset.fparam_get_pi(1).mid == Ncm.ModelMVND.id()
    assert mset.fparam_get_pi(1).pid == 1

    assert mset.fparam_get_pi(2).mid == Ncm.ModelMVND.id()
    assert mset.fparam_get_pi(2).pid == 2

    assert mset.fparam_get_pi(3).mid == Ncm.ModelRosenbrock.id()
    assert mset.fparam_get_pi(3).pid == 0

    assert mset.fparam_get_pi(4).mid == Ncm.ModelFunnel.id()
    assert mset.fparam_get_pi(4).pid == 3

    assert mset.fparam_get_pi(5).mid == Ncm.ModelFunnel.id()
    assert mset.fparam_get_pi(5).pid == 4

    assert mset.fparam_get_fpi(Ncm.ModelMVND.id(), 0) == 0
    assert mset.fparam_get_fpi(Ncm.ModelMVND.id(), 1) == 1
    assert mset.fparam_get_fpi(Ncm.ModelMVND.id(), 2) == 2
    assert mset.fparam_get_fpi(Ncm.ModelMVND.id(), 3) == -1

    assert mset.fparam_get_fpi(Ncm.ModelRosenbrock.id(), 0) == 3
    assert mset.fparam_get_fpi(Ncm.ModelRosenbrock.id(), 1) == -1

    assert mset.fparam_get_fpi(Ncm.ModelFunnel.id(), 0) == -1
    assert mset.fparam_get_fpi(Ncm.ModelFunnel.id(), 1) == -1
    assert mset.fparam_get_fpi(Ncm.ModelFunnel.id(), 2) == -1
    assert mset.fparam_get_fpi(Ncm.ModelFunnel.id(), 3) == 4
    assert mset.fparam_get_fpi(Ncm.ModelFunnel.id(), 4) == 5

    mset.peek(Ncm.ModelMVND.id()).param_set_ftype(2, Ncm.ParamType.FIXED)
    mset.prepare_fparam_map()

    assert mset.total_len() == 16
    assert mset.fparams_len() == 5

    assert mset.fparam_get_pi(0).mid == Ncm.ModelMVND.id()
    assert mset.fparam_get_pi(0).pid == 0

    assert mset.fparam_get_pi(1).mid == Ncm.ModelMVND.id()
    assert mset.fparam_get_pi(1).pid == 1

    assert mset.fparam_get_pi(2).mid == Ncm.ModelRosenbrock.id()
    assert mset.fparam_get_pi(2).pid == 0

    assert mset.fparam_get_pi(3).mid == Ncm.ModelFunnel.id()
    assert mset.fparam_get_pi(3).pid == 3

    assert mset.fparam_get_pi(4).mid == Ncm.ModelFunnel.id()
    assert mset.fparam_get_pi(4).pid == 4

    assert mset.fparam_get_fpi(Ncm.ModelMVND.id(), 0) == 0
    assert mset.fparam_get_fpi(Ncm.ModelMVND.id(), 1) == 1
    assert mset.fparam_get_fpi(Ncm.ModelMVND.id(), 2) == -1
    assert mset.fparam_get_fpi(Ncm.ModelMVND.id(), 3) == -1

    assert mset.fparam_get_fpi(Ncm.ModelRosenbrock.id(), 0) == 2
    assert mset.fparam_get_fpi(Ncm.ModelRosenbrock.id(), 1) == -1

    assert mset.fparam_get_fpi(Ncm.ModelFunnel.id(), 0) == -1
    assert mset.fparam_get_fpi(Ncm.ModelFunnel.id(), 1) == -1
    assert mset.fparam_get_fpi(Ncm.ModelFunnel.id(), 2) == -1
    assert mset.fparam_get_fpi(Ncm.ModelFunnel.id(), 3) == 3
    assert mset.fparam_get_fpi(Ncm.ModelFunnel.id(), 4) == 4


def test_mset_fparam_get_pi_by_name() -> None:
    """Test the NcmMSet free parameter map."""
    mset = Ncm.MSet.new_array(
        [Ncm.ModelFunnel.new(5), Ncm.ModelRosenbrock.new(), Ncm.ModelMVND.new(8)]
    )

    mvnd_id = Ncm.ModelMVND.id()
    funnel_id = Ncm.ModelFunnel.id()
    rosenbrock_id = Ncm.ModelRosenbrock.id()

    mset.peek(mvnd_id).param_set_ftype(0, Ncm.ParamType.FREE)
    mset.peek(mvnd_id).param_set_ftype(1, Ncm.ParamType.FREE)
    mset.peek(mvnd_id).param_set_ftype(2, Ncm.ParamType.FREE)

    mset.peek(rosenbrock_id).param_set_ftype(0, Ncm.ParamType.FREE)

    mset.peek(funnel_id).param_set_ftype(0, Ncm.ParamType.FREE)
    mset.peek(funnel_id).param_set_ftype(3, Ncm.ParamType.FREE)
    mset.peek(funnel_id).param_set_ftype(4, Ncm.ParamType.FREE)

    mset.prepare_fparam_map()

    assert mset.total_len() == 16
    assert mset.fparams_len() == 7

    assert mset.fparam_get_pi_by_name("NcmModelMVND:mu_0").mid == mvnd_id
    assert mset.fparam_get_pi_by_name("NcmModelMVND:mu_0").pid == 0

    assert mset.fparam_get_pi_by_name("NcmModelMVND:mu_1").mid == mvnd_id
    assert mset.fparam_get_pi_by_name("NcmModelMVND:mu_1").pid == 1

    assert mset.fparam_get_pi_by_name("NcmModelMVND:mu_2").mid == mvnd_id
    assert mset.fparam_get_pi_by_name("NcmModelMVND:mu_2").pid == 2

    assert mset.fparam_get_pi_by_name("NcmModelRosenbrock:x1").mid == rosenbrock_id
    assert mset.fparam_get_pi_by_name("NcmModelRosenbrock:x1").pid == 0

    assert mset.fparam_get_pi_by_name("NcmModelFunnel:nu").mid == funnel_id
    assert mset.fparam_get_pi_by_name("NcmModelFunnel:nu").pid == 0

    assert mset.fparam_get_pi_by_name("NcmModelFunnel:x_3").mid == funnel_id
    assert mset.fparam_get_pi_by_name("NcmModelFunnel:x_3").pid == 4

    assert mset.fparam_get_pi_by_name("NcmModelFunnel:Idontexist") is None
    assert mset.fparam_get_pi_by_name("NcmModelIdontexist:not") is None


def test_mset_fparam_get_pi_by_param_name_only() -> None:
    """Test the NcmMSet free parameter map."""
    mset = Ncm.MSet.new_array(
        [Ncm.ModelFunnel.new(5), Ncm.ModelRosenbrock.new(), Ncm.ModelMVND.new(8)]
    )

    mvnd_id = Ncm.ModelMVND.id()
    funnel_id = Ncm.ModelFunnel.id()
    rosenbrock_id = Ncm.ModelRosenbrock.id()

    mset.peek(mvnd_id).param_set_ftype(0, Ncm.ParamType.FREE)
    mset.peek(mvnd_id).param_set_ftype(1, Ncm.ParamType.FREE)
    mset.peek(mvnd_id).param_set_ftype(2, Ncm.ParamType.FREE)

    mset.peek(rosenbrock_id).param_set_ftype(0, Ncm.ParamType.FREE)

    mset.peek(funnel_id).param_set_ftype(0, Ncm.ParamType.FREE)
    mset.peek(funnel_id).param_set_ftype(3, Ncm.ParamType.FREE)
    mset.peek(funnel_id).param_set_ftype(4, Ncm.ParamType.FREE)

    mset.prepare_fparam_map()

    assert mset.total_len() == 16
    assert mset.fparams_len() == 7

    assert mset.fparam_get_pi_by_name("mu_0").mid == mvnd_id
    assert mset.fparam_get_pi_by_name("mu_0").pid == 0

    assert mset.fparam_get_pi_by_name("mu_1").mid == mvnd_id
    assert mset.fparam_get_pi_by_name("mu_1").pid == 1

    assert mset.fparam_get_pi_by_name("mu_2").mid == mvnd_id
    assert mset.fparam_get_pi_by_name("mu_2").pid == 2

    assert mset.fparam_get_pi_by_name("x1").mid == rosenbrock_id
    assert mset.fparam_get_pi_by_name("x1").pid == 0

    assert mset.fparam_get_pi_by_name("nu").mid == funnel_id
    assert mset.fparam_get_pi_by_name("nu").pid == 0

    assert mset.fparam_get_pi_by_name("x_3").mid == funnel_id
    assert mset.fparam_get_pi_by_name("x_3").pid == 4

    assert mset.fparam_get_pi_by_name("Idontexist") is None
    assert mset.fparam_get_pi_by_name("not") is None


def test_mset_fparam_get_pi_by_param_name_amiguity() -> None:
    """Test the NcmMSet free parameter map."""
    mset = Ncm.MSet.new_array(
        [Ncm.ModelFunnel.new(5), Ncm.ModelRosenbrock.new(), Ncm.ModelMVND.new(8)]
    )

    mvnd_id = Ncm.ModelMVND.id()
    funnel_id = Ncm.ModelFunnel.id()
    rosenbrock_id = Ncm.ModelRosenbrock.id()

    rosenbrock = mset.peek(rosenbrock_id)

    m = Ncm.Matrix.new(2, 2)
    m.set_identity()

    v = Ncm.Vector.new(2)
    v.set_all(0.0)

    reparam = Ncm.ReparamLinear.new(2, m, v)
    reparam.set_param_desc_full(
        0, "mu_0", "mu_0", -1.0, 1.0, 1.0, 0.0, 0.0, Ncm.ParamType.FREE
    )
    reparam.set_param_desc_full(
        1, "mu_1", "mu_1", -1.0, 1.0, 1.0, 0.0, 0.0, Ncm.ParamType.FREE
    )

    rosenbrock.set_reparam(reparam)
    mset.peek(mvnd_id).param_set_ftype(0, Ncm.ParamType.FREE)
    mset.peek(mvnd_id).param_set_ftype(1, Ncm.ParamType.FREE)
    mset.peek(mvnd_id).param_set_ftype(2, Ncm.ParamType.FREE)

    mset.peek(rosenbrock_id).param_set_ftype(0, Ncm.ParamType.FREE)

    mset.peek(funnel_id).param_set_ftype(0, Ncm.ParamType.FREE)
    mset.peek(funnel_id).param_set_ftype(3, Ncm.ParamType.FREE)
    mset.peek(funnel_id).param_set_ftype(4, Ncm.ParamType.FREE)

    mset.prepare_fparam_map()

    assert mset.total_len() == 16
    assert mset.fparams_len() == 7

    with pytest.raises(
        GLib.Error,
        match=re.compile(
            rf"^ncm-mset-error: ncm_mset_fparam_get_pi_by_name: more than one \[2\] "
            rf"parameters with the same name mu_0, use the full name to avoid "
            rf"ambiguities. "
            rf"\({int(Ncm.MSetError.PARAM_NAME_AMBIGUOUS)}\)$",
            re.DOTALL,
        ),
    ):
        _ = mset.fparam_get_pi_by_name("mu_0")


def test_mset_fparam_get_pi_by_param_full_stacking() -> None:
    """Test the NcmMSet free parameter map."""
    mset = Ncm.MSet.new_array(
        [
            Nc.ClusterPhotozGaussGlobal.new(0.0, 1.0, 1.234, 432.1),
            Nc.ClusterPhotozGaussGlobal.new(0.1, 1.1, 1.235, 432.2),
        ]
    )

    mset.peek_pos(Nc.ClusterPhotozGaussGlobal.id(), 0).param_set_ftype(
        Nc.ClusterPhotozGaussGlobalSParams.Z_BIAS, Ncm.ParamType.FREE
    )
    mset.peek_pos(Nc.ClusterPhotozGaussGlobal.id(), 1).param_set_ftype(
        Nc.ClusterPhotozGaussGlobalSParams.Z_BIAS, Ncm.ParamType.FREE
    )

    mset.prepare_fparam_map()

    assert (
        mset.fparam_get_pi_by_name("NcClusterRedshift:z-bias").mid
        == Nc.ClusterPhotozGaussGlobal.id()
    )
    assert (
        mset.fparam_get_pi_by_name("NcClusterRedshift:01:z-bias").pid
        == Nc.ClusterPhotozGaussGlobalSParams.Z_BIAS
    )


def test_mset_split_full_name_invalid_stackpos() -> None:
    """Test exception handling in NcmMSet.split_full_name."""
    with pytest.raises(
        GLib.Error,
        match=re.escape(
            rf"ncm-mset-error: ncm_mset_param_split_full_name: "
            rf"invalid stackpos number (132344 >= 1000). "
            rf"({int(Ncm.MSetError.FULLNAME_INVALID)})"
        ),
    ):
        _ = Ncm.MSet.split_full_name("NcHICosmo:132344:w")


def test_mset_split_full_name_invalid_name() -> None:
    """Test exception handling in NcmMSet.split_full_name."""
    found, _, _, _ = Ncm.MSet.split_full_name("ncHICosmo:132:w")
    assert not found

    found, _, _, _ = Ncm.MSet.split_full_name("NcHICosmo:1a32:w")
    assert not found

    found, _, _, _ = Ncm.MSet.split_full_name("NcHICosmo:132:w:")


def test_mset_peek_by_name_invalid_stackpos() -> None:
    """Test the NcmMSet peek_by_name function."""
    mset = Ncm.MSet.new_array(
        [Ncm.ModelFunnel.new(5), Ncm.ModelRosenbrock.new(), Ncm.ModelMVND.new(8)]
    )

    with pytest.raises(
        GLib.Error,
        match=re.escape(
            rf"ncm-mset-error: ncm_mset_peek_by_name: "
            rf"invalid stackpos number (sadd). ({int(Ncm.MSetError.NAMESPACE_INVALID)})"
        ),
    ):
        _ = mset.peek_by_name("NcmModelMVND:sadd")


def test_mset_set_pos_with_submodel() -> None:
    """Test the NcmMSet set_pos with submodel."""
    mset = Ncm.MSet.new_array(
        [Ncm.ModelFunnel.new(5), Ncm.ModelRosenbrock.new(), Ncm.ModelMVND.new(8)]
    )
    hireion = Nc.HIReionCamb.new()

    with pytest.raises(
        GLib.Error,
        match=re.compile(
            rf"^ncm-mset-error: ncm_mset_set_pos: cannot add "
            rf"model.*\({int(Ncm.MSetError.SUBMODEL)}\)$",
            re.DOTALL,
        ),
    ):
        mset.set_pos(hireion, 10)


def test_mset_pos_with_non_stackable_model() -> None:
    """Test the NcmMSet pos with non-stackable model."""
    mset = Ncm.MSet.new_array(
        [Ncm.ModelFunnel.new(5), Ncm.ModelRosenbrock.new(), Ncm.ModelMVND.new(8)]
    )
    cosmo = Nc.HICosmoLCDM.new()

    with pytest.raises(
        GLib.Error,
        match=re.compile(
            rf"^ncm-mset-error: ncm_mset_set_pos: cannot stack object in position "
            rf"10 NcmMSet, type `NcHICosmoLCDM' is not "
            rf"stackable.*\({int(Ncm.MSetError.MODEL_NOT_STACKABLE)}\)$",
            re.DOTALL,
        ),
    ):
        mset.set_pos(cosmo, 10)


def test_mset_set_fmap_invalid_parameter() -> None:
    """Test the NcmMSet set_fmap with invalid parameter."""
    mset = Ncm.MSet.new_array(
        [Ncm.ModelFunnel.new(5), Ncm.ModelRosenbrock.new(), Ncm.ModelMVND.new(8)]
    )

    with pytest.raises(
        GLib.Error,
        match=re.compile(
            rf"^ncm-mset-error: ncm_mset_set_fmap: cannot set fmap, invalid param "
            rf"`NcHICosmo:ble'..*\({int(Ncm.MSetError.FULLNAME_NOT_FOUND)}\)$",
            re.DOTALL,
        ),
    ):
        mset.set_fmap(["NcHICosmo:ble"], True)


def test_mset_param_get_ftype_invalid_mid() -> None:
    """Test the NcmMSet get_ftype with invalid mid."""
    mset = Ncm.MSet.new_array(
        [Ncm.ModelFunnel.new(5), Ncm.ModelRosenbrock.new(), Ncm.ModelMVND.new(8)]
    )

    with pytest.raises(
        GLib.Error,
        match=re.compile(
            rf"^ncm-mset-error: ncm_mset_param_get_ftype: cannot get ftype of "
            rf"model-id -1, model not set.*"
            rf"\({int(Ncm.MSetError.MODEL_NOT_SET)}\)$",
            re.DOTALL,
        ),
    ):
        mset.param_get_ftype(-1, 0)


def test_mset_param_get_by_full_name() -> None:
    """Test the NcmMSet param_get_by_full_name."""
    mset = Ncm.MSet.new_array(
        [Ncm.ModelFunnel.new(5), Ncm.ModelRosenbrock.new(), Ncm.ModelMVND.new(8)]
    )

    pi = mset.param_get_by_full_name("NcmModelFunnel:nu")
    assert pi.mid == Ncm.ModelFunnel.id()
    assert pi.pid == Ncm.ModelFunnelSParams.NU

    pi = mset.param_get_by_full_name("NcmModelFunnel:0:nu")
    assert pi.mid == Ncm.ModelFunnel.id()
    assert pi.pid == Ncm.ModelFunnelSParams.NU

    pi = mset.param_get_by_full_name(
        f"NcmModelFunnel:0:{int(Ncm.ModelFunnelSParams.NU)}"
    )
    assert pi.mid == Ncm.ModelFunnel.id()
    assert pi.pid == Ncm.ModelFunnelSParams.NU


def test_mset_param_get_by_full_name_invalid_ns() -> None:
    """Test the NcmMSet param_get_by_full_name."""
    mset = Ncm.MSet.new_array(
        [Ncm.ModelFunnel.new(5), Ncm.ModelRosenbrock.new(), Ncm.ModelMVND.new(8)]
    )

    with pytest.raises(
        GLib.Error,
        match=re.compile(
            rf"^ncm-mset-error: ncm_mset_param_get_by_full_name: namespace "
            rf"`NcmModelBla' not found. "
            rf"\({int(Ncm.MSetError.NAMESPACE_NOT_FOUND)}\)$",
            re.DOTALL,
        ),
    ):
        _ = mset.param_get_by_full_name("NcmModelBla:ble")


def test_mset_param_get_by_full_name_not_found() -> None:
    """Test the NcmMSet param_get_by_full_name."""
    mset = Ncm.MSet.new_array(
        [Ncm.ModelFunnel.new(5), Ncm.ModelRosenbrock.new(), Ncm.ModelMVND.new(8)]
    )

    assert mset.param_get_by_full_name("NcmModelFunnel:ble") is None


def test_mset_param_get_by_full_name_invalid() -> None:
    """Test the NcmMSet param_get_by_full_name."""
    mset = Ncm.MSet.new_array(
        [Ncm.ModelFunnel.new(5), Ncm.ModelRosenbrock.new(), Ncm.ModelMVND.new(8)]
    )

    with pytest.raises(
        GLib.Error,
        match=re.compile(
            rf"^ncm-mset-error: ncm_mset_param_get_by_full_name: invalid "
            rf"full name `NcmModelFunnelbleble'. "
            rf"\({int(Ncm.MSetError.FULLNAME_INVALID)}\)$",
            re.DOTALL,
        ),
    ):
        _ = mset.param_get_by_full_name("NcmModelFunnelbleble")


def test_mset_setitem() -> None:
    """Test the NcmMSet setitem."""
    mset = Ncm.MSet.new_array(
        [Ncm.ModelFunnel.new(5), Ncm.ModelRosenbrock.new(), Ncm.ModelMVND.new(8)]
    )

    mset[Ncm.ModelFunnel.id()] = Ncm.ModelFunnel.new(6)
    mset["NcmModelRosenbrock"] = Ncm.ModelRosenbrock.new()
    mset["NcmModelMVND"] = Ncm.ModelMVND.new(9)


def test_mset_getitem() -> None:
    """Test the NcmMSet getitem."""
    mset = Ncm.MSet.new_array(
        [Ncm.ModelFunnel.new(5), Ncm.ModelRosenbrock.new(), Ncm.ModelMVND.new(8)]
    )

    assert isinstance(mset[Ncm.ModelFunnel.id()], Ncm.ModelFunnel)
    assert isinstance(mset["NcmModelRosenbrock"], Ncm.ModelRosenbrock)
    assert isinstance(mset["NcmModelMVND"], Ncm.ModelMVND)


def test_mset_setitem_invalid_mid() -> None:
    """Test the NcmMSet setitem."""
    mset = Ncm.MSet.new_array(
        [Ncm.ModelFunnel.new(5), Ncm.ModelRosenbrock.new(), Ncm.ModelMVND.new(8)]
    )

    with pytest.raises(
        GLib.Error,
        match=re.compile(
            rf"^ncm-mset-error: ncm_mset___setitem__: model id mismatch, "
            rf"expected -1, got 2000. "
            rf"\({int(Ncm.MSetError.MODEL_ID_MISMATCH)}\)$",
            re.DOTALL,
        ),
    ):
        mset[-1] = Ncm.ModelFunnel.new(6)


def test_mset_setitem_mid_mismatch() -> None:
    """Test the NcmMSet setitem."""
    mset = Ncm.MSet.new_array(
        [Ncm.ModelFunnel.new(5), Ncm.ModelRosenbrock.new(), Ncm.ModelMVND.new(8)]
    )

    with pytest.raises(
        GLib.Error,
        match=re.compile(
            rf"^ncm-mset-error: ncm_mset___setitem__: model id mismatch, "
            rf"expected {Ncm.ModelRosenbrock.id()}, got {Ncm.ModelFunnel.id()}. "
            rf"\({int(Ncm.MSetError.MODEL_ID_MISMATCH)}\)$",
            re.DOTALL,
        ),
    ):
        mset[Ncm.ModelRosenbrock.id()] = Ncm.ModelFunnel.new(6)


def test_mset_setitem_invalid_arg() -> None:
    """Test the NcmMSet setitem."""
    mset = Ncm.MSet.new_array(
        [Ncm.ModelFunnel.new(5), Ncm.ModelRosenbrock.new(), Ncm.ModelMVND.new(8)]
    )

    with pytest.raises(
        GLib.Error,
        match=re.compile(
            rf"^ncm-mset-error: ncm_mset___setitem__: invalid argument type. "
            rf"\({int(Ncm.MSetError.MODEL_INVALID_ID)}\)$",
            re.DOTALL,
        ),
    ):
        mset[(1, 2, 4)] = Ncm.ModelFunnel.new(6)


def test_mset_setitem_invalid_ns() -> None:
    """Test the NcmMSet setitem."""
    mset = Ncm.MSet.new_array(
        [Ncm.ModelFunnel.new(5), Ncm.ModelRosenbrock.new(), Ncm.ModelMVND.new(8)]
    )

    with pytest.raises(
        GLib.Error,
        match=re.compile(
            rf"^ncm-mset-error: ncm_mset___setitem__: namespace "
            rf"`NcmModelBla' not found. "
            rf"\({int(Ncm.MSetError.NAMESPACE_NOT_FOUND)}\)$",
            re.DOTALL,
        ),
    ):
        mset["NcmModelBla"] = Ncm.ModelFunnel.new(6)


def test_mset_setitem_invalid_ns_stackpos() -> None:
    """Test the NcmMSet setitem."""
    mset = Ncm.MSet.new_array(
        [Ncm.ModelFunnel.new(5), Ncm.ModelRosenbrock.new(), Ncm.ModelMVND.new(8)]
    )

    with pytest.raises(
        GLib.Error,
        match=re.compile(
            rf"^ncm-mset-error: ncm_mset___setitem__: invalid namespace "
            rf"`NcmModelBla:zaa2'. "
            rf"\({int(Ncm.MSetError.NAMESPACE_INVALID)}\)$",
            re.DOTALL,
        ),
    ):
        mset["NcmModelBla:zaa2"] = Ncm.ModelFunnel.new(6)


def test_mset_getitem_invalid_ns() -> None:
    """Test the NcmMSet getitem."""
    mset = Ncm.MSet.new_array(
        [Ncm.ModelFunnel.new(5), Ncm.ModelRosenbrock.new(), Ncm.ModelMVND.new(8)]
    )

    with pytest.raises(
        GLib.Error,
        match=re.compile(
            rf"^ncm-mset-error: ncm_mset___getitem__: ncm_mset_fetch_by_name: "
            rf"model with name `NcmModelBla' not found. "
            rf"\({int(Ncm.MSetError.MODEL_NOT_SET)}\)$",
            re.DOTALL,
        ),
    ):
        _ = mset["NcmModelBla"]


def test_mset_getitem_invalid_ns_stackpos() -> None:
    """Test the NcmMSet getitem."""
    mset = Ncm.MSet.new_array(
        [Ncm.ModelFunnel.new(5), Ncm.ModelRosenbrock.new(), Ncm.ModelMVND.new(8)]
    )

    with pytest.raises(
        GLib.Error,
        match=re.compile(
            rf"^ncm-mset-error: ncm_mset___getitem__: ncm_mset_fetch_by_name: "
            rf"ncm_mset_peek_by_name: invalid stackpos number \(zaa2\). "
            rf"\({int(Ncm.MSetError.NAMESPACE_INVALID)}\)$",
            re.DOTALL,
        ),
    ):
        _ = mset["NcmModelBla:zaa2"]


def test_mset_save_no_comments(tmp_path):
    """Test the NcmMSet save function."""
    mset = Ncm.MSet.new_array(
        [
            Ncm.ModelRosenbrock.new(),
        ]
    )

    tmp_file = tmp_path / "simple_keyfile.ini"
    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
    mset.save(ser, tmp_file.as_posix(), False)

    assert tmp_file.exists()

    mset2 = Ncm.MSet.load(tmp_file.as_posix(), ser)
    assert mset2.total_len() == mset.total_len()


def test_mset_save_with_comments(tmp_path):
    """Test the NcmMSet save function."""
    mset = Ncm.MSet.new_array(
        [
            Ncm.ModelRosenbrock.new(),
        ]
    )

    tmp_file = tmp_path / "simple_keyfile.ini"
    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
    mset.save(ser, tmp_file.as_posix(), True)

    assert tmp_file.exists()

    mset2 = Ncm.MSet.load(tmp_file.as_posix(), ser)
    assert mset2.total_len() == mset.total_len()


def test_mset_load_empty(tmp_path):
    """Test the NcmMSet load function."""
    keyfile = GLib.KeyFile.new()

    tmp_file = tmp_path / "simple_keyfile.ini"
    keyfile.save_to_file(tmp_file.as_posix())

    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
    mset = Ncm.MSet.load(tmp_file.as_posix(), ser)

    assert mset is not None
    assert isinstance(mset, Ncm.MSet)


def test_mset_load_not_found():
    """Test the NcmMSet load function."""
    with pytest.raises(
        GLib.Error,
        match=re.compile(
            r"^g-file-error-quark: ncm_mset_load: Invalid mset configuration file: "
            r"file_not_found.ini: : No such file or directory \(\d*\)$",
            re.DOTALL,
        ),
    ):
        ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
        _ = Ncm.MSet.load("file_not_found.ini", ser)


def _save_lcdm_mset(tmp_path, name="cosmo.mset"):
    """Save an mset holding a plain NcHICosmoLCDM (carries its default NcBBN)."""
    mset = Ncm.MSet.new_array([Nc.HICosmoLCDM.new()])
    tmp_file = tmp_path / name
    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
    mset.save(ser, tmp_file.as_posix(), False)
    return tmp_file


def _mset_file_groups(path):
    """Split an mset keyfile into a header and a {group name: group text} dict."""
    text = path.read_text()
    parts = re.split(r"(?m)^(?=\[)", text)
    header = parts[0]
    groups = {p.split("]")[0][1:]: p for p in parts[1:]}
    return header, groups


def test_mset_load_duplicate_submodel_group(tmp_path):
    """A submodel type appearing in two groups is a clean error."""
    tmp_file = _save_lcdm_mset(tmp_path)
    header, groups = _mset_file_groups(tmp_file)
    dup = groups["NcBBN"].replace("[NcBBN]", "[NcBBN:1]")
    tmp_file.write_text(header + "".join(groups.values()) + dup)

    with pytest.raises(
        GLib.Error,
        match=r"^ncm-mset-error: ncm_mset_load: submodel `NcBBN' appears more than once",
    ):
        _ = Ncm.MSet.load(
            tmp_file.as_posix(), Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
        )


def test_mset_load_submodel_without_main_model(tmp_path):
    """A submodel group without its main model's group is a clean error."""
    tmp_file = _save_lcdm_mset(tmp_path)
    header, groups = _mset_file_groups(tmp_file)
    tmp_file.write_text(header + groups["NcBBN"])

    with pytest.raises(
        GLib.Error,
        match=(
            r"^ncm-mset-error: ncm_mset_load: cannot add submodel `NcBBNParthenope', "
            r"main model `NcHICosmo' not found"
        ),
    ):
        _ = Ncm.MSet.load(
            tmp_file.as_posix(), Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
        )


def test_mset_load_submodel_array_rejected(tmp_path):
    """The generic submodel-array property is forbidden in mset keyfiles."""
    tmp_file = _save_lcdm_mset(tmp_path)
    header, groups = _mset_file_groups(tmp_file)
    groups["NcHICosmo"] = groups["NcHICosmo"].rstrip() + "\nsubmodel-array=@ao []\n"
    tmp_file.write_text(header + "".join(groups.values()))

    with pytest.raises(
        GLib.Error,
        match=(
            r"^ncm-mset-error: ncm_mset_load: serialized version of mset models "
            r"cannot contain the submodel-array property"
        ),
    ):
        _ = Ncm.MSet.load(
            tmp_file.as_posix(), Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
        )


@pytest.mark.parametrize("group", ["NcHICosmo", "NcBBN"])
def test_mset_load_group_missing_type_key(tmp_path, group):
    """A group without its self-referencing type key is a clean error (host and submodel)."""
    tmp_file = _save_lcdm_mset(tmp_path)
    header, groups = _mset_file_groups(tmp_path / "cosmo.mset")
    groups[group] = re.sub(rf"(?m)^{group}=.*\n", "", groups[group])
    tmp_file.write_text(header + "".join(groups.values()))

    with pytest.raises(
        GLib.Error,
        match=r"^ncm-mset-error: ncm_mset_load: Every group must contain a key with same name",
    ):
        _ = Ncm.MSet.load(
            tmp_file.as_posix(), Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
        )


def test_mset_load_stackpos_qualified_submodel_group(tmp_path):
    """A stackpos-qualified submodel group name still resolves to its namespace."""
    tmp_file = _save_lcdm_mset(tmp_path)
    header, groups = _mset_file_groups(tmp_file)
    groups["NcBBN"] = groups["NcBBN"].replace("[NcBBN]", "[NcBBN:1]")
    tmp_file.write_text(header + "".join(groups.values()))

    mset = Ncm.MSet.load(
        tmp_file.as_posix(), Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
    )
    cosmo = mset.peek(Nc.HICosmo.id())
    assert type(cosmo.peek_bbn()).__name__ == "BBNParthenope"


def test_model_qualified_param_unknown_slot():
    """A qualified name naming a slot the model does not have is an error."""
    cosmo = Nc.HICosmoLCDM.new()

    with pytest.raises(
        GLib.Error,
        match=(
            r"^ncm-model-error: ncm_model_param_index_from_name_full: model "
            r"`NcHICosmoLCDM' has no attached submodel slot `nosuchslot'"
        ),
    ):
        _ = cosmo["nosuchslot:x"]


def test_model_qualified_param_slotless_model():
    """A qualified name on a model with no slots at all is an error."""
    prim = Nc.HIPrimPowerLaw.new()

    with pytest.raises(
        GLib.Error,
        match=(
            r"^ncm-model-error: ncm_model_param_index_from_name_full: model "
            r"`NcHIPrimPowerLaw' has no attached submodel slot `a'"
        ),
    ):
        _ = prim["a:b"]


def test_model_qualified_param_unknown_param_in_slot():
    """A qualified name with a valid slot but unknown parameter is an error."""
    cosmo = Nc.HICosmoLCDM(reion=Nc.HIReionCamb.new())

    with pytest.raises(
        GLib.Error,
        match=(
            r"^ncm-model-error: ncm_model_param_index_from_name_full: submodel slot "
            r"`reion' of model `NcHICosmoLCDM' does not have a parameter called `frobnicate'"
        ),
    ):
        _ = cosmo["reion:frobnicate"]


def test_model_qualified_param_empty_slot():
    """A qualified name for a declared but unfilled slot is an error."""
    cosmo = Nc.HICosmoLCDM.new()

    with pytest.raises(
        GLib.Error,
        match=(
            r"^ncm-model-error: ncm_model_param_index_from_name_full: model "
            r"`NcHICosmoLCDM' has no attached submodel slot `reion'"
        ),
    ):
        _ = cosmo["reion:z_re"]

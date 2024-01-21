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

from numcosmo_py import Ncm

Ncm.cfg_init()


def test_mset_new_empty():
    """Test the NcmMSet constructor."""
    mset = Ncm.MSet.empty_new()
    assert mset.total_len() == 0
    assert mset.fparam_len() == 0
    assert mset.nmodels() == 0

    mset.prepare_fparam_map()

    assert mset.fparams_len() == 0


def test_mset_new_array():
    """Test the NcmMSet constructor."""
    mset = Ncm.MSet.new_array(
        [Ncm.ModelFunnel.new(5), Ncm.ModelRosenbrock.new(), Ncm.ModelMVND.new(8)]
    )
    assert mset.total_len() == 16
    assert mset.fparam_len() == 0
    assert mset.nmodels() == 3

    mset.prepare_fparam_map()

    assert mset.fparams_len() == 0


def test_mset_peek_model():
    """Test the NcmMSet peek_model function."""
    mset = Ncm.MSet.new_array(
        [Ncm.ModelFunnel.new(5), Ncm.ModelRosenbrock.new(), Ncm.ModelMVND.new(8)]
    )
    assert isinstance(mset.peek(Ncm.ModelMVND.id()), Ncm.ModelMVND)
    assert isinstance(mset.peek(Ncm.ModelFunnel.id()), Ncm.ModelFunnel)
    assert isinstance(mset.peek(Ncm.ModelRosenbrock.id()), Ncm.ModelRosenbrock)


def test_mset_peek_model_by_name():
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


def test_mset_fparam_map():
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


def test_mset_fparam_map_change():
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


def test_mset_fparam_get_pi_by_name():
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

#!/usr/bin/env python
#
# test_py_data_funnel.py
#
# thu Dec 29 15:44:00 2023
# Copyright  2023  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_py_data_funnel.py
# Copyright (C) 2023 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

"""Tests on NcmDataFunnel class."""

import math
from numcosmo_py import Ncm

Ncm.cfg_init()


def test_constructor():
    """Test constructor."""

    funnel = Ncm.DataFunnel.new()
    assert funnel is not None
    assert isinstance(funnel, Ncm.DataFunnel)

    funnel2 = funnel.ref()
    assert funnel2 == funnel

    assert funnel.get_dof() > 0
    assert funnel.get_length() > 0


def test_eval():
    """Test eval."""

    funnel = Ncm.DataFunnel()
    funnel_model = Ncm.ModelFunnel.new(5)
    mset = Ncm.MSet.new_array([funnel_model])

    assert math.isfinite(funnel.m2lnL_val(mset))


def test_no_mean_vector():
    """Test no mean vector."""
    funnel = Ncm.DataFunnel.new()

    assert not funnel.has_mean_vector()


if __name__ == "__main__":
    test_constructor()
    test_eval()

#!/usr/bin/env python
#
# test_py_data_gaussmix2d.py
#
# thu Dec 30 19:02:32 2023
# Copyright  2023  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_py_data_gaussmix2d.py
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

"""Tests on NcmDataGaussMix2D class."""

import math
from numcosmo_py import Ncm

Ncm.cfg_init()


def test_constructor():
    """Test constructor."""

    gaussmix2d = Ncm.DataGaussMix2D.new()
    assert gaussmix2d is not None
    assert isinstance(gaussmix2d, Ncm.DataGaussMix2D)

    gaussmix2d2 = gaussmix2d.ref()
    assert gaussmix2d2 == gaussmix2d

    assert gaussmix2d.get_dof() > 0
    assert gaussmix2d.get_length() > 0


def test_eval():
    """Test eval."""

    gaussmix2d = Ncm.DataGaussMix2D()
    gaussmix2d_model = Ncm.ModelRosenbrock.new()
    mset = Ncm.MSet.new_array([gaussmix2d_model])

    assert math.isfinite(gaussmix2d.m2lnL_val(mset))


if __name__ == "__main__":
    test_constructor()
    test_eval()

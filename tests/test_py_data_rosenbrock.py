#!/usr/bin/env python
#
# test_py_data_rosenbrock.py
#
# thu Dec 31 11:29:00 2023
# Copyright  2023  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_py_data_rosenbrock.py
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

"""Tests on NcmDataRosenbrock class."""

import math
from numcosmo_py import Ncm

Ncm.cfg_init()


def test_constructor():
    """Test constructor."""

    rosenbrock = Ncm.DataRosenbrock.new()
    assert rosenbrock is not None
    assert isinstance(rosenbrock, Ncm.DataRosenbrock)

    rosenbrock2 = rosenbrock.ref()
    assert rosenbrock2 == rosenbrock

    assert rosenbrock.get_dof() > 0
    assert rosenbrock.get_length() > 0


def test_eval():
    """Test eval."""

    rosenbrock = Ncm.DataRosenbrock()
    rosenbrock_model = Ncm.ModelRosenbrock.new()
    mset = Ncm.MSet.new_array([rosenbrock_model])

    assert math.isfinite(rosenbrock.m2lnL_val(mset))


if __name__ == "__main__":
    test_constructor()
    test_eval()

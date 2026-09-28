#!/usr/bin/env python
#
# test_py_spline_vec.py
#
# Sat Mar 15 19:53:22 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_py_spline_vec.py
# Copyright (C) 2026 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

"""Python bindings of NcmSplineVec; the numerical tests are in test_ncm_spline_vec.c."""

import math

import numpy as np

from numcosmo_py import Ncm

Ncm.cfg_init()


def test_spline_vec_bindings() -> None:
    """The list constructor and the GArray outputs agree with the component splines."""
    x_arr = np.linspace(0.0, 10.0, 50)
    xv = Ncm.Vector.new_array(x_arr.tolist())
    yv_list = [
        Ncm.Vector.new_array((x_arr**2).tolist()),
        Ncm.Vector.new_array(np.cos(x_arr).tolist()),
    ]
    sv = Ncm.SplineVec.new_gpa(Ncm.SplineCubicNotaknot.new(), xv, yv_list, True)

    assert sv.get_len() == 2
    assert sv.get_nknots() == 50

    res = sv.eval_array(5.0)
    dres = sv.deriv_array(5.0)
    ires = sv.integ_array(2.0, 7.0)

    for i in range(2):
        s_i = sv.peek_spline(i)
        assert res[i] == s_i.eval(5.0)
        assert dres[i] == s_i.eval_deriv(5.0)
        assert ires[i] == s_i.eval_integ(2.0, 7.0)

    assert math.isclose(res[0], 25.0, rel_tol=1.0e-14)

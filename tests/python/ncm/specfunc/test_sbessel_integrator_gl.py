#!/usr/bin/env python
#
# test_py_sbessel_integrator_gl.py
#
# Thu Jan 09 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_py_sbessel_integrator_gl.py
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

"""Python bindings of NcmSBesselIntegratorGL.

Accuracy and the NcmSBesselIntegrator contract are tested in C
(tests/c/ncm/specfunc/test_ncm_sbessel_integrator.c); this checks what only Python
reaches: the callback path and the properties.
"""

import numpy as np

from numcosmo_py import Ncm

Ncm.cfg_init()


def test_properties() -> None:
    """Properties set at construction and by the setters read back."""
    integrator = Ncm.SBesselIntegratorGL(
        ell_range=Ncm.DTuple2.new(2, 7), npts=20, margin=3.0, nosc=4.0
    )

    assert integrator.get_ell_range() == (2, 7)
    assert integrator.props.npts == 20
    assert integrator.props.margin == 3.0
    assert integrator.props.nosc == 4.0

    integrator.set_npts(12)
    integrator.set_margin(6.0)
    integrator.set_nosc(8.0)

    assert integrator.get_npts() == 12
    assert integrator.get_margin() == 6.0
    assert integrator.get_nosc() == 8.0


def test_default_range() -> None:
    """Without ell_range the range is [0, 0]."""
    assert Ncm.SBesselIntegratorGL().get_ell_range() == (0, 0)


def test_python_callback_matches_c_shape() -> None:
    """A Python callback gives the numbers of the same shape evaluated in C."""
    center, std, a, b, k = 0.5, 0.05, 0.1, 0.9, 7.0
    integrator = Ncm.SBesselIntegratorGL.new(0, 6)

    def gauss(chi: float, _k: float) -> float:
        return np.exp(-0.5 * ((chi - center) / std) ** 2)

    res_py = Ncm.Vector.new(7)
    res_c = Ncm.Vector.new(7)
    integrator.integrate(gauss, a, b, k, res_py)
    integrator.integrate_gaussian(center, std, a, b, k, res_c)

    for ell in range(7):
        assert res_py.get(ell) == res_c.get(ell)
        assert integrator.integrate_ell(gauss, a, b, k, ell) == res_c.get(ell)

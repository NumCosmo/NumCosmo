#
# py_sline_mfunc.py
#
# Mon May 22 16:00:00 2023
# Copyright  2023  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# py_sline_mfunc.py
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

"""Example of defining a new MSetFunc to be used in a MCMC analysis."""


from scipy import integrate
from py_sline_model import PySLineModel

from numcosmo_py import Ncm, GObject


#
# Creating a new MSetFunc
#
class PyTestFunc(Ncm.MSetFunc1):
    """Example of defining a new MSetFunc to be used in a MCMC analysis."""

    def __init__(self):
        super().__init__(self, dimension=1, nvariables=0)
        self.symbol = r"\int_0^{10} f(x)\mathrm{d}x"
        self.name = r"intf"

    def do_eval1(self, mset, _):  # pylint: disable-msg=arguments-differ
        mid = mset.get_id_by_ns("NcPySLineModel")
        slm = mset.peek(mid)
        assert isinstance(slm, PySLineModel)

        res = integrate.quad(slm.f_x, 0.0, 10.0)

        return [res[0]]


#
# Register our new MSetFunc
#
GObject.type_register(PyTestFunc)

#
# py_sline_data.py
#
# Mon May 22 16:00:00 2023
# Copyright  2023  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# py_sline_data.py
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

"""Example of a Data class for the SLine model."""

import numpy as np

from py_sline_model import PySLineModel

from numcosmo_py import Ncm, GObject


#
# Creating a new class implementing our object Ncm.Data
#
class PySLineData(Ncm.Data):
    """A simple data class for the SLine model."""

    #
    # Creating three properties, len (length), dof (degrees of freedom) and data (NcmMatrix containing
    # the data).
    #
    len = GObject.Property(
        type=GObject.TYPE_UINT, default=0, flags=GObject.PARAM_READWRITE  # type: ignore
    )
    dof = GObject.Property(
        type=GObject.TYPE_UINT, default=0, flags=GObject.PARAM_READWRITE  # type: ignore
    )
    data = GObject.Property(type=Ncm.Matrix, flags=GObject.PARAM_READWRITE)

    #
    # The contructor assigns some default values and calls the father's constructor.
    #
    def __init__(self, length=600):
        Ncm.Data.__init__(self)
        self.len = length
        self.dof = self.len - 2
        self.data = Ncm.Matrix.new(self.len, 3)

    #
    # Implements the virtual method get_length.
    #
    def do_get_length(self) -> int:  # pylint: disable-msg=arguments-differ
        return self.len

    #
    # Implements the virtual method get_dof.
    #
    def do_get_dof(self) -> int:  # pylint: disable-msg=arguments-differ
        return self.dof

    #
    # Implements the virtual method `begin'.
    # This method usually do some groundwork in the data
    # before the actual calculations. For example, if the likelihood
    # involves the decomposition of a constant matrix, it can be done
    # during `begin' once and then used afterwards.
    #
    def do_begin(self):  # pylint: disable-msg=arguments-differ
        return

    #
    # Implements the virtual method `prepare'.
    # This method should do all the necessary calculations using mset
    # to be able to calculate the likelihood afterwards.
    #
    def do_prepare(self, mset):  # pylint: disable-msg=arguments-differ
        self.dof = self.len - mset.fparams_len()

    #
    # Implements the virtual method `resample'.
    # This method creates a new sample from the fiducial model.
    # It is necessary to the MC analysis but can be skipped if
    # doing only MCMC.
    #
    def do_resample(self, mset, rng):  # pylint: disable-msg=arguments-differ
        mid = mset.get_id_by_ns("NcPySLineModel")
        slm = mset.peek(mid)
        assert isinstance(slm, PySLineModel)

        my_data = []

        for _ in range(self.len):
            x = rng.uniform_gen(0.0, 10.0)
            s = x * rng.uniform_gen(0.4, 0.5)
            v = rng.gaussian_gen(slm.f_x(x), s)

            my_data.append([x, v, s])

        my_array = np.array(my_data).reshape(3 * self.len)
        self.data = Ncm.Matrix.new_array(my_array, 3)
        self.set_init(True)

    #
    # Implements the virtual method `m2lnL'.
    # This method should calculate the value of the likelihood for
    # the model set `mset'.
    #
    def do_m2lnL_val(self, mset):  # pylint: disable-msg=arguments-differ
        mid = mset.get_id_by_ns("NcPySLineModel")
        slm = mset.peek(mid)
        assert isinstance(slm, PySLineModel)

        m2lnL = 0.0
        for i in range(self.len):
            x = self.data.get(i, 0)
            v = self.data.get(i, 1)
            s = self.data.get(i, 2)
            m2lnL += ((slm.f_x(x) - v) / s) ** 2
        return m2lnL


#
# Register our new Python class PySLineData
#
GObject.type_register(PySLineData)

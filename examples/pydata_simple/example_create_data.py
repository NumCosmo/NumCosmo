#!/usr/bin/env python
#
# example_create_data.py
#
# Mon May 22 16:00:00 2023
# Copyright  2023  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# example_create_data.py
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

"""Create a new SLine data object and save it to a file."""

import sys
from typing import Union

import numpy as np
import matplotlib.pyplot as plt

from py_sline_model import PySLineModel
from py_sline_data import PySLineData
from py_sline_gauss import PySLineGauss

from numcosmo_py import Ncm


#
#  Initializing the library objects, this must be called before
#  any other library function.
#
Ncm.cfg_init()


def generate_data() -> None:
    """Create a new SLine data object and save it to a file."""

    #
    # Creating a new NcmRNG object with seed 123 and default algorithm
    #
    rng = Ncm.RNG.seeded_new(None, 123)

    #
    # Instantiating a new SLine model object and setting
    # some values for its parameters.
    #
    slm = PySLineModel()
    slm.props.alpha = 1.0
    slm.props.a = 0.5

    #
    # Instantiating a new empty SLine data object.
    #
    sld: Union[PySLineData, PySLineGauss]
    if (len(sys.argv) != 2) or (sys.argv[1] != "--plain" and sys.argv[1] != "--gauss"):
        print("usage: example_create_data.py --plain or --gauss")
        sys.exit(-1)
    elif sys.argv[1] == "--plain":
        sld = PySLineData(length=50)
    else:
        sld = PySLineGauss(length=50)
        sld.xv.set_array(np.linspace(0.0, 10.0, sld.get_size()))
        sld.create_random_cov(slm, rng)

    #
    # New Model set object including slm with parameters
    # set as free.
    #
    mset = Ncm.MSet.empty_new()
    mset.set(slm)
    mset.param_set_all_ftype(Ncm.ParamType.FREE)
    mset.prepare_fparam_map()

    #
    # Creating a new Serialization object, with a data
    # file does not exists, generate a new sample using
    # mset as fiducial model and save it to data_file.
    #
    # If data_file already exists, reload sld from it.
    #
    ser = Ncm.Serialize.new(Ncm.SerializeOpt.NONE)
    data_file = "example_data.obj"
    sld.resample(mset, rng)
    ser.to_binfile(sld, data_file)

    if isinstance(sld, PySLineGauss):
        #
        # Plotting the created data points
        #

        diag = sld.peek_mean().dup()

        sld.peek_cov().get_diag(diag)

        xa = sld.xv.dup_array()
        plt.errorbar(
            sld.xv.dup_array(),
            sld.peek_mean().dup_array(),
            yerr=diag.dup_array(),
            fmt="o",
            label=r"$y\pm\sigma$",
        )
        plt.plot(xa, [slm.f_x(x) for x in xa], lw=1.5, label="y(x)")

        plt.savefig("example_data_f.pdf")

        plt.show()


if __name__ == "__main__":
    generate_data()

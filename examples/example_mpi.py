#!/usr/bin/env python
#
# example_mpi_fit.py
#
# Mon May 22 16:00:00 2023
# Copyright  2023  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# example_mpi_fit.py
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

"""Example testing MPI fit."""

import sys
import timeit

import numpy as np

from numcosmo_py import Ncm

#
#  Initializing the library objects, this must be called before
#  any other library function.
#
sys.argv = Ncm.cfg_init_full(len(sys.argv), sys.argv)


def test_mpi() -> None:
    """Example testing MPI objects."""

    rng = Ncm.RNG.new(None)
    rng.set_random_seed(True)

    mj = Ncm.MPIJobTest.new()
    mj.set_rand_vector(12, rng)

    ser = Ncm.Serialize.new(Ncm.SerializeOpt.NONE)
    mj.init_all_slaves(ser)

    a = []
    b = []
    for t in np.arange(12.0):
        a.append(Ncm.Vector.new_array([t]))
        b.append(Ncm.Vector.new(1))

    # mj.run_array (a, b)
    t1 = timeit.Timer("mj.run_array (a, b)", "from __main__ import mj, a, b")
    print(f"Timming: {t1.timeit(1)}")

    for a_i, b_i in zip(a, b):
        i = int(a_i.get(0))
        print(
            f"{i:.3d} {b_i.get(0): 22.15g} {mj.props.vector.get(i): 22.15g} "
            f"{b_i.get(0) / mj.props.vector.get(i) - 1.0:e}"
        )

    # time.sleep (600)


if __name__ == "__main__":
    test_mpi()

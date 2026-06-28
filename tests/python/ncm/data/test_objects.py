#!/usr/bin/env python
#
# test_py_test_data_objects.py
#
# Fri May 19 10:00:00 2023
# Copyright  2023  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_py_test_data_objects.py
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

"""Test the data objects of the module py_test_data_objects."""

import os
import pytest
import numpy as np
from numcosmo_py import Ncm, Nc, GObject

Ncm.cfg_init()


def list_data_files():
    datadir = Ncm.cfg_get_data_directory()
    assert datadir is not None
    return [
        os.path.join(datadir, file)
        for file in os.listdir(datadir)
        if file.endswith(".obj")
    ]


@pytest.mark.parametrize("datafile", list_data_files())
def test_load_objects(datafile):
    """Test the loading of the data objects."""

    ser = Ncm.Serialize.new(Ncm.SerializeOpt.NONE)
    # Get the data directory

    # Load the data objects
    obj = ser.from_file(datafile)
    assert obj is not None
    assert isinstance(obj, GObject.Object)


@pytest.mark.parametrize("datafile", list_data_files())
def test_load_objects_run_data(datafile):
    """Test the loading of the data objects."""

    ser = Ncm.Serialize.new(Ncm.SerializeOpt.NONE)
    # Get the data directory

    cosmo = Nc.HICosmoDEXcdm()
    mset = Ncm.MSet.new_array([cosmo])
    dist = Nc.Distance.new(1000)

    dist.prepare(cosmo)

    # Load the data objects
    obj = ser.from_file(datafile)
    assert obj is not None
    assert isinstance(obj, GObject.Object)

    try:
        obj.set_dist(dist)
    except AttributeError:
        pass

    if isinstance(obj, Ncm.Data) and not isinstance(obj, Nc.DataClusterPseudoCounts):
        obj.prepare(mset)
        assert np.isfinite(obj.m2lnL_val(mset))

    if isinstance(obj, Ncm.Spline2d):
        xv = obj.peek_xv()
        yv = obj.peek_yv()
        assert xv is not None
        assert yv is not None
        assert np.all(np.isfinite(xv.dup_array()))
        assert np.all(np.isfinite(yv.dup_array()))

        assert np.isfinite(
            obj.eval(xv.get(0) + xv.get(xv.len() - 1), yv.get(0) + yv.get(yv.len() - 1))
        )

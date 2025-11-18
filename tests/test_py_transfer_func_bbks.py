#!/usr/bin/env python
#
# test_py_transfer_func_bbks.py
#
# Mon Nov 17 23:08:13 2025
# Copyright  2025  Mariana Penna-Lima
# <pennalima@unb.br>
#
# test_py_transfer_func_eh_bbks.py
# Copyright (C) 2025 Mariana Penna-Lima <pennalima@unb.br>
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

"""Unit tests for NcTransferFuncBBKS."""

import numpy as np
from numpy.testing import assert_allclose

from numcosmo_py import Ncm, Nc

Ncm.cfg_init()


def test_transfer_func_bbks_new():
    """Test NcTransferFuncBBKS instantiation."""
    tf = Nc.TransferFuncBBKS.new()
    assert isinstance(tf, Nc.TransferFunc)
    assert isinstance(tf, Nc.TransferFuncBBKS)


def test_transfer_func_bbks_eval():
    """Test NcTransferFuncBBKS evaluation."""
    tf = Nc.TransferFuncBBKS.new()
    cosmo = Nc.HICosmoDEXcdm()
    tf.prepare(cosmo)

    kh_array = np.geomspace(1e-3, 1e2, 50)

    for kh in kh_array:
        T = tf.eval(cosmo, kh)
        assert np.isfinite(T)
        assert T > 0.0
        assert T <= 1.0


def test_transfer_func_bbks_prepare():
    """Test NcTransferFuncBBKS prepare method."""
    tf = Nc.TransferFuncBBKS.new()
    cosmo = Nc.HICosmoDEXcdm()

    tf.prepare(cosmo)

    kh = 0.1
    T = tf.eval(cosmo, kh)
    assert np.isfinite(T)
    assert T > 0.0

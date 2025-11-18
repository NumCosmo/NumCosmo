#!/usr/bin/env python
#
# test_py_transfer_func_eh.py
#
# Mon Nov 17 22:57:27 2025
# Copyright  2025  Mariana Penna-Lima
# <pennalima@unb.br>
#
# test_py_transfer_func_eh.py
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

"""Unit tests for NcTransferFuncEH."""

import numpy as np
from numpy.testing import assert_allclose

from numcosmo_py import Ncm, Nc

Ncm.cfg_init()


def test_transfer_func_eh_new():
    """Test NcTransferFuncEH instantiation."""
    tf = Nc.TransferFuncEH.new()
    assert isinstance(tf, Nc.TransferFunc)
    assert isinstance(tf, Nc.TransferFuncEH)


def test_transfer_func_eh_eval():
    """Test NcTransferFuncEH evaluation."""
    tf = Nc.TransferFuncEH.new()
    cosmo = Nc.HICosmoDEXcdm()
    tf.prepare(cosmo)

    kh_array = np.geomspace(1e-3, 1e2, 50)

    for kh in kh_array:
        T = tf.eval(cosmo, kh)
        assert np.isfinite(T)
        assert T > 0.0
        assert T <= 1.0


def test_transfer_func_eh_prepare():
    """Test NcTransferFuncEH prepare method."""
    tf = Nc.TransferFuncEH.new()
    cosmo = Nc.HICosmoDEXcdm()

    tf.prepare(cosmo)

    kh = 0.1
    T = tf.eval(cosmo, kh)
    assert np.isfinite(T)
    assert T > 0.0


def test_transfer_func_eh_compare_colossus():
    """Test NcTransferFuncEH against Colossus."""

    T_colossus = np.array(
        [
            9.92285969e-01,
            9.88298611e-01,
            9.82470032e-01,
            9.74117687e-01,
            9.62416281e-01,
            9.46428025e-01,
            9.25161657e-01,
            8.97653184e-01,
            8.63049136e-01,
            8.20662543e-01,
            7.69966990e-01,
            7.10508039e-01,
            6.41805108e-01,
            5.63652847e-01,
            4.77770358e-01,
            3.90714369e-01,
            3.13267738e-01,
            2.53158863e-01,
            2.07737002e-01,
            1.61009854e-01,
            1.16475200e-01,
            9.09432713e-02,
            6.43466661e-02,
            4.68390224e-02,
            3.36138847e-02,
            2.37743706e-02,
            1.67008079e-02,
            1.16060302e-02,
            8.01101752e-03,
            5.49512475e-03,
            3.74748065e-03,
            2.54204928e-03,
            1.71592144e-03,
            1.15305372e-03,
            7.71612645e-04,
            5.14391782e-04,
            3.41720222e-04,
            2.26285093e-04,
            1.49405339e-04,
            9.83801541e-05,
            6.46214444e-05,
            4.23507494e-05,
            2.76974529e-05,
            1.80794045e-05,
            1.17803492e-05,
            7.66337017e-06,
            4.97761357e-06,
            3.22856342e-06,
            2.09135006e-06,
            1.35304381e-06,
        ]
    )

    tf = Nc.TransferFuncEH.new()
    cosmo = Nc.HICosmoDEXcdm()

    cosmo.omega_x2omega_k()
    cosmo["H0"] = 67.66
    cosmo["Omegab"] = 0.049
    cosmo["Omegac"] = 0.262
    cosmo["Omegak"] = 0.0
    cosmo["Tgamma0"] = 2.7255
    cosmo["w"] = -1.0

    tf.prepare(cosmo)

    kh_array = np.geomspace(1e-3, 1e2, 50)  # in h/Mpc

    T = [tf.eval(cosmo, kh) for kh in kh_array]
    assert_allclose(T, T_colossus, rtol=1e-8, atol=0.0)

#!/usr/bin/env python
#
# test_py_transfer_func_eh_no_baryon.py
#
# Mon Nov 17 16:19:27 2025
# Copyright  2025  Mariana Penna-Lima
# <pennalima@unb.br>
#
# test_py_transfer_func_eh_no_baryon.py
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

"""Unit tests for NcTransferFuncEHNoBaryon."""

import numpy as np
from numpy.testing import assert_allclose

from numcosmo_py import Ncm, Nc

Ncm.cfg_init()


def test_transfer_func_eh_no_baryon_new():
    """Test NcTransferFuncEHNoBaryon instantiation."""
    tf = Nc.TransferFuncEHNoBaryon.new()
    assert isinstance(tf, Nc.TransferFunc)
    assert isinstance(tf, Nc.TransferFuncEHNoBaryon)


def test_transfer_func_eh_no_baryon_eval():
    """Test NcTransferFuncEHNoBaryon evaluation."""
    tf = Nc.TransferFuncEHNoBaryon.new()
    cosmo = Nc.HICosmoDEXcdm()
    tf.prepare(cosmo)

    kh_array = np.geomspace(1e-3, 1e2, 50)

    for kh in kh_array:
        T = tf.eval(cosmo, kh)
        assert np.isfinite(T)
        assert T > 0.0
        assert T <= 1.0


def test_transfer_func_eh_no_baryon_prepare():
    """Test NcTransferFuncEHNoBaryon prepare method."""
    tf = Nc.TransferFuncEHNoBaryon.new()
    cosmo = Nc.HICosmoDEXcdm()

    tf.prepare(cosmo)

    kh = 0.1
    T = tf.eval(cosmo, kh)
    assert np.isfinite(T)
    assert T > 0.0


def test_transfer_func_eh_no_baryon_compare_colossus():
    """Test NcTransferFuncEHNoBaryon against Colossus."""

    T_colossus = np.array(
        [
            9.92101656e-01,
            9.88128882e-01,
            9.82379996e-01,
            9.74215842e-01,
            9.62860830e-01,
            9.47421895e-01,
            9.26931647e-01,
            9.00414354e-01,
            8.66960811e-01,
            8.25772811e-01,
            7.76093988e-01,
            7.16931088e-01,
            6.46845999e-01,
            5.65560990e-01,
            4.78422909e-01,
            3.95301041e-01,
            3.21777217e-01,
            2.58151363e-01,
            2.03604966e-01,
            1.57619869e-01,
            1.19749466e-01,
            8.93585687e-02,
            6.55862767e-02,
            4.74299697e-02,
            3.38573801e-02,
            2.38997741e-02,
            1.67103324e-02,
            1.15887822e-02,
            7.98085522e-03,
            5.46271571e-03,
            3.71887035e-03,
            2.51930755e-03,
            1.69900174e-03,
            1.14101701e-03,
            7.63309546e-04,
            5.08784368e-04,
            3.37986547e-04,
            2.23821256e-04,
            1.47787924e-04,
            9.73211063e-05,
            6.39285780e-05,
            4.18973531e-05,
            2.74005415e-05,
            1.78847959e-05,
            1.16526924e-05,
            7.57958369e-06,
            4.92260573e-06,
            3.19245045e-06,
            2.06764912e-06,
            1.33749766e-06,
        ]
    )

    tf = Nc.TransferFuncEHNoBaryon.new()
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

#!/usr/bin/env python
#
# test_py_hicosmo.py
#
# Sun Sep 08 14:33:22 2024
# Copyright  2024  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_py_hicosmo.py
# Copyright (C) 2024 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

"""Tests on NcHICosmoidem2 cosmological model."""

from numpy.testing import assert_allclose

from numcosmo_py import Ncm, Nc

Ncm.cfg_init()


def test_nc_hicosmo_idem2():
    """Test NcHICosmoIdem2 initialization."""
    idem2 = Nc.HICosmoIDEM2()
    assert idem2 is not None
    assert isinstance(idem2, Nc.HICosmoIDEM2)
    assert isinstance(idem2, Nc.HICosmo)
    assert isinstance(idem2, Ncm.Model)


def test_nc_hicosmo_idem2_reparam_cmb():
    """Test NcHICosmoIdem2 reparam_cmb."""
    idem2 = Nc.HICosmoIDEM2()
    assert idem2 is not None
    assert isinstance(idem2, Nc.HICosmoIDEM2)
    assert isinstance(idem2, Nc.HICosmo)
    assert isinstance(idem2, Ncm.Model)

    idem2.cmb_params()

    idem2["omegab"] = 0.022
    idem2["omegac"] = 0.12

    assert_allclose(idem2["omegab"], 0.022)
    assert_allclose(idem2["omegac"], 0.12)

    assert_allclose(
        idem2.orig_param_get_by_name("Omegab"), 0.022 / (idem2["H0"] / 100.0) ** 2
    )
    assert_allclose(
        idem2.orig_param_get_by_name("Omegac"), 0.12 / (idem2["H0"] / 100.0) ** 2
    )


def test_nc_hicosmo_idem2_reparam_omega_x2omega_k():
    """Test NcHICosmoIdem2 reparam_omega_x2omega_k."""
    idem2 = Nc.HICosmoIDEM2()
    assert idem2 is not None
    assert isinstance(idem2, Nc.HICosmoIDEM2)
    assert isinstance(idem2, Nc.HICosmo)
    assert isinstance(idem2, Ncm.Model)

    idem2.omega_x2omega_k()

    idem2["Omegak"] = 0.01

    assert_allclose(idem2["Omegak"], 0.01)


def test_nc_hicosmo_gcg():
    """Test NcHICosmoGCG initialization."""
    gcg = Nc.HICosmoGCG()
    assert gcg is not None
    assert isinstance(gcg, Nc.HICosmoGCG)
    assert isinstance(gcg, Nc.HICosmo)
    assert isinstance(gcg, Ncm.Model)


def test_nc_hicosmo_gcg_reparam_cmb():
    """Test NcHICosmoGCG reparam_cmb."""
    gcg = Nc.HICosmoGCG()
    assert gcg is not None
    assert isinstance(gcg, Nc.HICosmoGCG)
    assert isinstance(gcg, Nc.HICosmo)
    assert isinstance(gcg, Ncm.Model)

    gcg.cmb_params()

    gcg["omegab"] = 0.022
    gcg["omegac"] = 0.12

    assert_allclose(gcg["omegab"], 0.022)
    assert_allclose(gcg["omegac"], 0.12)

    assert_allclose(
        gcg.orig_param_get_by_name("Omegab"), 0.022 / (gcg["H0"] / 100.0) ** 2
    )
    assert_allclose(
        gcg.orig_param_get_by_name("Omegac"), 0.12 / (gcg["H0"] / 100.0) ** 2
    )


def test_nc_hicosmo_gcg_reparam_omega_x2omega_k():
    """Test NcHICosmoGCG reparam_omega_x2omega_k."""
    gcg = Nc.HICosmoGCG()
    assert gcg is not None
    assert isinstance(gcg, Nc.HICosmoGCG)
    assert isinstance(gcg, Nc.HICosmo)
    assert isinstance(gcg, Ncm.Model)

    gcg.omega_x2omega_k()

    gcg["Omegak"] = 0.01

    assert_allclose(gcg["Omegak"], 0.01)

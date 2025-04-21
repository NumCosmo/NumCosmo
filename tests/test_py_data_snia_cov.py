#!/usr/bin/env python
#
# test_py_data_snia_cov.py
#
# Sun Apr 21 11:15:07 2024
# Copyright  2024  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_py_data_snia_cov.py
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

"""Tests on NcmDataSNIACov class."""

import pytest

from numcosmo_py import Ncm, Nc

Ncm.cfg_init()


def test_constructor():
    """Test constructor."""
    snia_cov = Nc.DataSNIACov.new(False, 0)
    assert snia_cov is not None
    assert isinstance(snia_cov, Nc.DataSNIACov)

    snia_cov2 = snia_cov.ref()
    assert snia_cov2 == snia_cov


@pytest.mark.parametrize(
    "cov_id",
    [
        Nc.DataSNIAId.COV_DES_Y5_STATONLY,
        Nc.DataSNIAId.COV_DES_Y5_STAT_SYS,
        Nc.DataSNIAId.COV_JLA_SNLS3_SDSS_SYS_STAT,
        Nc.DataSNIAId.COV_JLA_SNLS3_SDSS_SYS_STAT_CMPL,
        Nc.DataSNIAId.COV_PANTHEON,
        Nc.DataSNIAId.COV_PANTHEON_PLUS_SH0ES_STAT,
        Nc.DataSNIAId.COV_PANTHEON_PLUS_SH0ES_SYS_STAT,
        Nc.DataSNIAId.COV_SNLS3_STAT_ONLY,
        Nc.DataSNIAId.COV_SNLS3_SYS_STAT,
    ],
)
def test_constructor_catalog_id(cov_id):
    """Test constructor with catalog id."""
    snia_cov = Nc.DataSNIACov.new_from_cat_id(cov_id, True)
    assert snia_cov is not None
    assert isinstance(snia_cov, Nc.DataSNIACov)

    assert snia_cov.get_dof() > 0
    assert snia_cov.get_length() > 0

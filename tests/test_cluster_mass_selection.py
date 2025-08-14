#!/usr/bin/env python
#
# test_cluster_mass_selection.py
#
# Mon Jan 15 10:19:40 2024
# Copyright  2023  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_py_serialize.py
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
#

"""Tests for the cluster mass with completeness and purity module."""

import pytest
import numpy as np
from numcosmo_py import Nc, Ncm

Ncm.cfg_init()

@pytest.fixture(name="cluster_m_selection")
def fixture_cluster_mass_selection() -> Nc.ClusterMassSelection:
    """Fixture for the NcClusterMassSelection."""
    cluster_m = Nc.ClusterMassSelection(lnR_min = np.log(5), lnR_max = np.log(200))
    return cluster_m

@pytest.fixture(name="cluster_m_selection", params=[1, 2, 3, 4, 5, 6])
def fixture_cluster_mass_selection(request):
    """Fixture for NcClusterMassSelection."""
    return request.param

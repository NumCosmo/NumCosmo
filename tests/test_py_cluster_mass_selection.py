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

lnM_bins_knots = np.linspace(13.0 * np.log(10), 16 * np.log(10) ,10)
z_bins_knots = np.array([0.1 , 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0, 1.1])

lnM_centers = 0.5 * (lnM_bins_knots[:-1] + lnM_bins_knots[1:])
z_centers = 0.5 * (z_bins_knots[:-1] + z_bins_knots[1:])

completeness_knots = Ncm.Matrix.new(len(z_centers),len(lnM_centers))
completeness_smooth_clipped = np.array([[0.00775574, 0.0257542, 0.03300902, 0.03199418, 0.02518367, 0.01505147, 0.00407155, 0.0, 0.0],
    [0.25260502, 0.27061681, 0.27605795, 0.2718841, 0.26105089, 0.24651397, 0.23122899, 0.21336439, 0.21044216],
    [0.40262404, 0.43122535, 0.4456371, 0.44914992, 0.44505445, 0.43664132, 0.42720117, 0.41831353, 0.42562491],
    [0.50266196, 0.54495916, 0.57169553, 0.58647906, 0.59291775, 0.59461957, 0.59519252, 0.60182781, 0.62621808],
    [0.59756795, 0.64919758, 0.68418232, 0.70655894, 0.72036417, 0.72963479, 0.73840753, 0.75946354, 0.80210609],
    [0.73219118, 0.78131995, 0.81304655, 0.83207695, 0.84311712, 0.85087303, 0.86005064, 0.88677706, 0.94317338],
    [0.95138081, 0.9787056, 0.98823728, 0.98572051, 0.97689998, 0.96752036, 0.96332631, 0.97932468, 1.0],
    [0.95138081, 0.9787056, 0.98823728, 0.98572051, 0.97689998, 0.96752036, 0.96332631, 0.97932468, 1.0],
    [0.95138081, 0.9787056, 0.98823728, 0.98572051, 0.97689998, 0.96752036, 0.96332631, 0.97932468, 1.0]])
 

for i in range(len(z_centers)):
    for j in range(len(lnM_centers)):
        completeness_knots.set(i,j, completeness_smooth_clipped.T[i][j])

completeness = Ncm.Spline2dBicubic(
spline=Ncm.SplineCubicNotaknot.new(),
x_vector=Ncm.Vector.new_array(npa_to_seq(lnM_centers)),
y_vector=Ncm.Vector.new_array(npa_to_seq(z_centers)),
z_matrix=completeness_knots
)

lnR_bins_knots = np.log(np.array([5, 10, 15, 20, 35, 70, 100, 200]))
z_bins_knots = np.array([0.1, 0.2 ,0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0, 1.1])

lnR_centers = 0.5 * (lnR_bins_knots[:-1] + lnR_bins_knots[1:])
z_centers = 0.5 * (z_bins_knots[:-1] + z_bins_knots[1:])

purity_knots = Ncm.Matrix.new(len(z_centers),len(lnR_centers))

purity_smooth_clipped = np.array([[1/0.56847321, 1/0.56553304, 1/0.56643402, 1/0.57076391, 1/0.57811045, 1/0.58806141,
  1/0.60020452, 1/0.6216277,  1/0.64566435],
 [1/0.67394648, 1/0.63268226, 1/0.61286432, 1/0.61046677, 1/0.62146372, 1/0.64182928,
  1/0.66753754, 1/0.70731085, 1/0.7364596 ],
 [1/0.70418518, 1/0.67834512, 1/0.66613058, 1/0.66496127, 1/0.67225691, 1/0.6854372,
  1/0.70192186, 1/0.72720015, 1/0.74539913],
 [1/0.71291281, 1/0.72903127, 1/0.73877806, 1/0.74339086, 1/0.7441073,  1/0.74216506,
  1/0.73880177, 1/0.73379976, 1/0.73256222],
 [1/0.69006831, 1/0.77528975, 1/0.82370733, 1/0.84235979, 1/0.83828587, 1/0.81852433,
  1/0.79011389, 1/0.74667884, 1/0.72337672],
 [1/0.66000857, 1/0.76869821, 1/0.83518768, 1/0.86739902, 1/0.87325424, 1/0.86067539,
  1/0.8375845,  1/0.80056751, 1/0.78445986],
 [1/0.64127851, 1/0.70230658, 1/0.75435678, 1/0.7986071,  1/0.83623552, 1/0.86842004,
  1/0.89633865, 1/0.93279484, 1/0.96627892]])
for i in range(len(z_centers)):
    for j in range(len(lnR_centers)):
        purity_knots.set(i,j, purity_smooth_clipped.T[i][j])

purity = Ncm.Spline2dBicubic(
spline=Ncm.SplineCubicNotaknot.new(),
x_vector=Ncm.Vector.new_array(npa_to_seq(lnR_centers)),
y_vector=Ncm.Vector.new_array(npa_to_seq(z_centers)),
z_matrix=purity_knots
)


@pytest.fixture(name="cluster_m_selection")
def fixture_cluster_mass_selection() -> Nc.ClusterMassSelection:
    """Fixture for the NcClusterMassSelection."""
    
    cluster_m = Nc.ClusterMassSelection(lnR_min = np.log(5), lnR_max = np.log(200))
    return cluster_m

@pytest.fixture(name="cluster_m_selection", params=[1, 2, 3, 4, 5, 6])
def fixture_cluster_mass_selection(request):
    """Fixture for NcClusterMassSelection."""
    
    cluster_m = Nc.ClusterMassSelection(lnR_min = np.log(5), lnR_max = np.log(200))
    return request.param

@pytest.fixture(name="cluster_m_completeness")
def test_completeness(cluster_m):
    """Fixture for NcClusterMassSelection."""
    cluster_m = Nc.ClusterMassSelection(lnR_min = np.log(5), lnR_max = np.log(200))
    
    assert(cluster_m.get_completeness() is None)
    assert(cluster_m.completeness(13.0 * np.log(10) , 0.1) == 1)

    
    cluster_m.set_completeness(completeness)
    assert(cluster_m.get_completeness() is not None)
    
    return 


@pytest.fixture(name="cluster_m_purity")
def test_purity(cluster_m):
    """Fixture for NcClusterMassSelection."""
    cluster_m = Nc.ClusterMassSelection(lnR_min = np.log(5), lnR_max = np.log(200))
    
    assert(cluster_m.get_purity() is None)
    assert(cluster_m.purity(np.log(5) , 0.1) == 1)

   
    cluster_m.set_purity(purity)
    assert(cluster_m.get_purity() is not None)
    
    return 


def test_moments(cluster_m):
     """Fixture for NcClusterMassSelection."""
    cluster_m = Nc.ClusterMassSelection(lnR_min = np.log(5), lnR_max = np.log(200))
    cluster_m.param_set_by_name("mup0", 4.12769558168741)
    cluster_m.param_set_by_name("mup1", 1.17476066603899)
    cluster_m.param_set_by_name("mup2", 0.393577193825473 )
    cluster_m.param_set_by_name("sigmap0", 0.408750324989284)
    cluster_m.param_set_by_name("sigmap1", -0.123232985316648)
    cluster_m.param_set_by_name("sigmap2", -0.0644996574273048 )
    cluster_m.set_completeness(completeness)
    cluster_m.set_purity(purity)

    assert(cluster_m.get_mean(13.5,0.1)  is not None)
    assert(cluster_m.get_std(13.5,0.1)  is not None)


def test_distribution(cluster_m):
     """Fixture for NcClusterMassSelection."""
    lnR_bins_test = np.linspace(np.log(5),np.log(200),100)
    lnR_centers_test = 0.5 * (lnR_bins_test[:-1] + lnR_bins_test[1:])
    cluster_m = Nc.ClusterMassSelection(lnR_min = np.log(5), lnR_max = np.log(200))
    cluster_m.param_set_by_name("mup0", 4.12769558168741)
    cluster_m.param_set_by_name("mup1", 1.17476066603899)
    cluster_m.param_set_by_name("mup2", 0.393577193825473 )
    cluster_m.param_set_by_name("sigmap0", 0.408750324989284)
    cluster_m.param_set_by_name("sigmap1", -0.123232985316648)
    cluster_m.param_set_by_name("sigmap2", -0.0644996574273048 )
    cluster_m.set_completeness(completeness)
    cluster_m.set_purity(purity)

    cluster_p = []
    time_p = []

    cluster_intp_bin = []
    time_intp_bin = []
    
    cluster_intp_init = time.time()
        cluster_intp_bin.append(cosmo, 13.5,0.1, None)
    cluster_intp_end = time.time()

    for i in range(len(lnR_bins_test)-1):
        cluster_p_init = time.time()
        cluster_p.append(cosmo, 13.5,0.1,lnR_centers_test[i], None)
        cluster_p_end = time.time()
        time_p.append(cluster_p_end - cluster_p_init)

        cluster_intp_bin_init = time.time()
        cluster_intp_bin.append(cosmo, 13.5,0.1,lnR_bins_test[i], lnR_bins_test[i+1], None)
        cluster_intp_bin_end = time.time()
        time_intp_bin.append(cluster_p_end - cluster_p_init)
        
        
        
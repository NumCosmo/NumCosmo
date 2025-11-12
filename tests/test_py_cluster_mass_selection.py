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
import timeit
import math
import numpy as np
from numcosmo_py import Nc, Ncm
from numcosmo_py.helper import npa_to_seq

Ncm.cfg_init()


@pytest.fixture(name="cluster_mass_selection")
def fixture_cluster_mass_selection():
    """"Fixture for the NcClusterMassSelection.""" ""

    # cosmological model
    cosmo = Nc.HICosmoDEXcdm()
    # cosmo.omega_x2omega_k()
    # cosmo.param_set_by_name("Omegak", 0.00)
    cosmo.param_set_by_name("Omegax", 1 - 0.2603)
    cosmo.param_set_by_name("H0", 71)
    cosmo.param_set_by_name("Omegab", 0.0406)
    cosmo.param_set_by_name("Omegac", 0.22)  # 0.2603
    cosmo.param_set_by_name("w", -1.0)  # -1.0

    prim = Nc.HIPrimPowerLaw.new()
    prim.props.n_SA = 0.967

    reion = Nc.HIReionCamb.new()

    cosmo.add_submodel(prim)
    cosmo.add_submodel(reion)

    tf = Nc.TransferFuncEH()

    psml = Nc.PowspecMLTransfer.new(tf)
    psml.require_kmin(1.0e-6)
    psml.require_kmax(1.0e3)

    psf = Ncm.PowspecFilter.new(psml, Ncm.PowspecFilterType.TOPHAT)
    psf.set_best_lnr0()
    psf.prepare(cosmo)

    old_amplitude = np.exp(prim.props.ln10e10ASA)
    prim.props.ln10e10ASA = np.log((0.8 / cosmo.sigma8(psf)) ** 2 * old_amplitude)

    cut = np.log(5)
    cluster_m = Nc.ClusterMassSelection(lnRichness_min=cut, lnRichness_max=np.log(200))
    cluster_m.param_set_by_name("mup0", 4.12769558168741)
    cluster_m.param_set_by_name("mup1", 1.17476066603899)
    cluster_m.param_set_by_name("mup2", 0.393577193825473)
    cluster_m.param_set_by_name("sigmap0", 0.408750324989284)
    cluster_m.param_set_by_name("sigmap1", -0.123232985316648)
    cluster_m.param_set_by_name("sigmap2", -0.0644996574273048)
    cluster_m.param_set_by_name("cut", cut)

    assert cluster_m.get_cut() == cut

    lnM_bins_knots = np.linspace(13.0 * np.log(10), 16 * np.log(10), 10)
    z_bins_knots = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0, 1.1])
    lnM_centers = 0.5 * (lnM_bins_knots[:-1] + lnM_bins_knots[1:])
    z_centers = 0.5 * (z_bins_knots[:-1] + z_bins_knots[1:])

    completeness_knots = Ncm.Matrix.new(len(z_centers), len(lnM_centers))

    completeness_smooth_clipped = np.array(
        [
            [
                0.00775574,
                0.0257542,
                0.03300902,
                0.03199418,
                0.02518367,
                0.01505147,
                0.00407155,
                0.0,
                0.0,
            ],
            [
                0.25260502,
                0.27061681,
                0.27605795,
                0.2718841,
                0.26105089,
                0.24651397,
                0.23122899,
                0.21336439,
                0.21044216,
            ],
            [
                0.40262404,
                0.43122535,
                0.4456371,
                0.44914992,
                0.44505445,
                0.43664132,
                0.42720117,
                0.41831353,
                0.42562491,
            ],
            [
                0.50266196,
                0.54495916,
                0.57169553,
                0.58647906,
                0.59291775,
                0.59461957,
                0.59519252,
                0.60182781,
                0.62621808,
            ],
            [
                0.59756795,
                0.64919758,
                0.68418232,
                0.70655894,
                0.72036417,
                0.72963479,
                0.73840753,
                0.75946354,
                0.80210609,
            ],
            [
                0.73219118,
                0.78131995,
                0.81304655,
                0.83207695,
                0.84311712,
                0.85087303,
                0.86005064,
                0.88677706,
                0.94317338,
            ],
            [
                0.95138081,
                0.9787056,
                0.98823728,
                0.98572051,
                0.97689998,
                0.96752036,
                0.96332631,
                0.97932468,
                1.0,
            ],
            [
                0.95138081,
                0.9787056,
                0.98823728,
                0.98572051,
                0.97689998,
                0.96752036,
                0.96332631,
                0.97932468,
                1.0,
            ],
            [
                0.95138081,
                0.9787056,
                0.98823728,
                0.98572051,
                0.97689998,
                0.96752036,
                0.96332631,
                0.97932468,
                1.0,
            ],
        ]
    )
    for i in range(len(z_centers)):
        for j in range(len(lnM_centers)):
            completeness_knots.set(i, j, completeness_smooth_clipped.T[i][j])

    completeness = Ncm.Spline2dBicubic(
        spline=Ncm.SplineCubicNotaknot.new(),
        x_vector=Ncm.Vector.new_array(npa_to_seq(lnM_centers)),
        y_vector=Ncm.Vector.new_array(npa_to_seq(z_centers)),
        z_matrix=completeness_knots,
    )

    lnR_bins_knots = np.log(np.array([5, 10, 15, 20, 35, 70, 100, 200]))
    z_bins_knots = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0, 1.1])

    lnR_centers = 0.5 * (lnR_bins_knots[:-1] + lnR_bins_knots[1:])
    z_centers = 0.5 * (z_bins_knots[:-1] + z_bins_knots[1:])

    ipurity_knots = Ncm.Matrix.new(len(z_centers), len(lnR_centers))

    ipurity_smooth_clipped = np.array(
        [
            [
                1 / 0.56847321,
                1 / 0.56553304,
                1 / 0.56643402,
                1 / 0.57076391,
                1 / 0.57811045,
                1 / 0.58806141,
                1 / 0.60020452,
                1 / 0.6216277,
                1 / 0.64566435,
            ],
            [
                1 / 0.67394648,
                1 / 0.63268226,
                1 / 0.61286432,
                1 / 0.61046677,
                1 / 0.62146372,
                1 / 0.64182928,
                1 / 0.66753754,
                1 / 0.70731085,
                1 / 0.7364596,
            ],
            [
                1 / 0.70418518,
                1 / 0.67834512,
                1 / 0.66613058,
                1 / 0.66496127,
                1 / 0.67225691,
                1 / 0.6854372,
                1 / 0.70192186,
                1 / 0.72720015,
                1 / 0.74539913,
            ],
            [
                1 / 0.71291281,
                1 / 0.72903127,
                1 / 0.73877806,
                1 / 0.74339086,
                1 / 0.7441073,
                1 / 0.74216506,
                1 / 0.73880177,
                1 / 0.73379976,
                1 / 0.73256222,
            ],
            [
                1 / 0.69006831,
                1 / 0.77528975,
                1 / 0.82370733,
                1 / 0.84235979,
                1 / 0.83828587,
                1 / 0.81852433,
                1 / 0.79011389,
                1 / 0.74667884,
                1 / 0.72337672,
            ],
            [
                1 / 0.66000857,
                1 / 0.76869821,
                1 / 0.83518768,
                1 / 0.86739902,
                1 / 0.87325424,
                1 / 0.86067539,
                1 / 0.8375845,
                1 / 0.80056751,
                1 / 0.78445986,
            ],
            [
                1 / 0.64127851,
                1 / 0.70230658,
                1 / 0.75435678,
                1 / 0.7986071,
                1 / 0.83623552,
                1 / 0.86842004,
                1 / 0.89633865,
                1 / 0.93279484,
                1 / 0.96627892,
            ],
        ]
    )
    for i in range(len(z_centers)):
        for j in range(len(lnR_centers)):
            ipurity_knots.set(i, j, ipurity_smooth_clipped.T[i][j])

    ipurity = Ncm.Spline2dBicubic(
        spline=Ncm.SplineCubicNotaknot.new(),
        x_vector=Ncm.Vector.new_array(npa_to_seq(lnR_centers)),
        y_vector=Ncm.Vector.new_array(npa_to_seq(z_centers)),
        z_matrix=ipurity_knots,
    )

    return cluster_m, completeness, ipurity, cosmo, psf


def test_cluster_mass_selection_completeness_and_purity(cluster_mass_selection):

    nsize = 500
    lnM = np.linspace(np.log(1e12), np.log(1e16), nsize)
    z = np.linspace(0, 2, nsize)
    lnR = np.linspace(np.log(1), np.log(1000), nsize)
    cluster_m, completeness, ipurity, cosmo, psf = cluster_mass_selection

    assert cluster_m.peek_completeness() is None
    assert cluster_m.peek_ipurity() is None

    for i in range(nsize):
        assert cluster_m.completeness(lnM[i], 0.5) == 1.0
        assert cluster_m.completeness(np.log(1e14), z[i]) == 1.0
        assert cluster_m.ipurity(np.log(5), z[i]) == 1.0
        assert cluster_m.ipurity(lnR[i], 0.5) == 1.0

    cluster_m.set_ipurity(ipurity)
    cluster_m.set_completeness(completeness)

    assert cluster_m.peek_completeness() is completeness
    assert cluster_m.peek_ipurity() is ipurity


def test_cluster_mass_selection_mean_std(cluster_mass_selection):

    nsize = 100
    cluster_m, completeness, ipurity, cosmo, psf = cluster_mass_selection
    cluster_m.set_ipurity(ipurity)
    cluster_m.set_completeness(completeness)

    lnM = np.linspace(np.log(1e13), np.log(1e16), nsize)
    z = np.linspace(0, 1.1, nsize)

    for i in range(nsize):
        assert cluster_m.get_mean(lnM[i], 0.5) >= cluster_m.get_mean_richness(
            lnM[i], 0.5
        )
        assert cluster_m.get_mean(np.log(1e14), z[i]) >= cluster_m.get_mean_richness(
            np.log(1e14), z[i]
        )

        assert cluster_m.get_std(lnM[i], 0.5) <= cluster_m.get_std_richness(lnM[i], 0.5)
        assert cluster_m.get_std(np.log(1e14), z[i]) <= cluster_m.get_std_richness(
            np.log(1e14), z[i]
        )

    cluster_m.param_set_by_name("cut", -1e20)
    for i in range(nsize):
        assert cluster_m.get_mean(lnM[i], 0.5) == cluster_m.get_mean_richness(
            lnM[i], 0.5
        )
        assert cluster_m.get_mean(np.log(1e14), z[i]) == cluster_m.get_mean_richness(
            np.log(1e14), z[i]
        )

        assert cluster_m.get_std(lnM[i], 0.5) == cluster_m.get_std_richness(lnM[i], 0.5)
        assert cluster_m.get_std(np.log(1e14), z[i]) == cluster_m.get_std_richness(
            np.log(1e14), z[i]
        )


def test_cluster_mass_selection_distribution(cluster_mass_selection):

    nsize = 100
    cluster_m, completeness, ipurity, cosmo, psf = cluster_mass_selection
    cluster_m.set_ipurity(ipurity)
    cluster_m.set_completeness(completeness)

    lnM = np.linspace(np.log(1e13), np.log(1e16), nsize)
    z = np.linspace(0, 1.1, nsize)
    lnR = np.linspace(np.log(5), np.log(200), nsize)

    my_globals = globals().copy()
    my_globals.update(
        {
            "nsize": nsize,
            "lnM": lnM,
            "z": z,
            "lnR": lnR,
            "np": np,
            "cluster_m": cluster_m,
            "cosmo": cosmo,
        }
    )

    cluster_p_mass_test = (
        """[cluster_m.p(cosmo, lnM[i], 0.5, [np.log(20)],None) for i in range(nsize)]"""
    )
    execution_time = timeit.timeit(cluster_p_mass_test, globals=my_globals, number=100)
    print(
        f"Average time per execution cluster_m.p mass: {execution_time / 100:.6f} seconds"
    )

    cluster_p_z_test = """[cluster_m.p(cosmo, np.log(1e14), z[i], [np.log(20)],None) for i in range(nsize)]"""
    execution_time = timeit.timeit(cluster_p_z_test, globals=my_globals, number=100)
    print(
        f"Average time per execution cluster_m.p redshift: {execution_time / 100:.6f} seconds"
    )

    cluster_p_richness_test = """[cluster_m.p(cosmo, np.log(1e14), 0.5, [lnR[i]],None) for i in range(nsize)]"""
    execution_time = timeit.timeit(
        cluster_p_richness_test, globals=my_globals, number=100
    )
    print(
        f"Average time per execution cluster_m.p richness: {execution_time / 100:.6f} seconds"
    )


def test_cluster_mass_selection_cumulative(cluster_mass_selection):

    nsize = 100
    cluster_m, completeness, ipurity, cosmo, psf = cluster_mass_selection
    cluster_m.set_ipurity(ipurity)
    cluster_m.set_completeness(completeness)

    lnM = np.linspace(np.log(1e13), np.log(1e16), nsize)
    z = np.linspace(0, 1.1, nsize)

    my_globals = globals().copy()
    my_globals.update(
        {
            "nsize": nsize,
            "lnM": lnM,
            "z": z,
            "np": np,
            "cluster_m": cluster_m,
            "cosmo": cosmo,
        }
    )

    cluster_intp_mass_test = (
        """[cluster_m.intp(cosmo, lnM[i], 0.5) for i in range(nsize)]"""
    )
    execution_time = timeit.timeit(
        cluster_intp_mass_test, globals=my_globals, number=100
    )
    print(
        f"Average time per execution cluster_m.intp mass: {execution_time / 100:.6f} seconds"
    )

    cluster_intp_z_test = (
        """[cluster_m.intp(cosmo, np.log(1e14), z[i]) for i in range(nsize)]"""
    )
    execution_time = timeit.timeit(cluster_intp_z_test, globals=my_globals, number=100)
    print(
        f"Average time per execution cluster_m.intp redshift: {execution_time / 100:.6f} seconds"
    )


def test_cluster_mass_selection_cumulative_bin(cluster_mass_selection):

    nsize = 100
    cluster_m, completeness, ipurity, cosmo, psf = cluster_mass_selection
    cluster_m.set_ipurity(ipurity)
    cluster_m.set_completeness(completeness)

    lnM = np.linspace(np.log(1e13), np.log(1e16), nsize)
    z = np.linspace(0, 1.1, nsize)
    lnR = np.linspace(np.log(5), np.log(200), nsize)

    my_globals = globals().copy()
    my_globals.update(
        {
            "nsize": nsize,
            "lnM": lnM,
            "z": z,
            "lnR": lnR,
            "np": np,
            "cluster_m": cluster_m,
            "cosmo": cosmo,
        }
    )

    cluster_intp_bin_mass_test = """[[cluster_m.intp_bin(cosmo, lnM[i], 0.5, [lnR[j]], [lnR[j+1]], None) for j in range(nsize-1)] for i in range(nsize)]"""
    execution_time = timeit.timeit(
        cluster_intp_bin_mass_test, globals=my_globals, number=100
    )
    print(
        f"Average time per execution cluster_m.intp_bin mass: {execution_time / 100:.6f} seconds"
    )

    cluster_intp_bin_redshift_test = """[[cluster_m.intp_bin(cosmo, np.log(1e14), z[i], [lnR[j]], [lnR[j+1]], None) for j in range(nsize-1)] for i in range(nsize)]"""
    execution_time = timeit.timeit(
        cluster_intp_bin_redshift_test, globals=my_globals, number=100
    )
    print(
        f"Average time per execution cluster_m.intp_bin redshift: {execution_time / 100:.6f} seconds"
    )


def test_cluster_mass_selection_limits(cluster_mass_selection):

    cluster_m, completeness, ipurity, cosmo, psf = cluster_mass_selection
    cluster_m.set_ipurity(ipurity)
    cluster_m.set_completeness(completeness)

    nsize = 100
    lnM = np.linspace(np.log(1e12), np.log(1e16), nsize)

    assert math.isclose(cluster_m.n_limits(cosmo)[0], lnM[0], rel_tol=1e-14)
    assert math.isclose(cluster_m.n_limits(cosmo)[1], lnM[-1], rel_tol=1e-14)

    assert math.isclose(
        cluster_m.p_limits(cosmo, [np.log(5)], [0])[0], lnM[0], rel_tol=1e-14
    )
    assert math.isclose(
        cluster_m.p_limits(cosmo, [np.log(5)], [0])[1], lnM[-1], rel_tol=1e-14
    )

    assert math.isclose(
        cluster_m.p_bin_limits(cosmo, [np.log(5)], [np.log(10)], [0])[0],
        lnM[0],
        rel_tol=1e-14,
    )
    assert math.isclose(
        cluster_m.p_bin_limits(cosmo, [np.log(5)], [np.log(10)], [0])[1],
        lnM[-1],
        rel_tol=1e-14,
    )


def test_cluster_mass_selection_resample(cluster_mass_selection):

    cluster_m, completeness, ipurity, cosmo, psf = cluster_mass_selection
    cluster_m.set_ipurity(ipurity)
    cluster_m.set_completeness(completeness)
    cluster_z = Nc.ClusterRedshiftNodist(z_min=0.1, z_max=1.1)

    nsize = 100
    lnM = np.linspace(np.log(1e13), np.log(1e16), nsize)
    z = np.linspace(0, 1.1, nsize)
    lnR = np.linspace(np.log(5), np.log(200), nsize)

    dist = Nc.Distance.new(2.0)
    dist.prepare(cosmo)

    mulf = Nc.MultiplicityFuncDespali.new()
    mulf.set_mdef(Nc.MultiplicityFuncMassDef.VIRIAL)
    hmf = Nc.HaloMassFunction.new(dist, psf, mulf)
    hmf.prepare(cosmo)
    hmf.set_area(439.790 * (np.pi / 180) ** 2)
    hbias = Nc.HaloBiasTinker.new(hmf)
    cad = Nc.ClusterAbundance.new(hmf, hbias)
    cad.set_area(439.790 * (np.pi / 180) ** 2)
    cad.prepare(cosmo, cluster_z, cluster_m)

    mset = Ncm.MSet.new_array([cosmo, cluster_z, cluster_m])
    rng = Ncm.RNG.seeded_new(None, 42)

    assert cluster_m.get_enable_rejection() == True
    ncount_rejection = Nc.DataClusterNCount.new(
        cad, "NcClusterRedshiftNodist", "NcClusterMassSelection"
    )
    ncount_rejection.init_from_sampling(mset, 439.790 * ((np.pi / 180) ** 2), rng)

    cluster_m.set_enable_rejection(False)
    assert cluster_m.get_enable_rejection() == False
    ncount_no_rejection = Nc.DataClusterNCount.new(
        cad, "NcClusterRedshiftNodist", "NcClusterMassSelection"
    )
    ncount_no_rejection.init_from_sampling(mset, 439.790 * ((np.pi / 180) ** 2), rng)

    assert (
        ncount_rejection.get_lnM_true().len() < ncount_no_rejection.get_lnM_true().len()
    )
    assert ncount_rejection.get_z_true().len() < ncount_no_rejection.get_z_true().len()

    assert (
        ncount_rejection.get_lnM_obs().col_len()
        < ncount_no_rejection.get_lnM_obs().col_len()
    )
    assert (
        ncount_rejection.get_z_obs().col_len()
        < ncount_no_rejection.get_z_obs().col_len()
    )


def test_cluster_mass_selection_hmf(cluster_mass_selection):

    cluster_m, completeness, ipurity, cosmo, psf = cluster_mass_selection
    cluster_m.set_ipurity(ipurity)
    cluster_m.set_completeness(completeness)
    cluster_z = Nc.ClusterRedshiftNodist(z_min=0.1, z_max=1.1)

    nsize = 100
    lnM = np.linspace(np.log(1e13), np.log(1e16), nsize)
    z = np.linspace(0, 1.1, nsize)
    lnR = np.linspace(np.log(5), np.log(200), nsize)

    dist = Nc.Distance.new(2.0)
    dist.prepare(cosmo)

    mulf = Nc.MultiplicityFuncDespali.new()
    mulf.set_mdef(Nc.MultiplicityFuncMassDef.VIRIAL)
    hmf = Nc.HaloMassFunction.new(dist, psf, mulf)
    hmf.prepare(cosmo)
    hmf.set_area(439.790 * (np.pi / 180) ** 2)
    hbias = Nc.HaloBiasTinker.new(hmf)
    cad = Nc.ClusterAbundance.new(hmf, hbias)
    cad.set_area(439.790 * (np.pi / 180) ** 2)
    cad.prepare(cosmo, cluster_z, cluster_m)

    my_globals = globals().copy()
    my_globals.update(
        {
            "nsize": nsize,
            "lnM": lnM,
            "z": z,
            "lnR": lnR,
            "np": np,
            "cluster_m": cluster_m,
            "cosmo": cosmo,
            "cluster_z": cluster_z,
            "cad": cad,
        }
    )

    cad_n_test = """cad.n(cosmo ,cluster_z , cluster_m)"""
    execution_time = timeit.timeit(cad_n_test, globals=my_globals, number=100)
    print(
        f"Average time per execution cad.n with selection: {execution_time / 100:.6f} seconds"
    )

    cad_d2n_richness_test = (
        """[cad.d2n(cosmo ,cluster_z , cluster_m,lnM[i],0.5) for i in range(nsize)]"""
    )
    execution_time = timeit.timeit(
        cad_d2n_richness_test, globals=my_globals, number=100
    )
    print(
        f"Average time per execution cad.d2n richness with selection : {execution_time / 100:.6f} seconds"
    )

    cad_d2n_redshift_test = """[cad.d2n(cosmo ,cluster_z , cluster_m,np.log(1e13),z[i]) for i in range(nsize)]"""
    execution_time = timeit.timeit(
        cad_d2n_redshift_test, globals=my_globals, number=100
    )
    print(
        f"Average time per execution cad.d2n redshift with selection: {execution_time / 100:.6f} seconds"
    )

    cad_intp_bin_d2n_richness_test = """[cad.intp_bin_d2n(cosmo ,cluster_z , cluster_m,[lnR[i]], [lnR[i+1]], None ,[0.1], [1.1], None) for i in range(nsize-1)]"""
    execution_time = timeit.timeit(
        cad_intp_bin_d2n_richness_test, globals=my_globals, number=1
    )
    print(
        f"Average time per execution cad.intp_bin_d2n richness with selection : {execution_time / 1:.6f} seconds"
    )

    cad_intp_bin_d2n_redshift_test = """[cad.intp_bin_d2n(cosmo ,cluster_z , cluster_m,[np.log(5)], [np.log(200)], None ,[z[i]], [z[i+1]], None)  for i in range(nsize-1)]"""
    execution_time = timeit.timeit(
        cad_intp_bin_d2n_redshift_test, globals=my_globals, number=1
    )
    print(
        f"Average time per execution cad.intp_bin_d2n redshift with selection: {execution_time / 1:.6f} seconds"
    )

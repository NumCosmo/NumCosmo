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

import timeit
import math
import pytest
import numpy as np
from numcosmo_py import Nc, Ncm
from numcosmo_py.helper import npa_to_seq

Ncm.cfg_init()

# Constants
SURVEY_AREA_DEG2 = 439.790
SURVEY_AREA_RAD2 = SURVEY_AREA_DEG2 * (np.pi / 180) ** 2
LN_RICHNESS_CUT = np.log(5)


@pytest.fixture(name="prim")
def fixture_prim() -> Nc.HIPrim:
    """Fixture for HIPrimPowerLaw."""
    prim = Nc.HIPrimPowerLaw.new()
    prim.props.n_SA = 0.967
    return prim


@pytest.fixture(name="reion")
def fixture_reion() -> Nc.HIReion:
    """Fixture for HIReionCamb."""
    return Nc.HIReionCamb.new()


@pytest.fixture(name="cosmo")
def fixture_cosmo(prim: Nc.HIPrim, reion: Nc.HIReion) -> Nc.HICosmo:
    """Create and configure a cosmological model."""
    cosmo = Nc.HICosmoDEXcdm()
    cosmo.param_set_by_name("Omegax", 1 - 0.2603)
    cosmo.param_set_by_name("H0", 71)
    cosmo.param_set_by_name("Omegab", 0.0406)
    cosmo.param_set_by_name("Omegac", 0.22)
    cosmo.param_set_by_name("w", -1.0)

    cosmo.add_submodel(prim)
    cosmo.add_submodel(reion)

    return cosmo


@pytest.fixture(name="psf")
def fixture_psf(cosmo: Nc.HICosmo, prim: Nc.HIPrim) -> Ncm.PowspecFilter:
    """Create and configure power spectrum filter."""
    tf = Nc.TransferFuncEH()
    psml = Nc.PowspecMLTransfer.new(tf)
    psml.require_kmin(1.0e-6)
    psml.require_kmax(1.0e3)

    psf = Ncm.PowspecFilter.new(psml, Ncm.PowspecFilterType.TOPHAT)
    psf.set_best_lnr0()
    psf.prepare(cosmo)

    # Normalize to sigma8 = 0.8
    old_amplitude = np.exp(prim["ln10e10ASA"])
    prim["ln10e10ASA"] = np.log((0.8 / cosmo.sigma8(psf)) ** 2 * old_amplitude)

    return psf


@pytest.fixture(name="cluster_m")
def fixture_cluster_m() -> Nc.ClusterMass:
    """Create and configure cluster mass selection model."""
    cluster_m = Nc.ClusterMassSelection(
        lnRichness_min=LN_RICHNESS_CUT, lnRichness_max=np.log(200)
    )
    cluster_m.param_set_by_name("mup0", 4.12769558168741)
    cluster_m.param_set_by_name("mup1", 1.17476066603899)
    cluster_m.param_set_by_name("mup2", 0.393577193825473)
    cluster_m.param_set_by_name("sigmap0", 0.408750324989284)
    cluster_m.param_set_by_name("sigmap1", -0.123232985316648)
    cluster_m.param_set_by_name("sigmap2", -0.0644996574273048)
    cluster_m.param_set_by_name("cut", LN_RICHNESS_CUT)

    return cluster_m


@pytest.fixture(name="completeness")
def fixture_completeness() -> Ncm.Spline2d:
    """Create completeness 2D spline from data."""
    lnM_bins_knots = np.linspace(13.0 * np.log(10), 16 * np.log(10), 10)
    z_bins_knots = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0, 1.1])
    lnM_centers = 0.5 * (lnM_bins_knots[:-1] + lnM_bins_knots[1:])
    z_centers = 0.5 * (z_bins_knots[:-1] + z_bins_knots[1:])

    completeness_data = np.array(
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

    completeness_knots = Ncm.Matrix.new(len(z_centers), len(lnM_centers))
    for i in range(len(z_centers)):
        for j in range(len(lnM_centers)):
            completeness_knots.set(i, j, completeness_data.T[i][j])

    return Ncm.Spline2dBicubic(
        spline=Ncm.SplineCubicNotaknot.new(),
        x_vector=Ncm.Vector.new_array(npa_to_seq(lnM_centers)),
        y_vector=Ncm.Vector.new_array(npa_to_seq(z_centers)),
        z_matrix=completeness_knots,
    )


@pytest.fixture(name="ipurity")
def fixture_ipurity() -> Ncm.Spline2d:
    """Create inverse purity 2D spline from data."""
    lnR_bins_knots = np.log(np.array([5, 10, 15, 20, 35, 70, 100, 200]))
    z_bins_knots = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0, 1.1])
    lnR_centers = 0.5 * (lnR_bins_knots[:-1] + lnR_bins_knots[1:])
    z_centers = 0.5 * (z_bins_knots[:-1] + z_bins_knots[1:])

    ipurity_data = np.array(
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

    ipurity_knots = Ncm.Matrix.new(len(z_centers), len(lnR_centers))
    for i in range(len(z_centers)):
        for j in range(len(lnR_centers)):
            ipurity_knots.set(i, j, ipurity_data.T[i][j])

    return Ncm.Spline2dBicubic(
        spline=Ncm.SplineCubicNotaknot.new(),
        x_vector=Ncm.Vector.new_array(npa_to_seq(lnR_centers)),
        y_vector=Ncm.Vector.new_array(npa_to_seq(z_centers)),
        z_matrix=ipurity_knots,
    )


def _benchmark_function(func_str: str, globals_dict: dict, number: int = 100):
    """Execute and time a function, printing the average execution time."""
    execution_time = timeit.timeit(func_str, globals=globals_dict, number=number)
    return execution_time / number


@pytest.fixture(name="cluster_z")
def fixture_cluster_z():
    """Create and configure cluster redshift object."""
    return Nc.ClusterRedshiftNodist(z_min=0.1, z_max=1.1)


@pytest.fixture(name="cluster_abundance")
def fixture_cluster_abundance(
    cosmo: Nc.HICosmo,
    psf: Ncm.PowspecFilter,
    cluster_z: Nc.ClusterRedshift,
    cluster_m: Nc.ClusterMass,
):
    """Create and configure cluster abundance object."""
    dist = Nc.Distance.new(2.0)
    dist.prepare(cosmo)

    mulf = Nc.MultiplicityFuncDespali.new()
    mulf.set_mdef(Nc.MultiplicityFuncMassDef.VIRIAL)

    hmf = Nc.HaloMassFunction.new(dist, psf, mulf)
    hmf.prepare(cosmo)
    hmf.set_area(SURVEY_AREA_RAD2)

    hbias = Nc.HaloBiasTinker.new(hmf)

    cad = Nc.ClusterAbundance.new(hmf, hbias)
    cad.set_area(SURVEY_AREA_RAD2)
    cad.prepare(cosmo, cluster_z, cluster_m)

    return cad


def test_get_cut(cluster_m: Nc.ClusterMass) -> None:
    """Fixture providing configured cluster mass selection model and dependencies."""
    assert cluster_m.get_cut() == LN_RICHNESS_CUT


def test_cluster_mass_selection_completeness_and_purity(
    cluster_m: Nc.ClusterMass, completeness: Ncm.Spline2d, ipurity: Ncm.Spline2d
) -> None:
    """Test completeness and purity functions with and without splines."""
    nsize = 500
    lnM = np.linspace(np.log(1e12), np.log(1e16), nsize)
    z = np.linspace(0, 2, nsize)
    lnR = np.linspace(np.log(1), np.log(1000), nsize)

    # Initially, no completeness/ipurity should be set
    assert cluster_m.peek_completeness() is None
    assert cluster_m.peek_ipurity() is None

    # Without splines, completeness and ipurity should be 1.0
    for i in range(nsize):
        assert cluster_m.completeness(lnM[i], 0.5) == 1.0
        assert cluster_m.completeness(np.log(1e14), z[i]) == 1.0
        assert cluster_m.ipurity(np.log(5), z[i]) == 1.0
        assert cluster_m.ipurity(lnR[i], 0.5) == 1.0

    # Set splines and verify they are properly assigned
    cluster_m.set_ipurity(ipurity)
    cluster_m.set_completeness(completeness)

    assert cluster_m.peek_completeness() is completeness
    assert cluster_m.peek_ipurity() is ipurity


def test_cluster_mass_selection_mean_std(
    cluster_m: Nc.ClusterMass, completeness: Ncm.Spline2d, ipurity: Ncm.Spline2d
) -> None:
    """Test mean and standard deviation with and without truncation."""
    nsize = 100
    cluster_m.set_ipurity(ipurity)
    cluster_m.set_completeness(completeness)

    lnM = np.linspace(np.log(1e13), np.log(1e16), nsize)
    z = np.linspace(0, 1.1, nsize)

    # With truncation, mean should be >= untruncated mean, std should be <= untruncated
    # std
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

    # Without truncation (cut at -infinity), mean and std should match
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


def test_cluster_mass_selection_distribution(
    cluster_m: Nc.ClusterMass,
    completeness: Ncm.Spline2d,
    ipurity: Ncm.Spline2d,
    cosmo: Nc.HICosmo,
) -> None:
    """Benchmark probability distribution function.

    Benchmark across mass, redshift, and richness."""
    nsize = 100
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

    # Benchmark p() function for different parameters
    tests = [
        (
            "mass",
            "[cluster_m.p(cosmo, lnM[i], 0.5, [np.log(20)], None) "
            "for i in range(nsize)]",
        ),
        (
            "redshift",
            "[cluster_m.p(cosmo, np.log(1e14), z[i], [np.log(20)], None) "
            "for i in range(nsize)]",
        ),
        (
            "richness",
            "[cluster_m.p(cosmo, np.log(1e14), 0.5, [lnR[i]], None) "
            "for i in range(nsize)]",
        ),
    ]

    for name, test_str in tests:
        avg_time = _benchmark_function(test_str, my_globals)
        print(f"Average time per execution cluster_m.p {name}: {avg_time:.6f} seconds")


def test_cluster_mass_selection_cumulative(
    cluster_m: Nc.ClusterMass,
    completeness: Ncm.Spline2d,
    ipurity: Ncm.Spline2d,
    cosmo: Nc.HICosmo,
) -> None:
    """Benchmark cumulative distribution function across mass and redshift."""
    nsize = 100
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

    # Benchmark intp() function for different parameters
    tests = [
        ("mass", "[cluster_m.intp(cosmo, lnM[i], 0.5) for i in range(nsize)]"),
        (
            "redshift",
            "[cluster_m.intp(cosmo, np.log(1e14), z[i]) for i in range(nsize)]",
        ),
    ]

    for name, test_str in tests:
        avg_time = _benchmark_function(test_str, my_globals)
        print(
            f"Average time per execution cluster_m.intp {name}: {avg_time:.6f} seconds"
        )


def test_cluster_mass_selection_cumulative_bin(
    cluster_m: Nc.ClusterMass,
    completeness: Ncm.Spline2d,
    ipurity: Ncm.Spline2d,
    cosmo: Nc.HICosmo,
) -> None:
    """Benchmark binned cumulative distribution function across mass and redshift."""
    nsize = 100
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

    # Benchmark intp_bin() function for different parameters
    tests = [
        (
            "mass",
            "[[cluster_m.intp_bin(cosmo, lnM[i], 0.5, [lnR[j]], [lnR[j+1]], None) "
            "for j in range(nsize-1)] for i in range(nsize)]",
        ),
        (
            "redshift",
            "[[cluster_m.intp_bin(cosmo, np.log(1e14), z[i], [lnR[j]], [lnR[j+1]], "
            "None) for j in range(nsize-1)] for i in range(nsize)]",
        ),
    ]

    for name, test_str in tests:
        avg_time = _benchmark_function(test_str, my_globals)
        print(
            f"Average time per execution cluster_m.intp_bin {name}: "
            f"{avg_time:.6f} seconds"
        )


def test_cluster_mass_selection_limits(
    cluster_m: Nc.ClusterMass,
    completeness: Ncm.Spline2d,
    ipurity: Ncm.Spline2d,
    cosmo: Nc.HICosmo,
) -> None:
    """Test mass limits for different selection functions."""
    cluster_m.set_ipurity(ipurity)
    cluster_m.set_completeness(completeness)

    nsize = 100
    lnM = np.linspace(np.log(1e12), np.log(1e16), nsize)

    # Test n_limits
    assert math.isclose(cluster_m.n_limits(cosmo)[0], lnM[0], rel_tol=1e-14)
    assert math.isclose(cluster_m.n_limits(cosmo)[1], lnM[-1], rel_tol=1e-14)

    # Test p_limits
    assert math.isclose(
        cluster_m.p_limits(cosmo, [np.log(5)], [0])[0], lnM[0], rel_tol=1e-14
    )
    assert math.isclose(
        cluster_m.p_limits(cosmo, [np.log(5)], [0])[1], lnM[-1], rel_tol=1e-14
    )

    # Test p_bin_limits
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


def test_cluster_mass_selection_resample(
    cluster_m: Nc.ClusterMass,
    completeness: Ncm.Spline2d,
    ipurity: Ncm.Spline2d,
    cosmo: Nc.HICosmo,
    cluster_abundance: Nc.ClusterAbundance,
) -> None:
    """Test resampling with and without rejection."""
    cluster_m.set_ipurity(ipurity)
    cluster_m.set_completeness(completeness)
    cluster_z = Nc.ClusterRedshiftNodist(z_min=0.1, z_max=1.1)

    cad = cluster_abundance
    mset = Ncm.MSet.new_array([cosmo, cluster_z, cluster_m])
    rng = Ncm.RNG.seeded_new(None, 42)

    # Test with rejection enabled
    assert cluster_m.get_enable_rejection()
    ncount_rejection = Nc.DataClusterNCount.new(
        cad, "NcClusterRedshiftNodist", "NcClusterMassSelection"
    )
    ncount_rejection.init_from_sampling(mset, SURVEY_AREA_RAD2, rng)

    # Test with rejection disabled
    cluster_m.set_enable_rejection(False)
    assert not cluster_m.get_enable_rejection()
    ncount_no_rejection = Nc.DataClusterNCount.new(
        cad, "NcClusterRedshiftNodist", "NcClusterMassSelection"
    )
    ncount_no_rejection.init_from_sampling(mset, SURVEY_AREA_RAD2, rng)

    # Rejection should produce fewer clusters
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


def test_cluster_mass_selection_hmf(
    cluster_m: Nc.ClusterMass,
    completeness: Ncm.Spline2d,
    ipurity: Ncm.Spline2d,
    cosmo: Nc.HICosmo,
    cluster_abundance: Nc.ClusterAbundance,
) -> None:
    """Benchmark cluster abundance calculations with selection effects."""
    cluster_m.set_ipurity(ipurity)
    cluster_m.set_completeness(completeness)
    cluster_z = Nc.ClusterRedshiftNodist(z_min=0.1, z_max=1.1)

    nsize = 100
    lnM = np.linspace(np.log(1e13), np.log(1e16), nsize)
    z = np.linspace(0, 1.1, nsize)
    lnR = np.linspace(np.log(5), np.log(200), nsize)

    cad = cluster_abundance

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

    # Benchmark cluster abundance functions
    tests = [
        ("cad.n", "cad.n(cosmo, cluster_z, cluster_m)", 100),
        (
            "cad.d2n mass",
            "[cad.d2n(cosmo, cluster_z, cluster_m, lnM[i], 0.5) "
            "for i in range(nsize)]",
            100,
        ),
        (
            "cad.d2n redshift",
            "[cad.d2n(cosmo, cluster_z, cluster_m, np.log(1e13), "
            "z[i]) for i in range(nsize)]",
            100,
        ),
        (
            "cad.intp_bin_d2n richness",
            "[cad.intp_bin_d2n(cosmo, cluster_z, cluster_m, [lnR[i]], "
            "[lnR[i+1]], None, [0.1], [1.1], None) for i in range(nsize-1)]",
            1,
        ),
        (
            "cad.intp_bin_d2n redshift",
            "[cad.intp_bin_d2n(cosmo, cluster_z, cluster_m, [np.log(5)], "
            "[np.log(200)], None, [z[i]], [z[i+1]], None) for i in range(nsize-1)]",
            1,
        ),
    ]

    for name, test_str, number in tests:
        avg_time = _benchmark_function(test_str, my_globals, number)
        print(
            f"Average time per execution {name} with selection: {avg_time:.6f} seconds"
        )

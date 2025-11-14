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

from timeit import timeit
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


@pytest.fixture(name="cluster_m_empty")
def fixture_cluster_m_empty() -> Nc.ClusterMassSelection:
    """Create cluster mass-richness relation model.

    Configures a selection model with richness cut at ln(5).
    """
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
    """Create completeness function as 2D spline.

    Returns a bicubic spline interpolating completeness as a function of
    cluster mass and redshift.
    """
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
    """Create inverse purity function as 2D spline.

    Returns a bicubic spline interpolating inverse purity as a function of
    richness and redshift.
    """
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


@pytest.fixture(name="cluster_m")
def fixture_cluster_m(
    cluster_m_empty: Nc.ClusterMassSelection,
    completeness: Ncm.Spline2d,
    ipurity: Ncm.Spline2d,
) -> Nc.ClusterMassSelection:
    """Create cluster mass-richness relation model.

    Configures a selection model with richness cut at ln(5). Includes
    completeness and ipurity splines.
    """
    assert isinstance(completeness, Ncm.Spline2dBicubic)
    assert isinstance(ipurity, Ncm.Spline2dBicubic)

    cluster_m_empty.set_completeness(completeness)
    cluster_m_empty.set_ipurity(ipurity)
    return cluster_m_empty


@pytest.fixture(name="cluster_z")
def fixture_cluster_z() -> Nc.ClusterRedshiftNodist:
    """Create cluster redshift distribution.

    Returns a no-distribution model with redshift range [0.1, 1.1].
    """
    return Nc.ClusterRedshiftNodist(z_min=0.1, z_max=1.1)


@pytest.fixture(name="cluster_abundance")
def fixture_cluster_abundance(
    cosmo: Nc.HICosmo,
    psf: Ncm.PowspecFilter,
    cluster_z: Nc.ClusterRedshiftNodist,
    cluster_m: Nc.ClusterMassSelection,
):
    """Create cluster abundance calculator.

    Configures halo mass function with Despali multiplicity function and
    Tinker bias for the survey area.
    """
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


def test_get_cut(cluster_m: Nc.ClusterMassSelection) -> None:
    """Test richness cut parameter."""
    assert cluster_m.get_cut() == LN_RICHNESS_CUT


def test_cluster_mass_selection_completeness_and_purity(
    cluster_m_empty: Nc.ClusterMassSelection,
    completeness: Ncm.Spline2d,
    ipurity: Ncm.Spline2d,
) -> None:
    """Test completeness and purity functions.

    Verifies that without splines, completeness and purity default to 1.0,
    and that splines can be properly set and retrieved.
    """
    nsize = 500
    lnM = np.linspace(np.log(1e12), np.log(1e16), nsize)
    z = np.linspace(0, 2, nsize)
    lnR = np.linspace(np.log(1), np.log(1000), nsize)

    # Initially, no completeness/ipurity should be set
    assert cluster_m_empty.peek_completeness() is None
    assert cluster_m_empty.peek_ipurity() is None

    # Without splines, completeness and ipurity should be 1.0
    for i in range(nsize):
        assert cluster_m_empty.completeness(lnM[i], 0.5) == 1.0
        assert cluster_m_empty.completeness(np.log(1e14), z[i]) == 1.0
        assert cluster_m_empty.ipurity(np.log(5), z[i]) == 1.0
        assert cluster_m_empty.ipurity(lnR[i], 0.5) == 1.0

    # Set splines and verify they are properly assigned
    assert isinstance(completeness, Ncm.Spline2dBicubic)
    assert isinstance(ipurity, Ncm.Spline2dBicubic)
    cluster_m_empty.set_ipurity(ipurity)
    cluster_m_empty.set_completeness(completeness)

    assert cluster_m_empty.peek_completeness() is completeness
    assert cluster_m_empty.peek_ipurity() is ipurity


def test_cluster_mass_selection_mean_std(cluster_m: Nc.ClusterMassSelection) -> None:
    """Test mean and standard deviation of truncated distribution.

    Verifies that truncation increases mean and decreases standard deviation
    compared to the untruncated Gaussian.
    """
    nsize = 100
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
    cluster_m: Nc.ClusterMassSelection, cosmo: Nc.HICosmo
) -> None:
    """Benchmark probability distribution function.

    Times the evaluation of p() across mass, redshift, and richness dimensions.
    """
    nsize = 100

    lnM = np.linspace(np.log(1e13), np.log(1e16), nsize)
    z = np.linspace(0, 1.1, nsize)
    lnR = np.linspace(np.log(5), np.log(200), nsize)

    tests = [
        (
            "mass",
            lambda: [
                cluster_m.p(cosmo, lnM[i], 0.5, [np.log(20)], None)
                for i in range(nsize)
            ],
        ),
        (
            "redshift",
            lambda: [
                cluster_m.p(cosmo, np.log(1e14), z[i], [np.log(20)], None)
                for i in range(nsize)
            ],
        ),
        (
            "richness",
            lambda: [
                cluster_m.p(cosmo, np.log(1e14), 0.5, [lnR[i]], None)
                for i in range(nsize)
            ],
        ),
    ]

    for name, func in tests:
        avg_time = timeit(func, number=100)
        print(f"Average time per execution cluster_m.p {name}: {avg_time:.6f} seconds")


def test_cluster_mass_selection_cumulative(
    cluster_m: Nc.ClusterMassSelection, cosmo: Nc.HICosmo
) -> None:
    """Benchmark cumulative distribution function across mass and redshift."""
    nsize = 100

    lnM = np.linspace(np.log(1e13), np.log(1e16), nsize)
    z = np.linspace(0, 1.1, nsize)

    tests = [
        ("mass", lambda: [cluster_m.intp(cosmo, lnM[i], 0.5) for i in range(nsize)]),
        (
            "redshift",
            lambda: [cluster_m.intp(cosmo, np.log(1e14), z[i]) for i in range(nsize)],
        ),
    ]

    for name, func in tests:
        avg_time = timeit(func, number=100)
        print(
            f"Average time per execution cluster_m.intp {name}: {avg_time:.6f} seconds"
        )


def test_cluster_mass_selection_cumulative_bin(
    cluster_m: Nc.ClusterMassSelection, cosmo: Nc.HICosmo
) -> None:
    """Benchmark binned cumulative distribution function across mass and redshift."""
    nsize = 100
    lnM = np.linspace(np.log(1e13), np.log(1e16), nsize)
    z = np.linspace(0, 1.1, nsize)
    lnR = np.linspace(np.log(5), np.log(200), nsize)

    tests = [
        (
            "mass",
            lambda: [
                [
                    cluster_m.intp_bin(cosmo, lnM[i], 0.5, [lnR[j]], [lnR[j + 1]], None)
                    for j in range(nsize - 1)
                ]
                for i in range(nsize)
            ],
        ),
        (
            "redshift",
            lambda: [
                [
                    cluster_m.intp_bin(
                        cosmo, np.log(1e14), z[i], [lnR[j]], [lnR[j + 1]], None
                    )
                    for j in range(nsize - 1)
                ]
                for i in range(nsize)
            ],
        ),
    ]

    for name, func in tests:
        avg_time = timeit(func, number=100)
        print(
            f"Average time per execution cluster_m.intp_bin {name}: "
            f"{avg_time:.6f} seconds"
        )


def test_cluster_mass_selection_limits(
    cluster_m: Nc.ClusterMassSelection, cosmo: Nc.HICosmo
) -> None:
    """Test mass limits for different selection functions."""
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
    cluster_m: Nc.ClusterMassSelection,
    cosmo: Nc.HICosmo,
    cluster_abundance: Nc.ClusterAbundance,
) -> None:
    """Test cluster resampling with rejection.

    Verifies that rejection sampling produces fewer clusters than sampling
    without rejection.
    """
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
    cluster_m: Nc.ClusterMassSelection,
    cosmo: Nc.HICosmo,
    cluster_abundance: Nc.ClusterAbundance,
) -> None:
    """Benchmark cluster abundance calculations.

    Times the evaluation of cluster number counts and differential number
    counts with selection effects applied.
    """
    assert isinstance(cluster_m, Nc.ClusterMassSelection)

    cluster_z = Nc.ClusterRedshiftNodist(z_min=0.1, z_max=1.1)

    nsize = 100
    lnM = np.linspace(np.log(1e13), np.log(1e16), nsize)
    z = np.linspace(0, 1.1, nsize)
    lnR = np.linspace(np.log(5), np.log(200), nsize)

    cad = cluster_abundance

    tests = [
        ("cad.n", lambda: cad.n(cosmo, cluster_z, cluster_m), 100),
        (
            "cad.d2n mass",
            lambda: [
                cad.d2n(cosmo, cluster_z, cluster_m, lnM[i], 0.5) for i in range(nsize)
            ],
            100,
        ),
        (
            "cad.d2n redshift",
            lambda: [
                cad.d2n(cosmo, cluster_z, cluster_m, np.log(1e13), z[i])
                for i in range(nsize)
            ],
            100,
        ),
        (
            "cad.intp_bin_d2n richness",
            lambda: [
                cad.intp_bin_d2n(
                    cosmo,
                    cluster_z,
                    cluster_m,
                    [lnR[i]],
                    [lnR[i + 1]],
                    None,
                    [0.1],
                    [1.1],
                    None,
                )
                for i in range(nsize - 1)
            ],
            1,
        ),
        (
            "cad.intp_bin_d2n redshift",
            lambda: [
                cad.intp_bin_d2n(
                    cosmo,
                    cluster_z,
                    cluster_m,
                    [np.log(5)],
                    [np.log(200)],
                    None,
                    [z[i]],
                    [z[i + 1]],
                    None,
                )
                for i in range(nsize - 1)
            ],
            1,
        ),
    ]

    for name, func, number in tests:
        avg_time = timeit(func, number=number)
        print(
            f"Average time per execution {name} with selection: {avg_time:.6f} seconds"
        )


def test_cluster_mass_selection_properties(
    cluster_m: Nc.ClusterMassSelection,
    completeness: Ncm.Spline2d,
    ipurity: Ncm.Spline2d,
) -> None:
    """Test setting properties via GObject interface."""
    assert isinstance(cluster_m, Nc.ClusterMassSelection)
    cluster_m.props.enable_rejection = False
    assert not cluster_m.props.enable_rejection

    cluster_m.props.enable_rejection = True
    assert cluster_m.props.enable_rejection

    cluster_m.props.ipurity = ipurity
    assert cluster_m.props.ipurity is ipurity

    cluster_m.props.completeness = completeness
    assert cluster_m.props.completeness is completeness


def test_serialization_deserialization(cluster_m: Nc.ClusterMassSelection) -> None:
    """Test serialization and deserialization."""
    assert isinstance(cluster_m, Nc.ClusterMassSelection)
    cluster_m.set_enable_rejection(True)

    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
    cluster_m_dup_obj = ser.dup_obj(cluster_m)
    assert isinstance(cluster_m_dup_obj, Nc.ClusterMassSelection)
    cluster_m_dup: Nc.ClusterMassSelection = cluster_m_dup_obj

    assert cluster_m_dup is not None
    assert isinstance(cluster_m_dup, Nc.ClusterMassSelection)
    assert cluster_m_dup is not cluster_m

    assert cluster_m_dup.get_enable_rejection() == cluster_m.get_enable_rejection()
    assert cluster_m_dup.peek_ipurity() is not None
    assert cluster_m_dup.peek_completeness() is not None


def test_cluster_mass_selection_p_basic(
    cluster_m: Nc.ClusterMassSelection, cosmo: Nc.HICosmo
) -> None:
    """Test probability distribution function.

    Verifies that p() returns positive probabilities for valid inputs,
    decreases with distance from mean richness, and integrates to the
    completeness factor.
    """
    lnM = np.log(1e14)
    z = 0.5
    mean_lnR = cluster_m.get_mean_richness(lnM, z)

    # Probability at mean should be positive
    p_mean = cluster_m.p(cosmo, lnM, z, [mean_lnR], None)
    assert p_mean > 0.0

    # Probability should decrease away from mean
    p_above = cluster_m.p(cosmo, lnM, z, [mean_lnR + 1.0], None)
    p_below = cluster_m.p(cosmo, lnM, z, [mean_lnR - 1.0], None)
    assert p_above < p_mean
    assert p_below < p_mean

    # Probability should be very small far from mean
    p_far = cluster_m.p(cosmo, lnM, z, [mean_lnR + 5.0], None)
    assert p_far < p_mean * 0.01


def test_cluster_mass_selection_p_integral(
    cluster_m: Nc.ClusterMassSelection, cosmo: Nc.HICosmo
) -> None:
    """Test probability distribution integrates to completeness.

    Verifies that integrating p() over richness equals the completeness
    factor using trapezoidal rule.
    """
    lnM = np.log(1e14)
    z = 0.5
    mean_lnR = cluster_m.get_mean_richness(lnM, z)
    std_lnR = cluster_m.get_std_richness(lnM, z)

    # Integrate over ±5 sigma range
    lnR_array = np.linspace(mean_lnR - 5 * std_lnR, mean_lnR + 5 * std_lnR, 200)
    p_array = np.array([cluster_m.p(cosmo, lnM, z, [lnR], None) for lnR in lnR_array])

    # Trapezoidal integration
    integral = np.trapezoid(p_array, lnR_array)
    assert integral < 1.0
    intp_val = cluster_m.intp(cosmo, lnM, z)

    assert np.isclose(integral, intp_val, rtol=1e-3)


def test_cluster_mass_selection_intp_basic(
    cluster_m: Nc.ClusterMassSelection, cosmo: Nc.HICosmo
) -> None:
    """Test cumulative distribution function.

    Verifies that intp() returns values between 0 and 1.
    """
    lnM = np.log(1e14)
    z = 0.5

    # intp should be between 0 and 1
    intp_val = cluster_m.intp(cosmo, lnM, z)
    assert 0.0 <= intp_val <= 1.0


def test_cluster_mass_selection_intp_bin_basic(
    cluster_m: Nc.ClusterMassSelection, cosmo: Nc.HICosmo
) -> None:
    """Test binned integration of probability distribution.

    Verifies that intp_bin() returns positive values and that
    summing bins recovers the total integral.
    """
    lnM = np.log(1e14)
    z = 0.5
    mean_lnR = cluster_m.get_mean_richness(lnM, z)
    std_lnR = cluster_m.get_std_richness(lnM, z)

    # Single bin should be positive
    bin_val = cluster_m.intp_bin(cosmo, lnM, z, [mean_lnR], [mean_lnR + std_lnR], None)
    assert bin_val > 0.0

    # Sum of bins should equal total integral
    lnR_edges = np.linspace(mean_lnR - 5 * std_lnR, mean_lnR + 5 * std_lnR, 11)
    bin_sum = sum(
        cluster_m.intp_bin(cosmo, lnM, z, [lnR_edges[i]], [lnR_edges[i + 1]], None)
        for i in range(len(lnR_edges) - 1)
    )
    total_integral = cluster_m.intp(cosmo, lnM, z)

    assert np.isclose(bin_sum, total_integral, rtol=1e-3)


def test_cluster_mass_selection_intp_bin_consistency(
    cluster_m: Nc.ClusterMassSelection, cosmo: Nc.HICosmo
) -> None:
    """Test intp_bin() matches manual integration of p().

    Verifies that intp_bin() over an interval matches trapezoidal
    integration of p() over the same interval.
    """
    lnM = np.log(1e14)
    z = 0.5
    mean_lnR = cluster_m.get_mean_richness(lnM, z)
    std_lnR = cluster_m.get_std_richness(lnM, z)

    lnR_lower = mean_lnR - std_lnR
    lnR_upper = mean_lnR + std_lnR

    # Manual integration
    lnR_array = np.linspace(lnR_lower, lnR_upper, 200)
    p_array = np.array([cluster_m.p(cosmo, lnM, z, [lnR], None) for lnR in lnR_array])
    manual_integral = np.trapezoid(p_array, lnR_array)

    # Internal integration via intp_bin
    intp_bin_val = cluster_m.intp_bin(cosmo, lnM, z, [lnR_lower], [lnR_upper], None)

    assert np.isclose(manual_integral, intp_bin_val, rtol=1e-3)


def test_cluster_mass_selection_truncated_gaussian(
    cluster_m_empty: Nc.ClusterMassSelection,
) -> None:
    """Test truncated Gaussian mean and std corrections.

    Verifies that get_mean() and get_std() correctly account for
    truncation at the richness cut, while get_mean_richness() and
    get_std_richness() return the untruncated Gaussian parameters.
    """
    lnM = np.log(1e14)
    z = 0.5
    cut = cluster_m_empty.get_cut()

    # Untruncated parameters
    mean_untrunc = cluster_m_empty.get_mean_richness(lnM, z)
    std_untrunc = cluster_m_empty.get_std_richness(lnM, z)

    # Truncated parameters
    mean_trunc = cluster_m_empty.get_mean(lnM, z)
    std_trunc = cluster_m_empty.get_std(lnM, z)

    # Truncation should increase mean (shift right)
    assert mean_trunc > mean_untrunc

    # Truncation should decrease std (reduce spread)
    assert std_trunc < std_untrunc

    # Verify truncation formulas
    alpha = (cut - mean_untrunc) / std_untrunc
    phi_alpha = np.exp(-0.5 * alpha**2) / np.sqrt(2 * np.pi)
    Phi_alpha = 0.5 * (1 + math.erf(alpha / np.sqrt(2)))

    # Mean correction: μ_trunc = μ + σ * φ(α) / (1 - Φ(α))
    expected_mean = mean_untrunc + std_untrunc * phi_alpha / (1 - Phi_alpha)
    assert np.isclose(mean_trunc, expected_mean, rtol=1e-10)

    # Std correction: σ_trunc = σ * sqrt(1 + α*φ(α)/(1-Φ(α)) - (φ(α)/(1-Φ(α)))²)
    lambda_alpha = phi_alpha / (1 - Phi_alpha)
    expected_std = std_untrunc * np.sqrt(1 + alpha * lambda_alpha - lambda_alpha**2)
    assert np.isclose(std_trunc, expected_std, rtol=1e-10)


def test_cluster_mass_selection_p_vec_z_lnMobs(
    cluster_m: Nc.ClusterMassSelection, cosmo: Nc.HICosmo
) -> None:
    """Test vectorized probability distribution function.

    Verifies that p_vec_z_lnMobs() matches individual p() calls
    for multiple redshifts and richness values.
    """
    lnM = np.log(1e14)
    z_array = np.linspace(0.1, 1.1, 100)
    lnR_array = np.log(np.geomspace(3.0, 1e3, 100))

    # Create Ncm.Vector for z and Ncm.Matrix for lnR (one column)
    z_vec = Ncm.Vector.new_array(z_array)
    lnR_mat = Ncm.Matrix.new_array(lnR_array, 1)

    # Allocate result array
    res = Ncm.Vector.new(len(z_array))

    # Vectorized call
    cluster_m.p_vec_z_lnMobs(cosmo, lnM, z_vec, lnR_mat, None, res)

    # Compare with individual calls
    for i, (z, lnR) in enumerate(zip(z_array, lnR_array)):
        p_scalar = cluster_m.p(cosmo, lnM, z, [lnR], None)
        p_vec = res.get(i)
        assert np.isclose(p_scalar, p_vec, rtol=1e-10)

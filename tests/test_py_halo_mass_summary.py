#!/usr/bin/env python
#
# test_py_halo_mass_summary.py
#
# Sun Nov 10 19:52:17 2024
# Copyright  2024  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_py_halo_mass_summary.py
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

"""Tests on NcHaloMassSummary model."""

import pytest
from numpy.testing import assert_allclose
import numpy as np

from numcosmo_py import Ncm, Nc

Ncm.cfg_init()


@pytest.fixture(
    name="mass_def",
    params=[
        Nc.HaloMassSummaryMassDef.CRITICAL,
        Nc.HaloMassSummaryMassDef.MEAN,
        Nc.HaloMassSummaryMassDef.VIRIAL,
    ],
    ids=["critical", "mean", "virial"],
)
def fixture_mass_def(request) -> Nc.HaloMassSummaryMassDef:
    """Fixture for HaloMassSummaryMassDef."""
    return request.param


@pytest.fixture(
    name="Delta",
    params=[200.0, 500.0],
    ids=["200", "500"],
)
def fixture_Delta(request) -> float:
    """Fixture for Delta."""
    return request.param


@pytest.fixture(name="prim")
def fixture_prim() -> Nc.HIPrim:
    """Create primordial power spectrum model."""
    prim = Nc.HIPrimPowerLaw.new()
    prim["n_SA"] = 0.96
    prim["ln10e10ASA"] = np.log(2.1e-9 * 1e10)
    return prim


@pytest.fixture(name="reion")
def fixture_reion() -> Nc.HIReion:
    """Create reionization model."""
    return Nc.HIReionCamb.new()


@pytest.fixture(name="cosmo")
def fixture_cosmo(prim: Nc.HIPrim, reion: Nc.HIReion) -> Nc.HICosmo:
    """Create cosmological model."""
    cosmo = Nc.HICosmoDEXcdm()
    cosmo.omega_x2omega_k()
    cosmo["Omegak"] = 0.0
    cosmo["H0"] = 67.66
    cosmo["Omegab"] = 0.049
    cosmo["Omegac"] = 0.262
    cosmo["w"] = -1.0
    cosmo["Tgamma0"] = 2.7255
    cosmo["ENnu"] = 3.046

    cosmo.add_submodel(prim)
    cosmo.add_submodel(reion)
    return cosmo


@pytest.fixture(name="dist")
def fixture_dist() -> Nc.Distance:
    """Fixture for Distance."""
    return Nc.Distance.new(5.0)


@pytest.fixture(name="psf")
def fixture_psf(cosmo: Nc.HICosmo, prim: Nc.HIPrim) -> Ncm.PowspecFilter:
    """Create power spectrum filter."""
    tf = Nc.TransferFuncEH()
    psml = Nc.PowspecMLTransfer.new(tf)
    psml.require_kmin(1.0e-6)
    psml.require_kmax(1.0e3)
    psf = Ncm.PowspecFilter.new(psml, Ncm.PowspecFilterType.TOPHAT)
    psf.set_best_lnr0()
    psf.prepare(cosmo)
    old_amplitude = np.exp(prim["ln10e10ASA"])
    prim["ln10e10ASA"] = np.log((0.8277 / cosmo.sigma8(psf)) ** 2 * old_amplitude)
    psml.prepare(cosmo)  # Re-prepare the power spectrum
    psf.prepare(cosmo)  # Re-prepare the filter after changing the amplitude
    return psf


@pytest.fixture(name="mfp", params=[Nc.MultiplicityFuncTinker])
def fixture_mass_function(
    request: pytest.FixtureRequest,
    dist: Nc.Distance,
    psf: Ncm.PowspecFilter,
    mass_def: Nc.HaloMassSummaryMassDef,
    Delta: float,
) -> Nc.HaloMassFunction:
    """Fixture for HaloMassFunction."""
    if mass_def == Nc.HaloMassSummaryMassDef.VIRIAL:
        pytest.skip("Tinker mass function does not support VIRIAL mass definition")
    mulf = request.param.new()
    if mass_def == Nc.HaloMassSummaryMassDef.MEAN:
        mulf.set_mdef(Nc.MultiplicityFuncMassDef.MEAN)
    elif mass_def == Nc.HaloMassSummaryMassDef.CRITICAL:
        mulf.set_mdef(Nc.MultiplicityFuncMassDef.CRITICAL)
    mulf.set_Delta(int(Delta))
    mulf.set_linear_interp(True)
    return Nc.HaloMassFunction.new(dist, psf, mulf)


@pytest.fixture(
    name="halo_mass_summary_simple",
    params=[Nc.HaloCMParam, Nc.HaloCMKlypin11],
    ids=["param", "klypin11"],
)
def fixture_halo_mass_summary_simple(
    request: pytest.FixtureRequest, mass_def: Nc.HaloMassSummaryMassDef, Delta: float
) -> Nc.HaloMassSummary:
    """Fixture for simple HaloMassSummary (no mass function needed)."""
    return request.param(mass_def=mass_def, Delta=Delta)


@pytest.fixture(
    name="halo_mass_summary_duffy08",
    params=[Nc.HaloCMDuffy08],
    ids=["duffy08"],
)
def fixture_halo_mass_summary_duffy08(
    request: pytest.FixtureRequest, mass_def: Nc.HaloMassSummaryMassDef
) -> Nc.HaloMassSummary:
    """Fixture for Duffy08 (only Delta=200)."""
    return request.param(mass_def=mass_def, Delta=200.0)


@pytest.fixture(
    name="halo_mass_summary_dutton14",
    params=[
        (Nc.HaloCMDutton14, Nc.HaloMassSummaryMassDef.CRITICAL),
        (Nc.HaloCMDutton14, Nc.HaloMassSummaryMassDef.VIRIAL),
    ],
    ids=["dutton14-critical", "dutton14-virial"],
)
def fixture_halo_mass_summary_dutton14(
    request: pytest.FixtureRequest,
) -> Nc.HaloMassSummary:
    """Fixture for Dutton14 (CRITICAL/200 or VIRIAL only)."""
    cls, mass_def = request.param
    return cls(mass_def=mass_def, Delta=200.0)


@pytest.fixture(
    name="halo_mass_summary_with_mfp",
    params=[Nc.HaloCMPrada12, Nc.HaloCMDiemer15],
    ids=["prada12", "diemer15"],
)
def fixture_halo_mass_summary_with_mfp(request, dist, psf) -> Nc.HaloMassSummary:
    """Fixture for HaloMassSummary that requires mass function (CRITICAL/200 only)."""
    mulf = Nc.MultiplicityFuncTinker.new()
    mulf.set_mdef(Nc.MultiplicityFuncMassDef.CRITICAL)
    mulf.set_Delta(200)
    mulf.set_linear_interp(True)
    mfp = Nc.HaloMassFunction.new(dist, psf, mulf)
    hms: Nc.HaloCMPrada12 | Nc.HaloCMDiemer15 = request.param(
        mass_def=Nc.HaloMassSummaryMassDef.CRITICAL, Delta=200.0, mass_function=mfp
    )

    return hms


def test_halo_mass_summary_simple_basic(
    halo_mass_summary_simple: Nc.HaloMassSummary, cosmo: Nc.HICosmo
) -> None:
    """Test simple HaloMassSummary basic properties."""
    assert halo_mass_summary_simple.concentration(cosmo, 0.0) > 0.0
    assert halo_mass_summary_simple.mass() > 0.0
    assert halo_mass_summary_simple.Delta(cosmo, 0.0) > 0.0
    assert halo_mass_summary_simple.Delta(cosmo, 1.0) > 0.0
    assert halo_mass_summary_simple.rho_bg(cosmo, 0.0) > 0.0
    assert halo_mass_summary_simple.rho_bg(cosmo, 1.0) > 0.0


def test_halo_mass_summary_duffy08_basic(
    halo_mass_summary_duffy08: Nc.HaloMassSummary, cosmo: Nc.HICosmo
) -> None:
    """Test Duffy08 HaloMassSummary basic properties."""
    assert halo_mass_summary_duffy08.concentration(cosmo, 0.0) > 0.0
    assert halo_mass_summary_duffy08.mass() > 0.0
    assert halo_mass_summary_duffy08.Delta(cosmo, 0.0) > 0.0
    assert halo_mass_summary_duffy08.Delta(cosmo, 1.0) > 0.0
    assert halo_mass_summary_duffy08.rho_bg(cosmo, 0.0) > 0.0
    assert halo_mass_summary_duffy08.rho_bg(cosmo, 1.0) > 0.0


def test_halo_mass_summary_dutton14_basic(
    halo_mass_summary_dutton14: Nc.HaloMassSummary, cosmo: Nc.HICosmo
) -> None:
    """Test Dutton14 HaloMassSummary basic properties."""
    assert halo_mass_summary_dutton14.concentration(cosmo, 0.0) > 0.0
    assert halo_mass_summary_dutton14.mass() > 0.0
    assert halo_mass_summary_dutton14.Delta(cosmo, 0.0) > 0.0
    assert halo_mass_summary_dutton14.Delta(cosmo, 1.0) > 0.0
    assert halo_mass_summary_dutton14.rho_bg(cosmo, 0.0) > 0.0
    assert halo_mass_summary_dutton14.rho_bg(cosmo, 1.0) > 0.0


def test_halo_mass_summary_with_mfp_basic(
    halo_mass_summary_with_mfp: Nc.HaloMassSummary, cosmo
):
    """Test HaloMassSummary with mass function basic properties."""
    assert halo_mass_summary_with_mfp.concentration(cosmo, 0.0) > 0.0
    assert halo_mass_summary_with_mfp.mass() > 0.0
    assert halo_mass_summary_with_mfp.Delta(cosmo, 0.0) > 0.0
    assert halo_mass_summary_with_mfp.Delta(cosmo, 1.0) > 0.0
    assert halo_mass_summary_with_mfp.rho_bg(cosmo, 0.0) > 0.0
    assert halo_mass_summary_with_mfp.rho_bg(cosmo, 1.0) > 0.0


def test_halo_mass_summary_simple_mass(
    halo_mass_summary_simple: Nc.HaloMassSummary,
) -> None:
    """Test simple HaloMassSummary mass."""
    log10MDelta = halo_mass_summary_simple["log10MDelta"]
    mass = 10.0**log10MDelta
    assert_allclose(halo_mass_summary_simple.mass(), mass, rtol=1e-5)


def test_halo_mass_summary_duffy08_mass(
    halo_mass_summary_duffy08: Nc.HaloMassSummary,
) -> None:
    """Test Duffy08 HaloMassSummary mass."""
    log10MDelta = halo_mass_summary_duffy08["log10MDelta"]
    mass = 10.0**log10MDelta
    assert_allclose(halo_mass_summary_duffy08.mass(), mass, rtol=1e-5)


def test_halo_mass_summary_dutton14_mass(
    halo_mass_summary_dutton14: Nc.HaloMassSummary,
) -> None:
    """Test Dutton14 HaloMassSummary mass."""
    log10MDelta = halo_mass_summary_dutton14["log10MDelta"]
    mass = 10.0**log10MDelta
    assert_allclose(halo_mass_summary_dutton14.mass(), mass, rtol=1e-5)


def test_halo_mass_summary_param_concentration(cosmo: Nc.HICosmo) -> None:
    """Test HaloCMParam concentration."""
    hms = Nc.HaloCMParam(mass_def=Nc.HaloMassSummaryMassDef.CRITICAL, Delta=200.0)
    cDelta = hms["cDelta"]
    assert_allclose(hms.concentration(cosmo, 0.0), cDelta, rtol=1e-5)


def test_halo_mass_summary_concentration_relations(
    halo_mass_summary_with_mfp: Nc.HaloMassSummary, cosmo
):
    """Test concentration-mass relation."""
    halo_mass_summary_with_mfp["log10MDelta"] = 12.5
    cDelta1 = halo_mass_summary_with_mfp.concentration(cosmo, 0.0)
    halo_mass_summary_with_mfp["log10MDelta"] = 13.5
    cDelta2 = halo_mass_summary_with_mfp.concentration(cosmo, 0.0)
    assert cDelta1 > cDelta2


def test_halo_cm_klypin11_compare_colossus(cosmo: Nc.HICosmo) -> None:
    """Test Klypin11 against Colossus reference data."""

    # Pre-generated from Colossus (shape: len(Z_ARRAY) x len(LOG10M_ARRAY))
    c_colossus = np.array(
        [
            12.48783463,
            11.51609042,
            10.61996274,
            9.79356748,
            9.03147839,
            8.32869147,
            7.6805921,
            7.08292476,
            6.53176507,
            6.02349401,
        ]
    )

    # Setup
    hms = Nc.HaloCMKlypin11(mass_def=Nc.HaloMassSummaryMassDef.VIRIAL, Delta=200.0)

    # NumCosmo computation
    log10m_array = np.linspace(
        np.log10(3.0) + 10.0, np.log10(5) + 14.0, 10
    )  # mass in Msun, colossus mass in Msun/h
    c_nc = np.zeros(len(log10m_array))
    for i, log10M in enumerate(log10m_array):
        hms["log10MDelta"] = log10M - np.log10(cosmo.h())
        c_nc[i] = hms.concentration(cosmo, 0.0)

    assert_allclose(c_nc, c_colossus, rtol=1e-5)


R_ARRAY = 10 ** np.linspace(-2.0, 2.0, 100)
SIGMA_ARRAY = [
    5.20970156,
    5.12178608,
    5.03424238,
    4.9470775,
    4.86029942,
    4.77391701,
    4.68794003,
    4.60237903,
    4.51724537,
    4.43255115,
    4.3483092,
    4.264533,
    4.18123667,
    4.09843492,
    4.01614297,
    3.93437124,
    3.85311995,
    3.77238973,
    3.69218305,
    3.61250413,
    3.53335887,
    3.4547548,
    3.37670097,
    3.2992079,
    3.22228747,
    3.14595287,
    3.07021851,
    2.99509992,
    2.92061366,
    2.84677727,
    2.77360914,
    2.70112843,
    2.62935499,
    2.55830928,
    2.48800918,
    2.4184611,
    2.34967033,
    2.28164444,
    2.21439309,
    2.14792784,
    2.08226206,
    2.0174107,
    1.95339017,
    1.89021813,
    1.82791335,
    1.76649556,
    1.70598523,
    1.64640342,
    1.58777163,
    1.53011159,
    1.47344472,
    1.41778312,
    1.36313117,
    1.30949534,
    1.25688427,
    1.20530848,
    1.15478008,
    1.10531247,
    1.05692005,
    1.0096179,
    0.96342154,
    0.91834662,
    0.87440864,
    0.83162273,
    0.79000311,
    0.74955364,
    0.71026628,
    0.67213536,
    0.63515813,
    0.59933419,
    0.56466495,
    0.53115305,
    0.49880193,
    0.46761523,
    0.43759643,
    0.40874833,
    0.38107269,
    0.35457191,
    0.32924973,
    0.30510658,
    0.28213908,
    0.26034002,
    0.23969838,
    0.22019944,
    0.2018249,
    0.18455303,
    0.16835887,
    0.15321526,
    0.13909655,
    0.12597412,
    0.11381496,
    0.10258239,
    0.09223676,
    0.08273617,
    0.07403709,
    0.06609493,
    0.0588635,
    0.0522927,
    0.0463341,
    0.04094235,
]


def test_sigma_agreement(psf: Ncm.PowspecFilter, cosmo: Nc.HICosmo) -> None:
    """Test that NumCosmo's sigma calculations are consistent and monotonic."""

    # Generate reference sigma values with current cosmology at z=0
    sigma_ref = np.zeros(len(R_ARRAY))
    for i, r in enumerate(R_ARRAY):
        sigma_ref[i] = psf.eval_sigma_lnr(0.0, np.log(r))

    # Test 1: Sigma should be monotonically decreasing with increasing R
    assert np.all(np.diff(sigma_ref) < 0), "Sigma should decrease with increasing R"

    # Test 2: Recompute and verify consistency
    sigma_test = np.zeros(len(R_ARRAY))
    for i, r in enumerate(R_ARRAY):
        sigma_test[i] = psf.eval_sigma_lnr(0.0, np.log(r))

    assert_allclose(sigma_test, sigma_ref, rtol=1e-10)

    # Test 3: Sigma8 should match the calibrated value (0.8)
    # NumCosmo computes sigma8 at R = 8 Mpc/h = 8.0 / h Mpc
    R_8Mpch = 8.0 / cosmo.h()
    sigma_8 = psf.eval_sigma_lnr(0.0, np.log(R_8Mpch))
    assert_allclose(sigma_8, 0.8277, rtol=1.0e-5)


def test_sigma_compare_colossus(psf: Ncm.PowspecFilter, cosmo: Nc.HICosmo) -> None:
    """Test that sigma8 matches between psf.eval_sigma_lnr and colossus."""
    sigma_a = [psf.eval_sigma_lnr(1.0, np.log(r / cosmo.h())) for r in R_ARRAY]

    assert_allclose(sigma_a, SIGMA_ARRAY, rtol=1.0e-2)


def test_halo_cm_prada12_compare_colossus(
    dist: Nc.Distance, psf: Ncm.PowspecFilter, cosmo: Nc.HICosmo
) -> None:
    """Test Prada12 against Colossus reference data."""

    # Pre-generated from Colossus (shape: len(Z_ARRAY) x len(log10m_array))
    c_colossus = np.array(
        [
            [11.79827496, 9.72828123, 7.92720837, 6.44725001, 5.3786824],
            [9.4485788, 7.93773446, 6.63831193, 5.60850411, 4.98398234],
            [7.62437844, 6.55155061, 5.65467982, 5.01202859, 4.85955869],
            [6.1623863, 5.44850512, 4.89954166, 4.64316421, 5.16257546],
            [5.19911526, 4.73116057, 4.4444199, 4.55281368, 6.04165491],
        ]
    )

    # Setup mass function for CRITICAL/200
    mulf = Nc.MultiplicityFuncTinker.new()
    mulf.set_mdef(Nc.MultiplicityFuncMassDef.CRITICAL)
    mulf.set_Delta(200)
    mulf.set_linear_interp(True)
    mfp = Nc.HaloMassFunction.new(dist, psf, mulf)

    # Setup
    hms = Nc.HaloCMPrada12(
        mass_def=Nc.HaloMassSummaryMassDef.CRITICAL, Delta=200.0, mass_function=mfp
    )

    # NumCosmo computation
    log10M_array = np.log10(10.0 ** np.arange(10.0, 15.0, 1.0) / cosmo.h())
    z_array = np.linspace(0.0, 2.0, 5)
    c_nc = np.zeros((len(z_array), len(log10M_array)))
    for i, z in enumerate(z_array):
        for j, log10M in enumerate(log10M_array):
            hms["log10MDelta"] = log10M
            c_nc[i, j] = hms.concentration(cosmo, z)

    assert_allclose(c_nc, c_colossus, rtol=3e-4)

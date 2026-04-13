"""Test fixtures for numcosmo_py."""

import pytest
import numpy as np
from numcosmo_py import Nc, Ncm
from numcosmo_py.cosmology import Cosmology

Ncm.cfg_init()


def pytest_addoption(parser):
    """Add custom command-line options for conditional test execution."""
    parser.addoption(
        "--run-mpi",
        action="store_true",
        default=False,
        help="Run tests marked with mpi",
    )
    parser.addoption(
        "--run-powspec",
        action="store_true",
        default=False,
        help="Run tests marked with powspec",
    )
    parser.addoption(
        "--run-xcor",
        action="store_true",
        default=False,
        help="Run tests marked with xcor",
    )


def pytest_configure(config):
    """Configure pytest marker behavior."""
    config.addinivalue_line("markers", "mpi: MPI-parallel tests (run with --run-mpi)")
    config.addinivalue_line(
        "markers", "powspec: Power spectrum tests (run with --run-powspec)"
    )
    config.addinivalue_line(
        "markers", "xcor: Cross-correlation tests (run with --run-xcor)"
    )


def pytest_collection_modifyitems(config, items):
    """Skip tests based on markers and command-line options."""
    run_mpi = config.getoption("--run-mpi")
    run_powspec = config.getoption("--run-powspec")
    run_xcor = config.getoption("--run-xcor")

    skip_mpi = pytest.mark.skip(reason="Need --run-mpi option to run")
    skip_powspec = pytest.mark.skip(reason="Need --run-powspec option to run")
    skip_xcor = pytest.mark.skip(reason="Need --run-xcor option to run")

    for item in items:
        if "mpi" in item.keywords and not run_mpi:
            item.add_marker(skip_mpi)
        if "powspec" in item.keywords and not run_powspec:
            item.add_marker(skip_powspec)
        if "xcor" in item.keywords and not run_xcor:
            item.add_marker(skip_xcor)


@pytest.fixture(name="prim")
def fixture_prim() -> Nc.HIPrim:
    """Create primordial power spectrum model."""
    prim = Nc.HIPrimPowerLaw.new()
    prim.props.n_SA = 0.967
    return prim


@pytest.fixture(name="reion")
def fixture_reion() -> Nc.HIReion:
    """Create reionization model."""
    return Nc.HIReionCamb.new()


@pytest.fixture(name="cosmo")
def fixture_cosmo(prim: Nc.HIPrim, reion: Nc.HIReion) -> Nc.HICosmo:
    """Create cosmological model with dark energy.

    Configures a flat LCDM cosmology with w = -1.
    """
    cosmo = Nc.HICosmoDEXcdm()
    cosmo.omega_x2omega_k()
    cosmo["Omegak"] = 0.0
    cosmo["H0"] = 71
    cosmo["Omegab"] = 0.0406
    cosmo["Omegac"] = 0.22
    cosmo["w"] = -1.0

    cosmo.add_submodel(prim)
    cosmo.add_submodel(reion)

    return cosmo


@pytest.fixture(name="dist")
def fixture_distribution() -> Nc.Distance:
    """Fixture for Distance."""
    return Nc.Distance.new(2000.0)


@pytest.fixture(name="psf")
def fixture_psf(cosmo: Nc.HICosmo, prim: Nc.HIPrim) -> Ncm.PowspecFilter:
    """Create power spectrum filter.

    Configures a tophat filter and normalizes to sigma8 = 0.8.
    """
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


@pytest.fixture(name="nc_cosmo_default", scope="module")
def fixture_nc_cosmo_default() -> Cosmology:
    """Fixture for default NumCosmo Cosmology."""
    return Cosmology.default()


@pytest.fixture(name="nc_cosmo_alt", scope="module")
def fixture_cosmology_alt() -> Cosmology:
    """Create a simple cosmology alternative for testing."""
    cosmology = Cosmology.default()
    cosmology.cosmo["H0"] = 75.0
    cosmology.cosmo["Omegab"] = 0.03
    cosmology.cosmo["Omegac"] = 0.22

    cosmology.prepare()

    return cosmology

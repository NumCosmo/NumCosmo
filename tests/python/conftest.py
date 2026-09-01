"""Test fixtures for numcosmo_py."""

import os
import sys
import warnings
import faulthandler
import pytest
import numpy as np
from gi import PyGIDeprecationWarning  # type: ignore
from gi.repository import GLib

warnings.filterwarnings(
    "ignore", message=".*unix_signal_add_full.*", category=PyGIDeprecationWarning
)
# flake8: noqa: E402
# pylint: disable=wrong-import-position
from numcosmo_py import Nc, Ncm
from numcosmo_py.cosmology import Cosmology

Ncm.cfg_init()

# A single slow-test dump only shows one snapshot in time, so there is no way to
# tell a genuine hang (identical stack next time) from a test that is merely slow
# but still making progress (advancing line/loop state next time) -- exactly the
# ambiguity that made a real CI timeout hard to diagnose. pytest's builtin
# faulthandler_timeout ini option can't help here: _pytest/faulthandler.py calls
# faulthandler.dump_traceback_later() without repeat=True, so it only ever fires
# once per test. This hook re-implements the same wrap (dup stderr's fd once so
# the dump survives pytest's output capturing, exactly as the builtin plugin
# does) but with repeat=True, so a slow test gets a fresh stack dump every
# _FAULTHANDLER_INTERVAL seconds for as long as it keeps running.
_FAULTHANDLER_INTERVAL = 300.0


def _stderr_fileno() -> int:
    try:
        fileno = sys.stderr.fileno()
        if fileno == -1:
            raise AttributeError()
        return fileno
    except (AttributeError, ValueError):
        # pytest-xdist monkeypatches sys.stderr with a non-file object.
        assert sys.__stderr__ is not None
        return sys.__stderr__.fileno()


# The library reports a fatal condition through its own GLib log handler
# (ncm_cfg.c installs one for G_LOG_LEVEL_ERROR), which writes to fd 2 and then
# aborts. pytest captures at fd level and typer's CliRunner adds a second layer,
# so that write lands in a buffer nobody ever prints: the process is gone before
# the capture is reported, and all CI shows is "Fatal Python error: Aborted"
# with no reason at all. Routing fatal messages to the duplicated fd -- the same
# one the faulthandler dumps use, and for the same reason -- puts them where
# they survive.
#
# ncm_cfg_set_error_log_handler() is the library's own hook for this;
# GLib.log_set_writer_func() does not see these messages, because a legacy
# g_log_set_handler() for the domain takes precedence over the structured
# writer.
def _install_fatal_log_mirror(fd: int) -> None:
    """Send the library's fatal messages to @fd, which capturing does not touch."""

    def logger(message: str) -> None:
        try:
            os.write(fd, message.encode())
        except OSError:
            pass

    Ncm.cfg_set_error_log_handler(logger)
    # Keep it alive: the C side stores the callback without a reference.
    _FATAL_LOGGERS.append(logger)


_FATAL_LOGGERS: list = []


def pytest_configure(config):
    """Dup stderr's fd once so later dumps survive per-test output capturing."""
    fd = os.dup(_stderr_fileno())
    config.stash["faulthandler_dup_fd"] = fd
    _install_fatal_log_mirror(fd)


def pytest_unconfigure(config):
    fd = config.stash.get("faulthandler_dup_fd", None)
    if fd is not None:
        os.close(fd)


@pytest.hookimpl(wrapper=True, trylast=True)
def pytest_runtest_protocol(item):
    """Re-arm a repeating faulthandler dump around every test (see above)."""
    fd = item.config.stash.get("faulthandler_dup_fd", None)
    if fd is not None:
        faulthandler.dump_traceback_later(
            _FAULTHANDLER_INTERVAL, repeat=True, file=fd, exit=False
        )
    try:
        return (yield)
    finally:
        if fd is not None:
            faulthandler.cancel_dump_traceback_later()


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
    parser.addoption(
        "--run-sphere-map",
        action="store_true",
        default=False,
        help="Run tests marked with sphere_map",
    )
    parser.addoption(
        "--run-app",
        action="store_true",
        default=False,
        help="Run tests marked with app",
    )
    parser.addoption(
        "--run-planck-data",
        action="store_true",
        default=False,
        help="Run tests marked with planck_data",
    )


#: Capability marker -> the option that opts into it. A test carrying the marker is
#: skipped unless its option is given.
_CAPABILITY_GATES = {
    "mpi": "--run-mpi",
    "powspec": "--run-powspec",
    "xcor": "--run-xcor",
    "sphere_map": "--run-sphere-map",
    "app": "--run-app",
    "planck_data": "--run-planck-data",
}


def pytest_collection_modifyitems(config, items):
    """Skip tests based on markers and command-line options.

    Gates on the marker, via ``get_closest_marker``, not on ``item.keywords``. A test's
    keywords include the names of the directories it lives under, and three of these
    gates -- app, xcor, powspec -- are also directory names, so keyword membership
    skipped every test under those directories whether or not it carried the marker.
    That is invisible while every file in them happens to be marked, and silently
    swallows the first one that is not.
    """
    for marker, option in _CAPABILITY_GATES.items():
        if config.getoption(option):
            continue

        skip = pytest.mark.skip(reason=f"Need {option} option to run")

        for item in items:
            if item.get_closest_marker(marker) is not None:
                item.add_marker(skip)


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
    cosmo = Nc.HICosmoDEXcdm(prim=prim, reion=reion)
    cosmo.omega_x2omega_k()
    cosmo["Omegak"] = 0.0
    cosmo["H0"] = 71
    cosmo["Omegab"] = 0.0406
    cosmo["Omegac"] = 0.22
    cosmo["w"] = -1.0

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

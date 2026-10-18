"""Test fixtures for numcosmo_py."""

import os
import sys
import warnings
import faulthandler
from pathlib import Path

import pytest
import numpy as np
from gi import PyGIDeprecationWarning  # type: ignore

# Installed before the GLib import, which is what emits it: pygobject's override
# machinery touches the deprecated attribute while building the namespace, so a filter
# added afterwards never sees the warning.
warnings.filterwarnings(
    "ignore", message=".*unix_signal_add_full.*", category=PyGIDeprecationWarning
)

from gi.repository import GLib

# flake8: noqa: E402
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


# The five header lines and the first six rows of the EDD export of All CF4 SNIa Samples
# (2026-10-17), for the NcDistMuObs importer and the generate bianchi command.
EDD_CF4_SNIA_EXCERPT = """HYPERLEDA|"All CF4 SNIa Samples"|"All CF4 SNIa Samples"|"All CF4 SNIa Samples"|"All CF4 SNIa Samples"|"All CF4 SNIa Samples"|"All CF4 SNIa Samples"|"All CF4 SNIa Samples"|"All CF4 SNIa Samples"|"All CF4 SNIa Samples"|"All CF4 SNIa Samples"|"All CF4 SNIa Samples"|"All CF4 SNIa Samples"|"All CF4 SNIa Samples"|"All CF4 SNIa Samples"|"All CF4 SNIa Samples"|"All CF4 SNIa Samples"|"All CF4 SNIa Samples"|"All CF4 SNIa Samples"|"All CF4 SNIa Samples"|"All CF4 SNIa Samples"|"All CF4 SNIa Samples"|"All CF4 SNIa Samples"|"All CF4 SNIa Samples"|"All CF4 SNIa Samples"|"All CF4 SNIa Samples"|"All CF4 SNIa Samples"|"All CF4 SNIa Samples"|"All CF4 SNIa Samples"|"All CF4 SNIa Samples"|"All CF4 SNIa Samples"|"All CF4 SNIa Samples"|"All CF4 SNIa Samples"|"All CF4 SNIa Samples"|"All CF4 SNIa Samples"|"All CF4 SNIa Samples"|"All CF4 SNIa Samples"|"All CF4 SNIa Samples"|"All CF4 SNIa Samples"|"All CF4 SNIa Samples"|"All CF4 SNIa Samples"|"All CF4 SNIa Samples"|"All CF4 SNIa Samples"|"All CF4 SNIa Samples"|"All CF4 SNIa Samples"
pgc|SNIa|IDT17g|1PGC|Vcmb|Vmod|DM_zp|DM_av|eDM|DM_pantheonplus|eDM_pantheonplus|DM_scolnic|eDM_scolnic|DM_csp1|eDM_csp1|DM_ganesh|eDM_ganesh|DM_rest|eDM_rest|DM_amanullah|eDM_amanullah|DM_jha|eDM_jha|DM_prieto|eDM_prieto|DM_hicken|eDM_hicken|DM_folatelli|eDM_folatelli|DM_walker|eDM_walker|DM_stahl|eDM_stahl|DM_sh0es|eDM_sh0es|DM_twins|eDM_twins|DM_avelino|eDM_avelino|RA|DE|glon|glat|sgl|sgb
%7s|%17s|%5d|%7d|%5d|%5d|%6.3f|%5.2f|%5.2f|%5.2f|%5.2f|%5.2f|%5.2f|%5.2f|%5.2f|%5.2f|%5.2f|%5.2f|%5.2f|%5.2f|%5.2f|%5.2f|%5.2f|%5.2f|%5.2f|%5.2f|%5.2f|%5.2f|%5.2f|%5.2f|%5.2f|%5.2f|%5.2f|%5.2f|%5.2f|%5.2f|%5.2f|%5.2f|%5.2f|%9.4f|%9.4f|%9.4f|%9.4f|%9.4f|%9.4f
---||||km/s|km/s|mag|mag|mag|mag|mag|mag|mag|mag|mag|mag|mag|mag|mag|mag|mag|mag|mag|mag|mag|mag|mag|mag|mag|mag|mag|mag|mag|mag|mag|mag|mag|mag|mag|deg|deg|deg|deg|deg|deg
"The Catalogue of Principal Galaxies (PGC) Number"|"SNIa name; multiple SN ID separated by /"|"Group ID in Tempel 2017 SDSS group catalog"|"PGC identification of brightest member of group"|"Radial velocity (cz) with respect to the CMB radiation"|"Velocity in CMB frame; cosmological adjustment described in CF2 data paper"|"Weighted average distance modulus shifted to absolute scale set by calibrators"|"Weighted average distance moduli (w=1/sigma^2) [excluding Walker+15]"|"Minimum error of all measured moduli"|"SNIa distance moduli in Scolnic Pantheonplus20"|"rms dispersion of the Pantheonplus20 distance measures (adding 0.15 mag in quadrature)"|"distance moduli in Scolnic SNIa 3 Tables on the Pantheonplus20 scales"|"rms dispersion of the Scolnic 3 Tables distance measures (adding 0.15 mag in quadrature)"|"distance moduli in CSP1 Burns18 (Burns+2018,ApJ,869:56) on the Pantheonplus20 scales"|"rms dispersion of the CSP1 Burns18 distance measures (adding 0.15 mag in quadrature)"|"distance moduli in LOSS Ganeshalingam13 (Ganeshalingam+2013,MNRAS,433,2240) on the Pantheonplus20 scales"|"rms dispersion of the LOSS Ganeshalingam13 distance measures (adding 0.15 mag in quadrature)"|"Rest+ (2014,ApJ,795:44) distance moduli in All CF3 SNIa on the Pantheonplus20 scales"|"rms dispersion of the Rest+ distance measures (adding 0.15 mag in quadrature)"|"Union2 (Amanullah et. 2010,ApJ,716,712) distance moduli in All CF3 SNIa on the Pantheonplus20 scales"|"rms dispersion of the Union2 distance measures (adding 0.15 mag in quadrature)"|"Jha et al. (2007,ApJ,659,122) distance moduli in All CF3 SNIa on the Pantheonplus20 scales"|"rms dispersion of the Jha"|"Prieto et al. (2006,ApJ,647,501) distance moduli in All CF3 SNIa on the Pantheonplus20 scales"|"rms dispersion of the Prieto"|"Hicken et al. (2009,ApJ,700,109) distance moduli in All CF3 SNIa on the Pantheonplus20 scales"|"rms dispersion of the Hicken"|"Folatelli et al. (2010,AJ,139,120) distance moduli in All CF3 SNIa on the Pantheonplus20 scales"|"rms dispersion of the Folatelli"|"Walker et al. (2015,ApJS,219:13) distance moduli on CF2 scale"|"rms dispersion of the Walker"|"Stahl et al. (2021,MNRAS,505,2300) distance moduli on Stahl"|"rms dispersion of the Stahl"|"Pantheon+SH0ES (2021,arXiv:2112.04510) distance moduli on SH0ES scale"|"Error estimates for Pantheonplus+SH0ES measures of corrected peak B magnitude"|"Boone et al. (2021,ApJ, 912, 21) distance moduli on scale assuming H0=75"|"Mean rms dispersion of twin sample distance measures"|"Avelino et al. (2019,ApJ, 887, 106); best of YJHK template or YJHK IRmax distance moduli on scale assuming H0=73.24"|"Approximated rms dispersion of Avelino distance measures"|"Right Ascension of host, J2000"|"Declination of host, J2000"|"Galactic longitude of host"|"Galactic latitude of host"|"Supergalactic longitude of host"|"Supergalactic latitude of host"
250|2006sr|0|250|6914|7042|34.912|35.05|0.13|35.07|0.25|35.03|0.20|||35.05|0.23|35.05|0.22|35.04|0.23|||||35.04|0.15|||||35.10|0.17|35.05|0.13|||||0.8958|23.2008|108.8693|-38.3606|316.3211|17.6720
282|2008hj|0|282|10921|11235|35.902|36.04|0.11|35.98|0.24|36.00|0.18|36.15|0.16|||36.07|0.26|||||||||||||36.12|0.13|35.96|0.11|||35.97|0.12|1.0064|-11.1742|85.7160|-70.5385|282.0792|8.9703
415|2003it|0|415|7165|7301|34.992|35.13|0.14|35.14|0.25|35.12|0.21|||||35.09|0.23|35.14|0.24|||||35.20|0.15|||||35.10|0.18|35.13|0.14|||||1.4517|27.4491|110.6597|-34.3329|320.8087|17.8681
465|2000dg|0|616|11092|11415|35.912|36.05|0.20|||||||36.00|0.24|||||||||||||||35.99|0.22|36.13|0.20|||||1.5588|8.8878|104.4513|-52.3545|301.9054|14.0485
515|2002hw|0|515|4796|4857|34.142|34.28|0.15|||||||34.23|0.25|||34.40|0.26|||||34.38|0.15|||||34.10|0.19|||||||1.7062|8.6285|104.5597|-52.6445|301.6830|13.8439
963|1998dk|0|963|3589|3623|33.312|33.45|0.16|33.37|0.28|33.31|0.24|||||||33.34|0.32|33.58|0.27|33.46|0.16|||||||33.63|0.27|||||||3.6328|-0.7377|102.8519|-62.1619|292.9641|9.5361
"""


@pytest.fixture(name="edd_cf4_excerpt_file")
def fixture_edd_cf4_excerpt_file(tmp_path: Path) -> Path:
    """Write the EDD excerpt of All CF4 SNIa Samples to a file."""
    path = tmp_path / "edd_cf4_snia.txt"
    path.write_text(EDD_CF4_SNIA_EXCERPT)
    return path

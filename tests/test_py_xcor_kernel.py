#
# test_py_xcor_kernel.py
#
# Tue Feb 18 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_py_xcor_kernel.py
# Copyright (C) 2026 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

"""Unit tests for NcXcorKernel.

Tests all kernel types using parametrized kernel_case fixture.
Kernels are collected with all relevant configurations (e.g., galaxy with/without magbias).
"""

from typing import Callable
from functools import cache
import pytest
from numpy.testing import assert_allclose
import numpy as np

from numcosmo_py import Ncm, Nc
from numcosmo_py.cosmology import Cosmology

Ncm.cfg_init()


@cache
def get_cosmology() -> Cosmology:
    """Create a default cosmology for testing."""
    return Cosmology.default()


@pytest.fixture(name="cosmology", scope="module")
def fixture_cosmology() -> Cosmology:
    """Create a default cosmology for testing."""
    return get_cosmology()


@pytest.fixture(name="integrator", scope="module")
def fixture_integrator() -> Ncm.SBesselIntegrator:
    """Create a spherical Bessel integrator."""
    return Ncm.SBesselIntegratorLevin.new(0, 2000)


def _create_galaxy_kernel(
    cosmology: Cosmology, integrator: Ncm.SBesselIntegrator, domagbias: bool = False
) -> Nc.XcorKernelGal:
    """Helper to create a galaxy kernel."""
    dn_dz = Ncm.SplineCubicNotaknot()
    z_array = np.linspace(0.1, 2.0, 50)
    dndz_array = np.exp(-((z_array - 0.8) ** 2) / (2.0 * 0.3**2))
    dndz_array /= np.trapezoid(dndz_array, z_array)

    xv = Ncm.Vector.new_array(z_array.tolist())
    yv = Ncm.Vector.new_array(dndz_array.tolist())
    dn_dz.set(xv, yv, True)

    gal = Nc.XcorKernelGal(
        dndz=dn_dz,
        dist=cosmology.dist,
        powspec=cosmology.ps_ml,
        integrator=integrator,
        domagbias=domagbias,
    )
    gal["bparam_0"] = 1.5
    return gal


def _create_cmb_lensing_kernel(
    cosmology: Cosmology, integrator: Ncm.SBesselIntegrator
) -> Nc.XcorKernelCMBLensing:
    """Helper to create a CMB lensing kernel."""
    lmax = 1000
    ell_vec = Ncm.Vector.new_array(np.arange(lmax + 1).tolist())
    return Nc.XcorKernelCMBLensing(
        dist=cosmology.dist,
        powspec=cosmology.ps_ml,
        recomb=cosmology.recomb,
        Nl=ell_vec,
        integrator=integrator,
    )


def _create_cmb_isw_kernel(
    cosmology: Cosmology, integrator: Ncm.SBesselIntegrator
) -> Nc.XcorKernelCMBISW:
    """Helper to create a CMB ISW kernel."""
    lmax = 1000
    ell_vec = Ncm.Vector.new_array(np.arange(lmax + 1).tolist())
    return Nc.XcorKernelCMBISW(
        dist=cosmology.dist,
        powspec=cosmology.ps_ml,
        recomb=cosmology.recomb,
        Nl=ell_vec,
        integrator=integrator,
    )


def _create_tsz_kernel(
    cosmology: Cosmology, integrator: Ncm.SBesselIntegrator
) -> Nc.XcorKerneltSZ:
    """Helper to create a tSZ kernel."""
    return Nc.XcorKerneltSZ(
        dist=cosmology.dist,
        powspec=cosmology.ps_ml,
        zmax=6.0,
        integrator=integrator,
    )


def _create_weak_lensing_kernel(
    cosmology: Cosmology, integrator: Ncm.SBesselIntegrator
) -> Nc.XcorKernelWeakLensing:
    """Helper to create a weak lensing kernel."""
    dn_dz = Ncm.SplineCubicNotaknot()
    z_array = np.linspace(0.1, 3.0, 50)
    dndz_array = (z_array / 1.0) ** 2 * np.exp(-((z_array / 1.0) ** 1.5))
    dndz_array /= np.trapezoid(dndz_array, z_array)

    xv = Ncm.Vector.new_array(z_array.tolist())
    yv = Ncm.Vector.new_array(dndz_array.tolist())
    dn_dz.set(xv, yv, True)

    return Nc.XcorKernelWeakLensing(
        dndz=dn_dz,
        dist=cosmology.dist,
        powspec=cosmology.ps_ml,
        integrator=integrator,
    )


def _collect_all_kernels(
    cosmology: Cosmology, integrator: Ncm.SBesselIntegrator
) -> list[tuple[str, Nc.XcorKernel]]:
    """Collect all kernel types with all configurations.

    Returns a list of tuples: (kernel_id, kernel)
    """
    cosmo = cosmology.cosmo
    all_kernels: list[tuple[str, Nc.XcorKernel]] = []

    # Define all kernel types with configurations
    kernels_to_create: list[tuple[str, Callable[[], Nc.XcorKernel]]] = [
        ("Gal", lambda: _create_galaxy_kernel(cosmology, integrator, domagbias=False)),
        (
            "Gal:MagBias",
            lambda: _create_galaxy_kernel(cosmology, integrator, domagbias=True),
        ),
        ("CMBLensing", lambda: _create_cmb_lensing_kernel(cosmology, integrator)),
        ("CMBISW", lambda: _create_cmb_isw_kernel(cosmology, integrator)),
        ("tSZ", lambda: _create_tsz_kernel(cosmology, integrator)),
        ("WeakLensing", lambda: _create_weak_lensing_kernel(cosmology, integrator)),
    ]

    # Create each kernel
    for kernel_id, kernel_factory in kernels_to_create:
        kernel = kernel_factory()
        kernel.prepare(cosmo)
        all_kernels.append((kernel_id, kernel))

    return all_kernels


@cache
def get_kernel_cases():
    """Get all kernel cases for parametrization."""
    cosmology = get_cosmology()
    integrator = Ncm.SBesselIntegratorLevin.new(0, 2000)
    return _collect_all_kernels(cosmology, integrator)


def pytest_generate_tests(metafunc: pytest.Metafunc) -> None:
    """Dynamically generate tests for all kernel types."""
    if "kernel_case" in metafunc.fixturenames:
        cases = get_kernel_cases()
        ids = [kernel_id for kernel_id, _ in cases]
        metafunc.parametrize("kernel_case", cases, ids=ids)


def test_limber_vs_limber_z(
    kernel_case: tuple[str, Nc.XcorKernel], cosmology: Cosmology
) -> None:
    """Test that limber and limber_z give consistent results."""
    _, kernel = kernel_case
    cosmo = cosmology.cosmo
    dist = cosmology.dist
    ps_ml = cosmology.ps_ml
    kernel.set_l_limber(0)  # Treat all ells as limber
    kernel.prepare(cosmo)

    # Test ell values
    ell_array = np.array([10, 50, 100, 500, 1000])

    # Get vectorized evaluation functions
    for ell in ell_array:
        limber_func = kernel.get_eval(cosmo, ell)
        nu = ell + 0.5

        # Evaluate both methods
        kmin, kmax = limber_func.get_range()

        k_array = np.geomspace(kmin, kmax, 100)
        for k in k_array:
            xi = nu / k
            k_Mpc = k / cosmo.RH_Mpc()
            z = dist.inv_comoving(cosmo, xi)
            pk = ps_ml.eval(cosmo, z, k_Mpc)
            result = limber_func.eval_array(k)[0]
            result_z = (
                kernel.eval_limber_z_full(cosmo, z, dist, ell)
                * np.sqrt(pk)
                * np.sqrt(np.pi / 2.0 / nu)
                / k
            )
            assert_allclose(
                result,
                result_z,
                rtol=1e-13,
                atol=0.0,
                err_msg=(
                    f"Limber and limber_z differ at ell={ell}, k={k:.3e}, z={z:.3f}"
                ),
            )

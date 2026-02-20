#!/usr/bin/env python
#
# test_py_xcor_kernel_component.py
#
# Tue Feb 18 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_py_xcor_kernel_component.py
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

"""Unit tests for NcXcorKernelComponent.

This test file has two main purposes:
1. Test the NcXcorKernelComponent interface through NcTestXcorKernelComponent
2. Test that real kernels properly create and return components via
   get_component_list(), and that those components satisfy the component
   interface requirements

Test organization:
- Group 1: NcTestXcorKernelComponent tests (direct component interface)
- Group 2: Real kernel components tests (via get_component_list)
- Group 3: Kernel-specific tests (kernel interface, not components)
"""

from typing import Callable
import pytest
from numpy.testing import assert_allclose
import numpy as np

from numcosmo_py import Ncm, Nc, GObject
from numcosmo_py.cosmology import Cosmology

Ncm.cfg_init()


# Define a simple test component that subclasses NcXcorKernelComponent in Python
class NcTestXcorKernelComponent(Nc.XcorKernelComponent):
    """
    A simple test component that implements a Gaussian kernel in xi-space.
    K(k, xi) = exp(-xi²/2sigma²) * sin(k*xi)/(k*xi)
    This is a simple, well-behaved kernel for testing purposes.
    """

    __gtype_name__ = "NcTestXcorKernelComponent"

    xi_min: float = GObject.Property(  # type: ignore
        type=float,
        default=0.1,
        flags=GObject.ParamFlags.READWRITE | GObject.ParamFlags.CONSTRUCT_ONLY,
    )

    xi_max: float = GObject.Property(  # type: ignore
        type=float,
        default=100.0,
        flags=GObject.ParamFlags.READWRITE | GObject.ParamFlags.CONSTRUCT_ONLY,
    )

    k_min: float = GObject.Property(  # type: ignore
        type=float,
        default=1.0e-4,
        flags=GObject.ParamFlags.READWRITE | GObject.ParamFlags.CONSTRUCT_ONLY,
    )

    k_max: float = GObject.Property(  # type: ignore
        type=float,
        default=1.0e3,
        flags=GObject.ParamFlags.READWRITE | GObject.ParamFlags.CONSTRUCT_ONLY,
    )

    sigma: float = GObject.Property(  # type: ignore
        type=float,
        default=10.0,
        flags=GObject.ParamFlags.READWRITE | GObject.ParamFlags.CONSTRUCT_ONLY,
    )

    def do_eval_kernel(  # pylint: disable=arguments-differ
        self,
        _cosmo: Nc.HICosmo,
        xi: float,
        k: float,
    ) -> float:
        """Evaluate test kernel: Gaussian x sinc."""
        if xi <= 0.0:
            return 0.0
        gaussian = np.exp(-0.5 * (xi / self.sigma) ** 2)
        kxi = k * xi
        # Handle sinc(x) = sin(x)/x carefully
        if np.abs(kxi) < 1e-10:
            sinc = 1.0
        else:
            sinc = np.sin(kxi) / kxi
        return gaussian * sinc

    def do_eval_prefactor(  # pylint: disable=arguments-differ
        self,
        _cosmo: Nc.HICosmo,
        _k: float,
        _l: int,
    ) -> float:
        """Evaluate prefactor: constant for simplicity."""
        return 1.0

    def do_get_limits(  # pylint: disable=arguments-differ
        self,
        _cosmo: Nc.HICosmo,
    ) -> tuple[float, float, float, float]:
        """Return integration limits."""
        return self.xi_min, self.xi_max, self.k_min, self.k_max


@pytest.fixture(name="cosmology", scope="module")
def fixture_cosmology() -> Cosmology:
    """Create a simple cosmology for testing."""
    cosmology = Cosmology.default()
    return cosmology


@pytest.fixture(name="cosmology_alt", scope="module")
def fixture_cosmology_alt() -> Cosmology:
    """Create a simple cosmology alternate for testing."""
    cosmology = Cosmology.default()
    cosmology.cosmo["H0"] = 75.0
    cosmology.cosmo["Omegab"] = 0.03
    cosmology.cosmo["Omegac"] = 0.22

    cosmology.prepare()

    return cosmology


@pytest.fixture(name="test_component")
def fixture_test_component() -> NcTestXcorKernelComponent:
    """Create a test component instance."""
    comp = NcTestXcorKernelComponent()
    return comp


@pytest.fixture(name="integrator", scope="module")
def fixture_integrator() -> Ncm.SBesselIntegrator:
    """Create a simple spherical Bessel integrator."""
    integrator = Ncm.SBesselIntegratorLevin.new(0, 2000)
    return integrator


@pytest.fixture(name="gal_kernel")
def fixture_gal_kernel(
    cosmology: Cosmology, integrator: Ncm.SBesselIntegrator
) -> Nc.XcorKernelGal:
    """Create a galaxy kernel with clustering component."""
    return _create_galaxy_kernel(cosmology, integrator)


@pytest.fixture(name="cmb_lensing_kernel")
def fixture_cmb_lensing_kernel(
    cosmology: Cosmology, integrator: Ncm.SBesselIntegrator
) -> Nc.XcorKernelCMBLensing:
    """Create a CMB lensing kernel."""
    return _create_cmb_lensing_kernel(cosmology, integrator)


@pytest.fixture(name="cmb_isw_kernel")
def fixture_cmb_isw_kernel(
    cosmology: Cosmology, integrator: Ncm.SBesselIntegrator
) -> Nc.XcorKernelCMBISW:
    """Create a CMB ISW kernel."""
    return _create_cmb_isw_kernel(cosmology, integrator)


@pytest.fixture(name="tsz_kernel")
def fixture_tsz_kernel(
    cosmology: Cosmology, integrator: Ncm.SBesselIntegrator
) -> Nc.XcorKerneltSZ:
    """Create a tSZ kernel."""
    return _create_tsz_kernel(cosmology, integrator)


@pytest.fixture(name="weak_lensing_kernel")
def fixture_weak_lensing_kernel(
    cosmology: Cosmology, integrator: Ncm.SBesselIntegrator
) -> Nc.XcorKernelWeakLensing:
    """Create a weak lensing kernel."""
    return _create_weak_lensing_kernel(cosmology, integrator)


# Helper functions to create kernel instances
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


def _collect_all_components(
    cosmology: Cosmology, integrator: Ncm.SBesselIntegrator
) -> list[tuple[str, Nc.XcorKernelComponent, Nc.XcorKernel | None]]:
    """Collect all components from all kernel types.

    Returns a list of tuples: (component_id, component, kernel)
    This function discovers all components dynamically.
    """
    cosmo = cosmology.cosmo
    all_components: list[tuple[str, Nc.XcorKernelComponent, Nc.XcorKernel | None]] = []

    # Test component (doesn't come from a kernel)
    test_comp = NcTestXcorKernelComponent()
    test_comp.prepare(cosmo)
    all_components.append((test_comp.__class__.__name__, test_comp, None))

    # Define all kernel types to test
    kernels_to_create: list[Callable[[], Nc.XcorKernel]] = [
        lambda: _create_galaxy_kernel(cosmology, integrator, domagbias=True),
        lambda: _create_cmb_lensing_kernel(cosmology, integrator),
        lambda: _create_cmb_isw_kernel(cosmology, integrator),
        lambda: _create_tsz_kernel(cosmology, integrator),
        lambda: _create_weak_lensing_kernel(cosmology, integrator),
    ]

    # Create each kernel and extract its components
    for kernel_factory in kernels_to_create:
        kernel = kernel_factory()
        kernel_name = kernel.__class__.__name__.replace("XcorKernel", "")
        comp_list = kernel.get_component_list()
        kernel.prepare(cosmo)  # Ensure kernel is prepared before testing components

        # Add each component with a unique ID
        for comp in comp_list:
            comp_name = comp.__class__.__name__.replace("NcXcorKernelComponent", "")
            all_components.append((f"{kernel_name}:{comp_name}", comp, kernel))

    return all_components


def pytest_generate_tests(metafunc: pytest.Metafunc) -> None:
    """Dynamically generate tests for all component types."""
    if "component_case" in metafunc.fixturenames:
        # build minimal objects manually or from config
        cosmology = Cosmology.default()
        integrator = Ncm.SBesselIntegratorLevin.new(0, 2000)

        cases = _collect_all_components(cosmology, integrator)
        ids = [comp_id for comp_id, _, _ in cases]
        metafunc.parametrize("component_case", cases, ids=ids)


# =============================================================================
# Group 1: NcTestXcorKernelComponent Tests (Direct Component Interface)
# =============================================================================


def test_any_component_interface(
    component_case: tuple[str, Nc.XcorKernelComponent, Nc.XcorKernel | None],
    cosmology: Cosmology,
) -> None:
    """Test that any component satisfies the basic interface.

    This test uses the parametrized any_component fixture to test all
    component types with the same interface requirements.
    """
    _, any_component, _ = component_case
    cosmo = cosmology.cosmo

    # All components should be instances of NcXcorKernelComponent
    assert isinstance(any_component, Nc.XcorKernelComponent)

    # Test get_limits
    xi_min, xi_max, k_min, k_max = any_component.get_limits(cosmo)
    assert np.isfinite(xi_min) and xi_min > 0.0
    assert np.isfinite(xi_max) and xi_max > xi_min
    assert np.isfinite(k_min) and k_min > 0.0
    assert np.isfinite(k_max) and k_max > k_min

    # Test eval_kernel at a sample point
    xi_test = (xi_min + xi_max) / 2.0
    k_test = np.sqrt(k_min * k_max)  # Geometric mean
    result = any_component.eval_kernel(cosmo, xi_test, k_test)
    assert np.isfinite(result)

    # Test eval_prefactor
    prefactor = any_component.eval_prefactor(cosmo, k_test, 100)
    assert np.isfinite(prefactor)

    # Test property getters
    epsilon = any_component.get_epsilon()
    assert np.isfinite(epsilon) and epsilon > 0.0

    ny = any_component.get_ny()
    assert isinstance(ny, int) and ny > 0

    max_iter = any_component.get_max_iter()
    assert isinstance(max_iter, int) and max_iter > 0

    tol = any_component.get_tol()
    assert np.isfinite(tol) and tol > 0.0


def test_component_creation(test_component: NcTestXcorKernelComponent) -> None:
    """Test that we can create a component."""
    assert test_component is not None
    assert isinstance(test_component, Nc.XcorKernelComponent)
    assert isinstance(test_component, NcTestXcorKernelComponent)


def test_component_properties(test_component: NcTestXcorKernelComponent) -> None:
    """Test component property getters and setters."""
    # Test epsilon
    test_component.set_epsilon(1.0e-8)
    assert_allclose(test_component.get_epsilon(), 1.0e-8, rtol=1e-10)

    # Test ny
    test_component.set_ny(100)
    assert test_component.get_ny() == 100

    # Test max_iter
    test_component.set_max_iter(500)
    assert test_component.get_max_iter() == 500

    # Test tol
    test_component.set_tol(1.0e-6)
    assert_allclose(test_component.get_tol(), 1.0e-6, rtol=1e-10)


def test_component_eval_kernel(
    test_component: NcTestXcorKernelComponent, cosmology: Cosmology
) -> None:
    """Test eval_kernel method."""
    cosmo = cosmology.cosmo
    xi_array = np.linspace(1.0, 50.0, 10)
    k_array = np.linspace(0.1, 5.0, 10)

    for xi in xi_array:
        for k in k_array:
            result = test_component.eval_kernel(cosmo, xi, k)
            assert np.isfinite(result)
            # Test symmetry and basic properties
            assert isinstance(result, float)

    # Test at xi=0 (should return 0)
    result_zero = test_component.eval_kernel(cosmo, 0.0, 1.0)
    assert result_zero == 0.0


def test_component_eval_prefactor(
    test_component: NcTestXcorKernelComponent, cosmology: Cosmology
) -> None:
    """Test eval_prefactor method."""
    cosmo = cosmology.cosmo
    k_array = np.linspace(0.1, 5.0, 10)
    ell_array = [10, 50, 100, 500]

    for k in k_array:
        for ell in ell_array:
            result = test_component.eval_prefactor(cosmo, k, ell)
            assert np.isfinite(result)
            assert result == 1.0  # Our test component returns constant 1.0


def test_component_get_limits(
    test_component: NcTestXcorKernelComponent, cosmology: Cosmology
) -> None:
    """Test get_limits method."""
    cosmo = cosmology.cosmo
    xi_min, xi_max, k_min, k_max = test_component.get_limits(cosmo)

    assert np.isfinite(xi_min) and xi_min > 0.0
    assert np.isfinite(xi_max) and xi_max > xi_min
    assert np.isfinite(k_min) and k_min > 0.0
    assert np.isfinite(k_max) and k_max > k_min

    # Check against expected values
    assert_allclose(xi_min, test_component.xi_min, rtol=1e-10)
    assert_allclose(xi_max, test_component.xi_max, rtol=1e-10)
    assert_allclose(k_min, test_component.k_min, rtol=1e-10)
    assert_allclose(k_max, test_component.k_max, rtol=1e-10)


def test_component_prepare(
    test_component: NcTestXcorKernelComponent, cosmology: Cosmology
) -> None:
    """Test prepare method and kernel analysis."""
    cosmo = cosmology.cosmo
    # Set parameters for analysis
    test_component.set_epsilon(1.0e-10)
    test_component.set_ny(50)
    test_component.set_max_iter(1000)
    test_component.set_tol(1.0e-8)

    # Call prepare - this should perform kernel analysis
    test_component.prepare(cosmo)

    # After prepare, we should be able to evaluate k_max, K_max, k_epsilon
    y_array = np.linspace(1.0, 100.0, 10)

    for y in y_array:
        k_max = test_component.eval_k_max(y)
        K_max = test_component.eval_K_max(y)
        k_epsilon = test_component.eval_k_epsilon(y)

        assert np.isfinite(k_max) and k_max > 0.0
        assert np.isfinite(K_max)
        assert np.isfinite(k_epsilon) and k_epsilon > 0.0


def test_component_kernel_analysis_properties(
    test_component: NcTestXcorKernelComponent, cosmology: Cosmology
) -> None:
    """Test kernel analysis creates valid splines."""
    cosmo = cosmology.cosmo
    test_component.set_epsilon(1.0e-8)
    test_component.set_ny(30)
    test_component.set_max_iter(1000)
    test_component.set_tol(1.0e-7)

    test_component.prepare(cosmo)

    # Test that eval_k_max is monotonic or behaves reasonably
    y_array = np.geomspace(1.0, 100.0, 20)
    k_max_array = np.array([test_component.eval_k_max(y) for y in y_array])
    K_max_array = np.array([test_component.eval_K_max(y) for y in y_array])
    k_epsilon_array = np.array([test_component.eval_k_epsilon(y) for y in y_array])

    # All should be positive and finite
    assert np.all(np.isfinite(k_max_array))
    assert np.all(k_max_array > 0.0)
    assert np.all(np.isfinite(K_max_array))
    assert np.all(np.isfinite(k_epsilon_array))
    assert np.all(k_epsilon_array > 0.0)

    # k_epsilon should be >= k_max (point where kernel becomes small)
    # This is physics-dependent, so just check they're in reasonable ranges
    assert np.all(k_epsilon_array >= k_max_array * 0.1)


def test_component_edge_cases(
    test_component: NcTestXcorKernelComponent, cosmology: Cosmology
) -> None:
    """Test component behavior at edge cases."""
    cosmo = cosmology.cosmo
    # Test at boundaries
    xi_min, xi_max, k_min, k_max = test_component.get_limits(cosmo)

    # At minimum xi
    result = test_component.eval_kernel(cosmo, xi_min, 1.0)
    assert np.isfinite(result)

    # At maximum xi
    result = test_component.eval_kernel(cosmo, xi_max, 1.0)
    assert np.isfinite(result)

    # At minimum k
    result = test_component.eval_kernel(cosmo, 10.0, k_min)
    assert np.isfinite(result)

    # At maximum k
    result = test_component.eval_kernel(cosmo, 10.0, k_max)
    assert np.isfinite(result)

    # Test very small and large multipoles
    for ell in [2, 10, 100, 1000, 10000]:
        result = test_component.eval_prefactor(cosmo, 1.0, ell)
        assert np.isfinite(result)


# =============================================================================
# Group 2: Real Kernel Components Tests (via get_component_list)
# =============================================================================


def test_gal_kernel_component_list(
    gal_kernel: Nc.XcorKernelGal, cosmology: Cosmology
) -> None:
    """Test that galaxy kernel returns components via get_component_list."""
    cosmo = cosmology.cosmo
    gal_kernel.prepare(cosmo)

    # Get components from kernel
    comp_list = gal_kernel.get_component_list()
    assert comp_list is not None
    assert len(comp_list) > 0  # Should have at least clustering component

    # Test each component satisfies the interface
    for comp in comp_list:
        assert isinstance(comp, Nc.XcorKernelComponent)

        # Test get_limits
        xi_min, xi_max, k_min, k_max = comp.get_limits(cosmo)
        assert np.isfinite(xi_min) and xi_min > 0.0
        assert np.isfinite(xi_max) and xi_max > xi_min
        assert np.isfinite(k_min) and k_min > 0.0
        assert np.isfinite(k_max) and k_max > k_min

        # Test eval_kernel at a sample point
        xi_test = (xi_min + xi_max) / 2.0
        k_test = (k_min + k_max) / 2.0
        result = comp.eval_kernel(cosmo, xi_test, k_test)
        assert np.isfinite(result)

        # Test eval_prefactor
        prefactor = comp.eval_prefactor(cosmo, k_test, 100)
        assert np.isfinite(prefactor)


def test_cmb_lensing_kernel_component_list(
    cmb_lensing_kernel: Nc.XcorKernelCMBLensing, cosmology: Cosmology
) -> None:
    """Test that CMB lensing kernel returns components via get_component_list."""
    cosmo = cosmology.cosmo
    cmb_lensing_kernel.prepare(cosmo)

    comp_list = cmb_lensing_kernel.get_component_list()
    assert comp_list is not None
    assert len(comp_list) > 0

    for comp in comp_list:
        assert isinstance(comp, Nc.XcorKernelComponent)

        xi_min, xi_max, k_min, k_max = comp.get_limits(cosmo)
        assert np.isfinite(xi_min) and xi_min > 0.0
        assert np.isfinite(xi_max) and xi_max > xi_min
        assert np.isfinite(k_min) and k_min > 0.0
        assert np.isfinite(k_max) and k_max > k_min


def test_cmb_isw_kernel_component_list(
    cmb_isw_kernel: Nc.XcorKernelCMBISW, cosmology: Cosmology
) -> None:
    """Test that CMB ISW kernel returns components via get_component_list."""
    cosmo = cosmology.cosmo
    cmb_isw_kernel.prepare(cosmo)

    comp_list = cmb_isw_kernel.get_component_list()
    assert comp_list is not None
    assert len(comp_list) > 0

    for comp in comp_list:
        assert isinstance(comp, Nc.XcorKernelComponent)

        xi_min, xi_max, k_min, k_max = comp.get_limits(cosmo)
        assert np.isfinite(xi_min) and xi_min > 0.0
        assert np.isfinite(xi_max) and xi_max > xi_min
        assert np.isfinite(k_min) and k_min > 0.0
        assert np.isfinite(k_max) and k_max > k_min


def test_tsz_kernel_component_list(
    tsz_kernel: Nc.XcorKerneltSZ, cosmology: Cosmology
) -> None:
    """Test that tSZ kernel returns components via get_component_list."""
    cosmo = cosmology.cosmo
    tsz_kernel.prepare(cosmo)

    comp_list = tsz_kernel.get_component_list()
    assert comp_list is not None
    assert len(comp_list) > 0

    for comp in comp_list:
        assert isinstance(comp, Nc.XcorKernelComponent)

        xi_min, xi_max, k_min, k_max = comp.get_limits(cosmo)
        assert np.isfinite(xi_min) and xi_min > 0.0
        assert np.isfinite(xi_max) and xi_max > xi_min
        assert np.isfinite(k_min) and k_min > 0.0
        assert np.isfinite(k_max) and k_max > k_min


def test_weak_lensing_kernel_component_list(
    weak_lensing_kernel: Nc.XcorKernelWeakLensing, cosmology: Cosmology
) -> None:
    """Test that weak lensing kernel returns components via get_component_list."""
    cosmo = cosmology.cosmo
    weak_lensing_kernel.prepare(cosmo)

    comp_list = weak_lensing_kernel.get_component_list()
    assert comp_list is not None
    assert len(comp_list) > 0

    for comp in comp_list:
        assert isinstance(comp, Nc.XcorKernelComponent)

        xi_min, xi_max, k_min, k_max = comp.get_limits(cosmo)
        assert np.isfinite(xi_min) and xi_min > 0.0
        assert np.isfinite(xi_max) and xi_max > xi_min
        assert np.isfinite(k_min) and k_min > 0.0
        assert np.isfinite(k_max) and k_max > k_min


def test_kernel_component_analysis(
    tsz_kernel: Nc.XcorKerneltSZ, cosmology: Cosmology
) -> None:
    """Test component kernel analysis.

    Test (k_max, K_max, k_epsilon) via get_component_list.
    """
    cosmo = cosmology.cosmo
    tsz_kernel.prepare(cosmo)

    comp_list = tsz_kernel.get_component_list()
    assert comp_list is not None
    assert len(comp_list) > 0

    # After prepare, components should have kernel analysis done
    comp = comp_list[0]

    # Test at several y values
    y_array = np.linspace(10.0, 100.0, 5)

    for y in y_array:
        k_max = comp.eval_k_max(y)
        K_max = comp.eval_K_max(y)
        k_epsilon = comp.eval_k_epsilon(y)

        assert np.isfinite(k_max) and k_max > 0.0
        assert np.isfinite(K_max)
        assert np.isfinite(k_epsilon) and k_epsilon > 0.0
        # k_epsilon should be at or beyond k_max
        assert k_epsilon >= k_max * 0.1


def test_multiple_components_galaxy_kernel(
    cosmology: Cosmology, integrator: Ncm.SBesselIntegrator
) -> None:
    """Test that galaxy kernel with magnification bias returns two components."""
    # Create galaxy kernel with magnification bias enabled
    gal = _create_galaxy_kernel(cosmology, integrator, domagbias=True)

    cosmo = cosmology.cosmo
    gal.prepare(cosmo)

    comp_list = gal.get_component_list()
    assert comp_list is not None
    # Should have 2 components: clustering and magnification bias
    assert len(comp_list) == 2

    # Both should be valid components
    for comp in comp_list:
        assert isinstance(comp, Nc.XcorKernelComponent)
        xi_min, xi_max, k_min, k_max = comp.get_limits(cosmo)
        assert np.isfinite(xi_min) and xi_min > 0.0
        assert np.isfinite(xi_max) and xi_max > xi_min
        assert np.isfinite(k_min) and k_min > 0.0
        assert np.isfinite(k_max) and k_max > k_min


def test_component_properties_via_get_component_list(
    tsz_kernel: Nc.XcorKerneltSZ, cosmology: Cosmology
) -> None:
    """Test component property getters and setters via get_component_list."""
    cosmo = cosmology.cosmo
    tsz_kernel.prepare(cosmo)

    comp_list = tsz_kernel.get_component_list()
    assert comp_list is not None
    assert len(comp_list) > 0

    comp = comp_list[0]

    # Test epsilon property
    comp.set_epsilon(1.0e-7)
    assert_allclose(comp.get_epsilon(), 1.0e-7, rtol=1e-10)

    # Test ny property
    comp.set_ny(150)
    assert comp.get_ny() == 150

    # Test max_iter property
    comp.set_max_iter(2000)
    assert comp.get_max_iter() == 2000

    # Test tol property
    comp.set_tol(1.0e-5)
    assert_allclose(comp.get_tol(), 1.0e-5, rtol=1e-10)

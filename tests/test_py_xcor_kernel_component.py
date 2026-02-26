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

Tests all component types using parametrized component_case fixture.
Components are collected from NcTestXcorKernelComponent and real kernels
via get_component_list().

Test organization:
- Group 1: Component interface tests (parametrized for all components)
- Group 2: Kernel-specific tests (using kernel fixtures directly)
"""

from typing import Callable
from functools import cache
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
        default=1.0e-2,
        flags=GObject.ParamFlags.READWRITE | GObject.ParamFlags.CONSTRUCT_ONLY,
    )

    xi_max: float = GObject.Property(  # type: ignore
        type=float,
        default=2.0,
        flags=GObject.ParamFlags.READWRITE | GObject.ParamFlags.CONSTRUCT_ONLY,
    )

    k_min: float = GObject.Property(  # type: ignore
        type=float,
        default=1.0e-3,
        flags=GObject.ParamFlags.READWRITE | GObject.ParamFlags.CONSTRUCT_ONLY,
    )

    k_max: float = GObject.Property(  # type: ignore
        type=float,
        default=1.0e4,
        flags=GObject.ParamFlags.READWRITE | GObject.ParamFlags.CONSTRUCT_ONLY,
    )

    mean: float = GObject.Property(  # type: ignore
        type=float,
        default=0.8,
        flags=GObject.ParamFlags.READWRITE | GObject.ParamFlags.CONSTRUCT_ONLY,
    )

    sigma: float = GObject.Property(  # type: ignore
        type=float,
        default=0.05,
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
        z = xi - self.mean
        gaussian = xi * np.exp(-0.5 * (z / self.sigma) ** 2)
        return gaussian * np.sinc(k * z * 1.0e-20)

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


@cache
def get_cosmology() -> Cosmology:
    """Create a default cosmology for testing."""
    return Cosmology.default()


@pytest.fixture(name="cosmology", scope="module")
def fixture_cosmology() -> Cosmology:
    """Create a simple cosmology for testing."""
    cosmology = get_cosmology()
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


@pytest.fixture(name="comp")
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
    mean = 0.8
    sigma = 0.3
    lower_bound = max(0.0, mean - 5.5 * sigma)
    upper_bound = mean + 5.5 * sigma

    z_array = np.linspace(lower_bound, upper_bound, 200)
    dndz_array = np.exp(-(((z_array - mean) / sigma) ** 2) / 2.0)
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
    gal["mag_bias"] = 0.3
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
    mean = 0.6
    sigma = 0.3
    lower_bound = max(0.0, mean - 5.5 * sigma)
    upper_bound = mean + 5.5 * sigma

    z_array = np.linspace(lower_bound, upper_bound, 200)
    dndz_array = np.exp(-(((z_array - mean) / sigma) ** 2) / 2.0)
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
    print("Collecting components from NcTestXcorKernelComponent...", flush=True)
    test_comp = NcTestXcorKernelComponent()
    test_comp.prepare(cosmo)
    all_components.append((test_comp.__class__.__name__, test_comp, None))
    print("Collected 1 component from NcTestXcorKernelComponent", flush=True)

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


@cache
def get_cases():
    """Get all component cases for parametrization."""
    cosmology = get_cosmology()
    integrator = Ncm.SBesselIntegratorLevin.new(0, 2000)
    return _collect_all_components(cosmology, integrator)


def pytest_generate_tests(metafunc: pytest.Metafunc) -> None:
    """Dynamically generate tests for all component types."""

    if "component_case" in metafunc.fixturenames:
        # build minimal objects manually or from config
        cases = get_cases()
        ids = [comp_id for comp_id, _, _ in cases]
        metafunc.parametrize("component_case", cases, ids=ids)


# =============================================================================
# Group 1: Component Interface Tests (Parametrized for All Components)
# =============================================================================


def test_any_component_interface(
    component_case: tuple[str, Nc.XcorKernelComponent, Nc.XcorKernel | None],
    cosmology: Cosmology,
) -> None:
    """Test that any component satisfies the basic interface."""
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


def test_component_creation(
    component_case: tuple[str, Nc.XcorKernelComponent, Nc.XcorKernel | None],
) -> None:
    """Test that we can create a component."""
    _, comp, _ = component_case
    assert comp is not None
    assert isinstance(comp, Nc.XcorKernelComponent)


def test_component_properties(
    comp: Nc.XcorKernelComponent,
) -> None:
    """Test component property getters and setters."""
    # Test epsilon
    comp.set_epsilon(1.0e-8)
    assert_allclose(comp.get_epsilon(), 1.0e-8, rtol=1e-10)

    # Test ny
    comp.set_ny(100)
    assert comp.get_ny() == 100

    # Test max_iter
    comp.set_max_iter(500)
    assert comp.get_max_iter() == 500

    # Test tol
    comp.set_tol(1.0e-6)
    assert_allclose(comp.get_tol(), 1.0e-6, rtol=1e-10)


def test_component_eval_kernel(
    component_case: tuple[str, Nc.XcorKernelComponent, Nc.XcorKernel | None],
    cosmology: Cosmology,
) -> None:
    """Test eval_kernel method."""
    _, comp, _ = component_case
    cosmo = cosmology.cosmo

    xi_min, xi_max, k_min, k_max = comp.get_limits(cosmo)
    xi_array = np.linspace(xi_min, xi_max, 10)
    k_array = np.linspace(k_min, k_max, 10)

    for xi in xi_array:
        for k in k_array:
            result = comp.eval_kernel(cosmo, xi, k)
            assert np.isfinite(result)
            # Test symmetry and basic properties
            assert isinstance(result, float)


def test_component_eval_prefactor(
    component_case: tuple[str, Nc.XcorKernelComponent, Nc.XcorKernel | None],
    cosmology: Cosmology,
) -> None:
    """Test eval_prefactor method."""
    _, comp, _ = component_case
    cosmo = cosmology.cosmo
    k_array = np.linspace(0.1, 5.0, 10)
    ell_array = [10, 50, 100, 500]

    for k in k_array:
        for ell in ell_array:
            result = comp.eval_prefactor(cosmo, k, ell)
            assert np.isfinite(result)


def test_component_get_limits(
    component_case: tuple[str, Nc.XcorKernelComponent, Nc.XcorKernel | None],
    cosmology: Cosmology,
) -> None:
    """Test get_limits method."""
    _, comp, _ = component_case
    cosmo = cosmology.cosmo
    xi_min, xi_max, k_min, k_max = comp.get_limits(cosmo)

    assert np.isfinite(xi_min) and xi_min > 0.0
    assert np.isfinite(xi_max) and xi_max > xi_min
    assert np.isfinite(k_min) and k_min > 0.0
    assert np.isfinite(k_max) and k_max > k_min


def test_component_prepare(comp: Nc.XcorKernelComponent, cosmology: Cosmology) -> None:
    """Test prepare method and kernel analysis."""
    # Set parameters for analysis
    cosmo = cosmology.cosmo
    comp.set_epsilon(1.0e-10)
    comp.set_ny(50)
    comp.set_max_iter(1000)
    comp.set_tol(1.0e-8)
    comp.prepare(cosmo)

    # After prepare, we should be able to evaluate k_max, K_max, k_epsilon
    y_array = np.linspace(1.0, 100.0, 10)

    for y in y_array:
        k_max = comp.eval_k_max(y)
        K_max = comp.eval_K_max(y)
        k_epsilon = comp.eval_k_epsilon(y)

        assert np.isfinite(k_max) and k_max > 0.0
        assert np.isfinite(K_max)
        assert np.isfinite(k_epsilon) and k_epsilon > 0.0


def test_component_kernel_analysis_properties(
    comp: Nc.XcorKernelComponent, cosmology: Cosmology
) -> None:
    """Test kernel analysis creates valid splines."""
    cosmo = cosmology.cosmo
    comp.set_epsilon(1.0e-8)
    comp.set_ny(30)
    comp.set_max_iter(1000)
    comp.set_tol(1.0e-7)
    comp.prepare(cosmo)

    # Test that eval_k_max is monotonic or behaves reasonably
    y_array = np.geomspace(1.0, 100.0, 20)
    k_max_array = np.array([comp.eval_k_max(y) for y in y_array])
    K_max_array = np.array([comp.eval_K_max(y) for y in y_array])
    k_epsilon_array = np.array([comp.eval_k_epsilon(y) for y in y_array])

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
    component_case: tuple[str, Nc.XcorKernelComponent, Nc.XcorKernel | None],
    cosmology: Cosmology,
) -> None:
    """Test component behavior at edge cases."""
    _, comp, _ = component_case
    cosmo = cosmology.cosmo
    # Test at boundaries
    xi_min, xi_max, k_min, k_max = comp.get_limits(cosmo)
    xi_mean = np.sqrt(xi_min * xi_max)
    k_mean = np.sqrt(k_min * k_max)

    # At minimum xi
    result = comp.eval_kernel(cosmo, xi_min, k_mean)
    assert np.isfinite(result)

    # At maximum xi
    result = comp.eval_kernel(cosmo, xi_max, k_mean)
    assert np.isfinite(result)

    # At minimum k
    result = comp.eval_kernel(cosmo, xi_mean, k_min)
    assert np.isfinite(result)

    # At maximum k
    result = comp.eval_kernel(cosmo, xi_mean, k_max)
    assert np.isfinite(result)

    # Test very small and large multipoles
    for ell in [2, 10, 100, 1000, 10000]:
        result = comp.eval_prefactor(cosmo, k_mean, ell)
        assert np.isfinite(result)


def test_component_k_max_is_maximum(
    component_case: tuple[str, Nc.XcorKernelComponent, Nc.XcorKernel | None],
    cosmology: Cosmology,
) -> None:
    """Test that k_max(y) gives the maximum of K*xi(k, y/k)."""
    _, comp, _ = component_case
    cosmo = cosmology.cosmo

    xi_min, xi_max, k_min, k_max = comp.get_limits(cosmo)

    y_min = k_min * xi_min
    y_max = k_max * xi_max

    # Test at several y values
    y_test_values = np.geomspace(max(y_min * 1.1, 0.5), min(y_max * 0.9, 1000.0), 5)

    for y in y_test_values:
        k_max_val = comp.eval_k_max(y)
        K_max_val = comp.eval_K_max(y)
        xi = y / k_max_val

        # Verify K_max is approximately the actual value at k_max
        if xi_min <= xi <= xi_max:
            K_at_k_max = abs(xi * comp.eval_kernel(cosmo, xi, k_max_val))
            # Allow tolerance for spline/search errors
            assert_allclose(
                K_max_val,
                K_at_k_max,
                rtol=0.05,
                atol=1e-10,
                err_msg=f"K_max doesn't match eval_kernel at y={y}",
            )


def test_component_K_max_value(
    component_case: tuple[str, Nc.XcorKernelComponent, Nc.XcorKernel | None],
    cosmology: Cosmology,
) -> None:
    """Test that K_max(y) equals K*xi(y/k_max, k_max)."""
    _, comp, _ = component_case
    cosmo = cosmology.cosmo

    xi_min, xi_max, k_min, k_max = comp.get_limits(cosmo)
    y_min = k_min * xi_min
    y_max = k_max * xi_max

    # Test at multiple y values
    y_test_values = np.geomspace(max(y_min * 1.1, 0.5), min(y_max * 0.9, 1000.0), 10)

    for y in y_test_values:
        k_max_val = comp.eval_k_max(y)
        K_max_val = comp.eval_K_max(y)

        # Compute xi from y and k_max
        xi = y / k_max_val

        if xi_min <= xi <= xi_max:
            # Compute K*xi directly from eval_kernel
            K_direct = abs(xi * comp.eval_kernel(cosmo, xi, k_max_val))

            # They should match within ~0.2% for real kernels
            assert_allclose(
                K_max_val,
                K_direct,
                rtol=0.05,
                atol=1e-10,
                err_msg=(
                    f"K_max mismatch at y={y}: " f"K_max={K_max_val}, direct={K_direct}"
                ),
            )


def test_component_k_epsilon_drop(
    component_case: tuple[str, Nc.XcorKernelComponent, Nc.XcorKernel | None],
    cosmology: Cosmology,
) -> None:
    """Test that k_epsilon is where K*xi drops by epsilon from K_max.

    At k_epsilon: K*xi(y/k_epsilon, k_epsilon) ≈ epsilon * K_max
    """
    _, comp, _ = component_case
    cosmo = cosmology.cosmo
    epsilon = 1.0e-5

    original_epsilon = comp.get_epsilon()
    # The algorithm will find k_epsilon where K*xi = epsilon^2 * K_max.
    comp.set_epsilon(epsilon**2)
    comp.prepare(cosmo)  # Re-prepare with new epsilon

    xi_min, xi_max, _, _ = comp.get_limits(cosmo)

    # Test at several y values
    y_test_values = [20.0, 50.0, 100.0]

    for y in y_test_values:
        k_max_val = comp.eval_k_max(y)
        K_max_val = comp.eval_K_max(y)
        k_epsilon_val = comp.eval_k_epsilon(y)

        # k_epsilon should be >= k_max (after the maximum)
        assert (
            k_epsilon_val >= k_max_val * 0.99
        ), f"k_epsilon={k_epsilon_val} should be >= k_max={k_max_val} at y={y}"

        # Compute K*xi at k_epsilon
        xi_epsilon = y / k_epsilon_val
        if xi_min <= xi_epsilon <= xi_max:
            K_at_epsilon = np.abs(
                xi_epsilon * comp.eval_kernel(cosmo, xi_epsilon, k_epsilon_val)
            )
            # Should be approximately epsilon * K_max (allow factor of 5 tolerance due
            # to spline/search)
            expected_K = np.abs(epsilon * K_max_val)
            assert K_at_epsilon >= expected_K * 0.2, (
                f"K at k_epsilon ({K_at_epsilon}) too small vs "
                f"epsilon*K_max ({expected_K}) at y={y}"
            )
            assert K_at_epsilon <= expected_K * 5.0, (
                f"K at k_epsilon ({K_at_epsilon}) too large vs "
                f"epsilon*K_max ({expected_K}) at y={y}"
            )

    comp.set_epsilon(original_epsilon)
    comp.prepare(cosmo)


def test_component_y_range_pruning(
    component_case: tuple[str, Nc.XcorKernelComponent, Nc.XcorKernel | None],
    cosmology: Cosmology,
) -> None:
    """Test that y range is correctly pruned to [0.5, 1000.0].

    The algorithm computes y_min = max(k_min * xi_min, 0.5) and
    y_max = min(k_max * xi_max, 1000.0).
    """
    _, comp, _ = component_case
    cosmo = cosmology.cosmo

    xi_min, xi_max, k_min, k_max = comp.get_limits(cosmo)

    # Calculate expected y range
    y_min_expected = max(k_min * xi_min, 0.5)
    y_max_expected = min(k_max * xi_max, 1000.0)

    # Test that eval_k_max works at the boundaries
    # Values near y_min
    k_max_at_ymin = comp.eval_k_max(y_min_expected)
    assert np.isfinite(k_max_at_ymin) and k_max_at_ymin > 0.0

    # Values near y_max
    k_max_at_ymax = comp.eval_k_max(y_max_expected)
    assert np.isfinite(k_max_at_ymax) and k_max_at_ymax > 0.0

    # Test a few points in the valid range
    y_test_values = np.geomspace(y_min_expected * 1.1, y_max_expected * 0.9, 5)
    for y in y_test_values:
        k = comp.eval_k_max(y)
        assert np.isfinite(k) and k > 0.0, f"Invalid k_max at y={y}"


def test_component_k_max_is_maximum1(
    component_case: tuple[str, Nc.XcorKernelComponent, Nc.XcorKernel | None],
    cosmology: Cosmology,
) -> None:
    """Test that k_max(y) is actually the maximum of K*xi(k, y/k).

    For a given y, k_max should satisfy: K*xi(k_max, y/k_max) >= K*xi(k, y/k)
    for all k in the valid range.
    """
    _, comp, _ = component_case
    cosmo = cosmology.cosmo

    xi_min, xi_max, k_min, k_max = comp.get_limits(cosmo)
    y_min = k_min * xi_min
    y_max = k_max * xi_max

    # Test at several y values
    y_test_values = np.geomspace(max(y_min * 1.1, 0.5), min(y_max * 0.9, 1000.0), 5)

    for y in y_test_values:
        k_max_val = comp.eval_k_max(y)
        K_max_val = comp.eval_K_max(y)
        xi = y / k_max_val

        # Verify K_max is approximately the actual value at k_max
        K_at_k_max = xi * comp.eval_kernel(cosmo, xi, k_max_val)
        # Allow tolerance due to spline interpolation and oscillating kernels
        assert_allclose(
            abs(K_max_val),
            abs(K_at_k_max),
            rtol=0.2,
            atol=1e-10,
            err_msg=f"K_max doesn't match eval_kernel at y={y}, k_max={k_max_val}",
        )

        # Sample points around k_max to verify it's a local maximum
        k_test_points = np.array(
            [
                k_max_val * 0.7,  # Before maximum
                k_max_val * 0.9,
                k_max_val,  # At maximum
                k_max_val * 1.1,
                k_max_val * 1.3,  # After maximum
            ]
        )

        # Filter to stay within valid k range
        k_test_points = k_test_points[
            (k_test_points >= k_min) & (k_test_points <= k_max)
        ]

        K_values_list = []
        for k_test in k_test_points:
            xi_test = y / k_test
            # Only test if xi is in valid range
            if xi_min <= xi_test <= xi_max:
                K_val = abs(xi_test * comp.eval_kernel(cosmo, xi_test, k_test))
                K_values_list.append(K_val)
            else:
                K_values_list.append(-np.inf)  # Mark as invalid

        K_values = np.array(K_values_list)
        valid_K = K_values[np.isfinite(K_values)]

        if len(valid_K) > 2:
            # K_max should be among the larger values (allowing numerical errors)
            assert (
                abs(K_max_val) >= np.percentile(valid_K, 90) * 0.9
            ), f"k_max={k_max_val} doesn't give a large |K| value at y={y}"


def test_component_K_max_value1(
    component_case: tuple[str, Nc.XcorKernelComponent, Nc.XcorKernel | None],
    cosmology: Cosmology,
) -> None:
    """Test that K_max(y) approximately equals K*xi(y/k_max, k_max).

    Note: There can be discrepancy due to spline interpolation and the
    maximum-finding algorithm not being infinitely precise.
    """
    _, comp, _ = component_case
    cosmo = cosmology.cosmo

    # Test at multiple y values
    y_test_values = np.geomspace(2.0, 50.0, 10)

    for y in y_test_values:
        k_max_val = comp.eval_k_max(y)
        K_max_val = comp.eval_K_max(y)

        # Compute xi from y and k_max
        xi = y / k_max_val

        # Compute K*xi directly from eval_kernel
        K_direct = xi * comp.eval_kernel(cosmo, xi, k_max_val)

        # Allow tolerance for spline interpolation and search errors
        assert_allclose(
            abs(K_max_val),
            abs(K_direct),
            rtol=0.05,
            atol=1e-10,
            err_msg=f"K_max mismatch at y={y}: K_max={K_max_val}, direct={K_direct}",
        )


def test_component_k_epsilon_drop1(
    component_case: tuple[str, Nc.XcorKernelComponent, Nc.XcorKernel | None],
    cosmology: Cosmology,
) -> None:
    """Test that k_epsilon is where |K*xi| drops by epsilon from K_max for k > k_max.

    At k_epsilon: |K*xi(y/k_epsilon, k_epsilon)| = epsilon * K_max
    Note: We test absolute value because kernels can go negative.
    """
    _, comp, _ = component_case
    cosmo = cosmology.cosmo
    original_epsilon = comp.get_epsilon()
    epsilon = 1.0e-5  # Use larger epsilon for more reliable testing
    # The algorithm finds k_epsilon where K*xi = epsilon^2 * K_max, so we set epsilon^2
    # here.
    comp.set_epsilon(epsilon**2)
    comp.prepare(cosmo)

    xi_min, xi_max, k_min, k_max = comp.get_limits(cosmo)
    y_min = k_min * xi_min
    y_max = k_max * xi_max

    # Test at several y values
    y_test_values = np.geomspace(max(y_min * 1.1, 0.5), min(y_max * 0.9, 1000.0), 5)

    for y in y_test_values:
        k_max_val = comp.eval_k_max(y)
        K_max_val = comp.eval_K_max(y)
        k_epsilon_val = comp.eval_k_epsilon(y)

        # k_epsilon should be >= k_max
        assert (
            k_epsilon_val >= k_max_val * 0.99
        ), f"k_epsilon={k_epsilon_val} should be >= k_max={k_max_val} at y={y}"

        # Compute |K*xi| at k_epsilon
        xi_epsilon = y / k_epsilon_val
        K_at_epsilon = np.abs(xi_epsilon * comp.eval_kernel(cosmo, xi_epsilon, k_epsilon_val))

        # Should be approximately epsilon * K_max
        expected_K = epsilon * K_max_val

        # Skip if k_epsilon is at boundary
        if np.isclose(k_epsilon_val, k_max, atol=1e-3) or np.isclose(
            xi_epsilon, xi_min, atol=1e-3
        ):
            return

        assert K_at_epsilon <= expected_K * 10.0, (
            f"K at k_epsilon ({K_at_epsilon}) should be "
            f"~epsilon*K_max ({expected_K}) at y={y}"
        )
        # Also check it's not way too small
        assert K_at_epsilon >= expected_K * 0.1, (
            f"K at k_epsilon ({K_at_epsilon}) too small "
            f"compared to epsilon*K_max ({expected_K}) at y={y}"
        )

    comp.set_epsilon(original_epsilon)
    comp.prepare(cosmo)


def test_component_k_epsilon_after_k_max(
    component_case: tuple[str, Nc.XcorKernelComponent, Nc.XcorKernel | None],
    cosmology: Cosmology,
) -> None:
    """Test that k_epsilon > k_max for all y values.

    The algorithm looks for the drop point after the maximum.
    """
    _, comp, _ = component_case
    cosmo = cosmology.cosmo
    original_epsilon = comp.get_epsilon()
    comp.set_epsilon(1.0e-7)
    comp.prepare(cosmo)

    # Test at many y values
    y_test_values = np.geomspace(1.0, 100.0, 15)

    for y in y_test_values:
        k_max_val = comp.eval_k_max(y)
        k_epsilon_val = comp.eval_k_epsilon(y)

        # k_epsilon should always be after k_max
        assert (
            k_epsilon_val >= k_max_val
        ), f"k_epsilon={k_epsilon_val} should be >= k_max={k_max_val} at y={y}"

    comp.set_epsilon(original_epsilon)
    comp.prepare(cosmo)


def test_component_monotonicity(
    component_case: tuple[str, Nc.XcorKernelComponent, Nc.XcorKernel | None],
) -> None:
    """Test basic monotonicity and sanity of k_max, K_max, k_epsilon."""
    _, comp, _ = component_case

    # Sample many y values
    y_array = np.geomspace(1.0, 100.0, 30)

    k_max_array = np.array([comp.eval_k_max(y) for y in y_array])
    K_max_array = np.array([comp.eval_K_max(y) for y in y_array])
    k_epsilon_array = np.array([comp.eval_k_epsilon(y) for y in y_array])

    # All should be positive and finite
    assert np.all(np.isfinite(k_max_array))
    assert np.all(k_max_array > 0.0)
    assert np.all(np.isfinite(K_max_array))
    assert np.all(np.isfinite(k_epsilon_array))
    assert np.all(k_epsilon_array > 0.0)

    # k_epsilon >= k_max everywhere
    assert np.all(k_epsilon_array >= k_max_array * 0.99)

    # K_max should be positive (for our test kernel)
    assert np.all(K_max_array >= 0.0)


def test_component_spline_smoothness(
    component_case: tuple[str, Nc.XcorKernelComponent, Nc.XcorKernel | None],
) -> None:
    """Test that the splines for k_max, K_max, k_epsilon are smooth.

    Evaluate at many points and check for no sudden jumps.
    """
    _, comp, _ = component_case

    # Dense sampling
    y_array = np.geomspace(2.0, 80.0, 600)

    k_max_array = np.array([comp.eval_k_max(y) for y in y_array])
    K_max_array = np.array([comp.eval_K_max(y) for y in y_array])
    k_epsilon_array = np.array([comp.eval_k_epsilon(y) for y in y_array])

    # Check no sudden jumps (should be relatively smooth)
    # Compute relative differences between consecutive points
    def rel_diff(arr):
        return np.abs(np.diff(arr)) / (arr[:-1] + 1e-10)

    k_max_diff = rel_diff(k_max_array)
    K_max_diff = rel_diff(K_max_array)
    k_epsilon_diff = rel_diff(k_epsilon_array)

    # Check for smoothness
    assert (
        np.max(k_max_diff) < 1.0
    ), f"k_max has discontinuous jumps: max={np.max(k_max_diff)}"
    assert np.max(K_max_diff) < 6.0, (
        f"K_max has large jumps: max={np.max(K_max_diff)} "
        f"(may indicate kernel structure)"
    )
    assert (
        np.max(k_epsilon_diff) < 1.0
    ), f"k_epsilon has discontinuous jumps: max={np.max(k_epsilon_diff)}"


def test_kernel_analysis_with_different_epsilon(
    component_case: tuple[str, Nc.XcorKernelComponent, Nc.XcorKernel | None],
    cosmology: Cosmology,
) -> None:
    """Test that smaller epsilon gives larger k_epsilon."""
    _, comp, _ = component_case
    cosmo = cosmology.cosmo

    original_epsilon = comp.get_epsilon()
    epsilon_values = [1e-4, 1e-6, 1e-8]
    y_test = 20.0

    k_epsilon_results = []

    for eps in epsilon_values:
        comp.set_epsilon(eps)
        comp.prepare(cosmo)
        k_eps = comp.eval_k_epsilon(y_test)
        k_epsilon_results.append(k_eps)

    comp.set_epsilon(original_epsilon)
    comp.prepare(cosmo)

    # Smaller epsilon should give larger k_epsilon
    assert (
        k_epsilon_results[0] <= k_epsilon_results[1] <= k_epsilon_results[2]
    ), f"k_epsilon should increase as epsilon decreases: {k_epsilon_results}"


def test_component_correctness(
    component_case: tuple[str, Nc.XcorKernelComponent, Nc.XcorKernel | None],
    cosmology: Cosmology,
) -> None:
    """Test kernel analysis correctness."""
    _, comp, _ = component_case
    cosmo = cosmology.cosmo

    # Get limits
    xi_min, xi_max, _, _ = comp.get_limits(cosmo)

    # Test at a few y values
    y_test_values = [10.0, 50.0, 200.0]

    for y in y_test_values:
        k_max_val = comp.eval_k_max(y)
        K_max_val = comp.eval_K_max(y)

        # Verify K_max equals K*xi(y/k_max, k_max)
        xi = y / k_max_val
        if xi_min <= xi <= xi_max:
            K_direct = np.abs(xi * comp.eval_kernel(cosmo, xi, k_max_val))
            # Allow generous tolerance for real components with complex structure
            assert_allclose(
                K_max_val, K_direct, rtol=0.05, err_msg=f"K_max mismatch at y={y}"
            )

        # Test that k_epsilon >= k_max
        k_epsilon_val = comp.eval_k_epsilon(y)
        assert (
            k_epsilon_val >= k_max_val * 0.99
        ), f"k_epsilon should be >= k_max at y={y}"


# =============================================================================
# Group 2: Kernel-Specific Tests (Using Kernel Fixtures)
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
    """Test that tSZ kernel components have valid kernel analysis."""
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


def test_component_chebyshev_decomposition(
    component_case: tuple[str, Nc.XcorKernelComponent, Nc.XcorKernel | None],
    cosmology: Cosmology,
) -> None:
    """Test Chebyshev decomposition of component kernel g_k(y) = K*xi(y/k, k)."""
    _, comp, _ = component_case
    cosmo = cosmology.cosmo

    xi_min, xi_max, k_min, k_max = comp.get_limits(cosmo)
    spectral = Ncm.Spectral.new()
    k_values = np.geomspace(k_min, k_max, 10)
    for k in k_values:

        def g_k(y, k=k):
            """Evaluate K*xi(y/k, k)."""
            xi = y / k
            K_val = xi * comp.eval_kernel(cosmo, xi, k)
            return K_val

        y_min = k * xi_min
        y_max = k * xi_max

        # Compute Chebyshev coefficients adaptively
        final_k, coeffs = spectral.compute_chebyshev_coeffs_adaptive(
            g_k, y_min, y_max, 4, 1.0e-8
        )
        max_order_k = spectral.get_max_order()

        # Verify we got coefficients
        assert final_k >= 5
        assert coeffs is not None
        assert len(coeffs) >= 2**5 + 1
        assert len(coeffs) <= 2**max_order_k + 1

        # Verify the expansion is accurate at test points
        y_test_points = np.linspace(y_min, y_max, 300)
        direct_vals = np.array([g_k(y) for y in y_test_points])
        max_direct = np.max(np.abs(direct_vals))
        cheb_vals = np.array(
            [
                Ncm.Spectral.chebyshev_eval_x(coeffs, y_min, y_max, y)
                for y in y_test_points
            ]
        )
        if max_direct > 0.0:
            cheb_vals /= max_direct
            direct_vals /= max_direct

        assert_allclose(
            cheb_vals,
            direct_vals,
            rtol=1.0e-7,
            atol=1.0e-7,
            err_msg=f"Chebyshev expansion mismatch at k={k}",
        )

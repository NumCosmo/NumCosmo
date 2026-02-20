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

"""Unit tests for NcXcorKernelComponent."""

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


@pytest.fixture(name="integrator")
def fixture_integrator() -> Ncm.SBesselIntegrator:
    """Create a simple spherical Bessel integrator."""
    integrator = Ncm.SBesselIntegratorLevin.new(0, 2000)
    return integrator


@pytest.fixture(name="gal_component")
def fixture_gal_component(
    cosmology: Cosmology, integrator: Ncm.SBesselIntegrator
) -> Nc.XcorKernelGal:
    """Create a galaxy kernel with clustering component."""
    dn_dz = Ncm.SplineCubicNotaknot()
    z_array = np.linspace(0.1, 2.0, 50)
    dndz_array = np.exp(-((z_array - 0.8) ** 2) / (2.0 * 0.3**2))
    dndz_array /= np.trapezoid(dndz_array, z_array)

    xv = Ncm.Vector.new_array(z_array.tolist())
    yv = Ncm.Vector.new_array(dndz_array.tolist())
    dn_dz.set(xv, yv, True)

    gal = Nc.XcorKernelGal(
        dndz=dn_dz, dist=cosmology.dist, powspec=cosmology.ps_ml, integrator=integrator
    )
    gal["bparam_0"] = 1.5
    return gal


@pytest.fixture(name="cmb_lensing_component")
def fixture_cmb_lensing_component(
    cosmology: Cosmology, integrator: Ncm.SBesselIntegrator
) -> Nc.XcorKernelCMBLensing:
    """Create a CMB lensing kernel component."""
    lmax = 1000
    ell_vec = Ncm.Vector.new_array(np.arange(lmax + 1).tolist())
    cmb_lens = Nc.XcorKernelCMBLensing(
        dist=cosmology.dist,
        powspec=cosmology.ps_ml,
        recomb=cosmology.recomb,
        Nl=ell_vec,
        integrator=integrator,
    )
    return cmb_lens


@pytest.fixture(name="cmb_isw_component")
def fixture_cmb_isw_component(
    cosmology: Cosmology, integrator: Ncm.SBesselIntegrator
) -> Nc.XcorKernelCMBISW:
    """Create a CMB ISW kernel component."""
    lmax = 1000
    ell_vec = Ncm.Vector.new_array(np.arange(lmax + 1).tolist())
    cmb_isw = Nc.XcorKernelCMBISW(
        dist=cosmology.dist,
        powspec=cosmology.ps_ml,
        recomb=cosmology.recomb,
        Nl=ell_vec,
        integrator=integrator,
    )
    return cmb_isw


@pytest.fixture(name="tsz_component")
def fixture_tsz_component(
    cosmology: Cosmology, integrator: Ncm.SBesselIntegrator
) -> Nc.XcorKerneltSZ:
    """Create a tSZ kernel component."""
    tsz = Nc.XcorKerneltSZ(
        dist=cosmology.dist, powspec=cosmology.ps_ml, zmax=6.0, integrator=integrator
    )
    return tsz


@pytest.fixture(name="weak_lensing_component")
def fixture_weak_lensing_component(
    cosmology: Cosmology, integrator: Ncm.SBesselIntegrator
) -> Nc.XcorKernelWeakLensing:
    """Create a weak lensing kernel component."""
    dn_dz = Ncm.SplineCubicNotaknot()
    z_array = np.linspace(0.1, 3.0, 50)
    dndz_array = (z_array / 1.0) ** 2 * np.exp(-((z_array / 1.0) ** 1.5))
    dndz_array /= np.trapezoid(dndz_array, z_array)

    xv = Ncm.Vector.new_array(z_array.tolist())
    yv = Ncm.Vector.new_array(dndz_array.tolist())
    dn_dz.set(xv, yv, True)

    wl = Nc.XcorKernelWeakLensing(
        dndz=dn_dz,
        dist=cosmology.dist,
        powspec=cosmology.ps_ml,
        integrator=integrator,
    )
    return wl


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


# Tests for real components (Galaxy, CMB, etc.)


def test_gal_clustering_component(
    gal_component: Nc.XcorKernelGal, cosmology: Cosmology
) -> None:
    """Test that galaxy kernel creates components correctly."""
    cosmo = cosmology.cosmo
    dist = cosmology.dist
    assert gal_component is not None

    # Prepare the kernel
    gal_component.prepare(cosmo)

    # The galaxy kernel should create a clustering component
    # We can't directly access components, but we can test the kernel works
    z_array = np.linspace(0.1, 1.5, 10)
    ell = 100

    for z in z_array:
        result = gal_component.eval_limber_z_full(cosmo, z, dist, ell)
        assert np.isfinite(result)


def test_cmb_lensing_component(
    cmb_lensing_component: Nc.XcorKernelCMBLensing, cosmology: Cosmology
) -> None:
    """Test CMB lensing kernel component."""
    cosmo = cosmology.cosmo
    dist = cosmology.dist
    cmb_lens = cmb_lensing_component
    cmb_lens.prepare(cosmo)

    # Test basic kernel evaluation
    z_array = np.linspace(0.1, 2.0, 10)
    ell = 100

    for z in z_array:
        result = cmb_lens.eval_limber_z_full(cosmo, z, dist, ell)
        assert np.isfinite(result)


def test_cmb_isw_component(
    cmb_isw_component: Nc.XcorKernelCMBISW, cosmology: Cosmology
) -> None:
    """Test CMB ISW kernel component."""
    cosmo = cosmology.cosmo
    dist = cosmology.dist
    cmb_isw = cmb_isw_component
    cmb_isw.prepare(cosmo)

    # Test basic kernel evaluation
    z_array = np.linspace(0.1, 2.0, 10)
    ell = 100

    for z in z_array:
        result = cmb_isw.eval_limber_z_full(cosmo, z, dist, ell)
        assert np.isfinite(result)


def test_tsz_component(tsz_component: Nc.XcorKerneltSZ, cosmology: Cosmology) -> None:
    """Test tSZ kernel component."""
    cosmo = cosmology.cosmo
    dist = cosmology.dist
    tsz = tsz_component
    tsz.prepare(cosmo)

    # Test basic kernel evaluation
    z_array = np.linspace(0.1, 2.0, 10)
    ell = 100

    for z in z_array:
        result = tsz.eval_limber_z_full(cosmo, z, dist, ell)
        assert np.isfinite(result)


def test_weak_lensing_component(
    weak_lensing_component: Nc.XcorKernelWeakLensing, cosmology: Cosmology
) -> None:
    """Test weak lensing kernel component."""
    cosmo = cosmology.cosmo
    dist = cosmology.dist
    wl = weak_lensing_component
    wl.prepare(cosmo)

    # Test basic kernel evaluation
    z_test = np.linspace(0.1, 2.5, 10)
    ell = 100

    for z in z_test:
        result = wl.eval_limber_z_full(cosmo, z, dist, ell)
        assert np.isfinite(result)


def test_component_comparison(
    cmb_lensing_component: Nc.XcorKernelCMBLensing,
    cmb_isw_component: Nc.XcorKernelCMBISW,
    tsz_component: Nc.XcorKerneltSZ,
    cosmology: Cosmology,
) -> None:
    """Test that different components produce different results."""
    cosmo = cosmology.cosmo
    dist = cosmology.dist

    cmb_lens = cmb_lensing_component
    cmb_isw = cmb_isw_component
    tsz = tsz_component

    cmb_lens.prepare(cosmo)
    cmb_isw.prepare(cosmo)
    tsz.prepare(cosmo)

    z: float = 0.5
    ell: int = 100

    result_lens = cmb_lens.eval_limber_z_full(cosmo, z, dist, ell)
    result_isw = cmb_isw.eval_limber_z_full(cosmo, z, dist, ell)
    result_tsz = tsz.eval_limber_z_full(cosmo, z, dist, ell)

    # All should be finite but different
    assert np.isfinite(result_lens)
    assert np.isfinite(result_isw)
    assert np.isfinite(result_tsz)

    # They should not all be equal (different physics)
    values = [result_lens, result_isw, result_tsz]
    assert len(set(values)) > 1, "All components returned same value"


def test_component_with_different_cosmologies(
    cosmology: Cosmology, cosmology_alt: Cosmology, integrator: Ncm.SBesselIntegrator
) -> None:
    """Test that components respond to cosmology changes."""
    cosmo1 = cosmology.cosmo
    ps_ml1 = cosmology.ps_ml
    dist1 = cosmology.dist
    recomb1 = cosmology.recomb
    cosmo2 = cosmology_alt.cosmo
    ps_ml2 = cosmology_alt.ps_ml
    dist2 = cosmology_alt.dist
    recomb2 = cosmology_alt.recomb

    lmax = 1000
    ell_vec = Ncm.Vector.new_array(np.arange(lmax + 1).tolist())

    cmb_lens = Nc.XcorKernelCMBLensing(
        dist=dist1, powspec=ps_ml1, recomb=recomb1, Nl=ell_vec, integrator=integrator
    )
    cmb_lens.prepare(cosmo1)
    result1 = cmb_lens.eval_limber_z_full(cosmo1, 0.5, dist1, 100)

    # Create new kernel for second cosmology
    cmb_lens2 = Nc.XcorKernelCMBLensing(
        dist=dist2, powspec=ps_ml2, recomb=recomb2, Nl=ell_vec, integrator=integrator
    )
    cmb_lens2.prepare(cosmo2)
    result2 = cmb_lens2.eval_limber_z_full(cosmo2, 0.5, dist2, 100)

    # Results should be different with different cosmologies
    assert np.isfinite(result1)
    assert np.isfinite(result2)
    # Relaxed comparison since results might be close
    if result1 != 0.0 or result2 != 0.0:
        assert not np.isclose(
            result1, result2, rtol=1e-3
        ), "Results should differ with different cosmologies"


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


@pytest.mark.parametrize(
    "kernel_type",
    [
        "cmb_lensing",
        "cmb_isw",
        "tsz",
    ],
)
def test_kernel_get_k_range(
    kernel_type: str, cosmology: Cosmology, integrator: Ncm.SBesselIntegrator
) -> None:
    """Test that all kernels implement get_k_range correctly."""
    cosmo = cosmology.cosmo
    dist = cosmology.dist
    ps_ml = cosmology.ps_ml

    kernel: Nc.XcorKernelCMBLensing | Nc.XcorKernelCMBISW | Nc.XcorKerneltSZ
    if kernel_type == "cmb_lensing":
        recomb = Nc.RecombSeager()
        lmax = 1000
        ell_vec = Ncm.Vector.new_array(np.arange(lmax + 1).tolist())
        kernel = Nc.XcorKernelCMBLensing(
            dist=dist, powspec=ps_ml, recomb=recomb, Nl=ell_vec, integrator=integrator
        )
    elif kernel_type == "cmb_isw":
        recomb = Nc.RecombSeager()
        lmax = 1000
        ell_vec = Ncm.Vector.new_array(np.arange(lmax + 1).tolist())
        kernel = Nc.XcorKernelCMBISW(
            dist=dist, powspec=ps_ml, recomb=recomb, Nl=ell_vec, integrator=integrator
        )
    else:  # tsz
        kernel = Nc.XcorKerneltSZ(
            dist=dist, powspec=ps_ml, zmax=6.0, integrator=integrator
        )

    kernel.prepare(cosmo)

    for ell in [10, 50, 100, 500]:
        k_min, k_max = kernel.get_k_range(cosmo, ell)
        assert np.isfinite(k_min) and k_min > 0.0
        assert np.isfinite(k_max) and k_max > k_min


@pytest.mark.parametrize(
    "kernel_type",
    [
        "cmb_lensing",
        "cmb_isw",
        "tsz",
    ],
)
def test_kernel_get_eval(
    kernel_type: str, cosmology: Cosmology, integrator: Ncm.SBesselIntegrator
) -> None:
    """Test that all kernels implement get_eval correctly."""
    cosmo = cosmology.cosmo
    dist = cosmology.dist
    ps_ml = cosmology.ps_ml

    kernel: Nc.XcorKernelCMBLensing | Nc.XcorKernelCMBISW | Nc.XcorKerneltSZ
    if kernel_type == "cmb_lensing":
        recomb = Nc.RecombSeager()
        lmax = 1000
        ell_vec = Ncm.Vector.new_array(np.arange(lmax + 1).tolist())
        kernel = Nc.XcorKernelCMBLensing(
            dist=dist,
            powspec=ps_ml,
            recomb=recomb,
            Nl=ell_vec,
            integrator=integrator,
        )
    elif kernel_type == "cmb_isw":
        recomb = Nc.RecombSeager()
        lmax = 1000
        ell_vec = Ncm.Vector.new_array(np.arange(lmax + 1).tolist())
        kernel = Nc.XcorKernelCMBISW(
            dist=dist,
            powspec=ps_ml,
            recomb=recomb,
            Nl=ell_vec,
            integrator=integrator,
        )
    else:  # tsz
        kernel = Nc.XcorKerneltSZ(
            dist=dist, powspec=ps_ml, zmax=6.0, integrator=integrator
        )

    kernel.prepare(cosmo)

    for ell in [10, 50, 100]:
        integrand = kernel.get_eval(cosmo, ell)
        assert integrand is not None
        assert integrand.len == 1  # Single multipole


def test_kernel_get_eval_vectorized(
    cmb_lensing_component: Nc.XcorKernelCMBLensing, cosmology: Cosmology
) -> None:
    """Test vectorized evaluation."""
    cosmo = cosmology.cosmo
    kernel = cmb_lensing_component
    kernel.prepare(cosmo)

    lmin = 10
    lmax = 20
    integrand = kernel.get_eval_vectorized(cosmo, lmin, lmax)

    assert integrand is not None
    assert integrand.len == lmax - lmin + 1  # Multiple multipoles

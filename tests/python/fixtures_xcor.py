#!/usr/bin/env python
#
# fixtures_xcor.py
#
# Thu Aug 01 11:44:12 2024
# Copyright  2024  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# fixtures_xcor.py
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

"""Fixtures for XCor tests."""

import itertools as it
from typing import Callable
import functools
import pytest
from pytest_lazy_fixtures import lf

import numpy as np
import numcosmo_py.cosmology as ncpy
from numcosmo_py import Ncm, Nc, GObject

Ncm.cfg_init()


# Define a simple test component that subclasses NcXcorKernelComponent in Python
class NcXcorKernelComponentTest(Nc.XcorKernelComponent):
    """
    A simple test component that implements a Gaussian kernel in xi-space.
    K(k, xi) = exp(-xi²/2sigma²) * sin(k*xi)/(k*xi)
    This is a simple, well-behaved kernel for testing purposes.
    """

    __gtype_name__ = "NcXcorKernelComponentTest"

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


@functools.cache
def _get_cosmology() -> ncpy.Cosmology:
    """Create a default cosmology for testing."""
    return ncpy.Cosmology.default()


@functools.cache
def _get_cosmology_alt() -> ncpy.Cosmology:
    """Create a default cosmology for testing."""
    cosmology = ncpy.Cosmology.default()
    cosmology.cosmo["H0"] = 75.0
    cosmology.cosmo["Omegab"] = 0.03
    cosmology.cosmo["Omegac"] = 0.22

    cosmology.prepare()

    return cosmology


@pytest.fixture(name="cosmology", scope="module")
def fixture_cosmology() -> ncpy.Cosmology:
    """Create a simple cosmology for testing."""
    return _get_cosmology()


@pytest.fixture(name="cosmology_alt", scope="module")
def fixture_cosmology_alt() -> ncpy.Cosmology:
    """Create a simple cosmology alternate for testing."""
    cosmology = _get_cosmology_alt()
    return cosmology


@functools.cache
def _lsst_y1_lens_bins() -> (
    tuple[list[Nc.GalaxySDObsRedshiftGauss], Nc.GalaxySDTrueRedshiftLSSTSRD]
):
    """Create LSST Y1 lens photo-z bins."""
    bins, gsdtr = Nc.GalaxySDObsRedshiftGauss.new_lsst_srd_bins(
        Nc.GalaxySDTrueRedshiftLSSTSRDType.Y1_LENS
    )
    return bins, gsdtr


@functools.cache
def _lsst_y1_source_bins() -> (
    tuple[list[Nc.GalaxySDObsRedshiftGauss], Nc.GalaxySDTrueRedshiftLSSTSRD]
):
    """Create LSST Y1 source photo-z bins."""
    bins, gsdtr = Nc.GalaxySDObsRedshiftGauss.new_lsst_srd_bins(
        Nc.GalaxySDTrueRedshiftLSSTSRDType.Y1_SOURCE
    )
    return bins, gsdtr


@functools.cache
def _lsst_y1_lens_bins_spec() -> list[tuple[str, int]]:
    """Generate specs for LSST Y1 lens bins."""
    bins, _ = _lsst_y1_lens_bins()
    return [(f"lsst_y1_lens_bin{i}", i) for i in range(len(bins))]


@functools.cache
def _lsst_y1_source_bins_spec() -> list[tuple[str, int]]:
    """Generate specs for LSST Y1 source bins."""
    bins, _ = _lsst_y1_source_bins()
    return [(f"lsst_y1_source_bin{i}", i) for i in range(len(bins))]


@pytest.fixture(
    name="lsst_y1_lens_bin", params=_lsst_y1_lens_bins_spec(), scope="module"
)
def fixture_lsst_y1_lens_bin(
    request: pytest.FixtureRequest,
) -> Nc.GalaxySDObsRedshiftGauss:
    """Fixture for LSST Y1 lens photo-z bins."""
    bins, _ = _lsst_y1_lens_bins()
    return bins[request.param[1]]


@pytest.fixture(
    name="lsst_y1_source_bin", params=_lsst_y1_source_bins_spec(), scope="module"
)
def fixture_lsst_y1_source_bin(
    request: pytest.FixtureRequest,
) -> Nc.GalaxySDObsRedshiftGauss:
    """Fixture for LSST Y1 source photo-z bins."""
    bins, _ = _lsst_y1_source_bins()
    return bins[request.param[1]]


@pytest.fixture(name="lsst_y1_lens_bins", scope="module")
def fixture_lsst_y1_lens_bins() -> tuple[list, Nc.GalaxySDTrueRedshiftLSSTSRD]:
    """Create LSST Y1 lens photo-z bins."""
    bins_array, gsdtr = _lsst_y1_lens_bins()
    return bins_array, gsdtr


@pytest.fixture(name="lsst_y1_source_bins", scope="module")
def fixture_lsst_y1_source_bins() -> tuple[list, Nc.GalaxySDTrueRedshiftLSSTSRD]:
    """Create LSST Y1 source photo-z bins."""
    bins_array, gsdtr = _lsst_y1_source_bins()
    return bins_array, gsdtr


@functools.cache
def _get_integrator() -> Ncm.SBesselIntegrator:
    """Create and cache a spherical Bessel integrator."""
    integrator = Ncm.SBesselIntegratorLevin(
        ell_min=0,
        ell_max=8,
        y_knots_min=1.0,
        y_knots_max=1.0e4,
        n_knots=10,
        ell_cache_max=1200,
        reltol=1.0e-13,
    )
    return integrator


@pytest.fixture(name="integrator", scope="module")
def fixture_integrator() -> Ncm.SBesselIntegrator:
    """Create a simple spherical Bessel integrator."""
    return _get_integrator()


@functools.cache
def _get_kernel_cmb_lens() -> Nc.XcorKernel:
    """Create and cache a CMB lensing kernel."""
    cosmology = _get_cosmology()
    integrator = _get_integrator()
    return Nc.XcorKernelCMBLensing(
        dist=cosmology.dist,
        powspec=cosmology.ps_ml,
        recomb=cosmology.recomb,
        Nl=Ncm.Vector.new_array(np.arange(3000 + 1)),
        integrator=integrator,
    )


@pytest.fixture(name="kernel_cmb_lens")
def fixture_kernel_cmb_lens() -> Nc.XcorKernel:
    """Fixture for NumCosmo CMB lensing tracer."""
    return _get_kernel_cmb_lens()


@functools.cache
def _get_kernel_cmb_isw() -> Nc.XcorKernel:
    """Create and cache a CMB ISW kernel."""
    cosmology = _get_cosmology()
    integrator = _get_integrator()
    return Nc.XcorKernelCMBISW(
        dist=cosmology.dist,
        powspec=cosmology.ps_ml,
        recomb=cosmology.recomb,
        Nl=Ncm.Vector.new_array(np.arange(3000 + 1)),
        integrator=integrator,
    )


@pytest.fixture(name="kernel_cmb_isw")
def fixture_kernel_cmb_isw() -> Nc.XcorKernel:
    """Fixture for NumCosmo CMB ISW tracer."""
    return _get_kernel_cmb_isw()


def _get_kernel_tsz() -> Nc.XcorKerneltSZ:
    """Create and cache a tSZ kernel."""
    cosmology = _get_cosmology()
    integrator = _get_integrator()
    return Nc.XcorKerneltSZ(
        dist=cosmology.dist,
        powspec=cosmology.ps_ml,
        zmax=6.0,
        integrator=integrator,
    )


@pytest.fixture(name="kernel_tsz")
def fixture_kernel_tsz() -> Nc.XcorKernel:
    """Fixture for NumCosmo tSZ tracer."""
    return _get_kernel_tsz()


def _get_kernel_gal(
    bin_idx: int,
    bias: float = 1.5,
    bias_nknots: int = 1,
    noise_bias: float = 1.234,
    domagbias: bool = True,
) -> Nc.XcorKernel:
    """Create XcorKernelGal from LSST photo-z bin."""
    bins, _ = _lsst_y1_lens_bins()
    dndz_spline = bins[bin_idx].compute_binned_dndz(None)
    cosmology = _get_cosmology()
    integrator = _get_integrator()
    kernel_gal = Nc.XcorKernelGal(
        dist=cosmology.dist,
        powspec=cosmology.ps_ml,
        bparam_length=bias_nknots,
        noise_bias=noise_bias,
        dndz=dndz_spline,
        domagbias=domagbias,
        integrator=integrator,
    )
    for i in range(bias_nknots):
        kernel_gal.orig_vparam_set(Nc.XcorKernelGalVParams.BIAS, i, bias)

    return kernel_gal


def _get_kernel_wl(
    bin_idx: int, nbar: float = 3.0, intr_shear: float = 7.0
) -> Nc.XcorKernel:
    """Create XcorKernelWeakLensing from LSST photo-z bin."""
    bins, _ = _lsst_y1_source_bins()
    dndz_spline = bins[bin_idx].compute_binned_dndz(None)
    cosmology = _get_cosmology()
    integrator = _get_integrator()
    nc_wl = Nc.XcorKernelWeakLensing(
        dist=cosmology.dist,
        powspec=cosmology.ps_ml,
        dndz=dndz_spline,
        nbar=nbar,
        intr_shear=intr_shear,
        integrator=integrator,
    )
    return nc_wl


@pytest.fixture(name="kernel_gal", params=_lsst_y1_lens_bins_spec(), ids=lambda x: x[1])
def fixture_kernel_gal(request: pytest.FixtureRequest) -> Nc.XcorKernel:
    """Fixture for NumCosmo galaxy tracer using LSST Y1 lens bins."""
    return _get_kernel_gal(request.param[1])


@pytest.fixture(
    name="kernel_wl", params=_lsst_y1_source_bins_spec(), ids=lambda x: x[1]
)
def fixture_kernel_wl(request: pytest.FixtureRequest) -> Nc.XcorKernel:
    """Fixture for NumCosmo weak lensing tracer using LSST Y1 source bins."""
    return _get_kernel_wl(request.param[1])


# Kernel fixture metadata for dynamic test generation
KERNEL_FIXTURES: list[tuple[str, int, Callable] | tuple[str, None, Callable]] = (
    [
        ("kernel_cmb_lens", None, _get_kernel_cmb_lens),
        ("kernel_cmb_isw", None, _get_kernel_cmb_isw),
        ("kernel_tsz", None, _get_kernel_tsz),
    ]
    + [
        (f"kernel_gal_bin{i}", i, lambda ii=i: _get_kernel_gal(ii))
        for _, i in _lsst_y1_lens_bins_spec()
    ]
    + [
        (f"kernel_wl_bin{i}", i, lambda ii=i: _get_kernel_wl(ii))
        for _, i in _lsst_y1_source_bins_spec()
    ]
)


def make_kernel_gal_fixture(kernel_name: str, bin_idx: int | None):
    """Create a fixture for a kernel."""
    assert bin_idx is not None

    @pytest.fixture(name=kernel_name)
    def _fixture() -> Nc.XcorKernel:
        return _get_kernel_gal(bin_idx)

    return _fixture


def make_kernel_wl_fixture(kernel_name: str, bin_idx: int | None):
    """Create a fixture for a kernel."""
    assert bin_idx is not None

    @pytest.fixture(name=kernel_name)
    def _fixture() -> Nc.XcorKernel:
        return _get_kernel_wl(bin_idx)

    return _fixture


for name, idx, _ in KERNEL_FIXTURES:
    if name.startswith("kernel_gal_bin"):
        globals()[name] = make_kernel_gal_fixture(name, idx)
    elif name.startswith("kernel_wl_bin"):
        globals()[name] = make_kernel_wl_fixture(name, idx)


@functools.cache
def _get_kernel_components() -> (
    list[tuple[str, Nc.XcorKernel | None, Nc.XcorKernelComponent]]
):
    """Create and cache a list of kernel components."""
    components: list[tuple[str, Nc.XcorKernel | None, Nc.XcorKernelComponent]] = [
        ("kernel_test_component", None, NcXcorKernelComponentTest())
    ]
    for kernel_name, _, kernel_func in KERNEL_FIXTURES:
        kernel = kernel_func()
        assert isinstance(kernel, Nc.XcorKernel)
        components += [
            (kernel_name, kernel, comp) for comp in kernel.get_component_list()
        ]
    return components


@functools.cache
def _get_kernel_components_spec() -> list[tuple[str, int]]:
    """Generate specs for kernel components."""
    components = _get_kernel_components()
    return [
        (
            (
                f"{kernel_name}_"
                f"{comp.__class__.__name__.replace('NcXcorKernelComponent', '')}"
            ),
            idx,
        )
        for idx, (kernel_name, _, comp) in enumerate(components)
    ]


@pytest.fixture(
    name="kernel_component", params=_get_kernel_components_spec(), ids=lambda x: x[0]
)
def fixture_kernel_component(
    request: pytest.FixtureRequest, cosmology: ncpy.Cosmology
) -> tuple[str, Nc.XcorKernel | None, Nc.XcorKernelComponent]:
    """Fixture for kernel components."""
    components = _get_kernel_components()
    comp_idx = request.param[1]
    kernel_name: str
    kernel: Nc.XcorKernel | None
    component: Nc.XcorKernelComponent
    kernel_name, kernel, component = components[comp_idx]
    if kernel is not None:
        kernel.prepare(cosmology.cosmo)
    return kernel_name, kernel, component


def pytest_generate_tests(metafunc: pytest.Metafunc) -> None:
    """
    Dynamically generate test cases based on test function signature.

    - Tests with 'k1' and 'k2' parameters: iterate over kernel pairs
    """

    if "k1" in metafunc.fixturenames and "k2" in metafunc.fixturenames:
        pairs: list[tuple] = []
        ids = []

        # Generate all combinations (can be customized)
        for (k1_name, _k1_idx, _k1_func), (
            k2_name,
            _k2_idx,
            _k2_func,
        ) in it.combinations_with_replacement(KERNEL_FIXTURES, 2):
            pairs.append((lf(k1_name), lf(k2_name)))
            ids.append(f"{k1_name}-{k2_name}")

        metafunc.parametrize("k1,k2", pairs, ids=ids)


@pytest.fixture(
    name="kernel_case",
    params=[(kernel_name, lf(kernel_name)) for kernel_name, _, _ in KERNEL_FIXTURES],
    ids=lambda x: x[0],
)
def fixture_kernel_case(request: pytest.FixtureRequest) -> tuple[str, Nc.XcorKernel]:
    """Fixture for parametrized kernels."""
    kernel_name, kernel = request.param
    assert isinstance(kernel, Nc.XcorKernel)
    return kernel_name, kernel


@pytest.fixture(name="kernel")
def fixture_kernel(kernel_case: tuple[str, Nc.XcorKernel]) -> Nc.XcorKernel:
    """Fixture for kernels."""
    _, kernel = kernel_case
    assert isinstance(kernel, Nc.XcorKernel)
    return kernel

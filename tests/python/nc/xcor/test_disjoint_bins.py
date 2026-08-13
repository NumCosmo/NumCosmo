#!/usr/bin/env python
#
# test_disjoint_bins.py
#
# Mon Aug 11 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_disjoint_bins.py
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

"""Non-Limber cross spectra of kernels with disjoint redshift support.

Two tomographic bins that do not overlap in redshift still correlate: they are
windows on the same 3D field, and the non-Limber kernel-space methods couple
them through the outer k integral rather than through a shared z integral. The
Limber-z tier is the only one for which a vanishing z overlap implies a
vanishing C_l, and applying that short-circuit to the kernel-space tiers
silently zeroed every disjoint-bin cross spectrum -- which is exactly the
super-sample-covariance use case (see numcosmo_py/ssc.py).

Also pins the radial-window normalization of NcXcorKernelClusterTophat: the
window is normalized in comoving volume alone, so that C_0/4pi is the
super-sample covariance S_ij of the bin.
"""

import numpy as np
import pytest
from numpy.testing import assert_allclose

from numcosmo_py import Nc, Ncm
from numcosmo_py.cosmology import Cosmology

pytest_plugins = [
    "python.fixtures_xcor",
]

pytestmark = pytest.mark.xcor

# Two disjoint top-hat bins, and one overlapping the first, all narrow enough
# that the Limber-z overlap of the disjoint pair is empty.
Z_BINS = [(0.1, 0.2), (0.6, 0.7)]


def _make_kernels(cosmology: Cosmology, bins=None) -> list[Nc.XcorKernel]:
    """Build non-Limber top-hat cluster kernels for the given redshift bins."""
    kernels = []

    for z_lower, z_upper in bins if bins is not None else Z_BINS:
        kernel = Nc.XcorKernelClusterTophat(
            dist=cosmology.dist,
            powspec=cosmology.ps_ml,
            z_lower=z_lower,
            z_upper=z_upper,
            integrator=Ncm.SBesselIntegratorLevin.new(0, 8),
        )
        kernel.set_l_limber(-1)
        kernel.prepare(cosmology.cosmo)
        kernels.append(kernel)

    return kernels


@pytest.mark.parametrize(
    "method", [Nc.XcorMethod.KERNEL_GSL, Nc.XcorMethod.KERNEL_CUBATURE]
)
def test_disjoint_bins_cross_is_nonzero(
    cosmology: Cosmology, method: Nc.XcorMethod
) -> None:
    """Non-Limber cross spectra of disjoint bins must not be identically zero."""
    k1, k2 = _make_kernels(cosmology)

    xcor = Nc.Xcor.new(cosmology.dist, cosmology.ps_ml, method)
    # The disjoint-bin cross spectrum is a small residual of a large
    # cancellation. KERNEL_GSL targets reltol * 1e-2 in the outer k integral and
    # aborts on GSL roundoff at the default 1e-6; the result is already stable
    # to five digits at 1e-4.
    if method == Nc.XcorMethod.KERNEL_GSL:
        xcor.set_reltol(1.0e-4)

    xcor.prepare(cosmology.cosmo)

    cross = Ncm.Vector.new(3)
    xcor.compute(k1, k2, cosmology.cosmo, 0, 2, cross)
    cross_values = np.array(cross.dup_array())

    assert np.all(np.isfinite(cross_values))
    assert np.all(cross_values != 0.0), (
        "disjoint tomographic bins have a non-zero non-Limber cross spectrum; "
        "the Limber-z redshift-overlap short-circuit must not be applied here"
    )

    # Cross-correlation of disjoint shells is bounded by the geometric mean of
    # the autos (Cauchy-Schwarz), and is negative at ell = 0 for shells that do
    # not touch.
    auto1, auto2 = Ncm.Vector.new(3), Ncm.Vector.new(3)
    xcor.compute(k1, k1, cosmology.cosmo, 0, 2, auto1)
    xcor.compute(k2, k2, cosmology.cosmo, 0, 2, auto2)

    bound = np.sqrt(np.array(auto1.dup_array()) * np.array(auto2.dup_array()))
    assert np.all(np.abs(cross_values) < bound)
    assert cross_values[0] < 0.0


def test_disjoint_bins_cross_matches_solver(cosmology: Cosmology) -> None:
    """The solver reproduces nc_xcor_compute() for disjoint-bin cross spectra."""
    k1, k2 = _make_kernels(cosmology)

    xcor = Nc.Xcor.new(cosmology.dist, cosmology.ps_ml, Nc.XcorMethod.KERNEL_CUBATURE)
    xcor.prepare(cosmology.cosmo)

    solver = Nc.XcorSolver.new()
    id1 = solver.register_kernel(k1)
    id2 = solver.register_kernel(k2)
    solver.request_cl(id1, id2, 0, 3)
    solver.plan_blocks(8)
    solver.solve(xcor, cosmology.cosmo)

    expected = Ncm.Vector.new(4)
    xcor.compute(k1, k2, cosmology.cosmo, 0, 3, expected)

    assert_allclose(
        np.array(solver.get_result(0).dup_array()),
        np.array(expected.dup_array()),
        rtol=1.0e-6,
        atol=1.0e-50,
    )


@pytest.mark.parametrize(
    "method",
    [
        Nc.XcorMethod.KERNEL_GSL,
        Nc.XcorMethod.KERNEL_CUBATURE,
        Nc.XcorMethod.KERNEL_FIXED,
    ],
)
def test_kernel_space_limber_disjoint_is_zero(
    cosmology: Cosmology, method: Nc.XcorMethod
) -> None:
    """Kernel-space methods still vanish for disjoint bins in the Limber tier.

    A Limber kernel is supported only where xi = (l + 1/2) / k lies inside its
    own radial range, so disjoint bins have disjoint support in k. Dropping the
    overlap short-circuit for every kernel-space method left this integrating
    the product of the two exponential extrapolation tails -- a numerical
    smoothing device, not physics -- which returned up to -48% of the auto
    spectrum instead of zero. The tier, not the method, decides.
    """
    kernels = []

    for z_lower, z_upper in Z_BINS:
        kernel = Nc.XcorKernelClusterTophat(
            dist=cosmology.dist,
            powspec=cosmology.ps_ml,
            z_lower=z_lower,
            z_upper=z_upper,
            integrator=Ncm.SBesselIntegratorLevin.new(0, 8),
        )
        kernel.set_l_limber(0)  # Limber everywhere, the class default
        kernel.prepare(cosmology.cosmo)
        kernels.append(kernel)

    xcor = Nc.Xcor.new(cosmology.dist, cosmology.ps_ml, method)
    xcor.prepare(cosmology.cosmo)

    cross = Ncm.Vector.new(3)
    xcor.compute(kernels[0], kernels[1], cosmology.cosmo, 10, 12, cross)

    assert_allclose(np.array(cross.dup_array()), 0.0, atol=0.0)

    # The auto spectrum of the same kernel is emphatically not zero, so this is
    # a real short-circuit and not an all-zero configuration.
    auto = Ncm.Vector.new(3)
    xcor.compute(kernels[0], kernels[0], cosmology.cosmo, 10, 12, auto)

    assert np.all(np.array(auto.dup_array()) > 0.0)


def test_limber_z_keeps_disjoint_short_circuit(cosmology: Cosmology) -> None:
    """The Limber-z tier still returns zero for disjoint bins, as it must."""
    k1, k2 = _make_kernels(cosmology)

    xcor = Nc.Xcor.new(cosmology.dist, cosmology.ps_ml, Nc.XcorMethod.LIMBER_Z_GSL)
    xcor.prepare(cosmology.cosmo)

    cross = Ncm.Vector.new(3)
    xcor.compute(k1, k2, cosmology.cosmo, 10, 12, cross)

    assert_allclose(np.array(cross.dup_array()), 0.0, atol=0.0)


def test_cluster_tophat_window_is_volume_normalized(cosmology: Cosmology) -> None:
    """The top-hat radial window integrates to one in comoving distance.

    The whole normalization is carried by 1/(V_upper - V_lower); an extra
    1/(z_upper - z_lower) would rescale every C_l by 1/(dz_i dz_j) and break
    the identification of C_0/4pi with the super-sample covariance.
    """
    cosmo = cosmology.cosmo
    dist = cosmology.dist
    dist.prepare(cosmo)

    z_lower, z_upper = 0.1, 0.2
    kernel = Nc.XcorKernelClusterTophat(
        dist=dist,
        powspec=cosmology.ps_ml,
        z_lower=z_lower,
        z_upper=z_upper,
        integrator=Ncm.SBesselIntegratorLevin.new(0, 8),
    )
    kernel.prepare(cosmo)

    # The Limber-z window is W(z) = xi_t^2, and its prefactor is 1/Delta V, so
    # int dz (dxi/dz) W(z) * prefactor must be exactly one.
    z_grid = np.linspace(z_lower, z_upper, 2001)
    xck = Nc.XcorKinetic()
    window = np.array([kernel.eval_limber_z(cosmo, float(z), xck, 0) for z in z_grid])
    dxi_dz = np.array([1.0 / cosmo.E(float(z)) for z in z_grid])
    prefactor = kernel.eval_limber_z_prefactor(cosmo, 0)

    assert_allclose(np.trapezoid(window * dxi_dz, z_grid) * prefactor, 1.0, rtol=1e-6)

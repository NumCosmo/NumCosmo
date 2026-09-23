#
# test_kernel_cmb_isw_source.py
#
# Thu Sep 11 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_kernel_cmb_isw_source.py
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

"""The source placement of the ISW kernel.

A photon that last scattered at chi' integrates the decay of the potential
over chi < chi' only, so with sources spread along the line of sight the ISW
kernel at chi is the thin-screen one times the fraction of photons that last
scatter beyond chi. Below the last-scattering shell that fraction is one and
the kernels coincide; inside the shell the visibility kernel goes to zero
smoothly where the thin screen stops; with reionization the fraction drops by
the rescattered tenth between the reionization bump and the shell.
"""

import numpy as np
import pytest
from numpy.testing import assert_allclose

from numcosmo_py import Nc, Ncm
from numcosmo_py.cosmology import Cosmology
from numcosmo_py.helper import duplicate_via_serialization

pytestmark = [pytest.mark.xcor, pytest.mark.xdist_group("cmb_isw_source")]

Ncm.cfg_init()

SOURCES = (
    Nc.XcorKernelCMBISWSource.THIN_SCREEN,
    Nc.XcorKernelCMBISWSource.VISIBILITY,
    Nc.XcorKernelCMBISWSource.VISIBILITY_REIONIZATION,
)


@pytest.fixture(name="cosmology", scope="module")
def fixture_cosmology() -> Cosmology:
    """Return the default cosmology with distances past recombination."""
    return Cosmology.default(dist_max_z=1500.0)


def _kernel(
    cosmology: Cosmology, source: Nc.XcorKernelCMBISWSource
) -> Nc.XcorKernelCMBISW:
    kernel = Nc.XcorKernelCMBISW(
        dist=cosmology.dist,
        powspec=cosmology.ps_ml,
        recomb=cosmology.recomb,
        Nl=Ncm.Vector.new_array([0.0] * 11),
        source=source,
        integrator=Ncm.SBesselIntegratorLevin.new(2, 9),
    )
    kernel.set_lmax(10)
    kernel.set_l_limber(-1)
    kernel.prepare(cosmology.cosmo)
    return kernel


@pytest.fixture(name="kernels", scope="module")
def fixture_kernels(cosmology: Cosmology) -> dict:
    """Return one prepared kernel per source."""
    return {source: _kernel(cosmology, source) for source in SOURCES}


def _eval(kernel: Nc.XcorKernelCMBISW, cosmo: Nc.HICosmo, z: float, k: float):
    comp = kernel.get_component_list()[0]
    chi = kernel.peek_dist().comoving(cosmo, z)
    return comp.eval_kernel(cosmo, chi, k)


def test_source_property_round_trips(cosmology: Cosmology) -> None:
    """The source is a property, settable and serialized."""
    for source in SOURCES:
        kernel = _kernel(cosmology, source)
        assert kernel.get_source() == source
        assert kernel.props.source == source
        copy = duplicate_via_serialization(kernel)
        assert copy.get_source() == source


def test_set_source_marks_the_kernel_outdated(cosmology: Cosmology) -> None:
    """Changing the source changes the support after the next prepare."""
    kernel = _kernel(cosmology, Nc.XcorKernelCMBISWSource.THIN_SCREEN)
    cosmo = cosmology.cosmo
    _, chi_max_thin, _, _ = kernel.get_component_list()[0].get_limits(cosmo)

    kernel.set_source(Nc.XcorKernelCMBISWSource.VISIBILITY)
    kernel.prepare(cosmo)
    _, chi_max_vis, _, _ = kernel.get_component_list()[0].get_limits(cosmo)

    assert kernel.get_source() == Nc.XcorKernelCMBISWSource.VISIBILITY
    assert chi_max_vis > chi_max_thin


def test_visibility_kernel_equals_thin_screen_below_the_shell(
    cosmology: Cosmology, kernels: dict
) -> None:
    """Below the shell every photon is still ahead: the survival fraction is one."""
    cosmo = cosmology.cosmo
    k = 0.05 * cosmo.RH_Mpc()
    thin = kernels[Nc.XcorKernelCMBISWSource.THIN_SCREEN]
    vis = kernels[Nc.XcorKernelCMBISWSource.VISIBILITY]
    for z in (0.1, 1.0, 10.0):
        assert_allclose(_eval(vis, cosmo, z, k), _eval(thin, cosmo, z, k), rtol=1.0e-12)
    # Residual ionization between the reionization bump and the shell scatters a
    # few tenths of a percent of the photons before z = 300.
    assert_allclose(
        _eval(vis, cosmo, 300.0, k), _eval(thin, cosmo, 300.0, k), rtol=1.0e-2
    )


def test_visibility_kernel_decays_across_the_shell(
    cosmology: Cosmology, kernels: dict
) -> None:
    """Inside the shell the kernel follows the cumulative visibility down to zero."""
    cosmo = cosmology.cosmo
    k = 0.05 * cosmo.RH_Mpc()
    thin = kernels[Nc.XcorKernelCMBISWSource.THIN_SCREEN]
    vis = kernels[Nc.XcorKernelCMBISWSource.VISIBILITY]
    _, z_max, _ = vis.get_z_range()
    z_lss = cosmology.dist.decoupling_redshift(cosmo)
    assert z_max > z_lss

    # The shell found by the visibility features starts near z = 720. The
    # residual ionization between the visibility minimum (z ~ 16) and the shell
    # scatters a few percent of the photons before z = 700.
    z_grid = np.linspace(700.0, z_max, 24)
    ratio = np.array(
        [_eval(vis, cosmo, z, k) / _eval(thin, cosmo, z, k) for z in z_grid]
    )
    assert np.all(np.diff(ratio) <= 1.0e-12)  # monotonically decreasing
    assert 0.95 < ratio[0] < 1.0
    assert 0.2 < ratio[np.argmin(np.abs(z_grid - z_lss))] < 0.8
    assert ratio[-1] < 1.0e-3


def test_reionization_lowers_the_kernel_by_the_rescattered_fraction(
    cosmology: Cosmology, kernels: dict
) -> None:
    """Between the reionization bump and the shell a tenth of the photons are gone."""
    cosmo = cosmology.cosmo
    recomb = cosmology.recomb
    k = 0.05 * cosmo.RH_Mpc()
    vis = kernels[Nc.XcorKernelCMBISWSource.VISIBILITY]
    reion = kernels[Nc.XcorKernelCMBISWSource.VISIBILITY_REIONIZATION]

    z_a = np.geomspace(1.0e-3, 50.0, 4000)
    g_a = np.array([recomb.v_tau(cosmo, -np.log1p(z)) / (1.0 + z) for z in z_a])
    frac = abs(np.trapezoid(g_a, z_a))
    assert 0.05 < frac < 0.15

    # Nothing has scattered yet at z = 0.1; by z = 50 the reionization bump is behind.
    assert_allclose(_eval(reion, cosmo, 0.1, k), _eval(vis, cosmo, 0.1, k), rtol=1.0e-3)
    ratio = _eval(reion, cosmo, 50.0, k) / _eval(vis, cosmo, 50.0, k)
    assert_allclose(ratio, 1.0 - frac, rtol=2.0e-2)


def test_spectra_are_finite_and_ordered(cosmology: Cosmology, kernels: dict) -> None:
    """Spreading the sources along the visibility lowers the ISW auto spectrum.

    The thin screen lets every photon collect the early ISW down to z_lss. Along
    the visibility, the photons that rescatter off residual electrons and the far
    half of the shell contribute less, so the auto spectrum, dominated here by the
    early ISW at recombination, drops by several percent; reionization removes
    another tenth of the photons before the shell.
    """
    cosmo = cosmology.cosmo
    xcor = Nc.Xcor.new(cosmology.dist, cosmology.ps_ml, Nc.XcorMethod.KERNEL_EXACT)
    xcor.set_closure_type(Nc.XcorKernelClosure.CHEBYSHEV)
    xcor.set_ell_batch_size(8)
    xcor.prepare(cosmo)
    cls = {}
    for source, kernel in kernels.items():
        vp = Ncm.Vector.new(8)
        xcor.compute(kernel, None, cosmo, 2, 9, vp)
        cls[source] = np.array(vp.dup_array())
        assert np.all(np.isfinite(cls[source]))
        assert np.all(cls[source] > 0.0)

    thin = cls[Nc.XcorKernelCMBISWSource.THIN_SCREEN]
    vis = cls[Nc.XcorKernelCMBISWSource.VISIBILITY]
    reion = cls[Nc.XcorKernelCMBISWSource.VISIBILITY_REIONIZATION]
    assert np.all(vis < thin)
    assert np.all(vis / thin > 0.85)
    assert np.all(reion <= vis)
    assert np.all(reion / vis > 0.7)

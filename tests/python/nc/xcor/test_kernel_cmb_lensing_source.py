#
# test_kernel_cmb_lensing_source.py
#
# Thu Sep 11 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_kernel_cmb_lensing_source.py
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

"""The source placement of the CMB lensing kernel.

NcXcorKernelCMBLensing places the CMB photons either on a thin screen at the
decoupling redshift or along the recombination visibility function of its
NcRecomb, with or without the reionization bump. The visibility sources are
integrated by NcXcorLensingEfficiency; at redshifts far below the shell they
must agree with the thin screen to the width of the shell over chi_*, and the
reionization variant must lower the kernel, between the reionization bump and
the shell, by the fraction of photons that last scatter at low redshift.
"""

import numpy as np
import pytest
from numpy.testing import assert_allclose

from numcosmo_py import Nc, Ncm
from numcosmo_py.cosmology import Cosmology
from numcosmo_py.helper import duplicate_via_serialization

pytestmark = [pytest.mark.xcor, pytest.mark.xdist_group("cmb_lensing_source")]

Ncm.cfg_init()

SOURCES = (
    Nc.XcorKernelCMBLensingSource.THIN_SCREEN,
    Nc.XcorKernelCMBLensingSource.VISIBILITY,
    Nc.XcorKernelCMBLensingSource.VISIBILITY_REIONIZATION,
)


@pytest.fixture(name="cosmology", scope="module")
def fixture_cosmology() -> Cosmology:
    """Return the default cosmology with distances past recombination."""
    return Cosmology.default(dist_max_z=1500.0)


def _kernel(
    cosmology: Cosmology, source: Nc.XcorKernelCMBLensingSource
) -> Nc.XcorKernelCMBLensing:
    kernel = Nc.XcorKernelCMBLensing(
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


def _eval(kernel: Nc.XcorKernelCMBLensing, cosmo: Nc.HICosmo, z: float, k: float):
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
    kernel = _kernel(cosmology, Nc.XcorKernelCMBLensingSource.THIN_SCREEN)
    cosmo = cosmology.cosmo
    _, chi_max_thin, _, _ = kernel.get_component_list()[0].get_limits(cosmo)

    kernel.set_source(Nc.XcorKernelCMBLensingSource.VISIBILITY)
    kernel.prepare(cosmo)
    _, chi_max_vis, _, _ = kernel.get_component_list()[0].get_limits(cosmo)

    assert kernel.get_source() == Nc.XcorKernelCMBLensingSource.VISIBILITY
    assert chi_max_vis > chi_max_thin


def test_visibility_support_reaches_past_decoupling(
    cosmology: Cosmology, kernels: dict
) -> None:
    """The visibility source extends the support to the far edge of the shell."""
    cosmo = cosmology.cosmo
    thin = kernels[Nc.XcorKernelCMBLensingSource.THIN_SCREEN]
    vis = kernels[Nc.XcorKernelCMBLensingSource.VISIBILITY]
    _, chi_thin, _, _ = thin.get_component_list()[0].get_limits(cosmo)
    _, chi_vis, _, _ = vis.get_component_list()[0].get_limits(cosmo)
    _, z_thin, _ = thin.get_z_range()
    _, z_vis, _ = vis.get_z_range()

    assert_allclose(z_thin, cosmology.dist.decoupling_redshift(cosmo), rtol=1.0e-12)
    assert z_vis > z_thin
    assert chi_vis > chi_thin
    # The shell is about a percent of chi_* wide.
    assert chi_vis / chi_thin - 1.0 < 0.02


def test_visibility_agrees_with_thin_screen_at_low_redshift(
    cosmology: Cosmology, kernels: dict
) -> None:
    """Far from the shell the two sources give the same kernel."""
    cosmo = cosmology.cosmo
    k = 0.05 * cosmo.RH_Mpc()
    thin = kernels[Nc.XcorKernelCMBLensingSource.THIN_SCREEN]
    vis = kernels[Nc.XcorKernelCMBLensingSource.VISIBILITY]

    for z, tol in ((0.1, 1.0e-4), (1.0, 1.0e-3), (5.0, 5.0e-3)):
        assert_allclose(_eval(vis, cosmo, z, k), _eval(thin, cosmo, z, k), rtol=tol)

    # Inside the shell the thin screen has already dropped to zero.
    z_in = 0.5 * (cosmology.dist.decoupling_redshift(cosmo) + vis.get_z_range()[1])
    assert _eval(vis, cosmo, z_in, k) > 0.0


def test_reionization_lowers_the_kernel_by_the_rescattered_fraction(
    cosmology: Cosmology, kernels: dict
) -> None:
    """Photons rescattered at reionization are lensed by less structure."""
    cosmo = cosmology.cosmo
    recomb = cosmology.recomb
    k = 0.05 * cosmo.RH_Mpc()
    vis = kernels[Nc.XcorKernelCMBLensingSource.VISIBILITY]
    reion = kernels[Nc.XcorKernelCMBLensingSource.VISIBILITY_REIONIZATION]

    # The fraction of photons that last scatter before z = 50.
    z_a = np.geomspace(1.0e-3, 50.0, 4000)
    g_a = np.array([recomb.v_tau(cosmo, -np.log1p(z)) / (1.0 + z) for z in z_a])
    frac = abs(np.trapezoid(g_a, z_a))
    assert 0.05 < frac < 0.15

    # At low z every source is far away and the efficiency does not care where:
    # the two kernels agree. Between the reionization bump and the shell only the
    # photons from the shell contribute, so the kernel is lowered by the
    # rescattered fraction.
    assert_allclose(_eval(reion, cosmo, 0.1, k), _eval(vis, cosmo, 0.1, k), rtol=5.0e-3)
    ratio = _eval(reion, cosmo, 50.0, k) / _eval(vis, cosmo, 50.0, k)
    assert_allclose(ratio, 1.0 - frac, rtol=2.0e-2)


def test_spectra_agree_between_thin_screen_and_visibility(
    cosmology: Cosmology, kernels: dict
) -> None:
    """The CMB lensing auto spectrum barely sees the width of the shell."""
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

    thin = cls[Nc.XcorKernelCMBLensingSource.THIN_SCREEN]
    vis = cls[Nc.XcorKernelCMBLensingSource.VISIBILITY]
    reion = cls[Nc.XcorKernelCMBLensingSource.VISIBILITY_REIONIZATION]
    assert_allclose(vis, thin, rtol=5.0e-3)
    assert np.all(reion < vis)
    assert np.all(reion / vis > 0.7)

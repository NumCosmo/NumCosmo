#
# test_kernel_gal_rsd_powspec.py
#
# Thu Sep 11 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_kernel_gal_rsd_powspec.py
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

"""The RSD term reads its growth rate from the power spectrum.

The linear growth rate f(k, z) of the redshift-space distortion component is
-(1 + z) dP/dz / (2 P), taken from the NcmPowspec the kernel was built with,
so a spectrum with scale-dependent growth carries it into the term and any
NcmPowspec, including a tabulated one, can drive it. These tests hold the
spectrum fixed and swap its representation: a tabulated copy of the
Eisenstein-Hu spectrum must give the same RSD spectra as the analytic one.
"""

import numpy as np
import pytest
from numpy.testing import assert_allclose

from numcosmo_py import Nc, Ncm
from numcosmo_py.cosmology import Cosmology

pytestmark = [pytest.mark.xcor, pytest.mark.xdist_group("rsd_powspec")]

Ncm.cfg_init()

LMIN, LMAX = 2, 12


@pytest.fixture(name="cosmology", scope="module")
def fixture_cosmology() -> Cosmology:
    """Return the default cosmology with the analytic linear spectrum."""
    return Cosmology.default(dist_max_z=20.0)


@pytest.fixture(name="ps_tabulated", scope="module")
def fixture_ps_tabulated(cosmology: Cosmology) -> Ncm.PowspecSpline2d:
    """Return the same linear spectrum tabulated in ln P on (z, ln k)."""
    cosmo, ps_ml = cosmology.cosmo, cosmology.ps_ml
    ps_ml.prepare(cosmo)
    za = np.linspace(0.0, 4.0, 81)
    ka = np.geomspace(ps_ml.get_kmin(), ps_ml.get_kmax(), 600)
    lnP = np.log(np.array([[ps_ml.eval(cosmo, z, k) for k in ka] for z in za]))
    sp2d = Ncm.Spline2dBicubic(
        spline=Ncm.SplineCubicNotaknot.new(),
        x_vector=Ncm.Vector.new_array(za.tolist()),
        y_vector=Ncm.Vector.new_array(np.log(ka).tolist()),
        z_matrix=Ncm.Matrix.new_array(lnP.T.flatten().tolist(), len(za)),
    )
    ps = Ncm.PowspecSpline2d.new(sp2d)
    ps.prepare(cosmo)
    return ps


def _dndz() -> Ncm.Spline:
    z_a = np.linspace(0.0, 1.5, 200)
    nz_a = np.exp(-0.5 * ((z_a - 0.5) / 0.1) ** 2)
    return Ncm.Spline.new_array(
        Ncm.SplineCubicNotaknot.new(), z_a.tolist(), nz_a.tolist(), True
    )


def _rsd_cl(cosmology: Cosmology, ps: Ncm.Powspec, *, dorsd: bool) -> np.ndarray:
    kernel = Nc.XcorKernelGal(
        dist=cosmology.dist,
        powspec=ps,
        bparam_length=1,
        dndz=_dndz(),
        domagbias=False,
        dorsd=dorsd,
        integrator=Ncm.SBesselIntegratorLevin.new(LMIN, LMAX),
    )
    kernel.orig_vparam_set(Nc.XcorKernelGalVParams.BIAS, 0, 1.0)
    kernel.set_l_limber(-1)
    kernel.prepare(cosmology.cosmo)
    xcor = Nc.Xcor.new(cosmology.dist, ps, Nc.XcorMethod.KERNEL_EXACT)
    xcor.prepare(cosmology.cosmo)
    res = Ncm.Vector.new(LMAX - LMIN + 1)
    xcor.compute(kernel, None, cosmology.cosmo, LMIN, LMAX, res)
    return np.array(res.dup_array())


def test_rsd_growth_rate_from_powspec_derivative(
    cosmology: Cosmology, ps_tabulated: Ncm.PowspecSpline2d
) -> None:
    """f(k, z) read from a tabulated spectrum reproduces the analytic one."""
    cosmo, ps_ml = cosmology.cosmo, cosmology.ps_ml
    for z, k in ((0.1, 0.02), (0.5, 0.1), (1.2, 0.5)):
        f_ml = -(1.0 + z) * ps_ml.deriv_z(cosmo, z, k) / (2.0 * ps_ml.eval(cosmo, z, k))
        f_tab = (
            -(1.0 + z)
            * ps_tabulated.deriv_z(cosmo, z, k)
            / (2.0 * ps_tabulated.eval(cosmo, z, k))
        )
        assert 0.4 < f_ml < 1.0
        assert_allclose(f_tab, f_ml, rtol=1.0e-4)


def test_rsd_spectrum_independent_of_powspec_representation(
    cosmology: Cosmology, ps_tabulated: Ncm.PowspecSpline2d
) -> None:
    """The RSD spectrum from the tabulated spectrum matches the analytic one."""
    cl_ml = _rsd_cl(cosmology, cosmology.ps_ml, dorsd=True)
    cl_tab = _rsd_cl(cosmology, ps_tabulated, dorsd=True)
    assert np.all(np.isfinite(cl_tab))
    assert_allclose(cl_tab, cl_ml, rtol=2.0e-4)


def test_rsd_term_changes_the_spectrum(cosmology: Cosmology) -> None:
    """The Kaiser term is not a no-op at low multipoles."""
    cl_rsd = _rsd_cl(cosmology, cosmology.ps_ml, dorsd=True)
    cl_no = _rsd_cl(cosmology, cosmology.ps_ml, dorsd=False)
    assert np.all(np.abs(cl_rsd / cl_no - 1.0) > 1.0e-3)


def test_rsd_runs_on_a_spectrum_without_closed_form_derivatives(
    cosmology: Cosmology,
) -> None:
    """A halofit spectrum has no deriv_z of its own; the default must serve."""
    ps_nl = Nc.PowspecMNLHaloFit.new(cosmology.ps_ml, 4.0, 1.0e-6)
    ps_nl.prepare(cosmology.cosmo)
    cl = _rsd_cl(cosmology, ps_nl, dorsd=True)
    assert np.all(np.isfinite(cl))

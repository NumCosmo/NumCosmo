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

Tests the kernel interface, particularly comparing limber_z and limber
implementations via get_eval_vectorized.
"""

import pytest
from numpy.testing import assert_allclose
import numpy as np

from numcosmo_py import Ncm, Nc
from numcosmo_py.cosmology import Cosmology

Ncm.cfg_init()


@pytest.fixture(name="cosmology", scope="module")
def fixture_cosmology() -> Cosmology:
    """Create a default cosmology for testing."""
    return Cosmology.default()


@pytest.fixture(name="integrator", scope="module")
def fixture_integrator() -> Ncm.SBesselIntegrator:
    """Create a spherical Bessel integrator."""
    return Ncm.SBesselIntegratorLevin.new(0, 2000)


@pytest.fixture(name="gal_kernel")
def fixture_gal_kernel(
    cosmology: Cosmology, integrator: Ncm.SBesselIntegrator
) -> Nc.XcorKernelGal:
    """Create a galaxy kernel."""
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
        domagbias=False,
    )
    gal["bparam_0"] = 1.5
    return gal


def test_limber_vs_limber_z(gal_kernel: Nc.XcorKernelGal, cosmology: Cosmology) -> None:
    """Test that limber and limber_z give consistent results via get_eval_vectorized."""
    cosmo = cosmology.cosmo
    dist = cosmology.dist
    ps_ml = cosmology.ps_ml
    gal_kernel.set_l_limber(0)  # All ells use limber
    gal_kernel.prepare(cosmo)

    # Test ell values
    ell_array = np.array([10, 50, 100, 500, 1000])

    # Get vectorized evaluation functions
    for ell in ell_array:
        limber_func = gal_kernel.get_eval(cosmo, ell)
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
                gal_kernel.eval_limber_z_full(cosmo, z, dist, ell)
                * np.sqrt(pk)
                * np.sqrt(np.pi / 2.0 / nu)
                / k
            )
            assert_allclose(
                result,
                result_z,
                rtol=1e-15,
                atol=0.0,
                err_msg=(
                    f"Limber and limber_z differ at ell={ell}, k={k:.3e}, z={z:.3f}"
                ),
            )

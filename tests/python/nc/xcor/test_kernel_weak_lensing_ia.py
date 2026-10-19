#!/usr/bin/env python
#
# test_kernel_weak_lensing_ia.py
#
# Mon Oct 19 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_kernel_weak_lensing_ia.py
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

"""NLA intrinsic alignment in NcXcorKernelWeakLensing against CCL.

The doia property adds the intrinsic shear as a component with kernel
-A_IA(z) C1 rho_crit Omega_m E(z) dn/dz / D(k, z) weighted by
j_l(k chi) / (k chi)^2 and the spin-2 factor, following CCL's
WeakLensingTracer with use_A_ia (transfer -A C1 rho_m / D, der_bessel = -1,
der_angles = 2). The growth factor is read from the power spectrum,
D(k, z) = sqrt(P_D(k, z) / P_D(k, 0)), P_D the growth-powspec property
or the kernel's power spectrum; for a linear spectrum it is CCL's D(z).
"""

import gc

import numpy as np
import pytest
from numpy.testing import assert_allclose

pytest.importorskip("pyccl")

import pyccl

from numcosmo_py import Nc, Ncm
from numcosmo_py.ccl.nc_ccl import create_nc_obj, CCLParams
from numcosmo_py.ccl import two_point

pytestmark = [pytest.mark.ccl, pytest.mark.xcor, pytest.mark.xdist_group("ccl_ia")]

Ncm.cfg_init()

Z_LEN = 1000
MU, SIGMA = 0.8, 0.15
IA_A, IA_ETA, IA_Z_PIV = 1.5, 1.1, 0.62


@pytest.fixture(name="ccl_cosmo", scope="module")
def fixture_ccl_cosmo():
    """CCL cosmology matching the NumCosmo one, linear power spectrum."""
    CCLParams.set_high_prec_params()
    return pyccl.Cosmology(
        Omega_c=0.25,
        Omega_b=0.05,
        h=0.7,
        n_s=0.96,
        sigma8=0.8,
        transfer_function="eisenstein_hu",
        matter_power_spectrum="linear",
    )


@pytest.fixture(name="cosmology", scope="module")
def fixture_cosmology(ccl_cosmo):
    """NumCosmo cosmology matching the CCL one."""
    return create_nc_obj(ccl_cosmo)


@pytest.fixture(name="cosmology_halofit", scope="module")
def fixture_cosmology_halofit():
    """NumCosmo cosmology with Takahashi halofit over the linear spectrum."""
    return create_nc_obj(
        pyccl.Cosmology(
            Omega_c=0.25,
            Omega_b=0.05,
            h=0.7,
            n_s=0.96,
            sigma8=0.8,
            transfer_function="eisenstein_hu",
            matter_power_spectrum="halofit",
        )
    )


@pytest.fixture(name="ccl_calculator_halofit", scope="module")
def fixture_ccl_calculator_halofit(cosmology_halofit):
    """CCL in calculator mode on NumCosmo's background, linear growth and halofit tables.

    CCL's NLA divides by the growth factor of the cosmology, here NumCosmo's
    linear one, and integrates with the nonlinear table.
    """
    CCLParams.set_high_prec_params()
    cosmo = cosmology_halofit.cosmo
    dist, ps_lin, ps_nl = (
        cosmology_halofit.dist,
        cosmology_halofit.ps_ml,
        cosmology_halofit.ps_mnl,
    )
    gf = ps_lin.peek_gf()
    z_bg = np.linspace(0.0, 10.0, 2000)
    a_bg = (1.0 / (1.0 + z_bg))[::-1]
    chi = np.array([dist.comoving(cosmo, z) for z in z_bg]) * cosmo.RH_Mpc()
    E = np.array([cosmo.E(z) for z in z_bg])
    D = np.array([gf.eval(cosmo, z) for z in z_bg]) / gf.eval(cosmo, 0.0)
    f = np.array([-(1.0 + z) * gf.eval_deriv(cosmo, z) for z in z_bg]) / (
        D * gf.eval(cosmo, 0.0)
    )
    k_pk = np.geomspace(1.0e-5, 1.0e2, 1000)
    z_pk = np.linspace(0.0, 3.0, 150)
    k_v = Ncm.Vector.new_array(k_pk.tolist())
    out = Ncm.Vector.new(len(k_pk))

    def table(ps):
        rows = []
        for z in z_pk[::-1]:
            ps.eval_vec(cosmo, float(z), k_v, out)
            rows.append(np.array(out.dup_array()))
        return np.array(rows)

    a_pk = (1.0 / (1.0 + z_pk))[::-1]
    return pyccl.CosmologyCalculator(
        Omega_c=0.25,
        Omega_b=0.05,
        h=0.7,
        n_s=0.96,
        sigma8=0.8,
        background={"a": a_bg, "chi": chi[::-1], "h_over_h0": E[::-1]},
        growth={"a": a_bg, "growth_factor": D[::-1], "growth_rate": f[::-1]},
        pk_linear={"a": a_pk, "k": k_pk, "delta_matter:delta_matter": table(ps_lin)},
        pk_nonlin={"a": a_pk, "k": k_pk, "delta_matter:delta_matter": table(ps_nl)},
    )


@pytest.fixture(name="dndz_arrays", scope="module")
def fixture_dndz_arrays():
    """Gaussian source bin shared by both codes."""
    z_a = np.linspace(0.0, 2.0, Z_LEN)
    nz_a = np.exp(-0.5 * ((z_a - MU) / SIGMA) ** 2)
    return z_a, nz_a


def _nc_wl_kernel(
    cosmology, dndz_arrays, *, doia: bool, lmax: int, ps=None, growth_ps=None
):
    z_a, nz_a = dndz_arrays
    dndz = Ncm.Spline.new_array(
        Ncm.SplineCubicNotaknot.new(), z_a.tolist(), nz_a.tolist(), True
    )
    kernel = Nc.XcorKernelWeakLensing(
        dist=cosmology.dist,
        powspec=cosmology.ps_ml if ps is None else ps,
        dndz=dndz,
        nbar=1.0,
        intr_shear=0.0,
        doia=doia,
        growth_powspec=growth_ps,
        integrator=Ncm.SBesselIntegratorLevin.new(2, lmax),
    )
    kernel["ia_A"] = IA_A
    kernel["ia_eta"] = IA_ETA
    kernel["ia_z_piv"] = IA_Z_PIV
    return kernel


def _ccl_wl_tracer(ccl_cosmo, dndz_arrays, *, has_ia: bool):
    z_a, nz_a = dndz_arrays
    ia_bias = None
    if has_ia:
        ia_bias = (z_a, IA_A * ((1.0 + z_a) / (1.0 + IA_Z_PIV)) ** IA_ETA)
    return pyccl.WeakLensingTracer(
        ccl_cosmo, dndz=(z_a, nz_a), ia_bias=ia_bias, use_A_ia=True, n_samples=Z_LEN
    )


def _nc_cl(
    cosmology, kernel, ells, l_limber, method=Nc.XcorMethod.KERNEL_EXACT, ps=None
):
    kernel.set_l_limber(l_limber)
    kernel.prepare(cosmology.cosmo)
    lmin, lmax = int(ells[0]), int(ells[-1])
    res = Ncm.Vector.new(lmax - lmin + 1)
    xcor = Nc.Xcor.new(cosmology.dist, cosmology.ps_ml if ps is None else ps, method)
    xcor.prepare(cosmology.cosmo)
    xcor.compute(kernel, kernel, cosmology.cosmo, lmin, lmax, res)
    return np.array(res.dup_array())[np.asarray(ells) - lmin]


def _ccl_limber(ccl_cosmo, tracer, ells):
    psp = ccl_cosmo.get_linear_power()
    return pyccl.angular_cl(
        ccl_cosmo, tracer, tracer, ells, p_of_k_a_lin=psp, p_of_k_a=psp
    )


def test_ia_limber_matches_ccl_limber(ccl_cosmo, cosmology, dndz_arrays) -> None:
    """Both Limber paths, with and without IA, against CCL's Limber.

    LIMBER_Z_GSL evaluates the kernel's eval_limber_z, the Limber tier of
    KERNEL_EXACT its components, so the two implementations of the IA term
    are compared with CCL separately. The IA part (total minus lensing only)
    is checked on its own, as at these multipoles it is a few percent of the
    total. Measured: total at most 4.2e-8 (z GSL) and 1.1e-6 (Limber tier),
    IA part at most 7.5e-7 on both.
    """
    ells = np.array([100, 300, 1000])

    plain = _nc_wl_kernel(cosmology, dndz_arrays, doia=False, lmax=int(ells[-1]))
    ia = _nc_wl_kernel(cosmology, dndz_arrays, doia=True, lmax=int(ells[-1]))
    ref_plain = _ccl_limber(
        ccl_cosmo, _ccl_wl_tracer(ccl_cosmo, dndz_arrays, has_ia=False), ells
    )
    ref_ia = _ccl_limber(
        ccl_cosmo, _ccl_wl_tracer(ccl_cosmo, dndz_arrays, has_ia=True), ells
    )

    for method, l_limber in (
        (Nc.XcorMethod.LIMBER_Z_GSL, 0),
        (Nc.XcorMethod.KERNEL_EXACT, 0),
    ):
        got_plain = _nc_cl(cosmology, plain, ells, l_limber, method)
        got_ia = _nc_cl(cosmology, ia, ells, l_limber, method)
        assert_allclose(got_ia, ref_ia, rtol=1.0e-5)
        assert_allclose(got_ia - got_plain, ref_ia - ref_plain, rtol=1.0e-5)


def test_ia_halofit_limber_matches_ccl_with_linear_growth(
    ccl_calculator_halofit, cosmology_halofit, dndz_arrays
) -> None:
    """Halofit spectrum with the linear spectrum as growth-powspec against CCL's NLA.

    CCL divides by the linear growth and integrates P_NL; growth-powspec set
    to the linear spectrum gives the same D. Left unset, D would be the ratio
    of the halofit spectra, which carries the nonlinear evolution. Both Limber
    paths, as in the linear test. Measured on the IA part with the linear
    growth-powspec: at most 4.6e-7 on both; with it unset 1.5e-3, 1.8e-2 and
    2.7e-3 at ell 100, 300, 1000, checked to fail.
    """
    ells = np.array([100, 300, 1000])
    ps_nl, ps_lin = cosmology_halofit.ps_mnl, cosmology_halofit.ps_ml
    calc = ccl_calculator_halofit

    def ccl(has_ia):
        tracer = _ccl_wl_tracer(calc, dndz_arrays, has_ia=has_ia)
        return pyccl.angular_cl(calc, tracer, tracer, ells)

    def nc(doia, growth_ps, method):
        kernel = _nc_wl_kernel(
            cosmology_halofit,
            dndz_arrays,
            doia=doia,
            lmax=int(ells[-1]),
            ps=ps_nl,
            growth_ps=growth_ps,
        )
        return _nc_cl(cosmology_halofit, kernel, ells, 0, method, ps=ps_nl)

    ref_ia = ccl(True) - ccl(False)
    for method in (Nc.XcorMethod.LIMBER_Z_GSL, Nc.XcorMethod.KERNEL_EXACT):
        plain = nc(False, None, method)
        got_ia = nc(True, ps_lin, method) - plain
        self_ia = nc(True, None, method) - plain

        assert_allclose(got_ia, ref_ia, rtol=1.0e-5)
        assert not np.allclose(self_ia, ref_ia, rtol=1.0e-3, atol=0.0)


def test_ia_exact_matches_ccl_kernel_bridge_at_low_ell(
    ccl_cosmo, cosmology, dndz_arrays
) -> None:
    """Exact-tier C_ell with IA against the CCL-kernel bridge where Limber fails.

    numcosmo_py.ccl.two_point.angular_cl integrates CCL's own lensing and IA
    kernels (der_bessel = -1, der_angles = 2) with the library, so the two
    sides share the integrator and nothing else. Measured deviations 9.8e-7,
    1.9e-6, 1.7e-8 at ell 2, 5, 10.
    """
    ells = np.array([2, 5, 10])

    kernel = _nc_wl_kernel(cosmology, dndz_arrays, doia=True, lmax=int(ells[-1]))
    got = _nc_cl(cosmology, kernel, ells, l_limber=-1)

    tracer = _ccl_wl_tracer(ccl_cosmo, dndz_arrays, has_ia=True)
    ref = two_point.angular_cl(cosmology, tracer, tracer, ells)

    assert_allclose(got, ref, rtol=1.0e-5)


def test_ia_component_list_and_properties(cosmology, dndz_arrays) -> None:
    """doia adds one component beside the lensing one; the parameters read back."""
    kernel = _nc_wl_kernel(cosmology, dndz_arrays, doia=True, lmax=10)
    assert kernel.get_property("doia") is True
    assert len(kernel.get_component_list()) == 2
    assert kernel["ia_A"] == IA_A
    assert kernel["ia_eta"] == IA_ETA
    assert kernel["ia_z_piv"] == IA_Z_PIV

    plain = Nc.XcorKernelWeakLensing(
        dist=cosmology.dist,
        powspec=cosmology.ps_ml,
        dndz=kernel.props.dndz,
        nbar=1.0,
        intr_shear=0.0,
    )
    assert plain.get_property("doia") is False
    assert len(plain.get_component_list()) == 1
    assert plain["ia_A"] == Nc.XCOR_KERNEL_WEAK_LENSING_DEFAULT_IA_A
    assert plain["ia_eta"] == Nc.XCOR_KERNEL_WEAK_LENSING_DEFAULT_IA_ETA
    assert plain["ia_z_piv"] == Nc.XCOR_KERNEL_WEAK_LENSING_DEFAULT_IA_Z_PIV
    assert plain.get_property("growth-powspec") is None

    with_growth = _nc_wl_kernel(
        cosmology, dndz_arrays, doia=True, lmax=10, growth_ps=cosmology.ps_ml
    )
    assert with_growth.get_property("growth-powspec") is cosmology.ps_ml

    # dropping the kernel disposes the IA component and its data
    del kernel, plain, with_growth
    gc.collect()


def test_ia_kernel_survives_serialization(cosmology, dndz_arrays) -> None:
    """A serialized copy keeps doia and the IA parameters and gives the same C_ell.

    NcXcorSolver hands each OpenMP thread NcmSerialize-duplicated objects, so
    a construct-only property that did not survive the round trip would
    silently drop the IA term there.
    """
    ells = np.array([2, 5, 10])
    kernel = _nc_wl_kernel(cosmology, dndz_arrays, doia=True, lmax=int(ells[-1]))

    ser = Ncm.Serialize.new(Ncm.SerializeOpt.KEEP_SHARED)
    copy = ser.dup_obj(kernel)

    assert copy.get_property("doia") is True
    assert len(copy.get_component_list()) == 2
    assert copy["ia_A"] == IA_A

    with_growth = _nc_wl_kernel(
        cosmology,
        dndz_arrays,
        doia=True,
        lmax=int(ells[-1]),
        growth_ps=cosmology.ps_ml,
    )
    assert ser.dup_obj(with_growth).get_property("growth-powspec") is not None

    original = _nc_cl(cosmology, kernel, ells, l_limber=-1)
    duplicate = _nc_cl(cosmology, copy, ells, l_limber=-1)

    assert_allclose(duplicate, original, rtol=1.0e-12)

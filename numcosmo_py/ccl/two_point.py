#
# two_point.py
#
# Tue Mar 4 11:14:17 2025
# Copyright  2025  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# two_point.py
# Copyright (C) 2025 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

"""Two-point correlation functions."""

import numpy as np
import pyccl

from numcosmo_py import Nc, Ncm
from numcosmo_py.cosmology import Cosmology


def compute_kernel(
    tracer: pyccl.Tracer,
    cosmology: Cosmology,
    ell: float,
    *,
    reltol: float = 0.0,
    lk: float = 0.0,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Compute the kernel for a given tracer.

    :param reltol: when positive, resample the tracer's components with a B-spline whose
        order is derived from this tolerance instead of a cubic spline. CCL hands kernels
        over as fixed tables, and a cubic reconstruction of a table stalls around 1e-9
        relative -- well above what the Bessel solver delivers -- so anything demanding
        better than that must not go through a cubic.
    :param lk: natural log of the wavenumber at which the transfer is evaluated. CCL's
        transfer is k-dependent; the default of 0.0 keeps the separable approximation.
    """
    cosmo = cosmology.cosmo
    dist = cosmology.dist
    Wchi_list, chi_list = tracer.get_kernel()
    assert chi_list is not None
    chi_a = np.array(chi_list[0])
    for chi_a_i, Wchi_a_i in zip(chi_list, Wchi_list):
        if reltol > 0.0:
            s: Ncm.Spline = Ncm.SplineBSpline.new_tol(reltol, 0.0)
            s.set_array(chi_a_i, Wchi_a_i, True)
        else:
            s = Ncm.Spline.new_array(
                Ncm.SplineCubicNotaknot.new(), chi_a_i, Wchi_a_i, True
            )
        Wchi_a_i[:] = np.array([s.eval(chi) for chi in chi_a])

    RH_Mpc = cosmo.RH_Mpc()
    nu = ell + 0.5

    bessel_factors_list = []

    for der_bessel in tracer.get_bessel_derivative():
        match der_bessel:
            case 0:
                bessel_factors_list.append(1.0)
            case -1:
                # j_l(x)/x^2, folded into the kernel as CCL does.
                bessel_factors_list.append(1.0 / nu**2)
            case 1 | 2:
                raise NotImplementedError(
                    f"der_bessel = {der_bessel} needs derivatives of j_l, which the "
                    "kernel-folding form used here cannot express."
                )
            case _:
                raise ValueError(f"unsupported der_bessel = {der_bessel}")

    z_a = np.array([dist.inv_comoving(cosmo, chi / RH_Mpc) for chi in chi_a])
    a_array = 1.0 / (1.0 + np.array(z_a))
    transfers_list = tracer.get_transfer(lk, a_array)
    ell_factors_list = tracer.get_f_ell(ell)

    Wtotal = np.zeros_like(chi_a)
    for Wchi_a, transfer, ell_factor, bessel_factor in zip(
        Wchi_list, transfers_list, ell_factors_list, bessel_factors_list
    ):
        assert Wchi_a is not None
        assert transfer is not None
        Wtotal += Wchi_a * transfer * ell_factor * bessel_factor

    return z_a[1:-1], chi_a[1:-1], Wtotal[1:-1]


def bessel_transform(
    tracer: pyccl.Tracer,
    cosmology: Cosmology,
    ell: int,
    k_a: np.ndarray,
    *,
    reltol: float = 1.0e-8,
    k_dependent: bool = False,
) -> np.ndarray:
    r"""Compute :math:`\Delta_\ell(k) = \int d\chi\, W(\chi)\, j_\ell(k\chi)`.

    The oscillatory integral is done with NumCosmo's Levin-type solver, which is told a
    tolerance rather than a resolution.

    :param k_a: wavenumbers in Mpc^-1.
    :param reltol: requested relative accuracy, applied to the kernel reconstruction and
        to the oscillatory solve alike.
    :param k_dependent: evaluate the tracer's transfer at each k instead of at k = 0.
        CCL's transfer is a function of (lk, a); the separable default discards that.
    """
    cosmo = cosmology.cosmo
    gf = Nc.GrowthFunc.new()
    gf.prepare_if_needed(cosmo)

    integ = Ncm.SBesselIntegratorLevin.new(ell, ell)
    integ = Ncm.SBesselIntegratorLevin.new_full(
        ell,
        ell,
        integ.get_y_knots_min(),
        integ.get_y_knots_max(),
        integ.get_n_knots(),
        integ.get_ell_cache_max(),
        reltol,
        integ.get_cheb_min_order(),
        reltol,
    )
    out = np.zeros_like(k_a)
    res = Ncm.Vector.new(1)
    cache: dict[float, tuple[Ncm.Spline, float, float]] = {}

    for i, k in enumerate(k_a):
        lk = float(np.log(k)) if k_dependent else 0.0
        key = lk if k_dependent else 0.0
        if key not in cache:
            z_a, chi_a, W_a = compute_kernel(
                tracer, cosmology, ell, reltol=reltol, lk=lk
            )
            # CCL's transfer carries no growth -- it evaluates the 2-D P(k, a) instead.
            # Pairing the kernel with P(k, 0) therefore requires folding D(z) in here.
            # Note this is the separable approximation P(k, z) = D(z)^2 P(k, 0); the
            # native NumCosmo path does not make that split.
            W_a = W_a * np.array([gf.eval(cosmo, float(zz)) for zz in z_a])
            sp = Ncm.SplineBSpline.new_tol(reltol, 0.0)
            sp.set_array(chi_a.tolist(), W_a.tolist(), True)
            cache[key] = (sp, float(chi_a[0]), float(chi_a[-1]))
        sp, chi_min, chi_max = cache[key]

        integ.integrate(
            lambda _p, chi, _k, _s=sp: _s.eval(chi),
            chi_min,
            chi_max,
            float(k),
            res,
            None,
        )
        out[i] = res.get(0)

    return out


def angular_cl(
    cosmology: Cosmology,
    tracer1: pyccl.Tracer,
    tracer2: pyccl.Tracer,
    ells: np.ndarray,
    *,
    reltol: float = 1.0e-8,
    n_k: int = 512,
    n_osc: int = 8,
    k_dependent: bool = False,
) -> np.ndarray:
    r"""Non-Limber angular power spectrum of two CCL tracers, solved by NumCosmo.

    .. math::
        C_\ell = \frac{2}{\pi}\int dk\, k^2 P_{\rm lin}(k)\,
                 \Delta^1_\ell(k)\,\Delta^2_\ell(k)

    The tracers, their kernels and the power spectrum all come from CCL; only the
    oscillatory integral and the reconstruction of the tabulated kernels are NumCosmo's.
    This makes the comparison against :func:`pyccl.angular_cl` an apples-to-apples one,
    and gives CCL users a non-Limber path driven by a tolerance.

    :param reltol: requested relative accuracy of the oscillatory solve.
    :param n_k: floor on the number of wavenumbers in the outer integral.
    :param n_osc: samples per oscillation of Delta_l(k) in the outer integral.
    :param k_dependent: keep the tracers' k-dependent transfer rather than evaluating it
        at k = 0. Costs one kernel rebuild per k.
    """
    ps_lin = cosmology.ps_ml
    cosmo = cosmology.cosmo
    ps_lin.prepare_if_needed(cosmo)
    out = np.zeros(len(ells), dtype=float)

    for i, ell in enumerate(ells):
        _, chi_a, W_a = compute_kernel(tracer1, cosmology, float(ell), reltol=reltol)
        nu = ell + 0.5

        # Effective support of the kernel. The tabulated chi range runs far past where
        # the kernel is non-negligible, and using its endpoints would put the turning
        # point nu/chi decades away from where the transform actually lives.
        w_abs = np.abs(W_a)
        keep = np.nonzero(w_abs > 1.0e-10 * w_abs.max())[0]
        chi_lo, chi_hi = float(chi_a[keep[0]]), float(chi_a[keep[-1]])

        # j_l(k chi) turns over at k chi = nu, and below that the transform is
        # evanescent; above nu/chi_lo it decays with the kernel's smoothness.
        k_lo = 0.2 * nu / chi_hi
        k_hi = 8.0 * nu / chi_lo

        # Delta_l(k) oscillates in k with period pi/chi_hi, so the outer grid has to
        # resolve that; a log grid spanning decades cannot. This is the one place the
        # method still asks for a resolution rather than deriving it -- see n_k_min.
        n_k_min = int(np.ceil(n_osc * (k_hi - k_lo) * chi_hi / np.pi))
        n_k_use = max(n_k, n_k_min)
        k_a = np.linspace(k_lo, k_hi, n_k_use)

        d1 = bessel_transform(
            tracer1, cosmology, int(ell), k_a, reltol=reltol, k_dependent=k_dependent
        )
        if tracer2 is tracer1:
            d2 = d1
        else:
            d2 = bessel_transform(
                tracer2,
                cosmology,
                int(ell),
                k_a,
                reltol=reltol,
                k_dependent=k_dependent,
            )

        # NumCosmo's power spectrum already takes k in Mpc^-1 and returns Mpc^3.
        pk = np.array([ps_lin.eval(cosmo, 0.0, float(k)) for k in k_a])
        integrand = k_a**2 * pk * d1 * d2

        sp = Ncm.SplineCubicNotaknot.new()
        sp.set_array(k_a.tolist(), integrand.tolist(), True)
        out[i] = (2.0 / np.pi) * sp.eval_integ(k_lo, k_hi)

    return out

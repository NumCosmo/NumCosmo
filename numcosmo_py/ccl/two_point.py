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

from typing import Callable

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
                # j_l(x)/x^2 at the Limber point x = nu. A single-ell kernel
                # cannot carry the exact 1/(k chi)^2; angular_cl does.
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


def tracer_components(
    tracer: pyccl.Tracer,
    cosmology: Cosmology,
    *,
    reltol: float = 1.0e-8,
    lk: float = 0.0,
) -> tuple[np.ndarray, list[np.ndarray], list[int]]:
    """Per-component kernels on a common chi grid, with the ell-dependent parts removed.

    The multipole enters a component only through scalars -- ``get_f_ell`` and the
    ``1/nu^2`` that accompanies ``der_bessel = -1`` -- so they can be stripped here and
    reapplied after the solve. That is what lets one factorisation serve a whole block
    of multipoles, which is the structure the solver is built around.

    Growth is folded in here: CCL's transfer carries none, because CCL evaluates the
    2-D P(k, a) rather than pairing a kernel with P(k, 0).
    """
    cosmo = cosmology.cosmo
    dist = cosmology.dist
    RH_Mpc = cosmo.RH_Mpc()
    gf = Nc.GrowthFunc.new()
    gf.prepare_if_needed(cosmo)

    Wchi_list, chi_list = tracer.get_kernel()
    assert chi_list is not None
    chi_a = np.array(chi_list[0])

    resampled = []
    for chi_a_i, Wchi_a_i in zip(chi_list, Wchi_list):
        sp: Ncm.Spline = Ncm.SplineBSpline.new_tol(reltol, 0.0)
        sp.set_array(list(chi_a_i), list(Wchi_a_i), True)
        resampled.append(np.array([sp.eval(float(c)) for c in chi_a]))

    z_a = np.array([dist.inv_comoving(cosmo, float(c) / RH_Mpc) for c in chi_a])
    a_array = 1.0 / (1.0 + z_a)
    D_a = np.array([gf.eval(cosmo, float(zz)) for zz in z_a])
    transfers_list = tracer.get_transfer(lk, a_array)
    der_bessel = list(tracer.get_bessel_derivative())

    comps = []
    for W_a, transfer, der in zip(resampled, transfers_list, der_bessel):
        if der not in (-1, 0, 1, 2):
            raise NotImplementedError(f"unsupported der_bessel = {der}")
        comps.append(W_a * np.asarray(transfer) * D_a)

    # Strip the endpoints, as compute_kernel does: CCL pads the grid.
    return chi_a[1:-1], [c[1:-1] for c in comps], der_bessel


def _kernel_callback(sp: Ncm.Spline) -> Callable[[object, float, float], float]:
    """Radial kernel callback for the integrator: the resampled table itself."""

    def kernel_cb(_p: object, chi: float, _k: float) -> float:
        return sp.eval(chi)

    return kernel_cb


def _spin2_kernel_callback(
    sp: Ncm.Spline, slope: float, chi0: float
) -> Callable[[object, float, float], float]:
    """Radial kernel callback carrying the spin-2 weight 1 / chi^2.

    Below the table's first sample ``chi0`` the kernel is continued linearly,
    ``W = slope * chi``, so the callback returns ``slope / chi`` there.
    """

    def kernel_cb(_p: object, chi: float, _k: float) -> float:
        if chi < chi0:
            return slope / chi

        return sp.eval(chi) / chi**2

    return kernel_cb


def block_transformer(
    tracer: pyccl.Tracer,
    cosmology: Cosmology,
    ell_min: int,
    ell_max: int,
    *,
    reltol: float = 1.0e-8,
) -> Callable[[np.ndarray], np.ndarray]:
    r"""Build the batched :math:`\Delta_\ell(k)` evaluator for one multipole block.

    The kernel splines and the integrator (one operator factorisation for the whole
    block) are set up once; the returned callable evaluates every multipole of the
    block on any wavenumber array, so an outer integral can grow its grid in pieces
    without paying the setup again.

    :returns: ``f(k_a) -> array`` of shape ``(ell_max - ell_min + 1, len(k_a))``.
    """
    n_ell = ell_max - ell_min + 1
    chi_a, comps, der_bessel = tracer_components(tracer, cosmology, reltol=reltol)
    chi_min, chi_max = float(chi_a[0]), float(chi_a[-1])

    proto = Ncm.SBesselIntegratorLevin.new(ell_min, ell_max)
    integ = Ncm.SBesselIntegratorLevin.new_full(
        ell_min,
        ell_max,
        proto.get_y_knots_min(),
        proto.get_y_knots_max(),
        proto.get_n_knots(),
        proto.get_ell_cache_max(),
        reltol,
        proto.get_cheb_min_order(),
        reltol,
    )
    res = Ncm.Vector.new(n_ell)
    ells = np.arange(ell_min, ell_max + 1)

    pieces: list[tuple[Callable, int, np.ndarray, float]] = []
    for comp, der in zip(comps, der_bessel):
        sp = Ncm.SplineBSpline.new_tol(reltol, 0.0)
        sp.set_array(chi_a.tolist(), comp.tolist(), True)
        # ell-dependent scalars, reapplied after the batched solve
        f_ell = np.array([tracer.get_f_ell(float(ell)) for ell in ells]).reshape(
            len(ells), -1
        )
        f_c = f_ell[:, 0]
        # der_bessel = -1 is CCL's spin-2 weight j_l(k chi) / (k chi)^2. It is
        # separable: 1/chi^2 multiplies the radial kernel inside the solve and
        # 1/k^2 comes out of the transform, so it is integrated exactly.
        # Replacing it by 1/nu^2 (its Limber value, x = nu) is 6% off at
        # ell = 2 on a lensing kernel. The weight is applied in the callback,
        # not to the table: W ~ chi near the origin, so W / chi^2 ~ 1 / chi
        # cannot be tabulated to the requested tolerance, while the solver's
        # panel fits in y = k chi handle it (NumCosmo's own lensing kernel
        # carries the same 1 / chi).
        if der == -1:
            # CCL's kernel table starts a few Mpc from the observer, where a
            # lensing kernel is already linear in chi. The 1 / chi left of the
            # first sample still contributes: at the ell = 2 peak the
            # transform is oscillation-suppressed and the omitted piece is
            # 0.2% of it (0.4% of C_2). Continue W linearly to the origin and
            # integrate from far inside the first cell.
            chi_from = 1.0e-3 * chi_min
            kernel_cb = _spin2_kernel_callback(sp, float(comp[0]) / chi_min, chi_min)
        else:
            chi_from = chi_min
            kernel_cb = _kernel_callback(sp)

        # positive orders are the Bessel-derivative weights (2 is CCL's RSD)
        pieces.append((kernel_cb, der, f_c, chi_from))

    def evaluate(k_a: np.ndarray) -> np.ndarray:
        out = np.zeros((n_ell, len(k_a)))
        for kernel_cb, der, f_c, chi_from in pieces:
            deriv = der if der > 0 else 0
            block = np.zeros((n_ell, len(k_a)))
            for ik, k in enumerate(k_a):
                integ.integrate_deriv(
                    kernel_cb, chi_from, chi_max, float(k), deriv, res, None
                )
                for j in range(n_ell):
                    block[j, ik] = res.get(j)

            if der == -1:
                block /= k_a[None, :] ** 2

            out += block * f_c[:, None]

        return out

    return evaluate


def bessel_transform_block(
    tracer: pyccl.Tracer,
    cosmology: Cosmology,
    ell_min: int,
    ell_max: int,
    k_a: np.ndarray,
    *,
    reltol: float = 1.0e-8,
) -> np.ndarray:
    r"""Batched :math:`\Delta_\ell(k)` for every multipole in ``[ell_min, ell_max]``.

    One operator factorisation serves the whole block, and one call per wavenumber
    returns every multipole in it. This is the reuse structure the method is designed
    around; solving multipole by multipole discards it.

    :returns: array of shape ``(ell_max - ell_min + 1, len(k_a))``.
    """
    return block_transformer(tracer, cosmology, ell_min, ell_max, reltol=reltol)(k_a)


def bessel_transform(
    tracer: pyccl.Tracer,
    cosmology: Cosmology,
    ell: int,
    k_a: np.ndarray,
    *,
    reltol: float = 1.0e-8,
    k_dependent: bool = False,
) -> np.ndarray:
    r"""Single-multipole :math:`\Delta_\ell(k)`.

    A convenience wrapper over :func:`bessel_transform_block` with a block of one. For
    more than a handful of multipoles use the batched form directly, which shares one
    factorisation across the block.
    """
    if k_dependent:
        raise NotImplementedError(
            "k-dependent transfers need one kernel per k; use tracer_components(lk=...)"
        )
    return bessel_transform_block(tracer, cosmology, ell, ell, k_a, reltol=reltol)[0]


_OUTER_CHUNK = 512


def _grow_outer_grid(
    integrand_of: Callable[[np.ndarray], np.ndarray],
    k_lo: float,
    k_hi: float,
    dk: float,
    *,
    chunk: int,
    tol: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Sample the outer integrand on a step-``dk`` grid whose ends the integrand sets.

    Both ends are placed by a bound on the *integrated* tail left outside,
    relative to the integral accumulated so far, per multipole. Upward from
    ``k_lo`` in chunks: past the peak the integrand falls at least as fast as
    ``k^-3`` (``P(k)/k^2`` for a spin-2 pair, faster for density kernels), so
    the tail beyond ``k_end`` is at most ``|I(k_end)| k_end / 2``; stop once that
    is below ``tol`` of the integral for every multipole, or at ``k_hi``.
    Downward from ``k_lo`` by halving: the integrand rises at least as
    ``k^(2 + n_s)``, so the piece below ``k_0`` is at most ``|I(k_0)| k_0 / 3``.
    A pointwise cut at ``tol`` of the peak is not enough here: on a lensing
    pair it left a tail worth 2e-4 of C_10.

    :returns: ``(k_a, vals)`` with ``vals`` of shape ``(n_ell, len(k_a))``.
    """

    def accumulated(k_a: np.ndarray, vals: np.ndarray) -> np.ndarray:
        return np.trapezoid(np.abs(vals), k_a, axis=1)

    k_a = np.empty(0)
    vals = np.empty((0, 0))
    k_start = k_lo
    while k_start < k_hi:
        k_stop = min(k_start + chunk * dk, k_hi)
        n_pts = max(2, int(np.round((k_stop - k_start) / dk)) + 1)
        k_c = np.linspace(k_start, k_stop, n_pts)
        if k_a.size:
            k_c = k_c[1:]  # shared endpoint
        v_c = integrand_of(k_c)
        k_a = np.concatenate([k_a, k_c])
        vals = v_c if vals.size == 0 else np.concatenate([vals, v_c], axis=1)
        k_start = k_stop

        peaked = np.all(np.argmax(np.abs(vals), axis=1) < vals.shape[1] - v_c.shape[1])
        tail = 0.5 * np.abs(vals[:, -1]) * k_a[-1]
        if peaked and np.all(tail < tol * accumulated(k_a, vals)):
            break

    # Downward the integrand is smooth (no oscillation below the first peak),
    # so each halving adds a handful of points rather than a step-dk chunk.
    k_floor = 1.0e-4 * k_lo
    while k_a[0] > k_floor:
        head = np.abs(vals[:, 0]) * k_a[0] / 3.0
        if np.all(head < tol * accumulated(k_a, vals)):
            break
        k_c = np.linspace(0.5 * k_a[0], k_a[0], 9)[:-1]
        v_c = integrand_of(k_c)
        k_a = np.concatenate([k_c, k_a])
        vals = np.concatenate([v_c, vals], axis=1)

    return k_a, vals


def angular_cl(
    cosmology: Cosmology,
    tracer1: pyccl.Tracer,
    tracer2: pyccl.Tracer,
    ells: np.ndarray,
    *,
    reltol: float = 1.0e-8,
    n_k: int = 512,
    n_osc: int = 8,
    block_size: int = 8,
    support_tol: float = 1.0e-6,
) -> np.ndarray:
    r"""Non-Limber angular power spectrum of two CCL tracers, solved by NumCosmo.

    .. math::
        C_\ell = \frac{2}{\pi}\int dk\, k^2 P_{\rm lin}(k)\,
                 \Delta^1_\ell(k)\,\Delta^2_\ell(k)

    The tracers, their kernels and the power spectrum are CCL's; only the oscillatory
    integral and the reconstruction of the tabulated kernels are NumCosmo's.

    Multipoles are processed in blocks that share one operator factorisation. Block
    width has an optimum rather than being as large as possible: the truncation order
    is chosen by a max-reduction across the block, so every member pays the order the
    hardest member needs, and the k-grid must span the whole block's range.

    :param reltol: requested relative accuracy of the oscillatory solve.
    :param n_k: floor on the number of wavenumbers in the outer integral.
    :param n_osc: samples per oscillation of Delta_l(k) in the outer integral.
    :param support_tol: fraction of the kernel's peak below which its tail is
        treated as outside the support, and the fraction of the outer integral
        allowed outside the k grid. This is the dominant cost control, not a
        safety margin: the sampling step and the starting window come from
        ``chi_lo`` and ``chi_hi``, so shrinking ``chi_lo`` raises the point
        count. On the Gaussian test tracer at ell = 32, 1e-10 (the old value)
        admits tail down to chi = 16.7 Mpc and needs 180,720 wavenumbers for
        293 s; 1e-6 needs 3,945 for 2.4 s, and the two C_ell agree to 2.6e-9.
        The grid then grows on the integrand itself, so a spin-2 tracer, whose
        transform tends to a constant at small k and to 1/k^2 at large k, is
        integrated to the same fraction (see :func:`_grow_outer_grid`).
        Tighten it only with a measurement in hand.
    :param block_size: number of consecutive multipoles sharing a factorisation.
        Clamped to what the solver will honour: NC_XCOR_KERNEL_MAX_ELL_BLOCK and the
        integrator's ell_cache_max. Wider is not better -- the truncation order is a
        max-reduction across the block, so every member pays the hardest member's
        order, and one k-grid must span the block's whole range. Measured on eight
        consecutive multipoles at fixed accuracy, the cost turns over at four:
        1 -> 214.5 s, 2 -> 128.0 s, 4 -> 98.9 s, 8 -> 117.1 s.
    """
    ps_lin = cosmology.ps_ml
    cosmo = cosmology.cosmo
    ps_lin.prepare_if_needed(cosmo)
    ells = np.atleast_1d(np.asarray(ells, dtype=int))
    out = np.zeros(len(ells), dtype=float)

    cache_max = Ncm.SBesselIntegratorLevin.new(0, 0).get_ell_cache_max()
    block_size = max(1, min(block_size, Nc.XCOR_KERNEL_MAX_ELL_BLOCK, cache_max))

    chi_a, comps, _ = tracer_components(tracer1, cosmology, reltol=reltol)
    w_abs = np.abs(sum(comps))
    # See support_tol: chi_lo sets the sampling step and the starting window.
    keep = np.nonzero(w_abs > support_tol * w_abs.max())[0]
    chi_lo, chi_hi = float(chi_a[keep[0]]), float(chi_a[keep[-1]])

    def pk_of(k_a: np.ndarray) -> np.ndarray:
        # NumCosmo's power spectrum already takes k in Mpc^-1 and returns Mpc^3.
        return np.array([ps_lin.eval(cosmo, 0.0, float(k)) for k in k_a])

    for start_i in range(0, len(ells), block_size):
        sel = np.arange(start_i, min(start_i + block_size, len(ells)))
        ell_min, ell_max = int(ells[sel[0]]), int(ells[sel[-1]])

        d1_of = block_transformer(tracer1, cosmology, ell_min, ell_max, reltol=reltol)
        d2_of = (
            d1_of
            if tracer2 is tracer1
            else block_transformer(tracer2, cosmology, ell_min, ell_max, reltol=reltol)
        )

        def integrand_of(k_a: np.ndarray) -> np.ndarray:
            return k_a[None, :] ** 2 * pk_of(k_a)[None, :] * d1_of(k_a) * d2_of(k_a)

        # The window a density kernel needs: Delta_l ~ (k chi_hi)^l below k_lo
        # and averaged out above k_hi. A spin-2 (der_bessel = -1) transform
        # tends to a constant at small k and falls as 1/k^2 at large k, so
        # neither bound holds for it; the grid is therefore grown on the
        # integrand itself from this window, downward until its low end is
        # below support_tol of the peak and upward until it stays there.
        nu_lo, nu_hi = ell_min + 0.5, ell_max + 0.5
        k_lo = 0.2 * nu_lo / chi_hi
        k_hi = 8.0 * nu_hi / chi_lo
        # sampling step from the fastest oscillation in k, period 2 pi / chi_hi
        n_k_min = int(np.ceil(n_osc * (k_hi - k_lo) * chi_hi / np.pi))
        dk = (k_hi - k_lo) / max(n_k, n_k_min)

        k_a, vals = _grow_outer_grid(
            integrand_of, k_lo, k_hi, dk, chunk=_OUTER_CHUNK, tol=support_tol
        )

        for j, i in enumerate(sel):
            row = int(ells[i]) - ell_min
            sp = Ncm.SplineCubicNotaknot.new()
            sp.set_array(k_a.tolist(), vals[row].tolist(), True)
            out[i] = (2.0 / np.pi) * sp.eval_integ(float(k_a[0]), float(k_a[-1]))

    return out

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


def _levin_integrator(
    ell_min: int, ell_max: int, reltol: float
) -> Ncm.SBesselIntegratorLevin:
    """A Levin integrator at the library's knot defaults with both tolerances set."""
    proto = Ncm.SBesselIntegratorLevin.new(ell_min, ell_max)

    return Ncm.SBesselIntegratorLevin.new_full(
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

    integ = _levin_integrator(ell_min, ell_max, reltol)
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


_KIND_BY_DERIVATIVES = {
    (0, 0): Nc.XcorKernelTableKind.DENSITY,
    (-1, 2): Nc.XcorKernelTableKind.SHEAR,
    (-1, 1): Nc.XcorKernelTableKind.CONVERGENCE,
    (2, 0): Nc.XcorKernelTableKind.RSD,
}


def tracer_component_tables(
    tracer: pyccl.Tracer,
    cosmology: Cosmology,
    *,
    order: int = 8,
) -> list[Nc.XcorComponentTable]:
    r"""One :class:`Nc.XcorComponentTable` per term of a CCL tracer.

    Each term keeps its own chi grid; a term whose window is identically zero
    is dropped. The window is ``W(chi) * transfer(a) * D(z)``
    -- CCL's kernel, its scale-factor transfer (bias, minus the growth rate for
    RSD, ...) and the linear growth, which CCL leaves out because it evaluates
    the two-dimensional ``P(k, a)`` while the radial family pairs the window
    with ``P(k, 0)``. The pair (``der_bessel``, ``der_angles``) selects the kind,
    which fixes the Bessel weight, the ``1/(k chi)^2`` factor and the ell
    prefactor on the NumCosmo side.

    The transfer is taken at ``lk = 0``: CCL's is formally k-dependent, but for
    the standard tracers it is a function of the scale factor alone. A
    k-dependent transfer would go through the radial family's scale-dependence
    slot, which this adapter does not populate.

    Do not oversample CCL's kernels. Its lensing kernel is a numerical integral
    with its own accuracy, and denser sampling exposes that noise rather than
    more shape: the relative fourth difference of a WeakLensingTracer table is
    2e-7 at 1000 samples and 1.7e-5 at 4000. The closure refinement resolves
    what it is given, so the 4000-sample table costs seven times the 1000-sample
    one for the same C_ell (2e-8 of peak from NumCosmo's own lensing kernel
    either way). CCL's analytic CMB lensing kernel has no such floor.

    :param order: B-spline order of the reconstruction (degree ``order - 1``).
    """
    cosmo = cosmology.cosmo
    dist = cosmology.dist
    RH_Mpc = cosmo.RH_Mpc()
    gf = Nc.GrowthFunc.new()
    gf.prepare_if_needed(cosmo)

    Wchi_list, chi_list = tracer.get_kernel()
    assert chi_list is not None
    # z(chi) comes from the distance object's inverse spline. Extend its reach
    # the way every NcXcorKernel does at construction, so a CMB lensing tracer
    # (z ~ 1100) inverts regardless of the cosmology's dist_z_max.
    dist.compute_inv_comoving(True)
    dist.require_zf(1.0e10)
    dist.prepare_if_needed(cosmo)

    der_bessel = [int(d) for d in tracer.get_bessel_derivative()]
    der_angles = [int(d) for d in tracer.get_angles_derivative()]

    comps = []
    for i, (chi_a, W_a, der, ang) in enumerate(
        zip(chi_list, Wchi_list, der_bessel, der_angles)
    ):
        kind = _KIND_BY_DERIVATIVES.get((der, ang))
        if kind is None:
            raise NotImplementedError(
                f"tracer term {i}: der_bessel = {der}, der_angles = {ang} has no "
                "NcXcorKernelTableKind"
            )

        W_a = np.asarray(W_a, dtype=float)
        if not np.any(W_a):
            # CCL emits an identically zero window for a term that cancels, e.g.
            # magnification at s = 0.4 where 5 s - 2 = 0; it carries nothing.
            continue

        chi_a = np.asarray(chi_a, dtype=float)
        z_a = np.array([dist.inv_comoving(cosmo, float(c) / RH_Mpc) for c in chi_a])
        a_a = 1.0 / (1.0 + z_a)
        transfer = np.asarray(tracer.get_transfer(0.0, a_a)[i], dtype=float)
        D_a = np.array([gf.eval(cosmo, float(zz)) for zz in z_a])
        W = W_a * transfer * D_a

        comps.append(
            Nc.XcorComponentTable.new_full(
                Ncm.Vector.new_array(chi_a.tolist()),
                Ncm.Vector.new_array(W.tolist()),
                kind,
                order,
                False,
            )
        )

    return comps


def tracer_kernel(
    tracer: pyccl.Tracer,
    cosmology: Cosmology,
    *,
    order: int = 8,
    kernel_reltol: float | None = None,
) -> Nc.XcorKernelTable:
    r"""A CCL tracer as a prepared :class:`Nc.XcorKernelTable`, never Limber.

    The kernel owns everything the solver needs: the components, their supports,
    the Levin integrator and the k-space closure. NumCosmo then discovers the k
    domain, places the knots and integrates the outer integral exactly, so the
    only things CCL contributes are the tabulated windows and the growth.

    The Levin integrator is the library default. Loosening its solve tolerance
    (1e-13) to the 1e-8 the raw transforms above use makes the closure
    refinement at 1e-6 chase solve noise on an RSD table and never finish; the
    accuracy knob here is ``kernel_reltol``, not the solve.

    :param order: B-spline order of the window reconstruction.
    :param kernel_reltol: when given, sets both the kernel's closure ``reltol``
        and ``scaled-abstol``; the library default (1e-4) applies otherwise.
    """
    comps = tracer_component_tables(tracer, cosmology, order=order)
    oa = Ncm.ObjArray.new()
    for comp in comps:
        oa.add(comp)

    kernel = Nc.XcorKernelTable.new_from_components(
        cosmology.dist, cosmology.ps_ml, oa, Ncm.SBesselIntegratorLevin.new(0, 8)
    )
    if kernel_reltol is not None:
        kernel.set_reltol(kernel_reltol)
        kernel.set_scaled_abstol(kernel_reltol)
    kernel.set_l_limber(-1)
    kernel.prepare(cosmology.cosmo)

    return kernel


def angular_cl(
    cosmology: Cosmology,
    tracer1: pyccl.Tracer,
    tracer2: pyccl.Tracer,
    ells: np.ndarray,
    *,
    order: int = 8,
    kernel_reltol: float | None = None,
    block_size: int = 8,
) -> np.ndarray:
    r"""Non-Limber angular power spectrum of two CCL tracers, solved by NumCosmo.

    .. math::
        C_\ell = \frac{2}{\pi}\int dk\, k^2 P_{\rm lin}(k)\,
                 \Delta^1_\ell(k)\,\Delta^2_\ell(k)

    The tracers are CCL's; everything else is :class:`Nc.XcorSolver` on
    :func:`tracer_kernel` kernels with the exact kernel-space method: the k
    domain is discovered, the knots are placed adaptively, each kernel's closure
    is shared across every pair and multipole block that needs it, and the
    outer integral is exact on the closure. Multipoles are requested one by one;
    contiguous runs merge into blocks of up to ``block_size``, which share one
    operator factorisation.

    :param order: B-spline order of the window reconstruction.
    :param kernel_reltol: closure tolerance for both kernels; see :func:`tracer_kernel`.
    :param block_size: largest number of consecutive multipoles per block. Clamped
        to what the solver honours: NC_XCOR_KERNEL_MAX_ELL_BLOCK and the
        integrator's ell_cache_max.
    """
    ells = np.atleast_1d(np.asarray(ells, dtype=int))
    cache_max = Ncm.SBesselIntegratorLevin.new(0, 0).get_ell_cache_max()
    block_size = max(1, min(block_size, Nc.XCOR_KERNEL_MAX_ELL_BLOCK, cache_max))

    kernel_1 = tracer_kernel(
        tracer1, cosmology, order=order, kernel_reltol=kernel_reltol
    )
    kernel_2 = (
        kernel_1
        if tracer2 is tracer1
        else tracer_kernel(tracer2, cosmology, order=order, kernel_reltol=kernel_reltol)
    )

    solver = Nc.XcorSolver.new()
    id_1 = solver.register_kernel(kernel_1)
    id_2 = id_1 if kernel_2 is kernel_1 else solver.register_kernel(kernel_2)
    for ell in ells:
        solver.request_cl(id_1, id_2, int(ell), int(ell))
    solver.plan_blocks(block_size)

    xcor = Nc.Xcor.new(cosmology.dist, cosmology.ps_ml, Nc.XcorMethod.KERNEL_EXACT)
    xcor.prepare(cosmology.cosmo)
    solver.solve(xcor, cosmology.cosmo)

    return np.array([solver.get_result(i).get(0) for i in range(len(ells))])

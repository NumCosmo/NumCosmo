#!/usr/bin/env python
#
# curvature_weight.py
#
# Sun Jun 22 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# curvature_weight.py
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

"""Data-driven local curvature weight for spline w(z)/q(z) reconstruction.

The global curvature prior penalizes the curvature norm uniformly along the
reconstruction interval. The data, however, constrains the reconstructed
function *non-locally*: a distance/expansion datum at redshift ``z*`` informs an
integral of ``w`` over ``[0, z*]`` (through ``d_L`` and ``rho_DE``), so "where
the data is" in redshift does not map pointwise onto the features of ``w(z)`` it
constrains. The exact object carrying that mapping is the Fisher information of
the spline knots, ``F_d = J^T C^{-1} J``, whose sensitivity Jacobian ``J``
already contains the full integral kernel.

This module turns ``F_d`` into a smooth local weight ``W(x)`` for the weighted
curvature prior (the ``wlp_*`` MSetFuncs). It uses the *function-space*
resolution: the local data information about the reconstructed value,

    I_d(x) = phi(x)^T F_d phi(x),

where ``phi(x)`` is the spline value basis (``d value / d knot``) at ``x``. The
weight is the resolution against a uniform reference precision ``I_ref``,

    W(x) = I_ref / (I_d(x) + I_ref) = 1 - I_d(x) / (I_d(x) + I_ref),

so ``W -> 1`` (full curvature prior) where the data is blind and ``W -> 0``
(prior relaxed) where the data resolves the function. ``I_ref`` is set from the
information scale of the problem via ``ref_factor`` (a dimensionless knob).

The per-knot Tikhonov resolution diagonal ``diag((F_d + F_p)^{-1} F_d)`` is the
discrete analogue, but on an over-resolved knot set it is rank-deficient and
its diagonal alternates (a checkerboard nullspace artifact); the function-space
information above is its smooth, knot-count-independent counterpart.
"""

from typing import Callable

import numpy as np
import numpy.typing as npt

from numcosmo_py import Nc, Ncm


def _matrix_to_numpy(matrix: Ncm.Matrix) -> npt.NDArray[np.float64]:
    """Copy an ``NcmMatrix`` into a dense, symmetrized NumPy array."""
    rows, cols = matrix.row_len(), matrix.col_len()
    out = np.array(matrix.dup_array(), dtype=float).reshape(rows, cols)
    return 0.5 * (out + out.T)


# Data families with an analytic mean vector + covariance, for which the expected
# Fisher F = J^T C^-1 J is defined. Empirical / non-Gaussian likelihoods (e.g.
# DataBaoEmpiricalFit) have no mean vector and cannot contribute an analytic
# Fisher, so they are skipped when building the weight (they remain in the
# likelihood used for the actual fit/MCMC).
_GAUSS_DATA_TYPES = (Ncm.DataGauss, Ncm.DataGaussCov, Ncm.DataGaussDiag)


def expected_fisher_dataset(dataset: Ncm.Dataset) -> Ncm.Dataset:
    """Sub-dataset of the Gaussian-family data, on which ``fisher_matrix`` works.

    Raises:
        ValueError: if no data object supports an analytic (expected) Fisher.
    """
    subset = Ncm.Dataset.new()
    for i in range(dataset.get_length()):
        data = dataset.peek_data(i)
        if isinstance(data, _GAUSS_DATA_TYPES):
            subset.append_data(data)
    if subset.get_length() == 0:
        raise ValueError(
            "no Gaussian-family data to build the expected Fisher weight from."
        )
    return subset


def knot_fparam_indices(mset: Ncm.MSet, knot_prefix: str) -> list[int]:
    """Free-parameter indices whose names start with ``knot_prefix``.

    The reconstruction knots (``w_i`` for the w-spline, ``qparam_i`` for the
    q-spline) must be set free before ``mset.prepare_fparam_map()``.
    """
    return [
        i
        for i in range(mset.fparam_len())
        if mset.fparam_name(i).startswith(knot_prefix)
    ]


def marginal_knot_fisher(
    fisher: npt.NDArray[np.float64], knot_idx: list[int]
) -> npt.NDArray[np.float64]:
    """Data Fisher on the knots, marginalized over the other free parameters.

    Uses the Schur complement ``F_kk - F_ko F_oo^{-1} F_ok`` so that degeneracies
    between the knots and nuisance/background parameters (``Omegac``, ``H0``,
    SNIa offsets, ...) are propagated, rather than held fixed.
    """
    other_idx = [i for i in range(fisher.shape[0]) if i not in knot_idx]
    f_kk = fisher[np.ix_(knot_idx, knot_idx)]
    if not other_idx:
        return f_kk
    f_ko = fisher[np.ix_(knot_idx, other_idx)]
    f_oo = fisher[np.ix_(other_idx, other_idx)]
    return f_kk - f_ko @ np.linalg.solve(f_oo, f_ko.T)


def value_basis(
    cosmo: Nc.HICosmo,
    knot_idx_in_model: list[int],
    eval_fn: Callable[[Nc.HICosmo, float], float],
    x_grid: npt.NDArray[np.float64],
) -> npt.NDArray[np.float64]:
    """Spline value basis ``phi[k, i] = d value(x_k) / d knot_i`` on ``x_grid``.

    The spline value at ``x`` is linear in its knots, so the basis is exact: it
    is read off by perturbing each knot by one unit from zero. The model's knot
    vector is restored on return.
    """
    n_model = cosmo.len()
    saved = np.array([cosmo.param_get(i) for i in knot_idx_in_model])

    def set_knots(values: npt.NDArray[np.float64]) -> None:
        for pid, val in zip(knot_idx_in_model, values):
            cosmo.param_set(pid, float(val))

    n_knots = len(knot_idx_in_model)
    set_knots(np.zeros(n_knots))
    base = np.array([eval_fn(cosmo, x) for x in x_grid])

    phi = np.zeros((len(x_grid), n_knots))
    unit = np.zeros(n_knots)
    for i in range(n_knots):
        unit[:] = 0.0
        unit[i] = 1.0
        set_knots(unit)
        phi[:, i] = np.array([eval_fn(cosmo, x) for x in x_grid]) - base

    set_knots(saved)
    assert n_model >= n_knots  # sanity: knots are a subset of the model params
    return phi


def information_weight_values(
    info: npt.NDArray[np.float64], ref_factor: float, floor: float
) -> npt.NDArray[np.float64]:
    """Map the local data information ``I_d(x)`` to ``W(x)`` in ``[floor, 1]``.

    ``I_ref = ref_factor * median(I_d > 0)`` sets the crossover scale, so
    ``ref_factor`` is the only dimensionless knob: smaller relaxes the prior more
    eagerly (data needs less information to take over), larger keeps it stronger.
    """
    positive = info[info > 0.0]
    info_ref = ref_factor * (np.median(positive) if positive.size else 1.0)
    weight = info_ref / (info + info_ref)
    return np.clip(weight, floor, 1.0)


def fisher_information_weight(
    dataset: Ncm.Dataset,
    mset: Ncm.MSet,
    cosmo: Nc.HICosmo,
    *,
    knot_prefix: str,
    eval_fn: Callable[[Nc.HICosmo, float], float],
    x_grid: npt.NDArray[np.float64],
    ref_factor: float = 1.0,
    floor: float = 1e-3,
) -> Ncm.Spline:
    """Build the data-driven local curvature weight ``W(x)`` as a spline.

    Args:
        dataset: the data whose Fisher information drives the weight.
        mset: the model set with the reconstruction knots set free and
            ``prepare_fparam_map()`` already called.
        cosmo: the spline cosmology (w-spline or q-spline) inside ``mset``.
        knot_prefix: parameter-name prefix selecting the knots (``"w_"`` or
            ``"qparam_"``).
        eval_fn: evaluates the reconstructed function at a curvature-coordinate
            point (``alpha`` for the w-spline, ``z`` for the q-spline).
        x_grid: the curvature-coordinate grid the weight spline is built on.
        ref_factor: crossover-scale knob (see ``information_weight_values``).
        floor: minimum weight, keeping ``W`` strictly positive for the norm.

    Returns:
        A not-a-knot ``NcmSpline`` ``W(x)`` on ``x_grid``, ready to feed a
        ``wlp_*`` curvature-prior function.
    """
    fisher = _matrix_to_numpy(expected_fisher_dataset(dataset).fisher_matrix(mset))
    knot_fpi = knot_fparam_indices(mset, knot_prefix)
    if not knot_fpi:
        raise ValueError(f"no free knots with prefix {knot_prefix!r} in mset")

    f_d = marginal_knot_fisher(fisher, knot_fpi)

    knot_pid = [
        i for i in range(cosmo.len()) if cosmo.param_name(i).startswith(knot_prefix)
    ]
    phi = value_basis(cosmo, knot_pid, eval_fn, x_grid)
    info = np.einsum("ki,ij,kj->k", phi, f_d, phi)

    weight = information_weight_values(info, ref_factor, floor)
    return Ncm.Spline.new(
        Ncm.SplineCubicNotaknot.new(),
        Ncm.Vector.new_array(x_grid.tolist()),
        Ncm.Vector.new_array(weight.tolist()),
        True,
    )


def wspline_curvature_weight(
    dataset: Ncm.Dataset,
    mset: Ncm.MSet,
    cosmo: Nc.HICosmoDEWSpline,
    *,
    n_grid: int = 64,
    ref_factor: float = 1.0,
    floor: float = 1e-3,
) -> Ncm.Spline:
    """Local curvature weight ``W(alpha)`` for a w-spline reconstruction.

    The curvature lives in ``alpha = ln(1 + z)``; the weight spline shares that
    abscissa, matching the ``NcHICosmoDEWSpline:wlp_*`` functions.
    """
    alpha_f = cosmo.get_alpha().dup_array()[-1]
    x_grid = np.linspace(0.0, alpha_f, n_grid)

    def eval_w(model: Nc.HICosmo, alpha: float) -> float:
        assert isinstance(model, Nc.HICosmoDE)
        return model.w_de(float(np.expm1(alpha)))

    return fisher_information_weight(
        dataset,
        mset,
        cosmo,
        knot_prefix="w_",
        eval_fn=eval_w,
        x_grid=x_grid,
        ref_factor=ref_factor,
        floor=floor,
    )


def qspline_curvature_weight(
    dataset: Ncm.Dataset,
    mset: Ncm.MSet,
    cosmo: Nc.HICosmoQSpline,
    *,
    z_max: float,
    n_grid: int = 64,
    ref_factor: float = 1.0,
    floor: float = 1e-3,
) -> Ncm.Spline:
    """Local curvature weight ``W(z)`` for a q-spline reconstruction.

    The curvature lives in redshift ``z``; the weight spline shares that
    abscissa, matching the ``NcHICosmoQSpline:wlp_*`` functions.
    """
    x_grid = np.linspace(0.0, z_max, n_grid)

    def eval_q(model: Nc.HICosmo, z: float) -> float:
        return model.q(float(z))

    return fisher_information_weight(
        dataset,
        mset,
        cosmo,
        knot_prefix="qparam_",
        eval_fn=eval_q,
        x_grid=x_grid,
        ref_factor=ref_factor,
        floor=floor,
    )

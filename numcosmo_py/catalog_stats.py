#
# catalog_stats.py
#
# Sun Jul 26 00:00:00 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# catalog_stats.py
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

"""Posterior-statistics helpers for catalog parameters and derived quantities.

Resolves catalog parameters and add-values by name or numeric index, parses
"name=parameter" variable bindings, and computes asymmetric confidence-interval
bounds around a median, mode, or arbitrary center. These are CLI-agnostic
helpers shared by the `catalog` CLI commands and by the corner-plot
derived-quantity support in `numcosmo_py.plotting.derived`.
"""

from __future__ import annotations

import enum
from typing import Callable, List, Optional

from . import Ncm

# 1, 2 and 3-sigma two-sided confidence levels.
SIGMA_LEVELS = (0.682689492137086, 0.954499736103642, 0.997300203936740)


class DerivedStat(str, enum.Enum):
    """Which posterior summary statistic to report."""

    MEDIAN = "median"
    MODE = "mode"
    BESTFIT = "bestfit"


def resolve_param(
    mset: Ncm.MSet, nadd_vals: int, name: str
) -> tuple[Optional[Ncm.MSetPIndex], int]:
    """Resolve a free-parameter name or numeric column index.

    Returns the (optional) `Ncm.MSetPIndex` and the full-row column index
    (add-values offset included).
    """
    pi = mset.fparam_get_pi_by_name(name)
    if pi is not None:
        return pi, mset.fparam_get_fpi(pi.mid, pi.pid) + nadd_vals

    if name.isnumeric():
        total_len = mset.fparams_len() + nadd_vals
        pindex = int(name)
        if pindex >= total_len or pindex < 0:
            raise ValueError(f'Invalid parameter index "{name}"=={pindex}.')
        if pindex >= nadd_vals:
            pi = mset.fparam_get_pi(pindex - nadd_vals)
        return pi, pindex

    raise ValueError(f"Parameter {name} not found.")


def parse_variable_bindings(
    mset: Ncm.MSet, nadd_vals: int, bindings: List[str], option_name: str
) -> dict[str, int]:
    """Parse repeated "name=parameter" (or bare "parameter") bindings.

    A binding without "=" is a shorthand for "parameter=parameter", i.e.
    binding a variable of the same name as the catalog parameter itself.
    """
    var_pindex: dict[str, int] = {}
    for binding in bindings:
        if "=" in binding:
            name, _, param = binding.partition("=")
            name, param = name.strip(), param.strip()
        else:
            name = param = binding.strip()

        if not name or not param:
            raise ValueError(
                f'Invalid {option_name} "{binding}", expected name=parameter '
                f"or a bare parameter name."
            )
        _, pindex = resolve_param(mset, nadd_vals, param)
        var_pindex[name] = pindex

    return var_pindex


def safe_eval_mode(sd1: Ncm.StatsDist1d, n_grid: int = 1000) -> float:
    """Return the mode, using a grid fallback at domain boundaries.

    The grid result is returned when the maximum density is at either
    boundary; otherwise the distribution's mode evaluator is used.

    `Ncm.StatsDist1d.eval_mode()` locates the mode with a GSL bracketing
    minimizer that requires the density maximum to be strictly interior to
    `[xi, xf]`; a skewed or prior-bounded posterior can put its maximum at the
    domain edge. `ncm_stats_dist1d_eval_mode()`
    (numcosmo/ncm/stats/ncm_stats_dist1d.c) detects that case and returns the
    coarse-grid estimate. A native build without that check violates the
    precondition instead and aborts the process with a fatal `g_error`, which
    Python cannot catch, so this function stands as a redundant guard for a
    numcosmo_py running against one. It repeats `eval_mode()`'s coarse linear
    pre-scan at the same resolution and delegates to `eval_mode()` only when
    the scan maximiser is strictly interior and above both endpoints.
    """
    xi = sd1.get_xi()
    xf = sd1.get_xf()
    if xi == xf:
        return xi

    xs = [xi + (xf - xi) / (n_grid - 1) * i for i in range(n_grid)]
    densities = [sd1.eval_p(x) for x in xs]
    argmax = max(range(n_grid), key=lambda i: densities[i])

    if densities[argmax] <= densities[0] or densities[argmax] <= densities[-1]:
        return xs[argmax]

    return sd1.eval_mode()


def median_bounds(
    sd1: Ncm.StatsDist1d, pa: float, median: float
) -> tuple[float, float]:
    """Asymmetric Pa-CI bounds (lower, upper errors) around the median."""
    v = (1.0 - pa) * 0.5
    lo = sd1.eval_inv_pdf(v)
    hi = sd1.eval_inv_pdf_tail(v)
    return median - lo, hi - median


def asymmetric_bounds(
    sd1: Ncm.StatsDist1d, pa: float, center: float
) -> tuple[float, float]:
    """Asymmetric Pa-CI bounds (lower, upper errors) around an arbitrary center.

    Port of mcat_analyze's `_nc_bestfit_error`: handles the case where
    `center` sits close enough to the domain edge that a symmetric-in-CDF
    interval around it would spill past `[xi, xf]`.
    """
    xi = sd1.get_xi()
    xf = sd1.get_xf()
    p_center = sd1.eval_pdf(center)
    pa_2 = pa * 0.5

    if p_center < pa_2:
        lo, hi = xi, sd1.eval_inv_pdf(pa)
    elif p_center > 1.0 - pa_2:
        lo, hi = sd1.eval_inv_pdf_tail(pa), xf
    else:
        lo = sd1.eval_inv_pdf(p_center - pa_2)
        hi = sd1.eval_inv_pdf(p_center + pa_2)

    return center - lo, hi - center


def stat_center_and_bounds(
    stat: DerivedStat, sd1: Ncm.StatsDist1d, bestfit_center: float
) -> tuple[float, List[tuple[float, float]]]:
    """Compute the (center, [(lo, hi) per SIGMA_LEVELS]) reported for `stat`.

    `bestfit_center` is the quantity evaluated at the catalog's best-fit row;
    it is only used when `stat` is `DerivedStat.BESTFIT`, but is required
    unconditionally to keep this a single, branch-free entry point for
    callers that report several statistics for the same distribution.
    """
    bound_fn: Callable[[Ncm.StatsDist1d, float, float], tuple[float, float]]

    if stat == DerivedStat.BESTFIT:
        center = bestfit_center
        bound_fn = asymmetric_bounds
    elif stat == DerivedStat.MODE:
        center = safe_eval_mode(sd1)
        bound_fn = asymmetric_bounds
    else:
        center = sd1.eval_inv_pdf(0.5)
        bound_fn = median_bounds

    bounds = [bound_fn(sd1, pa, center) for pa in SIGMA_LEVELS]
    return center, bounds

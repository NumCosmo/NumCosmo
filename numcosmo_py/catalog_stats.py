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

Utilities to resolve catalog parameters/add-values by name or numeric index,
parse "name=parameter" variable bindings, and compute asymmetric
confidence-interval bounds around a median, mode, or arbitrary center. These
are plain, CLI-agnostic helpers shared by the `catalog` CLI commands and by
the corner-plot derived-quantity support in `numcosmo_py.plotting.derived`.
"""

from __future__ import annotations

import enum
from typing import List, Optional

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
    """Parse repeated "name=parameter" strings into {name: column index}."""
    var_pindex: dict[str, int] = {}
    for binding in bindings:
        name, _, param = binding.partition("=")
        name, param = name.strip(), param.strip()
        if not name or not param:
            raise ValueError(
                f'Invalid {option_name} "{binding}", expected name=parameter.'
            )
        _, pindex = resolve_param(mset, nadd_vals, param)
        var_pindex[name] = pindex

    return var_pindex


def safe_eval_mode(sd1: Ncm.StatsDist1d, n_grid: int = 1000) -> float:
    """Mode of `sd1`, without risking `Ncm.StatsDist1d.eval_mode()`'s crash.

    `eval_mode()` locates the mode with a GSL bracketing minimizer that
    requires the density's maximum to be strictly interior to `[xi, xf]`
    (i.e. above the density at both endpoints). For a skewed or
    prior-bounded posterior whose true maximum sits at (or beyond,
    numerically, right at) the domain edge, that precondition fails and
    NumCosmo raises a fatal `g_error`, which aborts the whole process -- not
    something a Python `try/except` can catch.

    This replicates `eval_mode()`'s own coarse linear pre-scan (same
    resolution) to detect that edge case first. When the scan's maximiser is
    strictly interior and above both endpoints, it is safe to delegate to
    the precise `eval_mode()`. Otherwise, the coarse-grid maximiser is
    returned directly instead of risking the crash.
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

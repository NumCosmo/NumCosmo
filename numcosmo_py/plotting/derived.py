#
# derived.py
#
# Sun Jul 26 00:00:00 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# derived.py
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

"""Derived-quantity support for corner plots.

Computes an extra `CatalogData` column from a mathematical expression of
catalog parameters, so it can be plotted as an additional corner-plot
dimension alongside the catalog's own parameters.
"""

from __future__ import annotations

from typing import List, Optional

import numpy as np

from .. import Ncm
from ..catalog_stats import parse_variable_bindings
from ..safe_eval import compile_expr, SafeExprError
from . import CatalogData


def add_derived_column(
    cd: CatalogData,
    mset: Ncm.MSet,
    nadd_vals: int,
    mcat: Ncm.MSetCatalog,
    variable: List[str],
    expr: str,
    symbol: Optional[str],
    name: str,
) -> CatalogData:
    """Append a derived quantity, computed from `expr`, as a new column of `cd`.

    The underlying parameters bound by `variable` must be present among
    `cd.params_names` (i.e. not removed via --include/--exclude), since the
    values are read straight from the already thinned/burnt-in `cd.rows`
    rather than re-scanned from `mcat`.
    """
    var_pindex = parse_variable_bindings(
        mset, nadd_vals, variable, "--derived-variable"
    )

    try:
        evaluate = compile_expr(expr, var_pindex.keys())
    except SafeExprError as exc:
        raise ValueError(str(exc)) from exc

    col_index: dict[str, int] = {}
    for var_name, pindex in var_pindex.items():
        col_name = mcat.col_name(pindex)
        if col_name not in cd.params_names:
            raise ValueError(
                f'Parameter "{col_name}" (bound to "{var_name}") is not part of '
                f"this catalog's selection; adjust --include/--exclude."
            )
        col_index[var_name] = cd.params_names.index(col_name)

    values = np.array(
        [
            evaluate({var_name: row[idx] for var_name, idx in col_index.items()})
            for row in cd.rows
        ]
    )

    bestfit_value = None
    if cd.bestfit is not None:
        bestfit_value = evaluate(
            {var_name: cd.bestfit[idx] for var_name, idx in col_index.items()}
        )

    return cd.add_derived(
        name=name,
        symbol=symbol if symbol is not None else expr,
        values=values,
        bestfit_value=bestfit_value,
    )

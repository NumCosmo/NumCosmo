#
# getdist.py
#
# Wed Feb 8 10:00:00 2023
# Copyright  2023  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# getdist.py
# Copyright (C) 2023 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

"""NumCosmoPy getdist utilities."""

from typing import List, Tuple
import warnings

import re
import numpy as np
from getdist import MCSamples

from numcosmo_py import Ncm


def mcat_to_mcsamples(
    mcat: Ncm.MSetCatalog,
    name: str,
    asinh_transform: Tuple[int, ...] = (),
    burnin: int = 0,
    thin: int = 1,
    collapse: bool = False,
) -> Tuple[MCSamples, np.ndarray, np.ndarray]:
    """Converts a Ncm.MSetCatalog to a getdist.MCSamples object."""

    nchains: int = mcat.nchains()
    max_time: int = mcat.max_time()

    if burnin % nchains != 0:
        warnings.warn(
            f"burnin ({burnin}) is not a multiple of nchains ({nchains}). "
            f"burnin will be rounded down."
        )
    burnin_steps = burnin // nchains
    burnin = burnin_steps * nchains
    if burnin_steps >= max_time:
        raise ValueError("burnin is greater than the number of steps.")

    assert thin >= 1
    rows: np.ndarray = np.array(
        [
            mcat.peek_row(i * nchains + j).dup_array()
            for i in range(0, max_time, thin)
            for j in range(nchains)
        ]
    )

    params: List[str] = [mcat.col_symb(i) for i in range(mcat.ncols())]
    m2lnL: int = mcat.get_m2lnp_var()  # pylint:disable=invalid-name
    posterior: np.ndarray = 0.5 * rows[:, m2lnL]

    rows = np.delete(rows, m2lnL, 1)
    params = list(np.delete(params, m2lnL, 0))
    names = [re.sub("[^A-Za-z0-9_]", "", param) for param in params]

    weights = None
    if mcat.weighted():
        # Original index is nadd_vals - 1,
        # but since we removed m2lnL it is now nadd_vals - 2
        weight_index = mcat.nadd_vals() - 2
        assert weight_index >= 0
        weights = rows[:, weight_index]
        rows = np.delete(rows, weight_index, 1)
        params = list(np.delete(params, weight_index, 0))
        names = list(np.delete(names, weight_index, 0))

    if len(asinh_transform) > 0:
        rows[:, asinh_transform] = np.arcsinh(rows[:, asinh_transform])
        for i in asinh_transform:
            params[i] = f"\\mathrm{{sinh}}^{{-1}}({params[i]})"
            names[i] = f"asinh_{names[i]}"

    if not collapse:
        split_chains = np.array([rows[(burnin + n) :: nchains] for n in range(nchains)])
        split_posterior = np.array(
            [posterior[(burnin + n) :: nchains] for n in range(nchains)]
        )
    else:
        split_chains = rows[burnin::]
        split_posterior = posterior[burnin::]

    split_weights = None
    if weights is not None:
        if not collapse:
            split_weights = np.array(
                [weights[(burnin + n) :: nchains] for n in range(nchains)]
            )
        else:
            split_weights = weights[burnin::]

    mcsample = MCSamples(
        samples=split_chains,
        loglikes=split_posterior,
        names=names,
        labels=params,
        label=name,
        weights=split_weights,
    )
    return mcsample, rows, posterior

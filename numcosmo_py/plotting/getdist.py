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

from typing import List
import re
import numpy as np

from getdist import MCSamples
from numcosmo_py import Ncm


def mcat_to_mcsamples(
    mcat: Ncm.MSetCatalog, name: str, asinh_transform: List[int] = []
) -> MCSamples:
    """Converts a Ncm.MSetCatalog to a getdist.MCSamples object."""

    nchains: int = mcat.nchains()
    rows: np.ndarray = np.array(
        [mcat.peek_row(i).dup_array() for i in range(0, mcat.len())]
    )
    params: List[str] = [mcat.col_symb(i) for i in range(mcat.ncols())]
    m2lnL: int = mcat.get_m2lnp_var()  # pylint:disable=invalid-name
    posterior: np.ndarray = 0.5 * rows[:, m2lnL]

    rows = np.delete(rows, m2lnL, 1)
    params = list(np.delete(params, m2lnL, 0))
    names = [re.sub("[^A-Za-z0-9_]", "", param) for param in params]

    if len(asinh_transform) > 0:
        rows[:, asinh_transform] = np.arcsinh(rows[:, asinh_transform])
        for i in asinh_transform:
            params[i] = f"\\mathrm{{sinh}}^{{-1}}({params[i]})"
            names[i] = f"asinh_{names[i]}"

    split_chains = np.array([rows[n::nchains] for n in range(nchains)])
    split_posterior = np.array([posterior[n::nchains] for n in range(nchains)])

    mcsample = MCSamples(
        samples=split_chains,
        loglikes=split_posterior,
        names=names,
        labels=params,
        label=name,
    )

    return mcsample, rows, posterior

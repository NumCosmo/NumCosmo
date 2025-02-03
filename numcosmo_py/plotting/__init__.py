#
# __init__.py
#
# Wed Feb 8 10:00:00 2023
# Copyright  2023  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# __init__.py
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

"""NumCosmoPy plotting utilities."""

import warnings
import dataclasses

import numpy as np
import numpy.typing as npt
from getdist import MCSamples

from numcosmo_py import Ncm


@dataclasses.dataclass(frozen=True, kw_only=True)
class CatalogData:
    """Data for a catalog."""

    name: str
    nchains: int
    posterior: np.ndarray
    weights: np.ndarray | None
    rows: np.ndarray
    params_names: list[str]
    params_symbols: list[str]

    def __post_init__(self):
        """Post initialization.

        Checks that the data is valid.
        """
        assert len(self.rows) == len(self.posterior)
        if self.weights is not None:
            assert len(self.rows) == len(self.weights)

        assert len(self.params_names) == len(self.params_symbols)
        assert len(self.params_names) == self.rows.shape[1]

    def asinh_transform(self, indices: npt.NDArray[np.int64]) -> None:
        """Apply an asinh transformation to the catalog data."""
        assert len(indices) > 0
        assert len(indices) <= self.rows.shape[1]

        self.rows[:, indices] = np.arcsinh(self.rows[:, indices])
        for i in indices:
            self.params_symbols[i] = (
                f"\\mathrm{{sinh}}^{{-1}}({self.params_symbols[i]})"
            )
            self.params_names[i] = f"asinh_{self.params_names[i]}"

    def to_mcsamples(self, collapse: bool = False) -> MCSamples:
        """Convert the catalog data to a getdist.MCSamples object.

        :param collapse: If True, collapse the samples into a single chain.
        """
        if collapse:
            return MCSamples(
                samples=self.rows,
                loglikes=self.posterior,
                names=self.params_names,
                labels=self.params_symbols,
                label=self.name,
                weights=self.weights,
            )

        assert self.rows.shape[0] % self.nchains == 0
        ncols = self.rows.shape[1]

        return MCSamples(
            samples=self.rows.reshape(-1, self.nchains, ncols).transpose(1, 0, 2),
            loglikes=self.posterior.reshape(-1, self.nchains).T,
            names=self.params_names,
            labels=self.params_symbols,
            label=self.name,
            weights=(
                self.weights.reshape(-1, self.nchains).T
                if self.weights is not None
                else None
            ),
        )


def mcat_to_catalog_data(
    mcat: Ncm.MSetCatalog,
    name: str,
    burnin: int = 0,
    thin: int = 1,
    indices: npt.NDArray[np.int64] | None = None,
) -> CatalogData:
    """Convert a Ncm.MSetCatalog to a set of numpy arrays."""
    nchains: int = mcat.nchains()
    max_time: int = mcat.max_time()

    if burnin % nchains != 0:
        warnings.warn(
            f"burnin ({burnin}) is not a multiple of nchains ({nchains}). "
            f"burnin will be rounded down."
        )
    burnin_steps = burnin // nchains
    if burnin_steps >= max_time:
        raise ValueError("burnin is greater than the number of steps.")

    assert thin >= 1
    # We are thinning by thin, but in such a way that the last sample is included in the
    # final set of samples, that is why we start at burnin_steps + (max_time -
    # burnin_steps - 1) % thin, instead of burnin_steps.
    rows: np.ndarray = np.array(
        [
            mcat.peek_row(i * nchains + j).dup_array()
            for i in range(
                burnin_steps + (max_time - burnin_steps - 1) % thin, max_time, thin
            )
            for j in range(nchains)
        ]
    )

    if indices is not None:
        indices_array = np.array(indices)
    else:
        indices_array = np.arange(mcat.ncols())

    # Get the -2 log likelihood column
    m2lnL: int = mcat.get_m2lnp_var()  # pylint:disable=invalid-name
    posterior: np.ndarray = 0.5 * rows[:, m2lnL]
    indices_array = indices_array[indices_array != m2lnL]

    # Get the weights column
    weights = None
    if mcat.weighted():
        # Original index is nadd_vals - 1,
        # but since we removed m2lnL it is now nadd_vals - 2
        weight_index = mcat.nadd_vals() - 2
        assert weight_index >= 0
        weights = rows[:, weight_index]
        indices_array = indices_array[indices_array != weight_index]

    rows = rows[:, indices_array]
    param_symbols: list[str] = list(mcat.col_symb(i) for i in indices_array)
    param_names: list[str] = list(mcat.col_name(i) for i in indices_array)

    return CatalogData(
        name=name,
        nchains=nchains,
        posterior=posterior,
        weights=weights,
        rows=rows,
        params_names=param_names,
        params_symbols=param_symbols,
    )

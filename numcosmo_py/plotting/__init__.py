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

import os
import warnings
import dataclasses

import numpy as np
import numpy.typing as npt
import matplotlib.figure
import matplotlib.artist
import matplotlib.pyplot as plt

from numcosmo_py import Ncm
from .tools import set_rc_params_article

original_display = os.environ.get("DISPLAY", "")
if original_display is not None:
    os.environ["DISPLAY"] = ":0"
# pylint: disable=wrong-import-position
# pylint: disable=wrong-import-order
import getdist  # noqa: E402
import getdist.plots  # noqa: E402

# pylint: enable=wrong-import-order
# pylint: enable=wrong-import-position

if original_display is not None:
    os.environ["DISPLAY"] = original_display


@dataclasses.dataclass(frozen=True, kw_only=True)
class CatalogData:
    """Data for a catalog."""

    name: str
    nchains: int
    posterior: np.ndarray
    weights: np.ndarray | None
    rows: np.ndarray
    bestfit: np.ndarray | None
    params_names: list[str]
    params_symbols: list[str]

    def __post_init__(self):
        """Post initialization.

        Checks that the data is valid.
        """
        assert len(self.rows) == len(self.posterior)
        assert self.weights is None or (len(self.rows) == len(self.weights))
        assert len(self.params_names) == len(self.params_symbols)
        assert len(self.params_names) == self.rows.shape[1]
        assert self.bestfit is None or (len(self.bestfit) == self.rows.shape[1])

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

    def to_mcsamples(self, collapse: bool = False) -> getdist.MCSamples:
        """Convert the catalog data to a getdist.MCSamples object.

        :param collapse: If True, collapse the samples into a single chain.
        """
        if collapse:
            return getdist.MCSamples(
                samples=self.rows,
                loglikes=self.posterior,
                names=self.params_names,
                labels=self.params_symbols,
                label=self.name,
                weights=self.weights,
            )

        assert self.rows.shape[0] % self.nchains == 0
        ncols = self.rows.shape[1]

        return getdist.MCSamples(
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
            f"Burnin ({burnin}) is not a multiple of nchains ({nchains}). "
            f"Burnin will be rounded down."
        )
    burnin_steps = burnin // nchains
    if burnin_steps >= max_time:
        raise ValueError("Burnin is greater than the number of steps.")

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
        # weight_index is nadd_vals - 1,
        weight_index = mcat.nadd_vals() - 1
        assert weight_index >= 0
        weights = rows[:, weight_index]
        indices_array = indices_array[indices_array != weight_index]

    rows = rows[:, indices_array]
    param_symbols: list[str] = list(mcat.col_symb(i) for i in indices_array)
    param_names: list[str] = list(mcat.col_name(i) for i in indices_array)

    bestfit = np.array(mcat.get_bestfit_row().dup_array())[indices_array]

    return CatalogData(
        name=name,
        nchains=nchains,
        posterior=posterior,
        weights=weights,
        rows=rows,
        bestfit=bestfit,
        params_names=param_names,
        params_symbols=param_symbols,
    )


def _set_rasterized(element: matplotlib.artist.Artist):
    """Set rasterized to True for all elements in the figure."""
    if hasattr(element, "set_rasterized"):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            element.set_rasterized(True)
    if hasattr(element, "get_children"):
        for child in element.get_children():
            _set_rasterized(child)


def plot_mcsamples(
    mcsamples: list[getdist.MCSamples],
    markers: npt.NDArray[np.float64] | None = None,
    title_limit: int | None = None,
) -> matplotlib.figure.Figure:
    """Plot MCSamples using Matplotlib's object-oriented interface."""
    set_rc_params_article(ncol=2)
    g = getdist.plots.get_subplot_plotter(
        width_inch=plt.rcParams["figure.figsize"][0], rc_sizes=True
    )
    g.settings.linewidth = 0.01
    g.triangle_plot(
        mcsamples,
        shaded=False,
        filled=True,
        markers=markers,
        title_limit=title_limit,
    )

    _set_rasterized(g.fig)

    return g.fig

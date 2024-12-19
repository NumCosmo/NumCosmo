#
# catalog.py
#
# Min Apr 3 11:20:00 2023
# Copyright  2023  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# catalog.py
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

"""A simplified interface for the sampling catalogs."""

from typing import Optional
import numpy as np
import numpy.typing as npt

from getdist import MCSamples

from numcosmo_py import Ncm
from numcosmo_py.helper import npa_to_seq
from .model import build_mset
from ..plotting.getdist import mcat_to_mcsamples


class Catalog:
    """A simplified interface for the sampling catalogs."""

    def __init__(
        self,
        mcat: Optional[Ncm.MSetCatalog] = None,
        ndim: Optional[int] = None,
        nwalkers: Optional[int] = None,
        run_type: Optional[str] = None,
        weighted: bool = False,
    ):
        """Initialize the catalog."""
        if mcat is not None:
            self._catalog = mcat
            if ndim is not None or nwalkers is not None:
                raise ValueError("ndim and nwalkers must be None if mcat is not None")
            mset = mcat.peek_mset()
            if mcat.weighted() != weighted:
                raise ValueError("weighted must be the same as in mcat")
        elif ndim is not None and nwalkers is not None:
            mset = build_mset(ndim)

            self._catalog = Ncm.MSetCatalog.new_array(
                mset, 1, nwalkers, weighted, ["m2lnL"], [r"-2\ln(L)"]
            )

            self._catalog.set_m2lnp_var(0)
        else:
            raise ValueError("Either mcat or ndim and nwalkers must be provided.")

        self.ndim: int = mset.fparams_len()
        if run_type is not None:
            self._catalog.set_run_type(run_type)

    def add_samples(
        self,
        sample: npt.NDArray[np.float64],
        interweaved: bool = True,
    ):
        """Add a new sample to the catalog."""
        ncols = self._catalog.ncols()

        if len(sample.shape) == 1:
            if sample.shape[0] != ncols:
                raise ValueError("sample shape does not match catalog")

            self._catalog.add_from_vector(Ncm.Vector.new_array(npa_to_seq(sample)))

        elif len(sample.shape) == 2:
            if sample.shape[1] != ncols:
                raise ValueError("sample shape does not match catalog")
            if interweaved:
                for point in sample:
                    self._catalog.add_from_vector(Ncm.Vector.new_array(point))
            else:
                nwalkers = self._catalog.nchains()

                if sample.shape[0] % nwalkers != 0:
                    raise ValueError("sample shape does not match catalog")

                nsteps = sample.shape[0] // nwalkers
                for i in range(nsteps):
                    for j in range(nwalkers):
                        self._catalog.add_from_vector(
                            Ncm.Vector.new_array(sample[i + j * nsteps])
                        )

    def add_points_m2lnp(
        self,
        points: np.ndarray,
        m2lnp: np.ndarray,
        interweaved: bool = True,
        weights: Optional[np.ndarray] = None,
    ):
        """Add a new point to the catalog."""
        if self._catalog.weighted() and weights is None:
            raise ValueError("weights must be provided for a weighted catalog")

        if len(points.shape) != 2:
            raise ValueError("points must be a 2D array")

        if len(m2lnp.shape) != 1:
            raise ValueError("m2lnp must be a 1D array")

        if points.shape[0] != m2lnp.shape[0]:
            raise ValueError("points and m2lnp must have the same size")

        sample = np.insert(points, 0, m2lnp, axis=1)

        if weights is not None:
            if len(weights.shape) != 1:
                raise ValueError("weights must be a 1D array")

            if points.shape[0] != weights.shape[0]:
                raise ValueError("points and weights must have the same size")

            sample = np.insert(sample, 1, weights, axis=1)

        self.add_samples(sample, interweaved=interweaved)

    def print_status(self):
        """Print information about the MSet catalog.

        This includes the number of chains, the number of additional values,
        the number of parameters, the number of samples, the number of
        effective samples, the autocorrelation time, the maximum
        autocorrelation time, the maximum effective sample size, the
        maximum effective sample size time, the Heidelberger-Welch
        statistic, the maximum Heidelberger-Welch statistic.
        """
        mcat = self._catalog
        mcat.estimate_autocorrelation_tau(False)
        mcat.log_current_stats()

    def print_info(self, *, ntests=100):
        """Print information about the MSet catalog.

        This includes the number of chains, the number of additional values,
        the number of parameters, the number of samples, the number of
        effective samples, the autocorrelation time, the maximum
        autocorrelation time, the maximum effective sample size, the
        maximum effective sample size time, the Heidelberger-Welch
        statistic, the maximum Heidelberger-Welch statistic.
        """
        mcat = self._catalog
        mset = mcat.peek_mset()
        mcat.estimate_autocorrelation_tau(False)

        Ncm.cfg_msg_sepa()
        print(f"# Catalog run type: `{mcat.get_run_type()}`")
        print(f"# Catalog size:      {mcat.len()}.")
        print(f"# Catalog n-chains:  {mcat.nchains()}.")
        print(f"# Catalog nadd-vals: {mcat.nadd_vals()}.")
        print(f"# Catalog weighted:  {mcat.weighted()}")
        mcat.log_current_chain_stats()
        mcat.calc_max_ess_time(ntests, Ncm.FitRunMsgs.FULL)
        mcat.calc_heidel_diag(ntests, 0.0, Ncm.FitRunMsgs.FULL)
        mset.pretty_log()
        mcat.log_full_covar()
        mcat.log_current_stats()

    def trim(self, nsteps: int, thin: Optional[int] = None):
        """Trim the catalog to remove the first `nsteps` steps.

        Optionally, thin the catalog by a factor of `thin`.

        :param nsteps: number of steps to remove
        :param thin: thinning factor
        """
        if thin is None:
            thin = 1
        self._catalog.trim(nsteps, thin)

    def get_mean(self) -> np.ndarray:
        """Get the mean of the catalog."""
        return np.array(self._catalog.get_mean().dup_array())

    def get_covar(self) -> np.ndarray:
        """Get the covariance matrix of the catalog."""
        return np.array(self._catalog.get_covar().dup_array()).reshape(
            self.ndim, self.ndim
        )

    def get_mcsamples(self, name: str) -> MCSamples:
        """Get the MCSamples object from the catalog."""
        mcsample, _, _ = mcat_to_mcsamples(self._catalog, name)
        return mcsample

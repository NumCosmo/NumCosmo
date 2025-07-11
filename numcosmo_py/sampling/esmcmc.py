#
# esmcmc.py
#
# Wed Feb 8 10:00:00 2023
# Copyright  2023  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# esmcmc.py
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

"""Create a new ensemble sampler object."""

from typing import Optional, Union
import warnings
from enum import StrEnum, auto
from numcosmo_py import Ncm
from numcosmo_py.interpolation.stats_dist import (
    InterpolationMethod,
    InterpolationKernel,
)


class WalkerTypes(StrEnum):
    """Possible walkers for ensemble samplers."""

    APES = auto()
    STRETCH = auto()


def create_esmcmc(
    likelihood: Ncm.Likelihood,
    mset: Ncm.MSet,
    prefix: str,
    *,
    verbose: bool = True,
    fit_first: bool = False,
    robust: bool = False,
    use_apes_interpolation: bool = True,
    use_apes_threads: Optional[bool] = None,
    sampler: WalkerTypes = WalkerTypes.APES,
    interpolation_method: InterpolationMethod = InterpolationMethod.VKDE,
    interpolation_kernel: InterpolationKernel = InterpolationKernel.CAUCHY,
    nwalkers: int = 320,
    nthreads: int = 4,
    over_smooth: float = 1.0,
    local_fraction: Optional[float] = None,
    init_sampling_scale: float = 1.0e-1,
    start_mcat: Optional[Ncm.MSetCatalog] = None,
):
    """Create a new ensemble sampler object."""
    # New fit object using the likelihood.
    fit = Ncm.Fit.factory(
        Ncm.FitType.NLOPT,
        "ln-neldermead",
        likelihood,
        mset,
        Ncm.FitGradType.NUMDIFF_FORWARD,
    )

    # Printing fitting informations.
    if verbose:
        fit.log_info()
        message_level = Ncm.FitRunMsgs.SIMPLE
    else:
        message_level = Ncm.FitRunMsgs.NONE

    if fit_first:
        fit.run(message_level)

    init_sampler: Union[Ncm.MSetTransKernCat, Ncm.MSetTransKernGauss]
    if start_mcat is not None:
        init_sampler = Ncm.MSetTransKernCat.new(start_mcat, None)
        init_sampler.set_sampling(Ncm.MSetTransKernCatSampling.CHOOSE)
        init_sampler.set_mset(mset)
        init_sampler.set_prior_from_mset()
    else:
        init_sampler = Ncm.MSetTransKernGauss.new(0)
        init_sampler.set_mset(mset)
        init_sampler.set_prior_from_mset()
        init_sampler.set_cov_from_rescale(init_sampling_scale)

    #
    # Creates the ESMCMC walker object, this object is responsible
    # for moving the walkers in each interation, the stretch move
    # is affine invariant and therefore gives good results even for
    # very correlated parametric space.
    #

    walker: Union[Ncm.FitESMCMCWalkerAPES, Ncm.FitESMCMCWalkerStretch]
    if sampler == WalkerTypes.APES:
        walker = Ncm.FitESMCMCWalkerAPES.new(nwalkers, mset.fparams_len())
        # Sets the calibrated over-smoothing factor.
        walker.set_over_smooth(over_smooth)
        if local_fraction is not None:
            walker.set_local_frac(local_fraction)
        if robust:
            walker.set_cov_robust()
        walker.use_interp(use_apes_interpolation)
        walker.set_method(interpolation_method.genum)
        walker.set_k_type(interpolation_kernel.genum)
        if use_apes_threads is None:
            use_apes_threads = nwalkers >= 1000

        walker.set_use_threads(use_apes_threads)
        if use_apes_threads and nwalkers < 1000:
            warnings.warn(
                "Using threads with less than 1000 walkers can degrade performance."
            )
        elif not use_apes_threads and nwalkers >= 1000:
            warnings.warn(
                "Using threads with more than 1000 walkers can improve performance."
            )

    elif sampler == WalkerTypes.STRETCH:
        walker = Ncm.FitESMCMCWalkerStretch.new(nwalkers, mset.fparams_len())
    else:
        raise ValueError(f"Unknown sampler {sampler}")

    # Initialize the ESMCMC object using the objects above.
    esmcmc = Ncm.FitESMCMC.new(fit, nwalkers, init_sampler, walker, message_level)

    # Setting the number of threads to use.
    esmcmc.set_nthreads(nthreads)
    # Setting the file name to save the chains.
    filename = f"{prefix}_mcmc_catalog_{sampler}_{nwalkers}.fits"
    esmcmc.set_data_file(filename)

    return esmcmc


def mcat_print_info(mcat, *, ntests=100):
    """Print information about the MSet catalog.

    This includes the number of chains, the number of additional values,
    the number of parameters, the number of samples, the number of
    effective samples, the autocorrelation time, the maximum
    autocorrelation time, the maximum effective sample size, the
    maximum effective sample size time, the Heidelberger-Welch
    statistic, the maximum Heidelberger-Welch statistic.
    """
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

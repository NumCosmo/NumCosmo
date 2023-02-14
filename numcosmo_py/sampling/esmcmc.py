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

from enum import Enum
from numcosmo_py import Ncm, GEnum


class WalkerTypes(str, Enum):
    """Possible walkers for ensemble samplers."""

    APES = "apes"
    STRECTCH = "stretch"


class InterpolationMethod(GEnum):
    """Possible interpolation methods for APES walkers."""
    # pylint: disable=no-member
    KDE = Ncm.FitESMCMCWalkerAPESMethod.KDE
    VKDE = Ncm.FitESMCMCWalkerAPESMethod.VKDE


class InterpolationKernel(GEnum):
    """Possible interpolation kernels for APES walkers."""
    # pylint: disable=no-member
    CAUCHY = Ncm.FitESMCMCWalkerAPESKType.CAUCHY
    ST3 = Ncm.FitESMCMCWalkerAPESKType.ST3
    GAUSS = Ncm.FitESMCMCWalkerAPESKType.GAUSS


def create_esmcmc(
    likelihood: Ncm.Likelihood,
    mset: Ncm.MSet,
    prefix: str,
    *,
    verbose: bool = True,
    fit_first: bool = False,
    robust: bool = False,
    use_apes_interpolation: bool = True,
    sampler: WalkerTypes = WalkerTypes.APES,
    interpolation_method: InterpolationMethod = InterpolationMethod.VKDE,
    interpolation_kernel: InterpolationKernel = InterpolationKernel.CAUCHY,
    nwalkers: int = 320,
    nthreads: int = 4,
    over_smooth: float = 1.0,
    init_sampling_scale: float = 1.0e-1,
):
    """Create a new ensemble sampler object."""

    # New fit object using the likelihood.
    fit = Ncm.Fit.new(
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

    if sampler == WalkerTypes.APES:
        walker = Ncm.FitESMCMCWalkerAPES.new(nwalkers, mset.fparams_len())
        # Sets the calibrated over-smoothing factor.
        walker.set_over_smooth(over_smooth)
        if robust:
            walker.set_cov_robust()
        walker.use_interp(use_apes_interpolation)
        walker.set_method(interpolation_method.genum)
        walker.set_k_type(interpolation_kernel.genum)

    elif sampler == WalkerTypes.STRECTCH:
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

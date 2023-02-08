#
# rosenbrock_mcmc.py
#
# Wed Feb 8 10:00:00 2023
# Copyright  2023  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# rosenbrock_mcmc.py
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

"""Example of using the Rosenbrock function to test the MCMC sampler.
"""

import gi

gi.require_version("NumCosmo", "1.0")
gi.require_version("NumCosmoMath", "1.0")

# pylint:disable-next=wrong-import-position
from gi.repository import NumCosmoMath as Ncm  # noqa: E402


def run_rosenbrock_mcmc(
    sampler: str = "apes",
    nwalkers: int = 320,
    ssize: int = 5000000,
    verbose: bool = False,
    nthreads: int = 4,
):
    """Runs the Rosenbrock MCMC example."""

    # New Rosenbrock model object.
    mrb = Ncm.ModelRosenbrock()

    # New Model set object including slm with parameters
    # set as free.
    mset = Ncm.MSet.empty_new()
    mset.set(mrb)
    mset.param_set_all_ftype(Ncm.ParamType.FREE)
    mset.prepare_fparam_map()

    # Create a new data set object using the Rosenbrock function.
    drb = Ncm.DataRosenbrock.new()

    # Create a data set object and append the data.
    dset = Ncm.Dataset.new()
    dset.append_data(drb)

    # New likelihood object using dset.
    likelihood = Ncm.Likelihood.new(dset)

    # Creating a Fit object of type NLOPT using the fitting algorithm
    # ln-neldermead to fit the Modelset mset using the Likelihood lh
    # and using a numerical differentiation algorithm (NUMDIFF_FORWARD)
    # to obtain the gradient (if needed).
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

    #
    # Setting single thread calculation.
    #
    Ncm.func_eval_set_max_threads(nthreads)
    Ncm.func_eval_log_pool_stats()

    #
    # New Gaussian prior to provide the initial points for the chain.
    # It was created with size 0 (number of parameters), but once
    # initialized with mset the correct size is assigned.
    #
    # The initial sampler will use a diagonal covariance with the
    # diagonal terms being the parameters scale set by each model.
    #
    init_sampler = Ncm.MSetTransKernGauss.new(0)
    init_sampler.set_mset(mset)
    init_sampler.set_prior_from_mset()
    init_sampler.set_cov_from_rescale(100.0)

    #
    # Creates the ESMCMC walker object, this object is responsible
    # for moving the walkers in each interation, the stretch move
    # is affine invariant and therefore gives good results even for
    # very correlated parametric space.
    #

    if sampler == "apes":
        walker = Ncm.FitESMCMCWalkerAPES.new(nwalkers, mset.fparams_len())
        # Sets the calibrated over-smoothing factor.
        walker.set_over_smooth(0.2)
    elif sampler == "stretch":
        walker = Ncm.FitESMCMCWalkerStretch.new(nwalkers, mset.fparams_len())
    else:
        raise ValueError(f"Unknown sampler {sampler}")

    # Initialize the ESMCMC object using the objects above.
    if verbose:
        message_level = Ncm.FitRunMsgs.SIMPLE
    else:
        message_level = Ncm.FitRunMsgs.NONE

    esmcmc = Ncm.FitESMCMC.new(fit, nwalkers, init_sampler, walker, message_level)

    # Setting the number of threads to use.
    esmcmc.set_nthreads(nthreads)
    # Setting the file name to save the chains.
    esmcmc.set_data_file(f"rosenbrock_chains_{sampler}_{nwalkers}.fits")

    # Running the esmcmc.
    esmcmc.start_run()
    esmcmc.run(ssize // nwalkers)
    esmcmc.end_run()

    # Calculates the parameter means and covariance and set it into
    # the fit object and then print.
    if verbose:
        esmcmc.mean_covar()
        fit.log_covar()

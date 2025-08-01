#
# gaussmix2d.py
#
# Wed Feb 8 10:00:00 2023
# Copyright  2023  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# gaussmix2d.py
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

"""Example of using the Gauss Mixture 2d function to test the MCMC sampler."""

from numcosmo_py import Ncm
from numcosmo_py.sampling.esmcmc import (
    create_esmcmc,
    WalkerTypes,
    InterpolationMethod,
    InterpolationKernel,
)


def run_gaussmix2d_mcmc(
    *,
    ssize: int = 5000000,
    verbose: bool = True,
    fit_first: bool = True,
    robust: bool = False,
    use_apes_interpolation: bool = True,
    use_apes_threads: bool = True,
    sampler: WalkerTypes = WalkerTypes.APES,
    interpolation_method: InterpolationMethod = InterpolationMethod.VKDE,
    interpolation_kernel: InterpolationKernel = InterpolationKernel.CAUCHY,
    nwalkers: int = 320,
    nthreads: int = 4,
    over_smooth: float = 1.1,
    init_sampling_scale: float = 1.0e-1,
) -> str:
    """Runs the Gauss Mixture 2d MCMC example."""

    mrb = Ncm.ModelRosenbrock()
    mset = Ncm.MSet.empty_new()
    mset.set(mrb)
    mset.param_set_all_ftype(Ncm.ParamType.FREE)
    mset.prepare_fparam_map()

    drb = Ncm.DataGaussMix2D.new()
    dset = Ncm.Dataset.new()
    dset.append_data(drb)
    likelihood = Ncm.Likelihood.new(dset)

    esmcmc = create_esmcmc(
        likelihood,
        mset,
        "gaussmix2d",
        verbose=verbose,
        fit_first=fit_first,
        robust=robust,
        use_apes_interpolation=use_apes_interpolation,
        use_apes_threads=use_apes_threads,
        sampler=sampler,
        interpolation_method=interpolation_method,
        interpolation_kernel=interpolation_kernel,
        nwalkers=nwalkers,
        nthreads=nthreads,
        over_smooth=over_smooth,
        init_sampling_scale=init_sampling_scale,
    )

    esmcmc.start_run()
    esmcmc.run(ssize // nwalkers)
    esmcmc.end_run()

    fit = esmcmc.peek_fit()
    filename = esmcmc.peek_catalog().peek_filename()
    if verbose:
        esmcmc.mean_covar()
        fit.log_covar()

    return filename

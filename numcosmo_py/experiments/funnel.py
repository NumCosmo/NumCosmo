#
# funnel.py
#
# Wed Feb 8 10:00:00 2023
# Copyright  2023  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# funnel.py
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

from pathlib import Path
from typing import Optional

from numcosmo_py import Ncm
from numcosmo_py.sampling.esmcmc import (
    create_esmcmc,
    WalkerTypes,
    InterpolationKernel,
    InterpolationMethod,
)


def run_funnel_mcmc(
    *,
    dim: int = 10,
    ssize: int = 5000000,
    verbose: bool = True,
    fit_first: bool = False,
    robust: bool = False,
    use_apes_interpolation: bool = True,
    use_apes_threads: Optional[bool] = None,
    sampler: WalkerTypes = WalkerTypes.APES,
    interpolation_method: InterpolationMethod = InterpolationMethod.VKDE,
    interpolation_kernel: InterpolationKernel = InterpolationKernel.CAUCHY,
    nwalkers: int = 3000,
    nthreads: int = 4,
    over_smooth: float = 0.2,
    local_fraction: Optional[float] = None,
    init_sampling_scale: float = 1.0e2,
    start_catalog: Optional[Path] = None,
) -> str:
    """Runs the Funnel MCMC example."""

    mfu = Ncm.ModelFunnel.new(dim - 1)

    mset = Ncm.MSet.empty_new()
    mset.set(mfu)
    mset.param_set_all_ftype(Ncm.ParamType.FREE)
    mset.prepare_fparam_map()

    dfu = Ncm.DataFunnel.new()
    dset = Ncm.Dataset.new()
    dset.append_data(dfu)
    likelihood = Ncm.Likelihood.new(dset)

    start_mcat = None
    if start_catalog is not None:
        start_mcat = Ncm.MSetCatalog.new_from_file_ro(start_catalog.as_posix(), 0)

    esmcmc = create_esmcmc(
        likelihood,
        mset,
        f"funnel{dim}",
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
        local_fraction=local_fraction,
        init_sampling_scale=init_sampling_scale,
        start_mcat=start_mcat,
    )

    # Running the esmcmc.
    esmcmc.start_run()
    esmcmc.run(ssize // nwalkers)
    esmcmc.end_run()

    fit = esmcmc.peek_fit()
    filename = esmcmc.peek_catalog().peek_filename()
    if verbose:
        esmcmc.mean_covar()
        fit.log_covar()

    return filename

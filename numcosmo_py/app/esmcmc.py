#
# esmcmc.py
#
# Wed Feb 14 18:59:34 2024
# Copyright  2024  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# esmcmc.py
# Copyright (C) 2024 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

"""NumCosmo APP subcommand run the ESMCMC algorithm."""

import dataclasses
from enum import Enum
from pathlib import Path
from typing import Optional, Annotated, Union

import typer

from .. import Ncm
from ..interpolation.stats_dist import (
    InterpolationKernel,
    InterpolationMethod,
)
from .run_fit import RunCommonOptions


class IniSampler(str, Enum):
    """Initial sampler to use for the MCMC."""

    GAUSS_MSET = "gauss-mset"
    GAUSS_COV = "gauss-cov"
    FROM_CATALOG = "from-catalog"


class Parallezation(str, Enum):
    """Parallel sampler to use for the MCMC."""

    NONE = "none"
    MPI = "mpi"
    THREADS = "threads"


@dataclasses.dataclass
class RunMCMC(RunCommonOptions):
    """Computes the MCMC using APES."""

    nwalkers: Annotated[
        int,
        typer.Option(
            help="Number of walkers to use, 0 means automatic.",
            min=0,
        ),
    ] = 0

    nsamples: Annotated[
        int,
        typer.Option(
            help="Number of samples to compute.",
            min=1,
        ),
    ] = 10

    robust: Annotated[
        bool,
        typer.Option(
            help="Use robust covariance estimation. (Experimental)",
        ),
    ] = False

    interpolation_method: Annotated[
        InterpolationMethod,
        typer.Option(
            help="Interpolation method to use.",
        ),
    ] = InterpolationMethod.VKDE

    interpolation_kernel: Annotated[
        InterpolationKernel,
        typer.Option(
            help="Interpolation kernel to use.",
        ),
    ] = InterpolationKernel.CAUCHY

    over_smooth: Annotated[
        float,
        typer.Option(
            help="Over-smoothing parameter to use.",
            min=1.0e-2,
        ),
    ] = 1.0

    local_fraction: Annotated[
        Optional[float],
        typer.Option(
            help="Local fraction to use.",
            min=0.02,
        ),
    ] = None

    use_interpolation: Annotated[
        bool,
        typer.Option(
            help="Use interpolation to compute the weights of the APES approximation.",
        ),
    ] = True

    parallel: Annotated[
        Parallezation,
        typer.Option(
            help=(
                "Parallelization to use. Python likelihoods are not compatible with "
                "multi-threading."
            ),
        ),
    ] = Parallezation.NONE

    nthreads: Annotated[
        int,
        typer.Option(
            help="Number of threads to use when using multi-threading.",
            min=2,
        ),
    ] = 4

    initial_points_sampler: Annotated[
        Optional[IniSampler],
        typer.Option(
            help=(
                "Sampler to use for the initial points. "
                "Gaussian based samplers use a gaussian approximation of the "
                "posterior. The covariance matrix can be computed using the "
                "scales in the model-set or using the covariance matrix of the last "
                "fit. The catalog sampler uses a catalog of points to sample the "
                "initial points. The covariance can be specified using the "
                "--initial-sampler-covar option or using the covariance matrix in the "
                "product file."
            ),
        ),
    ] = IniSampler.GAUSS_MSET

    initial_sampler_rescale: Annotated[
        float,
        typer.Option(
            help=(
                "Rescale factor for the Gaussian based initial samplers. "
                "The rescale factor is applied to the covariance matrix or to the "
                "scales in the model-set."
            ),
            min=0.01,
        ),
    ] = 1.0

    initial_sampler_covar: Annotated[
        Optional[Path],
        typer.Option(
            help=(
                "Path to the covariance matrix file to use for the initial points "
                "sampler."
            ),
        ),
    ] = None

    initial_catalog: Annotated[
        Optional[Path],
        typer.Option(
            help="Path to the catalog file to use for the initial points sampler.",
        ),
    ] = None

    initial_catalog_burnin: Annotated[
        int,
        typer.Option(
            help="Number of points to discard from the catalog.",
            min=0,
        ),
    ] = 0

    exploration: Annotated[
        int,
        typer.Option(
            help=(
                "Number of samples to use for the exploration phase. The exploration "
                " phase should be discarded from the final samples."
            ),
        ),
    ] = 0

    def __post_init__(self) -> None:
        super().__post_init__()
        self.fit.log_info()

        fparams_len = self.mset.fparams_len()

        if self.nwalkers == 0:
            self.nwalkers = 100 * (fparams_len + 1)

        init_sampler: Union[Ncm.MSetTransKernGauss, Ncm.MSetTransKernCat]
        if self.initial_points_sampler == IniSampler.GAUSS_MSET:
            init_sampler = Ncm.MSetTransKernGauss.new(0)
            init_sampler.set_mset(self.mset)
            init_sampler.set_prior_from_mset()
            init_sampler.set_cov_from_rescale(self.initial_sampler_rescale)
        elif self.initial_points_sampler == IniSampler.GAUSS_COV:
            init_sampler = Ncm.MSetTransKernGauss.new(0)
            init_sampler.set_mset(self.mset)
            init_sampler.set_prior_from_mset()

            if self.initial_sampler_covar is None:
                if self.output_dict.get("covariance") is None:
                    raise RuntimeError(
                        "No covariance file given and the product file "
                        "does not contain a covariance matrix."
                    )
                cov = self.output_dict.get("covariance")
                if cov is None or not isinstance(cov, Ncm.Matrix):
                    raise RuntimeError(
                        "Covariance matrix cannot be found in the product file."
                    )
            else:
                ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
                cov_dict = ser.dict_str_from_yaml_file(
                    self.initial_sampler_covar.absolute().as_posix()
                )
                cov = cov_dict.get("covariance")
                if cov is None or not isinstance(cov, Ncm.Matrix):
                    raise RuntimeError(
                        f"Covariance matrix not found in file"
                        f"{self.initial_sampler_covar}"
                    )

            cov.scale(self.initial_sampler_rescale)
            init_sampler.set_cov(cov)
        elif self.initial_points_sampler == IniSampler.FROM_CATALOG:
            if self.initial_catalog is None:
                raise RuntimeError("No catalog file given.")
            if not self.initial_catalog.exists():
                raise RuntimeError(f"Catalog file {self.initial_catalog} not found.")

            mcat = Ncm.MSetCatalog.new_from_file_ro(
                self.initial_catalog.absolute().as_posix(),
                self.initial_catalog_burnin,
            )
            init_sampler = Ncm.MSetTransKernCat.new(mcat, None)

            init_sampler.set_sampling(Ncm.MSetTransKernCatSampling.CHOOSE)
            init_sampler.set_mset(self.mset)
            init_sampler.set_prior_from_mset()
        else:
            raise RuntimeError(
                f"Invalid initial points sampler {self.initial_points_sampler}."
            )

        apes_walker = Ncm.FitESMCMCWalkerAPES.new(self.nwalkers, fparams_len)
        apes_walker.set_over_smooth(self.over_smooth)
        if self.local_fraction is not None:
            apes_walker.set_local_frac(self.local_fraction)

        apes_walker.use_interp(self.use_interpolation)
        apes_walker.set_method(self.interpolation_method.genum)
        apes_walker.set_k_type(self.interpolation_kernel.genum)

        if self.parallel == Parallezation.THREADS.value:
            apes_walker.set_use_threads(True)
        elif self.parallel == Parallezation.MPI.value:
            if Ncm.cfg_mpi_nslaves() == 0:
                raise RuntimeError(
                    "MPI parallelization requested but MPI is not initialized."
                )
        else:
            apes_walker.set_use_threads(False)

        if self.exploration > 0:
            apes_walker.set_exploration(self.exploration)

        if self.functions is not None:
            esmcmc: Ncm.FitESMCMC = Ncm.FitESMCMC.new_funcs_array(
                self.fit,
                self.nwalkers,
                init_sampler,
                apes_walker,
                self.run_messages.genum,
                self.functions,
            )
        else:
            esmcmc = Ncm.FitESMCMC.new(
                self.fit,
                self.nwalkers,
                init_sampler,
                apes_walker,
                self.run_messages.genum,
            )

        if self.parallel == Parallezation.THREADS.value:
            esmcmc.set_nthreads(self.nthreads)

        if self.output is not None:
            esmcmc.set_data_file(
                self.output.absolute().with_suffix(".mcmc.fits").as_posix()
            )

        esmcmc.start_run()
        esmcmc.run(self.nsamples)
        esmcmc.end_run()

        self.end_experiment()

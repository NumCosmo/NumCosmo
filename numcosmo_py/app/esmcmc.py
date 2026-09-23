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
from enum import StrEnum, auto
from pathlib import Path
from typing import Optional, Annotated, Union

import typer

from .. import Ncm
from ..interpolation.stats_dist import (
    CrossValidationMethod,
    InterpolationKernel,
    InterpolationMethod,
)
from .run_fit import RunCommonOptions


class IniSampler(StrEnum):
    """Initial sampler to use for the MCMC."""

    GAUSS_MSET = "gauss-mset"
    GAUSS_COV = "gauss-cov"
    FROM_CATALOG = "from-catalog"


class Parallelization(StrEnum):
    """Parallel sampler to use for the MCMC."""

    NONE = auto()
    MPI = auto()
    THREADS = auto()


@dataclasses.dataclass(kw_only=True)
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

    cv_method: Annotated[
        CrossValidationMethod,
        typer.Option(
            help=(
                "Cross-validation used to choose the over-smoothing factor. NONE uses "
                "the value given by --over-smooth. SPLIT_NOFIT builds the approximation "
                "from a fraction --split-fraction of each block and chooses the "
                "over-smoothing factor on the remaining points of that block."
            ),
        ),
    ] = CrossValidationMethod.NONE

    split_fraction: Annotated[
        Optional[float],
        typer.Option(
            help="Fraction of each block used as kernel centres.",
            min=0.02,
            max=1.0,
        ),
    ] = None

    auto_kernel: Annotated[
        bool,
        typer.Option(
            help=(
                "Choose the interpolation kernel together with the over-smoothing "
                "factor, by the same out-of-sample objective. Requires --cv-method "
                "split-nofit and overrides --interpolation-kernel."
            ),
        ),
    ] = False

    center_shrink: Annotated[
        bool,
        typer.Option(
            help=(
                "Shrink the kernel centres toward the ensemble mean so that the APES "
                "approximation has the same covariance as the ensemble. Requires a "
                "kernel with a finite covariance, so not the Cauchy one."
            ),
        ),
    ] = False

    defensive_frac: Annotated[
        float,
        typer.Option(
            min=0.0,
            max=1.0,
            help=(
                "Weight of a wide Student-t component mixed into the APES proposal, "
                "centered on the ensemble mean with --defensive-scale times its "
                "covariance and --defensive-nu degrees of freedom. Keeps the proposal "
                "density positive where the kernels leave holes. Zero disables it."
            ),
        ),
    ] = 0.0

    defensive_scale: Annotated[
        float,
        typer.Option(min=1.0e-2, help="Covariance factor of the wide component."),
    ] = 4.0

    defensive_nu: Annotated[
        float,
        typer.Option(min=1.0, help="Degrees of freedom of the wide component."),
    ] = 3.0

    vkde_points_per_dim: Annotated[
        float,
        typer.Option(
            min=0.0,
            help=(
                "VKDE only: nearest neighbours per dimension for each local covariance, "
                "k = min(n, c d). Replaces --local-fraction when positive; a small "
                "ensemble then gives the KDE limit and a large one keeps the kernels local."
            ),
        ),
    ] = 0.0

    uniform_weights: Annotated[
        bool,
        typer.Option(
            help=(
                "Keep uniform kernel weights instead of the NNLS fit. The bandwidth and "
                "kernel cross-validation still run, including the methods that need the "
                "ensemble's -2lnL."
            ),
        ),
    ] = False

    parallel: Annotated[
        Parallelization,
        typer.Option(
            help=(
                "Parallelization to use. Python likelihoods are not compatible with "
                "multi-threading."
            ),
        ),
    ] = Parallelization.NONE

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
                "Length cap of the APES exploration phase, in iterations. With "
                "--exploration-qratio-floor 0 the phase accepts by the posterior ratio "
                "alone for exactly this many iterations; with a positive floor it ends "
                "earlier, after --exploration-stop-after consecutive iterations in which "
                "no acceptance was clipped. The phase runs only when the chain starts from "
                "its initial ensemble, and the catalog's markovian-id records where it ended."
            ),
            min=0,
        ),
    ] = 0

    exploration_qratio_floor: Annotated[
        float,
        typer.Option(
            help=(
                "Floor of the proposal-density ratio q(x)/q(x') in the acceptance during "
                "the exploration phase, so walkers where the proposal has almost no mass "
                "can leave. 0 disables the clip; 1 uses the posterior ratio alone for "
                "every blocked move."
            ),
            min=0.0,
            max=1.0,
        ),
    ] = 0.0

    exploration_stop_after: Annotated[
        int,
        typer.Option(
            help=(
                "Consecutive iterations in which no acceptance was clipped, after which "
                "the exploration phase ends and the clip is disarmed for the rest of the run."
            ),
            min=1,
        ),
    ] = 10

    skip_check: Annotated[
        bool,
        typer.Option(
            help="Skip the check of the last ensemble when continuing a run.",
        ),
    ] = True

    seed: Annotated[
        Optional[int],
        typer.Option(
            min=0,
            help=(
                "Seed of the random number generator used by the sampler, including "
                "the initial points. If not given, a seed is drawn and printed in the "
                "log."
            ),
        ),
    ] = None

    def __post_init__(self) -> None:
        """Run the ESMCMC algorithm."""
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
        # After the kernel, so that an incompatible pair is caught immediately.
        apes_walker.set_center_shrink(self.center_shrink)
        apes_walker.set_defensive_frac(self.defensive_frac)
        apes_walker.set_defensive_scale(self.defensive_scale)
        apes_walker.set_defensive_nu(self.defensive_nu)
        apes_walker.set_vkde_points_per_dim(self.vkde_points_per_dim)
        apes_walker.set_uniform_weights(self.uniform_weights)
        apes_walker.set_cv_type(self.cv_method.genum)
        apes_walker.set_auto_kernel(self.auto_kernel)
        if self.split_fraction is not None:
            apes_walker.set_split_frac(self.split_fraction)

        if self.parallel == Parallelization.THREADS.value:
            apes_walker.set_use_threads(True)
        elif self.parallel == Parallelization.MPI.value:
            if Ncm.cfg_mpi_nslaves() == 0:
                raise RuntimeError(
                    "MPI parallelization requested but MPI is not initialized."
                )
        else:
            apes_walker.set_use_threads(False)

        apes_walker.set_exploration(self.exploration)
        apes_walker.set_exploration_qratio_floor(self.exploration_qratio_floor)
        apes_walker.set_exploration_stop_after(self.exploration_stop_after)

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

        esmcmc.set_use_threads(self.parallel == Parallelization.THREADS.value)

        if self.output is not None:
            esmcmc.set_data_file(
                self.output.absolute().with_suffix(".mcmc.fits").as_posix()
            )

        esmcmc.set_skip_check(self.skip_check)

        if self.seed is not None:
            esmcmc.set_rng(Ncm.RNG.seeded_new(None, self.seed))

        esmcmc.start_run()
        esmcmc.run(self.nsamples)
        esmcmc.end_run()

        self.end_experiment()

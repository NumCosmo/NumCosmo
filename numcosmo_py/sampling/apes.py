#
# apes.py
#
# Sun Apr 2 15:40:00 2023
# Copyright  2023  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# apes.py
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

"""A simplified interface for the ensemble sampler NcmFitESMCMC using APES."""


from typing import Optional, Callable, Union, Tuple, List
import numpy as np

from numcosmo_py import Ncm
from numcosmo_py.interpolation.stats_dist import (
    InterpolationMethod,
    InterpolationKernel,
)

from .model import NcmModelGeneric, get_generic_model
from .catalog import Catalog


class APES:
    """A simplified interface for the ensemble sampler NcmFitESMCMC using APES."""

    def __init__(
        self,
        *,
        nwalkers: int,
        ndim: Optional[int],
        model: Optional[Ncm.Model],
        log_prob: Callable[[Union[np.ndarray, List[float]], Tuple], float],
        args: Tuple = (),
        verbose: bool = False,
        robust: bool = False,
        use_interpolation: bool = True,
        interpolation_method: InterpolationMethod = InterpolationMethod.VKDE,
        interpolation_kernel: InterpolationKernel = InterpolationKernel.CAUCHY,
        over_smooth: float = 0.2,
        local_fraction: Optional[float] = None,
    ):
        """Create a new APES sampler object."""

        self.log_prob = log_prob
        self.args = args
        self.verbose: bool = verbose
        self.nwalkers = nwalkers

        if ndim is None and model is None:
            raise ValueError("Either ndim or model must be provided.")
        if ndim is not None and model is not None:
            raise ValueError("Only one of ndim or model must be provided.")

        if ndim is not None:
            # pylint:disable-next=invalid-name
            model = NcmModelGeneric(theta_length=ndim)

        if model is None:
            raise ValueError("Either ndim or model must be provided.")

        self.mset = Ncm.MSet()
        self.mset.set(model)
        self.mset.param_set_all_ftype(Ncm.ParamType.FREE)
        self.mset.prepare_fparam_map()
        pself = self

        class APESLikehood(Ncm.Data):
            """APES likelihood."""

            def __init__(self) -> None:
                super().__init__()
                self.set_init(True)
                self.theta: Optional[np.ndarray] = None

            def do_get_length(self) -> int:  # pylint: disable-msg=arguments-differ
                return 100

            def do_begin(self) -> None:  # pylint: disable-msg=arguments-differ
                pass

            # pylint: disable-next=arguments-differ
            def do_prepare(self, mset: Ncm.MSet):
                model = get_generic_model(mset)
                self.theta = np.array(model.orig_params_peek_vector().dup_array())

            def do_m2lnL_val(self, _):  # pylint: disable-msg=arguments-differ
                return -2.0 * pself.log_prob(self.theta, pself.args)

        apes_likelihood = APESLikehood()

        dset = Ncm.Dataset()
        dset.append_data(apes_likelihood)

        likelihood = Ncm.Likelihood(dataset=dset)

        self.fit = Ncm.Fit.new(
            Ncm.FitType.NLOPT,
            "ln-neldermead",
            likelihood,
            self.mset,
            Ncm.FitGradType.NUMDIFF_FORWARD,
        )

        walker = Ncm.FitESMCMCWalkerAPES.new(nwalkers, self.mset.fparams_len())
        walker.set_over_smooth(over_smooth)
        if local_fraction is not None:
            walker.set_local_frac(local_fraction)
        if robust:
            walker.set_cov_robust()
        walker.use_interp(use_interpolation)
        walker.set_method(interpolation_method.genum)
        walker.set_k_type(interpolation_kernel.genum)

        init_sampler = Ncm.MSetTransKernGauss.new(0)
        init_sampler.set_mset(self.mset)
        init_sampler.set_prior_from_mset()
        init_sampler.set_cov_from_rescale(1.0)

        if verbose:
            msg_level = Ncm.FitRunMsgs.SIMPLE
        else:
            msg_level = Ncm.FitRunMsgs.NONE

        self.esmcmc = Ncm.FitESMCMC.new(
            self.fit, nwalkers, init_sampler, walker, msg_level
        )

        self.esmcmc.set_nthreads(1)

    def run_mcmc(self, initial_sample: np.ndarray, niter: int) -> None:
        """Run the sampler for niter iterations."""

        if initial_sample.shape != (self.nwalkers, self.mset.fparams_len()):
            raise ValueError("Initial sample has wrong shape.")

        mcat = self.esmcmc.peek_catalog()

        self.esmcmc.start_run()

        if mcat.len() == 0:
            for point in initial_sample:
                point_vector = Ncm.Vector.new_array(point)
                self.mset.fparams_set_vector(point_vector)
                m2lnL = self.fit.m2lnL_val()  # pylint:disable=invalid-name
                mcat.add_from_vector_array(point_vector, [m2lnL])
            assert mcat.len() == self.nwalkers

        assert mcat.len() >= self.nwalkers

        self.esmcmc.run(niter)
        self.esmcmc.end_run()

    def get_catalog(self) -> Catalog:
        """Get the catalog of the sampler."""

        return Catalog(mcat=self.esmcmc.peek_catalog())

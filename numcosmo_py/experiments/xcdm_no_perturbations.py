#
# xcdm_no_perturbations.py
#
# Sun Mar 8 12:12:00 2023
# Copyright  2023  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# xcdm_no_perturbations.py
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

""" Experiments the XCDM model using likelihoods that do not
depend on perturbations. """

import os
from typing import Optional

from numcosmo_py import Ncm, Nc
from numcosmo_py.sampling.esmcmc import (
    create_esmcmc,
    WalkerTypes,
    InterpolationMethod,
    InterpolationKernel,
)

from numcosmo_py.datasets.no_perturbations import (
    SNIaID,
    BAOID,
    HID,
    gen_snia_likelihood,
    gen_bao_likelihood,
    gen_h_likelihood,
)


def run_xcdm_nopert_mcmc(
    *,
    ssize: int = 5000000,
    verbose: bool = True,
    fit_first: bool = False,
    robust: bool = False,
    use_apes_interpolation: bool = True,
    sampler: WalkerTypes = WalkerTypes.APES,
    interpolation_method: InterpolationMethod = InterpolationMethod.VKDE,
    interpolation_kernel: InterpolationKernel = InterpolationKernel.CAUCHY,
    nwalkers: int = 2000,
    nthreads: int = 1,
    over_smooth: float = 0.2,
    init_sampling_scale: float = 1.0e0,
    flat: bool = False,
    use_neutrino: bool = False,
    z_f: float = 3.0,
    snia_id: Optional[SNIaID] = SNIaID.COV_PANTHEON_PLUS_SH0ES_SYS_STAT,
    bao_id: Optional[BAOID] = BAOID.ALL_COMBINED_JAN_2023,
    h_id: Optional[HID] = HID.ALL_COMBINED_JAN_2023,
) -> str:
    """Runs the XCDM model with no perturbations MCMC."""

    model_str = "xcdm"
    data_str = "data_local"

    if flat:
        model_str = f"{model_str}_flat"

    filename_base = f"{model_str}_{data_str}"
    progress_file = f"{filename_base}_progress.mset"

    ser = Ncm.Serialize.new(0)
    dset = Ncm.Dataset.new()

    if os.path.exists(progress_file):
        mset = Ncm.MSet.load(progress_file, ser)
    else:
        mset = Ncm.MSet.empty_new()
        
        if use_neutrino:
            cosmo = Nc.HICosmo.new_from_name (Nc.HICosmo, "NcHICosmoDEXcdm{'massnu-length':<1>}")
        else:
            cosmo = Nc.HICosmo.new_from_name(Nc.HICosmo, "NcHICosmoDEXcdm")

        cosmo.omega_x2omega_k()

        cosmo.param_set_by_name("H0", 70.0)
        cosmo.param_set_by_name("Omegab", 0.05)
        cosmo.param_set_by_name("Omegac", 0.25)
        cosmo.param_set_by_name("Omegak", 0.00)
        
        if use_neutrino:
            cosmo.orig_param_set(Nc.HICosmoDESParams.ENNU, 2.0328)    
            id = cosmo.vparam_index(Nc.HICosmoDEVParams.M, 0)
            cosmo.param_set_ftype(id, Ncm.ParamType.FREE)

        cosmo.props.H0_fit = True
        cosmo.props.Omegac_fit = True
        cosmo.props.Omegax_fit = not flat
        cosmo.props.w_fit = True

        mset.set(cosmo)

    dist = Nc.Distance(zf=z_f)

    if snia_id is not None:
        gen_snia_likelihood(dset, mset, dist, snia_id)

    if bao_id is not None:
        gen_bao_likelihood(dset, mset, dist, bao_id)

    if h_id is not None:
        gen_h_likelihood(dset, mset, h_id)

    mset.prepare_fparam_map()
    likelihood = Ncm.Likelihood.new(dset)

    esmcmc = create_esmcmc(
        likelihood,
        mset,
        "xcdm-no-perturbations",
        verbose=verbose,
        fit_first=fit_first,
        robust=robust,
        use_apes_interpolation=use_apes_interpolation,
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

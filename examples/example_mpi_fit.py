#!/usr/bin/env python
#
# example_mpi_fit.py
#
# Mon May 22 16:00:00 2023
# Copyright  2023  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# example_mpi_fit.py
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

"""Example testing MPI fit."""

import sys
import timeit

from numcosmo_py import Nc, Ncm

#
#  Initializing the library objects, this must be called before
#  any other library function.
#
sys.argv = Ncm.cfg_init_full(sys.argv)


def test_mpi_fit(test_union: bool = False) -> None:
    """Example testing MPI fit."""

    cosmo = Nc.HICosmoDEXcdm()

    cosmo.param_set_by_name("H0", 70.0)
    cosmo.param_set_by_name("Omegab", 0.05)
    cosmo.param_set_by_name("Omegac", 0.25)
    cosmo.param_set_by_name("Omegax", 0.70)
    cosmo.param_set_by_name("Tgamma0", 2.72)
    cosmo.param_set_by_name("w", -1.0)

    #
    #  Creating a new Modelset and set cosmo as the HICosmo model to be used.
    #
    mset = Ncm.MSet()
    mset.set(cosmo)

    #
    #  Setting parameters Omega_c and w to be fitted.
    #
    cosmo.props.Omegac_fit = True
    cosmo.props.w_fit = True

    #
    #  Creating a new Distance object optimized to redshift 2.
    #
    dist = Nc.Distance(zf=2.0)

    #
    #  Creating a new Data object from distance modulus catalogs.
    #
    if test_union:
        data_snia = Nc.DataDistMu.new_from_id(dist, Nc.DataSNIAId.SIMPLE_UNION2_1)
    else:
        snia = Nc.SNIADistCov.new(dist, 4)
        mset.set(snia)
        mset.param_set_mid_ftype(Nc.SNIADistCov.id(), Ncm.ParamType.FIXED)
        mset.param_set_ftype(
            Nc.SNIADistCov.id(), Nc.SNIADistCovSParams.ALPHA, Ncm.ParamType.FREE
        )
        mset.param_set_ftype(
            Nc.SNIADistCov.id(), Nc.SNIADistCovSParams.BETA, Ncm.ParamType.FREE
        )
        # mset.param_set_ftype (Nc.SNIADistCov.id (), Nc.SNIADistCovSParams.M1,
        # Ncm.ParamType.FREE)
        # mset.param_set_ftype (Nc.SNIADistCov.id (), Nc.SNIADistCovSParams.M2,
        # Ncm.ParamType.FREE)
        data_snia = Nc.DataSNIACov.new_from_cat_id(
            Nc.DataSNIAId.COV_JLA_SNLS3_SDSS_SYS_STAT_CMPL, False
        )

    #
    #  Creating a new Dataset and add snia to it.
    #
    dset = Ncm.Dataset()
    dset.append_data(data_snia)

    #
    #  Creating a Likelihood from the Dataset.
    #
    lh = Ncm.Likelihood(dataset=dset)

    #
    # Creating a Fit object of type NLOPT using the fitting algorithm ln-neldermead to
    # fit the Modelset mset using the Likelihood lh and using a numerical
    # differentiation algorithm (NUMDIFF_FORWARD) to obtain the gradient (if needed).
    #
    fit = Ncm.Fit.factory(
        Ncm.FitType.NLOPT, "ln-neldermead", lh, mset, Ncm.FitGradType.NUMDIFF_FORWARD
    )

    mj = Ncm.MPIJobFit.new(fit, None)

    ser = Ncm.Serialize.new(Ncm.SerializeOpt.NONE)
    mj.init_all_slaves(ser)

    # time.sleep (10000)

    rng = Ncm.RNG.new()
    rng.set_random_seed(False)

    param_a = []
    m2lnL_a = []
    for _ in range(24):
        w_i = rng.gaussian_gen(-1.0, 0.05)
        Oc_i = rng.gaussian_gen(0.25, 0.01)
        a_i = rng.gaussian_gen(0.141, 0.01)
        b_i = rng.gaussian_gen(3.101, 0.1)
        param_a.append(Ncm.Vector.new_array([Oc_i, w_i, a_i, b_i]))
        m2lnL_a.append(Ncm.Vector.new(1 + 4))

    # mj.run_array (param_a, m2lnL_a)

    t1 = timeit.Timer(
        "mj.run_array (param_a, m2lnL_a)", "from __main__ import mj, param_a, m2lnL_a"
    )
    print("Timming: ", t1.timeit(1))

    for param, m2lnL in zip(param_a, m2lnL_a):
        fit.params_set_vector(param)
        print(
            f"{m2lnL.get(0): -22.15g} {m2lnL.get(1): -22.15g} {m2lnL.get(2): -22.15g} "
            f"{m2lnL.get(3): -22.15g} {m2lnL.get(4): -22.15g}"
        )

    mj.free_all_slaves()


if __name__ == "__main__":
    test_mpi_fit()

#!/usr/bin/env python
#
# test_data_cluster_wl_factor_mc_bias.py
#
# Copyright (C) 2026 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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
# with this program. If not, see <http://www.gnu.org/licenses/>.

"""Tier 2 acceptance gate: from-model MC mass-recovery bias reproduction.

A previous Monte-Carlo study run through the legacy CLI
(``numcosmo_py/experiments/cluster_wl.py`` + ``NcDataClusterWL`` +
``NcGalaxySDShapeHSMGaussGlobal``) found a robust, highly significant
cluster-mass *underestimate* -- worse at lower mass, roughly -20% flat from
3e14-3e15 Msun and -28% to -32% at 1e14 Msun -- isolated to the
intrinsic-shape-noise marginalization approximation itself (the same
approximation ``NcGalaxyShapeFactorVarAdd`` reproduces bit-identically, see
``test_galaxy_shape_factor_var_add.py``).

This test reproduces that finding through the *new* orchestrator
(``NcDataClusterWLFactor`` + ``VarAdd``) via ``ncm_data_resample`` +
``Ncm.Fit``, at 1e14 Msun (where the bias was largest and most clearly
resolved in the original study). It is a from-model MC acceptance gate, not a
tight numerical parity check like Tier 1 -- it only needs to confirm the new
pipeline exhibits the *same known pathology* (a significant, negative bias
of comparable order of magnitude) before ``VarAdd`` is ever swapped for
``Quad``/``Laplace``/``FixedQuad`` to test whether that swap closes the gap.

A single 3000-galaxy cluster gives a *very* weak individual mass constraint
(the per-realization scatter in the fitted log10MDelta is order 0.5-1 dex),
and a meaningful minority of realizations (order 10%) have such a weak shear
signal that the fit runs away to ``log10MDelta``'s own lower bound (10.0)
rather than converging to a genuine (if noisy) estimate -- the same class of
methodological pitfall as the fit-bound-clamping bug found while first
measuring the legacy bias (see the plan/session notes): the bound itself
isn't too narrow (14 +/- a few dex of headroom), the realization is just
uninformative. These boundary-clamped realizations are excluded from the
statistics below as non-converged, exactly as one would exclude a
non-converged MCMC chain; with them excluded, this reproduction landed at
-20% to -30% (measured from a 300-realization run while writing this test),
squarely matching the legacy-measured range.

Runs a CI-sized 150-realization MC (a few minutes) rather than the original
study's full 3k-20k-galaxy sweep across masses, or the 300-realization run
used to first establish this test's acceptance band -- enough for the
(excluding-boundary-hits) mean bias to clear the significance bar reliably,
not to characterize the bias precisely.
"""

import numpy as np
import pytest
from numcosmo_py import Nc, Ncm

Ncm.cfg_init()

TRUE_LOG10M = 14.0
Z_CL = 0.2
R_MIN, R_MAX = 0.3, 5.0
N_GAL = 3000
SIGMA_INT = 0.28
STD_NOISE = 0.03
FIELD_HALF = 0.09  # degrees
N_REALIZATIONS = 150
SEED0 = 1000

# The original CLI-driven study found roughly -20% to -32% bias in this mass
# range. Accept anything significantly negative and within an order of
# magnitude of that -- a loose band, since this is a different survey
# realization/geometry, not a tight parity check (see module docstring).
MIN_ABS_REL_BIAS = 0.08
MAX_ABS_REL_BIAS = 0.50
MIN_SIGNIFICANCE_SIGMA = 1.5
MAX_NON_CONVERGED_FRACTION = 0.30


def _build():
    cosmo = Nc.HICosmoDEXcdm.new()
    dist = Nc.Distance.new(100.0)
    hms = Nc.HaloCMParam.new(Nc.HaloMassSummaryMassDef.MEAN, 200.0)
    dp = Nc.HaloDensityProfileNFW.new(hms)
    hp = Nc.HaloPosition.new(dist)
    smd = Nc.WLSurfaceMassDensity.new(dist)

    hms.param_set_by_name("cDelta", 5.0)
    hms.param_set_by_name("log10MDelta", TRUE_LOG10M)
    hms.param_set_ftype(1, Ncm.ParamType.FREE)  # only mass is fit; cDelta stays fixed
    hp.param_set_by_name("z", Z_CL)
    hp.prepare(cosmo)

    pop_shape = Nc.GalaxyShapePopGauss.new()
    pop_shape.param_set_by_name("sigma", SIGMA_INT)
    pop_z = Nc.GalaxyRedshiftPopLSSTSRD.new_y1_source()
    obs_z = Nc.GalaxyRedshiftObsGauss.new()

    mset = Ncm.MSet.empty_new()
    for model in (cosmo, dp, hp, smd, pop_shape, pop_z, obs_z):
        mset.set(model)
    mset.prepare_fparam_map()

    position_factor = Nc.GalaxyPositionFactorFlat.new(
        -FIELD_HALF, FIELD_HALF, -FIELD_HALF, FIELD_HALF
    )
    redshift_factor = Nc.GalaxyRedshiftFactorComposed.new(0.0, 5.0)
    shape_factor = Nc.GalaxyShapeFactorVarAdd.new(Nc.GalaxyWLObsEllipConv.TRACE_DET)

    pos_data0 = Nc.GalaxyPositionFactorData.new(position_factor, mset)
    z_data0 = Nc.GalaxyRedshiftFactorData.new(redshift_factor, mset)
    s_data0 = Nc.GalaxyShapeFactorData.new(shape_factor, mset, pos_data0, z_data0)
    cols = Nc.GalaxyShapeFactorData.required_columns(s_data0)

    obs = Nc.GalaxyWLObs.new(
        Nc.GalaxyWLObsEllipConv.TRACE_DET, Nc.WLEllipticityFrame.CELESTIAL, N_GAL, cols
    )
    for i in range(N_GAL):
        for c in cols:
            obs.set(c, i, 0.0)
        obs.set("std_noise", i, STD_NOISE)
        obs.set("sigma0", i, 0.05)  # fixed calibration input to gen(), never redrawn

    dcwlf = Nc.DataClusterWLFactor.new(position_factor, redshift_factor, shape_factor)
    dcwlf.set_obs(obs)
    dcwlf.set_cut(R_MIN, R_MAX)
    dcwlf.set_prec(1.0e-6)

    dset = Ncm.Dataset.new()
    dset.append_data(dcwlf)
    likelihood = Ncm.Likelihood.new(dset)
    fit = Ncm.Fit.factory(
        Ncm.FitType.NLOPT,
        "ln-neldermead",
        likelihood,
        mset,
        Ncm.FitGradType.NUMDIFF_FORWARD,
    )
    fit.set_params_reltol(1.0e-8)
    fit.set_m2lnL_reltol(1.0e-9)
    fit.set_maxiter(10000)

    return mset, hms, dcwlf, fit


def test_mass_recovery_bias_matches_known_pathology():
    """From-model MC: NcDataClusterWLFactor + VarAdd reproduces the known
    negative mass bias, within an order of magnitude of the legacy-measured
    -20% to -32% (see module docstring)."""
    mset, hms, dcwlf, fit = _build()
    lower_bound = hms.param_get_lower_bound(1)

    fitted_log10m = np.empty(N_REALIZATIONS)
    for i in range(N_REALIZATIONS):
        hms.param_set_by_name("log10MDelta", TRUE_LOG10M)
        rng = Ncm.RNG.seeded_new(None, SEED0 + i)
        dcwlf.resample(mset, rng)
        fit.run(Ncm.FitRunMsgs.NONE)
        fitted_log10m[i] = hms.param_get(1)

    # Exclude non-converged (boundary-clamped) realizations -- see module
    # docstring: a genuinely weak-signal realization, not a too-narrow bound.
    converged = fitted_log10m[fitted_log10m > lower_bound + 0.01]
    non_converged_fraction = 1.0 - len(converged) / N_REALIZATIONS
    assert non_converged_fraction < MAX_NON_CONVERGED_FRACTION, (
        f"too many non-converged (boundary-clamped) realizations: "
        f"{non_converged_fraction:.0%}"
    )

    mean_log10m = converged.mean()
    sem_log10m = converged.std(ddof=1) / np.sqrt(len(converged))
    mean_bias = 10.0**mean_log10m / 10.0**TRUE_LOG10M - 1.0
    significance = (TRUE_LOG10M - mean_log10m) / sem_log10m

    assert mean_bias < 0.0, f"expected a negative mass bias, got {mean_bias:+.1%}"
    assert significance > MIN_SIGNIFICANCE_SIGMA, (
        f"bias not significant: log10MDelta = {mean_log10m:.3f} +/- "
        f"{sem_log10m:.3f} vs true {TRUE_LOG10M} ({significance:.1f} sigma)"
    )
    assert MIN_ABS_REL_BIAS < -mean_bias < MAX_ABS_REL_BIAS, (
        f"bias {mean_bias:+.1%} outside the order-of-magnitude band matching "
        f"the legacy-measured -20% to -32%"
    )


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

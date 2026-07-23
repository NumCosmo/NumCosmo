#!/usr/bin/env python
#
# test_data_cluster_wl_factor_parity_combinations.py
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

"""Cross-check of the remaining ``NcDataClusterWLFactor`` vs legacy
``NcDataClusterWL`` combinations at the full orchestrator level.

``test_data_cluster_wl_factor_parity.py`` only ever exercises ONE
combination: shape=VarAdd+GaussGlobal x redshift=Composed+Gauss x
position=Flat. Legacy only ever provides an oracle for the *variance-add*
shape approximation (``VarAdd``), in either its Global (single shared
sigma, ``NcGalaxySDShapeHSMGaussGlobal``) or per-galaxy (``e_rms``/
``std_shape`` catalog column, ``NcGalaxySDShapeHSMGauss``) flavour, and for
the redshift scheme in either its Composed (population x observable,
``NcGalaxySDObsRedshiftGauss`` with ``use_true_z``) or Spline (per-galaxy
pre-tabulated p(z), ``NcGalaxySDObsRedshiftPz``) flavour -- position only
ever had one flavour (Flat). This file covers the three legacy-comparable
combinations left untested at the orchestrator level:

* shape=GaussGlobal x redshift=Spline   (redshift axis, new)
* shape=GaussLocal  x redshift=Spline   (both axes, new)
* shape=GaussLocal  x redshift=Composed (shape axis, new)

(GaussGlobal x Composed is the existing file's combination.) The no-legacy-
oracle schemes (Quad/SeriesLensed/FixedQuad/Laplace shape, Beta pop) are out
of scope here by design -- see the cross-check scoping discussion.

FROZEN REFERENCE VALUES: the dual-engine agreement documented above was
proven by running both engines live and is captured, not re-derived, in the
functions below. Values were captured from an actual passing run of this
file's original legacy-comparison code, at git rev ``77313f22``
(2026-07-16), then legacy (``NcDataClusterWL``) construction was removed so
these tests no longer depend on legacy at runtime -- legacy is slated for
deletion in a follow-up PR. Each frozen assertion keeps the tolerance
(``rtol``/``atol``) that the original live comparison used.
"""

import numpy as np
import pytest
from numpy.testing import assert_allclose

from numcosmo_py import Ncm, Nc

Ncm.cfg_init()

Z_CL = 0.2
ZP_MIN, ZP_MAX = 0.0, 5.0
ELLIP_CONV = Nc.GalaxyWLObsEllipConv.TRACE_DET
FRAME = Nc.WLEllipticityFrame.CELESTIAL
SIGMA_INT_GLOBAL = 0.3
PZ_HALF_WIDTH_SIGMAS = 6.0

# (ra, dec, zp, sigma0, e_rms, e1, e2, std_noise, c1, c2, m) -- includes
# nonzero multiplicative (m) and additive (c1, c2) calibration bias on every
# row, since those terms are exactly what a bug in the shear/bias wiring
# would show up in.
_GALAXIES = [
    (0.03, 0.02, 0.60, 0.030, 0.28, 0.05, -0.02, 0.03, 0.005, -0.003, 0.05),
    (-0.10, 0.15, 0.90, 0.040, 0.35, -0.04, 0.01, 0.05, -0.002, 0.004, -0.10),
    (0.05, -0.08, 0.15, 0.020, 0.22, 0.02, 0.03, 0.04, 0.008, -0.006, 0.08),
    (0.08, 0.05, 1.50, 0.050, 0.30, 0.03, -0.01, 0.04, 0.001, -0.001, 0.02),
    (-0.05, -0.03, 0.40, 0.030, 0.18, -0.01, 0.02, 0.035, -0.004, 0.007, -0.03),
]

_LOG10_MDELTA_GRID = [13.5, 14.0, 14.5]


def _build_lens_mset():
    cosmo = Nc.HICosmoDEXcdm.new()
    dist = Nc.Distance.new(100.0)
    hms = Nc.HaloCMParam.new(Nc.HaloMassSummaryMassDef.MEAN, 200.0)
    dp = Nc.HaloDensityProfileNFW.new(hms)
    hp = Nc.HaloPosition.new(dist)
    smd = Nc.WLSurfaceMassDensity.new(dist)

    hp.param_set_by_name("z", Z_CL)
    hp.prepare(cosmo)

    return cosmo, dp, hp, smd, hms


def _pz_spline(zp, sigma0, n=200):
    """A per-galaxy p(z) spline for the Spline/Pz combos; domain scaled to
    sigma0 to avoid a near-delta-function integrand (see
    adaptive-spline-degenerate-integrand-oom)."""
    z_min = max(1.0e-3, zp - PZ_HALF_WIDTH_SIGMAS * sigma0)
    z_max = zp + PZ_HALF_WIDTH_SIGMAS * sigma0
    zs = np.linspace(z_min, z_max, n)
    pz_vals = np.exp(-0.5 * ((zs - zp) / sigma0) ** 2) + 1.0e-6

    xv = Ncm.Vector.new_array(zs.tolist())
    yv = Ncm.Vector.new_array(pz_vals.tolist())
    spline = Ncm.SplineCubicNotaknot.new()
    spline.set(xv, yv, True)
    return spline


def _build_mset(shape_kind, z_kind):
    """One shared mset carrying the new engine's own population/observable
    models, keyed by MAIN id -- exactly the strangler-fig pattern used
    throughout this refactor (see calculators-must-not-hold-models).
    Legacy is no longer built in this file (see module docstring), so no
    legacy-side models are registered here."""
    cosmo, dp, hp, smd, hms = _build_lens_mset()
    mset = Ncm.MSet.empty_new()
    for model in (cosmo, dp, hp, smd):
        mset.set(model)

    if shape_kind == "global":
        pop_shape = Nc.GalaxyShapePopGauss.new()
        pop_shape.param_set_by_name("sigma", SIGMA_INT_GLOBAL)
    elif shape_kind == "local":
        pop_shape = Nc.GalaxyShapePopGaussLocal.new()
    else:
        raise ValueError(shape_kind)

    mset.set(pop_shape)

    if z_kind == "composed":
        pop_z = Nc.GalaxyRedshiftPopLSSTSRD.new_y1_source()
        obs_z = Nc.GalaxyRedshiftObsGauss.new()
        mset.set(pop_z)
        mset.set(obs_z)
    elif z_kind == "spline":
        pass  # NcGalaxyRedshiftFactorSpline needs no mset-registered model
    else:
        raise ValueError(z_kind)

    mset.prepare_fparam_map()
    return mset, hms


def _set_common_columns(obs, i, galaxy):
    ra, dec, _zp, _sigma0, _e_rms, e1, e2, std_noise, c1, c2, m = galaxy
    obs.set("ra", i, ra)
    obs.set("dec", i, dec)
    obs.set("z", i, 0.0)
    obs.set("epsilon_int_1", i, 0.0)
    obs.set("epsilon_int_2", i, 0.0)
    obs.set("epsilon_obs_1", i, e1)
    obs.set("epsilon_obs_2", i, e2)
    obs.set("std_noise", i, std_noise)
    obs.set("c1", i, c1)
    obs.set("c2", i, c2)
    obs.set("m", i, m)


def _build_new_obs(mset, shape_kind, z_kind):
    position_factor = Nc.GalaxyPositionFactorFlat.new(-0.2, 0.2, -0.2, 0.2)
    shape_factor = Nc.GalaxyShapeFactorVarAdd.new(ELLIP_CONV)

    if z_kind == "composed":
        redshift_factor = Nc.GalaxyRedshiftFactorComposed.new(ZP_MIN, ZP_MAX)
    elif z_kind == "spline":
        redshift_factor = Nc.GalaxyRedshiftFactorSpline.new()
    else:
        raise ValueError(z_kind)

    pos_data = Nc.GalaxyPositionFactorData.new(position_factor, mset)
    z_data = Nc.GalaxyRedshiftFactorData.new(redshift_factor, mset)
    s_data = Nc.GalaxyShapeFactorData.new(shape_factor, mset, pos_data, z_data)
    cols = Nc.GalaxyShapeFactorData.required_columns(s_data)

    obs = Nc.GalaxyWLObs.new(ELLIP_CONV, FRAME, len(_GALAXIES), cols)

    for i, galaxy in enumerate(_GALAXIES):
        _, _, zp, sigma0, e_rms, *_ = galaxy
        _set_common_columns(obs, i, galaxy)

        if z_kind == "composed":
            obs.set("zp", i, zp)
            obs.set("sigma0", i, sigma0)
        else:
            obs.set_pz(i, _pz_spline(zp, sigma0))

        if shape_kind == "local":
            obs.set("e_rms", i, e_rms)

    return position_factor, redshift_factor, shape_factor, obs


_COMBOS = [
    ("global", "spline"),
    ("local", "spline"),
    ("local", "composed"),
]
_COMBO_IDS = ["shape=global,z=spline", "shape=local,z=spline", "shape=local,z=composed"]

# Frozen legacy (LNINT) -2lnL values, keyed by (shape_kind, z_kind,
# log10_mdelta) -- see module docstring for provenance.
_M2LNL_PARITY_FROZEN = {
    ("global", "spline", 13.5): 1.303719131183008,
    ("global", "spline", 14.0): 1.3087951336465575,
    ("global", "spline", 14.5): 1.3246741139788694,
    ("local", "spline", 13.5): -1.2931531151570173,
    ("local", "spline", 14.0): -1.2863843595056863,
    ("local", "spline", 14.5): -1.265885887457485,
    ("local", "composed", 13.5): -20.819891846703342,
    ("local", "composed", 14.0): -20.813127069429562,
    ("local", "composed", 14.5): -20.79263171691778,
}


@pytest.mark.parametrize("shape_kind,z_kind", _COMBOS, ids=_COMBO_IDS)
@pytest.mark.parametrize("log10_mdelta", _LOG10_MDELTA_GRID)
def test_m2lnL_parity(shape_kind, z_kind, log10_mdelta):
    """NcDataClusterWLFactor's total -2lnL matches legacy's (LNINT) across
    masses, for each of the three previously-uncovered combinations. Now
    checked against frozen legacy values (see module docstring)."""
    mset, hms = _build_mset(shape_kind, z_kind)
    hms.param_set_by_name("log10MDelta", log10_mdelta)

    position_factor, redshift_factor, shape_factor, new_obs = _build_new_obs(
        mset, shape_kind, z_kind
    )

    dcwlf = Nc.DataClusterWLFactor.new(position_factor, redshift_factor, shape_factor)
    dcwlf.set_obs(new_obs)
    dcwlf.set_prec(1.0e-8)
    dcwlf.set_integ_method(Nc.DataClusterWLIntegMethod.LNINT)

    new_m2lnL = dcwlf.m2lnL_val(mset)

    assert_allclose(
        new_m2lnL, _M2LNL_PARITY_FROZEN[(shape_kind, z_kind, log10_mdelta)], rtol=1.0e-5
    )


# Frozen legacy per-galaxy -2lnP breakdown at log10MDelta=14.0, LNINT, keyed
# by (shape_kind, z_kind) -- see module docstring for provenance.
_GAL_PARITY_FROZEN = {
    ("global", "spline"): [
        0.4161808949411303,
        -0.13477623249471948,
        1.232379132917283,
        -0.6051355310897102,
        0.4001468693725738,
    ],
    ("local", "spline"): [
        0.16982133021347767,
        0.6554924901197818,
        0.03528566593842067,
        -0.557523303218451,
        -1.5894605425589154,
    ],
    ("local", "composed"): [
        -4.89848387332248,
        -3.00789079609498,
        -4.668946253982924,
        -1.4173938224390745,
        -6.820412323590107,
    ],
}


@pytest.mark.parametrize("shape_kind,z_kind", _COMBOS, ids=_COMBO_IDS)
def test_m2lnL_gal_parity(shape_kind, z_kind):
    """Per-galaxy breakdown localizes any mismatch to a single galaxy. Now
    checked against a frozen legacy breakdown (see module docstring)."""
    mset, hms = _build_mset(shape_kind, z_kind)
    hms.param_set_by_name("log10MDelta", 14.0)

    position_factor, redshift_factor, shape_factor, new_obs = _build_new_obs(
        mset, shape_kind, z_kind
    )

    dcwlf = Nc.DataClusterWLFactor.new(position_factor, redshift_factor, shape_factor)
    dcwlf.set_obs(new_obs)
    dcwlf.set_prec(1.0e-8)
    dcwlf.set_integ_method(Nc.DataClusterWLIntegMethod.LNINT)

    n = len(_GALAXIES)
    new_gal = Ncm.Vector.new(n)

    dcwlf.eval_m2lnP_gal(mset, new_gal)

    assert_allclose(
        new_gal.dup_array(), _GAL_PARITY_FROZEN[(shape_kind, z_kind)], rtol=1.0e-5
    )


# Frozen legacy FIXED_NODES (default for both engines) -2lnL values at
# log10MDelta=14.0, keyed by (shape_kind, z_kind).
_FIXED_NODES_PARITY_FROZEN = {
    ("global", "spline"): 1.308794969011216,
    ("local", "spline"): -1.2863845241404908,
    ("local", "composed"): -20.81312706848369,
}


@pytest.mark.parametrize("shape_kind,z_kind", _COMBOS, ids=_COMBO_IDS)
def test_m2lnL_parity_fixed_nodes(shape_kind, z_kind):
    """FIXED_NODES (both engines' default) agrees with legacy's own default
    exactly, for each combination -- checks the orchestrator's fixed-node
    wiring is correct for schemes other than Composed+GaussGlobal too. Now
    checked against frozen legacy values (see module docstring)."""
    mset, hms = _build_mset(shape_kind, z_kind)
    hms.param_set_by_name("log10MDelta", 14.0)

    position_factor, redshift_factor, shape_factor, new_obs = _build_new_obs(
        mset, shape_kind, z_kind
    )

    dcwlf = Nc.DataClusterWLFactor.new(position_factor, redshift_factor, shape_factor)
    dcwlf.set_obs(new_obs)
    dcwlf.set_prec(1.0e-8)
    dcwlf.set_integ_method(Nc.DataClusterWLIntegMethod.FIXED_NODES)

    new_m2lnL = dcwlf.m2lnL_val(mset)

    assert_allclose(
        new_m2lnL, _FIXED_NODES_PARITY_FROZEN[(shape_kind, z_kind)], rtol=1.0e-6
    )


# Frozen legacy resample() output, keyed by (shape_kind, z_kind) then seed
# (see module docstring for provenance). ra/dec (and zp when z=composed) are
# bit-identical (rtol=0, atol=0) between engines; epsilon_obs_1/2 agree to
# rtol=0 for shape=global and rtol=1e-8 for shape=local (see docstring below);
# m2lnL agrees to rtol=1e-5. All are frozen at the values legacy produced.
_RESAMPLE_FROZEN = {
    ("global", "spline"): {
        200: {
            "cols": {
                "ra": [
                    -0.10938103049993514,
                    0.09081182871013882,
                    -0.15151121113449334,
                    0.18835634412243962,
                    -0.032229752186685806,
                ],
                "dec": [
                    -0.002342543371666539,
                    0.16387784329891236,
                    -0.15068294679505903,
                    0.13041990375585624,
                    -0.006571839549873736,
                ],
                "epsilon_obs_1": [
                    0.45648033539105765,
                    0.1663836009567453,
                    -0.14695473345308083,
                    -0.13732948567043893,
                    -0.2500497729429979,
                ],
                "epsilon_obs_2": [
                    -0.11120904619770347,
                    0.2783205429969342,
                    0.17399840334584926,
                    0.6573851710590882,
                    0.07351834336953396,
                ],
            },
            "m2lnL": 11.018394045298528,
        },
        201: {
            "cols": {
                "ra": [
                    0.0991523265838623,
                    -0.08193042138591411,
                    -0.1613445422612131,
                    -0.1104948495514691,
                    0.19883551495149732,
                ],
                "dec": [
                    0.011815003330658916,
                    -0.14095807328119264,
                    -0.1591939722196891,
                    0.08045418165202675,
                    0.07062092834004918,
                ],
                "epsilon_obs_1": [
                    0.2400706797017566,
                    0.12772301607339823,
                    -0.45334167657545993,
                    0.442176612285069,
                    0.35228716992034903,
                ],
                "epsilon_obs_2": [
                    0.012534815790335677,
                    0.6295368994121021,
                    0.2117515767653088,
                    0.09270305453109257,
                    0.2935061754831782,
                ],
            },
            "m2lnL": 13.696895880397472,
        },
        202: {
            "cols": {
                "ra": [
                    -0.09437827505171298,
                    -0.10091683054342865,
                    -0.035511021222919215,
                    0.19020439470186828,
                    -0.062485269457101825,
                ],
                "dec": [
                    -0.16729319407118493,
                    0.196417405822625,
                    -0.009351967416939246,
                    -0.12252664668627275,
                    -0.11067259364819969,
                ],
                "epsilon_obs_1": [
                    0.3895059484949051,
                    0.05908098971770451,
                    0.2164032326223062,
                    -0.5579878644228332,
                    -0.05055960539221188,
                ],
                "epsilon_obs_2": [
                    -0.029635225039541573,
                    -0.37187112768756386,
                    0.13540266942023446,
                    0.32227311248560575,
                    -0.1425934220136959,
                ],
            },
            "m2lnL": 9.917228051964102,
        },
    },
    ("local", "spline"): {
        200: {
            "cols": {
                "ra": [
                    -0.10938103049993514,
                    0.09081182871013882,
                    -0.15151121113449334,
                    0.18835634412243962,
                    -0.032229752186685806,
                ],
                "dec": [
                    -0.002342543371666539,
                    0.16387784329891236,
                    -0.15068294679505903,
                    0.13041990375585624,
                    -0.006571839549873736,
                ],
                "epsilon_obs_1": [
                    0.4291072559243664,
                    0.18082875225268236,
                    -0.09979544586912373,
                    -0.13968376790912976,
                    -0.15383344922364683,
                ],
                "epsilon_obs_2": [
                    -0.1054227767664152,
                    0.33697758528964006,
                    0.11883368833322652,
                    0.6653166774948261,
                    0.05104050617141078,
                ],
            },
            "m2lnL": 8.258653268444114,
        },
        201: {
            "cols": {
                "ra": [
                    0.0991523265838623,
                    -0.08193042138591411,
                    -0.1613445422612131,
                    -0.1104948495514691,
                    0.19883551495149732,
                ],
                "dec": [
                    0.011815003330658916,
                    -0.14095807328119264,
                    -0.1591939722196891,
                    0.08045418165202675,
                    0.07062092834004918,
                ],
                "epsilon_obs_1": [
                    0.2260505469622028,
                    0.14810011589625438,
                    -0.32042062024330925,
                    0.4468670679673635,
                    0.17970738387510146,
                ],
                "epsilon_obs_2": [
                    0.013170848912011113,
                    0.7653632893752472,
                    0.16541137227461508,
                    0.09356092298968915,
                    0.17903879818465998,
                ],
            },
            "m2lnL": 10.536290272434375,
        },
        202: {
            "cols": {
                "ra": [
                    -0.09437827505171298,
                    -0.10091683054342865,
                    -0.035511021222919215,
                    0.19020439470186828,
                    -0.062485269457101825,
                ],
                "dec": [
                    -0.16729319407118493,
                    0.196417405822625,
                    -0.009351967416939246,
                    -0.12252664668627275,
                    -0.11067259364819969,
                ],
                "epsilon_obs_1": [
                    0.36722867630732325,
                    0.06761230493541316,
                    0.164123600178707,
                    -0.5650374389486277,
                    -0.02011099707880146,
                ],
                "epsilon_obs_2": [
                    -0.02648397024810383,
                    -0.43717546582244327,
                    0.0815714056480282,
                    0.32591027176316734,
                    -0.06426499382088673,
                ],
            },
            "m2lnL": 7.030637658467775,
        },
    },
    ("local", "composed"): {
        200: {
            "cols": {
                "ra": [
                    -0.09695609295740726,
                    -0.0696467756293714,
                    -0.05022352645173668,
                    0.04346251515671612,
                    0.18913167547434567,
                ],
                "dec": [
                    0.03776798395647362,
                    0.19441098322911052,
                    0.17739633871756494,
                    0.12772756282210726,
                    -0.16053321201176812,
                ],
                "zp": [
                    0.29510266643050737,
                    0.42202300798875003,
                    0.43394969000964734,
                    0.7478088064032384,
                    0.1843946257314928,
                ],
                "epsilon_obs_1": [
                    -0.192399396880827,
                    -0.05276509717334017,
                    -0.07503601723955886,
                    0.06285365069629259,
                    -0.10534080367151903,
                ],
                "epsilon_obs_2": [
                    0.08845010170029108,
                    -0.18958692209328404,
                    0.48462383285136706,
                    0.46095263514892826,
                    -0.23307743902836947,
                ],
            },
            "m2lnL": -14.942575662331928,
        },
        201: {
            "cols": {
                "ra": [
                    0.12387672355398538,
                    0.08033701749518513,
                    0.1345349228009582,
                    -0.10026516951620579,
                    0.07227679193019866,
                ],
                "dec": [
                    0.08354788023234028,
                    0.059893120837202295,
                    0.11473177750808512,
                    0.17497459146381938,
                    -0.13259673179593884,
                ],
                "zp": [
                    0.6548159847504184,
                    0.18953662452973885,
                    1.2915012975325568,
                    1.1885385129355623,
                    0.3187867346606224,
                ],
                "epsilon_obs_1": [
                    -0.25133825435297685,
                    -0.5665287545950531,
                    -0.09043737246226338,
                    -0.14585441472782618,
                    -0.1582377022520448,
                ],
                "epsilon_obs_2": [
                    -0.27507764870163026,
                    0.2719542005829586,
                    -0.4970742154970388,
                    -0.24370796143310608,
                    0.006387202401500203,
                ],
            },
            "m2lnL": -9.844361094347786,
        },
        202: {
            "cols": {
                "ra": [
                    -0.05359764909371734,
                    -0.16520272409543396,
                    -0.12252680212259293,
                    0.16697374572977425,
                    -0.07057789023965598,
                ],
                "dec": [
                    0.19094442645985876,
                    -0.08152155538366941,
                    0.17366129982221604,
                    -0.02985580512227511,
                    0.10222493124029813,
                ],
                "zp": [
                    0.7914676254190821,
                    1.5022120317884804,
                    1.399243632416395,
                    0.22081792386533694,
                    0.37007993719054366,
                ],
                "epsilon_obs_1": [
                    -0.11917460520568143,
                    0.09121427264944222,
                    -0.019209314609747438,
                    -0.14231676277603644,
                    -0.007509292045900447,
                ],
                "epsilon_obs_2": [
                    -0.1848270641810301,
                    -0.6816503091552719,
                    -0.18408346009791204,
                    -0.3326162082648963,
                    0.130892529767581,
                ],
            },
            "m2lnL": -12.604102574374945,
        },
    },
}


@pytest.mark.parametrize("shape_kind,z_kind", _COMBOS, ids=_COMBO_IDS)
def test_resample_matches_legacy(shape_kind, z_kind):
    """``resample()`` -- redshift draw, then position draw+radius rejection,
    then shape draw, all sharing one RNG stream -- must reproduce legacy
    seed-for-seed under the same seed for every combination, not just the
    already-covered GaussGlobal+Composed one. Now checks the raw regenerated
    columns (ra/dec/epsilon_obs_1/2, identically named on both sides for
    every combo; zp too when z=composed) and the resulting -2lnL against
    frozen legacy values (see module docstring).

    shape=local was not bit-identical (rtol=0): both sides draw the same
    underlying uniforms, but GaussLocal resolves its per-galaxy sigma from
    e_rms via an independent bisection (self-consistent with its own
    forward map), while legacy's HSMGauss resolves sigma from std_shape via
    a closed-form inversion of an algebraically equivalent but differently
    arranged formula -- the same rtol=1e-8 caveat already established at
    the single-galaxy level (test_galaxy_shape_pop_gauss_local.py). ra/dec
    (position, shape-independent) and shape=global were bit identical when
    this was captured (rtol=1e-12 kept here only to absorb cross-platform
    libm/ULP noise in the RNG's transcendental calls, not a real tolerance).
    """
    mset, hms = _build_mset(shape_kind, z_kind)
    hms.param_set_by_name("log10MDelta", 14.0)

    position_factor, redshift_factor, shape_factor, new_obs = _build_new_obs(
        mset, shape_kind, z_kind
    )

    dcwlf = Nc.DataClusterWLFactor.new(position_factor, redshift_factor, shape_factor)
    dcwlf.set_obs(new_obs)
    dcwlf.set_prec(1.0e-8)
    dcwlf.set_integ_method(Nc.DataClusterWLIntegMethod.LNINT)

    n = len(_GALAXIES)
    bit_exact_cols = ["ra", "dec"]
    shape_cols = ["epsilon_obs_1", "epsilon_obs_2"]
    if z_kind == "composed":
        bit_exact_cols.append("zp")

    for seed in (200, 201, 202):
        rng_new = Ncm.RNG.seeded_new(None, seed)

        dcwlf.resample(mset, rng_new)

        frozen = _RESAMPLE_FROZEN[(shape_kind, z_kind)][seed]

        for col in bit_exact_cols:
            new_vals = [new_obs.get(col, i) for i in range(n)]
            # rtol=1e-12 (not bit-exact): ra/dec/zp are RNG draws routed
            # through platform/libm-sensitive transcendentals (asin() in the
            # sky-footprint sampler, an ODE-spline inverse-CDF for composed
            # redshift) -- see the same reasoning in
            # test_galaxy_redshift_factor_spline_legacy_parity.py.
            assert_allclose(new_vals, frozen["cols"][col], rtol=1.0e-12, atol=0.0)

        for col in shape_cols:
            new_vals = [new_obs.get(col, i) for i in range(n)]
            assert_allclose(
                new_vals, frozen["cols"][col], rtol=1.0e-8, atol=0.0
            )

        new_m2lnL = dcwlf.m2lnL_val(mset)
        assert_allclose(new_m2lnL, frozen["m2lnL"], rtol=1.0e-5)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

#!/usr/bin/env python
#
# test_galaxy_redshift_factor_spline_legacy_parity.py
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

"""Golden-parity test: ``NcGalaxyRedshiftFactorSpline`` vs legacy
``NcGalaxySDObsRedshiftPz``.

``test_galaxy_redshift_factor_spline.py`` states there is "no legacy class
sharing identical math" to use as an oracle -- that is not quite right:
legacy's ``NcGalaxySDObsRedshiftPz`` and the new ``NcGalaxyRedshiftFactorSpline``
both do exactly the same thing (hand out a per-galaxy pre-tabulated p(z) as-is,
built via the identical ``-2*log(y+1e-5)`` inverse-CDF construction for
``gen``), reading/writing through the very same ``NcGalaxyWLObs`` dedicated
``pz`` slot (``nc_galaxy_wl_obs_{peek,set}_pz``). This test originally fed the
exact same ``NcmSpline`` into both and checked bit-for-bit (or seed-for-seed)
agreement, closing the one redshift scheme that had never actually been
cross-checked against its legacy counterpart.

FROZEN REFERENCE VALUES: the parity documented above was proven by running
both engines live and is captured, not re-derived, in the functions below.
Every comparison here was bit-for-bit exact (``rtol=0, atol=0``) between the
new and legacy engines. Values were captured from an actual passing run of
this file's original legacy-comparison code, at git rev ``77313f22``
(2026-07-16), then legacy (``NcGalaxySDObsRedshiftPz``) construction was
removed so these tests no longer depend on legacy at runtime -- legacy is
slated for deletion in a follow-up PR.
"""

import pytest
import numpy as np
from numpy.testing import assert_allclose

from numcosmo_py import Ncm, Nc

Ncm.cfg_init()

# (zp, sigma0, n); the spline domain is derived from (zp, sigma0) below --
# NOT set independently -- since a domain far wider than the peak turns the
# tabulated density into a near-delta-function integrand that breaks
# NcmStatsDist1dSpline's ODE-based inverse-CDF construction (both engines'
# `gen()` build one internally); see adaptive-spline-degenerate-integrand-oom.
_CASES = [
    (0.6, 0.05, 200),
    (0.9, 0.10, 150),
    (0.3, 0.02, 300),
]

_HALF_WIDTH_SIGMAS = 6.0


def _pz_bounds(zp, sigma0):
    z_min = max(1.0e-3, zp - _HALF_WIDTH_SIGMAS * sigma0)
    z_max = zp + _HALF_WIDTH_SIGMAS * sigma0
    return z_min, z_max


def _make_pz_spline(zp, sigma0, n):
    """A simple Gaussian-shaped p(z) spline on a domain proportional to
    sigma0 -- the exact physical shape does not matter here, only that both
    engines are fed the identical spline object."""
    z_min, z_max = _pz_bounds(zp, sigma0)
    zs = np.linspace(z_min, z_max, n)
    pz_vals = np.exp(-0.5 * ((zs - zp) / sigma0) ** 2) + 1.0e-6

    xv = Ncm.Vector.new_array(zs.tolist())
    yv = Ncm.Vector.new_array(pz_vals.tolist())
    spline = Ncm.SplineCubicNotaknot.new()
    spline.set(xv, yv, True)
    return spline


def _build_new(spline):
    gsdrs = Nc.GalaxyRedshiftFactorSpline.new()
    mset = Ncm.MSet.empty_new()
    data = Nc.GalaxyRedshiftFactorData.new(gsdrs, mset)
    gsdrs.data_set(data, spline)
    gsdrs.update_data(data)
    return gsdrs, mset, data


# Frozen legacy `integ()` output, keyed by (zp, sigma0, n, use_lnp), sampled
# on a 37-point z-grid spanning [z_min, z_max] (see module docstring).
_INTEG_FROZEN = {
    (0.6, 0.05, 200, False): [
        1.0152299797447125e-06,
        1.1064507465442652e-06,
        1.6658350524623818e-06,
        4.726598237171647e-06,
        1.9664424263284366e-05,
        8.464778722655721e-05,
        0.0003364619629935719,
        0.001204857086311452,
        0.0038669161072865644,
        0.01110999019706973,
        0.028566494094902153,
        0.06572952766350386,
        0.13533630294626703,
        0.24935321939544922,
        0.4111133483361653,
        0.6065316648974874,
        0.8007383874732328,
        0.9459604685496172,
        1.0000008962289033,
        0.9459604685496172,
        0.8007383874732328,
        0.606531664897486,
        0.4111133483361641,
        0.24935321939544922,
        0.13533630294626703,
        0.06572952766350382,
        0.02856649409490207,
        0.011109990197069737,
        0.003866916107286546,
        0.0012048570863114423,
        0.00033646196299357017,
        8.464778722655726e-05,
        1.9664424263284403e-05,
        4.726598237171627e-06,
        1.6658350524623818e-06,
        1.1064507465442644e-06,
        1.0152299797447125e-06,
    ],
    (0.6, 0.05, 200, True): [
        -13.800395390106855,
        -13.714353191284474,
        -13.305184027247073,
        -12.262304802959367,
        -10.836699429738461,
        -9.377011589991916,
        -7.997025452434879,
        -6.721394319644531,
        -5.555297961113495,
        -4.499910557682947,
        -3.555520783083354,
        -2.722207022743367,
        -1.999992465336545,
        -1.3888848359225647,
        -0.8888863158020672,
        -0.49999834273172344,
        -0.22222099266080061,
        -0.05555449880592045,
        8.962285017003174e-07,
        -0.05555449880592045,
        -0.22222099266080061,
        -0.49999834273172583,
        -0.8888863158020702,
        -1.3888848359225647,
        -1.999992465336545,
        -2.7222070227433677,
        -3.555520783083357,
        -4.499910557682947,
        -5.555297961113499,
        -6.721394319644539,
        -7.997025452434884,
        -9.377011589991916,
        -10.83669942973846,
        -12.262304802959372,
        -13.305184027247073,
        -13.714353191284474,
        -13.800395390106855,
    ],
    (0.9, 0.1, 150, False): [
        1.0152299797447125e-06,
        1.1064506900509827e-06,
        1.6658017926713652e-06,
        4.726463763195915e-06,
        1.9663812250980317e-05,
        8.464693014683659e-05,
        0.0003364611609602864,
        0.0012048599772707834,
        0.003866915524013956,
        0.011109975341010792,
        0.02856647024240528,
        0.06572952748832972,
        0.1353363359257702,
        0.24935326721891435,
        0.41111329970654803,
        0.606531666871117,
        0.8007383649574067,
        0.9459602514493405,
        1.0000006686476341,
        0.9459602514493405,
        0.8007383649574067,
        0.6065316668711169,
        0.4111132997065474,
        0.24935326721891343,
        0.1353363359257702,
        0.0657295274883297,
        0.02856647024240529,
        0.011109975341010759,
        0.0038669155240139545,
        0.001204859977270779,
        0.0003364611609602834,
        8.464693014683641e-05,
        1.966381225098038e-05,
        4.726463763195866e-06,
        1.665801792671364e-06,
        1.1064506900509827e-06,
        1.0152299797447127e-06,
    ],
    (0.9, 0.1, 150, True): [
        -13.800395390106855,
        -13.714353242342584,
        -13.305203993283326,
        -12.262333253841364,
        -10.836730553041123,
        -9.377021715289787,
        -7.997027836164184,
        -6.721391920226463,
        -5.55529811195014,
        -4.499911894864043,
        -3.5555216180653,
        -2.722207025408442,
        -1.9999922216510118,
        -1.3888846441325382,
        -0.8888864340896887,
        -0.4999983394777638,
        -0.2222210207796304,
        -0.05555472830843022,
        6.686474106052473e-07,
        -0.05555472830843022,
        -0.2222210207796304,
        -0.49999833947776395,
        -0.8888864340896903,
        -1.388884644132542,
        -1.9999922216510118,
        -2.7222070254084425,
        -3.5555216180653,
        -4.499911894864046,
        -5.55529811195014,
        -6.721391920226467,
        -7.997027836164193,
        -9.37702171528979,
        -10.83673055304112,
        -12.262333253841375,
        -13.305203993283328,
        -13.714353242342584,
        -13.800395390106855,
    ],
    (0.3, 0.02, 300, False): [
        1.0152299797447127e-06,
        1.1064532479092732e-06,
        1.6658336589007822e-06,
        4.726652381860588e-06,
        1.966444632238221e-05,
        8.46482121781428e-05,
        0.0003364625254051364,
        0.0012048597847490944,
        0.0038669186064227445,
        0.011109995287202249,
        0.0285665007047975,
        0.06572952871227786,
        0.13533628667168773,
        0.2493532088607586,
        0.4111132985656591,
        0.6065316676317574,
        0.8007384028001746,
        0.9459604626554236,
        1.0000009796902487,
        0.9459604626554236,
        0.8007384028001746,
        0.6065316676317591,
        0.41111329856565915,
        0.2493532088607586,
        0.13533628667168773,
        0.06572952871227784,
        0.028566500704797704,
        0.011109995287202253,
        0.0038669186064227436,
        0.0012048597847491003,
        0.00033646252540513633,
        8.464821217814287e-05,
        1.9664446322382207e-05,
        4.7266523818606074e-06,
        1.6658336589007883e-06,
        1.1064532479092732e-06,
        1.0152299797447127e-06,
    ],
    (0.3, 0.02, 300, True): [
        -13.800395390106855,
        -13.714350930576359,
        -13.305184863801797,
        -12.26229334770633,
        -10.836698307962141,
        -9.377006569771863,
        -7.997023780890525,
        -6.72139208001408,
        -5.555297314827073,
        -4.4999100995249055,
        -3.555520551697082,
        -2.7222070067874635,
        -1.9999925855894256,
        -1.3888848781706291,
        -0.8888864368648093,
        -0.49999833822368195,
        -0.22222097351979037,
        -0.055554505036829564,
        9.796897688074124e-07,
        -0.055554505036829564,
        -0.22222097351979037,
        -0.49999833822367906,
        -0.8888864368648091,
        -1.3888848781706291,
        -1.9999925855894256,
        -2.722207006787464,
        -3.5555205516970747,
        -4.4999100995249055,
        -5.555297314827074,
        -6.721392080014075,
        -7.997023780890525,
        -9.377006569771863,
        -10.836698307962143,
        -12.262293347706324,
        -13.305184863801793,
        -13.714350930576359,
        -13.800395390106855,
    ],
}


@pytest.mark.parametrize("zp,sigma0,n", _CASES)
@pytest.mark.parametrize("use_lnp", [False, True])
def test_integ_bit_parity(zp, sigma0, n, use_lnp):
    """The new engine's direct spline evaluation is checked against frozen
    legacy values, bit-for-bit (rtol=0, atol=0; see module docstring)."""
    z_min, z_max = _pz_bounds(zp, sigma0)
    spline = _make_pz_spline(zp, sigma0, n)
    gsdrs, mset, new_data = _build_new(spline)

    new_integ = gsdrs.integ(mset, use_lnp)

    zs = np.linspace(z_min, z_max, 37)
    got = np.array([new_integ.eval(z, new_data) for z in zs])

    assert_allclose(got, _INTEG_FROZEN[(zp, sigma0, n, use_lnp)], rtol=0.0, atol=0.0)


# Frozen legacy `get_integ_lim()` output, keyed by (zp, sigma0, n).
_INTEG_LIM_FROZEN = {
    (0.6, 0.05, 200): (0.29999999999999993, 0.9),
    (0.9, 0.1, 150): (0.29999999999999993, 1.5),
    (0.3, 0.02, 300): (0.18, 0.42),
}


@pytest.mark.parametrize("zp,sigma0,n", _CASES)
def test_get_integ_lim_bit_parity(zp, sigma0, n):
    """Now checked against frozen legacy values (see module docstring)."""
    spline = _make_pz_spline(zp, sigma0, n)
    gsdrs, mset, new_data = _build_new(spline)

    new_lim = gsdrs.get_integ_lim(mset, new_data)

    assert_allclose(new_lim, _INTEG_LIM_FROZEN[(zp, sigma0, n)], rtol=0.0, atol=0.0)


# Frozen legacy `norm()` output, keyed by (zp, sigma0, n).
_NORM_FROZEN = {
    (0.6, 0.05, 200): 0.12533201348428158,
    (0.9, 0.1, 150): 0.250664026968786,
    (0.3, 0.02, 300): 0.050132805393701345,
}


@pytest.mark.parametrize("zp,sigma0,n", _CASES)
def test_norm_bit_parity(zp, sigma0, n):
    """Now checked against frozen legacy values (see module docstring)."""
    spline = _make_pz_spline(zp, sigma0, n)
    gsdrs, mset, new_data = _build_new(spline)

    new_norm = gsdrs.norm(mset, new_data)

    assert_allclose(new_norm, _NORM_FROZEN[(zp, sigma0, n)], rtol=0.0, atol=0.0)


# Frozen legacy seed=7531 draw sequence (50 draws), keyed by (zp, sigma0, n).
_GEN_FROZEN = {
    (0.6, 0.05, 200): [
        0.594581595584194,
        0.555697229785205,
        0.6080573137431238,
        0.5420421054820126,
        0.5442755165336454,
        0.6342414966592744,
        0.6526172528995584,
        0.6729256798521954,
        0.5873666283318951,
        0.6050349987807394,
        0.6911438080791591,
        0.61492319779363,
        0.6048759913483023,
        0.5366277323783374,
        0.5951398256595205,
        0.571641726733605,
        0.5996558966460619,
        0.6288000912444185,
        0.5914830915519061,
        0.6128525340859615,
        0.5724799378132301,
        0.5932646998764642,
        0.6436724424540712,
        0.5480059816276426,
        0.5601399259288156,
        0.6286020741794702,
        0.5424803680622045,
        0.5402702715419719,
        0.5725392590573334,
        0.6000156349452272,
        0.6204410587111854,
        0.5758991396749499,
        0.5189874898570841,
        0.6151371551764132,
        0.6204972426831086,
        0.6268187972317772,
        0.6267119486599798,
        0.527478051928146,
        0.615728727061574,
        0.630043739523837,
        0.5905333860661707,
        0.6141786635836178,
        0.6947010364568685,
        0.5811397678253121,
        0.6958973070661899,
        0.6531835534204381,
        0.5187551551712194,
        0.6731026678219494,
        0.5086140840006861,
        0.5769115984670666,
    ],
    (0.9, 0.1, 150): [
        0.8891655219628054,
        0.811397226063801,
        0.9161172315392787,
        0.7840869031889874,
        0.7885532006777156,
        0.9684855519456753,
        1.0052373716873606,
        1.0458542002988687,
        0.8747360105104054,
        0.910072599547303,
        1.0822911270694402,
        0.9298489558785393,
        0.9097545825541659,
        0.7732580707105441,
        0.8902820150672441,
        0.8432862399849353,
        0.8993143089672494,
        0.9576027433825999,
        0.8829688105289208,
        0.925707660161412,
        0.844962484266227,
        0.8865317597120912,
        0.9873473082844904,
        0.7960142850088232,
        0.8202819591284111,
        0.957206705614617,
        0.7849633353592235,
        0.7805435771986611,
        0.845081113950735,
        0.9000338012212196,
        0.9408847125543883,
        0.8518004732159558,
        0.7379774620347033,
        0.9302768607979172,
        0.9409970805412717,
        0.9536401775812327,
        0.9534264838774646,
        0.7549591362274207,
        0.9314599758408789,
        0.9600900680524879,
        0.8810694825396358,
        0.9283599147457802,
        1.0894067182801999,
        0.8622820831115049,
        1.0917996203009512,
        1.006369986198676,
        0.7375128300192322,
        1.0462081702579005,
        0.7172306508756606,
        0.8538254508053371,
    ],
    (0.3, 0.02, 300): [
        0.2978319434738058,
        0.2822784401866003,
        0.30322225458461494,
        0.2768164127175716,
        0.2777096435773334,
        0.3136956911964138,
        0.32104569336220506,
        0.32916829376482853,
        0.29494609335557415,
        0.30201333080383913,
        0.3364539191787232,
        0.3059685055277113,
        0.3019497276913273,
        0.2746505543528892,
        0.2980552390359851,
        0.28865618854464525,
        0.29986168100624655,
        0.3115191892884721,
        0.296592621683406,
        0.3051403051578483,
        0.2889914363842377,
        0.29730519475747774,
        0.317467862575722,
        0.27920180164697206,
        0.2840553335628005,
        0.31143998322387595,
        0.27699169025908243,
        0.2761076924443385,
        0.2890151619158442,
        0.3000055790126009,
        0.30817562994429865,
        0.29035900653483193,
        0.2675945651877044,
        0.30605408106085147,
        0.3081981029536844,
        0.31072668534047915,
        0.31068394738184557,
        0.27099080591886615,
        0.306290692335714,
        0.3120166410274321,
        0.296212772390249,
        0.3056707190223793,
        0.33787633153381297,
        0.29245529620713767,
        0.3383546803327908,
        0.3212722059560783,
        0.26750164703068025,
        0.3292390783008585,
        0.2634452148624445,
        0.2907639913453698,
    ],
}


@pytest.mark.parametrize("zp,sigma0,n", _CASES)
def test_gen_matches_seed_for_seed(zp, sigma0, n):
    """Same seed -> same inverse-CDF construction -> identical draws.

    The new Spline scheme builds its lazy `dist` inside `gen()`, using the
    exact same -2*log(y+1e-5) transform, NcmStatsDist1dSpline with
    reltol=1e-5, and a do-while rejection loop against [z_min, z_max] that
    legacy's `NcGalaxySDObsRedshiftPz` used (which built its own lazily
    inside `prepare()`) -- so the RNG call sequence was identical, checked
    here against a frozen legacy draw sequence (see module docstring).
    """
    spline = _make_pz_spline(zp, sigma0, n)
    gsdrs, mset, new_data = _build_new(spline)

    n_draws = 50
    new_zs = np.empty(n_draws)

    rng_new = Ncm.RNG.seeded_new(None, 7531)

    for i in range(n_draws):
        gsdrs.gen(mset, new_data, rng_new)
        new_zs[i] = new_data.z

    # rtol=1e-12 (not bit-exact): gen() routes through NcmStatsDist1d's
    # inverse-CDF spline, built by a GSL adaptive-ODE solve
    # (ncm_ode_spline_prepare) and evaluated via atanh() -- both genuinely
    # sensitive to the platform's libm/compiler at the ULP level, unlike the
    # other frozen comparisons in this file which only evaluate @pz's
    # cubic spline directly (see module docstring). Observed on CI: 1/50
    # elements off by exactly 1 ULP on a runner different from the one that
    # captured these constants.
    assert_allclose(new_zs, _GEN_FROZEN[(zp, sigma0, n)], rtol=1.0e-12, atol=1.0e-12)


@pytest.mark.parametrize("zp,sigma0,n", _CASES)
def test_required_columns_bit_parity(zp, sigma0, n):
    """Neither scheme adds a scalar column: the spline lives on the
    NcGalaxyWLObs row's own dedicated pz slot for both -- legacy
    (``NcGalaxySDObsRedshiftPz``) was verified to return the identical
    ``["z"]`` list (see module docstring)."""
    spline = _make_pz_spline(zp, sigma0, n)
    _gsdrs, _mset, new_data = _build_new(spline)

    assert list(Nc.GalaxyRedshiftFactorData.required_columns(new_data)) == ["z"]


# Frozen legacy read-row/eval output, keyed by (zp, sigma0, n), sampled on an
# 11-point z-grid spanning [z_min, z_max].
_READ_ROW_FROZEN = {
    (0.6, 0.05, 200): [
        1.0152299797447125e-06,
        1.0929493004249398e-05,
        0.0015348088526178417,
        0.056135759784925236,
        0.48675331192013227,
        1.0000008962289033,
        0.48675331192013094,
        0.05613575978492524,
        0.0015348088526178397,
        1.0929493004249329e-05,
        1.0152299797447125e-06,
    ],
    (0.9, 0.1, 150): [
        1.0152299797447125e-06,
        1.0929475137298975e-05,
        0.0015348050756630628,
        0.05613575212777439,
        0.4867534353721662,
        1.0000006686476341,
        0.4867534353721662,
        0.056135752127774226,
        0.0015348050756630507,
        1.0929475137298907e-05,
        1.0152299797447127e-06,
    ],
    (0.3, 0.02, 300): [
        1.0152299797447127e-06,
        1.0929501674445409e-05,
        0.0015348103082824475,
        0.05613576229435426,
        0.4867532668697522,
        1.0000009796902487,
        0.48675326686975057,
        0.05613576229435427,
        0.0015348103082824475,
        1.0929501674445345e-05,
        1.0152299797447127e-06,
    ],
}


@pytest.mark.parametrize("zp,sigma0,n", _CASES)
def test_read_row_write_row_shared_obs(zp, sigma0, n):
    """The new scheme reads/writes through the NcGalaxyWLObs pz slot, so a
    catalog row round-trips through it identically to how it round-tripped
    through legacy's identical slot -- checked here against frozen legacy
    values (see module docstring)."""
    z_min, z_max = _pz_bounds(zp, sigma0)
    spline = _make_pz_spline(zp, sigma0, n)
    gsdrs = Nc.GalaxyRedshiftFactorSpline.new()
    mset = Ncm.MSet.empty_new()

    new_data = Nc.GalaxyRedshiftFactorData.new(gsdrs, mset)

    wlobs = Nc.GalaxyWLObs.new(
        Nc.GalaxyWLObsEllipConv.TRACE_DET, Nc.WLEllipticityFrame.CELESTIAL, 1, ["z"]
    )
    wlobs.set("z", 0, 0.0)
    wlobs.set_pz(0, spline)

    new_data.read_row(wlobs, 0)

    new_integ = gsdrs.integ(mset, False)

    zs = np.linspace(z_min, z_max, 11)
    got = np.array([new_integ.eval(z, new_data) for z in zs])
    assert_allclose(got, _READ_ROW_FROZEN[(zp, sigma0, n)], rtol=0.0, atol=0.0)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

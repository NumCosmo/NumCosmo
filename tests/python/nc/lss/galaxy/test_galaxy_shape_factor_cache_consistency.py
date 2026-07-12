#!/usr/bin/env python
#
# test_galaxy_shape_factor_cache_consistency.py
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

"""Direct regression test for NcGalaxyShapeFactor's hash/cache-refresh scheme.

``nc_galaxy_shape_factor_prepare()`` caches ``pr_prefactor`` (guarded by
``radius_hash``, driven by ``halo_position``+``cosmo``), ``lens_ctx``/``r_s``/
``rho_s`` (guarded by ``optzs_hash``, driven by ``cosmo``+``density_profile``+
its ``NcHaloMassSummary`` submodel+``surface_mass_density``+``halo_position``),
and refreshes ``pop_data`` unconditionally (guarded by ``pop_hash``, driven by
the ``NcGalaxyShapePop`` model) -- see the class documentation. Two real bugs
were found in this scheme by contrasting its output against a "no caching"
reference built the slow way (a brand-new factor object every call, so there
is nothing to go stale):

- ``optzs_hash`` didn't include the ``NcHaloMassSummary`` submodel's own
  pkey (a submodel's pkey never propagates into its parent's own top-level
  pkey outside ``NcmModelCtrl``, which this hash scheme deliberately doesn't
  use) -- revisiting a previously-seen mass value left ``r_s``/``rho_s``
  stale at whatever mass had last triggered a refresh.
- ``optzs_hash`` read ``z_cl`` from ``halo_position`` inside its refresh
  block but never hashed ``halo_position``'s own pkey -- changing only the
  lens redshift (holding cosmo/profile/mass fixed) left ``lens_ctx``/``r_s``/
  ``rho_s`` stale at the old ``z_cl``.

Both are now fixed (nc_galaxy_shape_factor.c's ``optzs_hash``). This test
exists so that any *future* omission in this hash scheme is caught directly
and cheaply here, rather than showing up only as an unexplained discrepancy
several layers up (that's how both bugs above were actually found: an
``NcDataClusterWLFactor`` Monte-Carlo mass-recovery run silently stopped
responding to mass changes after a resample() cycle).

Methodology: build one "cached" ``NcGalaxyShapeFactorVarAdd`` reused across a
sequence of mset configurations (deliberately *revisiting* earlier values,
the case that caught both bugs above -- a monotonic sweep would not), and one
"fresh" instance rebuilt from scratch at every step (so it can never hold
stale cached state by construction). If the integrand agrees at every step,
the caching is transparent; any disagreement is a cache-invalidation bug.
"""

import pytest
import numpy as np
from numpy.testing import assert_allclose

from numcosmo_py import Ncm, Nc

Ncm.cfg_init()

ELLIP_CONV = Nc.GalaxyWLObsEllipConv.TRACE_DET
GALAXY = dict(ra=0.03, dec=-0.02, z=0.6, e1=0.05, e2=-0.02, std_noise=0.03, c1=0.0, c2=0.0, m=0.0)


def _build_persistent_mset():
    """Build the mset's constituent NcmModels *once*. Real usage (and both
    bugs in the module docstring) only ever mutates a model's own parameter
    *values* over its lifetime -- ncm_model_state_get_pkey() is a per-object
    counter, not a value hash, so comparing pkeys across brand-new model
    instances built fresh each step (as opposed to mutating persistent ones)
    is meaningless and would produce spurious coincidental hash collisions.
    """
    cosmo = Nc.HICosmoDEXcdm.new()
    dist = Nc.Distance.new(100.0)
    hms = Nc.HaloCMParam.new(Nc.HaloMassSummaryMassDef.MEAN, 200.0)
    dp = Nc.HaloDensityProfileNFW.new(hms)
    hp = Nc.HaloPosition.new(dist)
    smd = Nc.WLSurfaceMassDensity.new(dist)
    pop = Nc.GalaxyShapePopGauss.new()
    pop_z = Nc.GalaxyRedshiftPopLSSTSRD.new_y1_source()
    obs_z = Nc.GalaxyRedshiftObsGauss.new()

    mset = Ncm.MSet.empty_new()
    for model in (cosmo, dp, hp, smd, pop, pop_z, obs_z):
        mset.set(model)

    return mset, cosmo, hms, hp, pop


def _set_config(mset, cosmo, hms, hp, pop, config, prev_config):
    """Apply `config`, touching only the parameters that actually differ from
    `prev_config` (None on the first call, meaning "set everything").

    This is not a cosmetic optimization: ncm_model_param_set bumps a model's
    pkey unconditionally on every call, even when re-assigning the value it
    already has. Unconditionally re-setting all four knobs every step (as an
    earlier version of this test did) makes every hash "change" every step
    regardless of which knob actually moved, silently defeating this test's
    entire purpose -- it would have reported both real bugs in the module
    docstring as passing. Real callers (e.g. a Fit engine varying one
    parameter across iterations) only ever touch the parameter(s) actually
    being varied, so this mirrors real usage, not just this test's needs.
    """
    log10_mdelta, c_delta, hp_z, pop_sigma = config
    plog10_mdelta, pc_delta, php_z, ppop_sigma = prev_config if prev_config is not None else (None,) * 4

    if log10_mdelta != plog10_mdelta:
        hms.param_set_by_name("log10MDelta", log10_mdelta)
    if c_delta != pc_delta:
        hms.param_set_by_name("cDelta", c_delta)
    if hp_z != php_z:
        hp.param_set_by_name("z", hp_z)
        hp.prepare(cosmo)
    if pop_sigma != ppop_sigma:
        pop.param_set_by_name("sigma", pop_sigma)
    return mset


def _build_mset(config):
    """A brand-new mset (new model instances) for the given config -- used
    as the "no caching" reference side, where fresh-instance pkeys are never
    compared against anything, only the resulting integrand values are."""
    mset, cosmo, hms, hp, pop = _build_persistent_mset()
    _set_config(mset, cosmo, hms, hp, pop, config, None)
    return mset


def _make_data(gsf, mset):
    posf = Nc.GalaxyPositionFactorFlat.new(-0.2, 0.2, -0.2, 0.2)
    pos_data = Nc.GalaxyPositionFactorData.new(posf, mset)
    zf = Nc.GalaxyRedshiftFactorComposed.new(0.0, 20.0)
    z_data = Nc.GalaxyRedshiftFactorData.new(zf, mset)
    data = Nc.GalaxyShapeFactorData.new(gsf, mset, pos_data, z_data)

    pos_data.ra = GALAXY["ra"]
    pos_data.dec = GALAXY["dec"]
    z_data.z = GALAXY["z"]
    gsf.data_set(
        data, GALAXY["e1"], GALAXY["e2"], GALAXY["std_noise"],
        GALAXY["c1"], GALAXY["c2"], GALAXY["m"], Nc.WLEllipticityFrame.CELESTIAL,
    )
    return data


# (log10_mdelta, c_delta, hp_z, pop_sigma) -- deliberately revisits earlier
# values (mass, concentration, lens redshift, pop sigma each get revisited at
# least once) since that is the case that caught both bugs in the module
# docstring; a monotonic sweep would not.
_SEQUENCE = [
    (14.0, 4.0, 0.2, 0.25),
    (14.5, 4.0, 0.2, 0.25),  # mass changes
    (14.0, 4.0, 0.2, 0.25),  # mass reverts to a previously-seen value
    (14.0, 5.0, 0.2, 0.25),  # concentration changes
    (14.0, 4.0, 0.2, 0.25),  # concentration reverts
    (14.0, 4.0, 0.35, 0.25),  # lens redshift changes (mass/profile/smd fixed)
    (14.0, 4.0, 0.2, 0.25),  # lens redshift reverts
    (14.0, 4.0, 0.2, 0.35),  # pop sigma changes
    (14.0, 4.0, 0.2, 0.25),  # pop sigma reverts
]

_Z_EVAL = np.linspace(0.05, 1.5, 15)


@pytest.mark.parametrize("use_lnp", [False, True])
def test_cached_matches_fresh_across_revisits(use_lnp):
    """A reused (cached) factor, its models mutated in place across the whole
    sequence (touching only the knob(s) that actually change each step, see
    _set_config), must match a freshly-built factor+models at every step."""
    gsf_cached = Nc.GalaxyShapeFactorVarAdd.new(ELLIP_CONV)
    mset_cached, cosmo, hms, hp, pop = _build_persistent_mset()

    prev_config = None
    for config in _SEQUENCE:
        _set_config(mset_cached, cosmo, hms, hp, pop, config, prev_config)
        prev_config = config

        data_cached = _make_data(gsf_cached, mset_cached)
        gsf_cached.prepare_data_array(mset_cached, [data_cached], True, True)
        integ_cached = gsf_cached.integ(mset_cached, use_lnp)

        mset_fresh = _build_mset(config)
        gsf_fresh = Nc.GalaxyShapeFactorVarAdd.new(ELLIP_CONV)
        data_fresh = _make_data(gsf_fresh, mset_fresh)
        gsf_fresh.prepare_data_array(mset_fresh, [data_fresh], True, True)
        integ_fresh = gsf_fresh.integ(mset_fresh, use_lnp)

        cached_vals = [integ_cached.eval(z, data_cached) for z in _Z_EVAL]
        fresh_vals = [integ_fresh.eval(z, data_fresh) for z in _Z_EVAL]

        assert_allclose(
            cached_vals, fresh_vals, rtol=1.0e-12, atol=0.0,
            err_msg=f"cached vs fresh mismatch at config={config}",
        )


def test_hash_changes_when_driving_model_changes():
    """Each hash must change whenever its own driving model(s) differ from
    the previous prepare() cycle -- including on a revisit (a hash should go
    back to reflecting the revisited state, not stay frozen). This only
    checks the necessary direction (a real change must never be missed); the
    converse doesn't hold in general, since a driving model's own prepare()
    can legitimately bump its pkey from lazy internal caching unrelated to
    any value actually visible to this hash -- that's spurious-but-harmless
    over-invalidation, not the under-invalidation bug this guards against
    (see test_cached_matches_fresh_across_revisits for the real end-to-end
    correctness proof, which does catch both directions)."""
    gsf = Nc.GalaxyShapeFactorVarAdd.new(ELLIP_CONV)
    mset, cosmo, hms, hp, pop = _build_persistent_mset()

    prev_radius_hash = None
    prev_optzs_hash = None
    prev_pop_hash = None
    prev_config = None

    for config in _SEQUENCE:
        log10_mdelta, c_delta, hp_z, pop_sigma = config
        _set_config(mset, cosmo, hms, hp, pop, config, prev_config)
        gsf.prepare(mset)

        radius_hash = gsf.get_radius_hash()
        optzs_hash = gsf.get_optzs_hash()
        pop_hash = gsf.get_pop_hash()

        if prev_config is not None:
            plog10_mdelta, pc_delta, php_z, ppop_sigma = prev_config

            radius_should_change = hp_z != php_z
            optzs_should_change = (log10_mdelta, c_delta, hp_z) != (plog10_mdelta, pc_delta, php_z)
            pop_should_change = pop_sigma != ppop_sigma

            if radius_should_change:
                assert radius_hash != prev_radius_hash, config
            if optzs_should_change:
                assert optzs_hash != prev_optzs_hash, config
            if pop_should_change:
                assert pop_hash != prev_pop_hash, config

        prev_radius_hash, prev_optzs_hash, prev_pop_hash = radius_hash, optzs_hash, pop_hash
        prev_config = config


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

#!/usr/bin/env python
#
# test_data_cluster_wl_factor_mask.py
#
# Mon Nov 16 2026
# Copyright  2026  Caio Lima de Oliveira
# <caiolimadeoliveira@pm.me>
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

"""Tests of the mask of NcDataClusterWLFactor: the rows of its catalog it uses.

The likelihood of a mask must equal, bit for bit, that of the catalog made of its rows
alone, whatever the redshift integration; the other rows are never read nor written.
"""

import subprocess
import sys

import numpy as np
import pytest

from numcosmo_py import Nc, Ncm
from numcosmo_py.helper import duplicate_via_serialization

Ncm.cfg_init()

ELLIP_CONV = Nc.GalaxyWLObsEllipConv.TRACE_DET
FRAME = Nc.WLEllipticityFrame.CELESTIAL

# ra, dec, zp, sigma0, e1, e2, std_noise
GALAXIES = [
    (0.03, 0.02, 0.60, 0.030, 0.05, -0.02, 0.03),
    (-0.10, 0.15, 0.90, 0.040, -0.04, 0.01, 0.05),
    (0.05, -0.08, 0.15, 0.020, 0.02, 0.03, 0.04),
    (0.08, 0.05, 1.50, 0.050, 0.03, -0.01, 0.04),
    (-0.05, -0.03, 0.35, 0.025, -0.01, 0.02, 0.035),
    (0.02, 0.09, 0.75, 0.035, 0.01, -0.03, 0.045),
    (-0.07, 0.04, 0.55, 0.030, -0.02, -0.02, 0.04),
    (0.06, -0.02, 1.10, 0.045, 0.04, 0.02, 0.05),
]
MASK = [0, 2, 3, 5, 7]

METHODS = [
    (Nc.DataClusterWLIntegMethod.LNINT, False),
    (Nc.DataClusterWLIntegMethod.FIXED_NODES, False),
    (Nc.DataClusterWLIntegMethod.FIXED_NODES, True),
    (Nc.DataClusterWLIntegMethod.CUBATURE, False),
]


def _mset() -> Ncm.MSet:
    """The models of the likelihood."""
    cosmo = Nc.HICosmoDEXcdm.new()
    dist = Nc.Distance.new(100.0)
    hms = Nc.HaloCMParam.new(Nc.HaloMassSummaryMassDef.MEAN, 200.0)
    hms.param_set_by_name("log10MDelta", 14.0)
    dp = Nc.HaloDensityProfileNFW.new(hms)
    hp = Nc.HaloPosition.new(dist)
    smd = Nc.WLSurfaceMassDensity.new(dist)
    hp.param_set_by_name("z", 0.2)
    hp.prepare(cosmo)
    pop_shape = Nc.GalaxyShapePopGauss.new()
    pop_shape.param_set_by_name("sigma", 0.3)
    pop_z = Nc.GalaxyRedshiftPopLSSTSRD.new_y1_source()
    obs_z = Nc.GalaxyRedshiftObsGauss.new()
    mset = Ncm.MSet.empty_new()
    for model in (cosmo, dp, hp, smd, pop_shape, pop_z, obs_z):
        mset.set(model)
    mset.prepare_fparam_map()
    return mset


def _data(rows: list[int], mset: Ncm.MSet, method=METHODS[0]) -> Nc.DataClusterWLFactor:
    """A likelihood of the galaxies rows of GALAXIES, as a catalog of those rows alone."""
    position_factor = Nc.GalaxyPositionFactorFlat.new(-0.2, 0.2, -0.2, 0.2)
    redshift_factor = Nc.GalaxyRedshiftFactorComposed.new(0.0, 5.0)
    shape_factor = Nc.GalaxyShapeFactorVarAdd.new(ELLIP_CONV)
    pos_data = Nc.GalaxyPositionFactorData.new(position_factor, mset)
    z_data = Nc.GalaxyRedshiftFactorData.new(redshift_factor, mset)
    s_data = Nc.GalaxyShapeFactorData.new(shape_factor, mset, pos_data, z_data)
    cols = Nc.GalaxyShapeFactorData.required_columns(s_data)

    obs = Nc.GalaxyWLObs.new(ELLIP_CONV, FRAME, len(rows), cols)
    for i, row in enumerate(rows):
        ra, dec, zp, sigma0, e1, e2, std_noise = GALAXIES[row]
        for name, value in (
            ("ra", ra),
            ("dec", dec),
            ("z", 0.0),
            ("zp", zp),
            ("sigma0", sigma0),
            ("epsilon_int_1", 0.0),
            ("epsilon_int_2", 0.0),
            ("epsilon_obs_1", e1),
            ("epsilon_obs_2", e2),
            ("std_noise", std_noise),
            ("c1", 0.0),
            ("c2", 0.0),
            ("m", 0.0),
        ):
            obs.set(name, i, value)

    data = Nc.DataClusterWLFactor.new(position_factor, redshift_factor, shape_factor)
    data.set_integ_method(method[0])
    data.set_auto_nodes(method[1])
    data.set_obs(obs)
    return data


ALL = list(range(len(GALAXIES)))


@pytest.mark.parametrize("method", METHODS)
def test_mask_equals_subset_catalog(method):
    """The likelihood of the masked catalog, and of each of its galaxies, is that of the
    catalog of the rows of the mask, bit for bit."""
    mset = _mset()
    masked = _data(ALL, mset, method)
    masked.set_mask(MASK)
    subset = _data(MASK, mset, method)

    assert masked.get_length() == len(MASK)
    assert masked.props.len == len(MASK)
    assert list(masked.get_mask()) == MASK
    assert masked.m2lnL_val(mset) == subset.m2lnL_val(mset)

    gal_masked = Ncm.Vector.new(len(MASK))
    gal_subset = Ncm.Vector.new(len(MASK))
    masked.eval_m2lnP_gal(mset, gal_masked)
    subset.eval_m2lnP_gal(mset, gal_subset)
    assert gal_masked.dup_array() == gal_subset.dup_array()

    # Removing the mask brings back every row
    masked.set_mask(None)
    assert not masked.get_mask()  # NULL reaches Python as []
    assert masked.get_length() == len(ALL)
    assert masked.m2lnL_val(mset) == _data(ALL, mset, method).m2lnL_val(mset)


def test_mask_serialized():
    """The mask is a property of the likelihood, kept by serialization."""
    mset = _mset()
    masked = _data(ALL, mset)
    masked.set_mask(MASK)
    dup = duplicate_via_serialization(masked)
    assert list(dup.get_mask()) == MASK
    assert dup.peek_obs().len() == len(ALL)
    assert dup.m2lnL_val(mset) == masked.m2lnL_val(mset)

    unmasked = duplicate_via_serialization(_data(ALL, mset))
    assert not unmasked.get_mask()
    assert unmasked.get_length() == len(ALL)


def test_mask_resample_writes_its_rows_only():
    """A resample draws the galaxies of the mask and writes their rows only."""
    mset = _mset()
    data = _data(ALL, mset)
    data.set_mask(MASK)
    obs = data.peek_obs()
    before = obs.peek_data().to_numpy().copy()

    rng = Ncm.RNG.seeded_new(None, 7)
    data.resample(mset, rng)

    after = obs.peek_data().to_numpy()
    others = [i for i in ALL if i not in MASK]
    np.testing.assert_array_equal(after[others], before[others])
    col = list(obs.peek_columns()).index("epsilon_obs_1")
    assert np.all(after[MASK, col] != before[MASK, col])


def test_mask_bootstrap():
    """A bootstrap resamples the galaxies of the mask, as it would those of the subset."""
    mset = _mset()
    masked = _data(ALL, mset, METHODS[1])
    masked.set_mask(MASK)
    subset = _data(MASK, mset, METHODS[1])
    for data in (masked, subset):
        data.bootstrap_create()
        data.bootstrap_resample(Ncm.RNG.seeded_new(None, 11))
    assert masked.m2lnL_val(mset) == subset.m2lnL_val(mset)


def test_shared_clones_with_other_masks():
    """Clones sharing the catalog, with masks of the same length, never take each other's
    node grid: each has the likelihood of its own rows."""
    mset = _mset()
    method = (Nc.DataClusterWLIntegMethod.FIXED_NODES, False)
    data = _data(ALL, mset, method)
    data.set_mask(MASK)
    dset = Ncm.Dataset.new()
    dset.append_data(data)
    ser = Ncm.Serialize.new(Ncm.SerializeOpt.KEEP_SHARED)
    dset.register_shared(ser)
    clone = duplicate_via_serialization(data, ser)
    assert clone.peek_obs() is data.peek_obs()
    other = [1, 2, 4, 6, 7]
    clone.set_mask(other)

    assert data.m2lnL_val(mset) == _data(MASK, mset, method).m2lnL_val(mset)
    assert clone.m2lnL_val(mset) == _data(other, mset, method).m2lnL_val(mset)
    assert data.m2lnL_val(mset) == _data(MASK, mset, method).m2lnL_val(mset)


_ABORT = """
import sys
from numcosmo_py import Nc, Ncm
sys.path.insert(0, sys.argv[1])
import test_data_cluster_wl_factor_mask as t
mset = t._mset()
data = t._data(t.ALL, mset)
data.set_mask([int(x) for x in sys.argv[2].split(",")])
data.m2lnL_val(mset)
"""


@pytest.mark.parametrize(
    "mask, message",
    [
        ("0,3,3", "strictly increasing"),
        ("2,1", "strictly increasing"),
        ("0,8", "row 8"),
    ],
)
def test_bad_mask_aborts(mask, message):
    """A mask that is not strictly increasing, or keeps a row the catalog does not have,
    aborts with the reason."""
    run = subprocess.run(
        [
            sys.executable,
            "-c",
            _ABORT,
            str(__import__("pathlib").Path(__file__).parent),
            mask,
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert run.returncode != 0
    assert message in run.stderr + run.stdout

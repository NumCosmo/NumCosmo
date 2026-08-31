#!/usr/bin/env python
#
# test_bbn_migration.py
#
# Sat August 30 11:00:00 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_bbn_migration.py
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
# with this program.  If not, see <http://www.gnu.org/licenses/>.

"""Tests that cosmologies written before NcBBN existed still load correctly."""

import json
import subprocess
import sys

import numpy as np
import pytest
from numpy.testing import assert_allclose

from gi.repository import GLib

from numcosmo_py import Nc, Ncm

Ncm.cfg_init()

# Yp used to be a parameter of the cosmology, with its fit type doubling as the
# switch between "use this value" and "predict it from BBN". The fixtures in
# data/truth_tables/bbn were written by that code; golden.json records what each
# one evaluated to then. The compat properties are inert sinks now: a fixed-Yp
# file lands on the default NcBBNParthenope (the same physics its fit type
# meant), while a free-Yp file requests a removed mode and must fail loudly --
# never load with silently different physics.
#
# Only lcdm_yp changes value, and deliberately: NcHICosmoLCDM had no BBN branch
# at all, so a fixed Yp was taken at face value. It now gets the same PArthENoPE
# prediction as everything else, which is the point of the migration.
EXPECTED = {
    "de_xcdm_yp_fixed": ("NcBBNParthenope", None),
    "lcdm_yp": ("NcBBNParthenope", 0.2452620852),
    "mset_de_xcdm": ("NcBBNParthenope", None),
}

FORMATS = ("obj", "bin", "yaml")


def load_golden():
    """The values each fixture evaluated to before the migration."""
    filename = Ncm.cfg_get_data_filename("truth_tables/bbn/golden.json", True)

    with open(filename, "r", encoding="utf-8") as f:
        return json.load(f)


def read_fixture(name, fmt):
    """Deserialize one fixture, returning its cosmology."""
    ser = Ncm.Serialize.new(0)
    filename = Ncm.cfg_get_data_filename(f"truth_tables/bbn/{name}.{fmt}", True)

    readers = {
        "obj": ser.from_file,
        "bin": ser.from_binfile,
        "yaml": ser.from_yaml_file,
    }
    obj = readers[fmt](filename)

    if isinstance(obj, Ncm.MSet):
        return obj.peek(Nc.HICosmo.id())

    return obj


@pytest.mark.parametrize("name", sorted(EXPECTED))
@pytest.mark.parametrize("fmt", FORMATS)
def test_old_file_gets_the_right_bbn_model(name, fmt):
    """A pre-NcBBN file lands on the NcBBN its Yp fit type implied."""
    cosmo = read_fixture(name, fmt)
    bbn = cosmo.peek_bbn()

    assert bbn is not None
    assert bbn.__gtype__.name == EXPECTED[name][0]


@pytest.mark.parametrize("name", sorted(EXPECTED))
@pytest.mark.parametrize("fmt", FORMATS)
def test_old_file_still_evaluates_the_same(name, fmt):
    """Yp is unchanged, except where the migration deliberately changes it."""
    golden = load_golden()[name]
    cosmo = read_fixture(name, fmt)

    expected = EXPECTED[name][1]

    if expected is None:
        expected = golden["Yp_4He"]

    assert_allclose(Nc.HICosmo.Yp_4He(cosmo), expected, rtol=1.0e-9)


@pytest.mark.parametrize("fmt", FORMATS)
def test_a_free_yp_file_fails_loudly(fmt):
    """The sampled-Yp compat mode was removed: loading such a file is fatal.

    The refusal is a g_error (the property setter has no error channel), so it
    has to be observed from a subprocess.
    """
    filename = Ncm.cfg_get_data_filename(
        f"truth_tables/bbn/de_xcdm_yp_free.{fmt}", True
    )
    script = (
        "from numcosmo_py import Ncm\n"
        "Ncm.cfg_init()\n"
        "ser = Ncm.Serialize.new(0)\n"
        f"readers = {{'obj': ser.from_file, 'bin': ser.from_binfile, 'yaml': ser.from_yaml_file}}\n"
        f"readers[{fmt!r}]({filename!r})\n"
    )
    result = subprocess.run(
        [sys.executable, "-c", script], capture_output=True, text=True, check=False
    )

    assert result.returncode != 0
    assert "removed sampled-Yp" in result.stderr


@pytest.mark.parametrize("fmt", FORMATS)
def test_other_submodels_survive(fmt):
    """The default NcBBN must not displace submodels the file carried."""
    cosmo = read_fixture("mset_de_xcdm", fmt)
    names = sorted(
        cosmo.peek_submodel(i).__gtype__.name for i in range(cosmo.get_submodel_len())
    )

    assert names == ["NcBBNParthenope", "NcHIPrimPowerLaw", "NcHIReionCamb"]


def test_a_new_cosmology_does_not_serialize_yp():
    """The deprecated properties are write-only, so they do not come back out."""
    ser = Ncm.Serialize.new(0)
    cosmo = Nc.HICosmoDEXcdm.new()

    assert "'Yp'" not in ser.to_string(cosmo, True)


def test_a_stale_reparam_still_attaches():
    """A reparametrization written before Yp left is one slot too long.

    NcHICosmoDE went from eight scalar parameters to seven, so every stored
    NcHICosmoDEReparamOk says length = 8. It carries no parameter *values* --
    only its length, its descriptors and its compatible type -- so the model
    rebuilds it at its own length instead of adopting a vector one too long,
    which used to abort inside ncm_vector_memcpy().
    """
    cosmo = Nc.HICosmoDEXcdm.new()

    assert cosmo.len() == 7

    cosmo.set_reparam(Nc.HICosmoDEReparamOk.new(8))

    assert cosmo.param_name(Nc.HICosmoDESParams.OMEGA_X) == "Omegak"
    assert cosmo.len() == 7


def test_a_reparam_of_a_gone_parameter_fails_loudly():
    """Resizing is only safe while the descriptors still fit.

    A descriptor pointing past the end of the model means the parameter it
    reparametrizes is gone, which no rebuild can repair -- so it must raise
    rather than abort or silently reparametrize the wrong slot.
    """
    cosmo = Nc.HICosmoDEXcdm.new()
    reparam = Nc.HICosmoDEReparamOk.new(8)

    # Slot 7 existed while the cosmology had eight parameters and does not now.
    reparam.set_param_desc_full(
        7, "Bogus", "B", 0.0, 1.0, 0.1, 0.0, 0.5, Ncm.ParamType.FIXED
    )

    with pytest.raises(GLib.Error, match="no longer exists"):
        cosmo.set_reparam(reparam)


# The three containers a catalog has ever used for its model-set, all written by
# tests/tools/make_reparam_catalog_fixtures.py with a pre-migration NcHICosmoDE:
# eight scalar parameters, an NcHICosmoDEReparamOk of length 8, a Yp property and
# no NcBBN submodel. `sidecar' keeps its model-set in a `.mset' GKeyFile beside
# the FITS instead of in HDU0.
CATALOGS = ("vardict", "object", "sidecar")

CATALOG_NROWS = 40
CATALOG_NCOLS = 3


def read_catalog(name):
    """Open one pre-migration catalog fixture read-only."""
    filename = Ncm.cfg_get_data_filename(
        f"truth_tables/bbn/old_reparam_{name}.mc.fits", True
    )

    return Ncm.MSetCatalog.new_from_file_ro(filename, 0)


@pytest.mark.parametrize("name", CATALOGS)
def test_an_old_catalog_still_opens(name):
    """Reading one used to abort inside ncm_vector_memcpy(): 8 into 7."""
    cosmo = read_catalog(name).peek_mset().peek(Nc.HICosmo.id())

    assert cosmo.__gtype__.name == "NcHICosmoDEXcdm"
    assert cosmo.len() == 7


@pytest.mark.parametrize("name", CATALOGS)
def test_an_old_catalog_keeps_its_reparametrization(name):
    """The stale length is bookkeeping; the reparametrization itself survives.

    NcHICosmoDEReparamOk names index 2 Omegak, and Yp was index 4, so nothing it
    names moved when Yp left. The value has to agree with the Omegax the file
    stores, which is what proves old2new() ran on the rebuilt reparametrization.
    """
    mset = read_catalog(name).peek_mset()
    cosmo = mset.peek(Nc.HICosmo.id())

    assert cosmo.peek_reparam() is not None
    assert cosmo.param_name(Nc.HICosmoDESParams.OMEGA_X) == "Omegak"

    expected = (
        1.0
        - cosmo.orig_param_get(Nc.HICosmoDESParams.OMEGA_X)
        - (
            cosmo.orig_param_get(Nc.HICosmoDESParams.OMEGA_C)
            + cosmo.orig_param_get(Nc.HICosmoDESParams.OMEGA_B)
        )
    )
    assert_allclose(cosmo.param_get(Nc.HICosmoDESParams.OMEGA_X), expected, atol=1.0e-3)


@pytest.mark.parametrize("name", CATALOGS)
def test_an_old_catalog_maps_its_free_parameters(name):
    """A catalog is unusable unless its columns still name its free parameters."""
    mcat = read_catalog(name)
    mset = mcat.peek_mset()

    assert mcat.len() == CATALOG_NROWS
    assert mcat.ncols() == CATALOG_NCOLS
    assert [mset.fparam_full_name(i) for i in range(mset.fparam_len())] == [
        "NcHICosmo:Omegak",
        "NcHICosmo:w",
    ]


@pytest.mark.parametrize("name", CATALOGS)
def test_an_old_catalog_gets_the_default_bbn(name):
    """The Yp the file carries is a fixed one, so it lands on PArthENoPE."""
    cosmo = read_catalog(name).peek_mset().peek(Nc.HICosmo.id())
    bbn = cosmo.peek_bbn()

    assert bbn is not None
    assert bbn.__gtype__.name == "NcBBNParthenope"


@pytest.mark.parametrize("name", CATALOGS)
def test_an_old_catalog_still_has_its_rows(name):
    """Loading the model-set must not disturb the chain itself."""
    mcat = read_catalog(name)
    rows = np.array([mcat.peek_row(i).dup_array() for i in range(mcat.len())])

    assert rows.shape == (CATALOG_NROWS, CATALOG_NCOLS)
    assert np.all(np.isfinite(rows))
    # Column 0 is m2lnL = Omegak^2 + (w + 1)^2, the fixture's stand-in likelihood.
    assert_allclose(rows[:, 0], rows[:, 1] ** 2 + (rows[:, 2] + 1.0) ** 2)


def test_a_catalog_whose_reparam_cannot_be_resized_fails_loudly():
    """Deserialization has no error channel, so the refusal is a g_error.

    ncm_model_set_reparam() is reached through the `reparam' property, which
    cannot report a GError -- the refusal must still be a clean abort with the
    reason in it, never a NULL dereference.
    """
    script = (
        "from numcosmo_py import Ncm\n"
        "Ncm.cfg_init()\n"
        "ser = Ncm.Serialize.new(0)\n"
        "matrix = Ncm.Matrix.new(3, 3)\n"
        "matrix.set_identity()\n"
        "vector = Ncm.Vector.new(3)\n"
        "vector.set_all(0.0)\n"
        "reparam = Ncm.ReparamLinear.new(3, matrix, vector)\n"
        "reparam.set_compat_type('NcmModelRosenbrock')\n"
        "text = ser.to_string(Ncm.ModelRosenbrock.new(), True)\n"
        "ser.from_string(\n"
        "    text[:-2] + \", 'reparam': <\" + ser.to_string(reparam, True) + '>})'\n"
        ")\n"
    )
    result = subprocess.run(
        [sys.executable, "-c", script], capture_output=True, text=True, check=False
    )

    assert result.returncode != 0
    assert "carries state of its own" in result.stderr

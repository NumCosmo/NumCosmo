#!/usr/bin/env python
#
# test_planck_simall.py
#
# Wed July 23 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_planck_simall.py
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

"""Tests on the native NcDataPlanckSimall (low-ell SimAll) likelihood."""

import os

import numpy as np
import pytest
from astropy.io import fits

from python.fixtures_planck import (
    SIMALL_LMIN,
    SIMALL_LMAX,
    SIMALL_STEP,
    FixedClBoltzmann,
    fixed_spectra,
    make_simall_cldf,
    planck_mset,
    theory_cls,
)
from numcosmo_py import Ncm, Nc
from numcosmo_py.experiments.planck_simall import (
    SIMALL_EE_RELPATH,
    SIMALL_BB_RELPATH,
    SIMALL_EEBB_RELPATH,
    build_simall,
)
from numcosmo_py.experiments.planck_lite import find_baseline_file

Ncm.cfg_init()

_EE = find_baseline_file(SIMALL_EE_RELPATH)
needs_ee = pytest.mark.skipif(_EE is None, reason="simall EE data not found")


def test_type_is_ncmdata():
    """The native type is an NcmData subclass."""
    assert issubclass(Nc.DataPlanckSimall, Ncm.Data)


@needs_ee
def test_construct_and_serialize_roundtrip():
    """Build from the simall EE clik data, serialize and reload intact."""
    simall = build_simall(_EE)
    assert simall.get_length() == 28

    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
    dset = Ncm.Dataset.new_array([simall])
    simall2 = ser.from_variant(ser.to_variant(dset)).peek_data(0)

    assert isinstance(simall2, Nc.DataPlanckSimall)
    assert simall2.get_length() == 28


@pytest.mark.planck_data
@pytest.mark.parametrize(
    "relpath,enum,length",
    [
        (SIMALL_EE_RELPATH, "BASELINE_18_LOWL_EE", 28),
        (SIMALL_BB_RELPATH, "BASELINE_18_LOWL_BB", 28),
        (SIMALL_EEBB_RELPATH, "BASELINE_18_LOWL_EEBB", 56),
    ],
)
def test_matches_clik_reference(relpath, enum, length):
    """Native simall m2lnL matches the clik low-ell reference (bit-exact).

    The tabulated log-probability with floor indexing reproduces clik exactly.
    """
    # pylint: disable=import-outside-toplevel
    from numcosmo_py.cosmology import create_cosmo, HIPrimModel
    from numcosmo_py.experiments.planck18 import mset_set_parameters, Planck18Types

    clik = find_baseline_file(relpath)
    if clik is None:
        pytest.skip(f"simall data not found ({relpath})")

    cbe = Nc.HIPertBoltzmannCBE.new()
    cosmo = create_cosmo(prim_model=HIPrimModel.POWER_LAW)
    planck = Nc.PlanckFICorTT()
    planck.params_set_default_ftype()
    mset = Ncm.MSet.new_array([planck, cosmo])
    mset_set_parameters(mset, Planck18Types.TT, HIPrimModel.POWER_LAW)
    mset.prepare_fparam_map()

    ref = Nc.DataPlanckLKL.full_new_id(getattr(Nc.DataPlanckLKLType, enum), cbe)
    native = build_simall(clik, cbe)
    assert native.get_length() == length

    ref.prepare(mset)
    native.prepare(mset)
    assert native.m2lnL_val(mset) == pytest.approx(ref.m2lnL_val(mset), rel=1.0e-13)


# -----------------------------------------------------------------------------
# Synthetic-data tests: no plc_3.0 tree needed, so they also run where the real
# Planck data is absent (CI included). See tests/python/fixtures_planck.py.
# -----------------------------------------------------------------------------


def test_synthetic_matches_table_lookup(tmp_path):
    """The native m2lnL is the tabulated log-probability sum, to the last bit.

    Builds a synthetic simall cldf tree, runs it through the production
    converter and checks the whole chain (mdb/has_cl parsing, table layout,
    Dl conversion, calibration and floor indexing) against a direct numpy
    lookup in the same table.
    """
    clik = make_simall_cldf(tmp_path, spectra=("EE", "BB"))
    pb = FixedClBoltzmann()
    calib = 1.003
    mset, _ = planck_mset(A_planck=calib)

    simall = build_simall(clik, pb)
    assert simall.get_length() == 2 * (SIMALL_LMAX - SIMALL_LMIN + 1)

    simall.prepare(mset)

    ell = np.arange(SIMALL_LMIN, SIMALL_LMAX + 1)
    nell = ell.size
    lnl = 0.0
    for tag in ("EE", "BB"):
        table = fits.getdata(os.path.join(clik, "clik", "lkl_0", f"prob{tag}"))
        cl = theory_cls(pb, tag, SIMALL_LMAX)[SIMALL_LMIN:]
        dl = cl / calib**2 * ell * (ell + 1.0) / (2.0 * np.pi)
        position = (dl / SIMALL_STEP).astype(int)
        assert np.all((position >= 0) & (position < table.shape[1]))
        # Accumulate the way the C does -- one running sum per spectrum, then
        # added together -- so "to the last bit" holds by construction instead
        # of depending on numpy's pairwise summation order.
        res = 0.0
        for i in range(nell):
            res += float(table[i, position[i]])
        lnl += res

    assert simall.m2lnL_val(mset) == -2.0 * lnl


def test_synthetic_out_of_range_table(tmp_path):
    """A theory Dl past the end of the table returns the rejected-point value."""
    # One step of 1e-12 puts every low-l EE bandpower beyond the last column.
    clik = make_simall_cldf(tmp_path, spectra=("EE",), nsteps=4, step=1.0e-12)
    pb = FixedClBoltzmann()
    mset, _ = planck_mset()

    simall = build_simall(clik, pb)
    simall.prepare(mset)

    assert simall.m2lnL_val(mset) == 1.0e30


def test_synthetic_te_spectrum(tmp_path):
    """A file with a TE table is read and evaluated through the TE code path."""
    clik = make_simall_cldf(tmp_path, spectra=("EE", "TE"))
    # The rejection below needs a TE that is below the tabulated range; ask for
    # it explicitly rather than relying on the stand-in spectrum's sign.
    spectra = fixed_spectra()
    spectra["TE"] = -np.abs(spectra["TE"]) - 1.0e-12
    pb = FixedClBoltzmann(spectra)
    mset, _ = planck_mset()

    simall = build_simall(clik, pb)
    assert simall.get_property("step-te") == SIMALL_STEP
    assert simall.get_property("prob-te") is not None

    simall.prepare(mset)
    # A negative TE is below the tabulated range: the likelihood must reject
    # the point rather than index the table out of bounds.
    assert simall.m2lnL_val(mset) == 1.0e30


def test_synthetic_serialize_roundtrip(tmp_path):
    """A synthetic simall survives a serialization round trip unchanged."""
    clik = make_simall_cldf(tmp_path, spectra=("EE", "BB"))
    pb = FixedClBoltzmann()
    mset, _ = planck_mset()

    simall = build_simall(clik, pb)
    simall.prepare(mset)
    m2lnl = simall.m2lnL_val(mset)

    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
    simall2 = ser.from_variant(ser.to_variant(simall))
    assert isinstance(simall2, Nc.DataPlanckSimall)
    assert simall2.get_length() == simall.get_length()

    simall2.set_hipert_boltzmann(pb)
    simall2.prepare(mset)
    assert simall2.m2lnL_val(mset) == m2lnl


def test_synthetic_data_only_build(tmp_path):
    """Built without a Cls source, the object holds the data and no Boltzmann.

    This is the form the release artifacts are serialized in; the Boltzmann is
    attached when the block is loaded into an experiment.
    """
    clik = make_simall_cldf(tmp_path)
    simall = build_simall(clik)

    assert simall.peek_hipert_boltzmann() is None
    assert simall.get_length() == SIMALL_LMAX - SIMALL_LMIN + 1


def test_synthetic_new_from_file(tmp_path):
    """The `new_from_file` loader reads back a serialized simall.

    Only the compact likelihoods are checked this way: the loader is shared
    across all five, and the text serialization of the high-l blocks is large.
    """
    clik = make_simall_cldf(tmp_path)
    data = build_simall(clik)

    path = str(tmp_path / "simall.obj")
    Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP).to_file(data, path)

    loaded = Nc.DataPlanckSimall.new_from_file(path)
    assert isinstance(loaded, Nc.DataPlanckSimall)
    assert loaded.get_length() == data.get_length()

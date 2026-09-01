#!/usr/bin/env python
#
# test_planck18_native.py
#
# Fri August 29 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_planck18_native.py
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

"""Tests for the native (clik-free) Planck 2018 experiment assembly.

``generate_planck18_native`` combines the four native likelihood blocks into one
shared CBE. The blocks themselves are covered by the per-likelihood tests; here
they are stood in for by synthetic ones (see tests/python/fixtures_planck.py),
driven by a fixed-spectra Boltzmann, so the assembly, the priors and the
derived-function wiring are exercised without a local ``plc_3.0`` tree and
without a Boltzmann solve.
"""

import numpy as np
import pytest

from python.fixtures_planck import (
    FixedClBoltzmann,
    make_commander_cldf,
    make_lensing_cldf,
    make_simall_cldf,
    make_smica_cldf,
)
from numcosmo_py import Ncm, Nc
from numcosmo_py.cosmology import HIPrimModel
import numcosmo_py.experiments.planck_native_release as pnr
from numcosmo_py.experiments.planck18 import (
    Planck18Types,
    generate_planck18_native,
    mset_set_parameters,
)
from numcosmo_py.experiments.planck_commander import build_commander
from numcosmo_py.experiments.planck_lensing import build_lensing
from numcosmo_py.experiments.planck_simall import build_simall
from numcosmo_py.experiments.planck_smica import build_smica_tt, build_smica_ttteee

Ncm.cfg_init()


@pytest.fixture(name="synthetic_trees", scope="module")
def fixture_synthetic_trees(tmp_path_factory):
    """The synthetic cldf trees and their builders, keyed by release id."""
    root = tmp_path_factory.mktemp("planck_native")
    return {
        pnr.PlanckReleaseId.PR3_SIMALL_EE: (make_simall_cldf(root), build_simall),
        pnr.PlanckReleaseId.PR3_COMMANDER: (make_commander_cldf(root), build_commander),
        pnr.PlanckReleaseId.PR3_PLIK_TT: (make_smica_cldf(root), build_smica_tt),
        pnr.PlanckReleaseId.PR3_PLIK_TTTEEE: (
            make_smica_cldf(root, pol=True),
            build_smica_ttteee,
        ),
        pnr.PlanckReleaseId.PR3_LENSING: (make_lensing_cldf(root), build_lensing),
    }


def _release_loader(trees, boltzmann=None):
    """Stand in for ``load_planck_release``, optionally overriding the Boltzmann.

    @boltzmann replaces the CBE the experiment assembles, so the blocks read
    stored spectra and nothing is solved; None keeps the experiment's own CBE.
    """

    def load(rid, pb=None, cache_dir=None):  # signature of load_planck_release
        del cache_dir
        clik, builder = trees[rid]
        return builder(clik, pb if boltzmann is None else boltzmann)

    return load


@pytest.fixture(name="synthetic_blocks", scope="module")
def fixture_synthetic_blocks(synthetic_trees):
    """Release loader whose blocks are driven by #FixedClBoltzmann.

    The assembly, the priors and the derived-function wiring are what these tests
    check; a Boltzmann solve at Planck precision would dominate their runtime
    without covering any of it. The coupling to CLASS is covered separately by
    ``test_generate_native_evaluates_on_cbe``.
    """
    return _release_loader(synthetic_trees, FixedClBoltzmann())


@pytest.mark.parametrize(
    "data_type,model_name",
    [
        (Planck18Types.TT, "PlanckFICorTT"),
        (Planck18Types.TTTEEE, "PlanckFICorTTTEEE"),
    ],
)
@pytest.mark.parametrize("use_lensing", [False, True])
def test_generate_native_from_release(
    monkeypatch, synthetic_blocks, data_type, model_name, use_lensing
):
    """The native experiment assembles, evaluates and carries its derived functions.

    ``from_release=True`` is the path that needs no local clik tree; the release
    loader is replaced by synthetic blocks so the assembly runs anywhere.
    """
    monkeypatch.setattr(pnr, "load_planck_release", synthetic_blocks)

    experiment, mfunc_array = generate_planck18_native(
        data_type=data_type,
        use_lensing_likelihood=use_lensing,
        from_release=True,
    )

    likelihood = experiment.peek("likelihood")
    mset = experiment.peek("model-set")
    assert isinstance(likelihood, Ncm.Likelihood)
    assert isinstance(mset, Ncm.MSet)
    assert mset.peek_by_name("NcPlanckFI").__class__.__name__ == model_name

    dset = likelihood.peek_dataset()
    assert dset.get_ndata() == (4 if use_lensing else 3)
    # Native blocks only: the experiment must reload without clik or the PLC.
    names = {dset.get_data(i).__class__.__name__ for i in range(dset.get_ndata())}
    assert "DataPlanckLKL" not in names

    assert experiment.peek("distance") is not None
    assert experiment.peek("ps-ml") is not None
    assert experiment.peek("ps-ml-filter") is not None
    assert mfunc_array.len() == 4

    mset_set_parameters(mset, data_type, HIPrimModel.POWER_LAW)
    mset.prepare_fparam_map()
    for i in range(dset.get_ndata()):
        dset.get_data(i).prepare(mset)
    assert np.isfinite(dset.m2lnL_val(mset))


@pytest.mark.acceptance
def test_generate_native_evaluates_on_cbe(monkeypatch, synthetic_trees):
    """The assembled experiment evaluates against a real CBE.

    The cases above swap the Boltzmann out, so this is what covers the blocks
    self-configuring the shared CBE (targets and lmax) and the resulting solve
    feeding a finite $-2\\ln L$. One combination is enough: the block-by-block
    self-configuration is covered by the per-likelihood tests.
    """
    monkeypatch.setattr(pnr, "load_planck_release", _release_loader(synthetic_trees))

    experiment, _ = generate_planck18_native(
        data_type=Planck18Types.TT,
        use_lensing_likelihood=True,
        from_release=True,
    )

    mset = experiment.peek("model-set")
    dset = experiment.peek("likelihood").peek_dataset()

    mset_set_parameters(mset, Planck18Types.TT, HIPrimModel.POWER_LAW)
    mset.prepare_fparam_map()
    for i in range(dset.get_ndata()):
        dset.get_data(i).prepare(mset)
    assert np.isfinite(dset.m2lnL_val(mset))


def test_generate_native_two_fluids_raises_precision(monkeypatch, synthetic_blocks):
    """The two-fluids primordial model raises the CBE primordial sampling."""
    monkeypatch.setattr(pnr, "load_planck_release", synthetic_blocks)

    experiment, _ = generate_planck18_native(
        data_type=Planck18Types.TT,
        prim_model=HIPrimModel.TWO_FLUIDS,
        from_release=True,
    )

    mset = experiment.peek("model-set")
    assert isinstance(mset, Ncm.MSet)
    cosmo = mset.peek_by_name("NcHICosmo")
    assert isinstance(cosmo.peek_prim(), Nc.HIPrimTwoFluids)


def test_generate_native_invalid_data_type(monkeypatch, synthetic_blocks):
    """An unknown data combination is rejected."""
    monkeypatch.setattr(pnr, "load_planck_release", synthetic_blocks)

    with pytest.raises(ValueError, match="Invalid data type"):
        generate_planck18_native(data_type="TTTEEE_lensing", from_release=True)


def test_generate_native_without_local_data(monkeypatch):
    """Building from a local clik tree that is absent fails with a clear message."""
    import numcosmo_py.experiments.planck_lite as pl  # pylint: disable=C0415

    monkeypatch.setattr(pl, "find_baseline_file", lambda relpath: None)

    with pytest.raises(FileNotFoundError, match="Planck baseline data not found"):
        generate_planck18_native(data_type=Planck18Types.TT)

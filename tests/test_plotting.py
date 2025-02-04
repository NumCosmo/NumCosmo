#!/usr/bin/env python
#
# test_plotting.py
#
# Mon Feb 03 17:37:17 2025
# Copyright  2025  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_py_serialize.py
# Copyright (C) 2025 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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
#

"""Tests for the plotting module."""

import re
import pytest
import numpy as np
from numpy.testing import assert_allclose
from getdist import MCSamples

from numcosmo_py import Ncm
from numcosmo_py.plotting import mcat_to_catalog_data, CatalogData

Ncm.cfg_init()


@pytest.fixture(name="mcat", params=[False, True], ids=["unweighted", "weighted"])
def fixture_mcat(request):
    """Return a cosmology object."""
    model = Ncm.ModelMVND.new(dim=10)
    mset = Ncm.MSet.new_array([model])
    mset.param_set_all_ftype(Ncm.ParamType.FREE)
    mset.prepare_fparam_map()
    nchains = 34

    mcat = Ncm.MSetCatalog(
        mset=mset,
        nadd_vals=2,
        nadd_val_names=["a", "b"],
        nadd_val_symbols=["c", "d"],
        nchains=nchains,
        m2lnp_var=0,
        weighted=request.param,
    )

    nsamples = 100 * nchains
    samples = np.random.rand(nsamples, 10) - 10.0
    add_vals = np.random.rand(nsamples, 3)
    for sample, add_val in zip(samples, add_vals):
        v = Ncm.Vector.new_array(sample)
        mset.param_set_vector(v)
        mcat.add_from_mset_array(mset, add_val)

    return mcat


def test_mcat_catalog_data(mcat):
    """Test log_current_stats()."""
    cd = mcat_to_catalog_data(mcat, "test")
    assert cd is not None
    assert isinstance(cd, CatalogData)

    # Check that the number of columns is correct, the posterior column is excluded
    assert cd.rows.shape[1] == 12 - 1
    assert cd.nchains == mcat.nchains()
    assert cd.rows.shape[0] == mcat.len()

    m2lnp_index = mcat.get_m2lnp_var()
    excluded_indices = [m2lnp_index]
    if mcat.weighted():
        excluded_indices.append(mcat.nadd_vals() - 1)

    # Check that the parameters are correct, the posterior column is excluded
    assert cd.params_names == [
        mcat.col_name(i) for i in range(mcat.ncols()) if i not in excluded_indices
    ]
    assert cd.params_symbols == [
        mcat.col_symb(i) for i in range(mcat.ncols()) if i not in excluded_indices
    ]

    rows = np.array([mcat.peek_row(i).dup_array() for i in range(mcat.len())])

    # Check that the posterior is correct (should be -\ln p)
    assert_allclose(cd.posterior, 0.5 * rows[:, m2lnp_index])
    assert_allclose(cd.rows, np.delete(rows, excluded_indices, axis=1))


def test_catalog_data_to_mcsamples(mcat):
    """Test log_current_stats()."""
    cd = mcat_to_catalog_data(mcat, "test")
    mcsample = cd.to_mcsamples()
    assert mcsample is not None
    assert isinstance(mcsample, MCSamples)

    assert mcsample.numrows == cd.rows.shape[0]
    assert mcsample.n == cd.rows.shape[1]
    for i, offset in enumerate(mcsample.chain_offsets):
        assert offset == i * 100
        if offset >= cd.rows.shape[0]:
            continue
        # We save it interweaved and MCSamples gets all chains and stacks them
        assert_allclose(
            mcsample.samples[offset : offset + 100, :], cd.rows[i :: cd.nchains, :]
        )


def test_catalog_data_bad_burnin(mcat):
    """Test the corret warning."""
    with pytest.warns(
        UserWarning,
        match=re.escape(
            "Burnin (13) is not a multiple of nchains (34). "
            "Burnin will be rounded down."
        ),
    ):
        _ = mcat_to_catalog_data(mcat, "test", burnin=13)


def test_catalog_data_too_many_burnin(mcat):
    """Test the corret warning."""
    with pytest.raises(ValueError, match="Burnin is greater than the number of steps."):
        _ = mcat_to_catalog_data(mcat, "test", burnin=34 * 10000)


def test_catalog_data_with_indices(mcat):
    """Test log_current_stats()."""
    indices = np.array([4, 6, 8])
    cd = mcat_to_catalog_data(mcat, "test", indices=indices)
    assert cd is not None
    assert isinstance(cd, CatalogData)
    assert cd.rows.shape[1] == 3

    # Check that the parameters are correct, the posterior column is excluded
    assert cd.params_names == [
        mcat.col_name(i) for i in range(mcat.ncols()) if i in indices
    ]
    assert cd.params_symbols == [
        mcat.col_symb(i) for i in range(mcat.ncols()) if i in indices
    ]
    rows = np.array(
        [np.array(mcat.peek_row(i).dup_array())[indices] for i in range(mcat.len())]
    )
    assert_allclose(cd.rows, rows)


def test_catalog_data_to_mcsamples_collapse(mcat):
    """Test log_current_stats()."""
    cd = mcat_to_catalog_data(mcat, "test")
    mcsample = cd.to_mcsamples(collapse=True)
    assert mcsample is not None
    assert isinstance(mcsample, MCSamples)

    assert mcsample.numrows == cd.rows.shape[0]
    assert mcsample.n == cd.rows.shape[1]
    assert_allclose(mcsample.samples, cd.rows)


def test_asinh_transform(mcat):
    """Test asinh_transform()."""
    cd = mcat_to_catalog_data(mcat, "test")
    cd_orig = mcat_to_catalog_data(mcat, "test")
    cd.asinh_transform([3])

    for i in range(cd.rows.shape[1]):
        if i == 3:
            continue
        assert_allclose(cd.rows[:, i], cd_orig.rows[:, i])
        assert cd.params_names[i] == cd_orig.params_names[i]
        assert cd.params_symbols[i] == cd_orig.params_symbols[i]

    assert_allclose(cd.posterior, cd_orig.posterior)
    assert_allclose(cd.rows[:, 3], np.arcsinh(cd_orig.rows[:, 3]))

    assert cd.params_names[3] == "asinh_mu_2"
    assert cd.params_symbols[3] == r"\mathrm{sinh}^{-1}({\mu}_2)"

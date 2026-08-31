#!/usr/bin/env python
#
# test_derived.py
#
# Sun Jul 26 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_derived.py
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

"""Unit tests for numcosmo_py.plotting.derived (corner-plot derived quantities)."""

import numpy as np
import pytest
from numpy.testing import assert_allclose

pytest.importorskip("getdist")
# flake8: noqa: E402
# pylint: disable=wrong-import-position

from numcosmo_py import Ncm
from numcosmo_py.plotting import mcat_to_catalog_data
from numcosmo_py.plotting.derived import add_derived_column

Ncm.cfg_init()

NCHAINS = 5
NSTEPS = 20


@pytest.fixture(name="mset_mcat")
def fixture_mset_mcat():
    """Return a (mset, mcat) pair for a small MVN model with 3 free parameters."""
    model = Ncm.ModelMVND.new(dim=3)
    mset = Ncm.MSet.new_array([model])
    mset.param_set_all_ftype(Ncm.ParamType.FREE)
    mset.prepare_fparam_map()

    mcat = Ncm.MSetCatalog(
        mset=mset,
        nadd_vals=1,
        nadd_val_names=["m2lnL"],
        nadd_val_symbols=["m2lnL"],
        nchains=NCHAINS,
        m2lnp_var=0,
        weighted=False,
    )

    rng = np.random.default_rng(42)
    nsamples = NSTEPS * NCHAINS
    samples = rng.normal(size=(nsamples, 3))
    for sample in samples:
        v = Ncm.Vector.new_array(sample)
        mset.param_set_vector(v)
        mcat.add_from_mset_array(mset, [float(np.sum(sample**2))])

    return mset, mcat


def test_add_derived_column_basic(mset_mcat):
    """The derived column matches evaluating the expression row by row."""
    mset, mcat = mset_mcat
    nadd_vals = mcat.nadd_vals()
    cd = mcat_to_catalog_data(mcat, "test")

    new_cd = add_derived_column(
        cd,
        mset,
        nadd_vals,
        mcat,
        ["x=mu_0", "y=mu_1"],
        "x + y",
        None,
        "derived",
    )

    assert new_cd.params_names == cd.params_names + ["derived"]
    # Symbol defaults to the expression itself.
    assert new_cd.params_symbols == cd.params_symbols + ["x + y"]

    i_x = cd.params_names.index(mcat.col_name(nadd_vals))
    i_y = cd.params_names.index(mcat.col_name(nadd_vals + 1))
    assert_allclose(new_cd.rows[:, -1], cd.rows[:, i_x] + cd.rows[:, i_y])


def test_add_derived_column_custom_symbol_and_name(mset_mcat):
    """--derived-symbol/--derived-name override the defaults."""
    mset, mcat = mset_mcat
    nadd_vals = mcat.nadd_vals()
    cd = mcat_to_catalog_data(mcat, "test")

    new_cd = add_derived_column(
        cd, mset, nadd_vals, mcat, ["x=mu_0"], "10**x", r"10^{x}", "pow10_x"
    )

    assert new_cd.params_names[-1] == "pow10_x"
    assert new_cd.params_symbols[-1] == r"10^{x}"

    i_x = cd.params_names.index(mcat.col_name(nadd_vals))
    assert_allclose(new_cd.rows[:, -1], 10.0 ** cd.rows[:, i_x])


def test_add_derived_column_bestfit(mset_mcat):
    """The best-fit entry is the expression evaluated at the catalog's best fit."""
    mset, mcat = mset_mcat
    nadd_vals = mcat.nadd_vals()
    cd = mcat_to_catalog_data(mcat, "test")
    assert cd.bestfit is not None

    new_cd = add_derived_column(
        cd, mset, nadd_vals, mcat, ["x=mu_0", "y=mu_1"], "x + y", None, "d"
    )

    i_x = cd.params_names.index(mcat.col_name(nadd_vals))
    i_y = cd.params_names.index(mcat.col_name(nadd_vals + 1))
    assert new_cd.bestfit[-1] == pytest.approx(cd.bestfit[i_x] + cd.bestfit[i_y])


def test_add_derived_column_bare_variable_shorthand(mset_mcat):
    """A bare "parameter" binding is shorthand for "parameter=parameter"."""
    mset, mcat = mset_mcat
    nadd_vals = mcat.nadd_vals()
    cd = mcat_to_catalog_data(mcat, "test")

    new_cd = add_derived_column(cd, mset, nadd_vals, mcat, ["mu_0"], "mu_0", None, "d")

    i_x = cd.params_names.index(mcat.col_name(nadd_vals))
    assert_allclose(new_cd.rows[:, -1], cd.rows[:, i_x])


def test_add_derived_column_invalid_binding_syntax(mset_mcat):
    """A --derived-variable string missing a name or a parameter raises ValueError."""
    mset, mcat = mset_mcat
    cd = mcat_to_catalog_data(mcat, "test")

    with pytest.raises(ValueError, match="expected name=parameter"):
        add_derived_column(cd, mset, mcat.nadd_vals(), mcat, ["=mu_0"], "x", None, "d")


def test_add_derived_column_unknown_parameter(mset_mcat):
    """Binding a name to a nonexistent parameter fails clearly."""
    mset, mcat = mset_mcat
    cd = mcat_to_catalog_data(mcat, "test")

    with pytest.raises(ValueError, match="not found"):
        add_derived_column(
            cd, mset, mcat.nadd_vals(), mcat, ["x=does_not_exist"], "x", None, "d"
        )


def test_add_derived_column_excluded_parameter(mset_mcat):
    """A parameter excluded from the catalog selection cannot be bound."""
    mset, mcat = mset_mcat
    nadd_vals = mcat.nadd_vals()
    # Keep only mu_1 and mu_2, dropping the column mu_0 is bound to.
    cd = mcat_to_catalog_data(
        mcat, "test", indices=np.array([nadd_vals + 1, nadd_vals + 2])
    )

    with pytest.raises(ValueError, match="is not part of this catalog's selection"):
        add_derived_column(cd, mset, nadd_vals, mcat, ["x=mu_0"], "x", None, "d")


def test_add_derived_column_unsafe_expression(mset_mcat):
    """Disallowed expression syntax is rejected, not silently evaluated."""
    mset, mcat = mset_mcat
    cd = mcat_to_catalog_data(mcat, "test")

    with pytest.raises(ValueError):
        add_derived_column(
            cd,
            mset,
            mcat.nadd_vals(),
            mcat,
            ["x=mu_0"],
            "__import__('os')",
            None,
            "d",
        )

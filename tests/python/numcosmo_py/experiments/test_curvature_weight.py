#!/usr/bin/env python
#
# test_curvature_weight.py
#
# Sun Jun 22 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_curvature_weight.py
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

"""Unit tests for the data-driven local curvature weight builder."""

import numpy as np
import pytest
from numpy.testing import assert_allclose

from numcosmo_py import Nc, Ncm
from numcosmo_py.datasets.hicosmo import add_snia_likelihood, SNIaID
from numcosmo_py.experiments import curvature_weight as cw

Ncm.cfg_init()


@pytest.fixture(name="wspline_experiment")
def fixture_wspline_experiment():
    """A w-spline + SNIa experiment with the knots and Omegac set free."""
    n_knots, z_max = 8, 2.0
    cosmo = Nc.HICosmoDEWSpline.new(n_knots, z_max)
    cosmo.omega_x2omega_k()
    cosmo["Omegak"] = 0.0
    for i in range(n_knots):
        cosmo.param_set_desc(f"w_{i}", {"fit": True})
    cosmo.param_set_desc("Omegac", {"fit": True})

    mset = Ncm.MSet.new_array([cosmo])
    dist = Nc.Distance.new(10.0)
    dataset = Ncm.Dataset.new()
    add_snia_likelihood(dataset, mset, dist, SNIaID.COV_DES_Y5_STATONLY)
    mset.prepare_fparam_map()
    return cosmo, mset, dataset, z_max


def test_information_weight_values_monotonic() -> None:
    """W(I) decreases with information and stays within [floor, 1]."""
    info = np.array([0.0, 0.1, 1.0, 10.0, 100.0])
    weight = cw.information_weight_values(info, ref_factor=1.0, floor=1e-3)

    assert np.all(weight >= 1e-3) and np.all(weight <= 1.0)
    assert np.all(np.diff(weight) < 0.0)  # strictly decreasing in information
    assert_allclose(weight[0], 1.0)  # zero information -> full prior


def test_marginal_knot_fisher_schur() -> None:
    """The marginal knot Fisher is the Schur complement over the other block."""
    rng = np.random.default_rng(0)
    a = rng.normal(size=(5, 5))
    fisher = a @ a.T + 5.0 * np.eye(5)  # SPD
    knot_idx = [0, 1, 2]
    other_idx = [3, 4]

    f_d = cw.marginal_knot_fisher(fisher, knot_idx)

    f_kk = fisher[np.ix_(knot_idx, knot_idx)]
    f_ko = fisher[np.ix_(knot_idx, other_idx)]
    f_oo = fisher[np.ix_(other_idx, other_idx)]
    expected = f_kk - f_ko @ np.linalg.solve(f_oo, f_ko.T)
    assert_allclose(f_d, expected)
    # Equivalent to the knot block of the inverse, inverted.
    cov = np.linalg.inv(fisher)
    assert_allclose(f_d, np.linalg.inv(cov[np.ix_(knot_idx, knot_idx)]))


def test_wspline_weight_relaxes_where_data_constrains(wspline_experiment) -> None:
    """W is small where SNIa constrain w(z) and rises to 1 where data is blind."""
    cosmo, mset, dataset, _ = wspline_experiment
    weight = cw.wspline_curvature_weight(
        dataset, mset, cosmo, n_grid=32, ref_factor=1.0
    )

    alpha_f = cosmo.get_alpha().dup_array()[-1]
    grid = np.linspace(0.0, alpha_f, 32)
    w_vals = np.array([weight.eval(x) for x in grid])

    assert np.all(w_vals >= 1e-3) and np.all(w_vals <= 1.0 + 1e-12)
    # High-z end is unconstrained by SNIa -> prior near full strength.
    assert w_vals[-1] > 0.9
    # The data-rich low/mid-z interior relaxes the prior well below the blind end.
    assert w_vals[: len(grid) // 2].min() < 0.3
    assert w_vals[: len(grid) // 2].min() < w_vals[-1]


def test_wspline_weight_changes_curvature_norm(wspline_experiment) -> None:
    """The data-driven weight actually localizes the curvature functional."""
    cosmo, mset, dataset, _ = wspline_experiment
    weight = cw.wspline_curvature_weight(
        dataset, mset, cosmo, n_grid=32, ref_factor=1.0
    )

    n_knots = cosmo.get_alpha().len()
    cosmo.props.w = Ncm.Vector.new_array(
        (-1.0 + 0.2 * np.sin(np.linspace(0.0, 3.0, n_knots))).tolist()
    )
    ctype = Ncm.SplineCurvatureType.D2
    plain = cosmo.lp_norm(ctype, 2.0)
    weighted = cosmo.weighted_lp_norm(ctype, 2.0, weight)

    assert np.isfinite(weighted) and weighted > 0.0
    assert not np.isclose(weighted, plain)


def test_wspline_weight_ref_factor_scaling(wspline_experiment) -> None:
    """A larger ref_factor keeps the prior stronger (weights closer to 1)."""
    cosmo, mset, dataset, _ = wspline_experiment
    grid = np.linspace(0.0, cosmo.get_alpha().dup_array()[-1], 32)

    soft = cw.wspline_curvature_weight(dataset, mset, cosmo, n_grid=32, ref_factor=0.1)
    hard = cw.wspline_curvature_weight(dataset, mset, cosmo, n_grid=32, ref_factor=10.0)
    w_soft = np.array([soft.eval(x) for x in grid])
    w_hard = np.array([hard.eval(x) for x in grid])

    assert np.all(w_hard >= w_soft - 1e-9)
    assert w_hard.mean() > w_soft.mean()


def test_wspline_weight_function_serializes(wspline_experiment) -> None:
    """A wlp_* prior carrying the data-driven weight survives serialization."""
    cosmo, mset, dataset, _ = wspline_experiment
    weight = cw.wspline_curvature_weight(
        dataset, mset, cosmo, n_grid=32, ref_factor=1.0
    )
    func = Ncm.MSetFuncList.new("NcHICosmoDEWSpline:wlp_w2", weight)
    direct = func.eval1(mset, 2.0)

    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
    reloaded = ser.from_yaml(ser.to_yaml(func))
    assert_allclose(reloaded.eval1(mset, 2.0), direct)

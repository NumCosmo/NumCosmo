#
# test_spectral_tensor.py
#
# Thu Nov 19 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_spectral_tensor.py
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

"""Tests of the NcmSpectralTensor bindings: arrays, tree, evaluation, serialization.

The numerics are tested in tests/c/ncm/algebra/test_ncm_spectral_tensor.c.
"""

import numpy as np
import pytest
from numpy.polynomial import chebyshev

from numcosmo_py import Ncm

Ncm.cfg_init()

N_POINTS = 200


class Panel:
    """One panel: its node, bounds, numbers of coefficients and coefficients."""

    def __init__(
        self,
        node: int,
        lower: list[float],
        upper: list[float],
        rng: np.random.Generator,
    ) -> None:
        self.node = node
        self.lower = np.array(lower)
        self.upper = np.array(upper)
        self.n_coeffs = [int(k) for k in rng.integers(1, 7, size=len(lower))]
        self.coeffs = rng.uniform(-1.0, 1.0, size=int(np.prod(self.n_coeffs)))

    def eval(self, x: np.ndarray) -> float:
        """The series at x, from the coefficients with the first variable fastest."""
        s = (2.0 * x - self.lower - self.upper) / (self.upper - self.lower)
        # Axis 0 is the last variable, the last axis the first one.
        c = self.coeffs.reshape(self.n_coeffs[::-1])
        for i, K in enumerate(self.n_coeffs):
            c = c @ chebyshev.chebvander(s[i], K - 1).ravel()
        return float(c)


def _make_table(
    n_vars: int, seed: int
) -> tuple[Ncm.SpectralTensor, list[Panel], np.ndarray, np.ndarray]:
    """A table on the product of [a_i, b_i], divided as a grid at m_i in every variable.

    The first leaf of the grid is split again at its midpoint in the last variable.
    Every leaf is a panel with random coefficients.
    """
    rng = np.random.default_rng(seed)
    a = -1.0 + 0.1 * np.arange(n_vars)
    m = 0.25 + 0.05 * np.arange(n_vars)
    b = 1.5 + 0.2 * np.arange(n_vars)
    st = Ncm.SpectralTensor.new(n_vars)
    panels: list[Panel] = []

    st.set_domain(a, b)
    for i in range(n_vars):
        assert st.split_empty(i, [a[i], m[i], b[i]]) == 2**i

    first = next(
        node
        for node in range(st.get_n_nodes())
        if st.get_node_kind(node) == Ncm.SpectralTensorNodeKind.EMPTY
    )
    lower, upper = st.get_node_bounds(first)
    child = st.split(first, n_vars - 1, [0.5 * (lower[-1] + upper[-1])])
    assert child == st.get_n_nodes() - 2

    for node in range(st.get_n_nodes()):
        if st.get_node_kind(node) != Ncm.SpectralTensorNodeKind.EMPTY:
            continue
        lower, upper = st.get_node_bounds(node)
        panel = Panel(node, lower, upper, rng)
        assert st.set_panel(node, panel.n_coeffs, panel.coeffs) == len(panels)
        panels.append(panel)

    return st, panels, a, b


@pytest.fixture(name="n_vars", params=[1, 2, 3])
def fixture_n_vars(request) -> int:
    """Number of variables of the table."""
    return request.param


def test_counts_and_range(n_vars: int) -> None:
    """Numbers of variables, nodes, empty nodes and panels, and the domain."""
    st, panels, a, b = _make_table(n_vars, 11)

    assert st.get_n_vars() == n_vars
    assert st.props.n_vars == n_vars
    assert st.get_n_panels() == len(panels) == 2**n_vars + 1
    assert st.get_n_nodes() == 2 ** (n_vars + 1) + 1
    assert st.get_n_empty() == 0
    for i in range(n_vars):
        assert st.get_range(i) == (a[i], b[i])


def test_tree(n_vars: int) -> None:
    """The root split node, the kinds of the nodes and the panel of each panel node."""
    st, panels, _, _ = _make_table(n_vars, 12)

    edges, var, first_child = st.get_node_split(0)
    assert (list(edges), var, first_child) == ([0.25], 0, 1)

    for p, panel in enumerate(panels):
        assert st.get_node_kind(panel.node) == Ncm.SpectralTensorNodeKind.PANEL
        assert st.get_node_panel(panel.node) == p

    kinds = [st.get_node_kind(node) for node in range(st.get_n_nodes())]
    assert kinds.count(Ncm.SpectralTensorNodeKind.PANEL) == len(panels)
    assert kinds.count(Ncm.SpectralTensorNodeKind.EMPTY) == 0


def test_get_panel(n_vars: int) -> None:
    """The arrays of every panel come back as set."""
    st, panels, _, _ = _make_table(n_vars, 13)

    for p, panel in enumerate(panels):
        coeffs, lower, upper, n_coeffs = st.get_panel(p)
        assert list(lower) == panel.lower.tolist()
        assert list(upper) == panel.upper.tolist()
        assert list(n_coeffs) == panel.n_coeffs
        assert list(coeffs) == panel.coeffs.tolist()


def test_eval(n_vars: int) -> None:
    """eval, eval_on and find_panel against the series evaluated by numpy."""
    st, panels, a, b = _make_table(n_vars, 14)
    rng = np.random.default_rng(15)

    for x in rng.uniform(a, b, size=(N_POINTS, n_vars)):
        p = st.find_panel(x)
        assert np.all(x >= panels[p].lower) and np.all(x < panels[p].upper)

        ref = panels[p].eval(x)
        l1 = np.sum(np.abs(panels[p].coeffs))
        val = st.eval(x)
        assert abs(val - ref) <= 1.0e-14 * l1
        assert st.eval_on(p, x) == val
        assert st.eval(x.tolist()) == val


def test_eval_array(n_vars: int) -> None:
    """eval_array on points in random order and grouped by panel, against the numpy series."""
    st, panels, a, b = _make_table(n_vars, 18)
    rng = np.random.default_rng(19)
    points = rng.uniform(a, b, size=(4 * N_POINTS, n_vars))
    found = np.array([st.find_panel(x) for x in points])

    for pts, ps in ((points, found), (points[np.argsort(found)], np.sort(found))):
        values = st.eval_array(pts.ravel())
        assert len(values) == len(pts)

        for x, p, val in zip(pts, ps, values):
            l1 = np.sum(np.abs(panels[p].coeffs))
            assert abs(val - panels[p].eval(x)) <= 1.0e-14 * l1

    assert len(st.eval_array([])) == 0


def test_serialize(n_vars: int) -> None:
    """A serialized copy has the same nodes and panels and gives the same values."""
    st, panels, a, b = _make_table(n_vars, 16)
    ser = Ncm.Serialize.new(Ncm.SerializeOpt.NONE)
    dup = ser.dup_obj(st)
    rng = np.random.default_rng(17)

    assert isinstance(dup, Ncm.SpectralTensor)
    assert dup.get_n_vars() == n_vars
    assert dup.get_n_nodes() == st.get_n_nodes()
    assert dup.get_n_panels() == len(panels)
    for node in range(st.get_n_nodes()):
        assert dup.get_node_kind(node) == st.get_node_kind(node)
        assert dup.get_node_bounds(node) == st.get_node_bounds(node)
    for p in range(len(panels)):
        assert dup.get_panel(p) == st.get_panel(p)

    for x in rng.uniform(a, b, size=(N_POINTS, n_vars)):
        assert dup.eval(x) == st.eval(x)


def test_reset() -> None:
    """reset leaves the domain with an empty root; a new domain replaces the old one."""
    st, _, a, b = _make_table(2, 18)

    st.reset()
    assert (st.get_n_nodes(), st.get_n_empty(), st.get_n_panels()) == (1, 1, 0)
    assert st.get_node_kind(0) == Ncm.SpectralTensorNodeKind.EMPTY
    assert st.get_range(0) == (a[0], b[0])
    assert st.get_range(1) == (a[1], b[1])

    st.set_domain([-2.0, 0.0], [3.0, 1.0])
    assert st.get_range(0) == (-2.0, 3.0)
    assert st.get_node_bounds(0) == ([-2.0, 0.0], [3.0, 1.0])


def test_node_kind_enum() -> None:
    """The node kinds are introspected with their nicks."""
    assert Ncm.SpectralTensorNodeKind.EMPTY.value_nick == "empty"
    assert Ncm.SpectralTensorNodeKind.SPLIT.value_nick == "split"
    assert Ncm.SpectralTensorNodeKind.PANEL.value_nick == "panel"

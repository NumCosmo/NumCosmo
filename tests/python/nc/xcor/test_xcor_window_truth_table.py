#!/usr/bin/env python
#
# test_xcor_window_truth_table.py
#
# Wed Aug 27 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_xcor_window_truth_table.py
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

"""Check that the certified radial-integral tables are references.

The radial integral

.. math:: I_\\ell^{(d)}(k) = \\int W(\\chi)\\, j_\\ell^{(d)}(k\\chi)\\, \\mathrm{d}\\chi

is where three quarters of an ``NcXcorSolver`` run is spent.
``data/truth_tables/xcor/xcor_window_ilk.json.gz`` and its ``_d1``, ``_d2``
companions hold it computed in Arb for every analytic window, on a grid that
follows the peak across multipoles. The library is compared with them in C: the
Levin integrator in ``tests/c/nc/xcor/test_nc_xcor_window_levin.c``, the closures
in ``tests/c/nc/xcor/test_nc_xcor_window_dense.c``, both reading the ``.bin``
``tests/tools/make_xcor_window_bin.py`` writes from these tables. What stays here
reads what only the JSON carries, the ball radii: each table must be certified far
below the tolerance the C test holds the integrator to.

Regenerating the tables needs FLINT; running this test does not. See
``tests/tools/make_xcor_window_truth_table.py``.
"""

import gzip
import json
import pathlib

import pytest

from numcosmo_py import Nc, Ncm

Ncm.cfg_init()

pytestmark = [pytest.mark.xcor, pytest.mark.xdist_group("window_truth")]

TRUTH_TABLE = "truth_tables/xcor/xcor_window_ilk.json.gz"

DERIV_TRUTH_TABLE = {
    1: "truth_tables/xcor/xcor_window_ilk_d1.json.gz",
    2: "truth_tables/xcor/xcor_window_ilk_d2.json.gz",
}

# The relative tolerances test_nc_xcor_window_levin.c holds the integrator to, at
# d = 0 and at d = 1, 2; the certified radii have to sit far below them.
RTOL = 1.0e-8
DERIV_RTOL = 1.0e-8


def _make_kernel(shape, dist, ps, sbi, ctor):
    """Build one analytic kernel from the constructor arguments in the table."""
    match shape:
        case "gauss":
            return Nc.XcorKernelAnalyticGauss.new_full(dist, ps, *ctor, sbi)
        case "tophat":
            return Nc.XcorKernelAnalyticTophat.new_full(dist, ps, *ctor, sbi)
        case "tophat_smooth":
            return Nc.XcorKernelAnalyticTophatSmooth.new_full(dist, ps, *ctor, sbi)
        case "student_t":
            return Nc.XcorKernelAnalyticStudentT.new_full(dist, ps, *ctor, sbi)
        case "power_exp":
            return Nc.XcorKernelAnalyticPowerExp.new_full(dist, ps, *ctor, sbi)
        case "lensing":
            return Nc.XcorKernelAnalyticLensing.new_full(dist, ps, *ctor, sbi)
        case "multi":
            mu, sigma, weight, n_sigma = ctor

            return Nc.XcorKernelAnalyticMulti.new_full(
                dist,
                ps,
                Ncm.Vector.new_array(mu),
                Ncm.Vector.new_array(sigma),
                Ncm.Vector.new_array(weight),
                n_sigma,
                sbi,
            )
        case _:
            raise ValueError(f"truth table names an unknown shape {shape!r}")


@pytest.fixture(name="truth_table", scope="module")
def fixture_truth_table() -> dict:
    """Load the certified table once for the module."""
    path = pathlib.Path(Ncm.cfg_get_data_filename(TRUTH_TABLE, True))

    with gzip.open(path, "rt") as f:
        return json.load(f)


def test_table_is_certified_far_below_the_tolerance(truth_table: dict) -> None:
    """The reference must be a reference: its own uncertainty cannot matter.

    Guards a regeneration that silently lowered the precision target -- the
    comparison would still pass while no longer checking anything. The target is per
    multipole, because near the turning point the precision needed grows steeply with
    ell: the high multipoles are certified to 1e-18 rather than 1e-25. Both checks
    below are what makes that safe, the second one being the one that matters.
    """
    targets = truth_table["target_rel"]
    targets = (
        targets if isinstance(targets, list) else [targets] * len(truth_table["ells"])
    )

    # Every target must sit far below the comparison this file makes, or the reference
    # stops being one.
    assert max(targets) < 1.0e-3 * RTOL

    for index, (ell, target) in enumerate(zip(truth_table["ells"], targets)):
        worst = max(
            (
                radius / abs(float(value))
                for shape in truth_table["shapes"].values()
                for value, radius in zip(shape["table"][index], shape["radius"][index])
                if float(value) != 0.0
            ),
            default=0.0,
        )

        assert worst < target, f"ell = {ell} certified only to {worst:.2e}"


@pytest.fixture(name="deriv_truth_tables", scope="module")
def fixture_deriv_truth_tables() -> dict:
    """Load the certified Bessel-derivative tables once for the module."""
    tables = {}

    for deriv, name in DERIV_TRUTH_TABLE.items():
        path = pathlib.Path(Ncm.cfg_get_data_filename(name, True))

        with gzip.open(path, "rt") as f:
            tables[deriv] = json.load(f)

    return tables


@pytest.mark.parametrize("deriv", sorted(DERIV_TRUTH_TABLE))
def test_deriv_table_is_certified_far_below_the_tolerance(
    deriv: int, deriv_truth_tables: dict
) -> None:
    """Each derivative table must be a reference, on the same terms as the d = 0 one."""
    table = deriv_truth_tables[deriv]

    assert table["bessel_deriv"] == deriv

    targets = table["target_rel"]
    targets = targets if isinstance(targets, list) else [targets] * len(table["ells"])

    assert max(targets) < 1.0e-3 * DERIV_RTOL

    for index, (ell, target) in enumerate(zip(table["ells"], targets)):
        worst = max(
            (
                radius / abs(float(value))
                for shape in table["shapes"].values()
                for value, radius in zip(shape["table"][index], shape["radius"][index])
                if float(value) != 0.0
            ),
            default=0.0,
        )

        assert (
            worst < target
        ), f"deriv {deriv}, ell = {ell} certified only to {worst:.2e}"

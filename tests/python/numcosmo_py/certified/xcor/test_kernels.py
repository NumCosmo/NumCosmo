#
# test_kernels.py
#
# Thu Sep 24 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_kernels.py
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

"""The library kernel of every vocabulary window.

Whether the kernels' integrals are right is tested against certified values in
``tests/python/nc/xcor/test_xcor_window_truth_table.py``; this file tests that the
builder constructs the window the vocabulary names.
"""

import gzip
import json
import pathlib

import numpy as np
import pytest

from numcosmo_py import Nc, Ncm
from numcosmo_py.certified.xcor import kernels
from numcosmo_py.certified.xcor import windows as W

Ncm.cfg_init()

pytestmark = [pytest.mark.xcor]

TRUTH_TABLE = "truth_tables/xcor/xcor_window_ilk.json.gz"


@pytest.fixture(name="cosmo_bits", scope="module")
def fixture_cosmo_bits() -> tuple:
    """Build the filler objects every analytic kernel is constructed with."""
    return kernels.make_cosmo_bits()


def _prepared(window, cosmo_bits: tuple) -> Nc.XcorKernelRadial:
    cosmo, dist, ps = cosmo_bits
    kernel = kernels.build_kernel(window, dist, ps)
    kernel.set_l_limber(-1)
    kernel.prepare(cosmo)

    return kernel


def test_cosmo_bits_distance_is_prepared(cosmo_bits: tuple) -> None:
    """The distance is usable as returned."""
    cosmo, dist, _ = cosmo_bits

    assert dist.comoving(cosmo, 1.0) > 0.0


@pytest.mark.parametrize("name", sorted(W.WINDOWS))
def test_every_window_builds(name: str, cosmo_bits: tuple) -> None:
    """Every window builds, prepares, and carries its weight and k-dependence."""
    window = W.get(name)
    kernel = _prepared(name, cosmo_bits)

    assert isinstance(kernel, kernels.FORM_CLASSES[window.form])
    assert kernel.get_n_comps() >= 1

    for comp in range(kernel.get_n_comps()):
        assert kernel.get_comp_bessel_deriv(comp) == window.deriv

    assert (kernel.peek_kdep() is not None) == bool(window.kdep)


def test_table_supports_match_the_vocabulary(cosmo_bits: tuple) -> None:
    """Each window of the certified table has the support the table was certified on."""
    path = pathlib.Path(Ncm.cfg_get_data_filename(TRUTH_TABLE, True))

    with gzip.open(path, "rt") as handle:
        table = json.load(handle)

    for case, entry in table["shapes"].items():
        kernel = _prepared(W.resolve(case), cosmo_bits)
        supports = [kernel.get_comp_support(c) for c in range(kernel.get_n_comps())]
        low = min(s[0] for s in supports)
        high = max(s[1] for s in supports)

        assert np.allclose([low, high], entry["support"], rtol=1.0e-15, atol=0.0), case


@pytest.mark.parametrize("alias", sorted(W.ALIASES))
def test_historical_names_build_the_current_window(
    alias: str, cosmo_bits: tuple
) -> None:
    """A retired name builds exactly the window its current name does."""
    old = _prepared(alias, cosmo_bits)
    new = _prepared(W.ALIASES[alias], cosmo_bits)

    assert old.get_n_comps() == new.get_n_comps()

    for comp in range(old.get_n_comps()):
        low, high = new.get_comp_support(comp)

        assert old.get_comp_support(comp) == (low, high)
        chis = np.linspace(low, high, 17)
        assert [old.eval_W_comp(comp, c) for c in chis] == [
            new.eval_W_comp(comp, c) for c in chis
        ]

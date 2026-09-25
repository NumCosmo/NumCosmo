#
# test_measure.py
#
# Thu Sep 24 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_measure.py
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

"""The radial campaign measures every cell against the right scale.

What it finds is recorded as campaign results, not gated here.
"""

import pathlib

import numpy as np
import pytest

from numcosmo_py import GLib, Ncm
from numcosmo_py.certified import tables
from numcosmo_py.certified.xcor import measure

Ncm.cfg_init()

pytestmark = [pytest.mark.xcor]


def _data(name: str) -> pathlib.Path:
    try:
        return pathlib.Path(Ncm.cfg_get_data_filename(name, True))
    except (RuntimeError, GLib.Error):
        pytest.skip(f"{name} is not installed")


def test_radial_campaign_scales_by_the_peak_table() -> None:
    """Every cell is measured, and scaled by its window's peak, not the tail's."""
    tail = _data("truth_tables/xcor/xcor_window_ilk_xnu_ext.bin")
    main = tables.convert(_data("truth_tables/xcor/xcor_window_ilk.json.gz"))
    peaks = measure.peaks(main)

    campaign = measure.radial(tail, "test", main, "xcor_window_ilk.json.gz")

    assert campaign.meta["config_scale_table"] == "xcor_window_ilk.json.gz"
    assert set(campaign.groups) == set(tables.groups(tables.read(tail)))

    for name, rows in campaign.groups.items():
        ells = np.asarray(rows["ell"])
        expected = [peaks[(name, int(e))] for e in ells]

        np.testing.assert_array_equal(np.asarray(rows["scale"]), expected)
        assert np.all(np.isfinite(np.asarray(rows["measured"])))
        assert np.all(np.abs(np.asarray(rows["reference"])) < np.asarray(expected))

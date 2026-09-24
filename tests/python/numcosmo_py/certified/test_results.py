#
# test_results.py
#
# Thu Sep 24 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_results.py
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

"""A campaign file records what was compared with what, and the deviations."""

import hashlib
import pathlib

import numpy as np
import pytest

from numcosmo_py import GLib, Ncm
from numcosmo_py.certified import results, tables

Ncm.cfg_init()

pytestmark = [pytest.mark.xcor]


@pytest.fixture(name="reference")
def fixture_reference(tmp_path: pathlib.Path) -> pathlib.Path:
    """Build a certified table in the binfile layout from the committed C_ell table."""
    try:
        source = pathlib.Path(
            Ncm.cfg_get_data_filename("truth_tables/xcor/xcor_kquad.json.gz", True)
        )
    except (RuntimeError, GLib.Error):
        pytest.skip("xcor_kquad.json.gz is not installed")

    path = tmp_path / "xcor_kquad.bin"
    tables.write(tables.convert(source), path)

    return path


def test_campaign_round_trip(reference: pathlib.Path, tmp_path: pathlib.Path) -> None:
    """Metadata, provenance and deviations survive the binfile."""
    _, cells = tables.groups(tables.read(reference))["A1"]
    measured = np.asarray(cells["value"]) * (1.0 + 1.0e-9)
    scale = np.full(len(cells), 2.0 * float(np.max(np.abs(cells["value"]))))

    campaign = results.Campaign(
        reference, "numcosmo", "test", {"closure": 1, "reltol": 1.0e-6}
    )
    campaign.add(
        "A1",
        {
            "ell": np.asarray(cells["ell"]),
            "reference": np.asarray(cells["value"]),
            "radius": np.asarray(cells["radius"]),
            "measured": measured,
            "scale": scale,
        },
    )
    path = tmp_path / "campaign.bin"
    campaign.write(path)

    back = tables.read(path)
    meta = tables.metadata(back)
    _, rows = tables.groups(back)["A1"]

    assert meta["kind"] == "results"
    assert (
        meta["reference_sha256"] == hashlib.sha256(reference.read_bytes()).hexdigest()
    )
    assert meta["reference_kind"] == "kquad"
    assert meta["config_reltol"] == 1.0e-6
    assert meta["env_host"] and meta["env_cpu"] and meta["env_libraries"]
    assert "numcosmo" in meta["env_libraries"]
    np.testing.assert_allclose(rows["dev_rel"], 1.0e-9, rtol=1.0e-6)
    np.testing.assert_allclose(
        rows["dev_scaled"],
        np.abs(measured - np.asarray(cells["value"])) / scale,
        rtol=1.0e-15,
    )


def test_missing_required_column_is_refused(reference: pathlib.Path) -> None:
    """A row without its scale cannot be stored: a deviation needs one."""
    campaign = results.Campaign(reference, "numcosmo", "test")

    with pytest.raises(ValueError, match="scale"):
        campaign.add("A1", {"reference": [1.0], "radius": [0.0], "measured": [1.0]})

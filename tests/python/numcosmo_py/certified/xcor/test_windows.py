#
# test_windows.py
#
# Thu Sep 24 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_windows.py
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

"""The window vocabulary is internally consistent.

Agreement with the generators' and the pair matrix's rosters is tested in
``tests/python/nc/xcor/test_window_vocabulary.py``; this file tests the vocabulary alone.
"""

import pytest

from numcosmo_py.certified.xcor import windows as W

pytestmark = [pytest.mark.xcor]


def test_every_vocabulary_name_is_reachable() -> None:
    """Aliases point at real windows, and no window is orphaned by a typo."""
    for alias, canonical in W.ALIASES.items():
        assert canonical in W.WINDOWS, f"alias {alias} points at nothing"
        assert alias not in W.WINDOWS, f"alias {alias} shadows a real window"

    for name, window in W.WINDOWS.items():
        assert window.form in W.FORM_PARAMS, f"{name} has an unknown form"
        assert len(window.params) == len(W.FORM_PARAMS[window.form]), (
            f"{name} carries {len(window.params)} parameters, "
            f"{window.form} takes {len(W.FORM_PARAMS[window.form])}"
        )
        assert 0 <= window.deriv <= 2, f"{name} has derivative order {window.deriv}"

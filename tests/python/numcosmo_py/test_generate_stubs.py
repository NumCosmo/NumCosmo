#!/usr/bin/env python
#
# test_generate_stubs.py
#
# Wed Oct 21 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_generate_stubs.py
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

"""Unit tests for numcosmo_py.generate_stubs.

The committed stubs are the generator output formatted by black (see
numcosmo_py/update_pyi.sh) and contain no comments, so comparing syntax trees
compares everything the stubs hold without depending on the black version.
"""

import ast
import difflib
import sys
from pathlib import Path

import pytest

import numcosmo_py
from numcosmo_py import generate_stubs

STUBS_DIR = Path(numcosmo_py.__file__).parent


def _normalize(source: str) -> str:
    return ast.unparse(ast.parse(source))


@pytest.mark.parametrize(
    "module, stub", [("NumCosmoMath", "ncm.pyi"), ("NumCosmo", "nc.pyi")]
)
def test_generate_stubs_matches_committed(
    module: str, stub: str, monkeypatch: pytest.MonkeyPatch, capsys
) -> None:
    """Generated stubs match the committed ones."""
    monkeypatch.setattr(sys, "argv", ["generate_stubs.py", module, "1.0"])
    generate_stubs.main_func()

    generated = _normalize(capsys.readouterr().out)
    committed = _normalize((STUBS_DIR / stub).read_text(encoding="utf-8"))

    if generated != committed:
        diff = difflib.unified_diff(
            committed.splitlines(),
            generated.splitlines(),
            f"committed {stub}",
            f"generated {stub}",
            lineterm="",
        )
        pytest.fail(
            f"{stub} differs from the generator output; run "
            "numcosmo_py/update_pyi.sh and commit the stubs:\n"
            + "\n".join(list(diff)[:60])
        )

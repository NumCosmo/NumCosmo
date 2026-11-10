#!/usr/bin/env python
#
# test_import_no_init.py
#
# Mon Nov 2 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_import_no_init.py
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

"""Test that importing a numcosmo_py module does not initialize NumCosmo.

Under an MPI launcher, Ncm.cfg_init() puts every rank other than 0 into the
worker loop. A call run at import time therefore stops those ranks before the
importing program registers its Python types or builds its models.
"""

import ast
from pathlib import Path

import numcosmo_py

PACKAGE_DIR = Path(numcosmo_py.__file__).parent


class _ImportTimeCfgInit(ast.NodeVisitor):
    """Collect the Ncm.cfg_init* calls that run when a module is imported.

    Function bodies run only when called and the body of
    ``if __name__ == "__main__"`` only when the file runs as a script, so both
    are skipped; module-level and class-body statements are visited.
    """

    def __init__(self) -> None:
        self.lines: list[int] = []

    def visit_FunctionDef(self, node: ast.FunctionDef) -> None:
        """Skip the function body."""

    def visit_AsyncFunctionDef(self, node: ast.AsyncFunctionDef) -> None:
        """Skip the function body."""

    def visit_Lambda(self, node: ast.Lambda) -> None:
        """Skip the lambda body."""

    def visit_If(self, node: ast.If) -> None:
        """Visit the statements of an if, skipping the body of the main guard."""
        if _is_main_guard(node.test):
            for stmt in node.orelse:
                self.visit(stmt)
        else:
            self.generic_visit(node)

    def visit_Call(self, node: ast.Call) -> None:
        """Record a call to an attribute named cfg_init*."""
        if isinstance(node.func, ast.Attribute) and node.func.attr.startswith(
            "cfg_init"
        ):
            self.lines.append(node.lineno)
        self.generic_visit(node)


def _is_main_guard(test: ast.expr) -> bool:
    """Return whether test is ``__name__ == "__main__"``."""
    return (
        isinstance(test, ast.Compare)
        and isinstance(test.left, ast.Name)
        and test.left.id == "__name__"
        and len(test.ops) == 1
        and isinstance(test.ops[0], ast.Eq)
        and len(test.comparators) == 1
        and isinstance(test.comparators[0], ast.Constant)
        and test.comparators[0].value == "__main__"
    )


def test_no_cfg_init_at_import() -> None:
    """No module of numcosmo_py calls Ncm.cfg_init* when it is imported."""
    offenders: list[str] = []
    for path in sorted(PACKAGE_DIR.rglob("*.py")):
        visitor = _ImportTimeCfgInit()
        visitor.visit(ast.parse(path.read_text(encoding="utf-8"), filename=str(path)))
        offenders.extend(
            f"{path.relative_to(PACKAGE_DIR)}:{line}" for line in visitor.lines
        )

    assert not offenders, f"cfg_init called at import time: {offenders}"

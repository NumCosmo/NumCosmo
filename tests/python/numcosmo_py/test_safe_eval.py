#!/usr/bin/env python
#
# test_safe_eval.py
#
# Sun Jul 26 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_safe_eval.py
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

"""Unit tests for numcosmo_py.safe_eval module."""

import math

import numpy as np
import pytest

from numcosmo_py.safe_eval import compile_expr, SafeExprError


def test_identity() -> None:
    """A bare variable name evaluates to itself."""
    evaluate = compile_expr("x", ["x"])
    assert evaluate({"x": 3.5}) == 3.5


def test_arithmetic_operators() -> None:
    """All whitelisted binary/unary operators compute the expected value."""
    evaluate = compile_expr("-x + 2 * y - x / y + x ** 2 + x // y + x % y", ["x", "y"])
    x, y = 7.0, 2.0
    expected = -x + 2 * y - x / y + x**2 + x // y + x % y
    assert evaluate({"x": x, "y": y}) == pytest.approx(expected)


def test_pow10_transform() -> None:
    """The motivating use case: 10**log10(M)."""
    evaluate = compile_expr("10**x", ["x"])
    assert evaluate({"x": 2.0}) == pytest.approx(100.0)


def test_multi_variable_expression() -> None:
    """Expression combining two bound variables, e.g. (y / 100)**2 * x."""
    evaluate = compile_expr("(y / 100)**2 * x", ["x", "y"])
    assert evaluate({"x": 0.3, "y": 70.0}) == pytest.approx((70.0 / 100.0) ** 2 * 0.3)


def test_whitelisted_functions() -> None:
    """Whitelisted math functions evaluate to their NumPy equivalents."""
    evaluate = compile_expr("log10(x) + sqrt(y) + exp(1)", ["x", "y"])
    x, y = 1000.0, 4.0
    expected = np.log10(x) + np.sqrt(y) + np.exp(1)
    assert evaluate({"x": x, "y": y}) == pytest.approx(expected)


def test_constants() -> None:
    """The pi/e constants are available without being declared as variables."""
    evaluate = compile_expr("pi + e", [])
    assert evaluate({}) == pytest.approx(math.pi + math.e)


@pytest.mark.parametrize(
    "expr",
    [
        "__import__('os').system('echo hi')",
        "().__class__",
        "x.__class__",
        "open('/etc/passwd')",
        "exec('1')",
        "[i for i in range(3)]",
        "lambda: 1",
        "1 if x > 0 else -1",
        "x[0]",
        "'a string'",
        "not_a_variable",
        "sinh_missing_paren",
        "cos(x, y)",  # extra positional arg is fine syntactically, rejected at runtime
        "unknown_func(x)",
        "os.system('echo hi')",
    ],
)
def test_rejects_unsafe_or_unknown_syntax(expr: str) -> None:
    """Anything outside the whitelist is rejected before it can execute."""
    with pytest.raises(SafeExprError):
        compile_expr(expr, ["x"])


def test_rejects_keyword_arguments() -> None:
    """Calls with keyword arguments are rejected even for whitelisted functions."""
    with pytest.raises(SafeExprError):
        compile_expr("log(x, base=10)", ["x"])


def test_rejects_syntax_error() -> None:
    """A malformed expression raises SafeExprError, not a raw SyntaxError."""
    with pytest.raises(SafeExprError):
        compile_expr("x +", ["x"])


def test_variable_shadowing_builtins_not_allowed_to_escape() -> None:
    """No builtins are reachable even though `__builtins__` name itself is blocked."""
    with pytest.raises(SafeExprError):
        compile_expr("__builtins__", [])

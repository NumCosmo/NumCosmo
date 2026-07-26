#
# safe_eval.py
#
# Sun Jul 26 00:00:00 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# safe_eval.py
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

"""Safe evaluation of simple mathematical expressions.

Restricted expression evaluator used by CLI options that let users combine
named values with basic arithmetic (e.g. "(y / 100)**2 * x"), without
exposing a general-purpose eval() to CLI input.
"""

from __future__ import annotations

import ast
import math
import operator
from typing import Callable, Collection, Mapping

import numpy as np


class SafeExprError(ValueError):
    """Raised when an expression contains disallowed or invalid syntax."""


_BINOPS: dict[type, Callable[[float, float], float]] = {
    ast.Add: operator.add,
    ast.Sub: operator.sub,
    ast.Mult: operator.mul,
    ast.Div: operator.truediv,
    ast.FloorDiv: operator.floordiv,
    ast.Mod: operator.mod,
    ast.Pow: operator.pow,
}

_UNARYOPS: dict[type, Callable[[float], float]] = {
    ast.UAdd: operator.pos,
    ast.USub: operator.neg,
}

_FUNCTIONS: dict[str, Callable[..., float]] = {
    "log": np.log,
    "log10": np.log10,
    "log2": np.log2,
    "exp": np.exp,
    "sqrt": np.sqrt,
    "abs": np.abs,
    "sin": np.sin,
    "cos": np.cos,
    "tan": np.tan,
    "arcsin": np.arcsin,
    "arccos": np.arccos,
    "arctan": np.arctan,
    "sinh": np.sinh,
    "cosh": np.cosh,
    "tanh": np.tanh,
}

_CONSTANTS: dict[str, float] = {
    "pi": math.pi,
    "e": math.e,
}


def _validate(node: ast.AST, variables: Collection[str]) -> None:
    """Recursively verify that `node` only uses whitelisted syntax."""
    if isinstance(node, ast.Expression):
        _validate(node.body, variables)
    elif isinstance(node, ast.BinOp):
        if type(node.op) not in _BINOPS:
            raise SafeExprError(f"Operator `{type(node.op).__name__}` is not allowed.")
        _validate(node.left, variables)
        _validate(node.right, variables)
    elif isinstance(node, ast.UnaryOp):
        if type(node.op) not in _UNARYOPS:
            raise SafeExprError(
                f"Unary operator `{type(node.op).__name__}` is not allowed."
            )
        _validate(node.operand, variables)
    elif isinstance(node, ast.Call):
        if (
            not isinstance(node.func, ast.Name)
            or node.func.id not in _FUNCTIONS
            or node.keywords
        ):
            raise SafeExprError(
                "Only calls to a fixed set of math functions are allowed: "
                f"{', '.join(sorted(_FUNCTIONS))}."
            )
        for arg in node.args:
            _validate(arg, variables)
    elif isinstance(node, ast.Name):
        if node.id not in variables and node.id not in _CONSTANTS:
            raise SafeExprError(f"Unknown name `{node.id}` in expression.")
    elif isinstance(node, ast.Constant):
        if not isinstance(node.value, (int, float)):
            raise SafeExprError("Only numeric constants are allowed.")
    else:
        raise SafeExprError(
            f"Expression syntax `{type(node).__name__}` is not allowed."
        )


def compile_expr(
    expr: str, variables: Collection[str]
) -> Callable[[Mapping[str, float]], float]:
    """Compile a restricted mathematical expression into a callable.

    The expression may use ``+ - * / // % **`` (and unary minus) with the
    given `variables` names, numeric literals, the constants ``pi``/``e``,
    and calls to a small whitelist of math functions (log, log10, log2,
    exp, sqrt, abs, sin, cos, tan, arcsin, arccos, arctan, sinh, cosh,
    tanh). Attribute access, subscripting, comprehensions, and calls to
    anything outside that whitelist are rejected by construction, so this
    is safe to use directly on CLI-provided strings.
    """
    try:
        tree = ast.parse(expr, mode="eval")
    except SyntaxError as exc:
        raise SafeExprError(f"Invalid expression `{expr}`: {exc}") from exc

    _validate(tree, variables)
    code = compile(tree, "<compute-expr>", "eval")

    def _evaluate(values: Mapping[str, float]) -> float:
        env = {**_FUNCTIONS, **_CONSTANTS, **values}
        return eval(code, {"__builtins__": {}}, env)  # AST-validated above: safe.

    return _evaluate

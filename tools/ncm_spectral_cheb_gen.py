#!/usr/bin/env python3
#
# ncm_spectral_cheb_gen.py
#
# Thu Nov 12 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# ncm_spectral_cheb_gen.py
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

"""Generate the Chebyshev transform kernels of NcmSpectral.

Writes the C header ``numcosmo/ncm/algebra/ncm_spectral_cheb_private.h``. For each
level ``l`` from 0 to :data:`MAX_LEVEL` the header has one function that maps the
values of a function at the ``2^l + 1`` Chebyshev-Lobatto nodes to its Chebyshev
coefficients, for a batch of columns. The body of each function is straight-line
code: one C statement per arithmetic operation.

Usage::

    python tools/ncm_spectral_cheb_gen.py            # write the header
    python tools/ncm_spectral_cheb_gen.py --check    # compare with the header

The build does not run this script. The header is committed, and ``--check``
exits with status 1 when the committed header differs from the generated text.

Method. Let ``n = 2^l``, let ``x_j`` be the value at the node ``cos(pi j / n)``
and let ``y_k`` be the coefficients, ``f = sum_k y_k T_k``. The nodes with even
``j`` are the nodes of level ``l - 1``, with coefficients ``c_k``. The odd nodes
enter through their type-II cosine transform::

    d_k = (2 / n) sum_i x_{2 i + 1} cos(pi (2 i + 1) k / n),   0 < k < n / 2,

and ``d_0`` by the same sum with the factor ``1 / n``. Then::

    y_k = c_k / 2 + d_k,   y_{n - k} = c_k / 2 - d_k,   0 <= k < n / 2,
    y_{n / 2} = c_{n / 2}.

The type-II transform of ``N`` values ``v_i`` splits in the same way. Its even
outputs are the type-II transform of the ``N / 2`` sums ``v_i + v_{N - 1 - i}``.
Its odd outputs ``O_k`` follow from the type-II transform ``Z_k`` of the scaled
differences ``2 cos(pi (2 i + 1) / (2 N)) (v_i - v_{N - 1 - i})`` by
``O_0 = Z_0 / 2`` and ``O_k = Z_k - O_{k - 1}``.

Every normalization factor is a power of two. The generator multiplies it into
the cosine factors, so the factors of the generated code are exact products.
"""

from __future__ import annotations

import argparse
import dataclasses
import math
import sys
from collections.abc import Sequence
from pathlib import Path

# Largest level with a generated kernel: 2^7 + 1 = 129 nodes.
MAX_LEVEL = 7

# File name of the generated header.
HEADER_NAME = "ncm_spectral_cheb_private.h"

# Path of the generated header in the source tree that holds this script.
DEFAULT_OUTPUT = (
    Path(__file__).resolve().parent.parent
    / "numcosmo"
    / "ncm"
    / "algebra"
    / HEADER_NAME
)

# Signature shared by every kernel, after the function name.
KERNEL_ARGUMENTS = "(const gdouble * restrict x, gdouble * restrict y, const gsize L)"


@dataclasses.dataclass(frozen=True)
class Input:
    """The value at one node: row ``row`` of the input array."""

    row: int


@dataclasses.dataclass(frozen=True)
class Temporary:
    """The result of statement ``index`` of a program."""

    index: int


Operand = Input | Temporary


@dataclasses.dataclass(frozen=True)
class Sum:
    """Statement that adds two operands."""

    left: Operand
    right: Operand


@dataclasses.dataclass(frozen=True)
class Difference:
    """Statement that subtracts ``right`` from ``left``."""

    left: Operand
    right: Operand


@dataclasses.dataclass(frozen=True)
class Product:
    """Statement that multiplies an operand by a constant factor."""

    operand: Operand
    factor: float


Statement = Sum | Difference | Product


@dataclasses.dataclass
class Program:
    """Straight-line program that computes the output rows from the input rows.

    The statements are in execution order. Statement ``i`` defines
    ``Temporary(i)`` and uses only inputs and earlier temporaries.
    """

    statements: list[Statement] = dataclasses.field(default_factory=list)
    outputs: list[Operand] = dataclasses.field(default_factory=list)

    def _append(self, statement: Statement) -> Temporary:
        """Append a statement and return the temporary it defines."""
        self.statements.append(statement)
        return Temporary(len(self.statements) - 1)

    def add(self, left: Operand, right: Operand) -> Temporary:
        """Append ``left + right``.

        :param left: first term
        :param right: second term
        :return: the sum
        """
        return self._append(Sum(left, right))

    def subtract(self, left: Operand, right: Operand) -> Temporary:
        """Append ``left - right``.

        :param left: minuend
        :param right: subtrahend
        :return: the difference
        """
        return self._append(Difference(left, right))

    def scale(self, operand: Operand, factor: float) -> Operand:
        """Append ``operand * factor``, or nothing when the factor is one.

        :param operand: value to multiply
        :param factor: constant factor
        :return: the product, or ``operand`` itself when ``factor`` is one
        """
        if factor == 1.0:
            return operand
        return self._append(Product(operand, factor))

    def count(self, kind: type[Statement]) -> int:
        """Count the statements of one kind.

        :param kind: statement class
        :return: number of statements of that class
        """
        return sum(isinstance(statement, kind) for statement in self.statements)


def odd_factors(n_points: int, scale: float) -> list[float]:
    """Compute the factors of the differences in a type-II cosine transform.

    :param n_points: number of points of the transform, a power of two above one
    :param scale: normalization of the odd outputs
    :return: ``scale * 2 cos(pi (2 i + 1) / (2 n_points))`` for ``i`` below
        ``n_points / 2``, or ``scale * cos(pi / 4)`` when ``n_points`` is two
    """
    if n_points == 2:
        return [scale * math.cos(math.pi / 4.0)]
    return [
        scale * 2.0 * math.cos(math.pi * (2 * i + 1) / (2.0 * n_points))
        for i in range(n_points // 2)
    ]


def odd_outputs(program: Program, scaled_differences: list[Operand]) -> list[Operand]:
    """Append the odd outputs of a type-II cosine transform.

    :param program: program that receives the statements
    :param scaled_differences: the differences ``v_i - v_{N - 1 - i}`` times
        their factors, ``N / 2`` values
    :return: the outputs ``O_k`` of odd index ``2 k + 1``, in the order of ``k``
    """
    if len(scaled_differences) == 1:
        return scaled_differences
    transformed = cosine_transform(program, scaled_differences, 1.0, 1.0)
    outputs = [program.scale(transformed[0], 0.5)]
    for value in transformed[1:]:
        outputs.append(program.subtract(value, outputs[-1]))
    return outputs


def cosine_transform(
    program: Program, values: list[Operand], scale: float, dc_scale: float
) -> list[Operand]:
    """Append the type-II cosine transform of ``values``.

    :param program: program that receives the statements
    :param values: the ``N`` input values ``v_i``, ``N`` a power of two
    :param scale: normalization of the outputs of index above zero
    :param dc_scale: normalization of the output of index zero
    :return: ``scale * sum_i v_i cos(pi (2 i + 1) k / (2 N))`` for ``k`` from
        zero to ``N - 1``, with ``dc_scale`` in place of ``scale`` at ``k = 0``
    """
    n_points = len(values)
    if n_points == 1:
        return [program.scale(values[0], dc_scale)]
    sums: list[Operand] = []
    scaled_differences: list[Operand] = []
    for i, factor in enumerate(odd_factors(n_points, scale)):
        mirror = values[n_points - 1 - i]
        sums.append(program.add(values[i], mirror))
        scaled_differences.append(
            program.scale(program.subtract(values[i], mirror), factor)
        )
    even = cosine_transform(program, sums, scale, dc_scale)
    odd = odd_outputs(program, scaled_differences)
    return [output for pair in zip(even, odd, strict=True) for output in pair]


def chebyshev_coefficients(
    program: Program, values: list[Operand], scale: float
) -> list[Operand]:
    """Append the Chebyshev coefficients of the node values ``values``.

    :param program: program that receives the statements
    :param values: the values at the ``n + 1`` nodes ``cos(pi j / n)``, ``n`` a
        power of two
    :param scale: factor applied to every coefficient
    :return: ``scale * y_k`` for ``k`` from zero to ``n``
    """
    n = len(values) - 1
    if n == 1:
        return [
            program.scale(program.add(values[0], values[1]), 0.5 * scale),
            program.scale(program.subtract(values[0], values[1]), 0.5 * scale),
        ]
    half_coarse = chebyshev_coefficients(program, values[0::2], 0.5 * scale)
    d = cosine_transform(program, values[1::2], scale * 2.0 / n, scale * 1.0 / n)
    lower: list[Operand] = []
    upper: list[Operand] = []
    for half_c_k, d_k in zip(half_coarse[: n // 2], d, strict=True):
        lower.append(program.add(half_c_k, d_k))
        upper.append(program.subtract(half_c_k, d_k))
    middle = program.add(half_coarse[n // 2], half_coarse[n // 2])
    return [*lower, middle, *reversed(upper)]


def kernel_program(level: int) -> Program:
    """Build the program of the kernel of one level.

    :param level: level ``l``, the kernel transforms ``2^l + 1`` node values
    :return: program whose output ``k`` is the coefficient ``y_k``
    """
    program = Program()
    inputs: list[Operand] = [Input(row) for row in range((1 << level) + 1)]
    program.outputs = chebyshev_coefficients(program, inputs, 1.0)
    return program


def kernel_name(level: int) -> str:
    """Name the C function of the kernel of one level.

    :param level: level of the kernel
    :return: C identifier of the function
    """
    return f"_ncm_spectral_cheb_kernel_{level}"


def render_operand(operand: Operand) -> str:
    """Render an operand as a C expression.

    :param operand: input row or temporary
    :return: element of the input array or name of the temporary
    """
    match operand:
        case Input(row):
            return f"x[{row} * L + c]"
        case Temporary(index):
            return f"t{index}"


def render_statement(index: int, statement: Statement) -> str:
    """Render a statement as the C definition of its temporary.

    :param index: position of the statement in its program
    :param statement: statement to render
    :return: one C declaration with initializer
    """
    match statement:
        case Sum(left, right):
            expression = f"{render_operand(left)} + {render_operand(right)}"
        case Difference(left, right):
            expression = f"{render_operand(left)} - {render_operand(right)}"
        case Product(operand, factor):
            expression = f"{render_operand(operand)} * {factor!r}"
    return f"const gdouble {render_operand(Temporary(index))} = {expression};"


def render_body(program: Program) -> list[str]:
    """Render the statements of a program, each output stored when defined.

    :param program: program to render
    :return: C statements of the loop body, without indentation
    :raises TypeError: if an output is an input row, which no kernel produces
    """
    stores: dict[int, list[int]] = {}
    for row, output in enumerate(program.outputs):
        if not isinstance(output, Temporary):
            raise TypeError(f"output {row} is not the result of a statement")
        stores.setdefault(output.index, []).append(row)
    lines: list[str] = []
    for index, statement in enumerate(program.statements):
        lines.append(render_statement(index, statement))
        lines.extend(
            f"y[{row} * L + c] = {render_operand(Temporary(index))};"
            for row in stores.get(index, [])
        )
    return lines


def render_kernel(level: int) -> str:
    """Render the C function of the kernel of one level.

    :param level: level of the kernel
    :return: C function definition preceded by its comment
    """
    program = kernel_program(level)
    n_mul = program.count(Product)
    n_add = program.count(Sum) + program.count(Difference)
    n_nodes = (1 << level) + 1
    body = "\n".join(f"    {line}" for line in render_body(program))
    return f"""\
/*
 * Level {level}, {n_nodes} nodes. Operations per column: {n_mul} multiplications,
 * {n_add} additions.
 */
static void
{kernel_name(level)} {KERNEL_ARGUMENTS}
{{
  gsize c;

  #pragma omp simd
  for (c = 0; c < L; c++)
  {{
{body}
  }}
}}
"""


PREAMBLE = f"""\
/***************************************************************************
 *            {HEADER_NAME}
 *
 *  Thu Nov 12 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2026 <vitenti@uel.br>
 * numcosmo is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * numcosmo is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/*
 * Generated by tools/ncm_spectral_cheb_gen.py. Do not edit: change the generator
 * and run it again.
 *
 * Kernels that compute Chebyshev coefficients from the values at the
 * Chebyshev-Lobatto nodes, for a batch of L columns. The kernel of level l reads
 * the 2^l + 1 rows of x, row j at x + j L holding the values at the node
 * cos (pi j / 2^l), and writes the 2^l + 1 rows of y, row k at y + k L holding
 * the coefficients y_k of f = sum_k y_k T_k. The arrays x and y must not overlap.
 * _ncm_spectral_cheb_kernels[l] is the kernel of level l, for l from 0 to
 * NCM_SPECTRAL_CHEB_MAX_LEVEL. The method is described in the generator.
 */

#ifndef _NCM_SPECTRAL_CHEB_PRIVATE_H_
#define _NCM_SPECTRAL_CHEB_PRIVATE_H_

#include <glib.h>

G_BEGIN_DECLS

/* *INDENT-OFF* */

#define NCM_SPECTRAL_CHEB_MAX_LEVEL ({MAX_LEVEL})

typedef void (*NcmSpectralChebKernel) {KERNEL_ARGUMENTS};
"""


def render_header() -> str:
    """Render the whole header.

    :return: text of the generated C header
    """
    kernels = "\n".join(render_kernel(level) for level in range(MAX_LEVEL + 1))
    entries = "\n".join(f"  {kernel_name(level)}," for level in range(MAX_LEVEL + 1))
    return f"""\
{PREAMBLE}
{kernels}
static const NcmSpectralChebKernel
_ncm_spectral_cheb_kernels[NCM_SPECTRAL_CHEB_MAX_LEVEL + 1] = {{
{entries}
}};

/* *INDENT-ON* */

G_END_DECLS

#endif /* _NCM_SPECTRAL_CHEB_PRIVATE_H_ */

"""


def main(argv: Sequence[str] | None = None) -> int:
    """Write the header, or compare it with the generated text.

    :param argv: command-line arguments, those of the process when ``None``
    :return: process exit status, 1 when ``--check`` finds a difference
    """
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
        help="header to write or to check (default: %(default)s)",
    )
    parser.add_argument(
        "--check",
        action="store_true",
        help="do not write; exit with status 1 if the header differs",
    )
    args = parser.parse_args(argv)
    header = render_header()
    if not args.check:
        args.output.write_text(header, encoding="utf-8")
        return 0
    if args.output.read_text(encoding="utf-8") == header:
        return 0
    print(f"{args.output} differs from the generated header", file=sys.stderr)
    return 1


if __name__ == "__main__":
    sys.exit(main())  # pragma: no cover

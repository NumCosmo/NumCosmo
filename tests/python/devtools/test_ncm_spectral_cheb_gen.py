"""Tests of the generator of the Chebyshev transform kernels of NcmSpectral."""

import math
import random
from pathlib import Path

import pytest

from tools import ncm_spectral_cheb_gen as gen

# The generated arithmetic agrees with the direct sums within this fraction of
# the l1 norm of the result; the largest value found in these tests is 1.8e-16.
TOLERANCE = 1.0e-15


def evaluate(program: gen.Program, inputs: list[float]) -> list[float]:
    """Run a program in double precision and return its outputs."""
    temporaries: list[float] = []

    def value(operand: gen.Operand) -> float:
        match operand:
            case gen.Input(row):
                return inputs[row]
            case gen.Temporary(index):
                return temporaries[index]

    for statement in program.statements:
        match statement:
            case gen.Sum(left, right):
                temporaries.append(value(left) + value(right))
            case gen.Difference(left, right):
                temporaries.append(value(left) - value(right))
            case gen.Product(operand, factor):
                temporaries.append(value(operand) * factor)
    return [value(output) for output in program.outputs]


def direct_cosine_transform(
    values: list[float], scale: float, dc_scale: float
) -> list[float]:
    """Evaluate the type-II cosine transform by its defining sum."""
    n_points = len(values)
    return [
        (dc_scale if k == 0 else scale)
        * math.fsum(
            v * math.cos(math.pi * ((2 * i + 1) * k % (4 * n_points)) / (2 * n_points))
            for i, v in enumerate(values)
        )
        for k in range(n_points)
    ]


def direct_chebyshev_coefficients(values: list[float]) -> list[float]:
    """Evaluate the Chebyshev coefficients of node values by the defining sum."""
    n = len(values) - 1
    coefficients = []
    for k in range(n + 1):
        terms = [values[0], values[n] if k % 2 == 0 else -values[n]]
        terms += [
            2.0 * values[j] * math.cos(math.pi * (j * k % (2 * n)) / n)
            for j in range(1, n)
        ]
        coefficients.append(math.fsum(terms) / (2 * n if k in (0, n) else n))
    return coefficients


def assert_close(result: list[float], expected: list[float]) -> None:
    """Assert agreement within TOLERANCE of the l1 norm of the expected values."""
    assert len(result) == len(expected)
    norm = math.fsum(abs(e) for e in expected)
    assert max(abs(r - e) for r, e in zip(result, expected)) <= TOLERANCE * norm


def node_values(level: int, seed: int) -> list[float]:
    """Sample a smooth function with random parameters at the nodes of a level."""
    rng = random.Random(seed)
    a, b, phi = rng.uniform(0.0, 3.0), rng.uniform(0.0, 5.0), rng.uniform(0.0, 6.0)
    n = 1 << level
    nodes = [math.cos(math.pi * j / n) for j in range(n + 1)]
    return [math.exp(a * s) * math.cos(b * s + phi) for s in nodes]


def test_program_statements() -> None:
    """Each builder method appends its statement and returns its temporary."""
    program = gen.Program()
    total = program.add(gen.Input(0), gen.Input(1))
    difference = program.subtract(gen.Input(0), total)
    product = program.scale(difference, 0.25)
    assert (total, difference, product) == (
        gen.Temporary(0),
        gen.Temporary(1),
        gen.Temporary(2),
    )
    assert program.statements == [
        gen.Sum(gen.Input(0), gen.Input(1)),
        gen.Difference(gen.Input(0), gen.Temporary(0)),
        gen.Product(gen.Temporary(1), 0.25),
    ]
    assert program.count(gen.Sum) == 1
    assert program.count(gen.Difference) == 1
    assert program.count(gen.Product) == 1


def test_program_scale_by_one() -> None:
    """A factor of one appends nothing and returns the operand."""
    program = gen.Program()
    assert program.scale(gen.Input(3), 1.0) == gen.Input(3)
    assert not program.statements


def test_odd_factors() -> None:
    """The factors are the scaled cosines of the odd quarter angles."""
    assert gen.odd_factors(2, 3.0) == [3.0 * gen.cos_rounded(math.pi / 4.0)]
    assert gen.odd_factors(4, 0.5) == [
        gen.cos_rounded(math.pi / 8.0),
        gen.cos_rounded(3.0 * math.pi / 8.0),
    ]


def test_cos_rounded() -> None:
    """The cosine is exact at zero and within one ulp of the libm cosine elsewhere."""
    assert gen.cos_rounded(0.0) == 1.0
    for n_points in [2**level for level in range(1, gen.MAX_LEVEL + 1)]:
        for i in range(n_points // 2):
            x = math.pi * (2 * i + 1) / (2.0 * n_points)
            assert abs(gen.cos_rounded(x) - math.cos(x)) <= math.ulp(math.cos(x))


@pytest.mark.parametrize("n_points", [1, 2, 4, 8, 16, 32, 64])
@pytest.mark.parametrize("scale, dc_scale", [(1.0, 1.0), (0.25, 0.125)])
def test_cosine_transform(n_points: int, scale: float, dc_scale: float) -> None:
    """The program agrees with the defining sum of the type-II transform."""
    rng = random.Random(n_points)
    values = [rng.uniform(-1.0, 1.0) for _ in range(n_points)]
    program = gen.Program()
    inputs: list[gen.Operand] = [gen.Input(i) for i in range(n_points)]
    program.outputs = gen.cosine_transform(program, inputs, scale, dc_scale)
    assert_close(
        evaluate(program, values), direct_cosine_transform(values, scale, dc_scale)
    )


@pytest.mark.parametrize("level", range(gen.MAX_LEVEL + 1))
@pytest.mark.parametrize("seed", [1, 2, 3])
def test_kernel_program(level: int, seed: int) -> None:
    """The kernel program agrees with the defining sum of the coefficients."""
    values = node_values(level, seed)
    program = gen.kernel_program(level)
    assert len(program.outputs) == (1 << level) + 1
    assert_close(evaluate(program, values), direct_chebyshev_coefficients(values))


@pytest.mark.parametrize("level", range(gen.MAX_LEVEL + 1))
def test_kernel_program_constant(level: int) -> None:
    """A constant has the coefficient y_0 alone, without rounding."""
    program = gen.kernel_program(level)
    n = 1 << level
    assert evaluate(program, [1.0] * (n + 1)) == [1.0] + [0.0] * n


def test_kernel_operation_counts() -> None:
    """The operation counts are those recorded for the measured kernels."""
    counts = []
    for level in range(3, gen.MAX_LEVEL + 1):
        program = gen.kernel_program(level)
        counts.append(
            (
                program.count(gen.Product),
                program.count(gen.Sum) + program.count(gen.Difference),
            )
        )
    assert counts == [(11, 30), (27, 76), (67, 190), (163, 464), (387, 1106)]


def test_render_statement() -> None:
    """Each statement kind renders as the definition of its temporary."""
    x0, t1 = gen.Input(0), gen.Temporary(1)
    assert gen.render_operand(x0) == "x[0 * L + c]"
    assert gen.render_operand(t1) == "t1"
    assert gen.render_statement(2, gen.Sum(x0, t1)) == (
        "const gdouble t2 = x[0 * L + c] + t1;"
    )
    assert gen.render_statement(3, gen.Difference(t1, x0)) == (
        "const gdouble t3 = t1 - x[0 * L + c];"
    )
    assert gen.render_statement(4, gen.Product(t1, 0.5)) == (
        "const gdouble t4 = t1 * 0.5;"
    )


def test_render_body_stores_at_definition() -> None:
    """Each output row is stored right after the statement that defines it."""
    program = gen.Program()
    total = program.add(gen.Input(0), gen.Input(1))
    difference = program.subtract(gen.Input(0), gen.Input(1))
    program.outputs = [difference, total, difference]
    assert gen.render_body(program) == [
        "const gdouble t0 = x[0 * L + c] + x[1 * L + c];",
        "y[1 * L + c] = t0;",
        "const gdouble t1 = x[0 * L + c] - x[1 * L + c];",
        "y[0 * L + c] = t1;",
        "y[2 * L + c] = t1;",
    ]


def test_render_body_rejects_input_output() -> None:
    """An output that is an input row is an error."""
    program = gen.Program(outputs=[gen.Input(0)])
    with pytest.raises(TypeError, match="output 0 is not the result"):
        gen.render_body(program)


def test_render_kernel_level_zero() -> None:
    """The kernel of level zero is the half sum and the half difference."""
    assert gen.render_kernel(0) == f"""\
/*
 * Level 0, 2 nodes. Operations per column: 2 multiplications,
 * 2 additions.
 */
static void
_ncm_spectral_cheb_kernel_0 {gen.KERNEL_ARGUMENTS}
{{
  gsize c;

  #pragma omp simd
  for (c = 0; c < L; c++)
  {{
    const gdouble t0 = x[0 * L + c] + x[1 * L + c];
    const gdouble t1 = t0 * 0.5;
    y[0 * L + c] = t1;
    const gdouble t2 = x[0 * L + c] - x[1 * L + c];
    const gdouble t3 = t2 * 0.5;
    y[1 * L + c] = t3;
  }}
}}
"""


def test_render_header() -> None:
    """The header has one kernel and one table entry per level."""
    header = gen.render_header()
    assert header == gen.render_header()
    assert header.count("#pragma omp simd") == gen.MAX_LEVEL + 1
    for level in range(gen.MAX_LEVEL + 1):
        name = gen.kernel_name(level)
        assert f"\n{name} (" in header
        assert f"\n  {name},\n" in header
    assert f"#define NCM_SPECTRAL_CHEB_MAX_LEVEL ({gen.MAX_LEVEL})" in header
    assert header.startswith("/*") and header.endswith("_PRIVATE_H_ */\n\n")


def test_main_writes_and_checks(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    """The script writes the header, and --check reports a modified one."""
    output = tmp_path / gen.HEADER_NAME
    assert gen.main(["--output", str(output)]) == 0
    assert output.read_text(encoding="utf-8") == gen.render_header()
    assert gen.main(["--output", str(output), "--check"]) == 0
    output.write_text(output.read_text(encoding="utf-8") + "\n", encoding="utf-8")
    assert gen.main(["--output", str(output), "--check"]) == 1
    assert "differs from the generated header" in capsys.readouterr().err


def test_committed_header_is_current() -> None:
    """The header in the source tree is the output of the generator."""
    assert gen.main(["--check"]) == 0

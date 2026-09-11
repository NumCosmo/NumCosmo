#!/usr/bin/env python
#
# test_xcor_cls_app.py
#
# Thu Sep 10 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_xcor_cls_app.py
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

"""Tests for the xcor cls CLI command."""

from pathlib import Path

import pytest
from typer.testing import CliRunner

import matplotlib

matplotlib.use("Agg")

# flake8: noqa: E402
# pylint: disable=wrong-import-position

import numpy as np
from matplotlib.figure import Figure
import matplotlib.pyplot as plt

from numcosmo_py import Nc, Ncm
from numcosmo_py.app import app
from numcosmo_py.app.xcor.cls import ComputeCls, EllSpacing, sample_ells
from numcosmo_py.app.xcor.common import XcorClosureOption, contiguous_runs
from numcosmo_py.app.xcor.plotting import percent_tick, style_ratio_axis

pytestmark = pytest.mark.app
runner = CliRunner()

Ncm.cfg_init()

# The analytic shapes carry a closed form and no photo-z machinery, so a whole
# spectrum over them costs about as much as one physical kernel's preparation.
TOPHAT_A = "radial-tophat chi_lower=600 chi_upper=700"
TOPHAT_B = "radial-tophat chi_lower=500 chi_upper=600"
TOPHAT_C = "radial-tophat chi_lower=400 chi_upper=500"


def _cls(*extra: str, kernels: tuple[str, ...] = (TOPHAT_A,)):
    """Run `xcor cls` over a short multipole range, without displaying."""
    args = ["xcor", "cls"]
    for kernel in kernels:
        args += ["--kernel", kernel]
    args += ["--ell-min", "10", "--ell-max", "20", "--n-ell", "3", "--no-show-plot"]
    return runner.invoke(app, [*args, *extra])


def test_cls_runs_and_writes_the_figure(tmp_path: Path) -> None:
    """The command computes a spectrum and saves it where it was told to."""
    output = tmp_path / "cls.png"
    result = _cls("--output", output.as_posix())

    assert result.exit_code == 0, result.output
    assert "1 kernel(s), 1 spectra" in result.output
    assert "Angular power spectra complete" in result.output
    assert output.exists()


def test_cls_defaults_to_the_exact_chebyshev_pair() -> None:
    """Unset, the command takes the quadrature that cannot fail to converge.

    'exact' integrates the closures on the common refinement of their knots and
    needs no tolerance, so an unattended run over a wide multipole range cannot
    abort partway for missing one.
    """
    result = _cls()

    assert result.exit_code == 0, result.output
    assert "method=exact, closure=chebyshev" in result.output


@pytest.mark.parametrize(
    ("kernels", "n_spectra"),
    [((TOPHAT_A,), 1), ((TOPHAT_A, TOPHAT_B), 3), ((TOPHAT_A, TOPHAT_B, TOPHAT_C), 6)],
)
def test_cls_pairs_every_kernel_by_default(
    kernels: tuple[str, ...], n_spectra: int
) -> None:
    """Every auto- and cross-pair is computed: n(n+1)/2 spectra for n kernels."""
    result = _cls(kernels=kernels)

    assert result.exit_code == 0, result.output
    assert f"{len(kernels)} kernel(s), {n_spectra} spectra" in result.output


@pytest.mark.parametrize(
    "kernels", [(TOPHAT_A,), (TOPHAT_A, TOPHAT_B), (TOPHAT_A, TOPHAT_B, TOPHAT_C)]
)
def test_cls_no_cross_keeps_the_auto_spectra_only(kernels: tuple[str, ...]) -> None:
    """--no-cross drops the cross-pairs, leaving one spectrum per kernel."""
    result = _cls("--no-cross", kernels=kernels)

    assert result.exit_code == 0, result.output
    assert f"{len(kernels)} kernel(s), {len(kernels)} spectra" in result.output


def test_view_kernel_no_cross_reaches_the_cls_run() -> None:
    """The pairing option is shared, so it steers `kernel view --cls` too."""
    result = runner.invoke(
        app,
        [
            "xcor",
            "kernel",
            "view",
            "--kernel",
            TOPHAT_A,
            "--kernel",
            TOPHAT_B,
            "--ell",
            "10",
            "--n-ell",
            "2",
            "--no-show-plot",
            "--cls",
            "--no-cross",
        ],
    )

    assert result.exit_code == 0, result.output
    assert "2 kernel(s), 2 spectra" in result.output


def test_cls_goes_past_the_kernel_block_cap() -> None:
    """The whole point: more multipoles than a kernel evaluates in one block.

    The kernel view hands its whole range to one get_eval_vectorized() call and
    is capped at NC_XCOR_KERNEL_MAX_ELL_BLOCK for it; this command goes through
    NcXcorSolver, which tiles the range into blocks of its own.
    """
    n_ell = Nc.XCOR_KERNEL_MAX_ELL_BLOCK + 20
    result = runner.invoke(
        app,
        [
            "xcor",
            "cls",
            "--kernel",
            TOPHAT_A,
            "--ell-min",
            "2",
            "--ell-max",
            str(n_ell + 1),
            "--ell-spacing",
            "all",
            "--no-show-plot",
        ],
    )

    assert result.exit_code == 0, result.output
    assert f"for {n_ell} multipole(s)" in result.output


def test_view_kernel_still_caps_its_multipole_block() -> None:
    """The kernel view keeps the cap the library imposes on it."""
    result = runner.invoke(
        app,
        [
            "xcor",
            "kernel",
            "view",
            "--kernel",
            TOPHAT_A,
            "--ell",
            "2",
            "--n-ell",
            str(Nc.XCOR_KERNEL_MAX_ELL_BLOCK + 1),
            "--no-show-plot",
        ],
    )

    assert result.exit_code != 0
    assert "Kernel evaluation complete" not in result.output


def test_cls_rejects_an_empty_range() -> None:
    """An upper bound below the lower one has nothing to compute."""
    result = runner.invoke(
        app,
        ["xcor", "cls", "--ell-min", "100", "--ell-max", "10", "--no-show-plot"],
    )

    assert result.exit_code != 0


def test_cls_rejects_log_spacing_from_zero() -> None:
    """Log spacing has no first sample at ell = 0."""
    result = runner.invoke(
        app,
        [
            "xcor",
            "cls",
            "--ell-min",
            "0",
            "--ell-max",
            "10",
            "--ell-spacing",
            "log",
            "--no-show-plot",
        ],
    )

    assert result.exit_code != 0


@pytest.mark.parametrize("spacing", ["log", "linear", "all"])
def test_cls_every_spacing_runs(spacing: str) -> None:
    """All three layouts reach the solver."""
    result = _cls("--ell-spacing", spacing)

    assert result.exit_code == 0, result.output
    assert "Angular power spectra complete" in result.output


@pytest.mark.parametrize(
    ("ells", "runs"),
    [
        ([5], [(5, 5)]),
        ([5, 6, 7], [(5, 7)]),
        ([5, 7, 9], [(5, 5), (7, 7), (9, 9)]),
        ([5, 6, 9, 10, 11, 20], [(5, 6), (9, 11), (20, 20)]),
    ],
)
def test_contiguous_runs(ells: list[int], runs: list[tuple[int, int]]) -> None:
    """A sampled multipole set is split into the ranges the solver can ask for."""
    assert contiguous_runs(np.array(ells)) == runs


def test_sample_ells_all_covers_the_range() -> None:
    """'all' spacing ignores --n-ell and takes every multipole."""
    ells = sample_ells(2, 12, 3, EllSpacing.ALL)

    assert ells.tolist() == list(range(2, 13))


def test_sample_ells_clamps_to_the_range_width() -> None:
    """Asking for more samples than the range holds gives the whole range."""
    ells = sample_ells(2, 6, 50, EllSpacing.LOG)

    assert ells.tolist() == [2, 3, 4, 5, 6]


@pytest.mark.parametrize("spacing", [EllSpacing.LOG, EllSpacing.LINEAR])
def test_sample_ells_is_sorted_unique_and_in_range(spacing: EllSpacing) -> None:
    """Rounding may collapse neighbours, but never reorders or escapes."""
    ells = sample_ells(2, 1000, 64, spacing)

    assert ells[0] == 2
    assert ells[-1] == 1000
    assert np.all(np.diff(ells) > 0)
    assert len(ells) <= 64


def test_sample_ells_log_is_denser_at_low_ell() -> None:
    """Log spacing is what makes a wide range affordable."""
    ells = sample_ells(2, 1000, 32, EllSpacing.LOG)
    half = len(ells) // 2

    assert ells[half] < np.sqrt(2 * 1000) * 2


def test_sampled_multipoles_match_the_contiguous_solve() -> None:
    """A gappy multipole set gets the same values as the full range does.

    Sampled multipoles are asked for one contiguous run at a time and the runs
    are stitched back together in order, so this covers both the request
    splitting and the reassembly against a solve that needs neither.
    """
    command = ComputeCls(
        kernel=[TOPHAT_A, TOPHAT_B],
        ell_min=10,
        ell_max=12,
        ell_spacing=EllSpacing.ALL,
        show_plot=False,
    )

    # Two kernels of one shape would otherwise share a label, and the index is
    # parenthesised rather than '#'-prefixed to survive a LaTeX-rendered legend.
    assert [label for label, _ in command.kernels] == ["Top-hat (1)", "Top-hat (2)"]

    dense_ells = np.arange(10, 17)
    dense = command._compute_cls(  # pylint: disable=protected-access
        dense_ells, -1, XcorClosureOption.SPLINE
    )

    sparse_ells = np.array([10, 13, 16])
    sparse = command._compute_cls(  # pylint: disable=protected-access
        sparse_ells, -1, XcorClosureOption.SPLINE
    )

    picked = np.searchsorted(dense_ells, sparse_ells)

    assert set(dense) == {(0, 0), (0, 1), (1, 1)}
    for pair, cl in sparse.items():
        assert cl.shape == sparse_ells.shape
        assert_allclose = np.testing.assert_allclose
        assert_allclose(cl, dense[pair][picked], rtol=1e-10)


def test_cls_compare_limber_draws_both_spectra_on_top(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """--compare-limber draws the Limber spectrum dashed beside the solid one.

    The ratio panel alone hides where the two curves sit, so the top panel
    carries both: one solid and one dashed line per pair, in the same color.
    """
    captured: list[Figure] = []
    original = Figure.savefig

    def _capture(self, *args, **kwargs):
        captured.append(self)
        return original(self, *args, **kwargs)

    monkeypatch.setattr(Figure, "savefig", _capture)

    output = tmp_path / "cls.png"
    result = _cls("--compare-limber", "--output", output.as_posix())

    assert result.exit_code == 0, result.output
    assert output.exists()
    assert len(captured) == 1
    top, bottom = captured[0].axes
    styles = [line.get_linestyle() for line in top.get_lines()]
    assert styles == ["-", "--"]
    labels = [str(line.get_label()) for line in top.get_lines()]
    assert labels[0].endswith("(Non-Limber)")
    assert labels[1].endswith("(Limber)")
    assert [line.get_color() for line in top.get_lines()][0] == (
        [line.get_color() for line in top.get_lines()][1]
    )
    # The ratio panel: one line per pair plus the zero guide; the percent
    # levels come from the symlog ticks, not from guide lines that would widen
    # the axis when the agreement is good.
    assert len(bottom.get_lines()) == 2
    assert bottom.get_yscale() == "symlog"
    fmt = bottom.yaxis.get_major_formatter()
    assert [fmt(v) for v in (0.0, 1.0e-4, 0.01, -0.1, 1.0)] == [
        "0",
        "0.01%",
        "1%",
        "-10%",
        "100%",
    ]
    # Limber over-predicts this top hat at ell 10-20, so every deviation is
    # positive and the negative side closes at the linear core.
    ratio = np.asarray(bottom.get_lines()[0].get_ydata(), dtype=float)
    assert np.all(ratio > 0.0)
    ymin, ymax = bottom.get_ylim()
    assert -1.0e-4 <= ymin < 0.0
    assert ratio.max() < ymax < 10.0 * ratio.max()


def test_style_ratio_axis_without_finite_deviations() -> None:
    """A ratio panel with no finite deviation keeps matplotlib's own limits.

    Every C_ell of a pair being zero makes the deviation NaN throughout; the
    axis is still styled, and the limits are left to matplotlib.
    """
    fig, ax = plt.subplots()
    ax.axhline(0.0)
    before = ax.get_ylim()
    style_ratio_axis(ax, [np.full(4, np.nan)])
    assert ax.get_yscale() == "symlog"
    assert ax.get_ylim() == before
    plt.close(fig)


def test_style_ratio_axis_limits_follow_both_signs() -> None:
    """Deviations of both signs open both sides of the axis by a quarter decade."""
    fig, ax = plt.subplots()
    ratio = np.array([-0.02, 1.0e-6, 0.3])
    ax.plot(ratio)
    style_ratio_axis(ax, [ratio])
    ymin, ymax = ax.get_ylim()
    assert ymin == pytest.approx(-0.02 * 10.0**0.25)
    assert ymax == pytest.approx(0.3 * 10.0**0.25)
    assert percent_tick(-0.02) == "-2%"
    plt.close(fig)

#
# plotting.py
#
# Fri Sep 11 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# plotting.py
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

"""Axis styling shared by the xcor plotting commands."""

from collections.abc import Iterable

import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter, SymmetricalLogLocator
import numpy as np


def percent_tick(value: float, _pos: int | None = None) -> str:
    """Format a fractional deviation as a percentage tick label.

    :param value: Deviation as a fraction (0.01 is one percent).
    :param _pos: Tick position, unused.
    :return: The label, "0" at zero and e.g. "-1%" or "0.01%" elsewhere.
    """
    if value == 0.0:
        return "0"
    return f"{value * 100.0:g}%"


def style_ratio_axis(
    ax: plt.Axes, ratios: Iterable[np.ndarray], linthresh: float = 1.0e-4
) -> None:
    """Put a fractional-deviation axis on a symmetric log scale with percent ticks.

    Deviations between two spectra span O(1) where Limber fails and 1e-3 or
    below where it holds, so the axis is symmetric log with a linear core of
    @linthresh (0.01% by default) that keeps the sign visible near agreement.
    Major ticks sit on the powers of ten, labelled in percent, and carry the
    darker grid lines: the 1% and 10% levels read off the panel without guide
    lines, which would stretch the axis when the agreement is better than that.

    The limits follow the data: a side with no deviation of its sign is closed
    just past zero, inside the linear core, so the zero line stays visible
    without a decade of empty space. The core itself is drawn at half a decade
    (``linscale``) so it does not dominate a panel whose deviations are small.

    :param ax: Axis to style; the ratio curves must already be drawn on it.
    :param ratios: The plotted deviations, used only for the limits.
    :param linthresh: Half-width of the linear core around zero.
    """
    finite = np.concatenate([np.asarray(r, dtype=float).ravel() for r in ratios])
    finite = finite[np.isfinite(finite)]
    ax.set_yscale("symlog", linthresh=linthresh, linscale=0.5)
    ax.yaxis.set_major_locator(SymmetricalLogLocator(linthresh=linthresh, base=10.0))
    ax.yaxis.set_minor_locator(
        SymmetricalLogLocator(
            linthresh=linthresh, base=10.0, subs=[float(i) for i in range(2, 10)]
        )
    )
    ax.yaxis.set_major_formatter(FuncFormatter(percent_tick))
    ax.grid(True, which="major", alpha=0.6)
    ax.grid(True, which="minor", alpha=0.15)
    if finite.size == 0:
        return
    margin = 10.0**0.25
    low, high = float(finite.min()), float(finite.max())
    sliver = 0.2 * linthresh
    bottom = -sliver if low >= -sliver else low * margin
    top = sliver if high <= sliver else high * margin
    ax.set_ylim(bottom, top)

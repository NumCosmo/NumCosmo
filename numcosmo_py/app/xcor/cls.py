#
# cls.py
#
# Thu Sep 10 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# cls.py
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

"""CLI command for computing cross-correlation angular power spectra."""

import dataclasses
import enum
from typing import Annotated

import numpy as np
import typer

from .common import XcorKernelCommon


class EllSpacing(str, enum.Enum):
    """How the multipoles are laid out inside the requested range."""

    ALL = "all"
    LINEAR = "linear"
    LOG = "log"


def sample_ells(lmin: int, lmax: int, n_ell: int, spacing: EllSpacing) -> np.ndarray:
    """Choose the multipoles to compute inside an inclusive range.

    Rounding two neighbouring samples onto the same integer collapses them, so
    the returned array can be shorter than @n_ell -- densely sampled low
    multipoles are the common case.

    :param lmin: First multipole of the range.
    :param lmax: Last multipole of the range, inclusive.
    :param n_ell: Number of samples wanted, ignored for 'all' spacing.
    :param spacing: Layout of the samples.
    :return: Strictly increasing integer array of multipoles.
    """
    if spacing is EllSpacing.ALL or n_ell >= lmax - lmin + 1:
        return np.arange(lmin, lmax + 1)

    if spacing is EllSpacing.LOG:
        samples = np.logspace(np.log10(lmin), np.log10(lmax), n_ell)
    else:
        samples = np.linspace(lmin, lmax, n_ell)

    return np.unique(np.rint(samples).astype(int))


@dataclasses.dataclass(kw_only=True)
class ComputeCls(XcorKernelCommon):
    r"""Compute and plot cross-correlation angular power spectra.

    The kernels are specified exactly as for 'numcosmo xcor kernel view', and
    every option that command shares with this one means the same thing. Where
    the kernel view draws the projected kernel and is bounded by the multipole
    block a kernel evaluates in one call, this command goes through
    NcXcorSolver, which tiles the requested range into blocks itself: any range
    can be asked for.

    All auto- and cross-spectra of the given kernels are computed; pass
    --no-cross for the auto-spectra alone.

    Examples:
        # 64 log-spaced multipoles from 2 to 1000, one kernel
        numcosmo xcor cls --kernel "cmb_lensing lmax=3000"

        # Every multipole from 2 to 200, auto-spectra only
        numcosmo xcor cls \\
            --kernel "number-counts survey=LSST-Y1 bin_idx=0 bias=1.5" \\
            --kernel "number-counts survey=LSST-Y1 bin_idx=1 bias=1.6" \\
            --ell-min 2 --ell-max 200 --ell-spacing all --no-cross

        # Non-Limber against Limber over a wide range
        numcosmo xcor cls \\
            --kernel "radial-tophat chi_lower=600 chi_upper=700" \\
            --kernel "radial-tophat chi_lower=500 chi_upper=600" \\
            --ell-min 10 --ell-max 2000 --n-ell 100 \\
            --compare-limber --output cls.png --no-show-plot
    """

    ell_min: Annotated[
        int,
        typer.Option(
            help="First multipole of the range.",
            show_default=True,
            min=0,
        ),
    ] = 2

    ell_max: Annotated[
        int,
        typer.Option(
            help="Last multipole of the range, inclusive.",
            show_default=True,
            min=0,
        ),
    ] = 1000

    n_ell: Annotated[
        int,
        typer.Option(
            help=(
                "Number of multipoles to sample inside the range. Ignored for "
                "--ell-spacing all, and clamped to the range's width. Samples "
                "rounding onto the same multipole collapse into one, so fewer "
                "may be computed than asked for."
            ),
            show_default=True,
            min=1,
        ),
    ] = 64

    ell_spacing: Annotated[
        EllSpacing,
        typer.Option(
            help=(
                "Layout of the sampled multipoles: 'log' or 'linear' for "
                "--n-ell samples, 'all' for every multipole in the range. "
                "Only consecutive multipoles share a solver block, so a "
                "sampled range costs about one block per sample."
            ),
            show_default=True,
        ),
    ] = EllSpacing.LOG

    def __post_init__(self) -> None:
        """Execute the angular power spectrum command.

        :raises typer.BadParameter: If the multipole range is empty or cannot
            be sampled the requested way.
        """
        if self.ell_max < self.ell_min:
            raise typer.BadParameter(
                f"--ell-max ({self.ell_max}) must not be below --ell-min "
                f"({self.ell_min})."
            )
        if self.ell_spacing is EllSpacing.LOG and self.ell_min < 1:
            raise typer.BadParameter(
                "--ell-spacing log needs --ell-min >= 1; use linear or all "
                "spacing to include ell = 0."
            )

        super().__post_init__()

        ells = sample_ells(self.ell_min, self.ell_max, self.n_ell, self.ell_spacing)
        self.compute_and_plot_cls(ells)

        print()
        print("[OK] Angular power spectra complete!")

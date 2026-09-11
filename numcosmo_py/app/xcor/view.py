#
# view.py
#
# Wed Mar 12 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# view.py
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

"""CLI command for viewing cross-correlation kernels."""

import dataclasses
from typing import Annotated, Optional
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D
import numpy as np
import typer

from numcosmo_py import Nc

from .common import XcorKernelCommon
from .kernels import get_kernel_help_text, get_kernel_registry_help_text


@dataclasses.dataclass(kw_only=True)
class KernelEvaluation:
    """Result of evaluating a kernel on a k-grid.

    :ivar name: Kernel name.
    :ivar method: Evaluation method (Limber or Non-Limber).
    :ivar kernel: NumCosmo XcorKernel object.
    :ivar evaluator: NumCosmo XcorKernelIntegrand evaluator.
    :ivar RH_Mpc: Hubble radius in Mpc.
    :ivar lmin: Minimum multipole value.
    """

    name: str
    method: str
    kernel: Nc.XcorKernel
    evaluator: Nc.XcorKernelIntegrand
    RH_Mpc: float
    lmin: int

    def range(self) -> tuple[float, float]:
        """Return k range in Mpc^-1.

        :return: Tuple of (k_min, k_max) in Mpc^-1.
        """
        kmin, kmax = self.evaluator.get_range()
        return kmin / self.RH_Mpc, kmax / self.RH_Mpc

    def evaluate(self, k_Mpc: np.ndarray) -> np.ndarray:
        """Evaluate kernel on a k grid.

        :param k_Mpc: Wave number array in Mpc^-1.
        :return: Kernel values evaluated at k_Mpc. Shape is (n_k, n_ell) if n_ell > 1,
            else (n_k,).
        """
        k = k_Mpc * self.RH_Mpc
        result = np.array([self.evaluator.eval_array(ki) for ki in k])
        # If single ell, flatten to 1D for backward compatibility
        if self.evaluator.len == 1:
            return result[:, 0]
        return result

    def range_intersection(self, *others: "KernelEvaluation") -> tuple[float, float]:
        """Get common k range intersection with other kernels.

        :param others: Other KernelEvaluation objects to intersect with.
        :return: Tuple of (k_min, k_max) representing the common range.
        :raises ValueError: If there is no overlapping k range between kernels.
        """
        kmin, kmax = self.range()
        for other in others:
            other_kmin, other_kmax = other.range()
            kmin = max(kmin, other_kmin)
            kmax = min(kmax, other_kmax)
        if kmin >= kmax:
            raise ValueError("No overlapping k range between kernels.")
        return kmin, kmax

    def plot(
        self,
        ax: plt.Axes,
        n_points: int,
        k_range: tuple[float, float] | None = None,
        color: str = "blue",
        linestyle: str = "-",
        linewidth: float = 2.0,
        alpha: float = 1.0,
    ) -> str:
        """Plot kernel on given axes.

        :param ax: Matplotlib axes object.
        :param k_range: Optional k range (k_min, k_max) in Mpc^-1 for plotting.
        :param color: Line color.
        :param linestyle: Line style.
        :param linewidth: Line width.
        :param alpha: Line transparency (0-1).
        :return: Label for the legend entry.
        """
        kmin, kmax = self.range()
        if k_range is not None:
            kmin = max(kmin, k_range[0])
            kmax = min(kmax, k_range[1])
        k_Mpc = np.logspace(np.log10(kmin), np.log10(kmax), n_points)
        kernel_values = self.evaluate(k_Mpc)
        k3_2 = k_Mpc ** (3 / 2)

        if self.evaluator.len == 1:
            # Single ell: plot one line
            ax.plot(
                k_Mpc,
                k3_2 * kernel_values,
                color=color,
                linestyle=linestyle,
                linewidth=linewidth,
                alpha=alpha,
            )
            label = f"{self.name} ({self.method}), $\\ell={self.lmin}$"
        else:
            # Multiple ells: use shades of the same color
            # Create color shades from darker to lighter
            base_rgb = mcolors.to_rgb(color)
            for i in range(self.evaluator.len):
                # Interpolate between 0.5 (darker) and 1.3 (lighter) of base color
                # Clamp to valid RGB range [0, 1]
                factor = 0.5 + (0.8 * i / max(1, self.evaluator.len - 1))
                shade_rgb = tuple(min(1.0, c * factor) for c in base_rgb)

                ax.plot(
                    k_Mpc,
                    k3_2 * kernel_values[:, i],
                    color=shade_rgb,
                    linestyle=linestyle,
                    linewidth=linewidth,
                    alpha=alpha,
                )
            # Single legend entry for all ells
            label = (
                f"{self.name} ({self.method}), "
                f"$\\ell={self.lmin}$-${self.lmin + self.evaluator.len - 1}$"
            )

        ax.set_xscale("log")
        ax.set_yscale("symlog", linthresh=1e-10)
        return label

    def plot_comparison(
        self,
        other: "KernelEvaluation",
        ax: plt.Axes,
        k_range: tuple[float, float] | None = None,
        color: str = "red",
        linestyle: str = "--",
    ) -> None:
        """Plot comparison of this kernel with another.

        :param other: Another KernelEvaluation object to compare with.
        :param ax: Matplotlib axes object.
        :param k_range: Optional k range (k_min, k_max) in Mpc^-1 for plotting.
        :param color: Line color.
        :param linestyle: Line style.
        """
        kmin, kmax = self.range_intersection(other)
        if k_range is not None:
            kmin = max(kmin, k_range[0])
            kmax = min(kmax, k_range[1])
        k_Mpc = np.logspace(np.log10(kmin), np.log10(kmax), 1000)
        kernel_values = self.evaluate(k_Mpc)
        other_values = other.evaluate(k_Mpc)

        if self.evaluator.len == 1:
            # Single ell comparison
            ratio = other_values / kernel_values - 1.0
            ax.loglog(k_Mpc, np.abs(ratio), color=color, linestyle=linestyle)
        else:
            # Multiple ells: use shades of the same color
            base_rgb = mcolors.to_rgb(color)
            for i in range(self.evaluator.len):
                # Create color shades from darker to lighter
                factor = 0.5 + (0.8 * i / max(1, self.evaluator.len - 1))
                shade_rgb = tuple(min(1.0, c * factor) for c in base_rgb)

                ratio = other_values[:, i] / kernel_values[:, i] - 1.0
                ax.loglog(
                    k_Mpc,
                    np.abs(ratio),
                    color=shade_rgb,
                    linestyle=linestyle,
                )


@dataclasses.dataclass(kw_only=True)
class KernelVariants:
    """Container for kernel evaluation results.

    :ivar main: Main kernel evaluation result.
    :ivar alternative: Alternative kernel evaluation result (e.g., Limber
        approximation).
    """

    main: KernelEvaluation
    alternative: KernelEvaluation | None = None


def ListKernels(
    kernel_type: Annotated[
        Optional[str],
        typer.Argument(
            help=(
                "Kernel type to document. Given one, every parameter it takes "
                "is listed with its type, default, range and meaning. Left out, "
                "all kernel types are summarised."
            ),
            show_default=False,
        ),
    ] = None,
) -> None:
    """List the available cross-correlation kernel types, or document one.

    Displays a formatted table showing all available kernel types, the NumCosmo
    object each builds, and its parameters -- or, given a kernel type, the full
    documentation of that one's parameters.

    :param kernel_type: A kernel type to document, or None for the summary.
    :raises typer.BadParameter: If the kernel type is not recognized.
    """
    if kernel_type is not None:
        try:
            print()
            print(get_kernel_help_text(kernel_type))
            print()
        except ValueError as error:
            raise typer.BadParameter(str(error)) from error
        return

    print("\nAvailable Cross-Correlation Kernel Types\n")
    print(get_kernel_registry_help_text())
    print(
        '\nUsage: numcosmo xcor kernel view --kernel "<kernel_type> param=value ..."'
        '\n       numcosmo xcor cls         --kernel "<kernel_type> param=value ..."'
        "\n\nRun 'numcosmo xcor kernel list <kernel_type>' for what each "
        "parameter means."
    )
    print()


@dataclasses.dataclass(kw_only=True)
class ViewKernel(XcorKernelCommon):
    """View cross-correlation kernels.

    This command visualizes cross-correlation kernels used in cosmological
    analyses. By default, kernels are evaluated using the non-Limber method
    and plotted with full lines. Use --compare-limber to also show the
    Limber approximation (with thinner, semi-transparent dashed lines) for
    comparison.

    The multipoles are handed to the kernel in one block, so --n-ell is capped
    at NC_XCOR_KERNEL_MAX_ELL_BLOCK. Use 'numcosmo xcor cls' for a spectrum
    over an arbitrarily wide multipole range.

    The kernel specification follows the format:

        "<kernel_name> key=value key=value ..."

    Available kernel types:
        - cmb_lensing: CMB lensing convergence kernel
        - cmb_isw: CMB Integrated Sachs-Wolfe kernel
        - tsz: Thermal Sunyaev-Zeldovich kernel
        - number-counts: Galaxy number counts kernel
        - weak-lensing: Weak lensing kernel

    Examples:
        # View CMB lensing kernel (non-Limber by default)
        numcosmo xcor kernel view --kernel "cmb_lensing lmax=3000" --ell 20

        # View Galaxy clustering kernel (LSST Y1 lens bin 0)
        numcosmo xcor kernel view \\
            --kernel "number-counts survey=LSST-Y1 bin_idx=0 bias=1.5" \\
            --ell 100

        # Compare non-Limber vs Limber for weak lensing
        numcosmo xcor kernel view \\
            --kernel "weak-lensing survey=LSST-Y10 bin_idx=2" \\
            --ell 50 \\
            --compare-limber \\
            --output wl_comparison.png
    """

    ell: Annotated[
        int,
        typer.Option(
            help="Minimum multipole value for kernel evaluation.",
            show_default=True,
            min=0,
        ),
    ] = 20

    n_ell: Annotated[
        int,
        typer.Option(
            help=(
                "Number of multipole values to evaluate (starting from ell). "
                "Capped at the largest block the kernel evaluates in one call; "
                "'numcosmo xcor cls' has no such cap."
            ),
            show_default=True,
            min=1,
            max=Nc.XCOR_KERNEL_MAX_ELL_BLOCK,
        ),
    ] = 1

    k_range: Annotated[
        tuple[float, float] | None,
        typer.Option(
            help="k value range for plotting (Mpc^-1).",
            show_default=True,
        ),
    ] = None

    n_points: Annotated[
        int,
        typer.Option(
            help="Number of k points to evaluate.",
            show_default=True,
            min=100,
            max=50000,
        ),
    ] = 5000

    cls: Annotated[
        bool,
        typer.Option(
            help=(
                "Also compute and plot the angular power spectra C_ell for every "
                "auto- and cross-pair of the requested kernels, over the multipole "
                "range set by --ell/--n-ell. Use 'numcosmo xcor cls' for a "
                "spectrum over a wider range than --n-ell allows."
            ),
            show_default=True,
        ),
    ] = False

    def __post_init__(self) -> None:
        """Execute the kernel view command.

        This method parses the kernel specification, creates the kernel
        objects, evaluates them, and generates visualization plots.

        :raises ValueError: If the kernel specification is invalid.
        :raises RuntimeError: If kernel evaluation fails.
        """
        if self.k_range is not None:
            if self.k_range[0] <= 0 or self.k_range[1] <= 0:
                raise ValueError("k range values must be positive")
            if self.k_range[0] >= self.k_range[1]:
                raise ValueError("k range min must be less than max")

        super().__post_init__()

        if self.n_ell > 1:
            print(
                f"  [OK] Evaluating {self.n_ell} multipoles: ell = {self.ell} "
                f"to {self.ell + self.n_ell - 1}"
            )
        else:
            print(f"  [OK] Evaluating single multipole: ell = {self.ell}")
        print()

        kernel_evals = []
        for kernel_label, kernel_obj in self.kernels:
            kernel_evals += self._evaluate_kernels(kernel_label, kernel_obj)

        # Plot results
        self._plot_results(kernel_evals)

        if self.cls:
            self.compute_and_plot_cls(np.arange(self.ell, self.ell + self.n_ell))

        print()
        print("[OK] Kernel visualization complete!")

    def _cls_output_path(self) -> Optional[Path]:
        """Where the C_ell figure is written.

        The kernel figure already owns --output, so the spectra go beside it
        under a suffixed name.

        :return: The suffixed output path, or None when nothing is to be saved.
        """
        if self.output is None:
            return None
        return self.output.with_name(f"{self.output.stem}_cls{self.output.suffix}")

    def _evaluate_kernels(
        self, kernel_label: str, kernel_obj: Nc.XcorKernel
    ) -> list[KernelVariants]:
        """Evaluate kernel(s) at specified multipole range.

        :param kernel_label: Label for the kernel.
        :param kernel_obj: NumCosmo XcorKernel object.
        :return: List of KernelVariants containing evaluation results.
        """
        lmax = self.ell + self.n_ell - 1
        if self.n_ell > 1:
            print(f"Evaluating kernels at ell = {self.ell} to {lmax}...")
        else:
            print(f"Evaluating kernels at ell = {self.ell}...")

        RH_Mpc = self.cosmo.RH_Mpc()
        kernel_obj.set_l_limber(self.l_limber)
        eval_kernel = kernel_obj.get_eval_vectorized(
            self.cosmo, self.ell, lmax, self.closure_type.to_nc()
        )

        primary_method = (
            "Non-Limber"
            if self.l_limber < 0
            else (
                "Kernel-Limber"
                if self.l_limber == 0
                else f"Limber(ell>={self.l_limber})"
            )
        )
        kernel_eval = KernelEvaluation(
            name=kernel_label,
            method=primary_method,
            kernel=kernel_obj,
            evaluator=eval_kernel,
            RH_Mpc=RH_Mpc,
            lmin=self.ell,
        )

        k_min, k_max = kernel_eval.range()

        print(f"  Non-Limber kernel k range: [{k_min:.2e}, {k_max:.2e}] Mpc^-1")

        kernel_eval_alt: KernelEvaluation | None = None
        if self._comparing:
            _, alt_label = self._curve_labels

            if self.compare_limber:
                # Same representation, the Limber approximation instead.
                kernel_obj.set_l_limber(0)
                alt_closure = self.closure_type
            else:
                # Same evaluation tier, the other representation of it.
                alt_closure = self._alt_closure_type

            eval_alt = kernel_obj.get_eval_vectorized(
                self.cosmo, self.ell, lmax, alt_closure.to_nc()
            )
            kernel_eval_alt = KernelEvaluation(
                name=kernel_label,
                kernel=kernel_obj,
                method=alt_label,
                evaluator=eval_alt,
                RH_Mpc=RH_Mpc,
                lmin=self.ell,
            )
            k_min_alt, k_max_alt = kernel_eval_alt.range()

            print(
                f"  {alt_label} kernel k range: [{k_min_alt:.2e}, "
                f"{k_max_alt:.2e}] Mpc^-1"
            )

        print("  [OK] Kernel evaluation complete")
        print()

        return [KernelVariants(main=kernel_eval, alternative=kernel_eval_alt)]

    def _plot_results(self, kernel_vars: list[KernelVariants]) -> None:
        """Plot kernel evaluation results.

        :param kernel_vars: List of KernelVariants to plot.
        """
        print("Plotting results...")

        # Define color palette for different kernels
        colors = plt.cm.tab10.colors  # type: ignore # pylint: disable=no-member
        ax1: plt.Axes
        ax2: plt.Axes

        if self._comparing:
            # Create figure with two subplots for comparison
            fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 12), sharex=True)
            if self.n_ell > 1:
                ell_range = f"$\\ell = {self.ell}$ to ${self.ell + self.n_ell - 1}$"
            else:
                ell_range = f"$\\ell = {self.ell}$"
            main_label, alt_label = self._curve_labels
            fig.suptitle(
                f"Kernel Comparison at {ell_range}\n"
                f"(solid = {main_label}; dashed/thin = {alt_label}, same color)",
                fontsize=14,
                fontweight="bold",
            )

            legend_handles = []
            legend_labels = []

            for idx, kernel_var in enumerate(kernel_vars):
                main_kernel = kernel_var.main
                alt_kernel = kernel_var.alternative

                assert alt_kernel is not None

                color = colors[idx % len(colors)]

                # Plot the main curve with full solid lines
                label_main = main_kernel.plot(
                    ax1,
                    n_points=self.n_points,
                    color=color,
                    linestyle="-",
                    linewidth=2.0,
                    alpha=1.0,
                    k_range=self.k_range,
                )
                # Plot the alternative with thinner dashed lines, reduced opacity
                alt_kernel.plot(
                    ax1,
                    n_points=self.n_points,
                    color=color,
                    linestyle="--",
                    linewidth=1.5,
                    alpha=0.6,
                    k_range=self.k_range,
                )

                # Plot the relative difference between the two
                main_kernel.plot_comparison(
                    alt_kernel, ax2, color=color, linestyle="-", k_range=self.k_range
                )

                # Add legend entries
                legend_handles.append(
                    Line2D([0], [0], color=color, linewidth=2, linestyle="-")
                )
                legend_labels.append(label_main)

            # Configure first subplot
            ax1.set_xlabel("$k$ [Mpc$^{-1}$]", fontsize=12)
            ax1.set_ylabel("$k^{3/2} |K(k)|$", fontsize=12)
            ax1.set_title("Kernel vs Wave Number", fontsize=11)
            ax1.legend(legend_handles, legend_labels, loc="best", fontsize=10)
            ax1.grid(True, alpha=0.3)

            # Configure second subplot
            ax2.axhline(y=0.0, color="k", linestyle="-", alpha=0.5, linewidth=0.8)
            ax2.set_xlabel("$k$ [Mpc$^{-1}$]", fontsize=12)
            ax2.set_ylabel(
                f"Relative Difference ({alt_label}/{main_label} - 1)", fontsize=12
            )
            ax2.set_title(f"{alt_label} vs {main_label} Comparison", fontsize=11)
            ax2.grid(True, alpha=0.3)
        else:
            # Create single plot
            fig, ax = plt.subplots(figsize=(10, 6))
            if self.n_ell > 1:
                ell_range = f"$\\ell = {self.ell}$ to ${self.ell + self.n_ell - 1}$"
            else:
                ell_range = f"$\\ell = {self.ell}$"
            fig.suptitle(
                f"Kernel at {ell_range}",
                fontsize=14,
                fontweight="bold",
            )

            legend_entries = []

            for idx, kernel_var in enumerate(kernel_vars):
                main_kernel = kernel_var.main
                color = colors[idx % len(colors)]

                label = main_kernel.plot(
                    ax,
                    n_points=self.n_points,
                    color=color,
                    linestyle="-",
                    k_range=self.k_range,
                )
                legend_entries.append(label)

            ax.set_xlabel("$k$ [Mpc$^{-1}$]", fontsize=12)
            ax.set_ylabel("$k^{3/2} |K(k)|$", fontsize=12)
            ax.set_title("Kernel vs Wave Number", fontsize=11)
            ax.legend(legend_entries, loc="best", fontsize=10)
            ax.grid(True, alpha=0.3)

        plt.tight_layout()

        # Save figure if output path specified
        if self.output:
            plt.savefig(self.output, dpi=150, bbox_inches="tight")
            print(f"  [OK] Figure saved to: {self.output}")

        # Show plot if requested
        if self.show_plot:
            print("  [OK] Displaying plot...")
            plt.show()
        else:
            print("  [OK] Plot generated (not displayed)")

        print()

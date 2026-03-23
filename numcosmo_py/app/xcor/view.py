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
from matplotlib.lines import Line2D
import numpy as np
import typer

from numcosmo_py import Nc, Ncm
from numcosmo_py.cosmology import Cosmology

from .kernels import (
    parse_kernel_spec,
    KERNEL_CONFIG_REGISTRY,
    KernelCMBLensingConfig,
    KernelCMBISWConfig,
    KernelTSZConfig,
    KernelGalaxyLSSTConfig,
    KernelWeakLensingLSSTConfig,
    KernelConfigTypes,
)

Ncm.cfg_init()


@dataclasses.dataclass(kw_only=True)
class KernelEvaluation:
    """Result of evaluating a kernel on a k-grid.

    :ivar name: Kernel name.
    :ivar method: Evaluation method (Limber or Non-Limber).
    :ivar kernel: NumCosmo XcorKernel object.
    :ivar evaluator: NumCosmo XcorKernelIntegrand evaluator.
    :ivar RH_Mpc: Hubble radius in Mpc.
    """

    name: str
    method: str
    kernel: Nc.XcorKernel
    evaluator: Nc.XcorKernelIntegrand
    RH_Mpc: float

    def range(self) -> tuple[float, float]:
        """Return k range in Mpc^-1.

        :return: Tuple of (k_min, k_max) in Mpc^-1.
        """
        kmin, kmax = self.evaluator.get_range()
        return kmin / self.RH_Mpc, kmax / self.RH_Mpc

    def evaluate(self, k_Mpc: np.ndarray) -> np.ndarray:
        """Evaluate kernel on a k grid.

        :param k_Mpc: Wave number array in Mpc^-1.
        :return: Kernel values evaluated at k_Mpc.
        """
        k = k_Mpc * self.RH_Mpc
        return np.array([self.evaluator.eval_array(ki)[0] for ki in k])

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
    ) -> None:
        """Plot kernel on given axes.

        :param ax: Matplotlib axes object.
        :param k_range: Optional k range (k_min, k_max) in Mpc^-1 for plotting.
        :param color: Line color.
        :param linestyle: Line style.
        """
        kmin, kmax = self.range()
        if k_range is not None:
            kmin = max(kmin, k_range[0])
            kmax = min(kmax, k_range[1])
        k_Mpc = np.logspace(np.log10(kmin), np.log10(kmax), n_points)
        kernel_values = self.evaluate(k_Mpc)
        k3_2 = k_Mpc ** (3 / 2)
        ax.plot(k_Mpc, k3_2 * kernel_values, color=color, linestyle=linestyle)
        ax.set_xscale("log")
        ax.set_yscale("symlog", linthresh=1e-10)

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
        ratio = other_values / kernel_values - 1.0

        ax.loglog(k_Mpc, np.abs(ratio), color=color, linestyle=linestyle)


@dataclasses.dataclass(kw_only=True)
class KernelVariants:
    """Container for kernel evaluation results.

    :ivar main: Main kernel evaluation result.
    :ivar alternative: Alternative kernel evaluation result (e.g., Limber
        approximation).
    """

    main: KernelEvaluation
    alternative: KernelEvaluation | None = None


@dataclasses.dataclass(kw_only=True)
class ViewKernel:
    """View cross-correlation kernels.

    This command visualizes cross-correlation kernels used in cosmological
    analyses. You can view single kernels or compare Limber vs non-Limber
    approximations.

    The kernel specification follows the format:

        "<kernel_name> key=value key=value ..."

    Available kernel types:
        - cmb_lensing: CMB lensing convergence kernel
        - cmb_isw: CMB Integrated Sachs-Wolfe kernel
        - tsz: Thermal Sunyaev-Zeldovich kernel
        - galaxy_lsst: Galaxy clustering kernel (LSST SRD bins)
        - weak_lensing_lsst: Weak lensing kernel (LSST SRD bins)

    Examples:
        # View CMB lensing kernel
        numcosmo xcor kernel view --kernel "cmb_lensing lmax=3000" --ell 20

        # View Galaxy clustering kernel (Y1 lens bin 0)
        numcosmo xcor kernel view \\
            --kernel "galaxy_lsst survey=y1 bin_idx=0 bias=1.5" \
            --ell 100

        # Compare Limber vs non-Limber for weak lensing
        numcosmo xcor kernel view \
            --kernel "weak_lensing_lsst survey=y10 bin_idx=2" \
            --ell 50 \\
            --compare-limber \\
            --output wl_comparison.png
    """

    kernel: Annotated[
        list[str],
        typer.Option(
            default_factory=lambda: ["cmb_lensing lmax=3000"],
            help=(
                "Kernel specification string. "
                "Format: '<kernel_name> key=value ...'. "
                f"Available types: {', '.join(KERNEL_CONFIG_REGISTRY.keys())}. "
                "Use --help to see kernel-specific parameters."
            ),
            show_default=True,
        ),
    ]

    ell: Annotated[
        int,
        typer.Option(
            help="Multipole value for kernel evaluation.",
            show_default=True,
            min=2,
        ),
    ] = 20

    k_range: Annotated[
        tuple[float, float] | None,
        typer.Option(
            help="k value range for plotting (Mpc^-1).",
            show_default=True,
        ),
    ] = None

    compare_limber: Annotated[
        bool,
        typer.Option(
            help="Compare Limber vs non-Limber approximations.",
            show_default=True,
        ),
    ] = False

    n_points: Annotated[
        int,
        typer.Option(
            help="Number of k points to evaluate.",
            show_default=True,
            min=100,
            max=50000,
        ),
    ] = 5000

    output: Annotated[
        Optional[Path],
        typer.Option(
            help="Output file path for plot (e.g., kernel_plot.png).",
        ),
    ] = None

    show_plot: Annotated[
        bool,
        typer.Option(
            help="Display plot interactively.",
            show_default=True,
        ),
    ] = True

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

        # Create cosmology with appropriate maximum redshift
        # Use larger dist_max_z for non-Limber calculations
        print("Creating cosmology...")
        dist_max_z = 1000.0 if self.compare_limber else 10.0
        nc_cosmo = Cosmology.default(dist_max_z=dist_max_z)
        self.cosmo = nc_cosmo.cosmo
        self.dist = nc_cosmo.dist
        self.ps_ml = nc_cosmo.ps_ml
        self.recomb = nc_cosmo.recomb
        print(f"  ✓ H0 = {self.cosmo['H0']:.2f}, Omega_b = {self.cosmo['Omegab']:.4f}")
        print(f"  ✓ Maximum redshift: {dist_max_z}")
        print()

        # Create integrator (matching fixtures pattern)
        print("Creating integrator...")
        self.integrator = Ncm.SBesselIntegratorLevin.new(0, 8)
        print("  ✓ Levin integrator created")
        print()

        print("Parsing kernel specification...")
        kernel_evals = []
        for k0 in self.kernel:
            kernel_name, kernel_config = parse_kernel_spec(k0)
            print(f"  ✓ Kernel type: {kernel_name}")
            print(f"  ✓ Configuration: {kernel_config}")
            print()

            # Create kernel(s)
            kernel_label, kernel_obj = self._create_kernels(kernel_config)

            # Evaluate kernels
            kernel_evals += self._evaluate_kernels(kernel_label, kernel_obj)

        # Plot results
        self._plot_results(kernel_evals)

        print()
        print("✓ Kernel visualization complete!")

    def _create_kernels(
        self, kernel_config: KernelConfigTypes
    ) -> tuple[str, Nc.XcorKernel]:
        """Create kernel objects based on configuration.

        :param kernel_config: Kernel configuration object.
        :return: Tuple of (kernel_label, kernel_object).
        """
        print("Creating kernel(s)...")

        match kernel_config:
            case KernelCMBLensingConfig():
                return self._create_cmb_lensing_kernels(kernel_config)
            case KernelCMBISWConfig():
                return self._create_cmb_isw_kernels(kernel_config)
            case KernelTSZConfig():
                return self._create_tsz_kernels(kernel_config)
            case KernelGalaxyLSSTConfig():
                return self._create_galaxy_lsst_kernels(kernel_config)
            case KernelWeakLensingLSSTConfig():
                return self._create_weak_lensing_lsst_kernels(kernel_config)
            case _:
                raise ValueError(f"Unknown kernel type: {type(kernel_config)}")

        print("  ✓ Kernels created and prepared")
        print()

    def _create_cmb_lensing_kernels(
        self, config: KernelCMBLensingConfig
    ) -> tuple[str, Nc.XcorKernelCMBLensing]:
        """Create CMB lensing kernel.

        :param config: CMB lensing configuration.
        :return: Tuple of (kernel_label, kernel_object).
        """
        assert isinstance(config, KernelCMBLensingConfig)

        lmax = config.lmax
        Nl = Ncm.Vector.new_array(np.arange(lmax + 1))
        Nl.set_zero()

        # Create primary kernel (non-Limber if compare_limber, Limber otherwise)
        kernel_obj = Nc.XcorKernelCMBLensing(
            dist=self.dist,
            powspec=self.ps_ml,
            recomb=self.recomb,
            Nl=Nl,
            lmax=lmax,
            integrator=self.integrator,
        )
        kernel_obj.set_lmax(lmax)
        kernel_obj.prepare(self.cosmo)

        kernel_label = "CMB Lensing"

        return kernel_label, kernel_obj

    def _create_cmb_isw_kernels(
        self, config: KernelCMBISWConfig
    ) -> tuple[str, Nc.XcorKernelCMBISW]:
        """Create CMB ISW kernel.

        :param config: CMB ISW configuration.
        :return: Tuple of (kernel_label, kernel_object).
        """
        assert isinstance(config, KernelCMBISWConfig)

        lmax = config.lmax
        Nl = Ncm.Vector.new_array(np.arange(lmax + 1))
        Nl.set_zero()

        # Create primary kernel
        kernel_obj = Nc.XcorKernelCMBISW(
            dist=self.dist,
            powspec=self.ps_ml,
            recomb=self.recomb,
            Nl=Nl,
            lmax=lmax,
            integrator=self.integrator,
        )
        kernel_obj.set_lmax(lmax)
        kernel_obj.prepare(self.cosmo)

        kernel_label = "CMB ISW"

        return kernel_label, kernel_obj

    def _create_tsz_kernels(
        self, config: KernelTSZConfig
    ) -> tuple[str, Nc.XcorKerneltSZ]:
        """Create tSZ kernel.

        :param config: tSZ configuration.
        :return: Tuple of (kernel_label, kernel_object).
        """
        assert isinstance(config, KernelTSZConfig)

        # Create primary kernel
        kernel_obj = Nc.XcorKerneltSZ(
            dist=self.dist,
            powspec=self.ps_ml,
            zmax=config.zmax,
            integrator=self.integrator,
        )
        kernel_obj.prepare(self.cosmo)
        kernel_label = "tSZ"

        return kernel_label, kernel_obj

    def _create_galaxy_lsst_kernels(
        self, config: KernelGalaxyLSSTConfig
    ) -> tuple[str, Nc.XcorKernelGal]:
        """Create Galaxy LSST kernel.

        :param config: Galaxy LSST configuration.
        :return: Tuple of (kernel_label, kernel_object).
        """
        assert isinstance(config, KernelGalaxyLSSTConfig)

        # Get LSST bins
        bins, _ = Nc.GalaxySDObsRedshiftGauss.new_lsst_srd_bins(config.bin_type.genum)
        dndz_spline = bins[config.bin_idx].compute_binned_dndz(None)

        # Create primary kernel
        kernel_obj = Nc.XcorKernelGal(
            dist=self.dist,
            powspec=self.ps_ml,
            dndz=dndz_spline,
            domagbias=config.domagbias,
            integrator=self.integrator,
        )
        kernel_obj.orig_vparam_set(Nc.XcorKernelGalVParams.BIAS, 0, config.bias)
        kernel_obj.orig_param_set(Nc.XcorKernelGalSParams.MAG_BIAS, config.mag_bias)
        kernel_obj.prepare(self.cosmo)

        kernel_label = (
            f"Galaxy LSST ({config.survey.value.upper()} bin {config.bin_idx})"
        )

        return kernel_label, kernel_obj

    def _create_weak_lensing_lsst_kernels(
        self, config: KernelWeakLensingLSSTConfig
    ) -> tuple[str, Nc.XcorKernelWeakLensing]:
        """Create Weak Lensing LSST kernel.

        :param config: Weak lensing LSST configuration.
        :return: Tuple of (kernel_label, kernel_object).
        """
        assert isinstance(config, KernelWeakLensingLSSTConfig)

        # Get LSST bins
        bins, _ = Nc.GalaxySDObsRedshiftGauss.new_lsst_srd_bins(config.bin_type.genum)
        dndz_spline = bins[config.bin_idx].compute_binned_dndz(None)

        # Create primary kernel
        kernel_obj = Nc.XcorKernelWeakLensing(
            dist=self.dist,
            powspec=self.ps_ml,
            dndz=dndz_spline,
            nbar=config.nbar,
            intr_shear=config.intr_shear,
            integrator=self.integrator,
        )
        kernel_obj.prepare(self.cosmo)

        kernel_label = (
            f"Weak Lensing LSST ({config.survey.value.upper()} bin {config.bin_idx})"
        )

        return kernel_label, kernel_obj

    def _evaluate_kernels(
        self, kernel_label: str, kernel_obj: Nc.XcorKernel
    ) -> list[KernelVariants]:
        """Evaluate kernel(s) at specified multipole.

        :param kernel_label: Label for the kernel.
        :param kernel_obj: NumCosmo XcorKernel object.
        :return: List of KernelVariants containing evaluation results.
        """
        print(f"Evaluating kernels at ell = {self.ell}...")

        RH_Mpc = self.cosmo.RH_Mpc()
        # Get evaluator for primary kernel and its range
        kernel_obj.set_l_limber(-1)
        eval_kernel = kernel_obj.get_eval_vectorized(self.cosmo, self.ell, self.ell)

        kernel_eval = KernelEvaluation(
            name=kernel_label,
            method="Limber" if not self.compare_limber else "Non-Limber",
            kernel=kernel_obj,
            evaluator=eval_kernel,
            RH_Mpc=RH_Mpc,
        )

        k_min, k_max = kernel_eval.range()

        print(f"  Primary kernel k range: [{k_min:.2e}, {k_max:.2e}] Mpc^-1")

        kernel_eval_limber: KernelEvaluation | None = None
        if self.compare_limber:
            kernel_obj.set_l_limber(0)
            eval_limber = kernel_obj.get_eval_vectorized(self.cosmo, self.ell, self.ell)
            kernel_eval_limber = KernelEvaluation(
                name=kernel_label,
                kernel=kernel_obj,
                method="Limber",
                evaluator=eval_limber,
                RH_Mpc=RH_Mpc,
            )
            k_min_limber, k_max_limber = kernel_eval_limber.range()

            print(
                f"  Limber kernel k range: [{k_min_limber:.2e}, "
                f"{k_max_limber:.2e}] Mpc^-1"
            )

        print("  ✓ Kernel evaluation complete")
        print()

        return [KernelVariants(main=kernel_eval, alternative=kernel_eval_limber)]

    def _plot_results(self, kernel_vars: list[KernelVariants]) -> None:
        """Plot kernel evaluation results.

        :param kernel_vars: List of KernelVariants to plot.
        """
        print("Plotting results...")

        # Define color palette for different kernels
        colors = plt.cm.tab10.colors  # type: ignore # pylint: disable=no-member
        ax1: plt.Axes
        ax2: plt.Axes

        if self.compare_limber:
            # Create figure with two subplots for comparison
            fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 12), sharex=True)
            fig.suptitle(
                f"Kernel Comparison at $\\ell = {self.ell}$ "
                "(dashed = Non-Limber, solid = Limber)",
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

                # Plot both kernels with same color, different line styles
                main_kernel.plot(
                    ax1,
                    n_points=self.n_points,
                    color=color,
                    linestyle="--",
                    k_range=self.k_range,
                )
                alt_kernel.plot(
                    ax1,
                    n_points=self.n_points,
                    color=color,
                    linestyle="-",
                    k_range=self.k_range,
                )

                # Plot comparison
                main_kernel.plot_comparison(
                    alt_kernel, ax2, color=color, linestyle="-", k_range=self.k_range
                )

                # Add single legend entry per kernel (with both line styles shown)
                legend_handles.append(
                    Line2D([0], [0], color=color, linewidth=2, linestyle="-")
                )
                legend_labels.append(main_kernel.name)

            # Configure first subplot
            ax1.set_xlabel("$k$ [Mpc$^{-1}$]", fontsize=12)
            ax1.set_ylabel("$k^{3/2} |K(k)|$", fontsize=12)
            ax1.set_title("Kernel vs Wave Number", fontsize=11)
            ax1.legend(legend_handles, legend_labels, fontsize=10)
            ax1.grid(True, alpha=0.3)

            # Configure second subplot
            ax2.axhline(y=1.0, color="k", linestyle="--", alpha=0.5)
            ax2.set_xlabel("$k$ [Mpc$^{-1}$]", fontsize=12)
            ax2.set_ylabel("Relative Difference", fontsize=12)
            ax2.set_title("Limber vs Non-Limber Comparison", fontsize=11)
            ax2.grid(True, alpha=0.3)
        else:
            # Create single plot
            fig, ax = plt.subplots(figsize=(10, 6))
            fig.suptitle(
                f"Kernel at $\\ell = {self.ell}$",
                fontsize=14,
                fontweight="bold",
            )

            legend_entries = []
            for idx, kernel_var in enumerate(kernel_vars):
                main_kernel = kernel_var.main
                color = colors[idx % len(colors)]

                main_kernel.plot(
                    ax,
                    n_points=self.n_points,
                    color=color,
                    linestyle="-",
                    k_range=self.k_range,
                )
                legend_entries.append(f"{main_kernel.name} ({main_kernel.method})")

            ax.set_xlabel("$k$ [Mpc$^{-1}$]", fontsize=12)
            ax.set_ylabel("$k^{3/2} |K(k)|$", fontsize=12)
            ax.set_title("Kernel vs Wave Number", fontsize=11)
            ax.legend(legend_entries, fontsize=10)
            ax.grid(True, alpha=0.3)

        plt.tight_layout()

        # Save figure if output path specified
        if self.output:
            plt.savefig(self.output, dpi=150, bbox_inches="tight")
            print(f"  ✓ Figure saved to: {self.output}")

        # Show plot if requested
        if self.show_plot:
            print("  ✓ Displaying plot...")
            plt.show()
        else:
            print("  ✓ Plot generated (not displayed)")

        print()

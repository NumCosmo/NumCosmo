#
# common.py
#
# Thu Sep 10 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# common.py
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

r"""Common block for the cross-correlation kernel commands.

Everything the kernel view and the angular power spectrum commands share --
the kernel specification, the cosmology, the integrator, the comparison mode
and the $C_\\ell$ solve itself -- lives here, so the two commands stay in sync
through inheritance rather than through duplicated option blocks.
"""

import dataclasses
import enum
import time
from typing import Annotated, Any, Optional
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import typer

from numcosmo_py import Nc, Ncm
from numcosmo_py.cosmology import Cosmology

from .kernels import (
    CMBLensingSource,
    _KernelRadialConfig,
    parse_kernel_spec,
    LSSTBinType,
    KernelCMBLensingConfig,
    KernelCMBISWConfig,
    KernelTSZConfig,
    KernelNumberCountsConfig,
    KernelWeakLensingConfig,
    KernelClusterTophatConfig,
    KernelRadialGaussConfig,
    KernelRadialTophatConfig,
    KernelRadialTophatSmoothConfig,
    KernelRadialStudentTConfig,
    KernelRadialPowerExpConfig,
    KernelRadialLensingConfig,
    KernelRadialMultiConfig,
    KernelConfigTypes,
)

Ncm.cfg_init()


class XcorMethodOption(str, enum.Enum):
    """Quadrature methods available for the C_ell computation."""

    CUBATURE = "cubature"
    GSL = "gsl"
    EXACT = "exact"

    def to_nc(self) -> Nc.XcorMethod:
        """Convert to the corresponding #NcXcorMethod value.

        :return: The NumCosmo enumeration value.
        """
        match self:
            case XcorMethodOption.CUBATURE:
                return Nc.XcorMethod.KERNEL_CUBATURE
            case XcorMethodOption.GSL:
                return Nc.XcorMethod.KERNEL_GSL
            case XcorMethodOption.EXACT:
                return Nc.XcorMethod.KERNEL_EXACT
        raise ValueError(f"Unknown method: {self}")


class XcorClosureOption(str, enum.Enum):
    """Representations available for the k-space closure."""

    SPLINE = "spline"
    CHEBYSHEV = "chebyshev"

    def to_nc(self) -> Nc.XcorKernelClosure:
        """Convert to the corresponding #NcXcorKernelClosure value.

        :return: The NumCosmo enumeration value.
        """
        match self:
            case XcorClosureOption.SPLINE:
                return Nc.XcorKernelClosure.SPLINE
            case XcorClosureOption.CHEBYSHEV:
                return Nc.XcorKernelClosure.CHEBYSHEV
        raise ValueError(f"Unknown closure type: {self}")


def contiguous_runs(ells: np.ndarray) -> list[tuple[int, int]]:
    r"""Split a sorted multipole array into its maximal contiguous runs.

    The solver requests one $\\ell$-range at a time, so a sampled (and thus
    gappy) multipole set is asked for as one request per run rather than as one
    request per multipole: consecutive multipoles then still share a block.

    :param ells: Strictly increasing array of multipoles.
    :return: List of inclusive (lmin, lmax) pairs, in increasing order.
    """
    runs: list[tuple[int, int]] = []
    for ell in (int(ell) for ell in ells):
        if runs and ell == runs[-1][1] + 1:
            runs[-1] = (runs[-1][0], ell)
        else:
            runs.append((ell, ell))
    return runs


@dataclasses.dataclass(kw_only=True)
class XcorKernelCommon:
    r"""Common options for the cross-correlation kernel commands.

    Builds the cosmology, the Levin integrator and every requested kernel, and
    carries the $C_\\ell$ solve shared by :class:`ViewKernel` and
    :class:`ComputeCls`. All commands that work with cross-correlation kernels
    should inherit from this class.
    """

    kernel: Annotated[
        list[str],
        typer.Option(
            default_factory=lambda: ["cmb_lensing lmax=3000"],
            help=(
                "Kernel specification string. "
                "Format: '<kernel_name> key=value ...'. "
                "Use 'numcosmo xcor kernel list' to see the available kernel "
                "types, and 'numcosmo xcor kernel list <kernel_type>' for what "
                "each of that one's parameters means."
            ),
            show_default=True,
        ),
    ]

    l_limber: Annotated[
        int,
        typer.Option(
            help=(
                "Limber threshold for the primary evaluation "
                "(-1: never [true non-Limber], 0: always [kernel-Limber], "
                "N>0: Limber for ell>=N). See dev-notes/"
                "xcor_ultralevin_batching_plan.md for tier semantics."
            ),
            show_default=True,
        ),
    ] = -1

    closure_type: Annotated[
        XcorClosureOption,
        typer.Option(
            help=(
                "Representation fitted to the sampled kernel. 'spline' bisects "
                "until it meets a tolerance; 'chebyshev' expands on panels of a "
                "prescribed order. Both plot and both compute C_ell, so the two "
                "can be compared directly. Limber multipoles keep the spline "
                "whatever this is set to."
            ),
            show_default=True,
        ),
    ] = XcorClosureOption.CHEBYSHEV

    compare_limber: Annotated[
        bool,
        typer.Option(
            help=(
                "Also show Limber approximation for comparison "
                "(with thinner dashed lines)."
            ),
            show_default=True,
        ),
    ] = False

    compare_closure: Annotated[
        bool,
        typer.Option(
            help=(
                "Also show the other closure representation for comparison "
                "(with thinner dashed lines). Mutually exclusive with "
                "--compare-limber: there is one alternative curve. Under Limber "
                "both representations are the spline, so pair this with "
                "--l-limber -1 for it to show anything."
            ),
            show_default=True,
        ),
    ] = False

    integrator_reltol: Annotated[
        Optional[float],
        typer.Option(
            min=0.0,
            max=1.0,
            help=(
                "NcmSBesselIntegratorLevin ODE solve relative tolerance "
                "(library default: 1e-13, near machine precision). "
                "This is the dominant cost/precision knob for tier 3 -- "
                "see dev-notes/xcor_ultralevin_batching_plan.md sec 9.4. "
                "Leave unset to keep the library default."
            ),
        ),
    ] = None

    integrator_cheb_reltol: Annotated[
        Optional[float],
        typer.Option(
            min=0.0,
            max=1.0,
            help=(
                "NcmSBesselIntegratorLevin integrand Chebyshev-fit relative "
                "tolerance (library default: 1e-8). The looser of this and "
                "--integrator-reltol bounds the result. Leave unset to keep "
                "the library default."
            ),
        ),
    ] = None

    integrator_max_order: Annotated[
        Optional[int],
        typer.Option(
            help=(
                "NcmSBesselIntegratorLevin maximum spectral order "
                "(library default: 16384). Leave unset to keep the library "
                "default."
            ),
        ),
    ] = None

    cls_method: Annotated[
        XcorMethodOption,
        typer.Option(
            help=(
                "Quadrature used for the C_ell computation. 'cubature' and 'gsl' "
                "target a tolerance and abort if they cannot reach it; 'exact' "
                "integrates the closures exactly, on the common refinement of "
                "their knots, so it needs no tolerance and cannot fail to "
                "converge."
            ),
            show_default=True,
        ),
    ] = XcorMethodOption.EXACT

    cls_block_size: Annotated[
        int,
        typer.Option(
            min=1,
            help=(
                "Multipole block size handed to NcXcorSolver.plan_blocks(). "
                "Eight was the fastest block size measured; see "
                "dev-notes/xcor_ultralevin_batching_plan.md section 1.3."
            ),
            show_default=True,
        ),
    ] = 8

    cross: Annotated[
        bool,
        typer.Option(
            help=(
                "Include the cross spectra of every kernel pair. With --no-cross "
                "only the auto spectra are computed, which is n(n+1)/2 -> n "
                "solves for n kernels."
            ),
            show_default=True,
        ),
    ] = True

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
        """Build the cosmology, the integrator and every requested kernel.

        :raises ValueError: If the kernel specification is invalid.
        :raises typer.BadParameter: If two mutually exclusive options are given.
        """
        # One alternative curve, so one thing to compare against.
        if self.compare_limber and self.compare_closure:
            raise ValueError(
                "--compare-limber and --compare-closure both draw the "
                "alternative curve; pick one."
            )

        # typer's min is inclusive, but zero tolerance aborts in the library.
        for name, tol in (
            ("--integrator-reltol", self.integrator_reltol),
            ("--integrator-cheb-reltol", self.integrator_cheb_reltol),
        ):
            if tol is not None and tol <= 0.0:
                raise typer.BadParameter(f"{name} must be positive, got {tol}.")

        self._create_cosmology()
        self._create_integrator()

        print("Parsing kernel specification...")
        self._solver_cache: Optional[tuple[Nc.XcorSolver, list[int]]] = None
        self.kernels: list[tuple[str, Nc.XcorKernel]] = []
        for spec in self.kernel:
            kernel_name, kernel_config = parse_kernel_spec(spec)
            print(f"  [OK] Kernel type: {kernel_name}")
            print(f"  [OK] Configuration: {kernel_config}")
            print()

            self.kernels.append(self._create_kernels(kernel_config))

        self._disambiguate_labels()

    def _disambiguate_labels(self) -> None:
        """Add an index to kernels whose labels came out the same.

        A label names the shape, not the parameters, so two bins of one shape
        are indistinguishable in a legend -- which is precisely the case a
        cross-spectrum is asked for in.
        """
        counts: dict[str, int] = {}
        for label, _ in self.kernels:
            counts[label] = counts.get(label, 0) + 1

        seen: dict[str, int] = {}
        for idx, (label, kernel_obj) in enumerate(self.kernels):
            if counts[label] == 1:
                continue
            seen[label] = seen.get(label, 0) + 1
            # Parentheses, not '#': a label is drawn as LaTeX when the user has
            # usetex on, and '#' is a macro parameter character there.
            self.kernels[idx] = (f"{label} ({seen[label]})", kernel_obj)

    def _create_cosmology(self) -> None:
        """Create the cosmology every kernel is built against."""
        # Always use larger dist_max_z for non-Limber calculations (default behavior)
        print("Creating cosmology...")
        dist_max_z = 1000.0
        nc_cosmo = Cosmology.default(dist_max_z=dist_max_z)
        self.cosmo = nc_cosmo.cosmo
        self.dist = nc_cosmo.dist
        self.ps_ml = nc_cosmo.ps_ml
        self.recomb = nc_cosmo.recomb
        print(
            f"  [OK] H0 = {self.cosmo['H0']:.2f}, Omega_b = {self.cosmo['Omegab']:.4f}"
        )
        print(f"  [OK] Maximum redshift: {dist_max_z}")
        print()

    def _create_integrator(self) -> None:
        """Create the Levin integrator shared by every kernel.

        Pass the tolerances at construction: set_reltol() would rebuild every
        operator to apply them, which is pure waste when nothing has run yet.
        """
        print("Creating integrator...")
        defaults = Ncm.SBesselIntegratorLevin.new(0, 8)
        self.integrator = Ncm.SBesselIntegratorLevin.new_full(
            0,
            8,
            defaults.get_x_knots_min(),
            defaults.get_x_knots_max(),
            defaults.get_n_knots(),
            defaults.get_ell_cache_max(),
            (
                self.integrator_reltol
                if self.integrator_reltol is not None
                else defaults.get_reltol()
            ),
            defaults.get_cheb_min_order(),
            (
                self.integrator_cheb_reltol
                if self.integrator_cheb_reltol is not None
                else defaults.get_cheb_reltol()
            ),
        )
        if self.integrator_max_order is not None:
            self.integrator.set_max_order(self.integrator_max_order)
        # A closure cannot be fitted to more precision than the integrator
        # samples, and the library refuses the pairing. The integrator
        # tolerances are the ones this command exposes, so a request to compute
        # loosely is honoured by loosening the fit to match rather than by
        # failing.
        self.closure_tol_floor = max(
            self.integrator.get_reltol(), self.integrator.get_cheb_reltol()
        )
        print(
            f"  [OK] Levin integrator created "
            f"(reltol={self.integrator.get_reltol():.1e}, "
            f"cheb_reltol={self.integrator.get_cheb_reltol():.1e}, "
            f"max_order={self.integrator.get_max_order()})"
        )
        print()

    @property
    def _alt_closure_type(self) -> XcorClosureOption:
        """The representation the comparison curve uses.

        :return: The option other than the one --closure-type selected.
        """
        if self.closure_type is XcorClosureOption.SPLINE:
            return XcorClosureOption.CHEBYSHEV
        return XcorClosureOption.SPLINE

    @property
    def _comparing(self) -> bool:
        """Whether a second, alternative curve is drawn beside the main one.

        :return: True when either comparison mode is on.
        """
        return self.compare_limber or self.compare_closure

    @property
    def _curve_labels(self) -> tuple[str, str]:
        """Names for the main and the alternative curve, for titles and axes.

        :return: Tuple of (main label, alternative label).
        """
        if self.compare_closure:
            return self.closure_type.value.capitalize(), (
                self._alt_closure_type.value.capitalize()
            )
        main = (
            "Non-Limber"
            if self.l_limber < 0
            else ("Kernel-Limber" if self.l_limber == 0 else "Limber")
        )
        return main, "Limber"

    def _create_kernels(
        self, kernel_config: KernelConfigTypes
    ) -> tuple[str, Nc.XcorKernel]:
        """Create kernel objects based on configuration.

        :param kernel_config: Kernel configuration object.
        :return: Tuple of (kernel_label, kernel_object).
        """
        print("Creating kernel(s)...")

        # Each branch returns its own kernel subclass; widen to the common base so
        # the first branch does not fix the type of the others.
        result: tuple[str, Nc.XcorKernel]

        match kernel_config:
            case KernelCMBLensingConfig():
                result = self._create_cmb_lensing_kernels(kernel_config)
            case KernelCMBISWConfig():
                result = self._create_cmb_isw_kernels(kernel_config)
            case KernelTSZConfig():
                result = self._create_tsz_kernels(kernel_config)
            case KernelNumberCountsConfig():
                result = self._create_number_counts_kernels(kernel_config)
            case KernelWeakLensingConfig():
                result = self._create_weak_lensing_kernels(kernel_config)
            case KernelClusterTophatConfig():
                result = self._create_cluster_tophat_kernels(kernel_config)
            case (
                KernelRadialGaussConfig()
                | KernelRadialTophatConfig()
                | KernelRadialTophatSmoothConfig()
                | KernelRadialStudentTConfig()
                | KernelRadialPowerExpConfig()
                | KernelRadialLensingConfig()
                | KernelRadialMultiConfig()
            ):
                result = self._create_radial_kernels(kernel_config)
            case _:
                raise ValueError(f"Unknown kernel type: {type(kernel_config)}")

        self._apply_closure_tolerance_floor(result[1])

        print("  [OK] Kernels created and prepared")
        print()

        return result

    def _apply_closure_tolerance_floor(self, kernel: Nc.XcorKernel) -> None:
        """Keep the closure fit no tighter than the integrator samples.

        The library refuses a kernel that fits its k-space closure tighter than
        its integrator carries -- below that the sampled window is not a smooth
        function and no representation converges on it. This command exposes the
        integrator tolerances and not the kernel's, so a request to compute
        loosely is honoured by loosening the fit to match.

        :param kernel: Kernel whose fit tolerances may need loosening.
        """
        floor = self.closure_tol_floor

        if kernel.get_reltol() >= floor and kernel.get_peak_epsilon() >= floor:
            return

        kernel.set_reltol(max(kernel.get_reltol(), floor))
        kernel.set_peak_epsilon(max(kernel.get_peak_epsilon(), floor))
        print(
            f"  [OK] Closure fit tolerances raised to {floor:.1e} to match the "
            f"integrator"
        )

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
            source=config.source.genum,
            integrator=self.integrator,
        )
        kernel_obj.set_lmax(lmax)
        kernel_obj.prepare(self.cosmo)

        kernel_label = (
            config.label
            if config.source is CMBLensingSource.THIN_SCREEN
            else f"{config.label} ({config.source.value})"
        )

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

        kernel_label = config.label

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
        kernel_label = config.label

        return kernel_label, kernel_obj

    def _lsst_srd_bin_dndz(
        self, bin_type: LSSTBinType, bin_idx: int, survey: str
    ) -> Ncm.Spline:
        """Compute the dN/dz spline for one LSST-SRD photo-z bin.

        :param bin_type: LSST year/sample bin type.
        :param bin_idx: Index of the bin within that type's edges.
        :param survey: Survey label, used only for the error message.
        :return: The binned dN/dz spline.
        """
        edges, population, observable_population = (
            Nc.GalaxyRedshiftBinning.lsst_srd_edges(bin_type.genum)
        )
        n_bins = edges.len() - 1

        if bin_idx >= n_bins:
            raise ValueError(
                f"Bin index {bin_idx} is out of range for survey '{survey}'"
            )

        binning = Nc.GalaxyRedshiftBinning.new()

        return binning.compute_dndz(
            population,
            observable_population,
            edges.get(bin_idx),
            edges.get(bin_idx + 1),
        )

    def _create_number_counts_kernels(
        self, config: KernelNumberCountsConfig
    ) -> tuple[str, Nc.XcorKernelGal]:
        """Create number counts kernel.

        :param config: Number counts configuration.
        :return: Tuple of (kernel_label, kernel_object).
        """
        assert isinstance(config, KernelNumberCountsConfig)

        dndz_spline = self._lsst_srd_bin_dndz(
            config.bin_type, config.bin_idx, config.survey
        )

        # Create primary kernel
        kernel_obj = Nc.XcorKernelGal(
            dist=self.dist,
            powspec=self.ps_ml,
            dndz=dndz_spline,
            domagbias=config.domagbias,
            dorsd=config.dorsd,
            integrator=self.integrator,
        )
        kernel_obj.orig_vparam_set(Nc.XcorKernelGalVParams.BIAS, 0, config.bias)
        kernel_obj.orig_param_set(Nc.XcorKernelGalSParams.MAG_BIAS, config.mag_bias)
        kernel_obj.prepare(self.cosmo)

        kernel_label = f"{config.label} ({config.survey} bin {config.bin_idx})"

        return kernel_label, kernel_obj

    def _create_weak_lensing_kernels(
        self, config: KernelWeakLensingConfig
    ) -> tuple[str, Nc.XcorKernelWeakLensing]:
        """Create weak lensing kernel.

        :param config: Weak lensing configuration.
        :return: Tuple of (kernel_label, kernel_object).
        """
        assert isinstance(config, KernelWeakLensingConfig)

        dndz_spline = self._lsst_srd_bin_dndz(
            config.bin_type, config.bin_idx, config.survey
        )

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

        kernel_label = f"{config.label} ({config.survey} bin {config.bin_idx})"

        return kernel_label, kernel_obj

    def _create_radial_kernels(
        self, config: KernelConfigTypes
    ) -> tuple[str, Nc.XcorKernelRadial]:
        """Create one of the analytic radial windows.

        One builder for all seven shapes: they differ only in which closed form
        they carry, and every one takes its parameters straight through as
        construct properties. The vector-valued ones are wrapped on the way in.

        These are the shapes the Arb truth tables certify, so a curve drawn from
        one of them can be shown against proven values -- which no physical kernel
        above can offer, since none has a closed form.

        :param config: One of the radial window configurations.
        :return: Tuple of (kernel_label, kernel_object).
        """
        props: dict[str, Any] = {}
        for field, value in config.model_dump().items():
            props[field] = (
                Ncm.Vector.new_array(value) if isinstance(value, list) else value
            )

        kernel_obj = config.nc_type(
            dist=self.dist,
            powspec=self.ps_ml,
            integrator=self.integrator,
            **props,
        )
        kernel_obj.prepare(self.cosmo)

        # Only the radial shapes carry the Bessel-derivative order; the dispatch
        # table above already rejected everything else.
        assert isinstance(config, _KernelRadialConfig)
        deriv = config.bessel_deriv
        weight = "" if deriv == 0 else f", $j_\\ell^{{({deriv})}}$"
        kernel_label = f"{config.label}{weight}"

        return kernel_label, kernel_obj

    def _create_cluster_tophat_kernels(
        self, config: KernelClusterTophatConfig
    ) -> tuple[str, Nc.XcorKernelClusterTophat]:
        """Create cluster tophat kernel.

        :param config: Cluster tophat configuration.
        :return: Tuple of (kernel_label, kernel_object).
        """
        assert isinstance(config, KernelClusterTophatConfig)

        # Create primary kernel
        kernel_obj = Nc.XcorKernelClusterTophat(
            dist=self.dist,
            powspec=self.ps_ml,
            z_lower=config.z_lower,
            z_upper=config.z_upper,
            integrator=self.integrator,
        )
        kernel_obj.prepare(self.cosmo)

        kernel_label = (
            f"{config.label} (z=[{config.z_lower:.2f}, {config.z_upper:.2f}])"
        )

        return kernel_label, kernel_obj

    @property
    def _pairs(self) -> list[tuple[int, int]]:
        """The kernel index pairs whose spectra are computed.

        :return: Every auto- and cross-pair, or only the auto-pairs under
            --no-cross.
        """
        n_kernels = len(self.kernels)
        if not self.cross:
            return [(i, i) for i in range(n_kernels)]
        return [(i, j) for i in range(n_kernels) for j in range(i, n_kernels)]

    def _cls_output_path(self) -> Optional[Path]:
        """Where the C_ell figure is written.

        :return: The output path, or None when nothing is to be saved.
        """
        return self.output

    def compute_and_plot_cls(self, ells: np.ndarray) -> None:
        """Compute the requested spectra over @ells and plot them.

        :param ells: Strictly increasing array of multipoles.
        """
        cls_main = self._compute_cls(ells, self.l_limber, self.closure_type)
        cls_alt = None
        if self.compare_limber:
            cls_alt = self._compute_cls(ells, 0, self.closure_type)
        elif self.compare_closure:
            cls_alt = self._compute_cls(ells, self.l_limber, self._alt_closure_type)
        self._plot_cls(ells, cls_main, cls_alt)

    def _solver(self) -> tuple[Nc.XcorSolver, list[int]]:
        """Return the solver the C_ell runs share, registering kernels once.

        :return: Tuple of (solver, kernel ids in :attr:`kernels` order).
        """
        if self._solver_cache is None:
            solver = Nc.XcorSolver.new()
            ids = [solver.register_kernel(kernel) for _, kernel in self.kernels]
            self._solver_cache = (solver, ids)
        return self._solver_cache

    def _compute_cls(
        self,
        ells: np.ndarray,
        l_limber: int,
        closure_type: XcorClosureOption,
    ) -> dict[tuple[int, int], np.ndarray]:
        """Compute C_ell for every requested pair of the created kernels.

        Every kernel is put in the requested Limber mode first: the kernel view
        leaves the kernels in Limber mode when ``--compare-limber`` is used, so
        the mode must be set explicitly here rather than assumed.

        Multipoles are asked for one contiguous run at a time. A run is the unit
        the solver can batch; a sampled range is a sequence of short runs, and
        asking for each separately is what keeps the solve to the multipoles
        actually wanted.

        :param ells: Strictly increasing array of multipoles.
        :param l_limber: Limber threshold to apply to every kernel.
        :param closure_type: Representation to fit to every kernel.
        :return: Mapping from (i, j) kernel index pair to the C_ell array, one
            entry per multipole in @ells.
        """
        method_label = (
            "non-Limber"
            if l_limber < 0
            else ("Limber" if l_limber == 0 else f"Limber(ell>={l_limber})")
        )
        pairs = self._pairs
        runs = contiguous_runs(ells)

        print(
            f"Computing C_ell ({method_label}) for {len(ells)} multipole(s), "
            f"ell = {ells[0]} to {ells[-1]}..."
        )
        print(
            f"  {len(self.kernels)} kernel(s), {len(pairs)} spectra, "
            f"method={self.cls_method.value}, closure={closure_type.value}, "
            f"block size={self.cls_block_size}"
        )

        self.dist.prepare_if_needed(self.cosmo)
        self.ps_ml.prepare_if_needed(self.cosmo)

        for _, kernel_obj in self.kernels:
            kernel_obj.set_l_limber(l_limber)

        xcor = Nc.Xcor.new(self.dist, self.ps_ml, self.cls_method.to_nc())
        xcor.set_closure_type(closure_type.to_nc())

        # One solver for the whole command: it keeps a factorised integrator per
        # block, and plan_blocks() hands those back whenever the plan comes out
        # the same. The comparison run tiles the same multipoles as the main one,
        # so it starts warm.
        solver, ids = self._solver()
        solver.clear_requests()
        for i, j in pairs:
            for lmin, lmax in runs:
                solver.request_cl(ids[i], ids[j], lmin, lmax)
        solver.plan_blocks(self.cls_block_size)

        start = time.monotonic()
        solver.solve(xcor, self.cosmo)
        elapsed = time.monotonic() - start

        # One request per (pair, run), in that order: a pair's runs concatenate
        # back into one array aligned with @ells.
        result = {}
        for pair_idx, (i, j) in enumerate(pairs):
            offset = pair_idx * len(runs)
            result[(i, j)] = np.concatenate(
                [
                    np.array(solver.get_result(offset + run_idx).dup_array())
                    for run_idx in range(len(runs))
                ]
            )

        print(
            f"  [OK] {len(pairs)} spectra in {elapsed:.2f} s "
            f"({solver.get_n_blocks()} multipole block(s))"
        )
        print()

        return result

    def _plot_cls(
        self,
        ells: np.ndarray,
        cls_main: dict[tuple[int, int], np.ndarray],
        cls_alt: dict[tuple[int, int], np.ndarray] | None,
    ) -> None:
        """Plot the angular power spectra, optionally against the alternative.

        :param ells: Multipoles the spectra were computed at.
        :param cls_main: C_ell computed with the primary mode and representation.
        :param cls_alt: C_ell from the comparison run, or None.
        """
        print("Plotting C_ell...")

        colors = plt.cm.tab10.colors  # type: ignore # pylint: disable=no-member
        ax1: plt.Axes

        if cls_alt is not None:
            fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10), sharex=True)
        else:
            fig, ax1 = plt.subplots(1, 1, figsize=(10, 6))
            ax2 = None

        for idx, ((i, j), cl) in enumerate(cls_main.items()):
            color = colors[idx % len(colors)]
            label = (
                self.kernels[i][0]
                if i == j
                else f"{self.kernels[i][0]} x {self.kernels[j][0]}"
            )
            ax1.plot(ells, np.abs(cl), color=color, label=label)
            if cls_alt is not None and ax2 is not None:
                with np.errstate(divide="ignore", invalid="ignore"):
                    ratio = np.where(cl != 0.0, cls_alt[(i, j)] / cl - 1.0, np.nan)
                ax2.plot(ells, ratio, color=color, label=label)

        ax1.set_ylabel(r"$|C_\ell|$")
        ax1.set_yscale("log")
        if len(ells) > 1:
            ax1.set_xscale("log")
        ax1.grid(True, alpha=0.3)
        ax1.legend(fontsize=8)
        ax1.set_title("Angular power spectra", fontweight="bold")

        if ax2 is not None:
            ax2.axhline(0.0, color="black", lw=0.8)
            main_label, alt_label = self._curve_labels
            ax2.set_ylabel(
                rf"$C_\ell^{{\rm {alt_label}}}/C_\ell^{{\rm {main_label}}} - 1$"
            )
            ax2.set_xlabel(r"$\ell$")
            if len(ells) > 1:
                ax2.set_xscale("log")
            ax2.grid(True, alpha=0.3)
        else:
            ax1.set_xlabel(r"$\ell$")

        fig.tight_layout()

        cls_output = self._cls_output_path()
        if cls_output is not None:
            fig.savefig(cls_output, dpi=150, bbox_inches="tight")
            print(f"  [OK] Saved to {cls_output}")

        if self.show_plot:
            plt.show()
        else:
            plt.close(fig)

        print()

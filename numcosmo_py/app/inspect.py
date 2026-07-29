#
# inspect.py
#
# Tue May 13 12:00:00 2026
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
# with this program. If not, see <http://www.gnu.org/licenses/>.

"""NumCosmo APP commands for read-only experiment inspection."""

from __future__ import annotations

import dataclasses
import math
from enum import StrEnum, auto
from pathlib import Path
from typing import Annotated

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import typer
from matplotlib.patches import Circle
from rich.table import Table

from numcosmo_py import Ncm, Nc
from .loading import LoadExperiment


@dataclasses.dataclass(kw_only=True)
class InspectExperiment(LoadExperiment):
    """Load experiment and helper objects for inspection commands."""

    def __post_init__(self) -> None:
        """Load experiment objects using the standard app loading order."""
        super().__post_init__()

        # Keep dataset loading behavior consistent with the rest of the app while
        # exposing it as an attribute for inspection commands.
        ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
        dataset_file = self.experiment.with_suffix(".dataset.gvar")
        if dataset_file.exists():
            dataset = ser.from_binfile(dataset_file.absolute().as_posix())
            if not isinstance(dataset, Ncm.Dataset):
                raise RuntimeError(f"Invalid dataset file {dataset_file}.")
            self.dataset = dataset
        else:
            self.dataset = self.likelihood.peek_dataset()

    def _find_cluster_ncounts(self) -> Nc.DataClusterNCountsGauss:
        """Return the first DataClusterNCountsGauss in the dataset."""
        for i in range(self.dataset.get_length()):
            data_i = self.dataset.peek_data(i)
            if isinstance(data_i, Nc.DataClusterNCountsGauss):
                return data_i

        raise RuntimeError(
            "No Nc.DataClusterNCountsGauss data object found in dataset."
        )

    def _find_wl_factor(self) -> Nc.DataClusterWLFactor:
        """Return the first DataClusterWLFactor in the dataset."""
        for i in range(self.dataset.get_length()):
            data_i = self.dataset.peek_data(i)
            if isinstance(data_i, Nc.DataClusterWLFactor):
                return data_i

        raise RuntimeError("No Nc.DataClusterWLFactor data object found in dataset.")

    def _find_galaxy_shape_pop(self) -> Nc.GalaxyShapePop:
        """Return the first GalaxyShapePop model in the model set."""
        for i in range(self.mset.nmodels()):
            model = self.mset.peek_array_pos(i)
            if isinstance(model, Nc.GalaxyShapePop):
                return model

        raise RuntimeError("No Nc.GalaxyShapePop model found in the model set.")


@dataclasses.dataclass(kw_only=True)
class InspectSummary(InspectExperiment):
    """Print a concise and general summary of experiment contents."""

    show_model_params: Annotated[
        bool,
        typer.Option(
            "--show-model-params/--no-show-model-params",
            help=(
                "Show detailed model parameter table including values, bounds, "
                "scale, and fit type."
            ),
        ),
    ] = False

    def __post_init__(self) -> None:
        """Print high-level dataset and model inventory diagnostics."""
        super().__post_init__()

        dset = self.dataset
        m2lnl = self.likelihood.m2lnL_val(self.mset)

        main_table = Table(title="Experiment summary", show_header=False, expand=False)
        main_table.add_column(style="bold bright_cyan")
        main_table.add_column(style="bold bright_green")
        main_table.add_row("Likelihood m2lnL", f"{m2lnl:.8g}")
        main_table.add_row("Dataset entries", f"{dset.get_length()}")
        main_table.add_row("Models present", f"{self.mset.nmodels()}")
        self.console.print(main_table)

        data_table = Table(title="Dataset likelihood terms", expand=False)
        data_table.add_column("Index", justify="right", style="bold bright_cyan")
        data_table.add_column("Type", style="bold")
        data_table.add_column("Length", justify="right")
        data_table.add_column("Description")
        for i in range(dset.get_length()):
            data_i = dset.peek_data(i)
            data_type = data_i.__class__.__name__
            data_desc = data_i.get_desc()
            data_len = data_i.get_length()
            data_table.add_row(str(i), data_type, str(data_len), data_desc)
        self.console.print(data_table)

        models_table = Table(title="Model set", expand=False)
        models_table.add_column("Index", justify="right", style="bold bright_cyan")
        models_table.add_column("Namespace")
        models_table.add_column("Name")
        models_table.add_column("Nick")
        for i in range(self.mset.nmodels()):
            model = self.mset.peek_array_pos(i)
            model_mid = model.id()
            model_ns = Ncm.MSet.get_ns_by_id(model_mid)
            models_table.add_row(str(i), model_ns, model.name(), model.nick())
        self.console.print(models_table)

        if self.show_model_params:
            for i in range(self.mset.nmodels()):
                model = self.mset.peek_array_pos(i)
                model_mid = model.id()
                model_ns = Ncm.MSet.get_ns_by_id(model_mid)

                param_table = Table(
                    title=f"Parameters: {model_ns} ({model.nick()})",
                    expand=False,
                )
                param_table.add_column("#", justify="right", style="bold bright_cyan")
                param_table.add_column("Name")
                param_table.add_column("Symbol")
                param_table.add_column("Value", justify="right")
                param_table.add_column("Lower", justify="right")
                param_table.add_column("Upper", justify="right")
                param_table.add_column("Scale", justify="right")
                param_table.add_column("Fit type")

                for p in range(model.len()):
                    ftype = model.param_get_ftype(p)
                    ftype_str = ftype.name if hasattr(ftype, "name") else str(ftype)
                    param_table.add_row(
                        str(p),
                        model.param_name(p),
                        model.param_symbol(p),
                        f"{model.param_get(p):.6g}",
                        f"{model.param_get_lower_bound(p):.6g}",
                        f"{model.param_get_upper_bound(p):.6g}",
                        f"{model.param_get_scale(p):.6g}",
                        ftype_str,
                    )

                self.console.print(param_table)

        self.close_logging()


@dataclasses.dataclass(kw_only=True)
class InspectClusterNCounts(InspectExperiment):
    """Plot cluster-count data vector and covariance diagnostics."""

    output_prefix: Annotated[
        Path | None,
        typer.Option(
            "--output-prefix",
            help=(
                "Prefix for output image files. "
                "If set to path/to/foo, figures are saved as "
                "foo_data_vector.png, foo_covariance.png, and foo_sij.png."
            ),
        ),
    ] = None

    show_plot: Annotated[
        bool,
        typer.Option(
            "--show-plot/--no-show-plot",
            help="Display figures interactively.",
        ),
    ] = True

    cmap: Annotated[
        str,
        typer.Option(help="Matplotlib colormap name for heatmaps."),
    ] = "viridis"

    dpi: Annotated[
        int,
        typer.Option(help="Figure DPI when saving files.", min=50),
    ] = 150

    log_data: Annotated[
        bool,
        typer.Option(
            "--log-data/--linear-data",
            help="Use log10 scale for data-vector heatmap color values.",
        ),
    ] = False

    show_sij: Annotated[
        bool,
        typer.Option(
            "--show-sij/--no-show-sij",
            help="Plot fitting and resample S_ij matrices when present.",
        ),
    ] = True

    def __post_init__(self) -> None:
        """Generate cluster-count diagnostic plots."""
        super().__post_init__()

        ncounts = self._find_cluster_ncounts()

        data_vec = ncounts.peek_mean().to_numpy()
        cov_matrix, _ = ncounts.compute_cov(self.mset)
        cov = cov_matrix.to_numpy()

        z_obs = ncounts.get_z_obs().to_numpy()
        lnm_obs = ncounts.get_lnM_obs().to_numpy()

        n_z_bins = max(z_obs.size - 1, 0)
        n_m_bins = max(lnm_obs.size - 1, 0)

        if n_z_bins * n_m_bins == data_vec.size and n_z_bins > 0 and n_m_bins > 0:
            data_grid = data_vec.reshape((n_z_bins, n_m_bins))
        else:
            self.console.print(
                "Could not infer (z, lnM) 2D grid from bins; plotting flattened vector."
            )
            data_grid = data_vec.reshape((1, data_vec.size))

        self._plot_data_vector(data_grid)
        self._plot_covariance(cov)

        if self.show_sij:
            fit_sij = ncounts.get_s_matrix()
            resample_sij = ncounts.get_resample_s_matrix()
            if fit_sij is not None or resample_sij is not None:
                self._plot_sij(fit_sij, resample_sij)
            else:
                self.console.print("No S_ij matrices available to plot.")

        if self.show_plot:
            plt.show()

        self.close_logging()

    def _save_or_note(self, fig: plt.Figure, suffix: str) -> None:
        """Save a figure when output prefix is provided."""
        if self.output_prefix is None:
            return

        out_file = self.output_prefix.parent / f"{self.output_prefix.name}_{suffix}.png"
        fig.savefig(out_file, dpi=self.dpi, bbox_inches="tight")
        self.console.print(f"Saved figure: {out_file}")

    def _plot_data_vector(self, data_grid: np.ndarray) -> None:
        """Plot data-vector heatmap."""
        plot_data = data_grid
        if self.log_data:
            # Clip at tiny positive values to keep finite log colors.
            plot_data = np.log10(np.clip(data_grid, 1.0e-16, None))

        fig, ax = plt.subplots(figsize=(8.0, 5.5))
        im = ax.imshow(plot_data, origin="lower", aspect="auto", cmap=self.cmap)
        cbar = fig.colorbar(im, ax=ax)
        cbar.set_label("log10(counts)" if self.log_data else "counts")
        ax.set_xlabel("lnM_obs bin index")
        ax.set_ylabel("z bin index")
        ax.set_title("Cluster-count data vector")
        fig.tight_layout()
        self._save_or_note(fig, "data_vector")

    def _plot_covariance(self, cov: np.ndarray) -> None:
        """Plot covariance and correlation heatmaps."""
        std = np.sqrt(np.clip(np.diag(cov), 1.0e-300, None))
        corr = cov / np.outer(std, std)

        fig, (ax_cov, ax_corr) = plt.subplots(1, 2, figsize=(12.5, 5.2))

        im_cov = ax_cov.imshow(cov, origin="lower", aspect="auto", cmap=self.cmap)
        fig.colorbar(im_cov, ax=ax_cov, fraction=0.046, pad=0.04)
        ax_cov.set_title("Covariance")
        ax_cov.set_xlabel("bin index")
        ax_cov.set_ylabel("bin index")

        im_corr = ax_corr.imshow(
            corr,
            origin="lower",
            aspect="auto",
            cmap="coolwarm",
            vmin=-1.0,
            vmax=1.0,
        )
        fig.colorbar(im_corr, ax=ax_corr, fraction=0.046, pad=0.04)
        ax_corr.set_title("Correlation")
        ax_corr.set_xlabel("bin index")
        ax_corr.set_ylabel("bin index")

        fig.tight_layout()
        self._save_or_note(fig, "covariance")

    def _plot_sij(
        self,
        fit_sij: Ncm.Matrix | None,
        resample_sij: Ncm.Matrix | None,
    ) -> None:
        """Plot available S_ij matrices."""
        fit_arr = fit_sij.to_numpy() if fit_sij is not None else None
        resample_arr = resample_sij.to_numpy() if resample_sij is not None else None

        n_panels = int(fit_arr is not None) + int(resample_arr is not None)
        fig, axes = plt.subplots(1, n_panels, figsize=(6.0 * n_panels, 5.0))

        if n_panels == 1:
            axes = [axes]

        i = 0
        if fit_arr is not None:
            im = axes[i].imshow(fit_arr, origin="lower", aspect="auto", cmap=self.cmap)
            fig.colorbar(im, ax=axes[i], fraction=0.046, pad=0.04)
            axes[i].set_title("Fitting S_ij")
            axes[i].set_xlabel("z bin j")
            axes[i].set_ylabel("z bin i")
            i += 1

        if resample_arr is not None:
            im = axes[i].imshow(
                resample_arr, origin="lower", aspect="auto", cmap=self.cmap
            )
            fig.colorbar(im, ax=axes[i], fraction=0.046, pad=0.04)
            axes[i].set_title("Resample S_ij")
            axes[i].set_xlabel("z bin j")
            axes[i].set_ylabel("z bin i")

        fig.tight_layout()
        self._save_or_note(fig, "sij")


class GalaxyShapeIntegrandComponent(StrEnum):
    """Which term of the integrand to plot.

    The noise kernel lives in chi_L-space alone and never depends on g (see
    NcGalaxyShapeFactorFixedQuad's own class docs); only the population/
    Jacobian term does. Its dynamic range is typically ~10x smaller than the
    noise term's, so under TOTAL's shared linear color scale, moving g can
    look like it does nothing even though it visibly does once isolated.
    """

    TOTAL = auto()
    NOISE = auto()
    POPULATION = auto()
    ALL = auto()


@dataclasses.dataclass(kw_only=True)
class InspectGalaxyShapeIntegrand(InspectExperiment):
    """Plot a heatmap of the shear-marginalization integrand for one galaxy.

    Evaluates, over a grid of the lens-plane variable chi_L,
    ln[ N_2(eps_obs - chi_L; sigma_noise^2) * |det J| * P_pop(r) / (2*pi*r) ]
    (r = |chi_I(chi_L, g)|), i.e. the exact same quantity
    NcGalaxyShapeFactorFixedQuad/Quad integrate over chi_L to get the
    marginal likelihood -- without going through any engine's own
    quadrature, so it cannot be biased by whichever grid that engine
    happens to use. Two markers show where each of the two competing pulls
    is centered: eps_obs (the noise kernel's own center) and the chi_L image
    of chi_I=0 (where the population's own density is largest, mapped
    through the shear at the given g) -- how far apart they are, and how
    sharply peaked the population term is between them, is exactly what
    determines how hard this galaxy is for any fixed quadrature grid to
    resolve. Useful for judging quadrature-design choices (or picking g/eps
    combinations to stress-test) by looking at the actual integrand instead
    of guessing from aggregate -2lnL numbers.
    """

    galaxy_index: Annotated[
        int,
        typer.Option(
            help=("Index of the galaxy (row in the WL observation table) to " "plot."),
            min=0,
        ),
    ] = 0

    g1: Annotated[
        float,
        typer.Option(help="Reduced shear g1 (real part) to evaluate at."),
    ] = 0.0

    g2: Annotated[
        float,
        typer.Option(help="Reduced shear g2 (imaginary part) to evaluate at."),
    ] = 0.0

    n_grid: Annotated[
        int,
        typer.Option(help="Heatmap grid resolution per axis.", min=20),
    ] = 400

    component: Annotated[
        GalaxyShapeIntegrandComponent,
        typer.Option(
            help=(
                "Which term to plot: the full integrand (total), just the "
                "g-INDEPENDENT noise kernel (noise), or just the "
                "g-dependent population/Jacobian term (population) -- the "
                "one that actually moves when you change --g1/--g2, easy "
                "to miss under 'total' since the noise term's dynamic "
                "range is usually much larger. 'all' plots all three side "
                "by side, each with its own color scale."
            ),
        ),
    ] = GalaxyShapeIntegrandComponent.TOTAL

    output_prefix: Annotated[
        Path | None,
        typer.Option(
            "--output-prefix",
            help=(
                "Prefix for the output image file. If set to path/to/foo, "
                "the figure is saved as foo_shape_integrand.png."
            ),
        ),
    ] = None

    show_plot: Annotated[
        bool,
        typer.Option(
            "--show-plot/--no-show-plot",
            help="Display the figure interactively.",
        ),
    ] = True

    cmap: Annotated[
        str,
        typer.Option(help="Matplotlib colormap name for the heatmap."),
    ] = "viridis"

    dpi: Annotated[
        int,
        typer.Option(help="Figure DPI when saving the file.", min=50),
    ] = 150

    show_nodes: Annotated[
        bool,
        typer.Option(
            "--show-nodes/--no-show-nodes",
            help=(
                "Overlay NcGalaxyShapeFactorFixedQuad's own quadrature "
                "domain (nodes + boundary) for this galaxy, when the "
                "engine is FixedQuad. Replicates _regen_domain()'s branch "
                "selection and node placement in Python (same "
                "Gauss-Legendre rule, just via numpy instead of GSL) since "
                "the engine keeps its per-galaxy node array private; if "
                "auto-lens-nodes is on, this shows the configured n-lens "
                "upper bound, not the (possibly smaller) per-galaxy "
                "calibrated count actually used."
            ),
        ),
    ] = True

    def __post_init__(self) -> None:
        """Build and plot the integrand heatmap for the requested galaxy."""
        super().__post_init__()

        dcwlf = self._find_wl_factor()
        pop = self._find_galaxy_shape_pop()
        shape_factor = dcwlf.props.shape_factor
        ellip_conv = shape_factor.get_property("ellip-conv")

        # prepare() resolves the shape factor's cached pop/halo_position
        # references from the current mset; update_data_pop() then syncs
        # this specific data object's pop_data fragment (e.g. the Beta
        # population's cached ln-normalization) to that pop's CURRENT
        # alpha/beta -- skipping either leaves pop_data with whatever
        # ldata it was default-constructed with (NOT the experiment's
        # actual population), silently giving a wrong, alpha/beta-blind
        # eval_p() (caught 2026-07-29: matched 1/(r*(1-r)) exactly instead
        # of the real Beta(1.55,1.62) density).
        shape_factor.prepare(self.mset)

        posf = dcwlf.props.position_factor
        zf = dcwlf.props.redshift_factor
        pos_data = Nc.GalaxyPositionFactorData.new(posf, self.mset)
        z_data = Nc.GalaxyRedshiftFactorData.new(zf, self.mset)
        data = Nc.GalaxyShapeFactorData.new(shape_factor, self.mset, pos_data, z_data)
        obs = dcwlf.props.obs
        if not (0 <= self.galaxy_index < obs.len()):
            raise RuntimeError(
                f"--galaxy-index {self.galaxy_index} out of range " f"[0,{obs.len()})."
            )
        data.read_row(obs, self.galaxy_index)
        shape_factor.update_data_pop(data)

        eps_obs = complex(data.epsilon_obs_1, data.epsilon_obs_2)
        sig2 = data.std_noise**2
        g_ncm = Ncm.Complex.new()
        g_ncm.set(self.g1, self.g2)

        if ellip_conv == Nc.GalaxyWLObsEllipConv.TRACE:
            apply_shear = Nc.wl_ellipticity_apply_shear_trace_ptr
            apply_shear_inv = Nc.wl_ellipticity_apply_shear_inv_trace_ptr
            lndet_jac = Nc.wl_ellipticity_lndet_jac_trace_ptr
        else:
            apply_shear = Nc.wl_ellipticity_apply_shear_trace_det_ptr
            apply_shear_inv = Nc.wl_ellipticity_apply_shear_inv_trace_det_ptr
            lndet_jac = Nc.wl_ellipticity_lndet_jac_trace_det_ptr

        zero = Ncm.Complex.new()
        zero.set(0.0, 0.0)
        chi_l0 = Ncm.Complex.new()
        apply_shear(g_ncm, zero, chi_l0)
        pop_pull = complex(chi_l0.Re(), chi_l0.Im())

        xs = np.linspace(-1.05, 1.05, self.n_grid)
        xx, yy = np.meshgrid(xs, xs)
        chi_l_grid = xx + 1j * yy
        mask = np.abs(chi_l_grid) <= 1.0
        chi_l_flat = chi_l_grid[mask]

        chi_obs_ncm = Ncm.Complex.new()
        chi_i_ncm = Ncm.Complex.new()
        r_arr = np.empty(chi_l_flat.size)
        lnjac_arr = np.empty(chi_l_flat.size)

        for k, chi_l in enumerate(chi_l_flat):
            chi_obs_ncm.set(chi_l.real, chi_l.imag)
            apply_shear_inv(g_ncm, chi_obs_ncm, chi_i_ncm)
            r_arr[k] = math.hypot(chi_i_ncm.Re(), chi_i_ncm.Im())
            lnjac_arr[k] = lndet_jac(g_ncm, chi_obs_ncm)

        r_safe = np.clip(r_arr, 1.0e-12, None)
        ln_p_pop = np.log(
            np.clip(
                [pop.eval_p(data.pop_data, float(r)) for r in r_safe],
                1.0e-300,
                None,
            )
        )
        ln_area = np.log(2.0 * np.pi * r_safe)
        ln_pop_term_flat = ln_p_pop - ln_area + lnjac_arr
        ln_noise_flat = -np.abs(eps_obs - chi_l_flat) ** 2 / (2.0 * sig2) - math.log(
            2.0 * np.pi * sig2
        )

        values_and_labels = {
            GalaxyShapeIntegrandComponent.NOISE: (
                ln_noise_flat,
                "ln(noise kernel) [g-independent]",
            ),
            GalaxyShapeIntegrandComponent.POPULATION: (
                ln_pop_term_flat,
                "ln(population * |det J|) [the g-dependent term]",
            ),
            GalaxyShapeIntegrandComponent.TOTAL: (
                ln_pop_term_flat + ln_noise_flat,
                "ln(integrand) [sum]",
            ),
        }

        if self.show_nodes and isinstance(shape_factor, Nc.GalaxyShapeFactorFixedQuad):
            chi_l_mat, node_weights = shape_factor.peek_domain(pop, data)
            nodes = chi_l_mat.to_numpy()
            node_w = node_weights.to_numpy()
        elif self.show_nodes:
            self.console.print(
                "--show-nodes requested but the shape factor is "
                f"{type(shape_factor).__name__}, not FixedQuad -- skipping."
            )
            nodes = None
            node_w = None
        else:
            nodes = None
            node_w = None

        def _plot_panel(
            ax: plt.Axes,
            component: GalaxyShapeIntegrandComponent,
            norm: mcolors.Normalize | None = None,
            draw_colorbar: bool = True,
        ) -> plt.cm.ScalarMappable:
            values_flat, label = values_and_labels[component]
            plot_grid = np.full(xx.shape, np.nan)
            plot_grid[mask] = values_flat

            im = ax.imshow(
                plot_grid,
                origin="lower",
                extent=(xs[0], xs[-1], xs[0], xs[-1]),
                cmap=self.cmap,
                norm=norm,
            )
            if draw_colorbar:
                ax.figure.colorbar(im, ax=ax, label=label, fraction=0.046, pad=0.04)
            ax.add_patch(
                Circle((0.0, 0.0), 1.0, fill=False, edgecolor="white", linewidth=1.0)
            )
            ax.scatter(
                [eps_obs.real],
                [eps_obs.imag],
                marker="x",
                color="red",
                s=80,
                label="eps_obs (noise pull)",
            )
            ax.scatter(
                [pop_pull.real],
                [pop_pull.imag],
                marker="*",
                color="cyan",
                s=120,
                label="chi_I=0 image (pop pull)",
            )

            if nodes is not None:
                w_max = node_w.max() if node_w.size and node_w.max() > 0 else 1.0
                ax.scatter(
                    nodes[:, 0],
                    nodes[:, 1],
                    marker=".",
                    color="white",
                    edgecolors="black",
                    linewidths=0.3,
                    s=10.0 + 60.0 * (node_w / w_max),
                    alpha=0.7,
                    label=f"FixedQuad nodes (n={nodes.shape[0]})",
                )

            ax.set_xlabel("Re(chi_L)")
            ax.set_ylabel("Im(chi_L)")
            ax.set_title(str(component))
            ax.set_xlim(xs[0], xs[-1])
            ax.set_ylim(xs[0], xs[-1])
            ax.set_aspect("equal")

            return im

        suptitle = (
            f"Galaxy #{self.galaxy_index}: g=({self.g1:+.3f},{self.g2:+.3f}), "
            f"eps_obs=({eps_obs.real:+.3f},{eps_obs.imag:+.3f}), "
            f"std_noise={data.std_noise:.3f}"
        )

        if self.component == GalaxyShapeIntegrandComponent.ALL:
            order = (
                GalaxyShapeIntegrandComponent.NOISE,
                GalaxyShapeIntegrandComponent.TOTAL,
                GalaxyShapeIntegrandComponent.POPULATION,
            )
            all_values = [values_and_labels[c][0] for c in order]
            all_concat = np.concatenate(all_values)
            vmax = all_concat.max()
            # The true min is dominated by the noise term's far corners --
            # for a small std_noise the kernel falls off so steeply with
            # distance that most of the disc's AREA (not its probability
            # mass) sits far below the peak, stretching the shared scale
            # down to ln(value)~-80+ and washing out the informative region
            # near the peak. A percentile of pixel values doesn't fix this
            # (it reflects area, not relevance), so clip to a fixed range
            # below the peak instead: exp(-40) is already ~4e-18 of the
            # peak, well past anything visually or physically relevant.
            vmin = max(all_concat.min(), vmax - 40.0)
            # A single shared Normalize instance -- not just equal vmin/vmax
            # passed to three separate imshow calls -- so that interactively
            # rescaling the colorbar (e.g. scrolling/zooming on it) mutates
            # the one norm object all three images reference, instead of
            # only the first panel's own private norm.
            shared_norm = mcolors.Normalize(vmin=vmin, vmax=vmax)

            fig, axes = plt.subplots(
                1, 3, figsize=(19.5, 6.5), sharex=True, sharey=True
            )
            ims = [
                _plot_panel(ax, component, norm=shared_norm, draw_colorbar=False)
                for ax, component in zip(axes, order)
            ]
            fig.colorbar(
                ims[0],
                ax=axes.tolist(),
                label="ln(value) [shared scale across all three panels]",
                fraction=0.025,
                pad=0.01,
            )
            axes[0].legend(loc="upper right", framealpha=0.7, fontsize=8)
            fig.suptitle(suptitle, fontsize=10)
        else:
            fig, ax = plt.subplots(figsize=(7.5, 6.5))
            _plot_panel(ax, self.component)
            ax.legend(loc="upper right", framealpha=0.7)
            ax.set_title(f"{suptitle}\n[{self.component}]", fontsize=10)
            fig.tight_layout()

        if self.output_prefix is not None:
            suffix = (
                "shape_integrand"
                if self.component == GalaxyShapeIntegrandComponent.TOTAL
                else f"shape_integrand_{self.component}"
            )
            out_file = (
                self.output_prefix.parent / f"{self.output_prefix.name}_{suffix}.png"
            )
            fig.savefig(out_file, dpi=self.dpi, bbox_inches="tight")
            self.console.print(f"Saved figure: {out_file}")

        if self.show_plot:
            plt.show()

        self.close_logging()

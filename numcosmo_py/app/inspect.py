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
from pathlib import Path
from typing import Annotated

import matplotlib.pyplot as plt
import numpy as np
import typer
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

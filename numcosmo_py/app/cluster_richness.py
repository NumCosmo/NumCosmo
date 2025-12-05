#
# cluster_richness.py
#
# Copyright (C) 2025 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

"""NumCosmo APP subcommand for cluster richness analysis.

This module provides CLI commands for analyzing cluster mass-richness
scaling relations using the cluster_richness analysis package.
"""

import dataclasses
from pathlib import Path
from typing import Annotated, Optional

import numpy as np
import typer
from astropy.table import Table
from rich.panel import Panel
from rich.table import Table as RichTable


from numcosmo_py import Nc
from numcosmo_py.analysis.cluster_richness import (
    RichnessModelType,
    CutAnalyzer,
    ClusterData,
    MockStudy,
    get_model_param_names,
    model_params_as_list,
    PARAM_FORMAT,
    compute_binned_statistics,
    plot_diagnostic_summary,
    mean_lnR_truncated,
    std_lnR_truncated,
    CutAnalysisResult,
)
from .logging import AppLogging


@dataclasses.dataclass(kw_only=True)
class RunClusterRichnessAnalysis(AppLogging):
    """Run cluster mass-richness scaling relation analysis.

    Analyzes cluster data from a FITS file to extract mass-richness
    scaling relation parameters. Supports mock studies for bias assessment.
    """

    # Required argument: FITS file path
    fits_file: Annotated[
        Path,
        typer.Argument(
            help="Path to FITS file containing cluster data.",
            exists=True,
            readable=True,
        ),
    ]

    # Column names with defaults
    mass_column: Annotated[
        str,
        typer.Option(
            "--mass-col",
            "-m",
            help="Name of the halo mass column in the FITS file.",
        ),
    ] = "halo_mass"

    redshift_column: Annotated[
        str,
        typer.Option(
            "--redshift-col",
            "-z",
            help="Name of the redshift column in the FITS file.",
        ),
    ] = "redshift"

    richness_column: Annotated[
        str,
        typer.Option(
            "--richness-col",
            "-r",
            help="Name of the richness column in the FITS file.",
        ),
    ] = "richness"

    sigma_lnR_column: Annotated[
        Optional[str],
        typer.Option(
            "--sigma-lnR-col",
            help=(
                "Name of the richness uncertainty column in the FITS file. "
                "If not provided, will look for '{richness_column}_err'."
            ),
        ),
    ] = None

    ignore_noise: Annotated[
        bool,
        typer.Option(
            "--ignore-noise",
            help=(
                "Ignore richness measurement uncertainties in the analysis. "
                "By default, uncertainties are included."
            ),
            is_flag=True,
        ),
    ] = False

    # Model configuration
    model_type: Annotated[
        RichnessModelType,
        typer.Option(
            "--model-type",
            "-t",
            help="Type of richness model to use.",
        ),
    ] = RichnessModelType.ASCASO

    # Analysis options
    cuts: Annotated[
        Optional[str],
        typer.Option(
            "--cuts",
            "-c",
            help=(
                "Comma-separated richness cut values (not log). "
                "Default: '5,10,15,20,30'"
            ),
        ),
    ] = None

    n_mocks: Annotated[
        int,
        typer.Option(
            "--n-mocks",
            "-n",
            help="Number of mock realizations for bias study.",
            min=0,
        ),
    ] = 10000

    n_bootstrap: Annotated[
        int,
        typer.Option(
            "--n-bootstrap",
            "-b",
            help="Number of bootstrap samples per analysis.",
            min=0,
        ),
    ] = 200

    seed: Annotated[
        int,
        typer.Option(
            "--seed",
            "-s",
            help="Random seed for reproducibility.",
        ),
    ] = 42

    # FITS HDU option
    hdu: Annotated[
        int,
        typer.Option(
            "--hdu",
            help="HDU number to read from FITS file.",
        ),
    ] = 1

    # Output prefix
    output_prefix: Annotated[
        str,
        typer.Option(
            "--output-prefix",
            "-o",
            help="Prefix for output files (database, catalogs).",
        ),
    ] = "cluster_richness"

    run_mocks: Annotated[
        bool,
        typer.Option(
            "--run-mocks",
            help="Run mock study for bias assessment (requires --run-analysis).",
            is_flag=True,
        ),
    ] = False

    run_diagnostics: Annotated[
        bool,
        typer.Option(
            "--run-diagnostics",
            "-d",
            help="Run diagnostic analysis on the results.",
            is_flag=True,
        ),
    ] = False

    n_bins: Annotated[
        int,
        typer.Option(
            "--n-bins",
            help="Number of bins for the diagnostic analysis.",
        ),
    ] = 20

    show_plots: Annotated[
        bool,
        typer.Option(
            "--show-plots",
            "-p",
            help="Display diagnostic plots (requires --run-analysis).",
            is_flag=True,
        ),
    ] = False

    compute_mcmc: Annotated[
        bool,
        typer.Option(
            "--mcmc/--no-mcmc",
            help="Compute MCMC for uncertainty estimation.",
        ),
    ] = False

    compute_bootstrap: Annotated[
        bool,
        typer.Option(
            "--bootstrap/--no-bootstrap",
            help="Compute bootstrap for uncertainty estimation.",
        ),
    ] = False

    def __post_init__(self) -> None:
        """Run the cluster richness analysis."""
        super().__post_init__()

        # Load and display data summary
        data, cuts_array = self._load_data()

        # Print summary
        self._print_summary_only(data, cuts_array)

        # Run real data analysis if requested
        real_results = None
        model_fiducial = None

        real_results, model_fiducial = self._run_real_analysis(data, cuts_array)

        # Run mock study if requested
        if self.run_mocks and real_results is not None and model_fiducial is not None:
            self._run_mock_study(data, cuts_array, real_results, model_fiducial)

        # Run diagnostics if requested
        if (self.run_diagnostics or self.show_plots) and real_results is not None:
            self._run_diagnostics(data, cuts_array, real_results)

    def _load_data(
        self,
    ) -> tuple[ClusterData, np.ndarray]:
        """Load data from FITS file and return ClusterData.

        :return: Tuple of (ClusterData, cuts_array)
        """
        # Parse cuts
        if self.cuts is None:
            cuts_array = np.log(np.array([5, 10, 15, 20, 30]))
        else:
            cuts_list = [float(x.strip()) for x in self.cuts.split(",")]
            cuts_array = np.log(np.array(cuts_list))

        # Load data from FITS file
        self.console.print(
            f"[bold cyan]Loading data from {self.fits_file}...[/bold cyan]"
        )
        table_halos = Table.read(self.fits_file, hdu=self.hdu)

        # Determine sigma_lnR column name
        sigma_col = (
            self.sigma_lnR_column
            if self.sigma_lnR_column is not None
            else f"{self.richness_column}_err"
        )

        # Extract columns
        try:
            mass = np.array(table_halos[self.mass_column])
            z = np.array(table_halos[self.redshift_column])
            richness = np.array(table_halos[self.richness_column])
            richness_err = np.array(table_halos[sigma_col])
        except KeyError as e:
            available_cols = ", ".join(table_halos.colnames)
            self.console.print(f"[bold red]Error: Column {e} not found.[/bold red]")
            self.console.print(f"[yellow]Available columns: {available_cols}[/yellow]")
            raise typer.Exit(code=1) from e

        # Convert to log values
        lnM = np.log(mass)
        lnR = np.log(richness)
        if self.ignore_noise:
            sigma_lnR = np.zeros_like(lnR)
        else:
            sigma_lnR = richness_err / richness  # Propagate error to log space

        data = ClusterData(lnM=lnM, z=z, lnR=lnR, sigma_lnR=sigma_lnR)
        return data, cuts_array

    def _print_summary_only(
        self,
        data: ClusterData,
        cuts_array: np.ndarray,
    ) -> None:
        """Print data summary without running analysis.

        :param data: ClusterData with lnM, z, lnR, sigma_lnR
        :param cuts_array: Array of richness cuts (log values)
        """
        self.console.print("\n[bold cyan]Data Summary[/bold cyan]")
        self.console.print(f"  Total clusters: {len(data)}")
        self.console.print(
            f"  Mass range: [{np.exp(data.lnM.min()):.2e}, "
            f"{np.exp(data.lnM.max()):.2e}]"
        )
        self.console.print(
            f"  Redshift range: [{data.z.min():.3f}, {data.z.max():.3f}]"
        )
        self.console.print(
            f"  Richness range: [{np.exp(data.lnR.min()):.1f}, "
            f"{np.exp(data.lnR.max()):.1f}]"
        )
        self.console.print(
            f"  σ(ln R) range: [{data.sigma_lnR.min():.4f}, "
            f"{data.sigma_lnR.max():.4f}]"
        )

        # Show cluster counts per cut
        table = RichTable(
            title="Clusters per Richness Cut",
            show_header=True,
            header_style="bold magenta",
        )
        table.add_column("Richness Cut", style="green")
        table.add_column("N clusters", style="cyan")
        table.add_column("Fraction", style="yellow")

        for cut in cuts_array:
            n_above = np.sum(data.lnR >= cut)
            fraction = n_above / len(data) * 100
            table.add_row(
                f"λ ≥ {np.exp(cut):.1f}",
                f"{n_above}",
                f"{fraction:.1f}%",
            )

        self.console.print(Panel(table))

    def _run_real_analysis(
        self,
        data: ClusterData,
        cuts_array: np.ndarray,
    ) -> tuple[dict, Nc.ClusterMassRichness]:
        """Run real data analysis.

        :param data: ClusterData with lnM, z, lnR, sigma_lnR
        :param cuts_array: Array of richness cuts (log values)
        :return: Tuple of (results dict, fiducial model)
        """
        self.console.print(f"  Loaded {len(data)} clusters")
        self.console.print(
            f"  Mass range: [{np.exp(data.lnM.min()):.2e}, "
            f"{np.exp(data.lnM.max()):.2e}]"
        )
        self.console.print(
            f"  Redshift range: [{data.z.min():.3f}, {data.z.max():.3f}]"
        )
        self.console.print(
            f"  Richness range: [{np.exp(data.lnR.min()):.1f}, "
            f"{np.exp(data.lnR.max()):.1f}]"
        )

        # Create initial model
        model_init = self._create_initial_model()

        # Analyze real data
        self.console.print("\n[bold magenta]PHASE 1: Real Data Analysis[/bold magenta]")
        real_analyzer = CutAnalyzer(
            data,
            cuts_array.tolist(),
            n_bootstrap=self.n_bootstrap,
            compute_mcmc=self.compute_mcmc,
            compute_bootstrap=self.compute_bootstrap,
            file_prefix=f"{self.output_prefix}_real",
            console=self.console,
        )
        real_results = real_analyzer.analyze(model_init=model_init)

        # Extract best-fit from first cut as fiducial for mocks
        first_cut = cuts_array[0]
        model_fiducial = real_results[first_cut].bestfit
        param_names = get_model_param_names(model_fiducial)
        param_values = model_params_as_list(model_fiducial)

        self.console.print(
            f"\n[green]Best-fit from cut {np.exp(first_cut):.1f} "
            f"(fiducial for mocks):[/green]"
        )
        for name, value in zip(param_names, param_values):
            self.console.print(f"  {name}={value:{PARAM_FORMAT}}")

        return real_results, model_fiducial

    def _run_mock_study(
        self,
        data: ClusterData,
        cuts_array: np.ndarray,
        real_results: dict,
        model_fiducial: Nc.ClusterMassRichness,
    ) -> None:
        """Run mock study for bias assessment.

        :param data: ClusterData with lnM, z, lnR, sigma_lnR
        :param cuts_array: Array of richness cuts (log values)
        :param real_results: Results from real data analysis
        :param model_fiducial: Fiducial model from real data analysis
        """
        if self.n_mocks <= 0:
            self.console.print("[yellow]Skipping mock study: n_mocks is 0.[/yellow]")
            return

        self.console.print("\n[bold magenta]PHASE 2: Mock Study[/bold magenta]")
        mock_study = MockStudy(
            model_fiducial,
            data,
            cuts_array.tolist(),
            n_mocks=self.n_mocks,
            n_bootstrap=self.n_bootstrap,
            file_prefix=f"{self.output_prefix}_mock",
            fiducial_results=real_results,
            console=self.console,
        )
        mock_study.run(seed=self.seed)

        # GOF Analysis
        self.console.print(
            "\n[bold magenta]PHASE 3: Goodness-of-Fit Analysis[/bold magenta]"
        )
        gof_stats = mock_study.compute_gof_statistics()
        mock_study.display_gof_results(gof_stats)

        self.console.print("\n[bold green]Mock study complete![/bold green]")
        self.console.print(
            f"[cyan]Mock study: {len(mock_study.mock_results)} mocks analyzed[/cyan]"
        )

    def _run_diagnostics(
        self,
        data: ClusterData,
        cuts_array: np.ndarray,
        real_results: dict[float, CutAnalysisResult],
    ) -> None:
        """Run diagnostic analysis and optionally show plots.

        :param data: ClusterData with lnM, z, lnR, sigma_lnR
        :param cuts_array: Array of richness cuts (log values)
        :param real_results: Results from real data analysis
        """
        self.console.print("\n[bold magenta]Diagnostic Analysis[/bold magenta]")

        for cut in cuts_array:
            if cut not in real_results:
                continue

            result = real_results[cut]
            model = result.bestfit

            # Filter data for this cut
            data_cut = data.apply_cut(cut)

            # Compute model predictions
            mu = np.array(
                [model.mu(data_cut.lnM[i], data_cut.z[i]) for i in range(len(data_cut))]
            )
            sigma = np.array(
                [
                    model.sigma(data_cut.lnM[i], data_cut.z[i])
                    for i in range(len(data_cut))
                ]
            )

            # Compute total scatter: intrinsic model + catalog uncertainty
            sigma_total = np.sqrt(sigma**2 + data_cut.sigma_lnR**2)

            # Compute predicted mean and std of truncated distribution using total scatter
            mean_pred = mean_lnR_truncated(mu, sigma_total, cut)
            std_pred = std_lnR_truncated(mu, sigma_total, cut)
            assert isinstance(mean_pred, np.ndarray)
            assert isinstance(std_pred, np.ndarray)

            # Compute statistics
            stats = compute_binned_statistics(
                mean_pred,
                std_pred,
                data_cut.lnR,
                mu,
                sigma,
                data_cut.sigma_lnR,
                cut,
                nbins=self.n_bins,
            )

            if self.run_diagnostics:
                self.console.print(f"\n[cyan]Cut λ ≥ {np.exp(cut):.1f}:[/cyan]")
                self.console.print(f"  N clusters: {len(data_cut)}")

                # Print some diagnostic statistics
                sigma_emp = np.nanmean(stats["sample_std"])
                sigma_mod = np.nanmean(stats["bin_std_pred"])
                self.console.print(f"  Mean empirical σ: {sigma_emp:.3f}")
                self.console.print(f"  Mean model σ: {sigma_mod:.3f}")

            if self.show_plots:
                plot_diagnostic_summary(stats, cut, show=True)

    def _create_initial_model(self) -> Nc.ClusterMassRichness:
        """Create the initial model based on model_type.

        :return: Configured NcClusterMassRichness model
        """
        model_init: Nc.ClusterMassRichness
        if self.model_type == RichnessModelType.ASCASO:
            model_init = Nc.ClusterMassAscaso(lnRichness_min=0.0, lnRichness_max=20.0)
            model_init["mup0"] = 4.0
            model_init["mup1"] = 1.0
            model_init["mup2"] = 0.2
            model_init["sigmap0"] = 0.5
            model_init["sigmap1"] = 0.03
            model_init["sigmap2"] = 0.15
        elif self.model_type == RichnessModelType.EXT:
            model_init = Nc.ClusterMassExt(lnRichness_min=0.0, lnRichness_max=20.0)
            model_init["mup0"] = 4.0
            model_init["mup1"] = 1.0
            model_init["mup2"] = 0.1
            model_init["mup3"] = 0.01
            model_init["sigmap0"] = -0.3
            model_init["sigmap1"] = -0.08
            model_init["sigmap2"] = 0.005
        else:
            raise ValueError(f"Unknown model type: {self.model_type}")

        return model_init

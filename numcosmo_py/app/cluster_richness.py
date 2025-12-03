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

from numcosmo_py import Nc
from numcosmo_py.analysis.cluster_richness import (
    RichnessModelType,
    CutAnalyzer,
    MockStudy,
    COMPUTE_MCMC,
    COMPUTE_BOOTSTRAP,
    get_model_param_names,
    model_params_as_list,
    PARAM_FORMAT,
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

    # Optional: skip phases
    skip_mocks: Annotated[
        bool,
        typer.Option(
            "--skip-mocks",
            help="Skip mock study phase (only analyze real data).",
        ),
    ] = False

    compute_mcmc: Annotated[
        bool,
        typer.Option(
            "--mcmc/--no-mcmc",
            help="Compute MCMC for uncertainty estimation.",
        ),
    ] = COMPUTE_MCMC

    compute_bootstrap: Annotated[
        bool,
        typer.Option(
            "--bootstrap/--no-bootstrap",
            help="Compute bootstrap for uncertainty estimation.",
        ),
    ] = COMPUTE_BOOTSTRAP

    def __post_init__(self) -> None:
        """Run the cluster richness analysis."""
        super().__post_init__()
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

        # Extract columns
        try:
            mass = np.array(table_halos[self.mass_column])
            z = np.array(table_halos[self.redshift_column])
            richness = np.array(table_halos[self.richness_column])
        except KeyError as e:
            available_cols = ", ".join(table_halos.colnames)
            self.console.print(f"[bold red]Error: Column {e} not found.[/bold red]")
            self.console.print(f"[yellow]Available columns: {available_cols}[/yellow]")
            raise typer.Exit(code=1) from e

        # Convert to log values
        lnM = np.log(mass)
        lnR = np.log(richness)

        self.console.print(f"  Loaded {len(lnM)} clusters")
        self.console.print(
            f"  Mass range: [{np.exp(lnM.min()):.2e}, {np.exp(lnM.max()):.2e}]"
        )
        self.console.print(f"  Redshift range: [{z.min():.3f}, {z.max():.3f}]")
        self.console.print(
            f"  Richness range: [{np.exp(lnR.min()):.1f}, {np.exp(lnR.max()):.1f}]"
        )

        # Create initial model
        model_init = self._create_initial_model()

        # Phase 1: Analyze real data
        self.console.print("\n[bold magenta]PHASE 1: Real Data Analysis[/bold magenta]")
        real_analyzer = CutAnalyzer(
            lnM,
            z,
            lnR,
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

        # Phase 2: Mock analysis (if not skipped)
        if not self.skip_mocks and self.n_mocks > 0:
            self.console.print("\n[bold magenta]PHASE 2: Mock Study[/bold magenta]")
            mock_study = MockStudy(
                model_fiducial,
                lnM,
                z,
                cuts_array.tolist(),
                n_mocks=self.n_mocks,
                n_bootstrap=self.n_bootstrap,
                file_prefix=f"{self.output_prefix}_mock",
                fiducial_results=real_results,
                console=self.console,
            )
            mock_study.run(seed=self.seed)

            # Phase 3: GOF Analysis
            self.console.print(
                "\n[bold magenta]PHASE 3: Goodness-of-Fit Analysis[/bold magenta]"
            )
            gof_stats = mock_study.compute_gof_statistics()
            mock_study.display_gof_results(gof_stats)

            self.console.print("\n[bold green]Analysis complete![/bold green]")
            self.console.print(
                f"[cyan]Real data: {len(real_results)} cuts analyzed[/cyan]"
            )
            self.console.print(
                f"[cyan]Mock study: {len(mock_study.mock_results)} mocks analyzed"
                "[/cyan]"
            )
        else:
            self.console.print(
                "\n[bold green]Real data analysis complete![/bold green]"
            )
            self.console.print(f"[cyan]{len(real_results)} cuts analyzed[/cyan]")

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

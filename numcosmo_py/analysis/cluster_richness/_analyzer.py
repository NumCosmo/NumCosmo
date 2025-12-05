#
# _analyzer.py
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

"""CutAnalyzer class for cluster mass-richness analysis.

This module provides the main analysis pipeline for studying cluster
mass-richness relations with progressive richness cuts.
"""

from dataclasses import dataclass
from typing import Optional

import numpy as np
from rich.console import Console
from rich.panel import Panel
from rich.table import Table as RichTable

from numcosmo_py import Nc, Ncm

from ._parameters import (
    CutAnalysisResult,
    dup_model,
    model_params_as_list,
    model_params_from_list,
    get_model_param_names,
)
from ._utils import setup_model_fit_params, PARAM_FORMAT


@dataclass
class ClusterData:
    """Container for cluster observable data.

    This dataclass holds all the observables for a cluster sample:
    - Independent variables: lnM (log mass), z (redshift)
    - Dependent variable: lnR (log richness)
    - Observational uncertainty: sigma_lnR (catalog uncertainty on log richness)

    All arrays must have the same length.
    """

    lnM: np.ndarray
    """Log mass array (natural log of solar masses)"""

    z: np.ndarray
    """Redshift array"""

    lnR: np.ndarray
    """Log richness array (natural log)"""

    sigma_lnR: np.ndarray
    """Catalog uncertainty on log richness (natural log units)"""

    def __post_init__(self):
        """Validate that all arrays have the same length."""
        lengths = {
            "lnM": len(self.lnM),
            "z": len(self.z),
            "lnR": len(self.lnR),
            "sigma_lnR": len(self.sigma_lnR),
        }
        if len(set(lengths.values())) != 1:
            raise ValueError(f"All arrays must have the same length. Got: {lengths}")

    def __len__(self) -> int:
        """Return the number of clusters in the dataset."""
        return len(self.lnM)

    def apply_cut(self, lnR_cut: float) -> "ClusterData":
        """Apply a richness cut and return a new ClusterData instance.

        :param lnR_cut: Richness cut in log units
        :return: New ClusterData with only clusters above the cut
        """
        mask = self.lnR >= lnR_cut
        return ClusterData(
            lnM=self.lnM[mask],
            z=self.z[mask],
            lnR=self.lnR[mask],
            sigma_lnR=self.sigma_lnR[mask],
        )


def _get_default_console() -> Console:
    """Get or create the default console instance.

    :return: Default Console instance
    """
    return Console()


class CutAnalyzer:
    """Analyze cluster mass-richness relation with progressive cuts.

    This class handles:
    - Setting up the mass-richness model (any NcClusterMassRichness subclass)
    - Applying richness cuts to the dataset
    - Running MCMC and bootstrap analyses
    - Storing and displaying results
    """

    def __init__(
        self,
        data: ClusterData,
        cuts: list[float],
        n_bootstrap: int = 100,
        compute_mcmc: bool = False,
        compute_bootstrap: bool = False,
        file_prefix: str | None = None,
        sample_desc: str = "Sample",
        verbose: bool = True,
        console: Optional[Console] = None,
    ):
        """Initialize the analyzer.

        :param data: ClusterData object containing lnM, z, lnR, and sigma_lnR
        :param cuts: List of richness cuts (in log units)
        :param n_bootstrap: Number of bootstrap resamples (default: 100)
        :param compute_mcmc: Whether to run MCMC (default: False)
        :param compute_bootstrap: Whether to run bootstrap (default: False)
        :param file_prefix: Prefix for output files (default: None)
        :param sample_desc: Description of the sample for display (default: "Sample")
        :param verbose: Whether to print progress (default: True)
        :param console: Rich Console for output (default: creates new Console)
        """
        self.data = data
        self.cuts = cuts
        self.n_bootstrap = n_bootstrap
        self.compute_mcmc = compute_mcmc
        self.compute_bootstrap = compute_bootstrap
        self.file_prefix = file_prefix
        self.sample_desc = sample_desc
        self.verbose = verbose
        self.console = console if console is not None else _get_default_console()
        self.results: dict[float, CutAnalysisResult] = {}

        # Will be initialized in _initialize_model
        self.cluster_m: Nc.ClusterMassRichness | None = None

    def _initialize_model(self, model_init: Nc.ClusterMassRichness) -> None:
        """Initialize the cluster mass model from a template.

        :param model_init: Template model with initial parameters
        """
        # Duplicate the model to avoid modifying the original
        self.cluster_m = dup_model(model_init)
        setup_model_fit_params(self.cluster_m)

    def _setup_real_data(self) -> Nc.DataClusterMassRich:
        """Create NumCosmo data object from ClusterData.

        :return: DataClusterMassRich object
        """
        nc_data = Nc.DataClusterMassRich.new()

        lnM_v = Ncm.Vector.new_array(self.data.lnM)
        z_v = Ncm.Vector.new_array(self.data.z)
        lnR_v = Ncm.Vector.new_array(self.data.lnR)
        sigma_lnR_v = Ncm.Vector.new_array(self.data.sigma_lnR)

        nc_data.set_data(lnM_v, z_v, lnR_v, sigma_lnR_v)
        return nc_data

    def _analyze_cut(
        self, data: Nc.DataClusterMassRich, cut: float
    ) -> CutAnalysisResult:
        """Analyze a single cut.

        :param data: The cluster data
        :param mset: The model set
        :param cut: The richness cut value (log units)
        :return: Analysis result for this cut
        """
        assert self.cluster_m is not None
        cluster_m = dup_model(self.cluster_m)
        mset = Ncm.MSet.new_array([cluster_m])
        mset.prepare_fparam_map()

        # Apply cut
        data.apply_cut(cut)
        cluster_m["cut"] = cut
        n_clusters = int(np.sum(self.data.lnR >= cut))

        # Best-fit
        if self.verbose:
            self.console.print("    [yellow]→ Running fit...[/yellow]", end="")
        dset = Ncm.Dataset.new()
        dset.append_data(data)
        likelihood = Ncm.Likelihood.new(dset)

        fit = Ncm.Fit.factory(
            Ncm.FitType.NLOPT,
            "ln-neldermead",
            likelihood,
            mset,
            Ncm.FitGradType.NUMDIFF_CENTRAL,
        )
        fit.set_params_reltol(1.0e-11)
        fit.set_m2lnL_reltol(1.0e-11)
        fit.run_restart(Ncm.FitRunMsgs.NONE, 1.0e-3, 0.0, None, None)
        bestfit_model = dup_model(cluster_m)
        m2lnL = fit.peek_state().get_m2lnL_curval()

        if self.verbose:
            self.console.print(" [green]✓[/green]")

        # MCMC (if requested)
        if self.compute_mcmc:
            if self.verbose:
                self.console.print(
                    "    [yellow]→ Running MCMC (1000 walkers)...[/yellow]", end=""
                )
            nwalkers = 1000
            init_sampler = Ncm.MSetTransKernGauss.new(0)
            init_sampler.set_mset(mset)
            init_sampler.set_prior_from_mset()
            init_sampler.set_cov_from_rescale(1.0)
            nparams = mset.fparam_len()
            apes_walker = Ncm.FitESMCMCWalkerAPES.new(nwalkers, nparams)
            esmcmc = Ncm.FitESMCMC.new(
                fit, nwalkers, init_sampler, apes_walker, Ncm.FitRunMsgs.NONE
            )

            esmcmc.set_nthreads(12)
            # Set data file if prefix is provided
            if self.file_prefix is not None:
                mcmc_file = f"{self.file_prefix}_{np.exp(cut):.2f}_mcmc.fits"
                esmcmc.set_data_file(mcmc_file)

            esmcmc.start_run()
            esmcmc.run(400)
            esmcmc.end_run()

            mcat = esmcmc.get_catalog()
            rows = np.array(
                [mcat.peek_row(i).dup_array()[1:] for i in range(150, mcat.len())]
            )
            # Compute mean from rows
            mcmc_mean_model = dup_model(cluster_m)
            model_params_from_list(mcmc_mean_model, np.mean(rows, axis=0).tolist())

            mcmc_median_model = dup_model(cluster_m)
            model_params_from_list(mcmc_median_model, np.median(rows, axis=0).tolist())

            if self.verbose:
                self.console.print(" [green]✓[/green]")
        else:
            mcmc_mean_model = dup_model(bestfit_model)
            mcmc_median_model = dup_model(bestfit_model)

        # Bootstrap
        if self.verbose:
            self.console.print(
                (
                    f"    [yellow]→ Running bootstrap "
                    f"({self.n_bootstrap} samples)...[/yellow]"
                ),
                end="",
            )
        if self.compute_bootstrap:
            fitmc = Ncm.FitMC.new(
                fit, Ncm.FitMCResampleType.BOOTSTRAP_NOMIX, Ncm.FitRunMsgs.NONE
            )
            fitmc.set_nthreads(12)

            # Set data file if prefix is provided
            if self.file_prefix is not None:
                bs_file = f"{self.file_prefix}_{np.exp(cut):.2f}_bs.fits"
                fitmc.set_data_file(bs_file)

            fitmc.start_run()
            fitmc.run(self.n_bootstrap)
            fitmc.end_run()
            fitmc.mean_covar()

            bootstrap_model = dup_model(cluster_m)
            model_params_from_list(
                bootstrap_model,
                fitmc.get_catalog().get_mean().dup_array(),
            )
        else:
            bootstrap_model = dup_model(bestfit_model)
        if self.verbose:
            self.console.print(" [green]✓[/green]")

        return CutAnalysisResult(
            cut=cut,
            n_clusters=n_clusters,
            bestfit=bestfit_model,
            mcmc_mean=mcmc_mean_model,
            mcmc_median=mcmc_median_model,
            bootstrap_mean=bootstrap_model,
            m2lnL=m2lnL,
        )

    def analyze(
        self, model_init: Nc.ClusterMassRichness
    ) -> dict[float, CutAnalysisResult]:
        """Analyze dataset with progressive cuts.

        Iterates through all cuts, performing the complete analysis pipeline on each.

        :param model_init: Initial model with starting parameters.
        :return: Results indexed by cut value
        """
        self._initialize_model(model_init)
        data = self._setup_real_data()

        if self.verbose:
            self.console.print(
                f"[cyan]Analyzing {len(self.data)} objects with "
                f"{len(self.cuts)} cuts...[/cyan]"
            )

        for cut in self.cuts:
            if self.verbose:
                self.console.print(f"  Processing cut {np.exp(cut):.2f}...")
            result = self._analyze_cut(data, cut)
            self.results[cut] = result
            if self.verbose:
                self.console.print(
                    f"    Found {result.n_clusters} objects. "
                    f"μ₀={result.bestfit['mup0']:.3f}"
                )

        # Display results table
        self.display_results(
            f"[bold cyan]{self.sample_desc}: Analysis Results[/bold cyan]"
        )
        return self.results

    def display_results(self, title: str = "Analysis Results"):
        """Display analysis results in a formatted table.

        :param title: Title for the results table
        """
        if not self.results:
            return

        # Get parameter names from first result
        first_result = next(iter(self.results.values()))
        param_names = get_model_param_names(first_result.bestfit)
        labels = param_names  # Use actual parameter names

        table = RichTable(title=title, show_header=True, header_style="bold cyan")
        table.add_column("Cut", style="green", width=4)
        table.add_column("N_clusters", style="cyan", width=7)
        table.add_column("Best-fit: " + ", ".join(labels), style="yellow")
        table.add_column("MCMC (mean): " + ", ".join(labels), style="magenta")
        table.add_column("Bootstrap: " + ", ".join(labels), style="blue")

        for cut in sorted(self.results.keys()):
            result = self.results[cut]
            table.add_row(
                f"{np.exp(cut):.1f}",
                str(result.n_clusters),
                ", ".join(
                    [
                        f"{val: {PARAM_FORMAT}}"
                        for val in model_params_as_list(result.bestfit)
                    ]
                ),
                ", ".join(
                    [
                        f"{val: {PARAM_FORMAT}}"
                        for val in model_params_as_list(result.mcmc_mean)
                    ]
                ),
                ", ".join(
                    [
                        f"{val: {PARAM_FORMAT}}"
                        for val in model_params_as_list(result.bootstrap_mean)
                    ]
                ),
            )

        self.console.print(Panel(table))

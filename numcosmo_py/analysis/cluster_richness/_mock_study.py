#
# _mock_study.py
#
# Copyright (C) 2024 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

"""MockStudy class for mock data generation and analysis.

This module provides tools for generating mock realizations from a fiducial
model and analyzing them to assess parameter biases and uncertainties.
"""

import numpy as np
from rich.console import Console, Group
from rich.live import Live
from rich.panel import Panel
from rich.table import Table as RichTable

from numcosmo_py import Nc, Ncm

from ._parameters import (
    CutAnalysisResult,
    dup_model,
    model_params_as_list,
    get_model_param_names,
)
from ._utils import PARAM_FORMAT
from ._analyzer import CutAnalyzer, COMPUTE_MCMC, COMPUTE_BOOTSTRAP
from ._database import BestfitDatabase


#: Rich console for output
console = Console()


class MockStudy:
    """Generate mocks from a fiducial model and analyze each with CutAnalyzer.

    This class orchestrates mock data generation and analysis:
    1. Generates mock realizations from fiducial parameters
    2. Analyzes each mock using CutAnalyzer
    3. Aggregates results across mocks

    Works with any NcClusterMassRichness subclass (Ascaso, Ext, etc.).
    """

    def __init__(
        self,
        model_fiducial: Nc.ClusterMassRichness,
        lnM: np.ndarray,
        z: np.ndarray,
        cuts: list[float],
        n_mocks: int = 100,
        n_bootstrap: int = 100,
        file_prefix: str | None = None,
        verbose: bool = True,
        fiducial_results: dict[float, CutAnalysisResult] | None = None,
        db_path: str = "bestfits.db",
        recompute: bool = False,
    ):
        """Initialize mock study.

        :param model_fiducial: Fiducial model to generate mocks from
        :param lnM: Log mass array (template)
        :param z: Redshift array (template)
        :param cuts: List of richness cuts to analyze
        :param n_mocks: Number of mocks to generate (default: 100)
        :param n_bootstrap: Number of bootstrap resamples per analysis (default: 100)
        :param file_prefix: If provided, pass seed-based prefix to CutAnalyzer for
            saving data (default: None)
        :param verbose: Whether to print progress (default: True)
        :param fiducial_results: Results from real data analysis for GOF comparison
            (default: None)
        :param db_path: Path to SQLite database for bestfit caching (default:
            "bestfits.db")
        :param recompute: Force recomputation of all bestfits, ignoring cached values
            (default: False)
        """
        self.model_fiducial = model_fiducial
        self.lnM = lnM
        self.z = z
        self.cuts = cuts
        self.n_mocks = n_mocks
        self.n_bootstrap = n_bootstrap
        self.file_prefix = file_prefix
        self.mock_results: list[dict[float, CutAnalysisResult]] = []
        self.fiducial_results = fiducial_results
        self.verbose = verbose
        self.db = BestfitDatabase(db_path)
        self.recompute = recompute
        self.cached_bestfits: dict[tuple[int, float], Nc.ClusterMassRichness] = {}

    def _generate_mock(
        self, mock_seed: int
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Generate a single mock realization.

        :param mock_seed: Random seed for reproducible generation
        :return: Tuple of (lnM, z, lnR_mock) for the generated mock
        """
        # Duplicate fiducial model to avoid modifying it
        cluster_m = dup_model(self.model_fiducial)
        cluster_m["cut"] = np.log(5.0)  # Same initial cut of the real data

        mset = Ncm.MSet.new_array([cluster_m])
        mset.prepare_fparam_map()

        rng = Ncm.RNG.seeded_new(None, int(mock_seed))
        data = Nc.DataClusterMassRich.new()

        lnM_v = Ncm.Vector.new_array(self.lnM)
        z_v = Ncm.Vector.new_array(self.z)
        rich = Ncm.Vector.new(len(self.lnM))
        rich.set_zero()

        data.set_data(lnM_v, z_v, rich)
        data.resample(mset, rng)

        # Extract resampled data
        lnR_mock = np.array(
            [data.peek_lnR().get(i) for i in range(data.peek_lnR().len())]
        )
        return self.lnM, self.z, lnR_mock

    def run(self, seed: int | None = None):
        """Run the mock study.

        Generates mocks and analyzes each one, storing results in mock_results.

        :param seed: Random seed for reproducibility (default: None)
        """
        rng = np.random.default_rng(seed)
        mock_seeds = [rng.integers(0, 2**31) for _ in range(self.n_mocks)]

        console.print(
            f"[magenta]Running mock study: {self.n_mocks} mocks "
            f"× {len(self.cuts)} cuts[/magenta]"
        )

        # Get parameter names for display
        param_names = get_model_param_names(self.model_fiducial)

        # Initialize tracking for running means
        running_bestfit: dict[float, list[Nc.ClusterMassRichness]] = {
            cut: [] for cut in self.cuts
        }
        running_mcmc: dict[float, list[Nc.ClusterMassRichness]] = {
            cut: [] for cut in self.cuts
        }
        running_mcmc_median: dict[float, list[Nc.ClusterMassRichness]] = {
            cut: [] for cut in self.cuts
        }
        running_bs: dict[float, list[Nc.ClusterMassRichness]] = {
            cut: [] for cut in self.cuts
        }
        running_m2lnL: dict[float, list[float]] = {cut: [] for cut in self.cuts}

        def create_gof_table():
            """Create GOF statistics table updated in real-time."""
            gof_table = RichTable(
                title="GOF: Fiducial m2lnL vs Mock Distribution",
                show_header=True,
                header_style="bold magenta",
            )
            gof_table.add_column("Cut", style="green", width=6)
            gof_table.add_column("Fiducial m2lnL", style="yellow", width=14)
            gof_table.add_column("Mock Mean ± Std", style="cyan", width=18)
            gof_table.add_column("p-value", style="blue", width=10)
            gof_table.add_column("z-score", style="magenta", width=10)
            gof_table.add_column("n_mocks", style="white", width=8)

            for cut in sorted(self.cuts):
                if self.fiducial_results is None or cut not in self.fiducial_results:
                    continue
                if len(running_m2lnL[cut]) == 0:
                    continue

                fiducial_m2lnL = self.fiducial_results[cut].m2lnL
                mock_m2lnLs = np.array(running_m2lnL[cut])
                mock_mean = np.mean(mock_m2lnLs)
                mock_std = np.std(mock_m2lnLs)
                p_value = np.sum(mock_m2lnLs >= fiducial_m2lnL) / len(mock_m2lnLs)
                z_score = (
                    (fiducial_m2lnL - mock_mean) / mock_std if mock_std > 0 else 0.0
                )

                gof_table.add_row(
                    f"{np.exp(cut):.1f}",
                    f"{fiducial_m2lnL:.2f}",
                    f"{mock_mean:.2f} ± {mock_std:.2f}",
                    f"{p_value:.4f}",
                    f"{z_score:+.2f}",
                    f"{len(mock_m2lnLs)}",
                )

            return gof_table

        def create_status_table():
            """Create the status table.

            Contains with latest and running mean values for all cuts.
            """
            table = RichTable(
                title="Mock Study Progress",
                show_header=True,
                header_style="bold cyan",
            )
            table.add_column("Cut", style="green", width=4)
            labels = param_names
            table.add_column("Latest: " + ", ".join(labels), style="white")
            table.add_column("Mean: " + ", ".join(labels), style="white")

            # Iterate through all cuts
            for cut in sorted(self.cuts):
                cut_label = f"{np.exp(cut):.1f}"

                if cut in running_bestfit and len(running_bestfit[cut]) > 0:
                    latest_bf = running_bestfit[cut][-1]
                    bf_list = [model_params_as_list(p) for p in running_bestfit[cut]]
                    mean_bf = np.mean(bf_list, axis=0)
                    std_bf = np.std(bf_list, axis=0)
                    table.add_row(
                        cut_label,
                        ", ".join(
                            [
                                f"{val: {PARAM_FORMAT}}"
                                for val in model_params_as_list(latest_bf)
                            ]
                        ),
                        ", ".join(
                            [
                                f"{val: {PARAM_FORMAT}} ± {err: {PARAM_FORMAT}}"
                                for val, err in zip(mean_bf, std_bf)
                            ]
                        ),
                    )

                if cut in running_mcmc and len(running_mcmc[cut]) > 0:
                    latest_mc = running_mcmc[cut][-1]
                    mc_list = [model_params_as_list(p) for p in running_mcmc[cut]]
                    mean_mc = np.mean(mc_list, axis=0)
                    std_mc = np.std(mc_list, axis=0)
                    table.add_row(
                        "",
                        ", ".join(
                            [
                                f"{val: {PARAM_FORMAT}}"
                                for val in model_params_as_list(latest_mc)
                            ]
                        ),
                        ", ".join(
                            [
                                f"{val: {PARAM_FORMAT}} ± {err: {PARAM_FORMAT}}"
                                for val, err in zip(mean_mc, std_mc)
                            ]
                        ),
                    )

                if cut in running_mcmc_median and len(running_mcmc_median[cut]) > 0:
                    latest_mc = running_mcmc_median[cut][-1]
                    mc_list = [
                        model_params_as_list(p) for p in running_mcmc_median[cut]
                    ]
                    mean_mc = np.mean(mc_list, axis=0)
                    std_mc = np.std(mc_list, axis=0)
                    table.add_row(
                        "",
                        ", ".join(
                            [
                                f"{val: {PARAM_FORMAT}}"
                                for val in model_params_as_list(latest_mc)
                            ]
                        ),
                        ", ".join(
                            [
                                f"{val: {PARAM_FORMAT}} ± {err: {PARAM_FORMAT}}"
                                for val, err in zip(mean_mc, std_mc)
                            ]
                        ),
                    )

                if cut in running_bs and len(running_bs[cut]) > 0:
                    latest_bs = running_bs[cut][-1]
                    bs_list = [model_params_as_list(p) for p in running_bs[cut]]
                    mean_bs = np.mean(bs_list, axis=0)
                    std_bs = np.std(bs_list, axis=0)
                    table.add_row(
                        "",
                        ", ".join(
                            [
                                f"{val: {PARAM_FORMAT}}"
                                for val in model_params_as_list(latest_bs)
                            ]
                        ),
                        ", ".join(
                            [
                                f"{val: {PARAM_FORMAT}} ± {err: {PARAM_FORMAT}}"
                                for val, err in zip(mean_bs, std_bs)
                            ]
                        ),
                    )

            return table

        with Live(console=console, refresh_per_second=1) as live:
            # Load existing bestfits from database if not recomputing
            if not self.recompute:
                all_computed_seeds = self.db.get_all_computed_seeds()
                for mock_seed in all_computed_seeds:
                    for cut in self.cuts:
                        bestfit = self.db.get_bestfit(mock_seed, cut)
                        if bestfit:
                            self.cached_bestfits[(mock_seed, cut)] = bestfit
                if all_computed_seeds:
                    console.print(
                        f"[cyan]Loaded {len(all_computed_seeds)} seeds from database "
                        f"({self.db.count_entries()} total entries)[/cyan]"
                    )

            for i, mock_seed_int in enumerate(mock_seeds):
                mock_seed = int(mock_seed_int)  # type: ignore
                # Generate mock
                lnM_mock, z_mock, lnR_mock = self._generate_mock(mock_seed)

                # Create seed-based prefix if file_prefix is provided
                mock_file_prefix = None
                if self.file_prefix is not None:
                    mock_file_prefix = f"{self.file_prefix}_seed{int(mock_seed)}"

                # Analyze with CutAnalyzer
                analyzer = CutAnalyzer(
                    lnM_mock,
                    z_mock,
                    lnR_mock,
                    self.cuts,
                    n_bootstrap=self.n_bootstrap,
                    compute_mcmc=COMPUTE_MCMC,
                    compute_bootstrap=COMPUTE_BOOTSTRAP,
                    file_prefix=mock_file_prefix,
                    sample_desc=f"Mock Seed: {int(mock_seed)}",
                    verbose=False,
                )
                results = analyzer.analyze(model_init=self.model_fiducial)
                self.mock_results.append(results)

                # Save bestfits to database with transaction wrapping
                for cut in self.cuts:
                    if cut in results:
                        self.db.insert_bestfit(mock_seed, results[cut])

                # Update running means for all cuts
                for cut in self.cuts:
                    if cut in results:
                        running_bestfit[cut].append(results[cut].bestfit)
                        running_mcmc[cut].append(results[cut].mcmc_mean)
                        running_mcmc_median[cut].append(results[cut].mcmc_median)
                        running_bs[cut].append(results[cut].bootstrap_mean)
                        running_m2lnL[cut].append(results[cut].m2lnL)

                # Create display content
                status_table = create_status_table()
                gof_table = create_gof_table()
                progress_text = f"[cyan]Processed {i + 1}/{self.n_mocks} mocks[/cyan]"

                # Combine tables for display
                display_content = Group(
                    Panel(status_table, subtitle=progress_text),
                    Panel(gof_table),
                )

                # Update live display
                live.update(display_content)

        console.print("[green]Mock study complete![/green]")

    def compute_gof_statistics(self) -> dict[float, dict[str, float]]:
        """Compute goodness-of-fit (GOF) statistics using parametric bootstrap.

        Compares the fiducial m2lnL values against the distribution from mocks using
        a parametric bootstrap with refit approach. For each cut, computes:
        - p-value: fraction of mocks with m2lnL >= fiducial m2lnL
        - z-score: standardized deviation from mock mean
        - Mean and std of mock m2lnL distribution

        :return: dictionary mapping cuts to GOF statistics
        """
        if not self.mock_results:
            console.print("[red]No mock results available for GOF computation[/red]")
            return {}

        if self.fiducial_results is None:
            console.print(
                "[red]Fiducial results not provided for GOF computation[/red]"
            )
            return {}

        gof_stats: dict[float, dict[str, float]] = {}

        for cut in sorted(self.cuts):
            if cut not in self.fiducial_results:
                continue

            fiducial_m2lnL = self.fiducial_results[cut].m2lnL
            mock_m2lnLs_list = [
                mock_result[cut].m2lnL
                for mock_result in self.mock_results
                if cut in mock_result
            ]

            if len(mock_m2lnLs_list) == 0:
                continue

            mock_m2lnLs = np.array(mock_m2lnLs_list)
            mock_mean = np.mean(mock_m2lnLs)
            mock_std = np.std(mock_m2lnLs)

            # Compute p-value: fraction of mocks with m2lnL >= fiducial
            p_value = np.sum(mock_m2lnLs >= fiducial_m2lnL) / len(mock_m2lnLs)

            # Compute z-score
            z_score = (fiducial_m2lnL - mock_mean) / mock_std if mock_std > 0 else 0.0

            gof_stats[cut] = {
                "fiducial_m2lnL": float(fiducial_m2lnL),
                "mock_mean": float(mock_mean),
                "mock_std": float(mock_std),
                "p_value": float(p_value),
                "z_score": float(z_score),
                "n_mocks": len(mock_m2lnLs),
            }

        return gof_stats

    def display_gof_results(self, gof_stats: dict[float, dict[str, float]]):
        """Display GOF statistics in a formatted panel.

        :param gof_stats: dictionary of GOF statistics from compute_gof_statistics
        """
        if not gof_stats:
            return

        table = RichTable(
            title="Goodness-of-Fit Analysis (Parametric Bootstrap with Refit)",
            show_header=True,
            header_style="bold magenta",
        )
        table.add_column("Cut", style="green", width=6)
        table.add_column("Fiducial m2lnL", style="yellow", width=14)
        table.add_column("Mock Mean ± Std", style="cyan", width=18)
        table.add_column("p-value", style="blue", width=10)
        table.add_column("z-score", style="magenta", width=10)
        table.add_column("n_mocks", style="white", width=8)

        for cut in sorted(gof_stats.keys()):
            stats = gof_stats[cut]
            table.add_row(
                f"{np.exp(cut):.1f}",
                f"{stats['fiducial_m2lnL']:.2f}",
                f"{stats['mock_mean']:.2f} ± {stats['mock_std']:.2f}",
                f"{stats['p_value']:.4f}",
                f"{stats['z_score']:+.2f}",
                f"{stats['n_mocks']}",
            )

        console.print(Panel(table))

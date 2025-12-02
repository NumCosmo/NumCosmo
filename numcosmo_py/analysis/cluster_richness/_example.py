#
# _example.py
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

"""Example usage of the cluster richness analysis package.

This module provides example functions demonstrating the complete
analysis workflow.
"""

import numpy as np
from astropy.table import Table
from rich.console import Console

from numcosmo_py import Nc

from ._parameters import model_params_as_list, get_model_param_names
from ._utils import PARAM_FORMAT
from ._analyzer import CutAnalyzer, COMPUTE_MCMC, COMPUTE_BOOTSTRAP
from ._mock_study import MockStudy


#: Rich console for output
console = Console()


def run_analysis_example(
    fits_file: str = "match_ID_v2.fits",
    n_mocks: int = 10000,
    n_bootstrap: int = 200,
    seed: int = 42,
    model_type: str = "ascaso",
) -> None:
    """Run the full analysis: real data -> mocks from best-fit.

    This function demonstrates the complete workflow:
    1. Load real cluster data from FITS file
    2. Analyze real data with MCMC to extract best-fit parameters
    3. Use best-fit as fiducial for mock generation
    4. Run mock study to assess parameter bias and uncertainties

    :param fits_file: Path to FITS file with cluster data (default: "match_ID_v2.fits")
    :param n_mocks: Number of mock realizations (default: 10000)
    :param n_bootstrap: Number of bootstrap samples per analysis (default: 200)
    :param seed: Random seed for reproducibility (default: 42)
    :param model_type: Type of richness model to use ("ascaso" or "ext")
        (default: "ascaso")
    """
    cuts = np.log(np.array([5, 10, 15, 20, 30]))

    # Load real data
    console.print("[bold cyan]Loading real data...[/bold cyan]")
    table_halos = Table.read(fits_file, hdu=1)
    table_halos["lnM"] = np.log(table_halos["halo_mass"])
    lnM = np.array(table_halos["lnM"])
    z = np.array(table_halos["redshift"])
    lnR = np.log(table_halos["richness"])
    console.print(f"  Loaded {len(lnM)} clusters")

    # Create initial model based on type
    model_init: Nc.ClusterMassRichness
    if model_type.lower() == "ascaso":
        model_init = Nc.ClusterMassAscaso(lnRichness_min=0.0, lnRichness_max=20.0)
        model_init["mup0"] = 4.0
        model_init["mup1"] = 1.0
        model_init["mup2"] = 0.2
        model_init["mup3"] = 0.0
        model_init["sigmap0"] = 0.5
        model_init["sigmap1"] = 0.03
        model_init["sigmap2"] = 0.15
    elif model_type.lower() == "ext":
        model_init = Nc.ClusterMassExt(lnRichness_min=0.0, lnRichness_max=20.0)
        model_init["mup0"] = 4.0
        model_init["mup1"] = 1.0
        model_init["mup2"] = 0.1
        model_init["mup3"] = 0.01
        model_init["sigmap0"] = -0.3
        model_init["sigmap1"] = -0.08
        model_init["sigmap2"] = 0.005
    else:
        raise ValueError(f"Unknown model type: {model_type}. Use 'ascaso' or 'ext'.")

    # Phase 1: Analyze real data
    console.print("\n[bold magenta]PHASE 1: Real Data Analysis[/bold magenta]")
    real_analyzer = CutAnalyzer(
        lnM,
        z,
        lnR,
        cuts.tolist(),
        n_bootstrap=n_bootstrap,
        compute_mcmc=COMPUTE_MCMC,
        compute_bootstrap=COMPUTE_BOOTSTRAP,
        file_prefix="real",
    )
    real_results = real_analyzer.analyze(model_init=model_init)

    # Extract best-fit from first cut as fiducial for mocks
    model_fiducial = real_results[cuts[0]].bestfit
    param_names = get_model_param_names(model_fiducial)
    param_values = model_params_as_list(model_fiducial)

    console.print(
        f"\n[green]Using best-fit from cut {np.exp(cuts[0]):.1f} as fiducial:[/green]"
    )
    for name, value in zip(param_names, param_values):
        console.print(f"  {name}={value:{PARAM_FORMAT}}")

    # Phase 2: Mock analysis from best-fit
    console.print("\n[bold magenta]PHASE 2: Mock Study[/bold magenta]")
    mock_study = MockStudy(
        model_fiducial,
        lnM,
        z,
        cuts.tolist(),
        n_mocks=n_mocks,
        n_bootstrap=n_bootstrap,
        file_prefix="mock",
        fiducial_results=real_results,
    )
    mock_study.run(seed=seed)

    # Phase 3: GOF Analysis
    console.print("\n[bold magenta]PHASE 3: Goodness-of-Fit Analysis[/bold magenta]")
    gof_stats = mock_study.compute_gof_statistics()
    mock_study.display_gof_results(gof_stats)

    console.print("\n[bold green]Analysis complete![/bold green]")
    console.print(f"\n[cyan]Real data: {len(real_results)} cuts analyzed[/cyan]")
    console.print(
        f"[cyan]Mock study: {len(mock_study.mock_results)} mocks analyzed[/cyan]"
    )


if __name__ == "__main__":
    run_analysis_example()

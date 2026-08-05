#
# catalog.py
#
# Wed Feb 14 18:59:34 2024
# Copyright  2024  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# catalog.py
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

"""NumCosmo APP subcommands to analyze catalogs."""

import math
import dataclasses
from typing import Optional, Annotated, List
from pathlib import Path
import typer
from rich.table import Table
from rich.text import Text
from rich.progress import track
import numpy as np

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

from .. import Ncm
from ..interpolation.stats_dist import (
    create_stats_dist,
    CrossValidationMethod,
    InterpolationKernel,
    InterpolationMethod,
)
from ..plotting.tools import set_rc_params_article, confidence_ellipse
from ..safe_eval import compile_expr, SafeExprError
from ..catalog_stats import (
    DerivedStat,
    parse_variable_bindings,
    resolve_param,
    stat_center_and_bounds,
)
from .loading import LoadCatalog, LoadExperiment, LoadedCatalog, load_catalog
from .logging import AppLogging
from ..plotting import mcat_to_catalog_data, plot_mcsamples
from ..plotting.derived import add_derived_column


@dataclasses.dataclass(kw_only=True)
class AnalyzeMCMC(LoadCatalog):
    """Analyzes the results of a MCMC run."""

    evidence: Annotated[
        bool,
        typer.Option(
            help=(
                "Computes the ln-Bayesian evidence and the 1sigma parameter "
                "space ln-volume."
            ),
        ),
    ] = False

    def __post_init__(self) -> None:
        """Analyzes the results of a MCMC run."""
        super().__post_init__()

        mcat = self.mcat
        fs = self.full_stats
        if self.nitems >= 10:
            mcat.estimate_autocorrelation_tau(False)

        desc_color = "bold bright_cyan"
        values_color = "bold bright_green"
        main_table = Table(title="Catalog information")
        main_table.show_header = False

        main_table.add_column(justify="left")

        details = Table(title="Run details", expand=False)
        details.show_header = False
        details.add_column(justify="left", style=desc_color)
        details.add_column(justify="right", style=values_color)

        details.add_row("Run type", mcat.get_run_type())
        details.add_row("Size", f"{mcat.len()}")
        details.add_row("Number of Iterations", f"{mcat.max_time()}")
        details.add_row("Number of chains", f"{self.nchains}")
        details.add_row("Number of parameters", f"{self.fparams_len}")
        details.add_row("Number of extra columns", f"{self.nadd_vals}")
        details.add_row("Weighted", f"{mcat.weighted()}")
        main_table.add_row(details)

        if self.nitems == 0:
            self.console.print(main_table)
            self.console.print("#  Empty catalog!")

            self.close_logging()
            return

        # Global diagnostics

        global_diag = Table(
            title="Global Convergence Diagnostics",
            expand=False,
        )
        global_diag.add_column("Diagnostic Statistic", justify="left", style=desc_color)
        global_diag.add_column("Suggested cut-off", justify="left", style=values_color)
        global_diag.add_column("Worst parameter", justify="left", style=values_color)
        global_diag.add_column("AR model order", justify="left", style=values_color)
        global_diag.add_column("Value", justify="left", style=values_color)

        param_diag = Table(title="Parameters", expand=False, show_lines=True)
        param_diag_matrix = []

        # Parameter names
        param_diag.add_column(
            "Parameter", justify="left", style=desc_color, vertical="middle"
        )
        param_diag_matrix.append([mcat.col_full_name(i) for i in self.indices])

        # Values color
        val_color = values_color
        # Parameter best fit
        best_fit_vec = mcat.get_bestfit_row()
        if best_fit_vec is not None:
            param_diag.add_column(
                "Best-fit", justify="left", style=val_color, vertical="middle"
            )
            param_diag_matrix.append(
                [f"{best_fit_vec.get(i): .6g}" for i in self.indices]
            )

        # Parameter mean
        param_diag.add_column(
            "Mean", justify="left", style=val_color, vertical="middle"
        )
        param_diag_matrix.append([f"{fs.get_mean(i): .6g}" for i in self.indices])

        # Standard Deviation

        param_diag.add_column(
            "Standard Deviation", justify="left", style=val_color, vertical="middle"
        )
        param_diag_matrix.append([f"{fs.get_sd(i): .6g}" for i in self.indices])

        if self.nitems >= 10:
            # Mean Standard Deviation
            param_diag.add_column(
                "Mean Standard Deviation",
                justify="left",
                style=val_color,
                vertical="middle",
            )
            tau_vec = mcat.peek_autocorrelation_tau()

            mean_sd_array = [
                np.sqrt(fs.get_var(i) * tau_vec.get(i) / fs.nitens())
                for i in self.indices
            ]
            param_diag_matrix.append([f"{mean_sd: .6g}" for mean_sd in mean_sd_array])

            # Autocorrelation Time
            tau_row = []
            tau_row.append("Autocorrelation time (tau)")
            tau_row.append("NA")
            tau_row.append(
                f"{tau_vec.get_max():.0f} "
                f"({mcat.col_full_name(tau_vec.get_max_index())})"
            )
            tau_row.append("NA")
            tau_row.append(f"{tau_vec.get_max():.3f}")
            global_diag.add_row(*tau_row)

            param_diag.add_column(
                "tau", justify="left", style=val_color, vertical="middle"
            )
            param_diag_matrix.append([f"{tau_vec.get(i): .6g}" for i in self.indices])

        if self.nchains > 1:
            # Gelman Rubin
            gelman_rubin_row = []
            gelman_rubin_row.append("Gelman-Rubin (G&B) Shrink Factor (R-1)")
            skf = [mcat.get_param_shrink_factor(i) - 1 for i in self.indices]
            gelman_rubin_row.append("NA")
            gr_worst = int(np.argmin(skf))
            gelman_rubin_row.append(
                f"{skf[gr_worst]:.3f} ({mcat.col_full_name(gr_worst)})"
            )
            gelman_rubin_row.append("NA")
            gelman_rubin_row.append(f"{mcat.get_shrink_factor() - 1:.3f}")
            global_diag.add_row(*gelman_rubin_row)

            param_diag.add_column(
                "G&R", justify="left", style=val_color, vertical="middle"
            )
            param_diag_matrix.append([f"{skf_i:.3f}" for skf_i in skf])

        # Constant Break

        cb = [self.stats.estimate_const_break(i) for i in self.indices]
        cb_worst = int(np.argmax(cb))
        const_break_row = []
        const_break_row.append("Constant Break (CB) (iterations, points)")
        const_break_row.append(f"{cb[cb_worst]:.0f}")
        const_break_row.append(f"{cb[cb_worst]:.0f} ({mcat.col_full_name(cb_worst)})")
        const_break_row.append("NA")
        const_break_row.append(f"{cb[cb_worst]:.0f}")
        global_diag.add_row(*const_break_row)

        param_diag.add_column("CB", justify="left", style=val_color)
        param_diag_matrix.append(
            [f"{cb_i:.0f} {cb_i * self.nchains:.0f}" for cb_i in cb]
        )

        if self.nitems >= 10:
            # Effective sample size
            (
                ess_vec,
                ess_best_cutoff,
                ess_worst_index,
                ess_worst_order,
                ess_worst_ess,
            ) = self.stats.max_ess_time(100)
            ess_row = []
            ess_row.append("Effective Sample Size (ESS) (ensembles, points)")
            ess_row.append(f"{ess_best_cutoff}")
            ess_row.append(
                f"{ess_vec.get(ess_worst_index):.0f} "
                f"({mcat.col_full_name(ess_worst_index)})"
            )
            ess_row.append(f"{ess_worst_order}")
            ess_row.append(f"{ess_worst_ess:.0f}")
            global_diag.add_row(*ess_row)

            param_diag.add_column("ESS", justify="left", style=val_color)
            param_diag_matrix.append(
                [
                    f"{ess_vec.get(i):.0f} {ess_vec.get(i) * self.nchains:.0f}"
                    for i in self.indices
                ]
            )

            # Heidelberger and Welch

            hw_pvalue = 1.0 - 0.95 ** (1.0 / self.fparams_len)
            (
                hw_vec,
                hw_best_cutoff,
                hw_worst_index,
                hw_worst_order,
                hw_worst_pvalue,
            ) = self.stats.heidel_diag(100, hw_pvalue)

            hw_row = []
            hw_row.append(f"Heidelberger and Welch p-value (>{hw_pvalue * 100.0:.1f}%)")

            if hw_best_cutoff >= 0:
                hw_row.append(f"{hw_best_cutoff}")
            else:
                hw_row.append("All parameters fail")
            hw_row.append(
                f"{(1.0 - hw_worst_pvalue) * 100.0:.1f}% "
                f"({mcat.col_full_name(hw_worst_index)})"
            )
            hw_row.append(f"{hw_worst_order}")
            hw_row.append(f"{(1.0 - hw_worst_pvalue) * 100.0:.1f}%")
            global_diag.add_row(*hw_row)

            param_diag.add_column(
                "H&W",
                justify="left",
                style=val_color,
            )
            param_diag_matrix.append(
                [f"{(1.0 - hw_vec.get(i)) * 100.0:.1f}" for i in self.indices]
            )

        for row in np.array(param_diag_matrix).T:
            param_diag.add_row(*row)

        # Add the global diagnostics to the main table
        main_table.add_row(global_diag)
        main_table.add_row(param_diag)

        covariance_matrix = Table(title="Covariance Matrix", expand=False)
        covariance_matrix.add_column("Parameter", justify="right", style="bold")
        for i in self.indices:
            covariance_matrix.add_column(
                mcat.col_name(i).split(":")[-1], justify="right"
            )

        for i in self.indices:
            row = [mcat.col_name(i).split(":")[-1]]
            for j in self.indices:
                cov_ij = fs.get_cor(i, j)
                cor_ij_string = f"{cov_ij * 100.0: 3.0f}%"
                styles_array = [
                    "bold bright_red",
                    "bright_red",
                    "dim bright_red",
                    "dim bright_green",
                    "bright_green",
                    "bold bright_green",
                ]
                cov_color_index = int(np.round((cov_ij + 1.0) * 2.5))

                row.append(Text(cor_ij_string, style=styles_array[cov_color_index]))
            covariance_matrix.add_row(*row)

        main_table.add_row(covariance_matrix)

        if self.evidence:
            evidence_table = Table(title="Posterior Analysis", expand=False)
            evidence_table.add_column("Evidence Type", justify="left", style=desc_color)
            evidence_table.add_column(
                "Value +/- 1-sigma", justify="left", style=values_color
            )
            evidence_table.add_column(
                "1-sigma volume", justify="left", style=values_color
            )
            evidence_table.add_column(
                "2-sigma volume", justify="left", style=values_color
            )

            be, be_sd = mcat.get_post_lnnorm()
            lnevol_1s, glnvol_1s = mcat.get_post_lnvol(0.682689492137086)
            lnevol_2s, glnvol_2s = mcat.get_post_lnvol(0.954499736103642)

            evidence_table.add_row(
                "Bayesian ln-Evidence", f"{be:.5g} +/- {be_sd:.5g}", "--", "--"
            )

            evidence_table.add_row(
                "Posterior ln-volume",
                "--",
                f"{lnevol_1s:.5g}",
                f"{lnevol_2s:.5g}",
            )
            evidence_table.add_row(
                "Posterior ln-volume (Gaussian approx.)",
                "--",
                f"{glnvol_1s:.5g}",
                f"{glnvol_2s:.5g}",
            )

            main_table.add_row(evidence_table)

        self.console.print(main_table)

        self.close_logging()


@dataclasses.dataclass(kw_only=True)
class CalibrateCatalog(LoadCatalog):
    """Calibrate the APES sampler using a given catalog."""

    robust: Annotated[
        bool,
        typer.Option(
            help="Use robust covariance estimation.",
        ),
    ] = False

    interpolation_method: Annotated[
        InterpolationMethod,
        typer.Option(
            help="Interpolation method to use.",
        ),
    ] = InterpolationMethod.VKDE

    interpolation_kernel: Annotated[
        InterpolationKernel,
        typer.Option(
            help="Interpolation kernel to use.",
        ),
    ] = InterpolationKernel.CAUCHY

    cv_method: Annotated[
        CrossValidationMethod,
        typer.Option(
            help=(
                "Cross-validation method to use. If NONE, no cross-validation is "
                "used and only weights information is printed. If SPLIT, the sample "
                "is split into two parts, one for training and the other for testing. "
                "If SPLIT_NOFIT, the sample is split into two parts, one for training "
                "and the other for testing, but equal weights are used for both parts."
            ),
        ),
    ] = CrossValidationMethod.SPLIT_NOFIT

    over_smooth: Annotated[
        float,
        typer.Option(
            help="Over-smoothing parameter to use.",
            min=1.0e-2,
        ),
    ] = 1.0

    split_fraction: Annotated[
        Optional[float],
        typer.Option(
            help="Split fraction to use.",
            min=0.02,
        ),
    ] = None

    local_fraction: Annotated[
        Optional[float],
        typer.Option(
            help="Local fraction to use.",
            min=0.02,
        ),
    ] = None

    interpolate: Annotated[
        bool,
        typer.Option(
            help="Use interpolation to compute the weights of the APES approximation.",
        ),
    ] = True

    ntries: Annotated[
        int,
        typer.Option(
            help="Number of tries to sample from the calibrated distribution.",
            min=1,
        ),
    ] = 100

    use_half: Annotated[
        bool,
        typer.Option(
            help="Use half of the walkers to calibrate the sampler.",
        ),
    ] = True

    verbose: Annotated[
        bool,
        typer.Option(
            help="Prints verbose information.",
        ),
    ] = False

    plot_2d: Annotated[
        bool,
        typer.Option(
            help="Plots 2D confidence ellipses.",
        ),
    ] = False

    plot_2d_scatter: Annotated[
        bool,
        typer.Option(
            help="Plots 2D scatter plots.",
        ),
    ] = False

    def __post_init__(self) -> None:
        """Calibrate the APES sampler using a given catalog."""
        super().__post_init__()

        mcat = self.mcat
        m2lnL_id = mcat.get_m2lnp_var()  # pylint: disable-msg=invalid-name
        mcat_len = mcat.len()

        nwalkers = self.nchains
        if self.use_half:
            nwalkers = nwalkers // 2

        last_e = [mcat.peek_row(mcat_len - nwalkers + i) for i in range(nwalkers)]
        ncols = mcat.ncols()
        nvar = ncols - self.nadd_vals
        params = [
            "$" + mcat.col_symb(i) + "$" for i in range(self.nadd_vals, mcat.ncols())
        ]

        sdist = create_stats_dist(
            robust=self.robust,
            interpolation_method=self.interpolation_method,
            interpolation_kernel=self.interpolation_kernel,
            cv_method=self.cv_method,
            dim=nvar,
            over_smooth=math.fabs(self.over_smooth),
            split_fraction=self.split_fraction,
            local_fraction=self.local_fraction,
            verbose=self.verbose,
        )

        m2lnL = []
        for row in last_e:
            m2lnL.append(row.get(m2lnL_id))
            sdist.add_obs(row.get_subvector(self.nadd_vals, nvar))

        m2lnL_v = Ncm.Vector.new_array(m2lnL)
        if self.interpolate:
            sdist.prepare_interp(m2lnL_v)
        else:
            sdist.prepare()

        ovs = sdist.get_over_smooth()

        main_table = Table(
            title="Catalog calibration information",
            caption=(
                "APES approximation of the posterior distribution. The calibration "
                "information shows how well the APES approximation fits the last "
                "nwalkers sample of the MCMC chain. Too concentrated weights "
                "indicate that the APES approximation is not a good fit."
            ),
            min_width=88,
        )
        main_table.show_header = False

        main_table.add_column(justify="left", style="bold bright_cyan")
        main_table.add_column(justify="right", style="bold bright_green")

        main_table.add_row("Number of walkers", f"{nwalkers}")
        main_table.add_row("Number of parameters", f"{nvar}")
        main_table.add_row("Number of tries", f"{self.ntries}")
        main_table.add_row("Over-smoothing parameter", f"{ovs:.2f}")
        main_table.add_row("Interpolation method", f"{self.interpolation_method.value}")
        main_table.add_row("Interpolation kernel", f"{self.interpolation_kernel.value}")
        main_table.add_row("Cross-validation method", f"{self.cv_method.value}")
        main_table.add_row("Split fraction", f"{self.split_fraction}")
        main_table.add_row("Local fraction", f"{self.local_fraction}")
        main_table.add_row("Use interpolation", f"{self.interpolate}")
        main_table.add_row("Use half of the walkers", f"{self.use_half}")

        rng = Ncm.RNG.new()
        var_vector = Ncm.Vector.new(nvar)

        try_sample_array = []
        for _ in range(self.ntries):
            sdist.sample(var_vector, rng)
            try_sample_array.append(var_vector.dup_array())

        try_sample = np.array(try_sample_array)

        weights = np.array(sdist.peek_weights().dup_array())
        weights = weights / np.sum(weights)
        max_w = np.max(weights[np.nonzero(weights)])
        min_w = np.min(weights[np.nonzero(weights)])

        main_table.add_row("Min weight", f"{min_w:.2e}")
        main_table.add_row("Max weight", f"{max_w:.2e}")
        main_table.add_row(
            "Mean weight", f"{np.mean(weights[np.nonzero(weights)]):.2e}"
        )
        main_table.add_row(
            "Weight standard deviation", f"{np.std(weights[np.nonzero(weights)]):.2e}"
        )
        main_table.add_row(
            "Median weight", f"{np.median(weights[np.nonzero(weights)]):.2e}"
        )
        main_table.add_row("Non-zero weights", f"{np.count_nonzero(weights)}")
        main_table.add_row("Final bandwidth", f"{sdist.get_href():.2f}")

        self.console.print(main_table)

        if self.plot_2d:
            for a in range(nvar):  # pylint: disable-msg=invalid-name
                for b in range(a + 1, nvar):  # pylint: disable-msg=invalid-name
                    indices = np.array([a, b])
                    print(f"# {indices}")

                    _, axis = plt.subplots(1, 1, figsize=(16, 8))

                    # pylint: disable-next=invalid-name
                    for ii in range(0, int(sdist.get_n_kernels())):
                        y_i, cov_i, _, w_i = sdist.get_Ki(ii)
                        mean = np.array(y_i.dup_array())
                        cov = np.array(
                            [[cov_i.get(i, j) for j in indices] for i in indices]
                        )
                        cov = cov * 1.0
                        w_i = weights[ii]

                        if w_i > 0.0:
                            confidence_ellipse(
                                mean[indices],
                                cov,
                                axis,
                                edgecolor="red",
                                facecolor="red",
                            )
                    if self.plot_2d_scatter:
                        axis.scatter(try_sample[:, a], try_sample[:, b])
                    plt.axis("auto")
                    plt.xlabel(params[a])
                    plt.ylabel(params[b])
                    plt.grid()
                    plt.show()


@dataclasses.dataclass(kw_only=True)
class PlotCorner(AppLogging):
    """Plots the corner plot of one or more catalogs.

    Each catalog file is self-sufficient (model-set and, if used,
    --derived-* parameter resolution all come from the file itself), so no
    experiment file is needed. Give more than one catalog to overlay them
    in the same corner plot, e.g. to compare two runs.
    """

    mcmc_file: Annotated[
        List[Path],
        typer.Argument(
            help=(
                "Path(s) to the MCMC catalog file(s). Give more than one to "
                "overlay them in the same corner plot."
            ),
        ),
    ]

    burnin: Annotated[
        int,
        typer.Option(
            help=(
                "Number of iterations (ensemble steps) to discard as burnin, "
                "applied to every catalog."
            ),
            min=0,
        ),
    ] = 0

    tail: Annotated[
        Optional[int],
        typer.Option(
            help=(
                "Keep only the last N iterations (ensemble steps) of every "
                "catalog instead of burning in from the start. Give at most "
                "one of --burnin/--tail."
            ),
            min=0,
        ),
    ] = None

    include: Annotated[
        Optional[List[str]],
        typer.Option(
            help="List of parameters and or model names to include in the analysis.",
        ),
    ] = None

    exclude: Annotated[
        Optional[List[str]],
        typer.Option(
            help="List of parameters and or model names to exclude from the analysis.",
        ),
    ] = None

    output: Annotated[
        Optional[Path],
        typer.Option(
            "--output",
            "-o",
            help="Path to the output file, if given.",
        ),
    ] = None

    plot_name: Annotated[
        Optional[List[str]],
        typer.Option(
            help=(
                "Legend name(s) for each catalog, matched by position. May be "
                "given fewer times than there are catalogs; missing entries "
                "default to the corresponding file's stem."
            ),
        ),
    ] = None

    remove_index: Annotated[
        Optional[List[int]],
        typer.Option(
            help="Index of the parameter to remove.",
            min=0,
        ),
    ] = None

    mark_bestfit: Annotated[
        bool,
        typer.Option(
            help="Mark the best-fit parameters of the first catalog.",
        ),
    ] = False

    title_limit: Annotated[
        int,
        typer.Option(
            help="Include the n-sigma limit in the title.",
            min=0,
        ),
    ] = 0

    auto_thin: Annotated[
        bool,
        typer.Option(
            help=(
                "Automatically thin the data using the estimate of "
                "the autocorrelation time."
            ),
        ),
    ] = True

    show: Annotated[
        bool,
        typer.Option(
            help="Show the figure.",
        ),
    ] = True

    derived_variable: Annotated[
        Optional[List[str]],
        typer.Option(
            "--derived-variable",
            help=(
                "Bind a name usable in --derived-expr to a catalog parameter "
                "or add-value, as name=parameter (or just parameter, short "
                "for parameter=parameter). May be given multiple times. "
                "Requires --derived-expr."
            ),
        ),
    ] = None

    derived_expr: Annotated[
        Optional[str],
        typer.Option(
            help=(
                "Add an extra corner-plot dimension computed from this "
                "expression of the --derived-variable bindings, e.g. "
                "'10**x'. See `catalog derived-error --help` for the "
                "supported syntax."
            ),
        ),
    ] = None

    derived_symbol: Annotated[
        Optional[str],
        typer.Option(
            help=(
                "Axis label for --derived-expr. Defaults to the expression " "itself."
            ),
        ),
    ] = None

    derived_name: Annotated[
        str,
        typer.Option(
            help="Internal column name for --derived-expr.",
        ),
    ] = "derived"

    def __post_init__(self) -> None:
        """Corner plot of one or more catalogs."""
        super().__post_init__()

        if not self.mcmc_file:
            raise typer.BadParameter("At least one MCMC catalog file is required.")

        derived_variable: List[str] = []
        if self.derived_expr is not None:
            if not self.derived_variable:
                raise typer.BadParameter(
                    "--derived-expr requires at least one --derived-variable."
                )
            derived_variable = self.derived_variable

        plot_names = self.plot_name or []
        if len(plot_names) > len(self.mcmc_file):
            raise typer.BadParameter("More --plot-name values than catalog files.")

        mcsamples = []
        bestfit = None
        for i, mcmc_file in enumerate(self.mcmc_file):
            name = plot_names[i] if i < len(plot_names) else mcmc_file.stem
            loaded = load_catalog(
                mcmc_file, self.burnin, self.tail, self.include, self.exclude
            )

            thin = 1
            if self.auto_thin:
                loaded.mcat.estimate_autocorrelation_tau(False)
                thin = int(np.ceil(loaded.mcat.peek_autocorrelation_tau().get_max()))

            cd = mcat_to_catalog_data(
                loaded.mcat, name, indices=loaded.indices, thin=thin
            )
            if self.derived_expr is not None:
                cd = add_derived_column(
                    cd,
                    loaded.mset,
                    loaded.nadd_vals,
                    loaded.mcat,
                    derived_variable,
                    self.derived_expr,
                    self.derived_symbol,
                    self.derived_name,
                )

            if i == 0 and self.mark_bestfit:
                bestfit = cd.bestfit

            mcsamples.append(cd.to_mcsamples(collapse=True))

        fig = plot_mcsamples(mcsamples, markers=bestfit, title_limit=self.title_limit)
        if self.output is not None:
            filename = self.output.with_suffix(".corner.pdf").absolute().as_posix()
            fig.savefig(filename, bbox_inches="tight", dpi=300)

        if self.show:
            plt.show()

        self.close_logging()


@dataclasses.dataclass(kw_only=True)
class ParameterAnalysis(LoadCatalog):
    """Plots the corner plot of the catalog."""

    plot_name: Annotated[
        Optional[str],
        typer.Option(help="Name of the plot file."),
    ] = None

    param_name: Annotated[
        str,
        typer.Option(help="Name of the plot file."),
    ]

    def __post_init__(self) -> None:
        """Parameter analysis."""
        super().__post_init__()

        self.pi, self.pindex = resolve_param(self.mset, self.nadd_vals, self.param_name)
        self.symbol: str = self.mcat.col_symb(self.pindex)


@dataclasses.dataclass(kw_only=True)
class VisualHW(ParameterAnalysis):
    """Visual Heidelberger and Welch."""

    def __post_init__(self) -> None:
        """Visual Heidelberger and Welch."""
        super().__post_init__()

        cumsum_vec, mean, var = self.stats.visual_heidel_diag(self.pindex, 0)
        cumsum = np.array(cumsum_vec.dup_array())
        mean_a = (np.array(range(len(cumsum))) + 1.0) * mean

        set_rc_params_article(ncol=2)
        _, ax = plt.subplots()

        ax.title.set_text("Visual Heidelberger and Welch.")
        ax.plot(cumsum, label=f"Cumulative sum -- ${self.symbol}$")
        ax.plot(mean_a, label=f"Mean, standard deviation = {np.sqrt(var):.2f}")

        ax.set_xlabel("Iterations")
        ax.set_ylabel("Cumulative sum")
        ax.legend(loc="best")
        if self.plot_name is not None:
            plt.savefig(self.plot_name, bbox_inches="tight")

        plt.show()

        self.close_logging()


@dataclasses.dataclass(kw_only=True)
class ParameterEvolution(ParameterAnalysis):
    """Parameter evolution."""

    grid_size: Annotated[
        int,
        typer.Option(help="Grid size."),
    ] = 200

    def __post_init__(self) -> None:
        """Parameter evolution."""
        super().__post_init__()

        mcat = self.mcat
        if self.pi is not None:
            param_vec, evol_matrix = mcat.calc_param_ensemble_evol(
                self.pi, self.grid_size, Ncm.FitRunMsgs.NONE
            )
        else:
            param_vec, evol_matrix = mcat.calc_add_param_ensemble_evol(
                self.pindex, self.grid_size, Ncm.FitRunMsgs.NONE
            )

        param = np.array(param_vec.dup_array())
        evol_a = np.abs(evol_matrix.dup_array())
        min_evol_a = min(evol_a[evol_a > 0.0])
        evol_a[evol_a == 0.0] = min_evol_a
        evol = evol_a.reshape((-1, self.grid_size))

        set_rc_params_article(ncol=2)
        _, ax = plt.subplots()

        ax.title.set_text(f"Parameter evolution -- ${self.symbol}$")
        cax = ax.matshow(
            evol.T,
            aspect="auto",
            origin="lower",
            cmap="viridis",
            extent=[0, evol.shape[0], param[0], param[-1]],
            norm=LogNorm(),
        )
        ax.set_xlabel("Iterations")
        ax.set_ylabel(f"${self.symbol}$")
        plt.colorbar(cax, label=rf"$p_t\left({self.symbol}\right)$")

        if self.plot_name is not None:
            plt.savefig(self.plot_name, bbox_inches="tight")
        plt.show()

        self.close_logging()


@dataclasses.dataclass(kw_only=True)
class DerivedQuantityError(LoadCatalog):
    """Posterior statistic and asymmetric error bars for derived quantities.

    Builds the posterior distribution of one or more arbitrary expressions
    of catalog parameters (e.g. ``(y / 100)**2 * x``), sharing the same
    --variable bindings, by evaluating each expression on every sample in
    the catalog in a single pass. Reports each quantity's median, mode,
    and/or best-fit value together with 1, 2 and 3-sigma asymmetric error
    bars in one table.
    """

    variable: Annotated[
        List[str],
        typer.Option(
            "--variable",
            "-x",
            help=(
                "Bind a name usable in --expr to a catalog parameter or "
                "add-value, as name=parameter (or just parameter, short for "
                "parameter=parameter). May be given multiple times, e.g. "
                "-x x=Omega_m -x y=H0, or -x log10MDelta."
            ),
        ),
    ]

    expr: Annotated[
        List[str],
        typer.Option(
            "--expr",
            help=(
                "Mathematical expression combining the bound variables, e.g. "
                "'(y / 100)**2 * x'. May be given multiple times to report "
                "several derived quantities in one table, all sharing the "
                "same --variable bindings. Supports + - * / // % ** and "
                "unary minus, the constants pi/e, and the functions log, "
                "log10, log2, exp, sqrt, abs, sin, cos, tan, arcsin, arccos, "
                "arctan, sinh, cosh, tanh."
            ),
        ),
    ]

    stat: Annotated[
        Optional[List[DerivedStat]],
        typer.Option(
            help=(
                "Statistic(s) to report. Unlike median/bestfit, mode is NOT "
                "invariant under a nonlinear --expr: mode(g(x)) != g(mode(x)) "
                "in general, and can land far from g(median(x))/g(bestfit(x)) "
                "for a strongly nonlinear g (e.g. 10**x) near a prior bound."
            ),
        ),
    ] = None

    symbol: Annotated[
        Optional[List[str]],
        typer.Option(
            "--symbol",
            help=(
                "Display symbol(s) for each --expr, matched by position. "
                "May be given fewer times than --expr; missing entries "
                "default to the corresponding --expr string itself."
            ),
        ),
    ] = None

    def __post_init__(self) -> None:
        """Compute the posterior statistic(s) for the derived quantities."""
        super().__post_init__()

        if not self.variable:
            raise typer.BadParameter("At least one --variable binding is required.")
        if not self.expr:
            raise typer.BadParameter("At least one --expr is required.")
        if not self.stat:
            self.stat = [DerivedStat.MEDIAN]
        if self.symbol and len(self.symbol) > len(self.expr):
            raise typer.BadParameter("More --symbol values than --expr values.")

        symbols = self.symbol or []
        displays = [
            symbols[i] if i < len(symbols) else self.expr[i]
            for i in range(len(self.expr))
        ]

        try:
            var_pindex = parse_variable_bindings(
                self.mset, self.nadd_vals, self.variable, "--variable"
            )
        except ValueError as exc:
            raise typer.BadParameter(str(exc)) from exc

        try:
            evaluators = [compile_expr(e, var_pindex.keys()) for e in self.expr]
        except SafeExprError as exc:
            raise typer.BadParameter(str(exc)) from exc

        mcat = self.mcat
        epdfs = [Ncm.StatsDist1dEPDF.new(1.0e-3) for _ in self.expr]
        for i in range(mcat.len()):
            row = mcat.peek_row(i)
            values = {name: row.get(pindex) for name, pindex in var_pindex.items()}
            for epdf, evaluate in zip(epdfs, evaluators):
                epdf.add_obs(evaluate(values))
        for epdf in epdfs:
            epdf.prepare()

        bestfit_row = mcat.get_bestfit_row()
        bf_values = {
            name: bestfit_row.get(pindex) for name, pindex in var_pindex.items()
        }

        bindings_desc = ", ".join(self.variable)
        table = Table(title=f"Derived quantities (where {bindings_desc})", expand=False)
        table.add_column("Quantity", justify="left", style="bold bright_cyan")
        table.add_column("Statistic", justify="left", style="bold bright_cyan")
        table.add_column("Value", justify="right", style="bold bright_green")
        for n_sigma in (1, 2, 3):
            table.add_column(
                f"{n_sigma}-sigma [-,+]", justify="right", style="bold bright_green"
            )

        for evaluate, display, sd1 in zip(evaluators, displays, epdfs):
            bestfit_center = evaluate(bf_values)
            for stat in self.stat:
                center, bounds = stat_center_and_bounds(stat, sd1, bestfit_center)
                row_cells = [display, stat.value, f"{center: .6g}"]
                row_cells.extend(f"-{lo:.4g} / +{hi:.4g}" for lo, hi in bounds)
                table.add_row(*row_cells)

        self.console.print(table)

        self.close_logging()


@dataclasses.dataclass(kw_only=True)
class GetBestFit(LoadCatalog):
    """Get best-fit parameters."""

    def __post_init__(self) -> None:
        """Get best-fit parameters."""
        super().__post_init__()

        best_fit = np.array(self.mcat.get_bestfit_row().dup_array(), dtype=np.float64)[
            self.nadd_vals :
        ]
        self.mset.fparams_set_array(best_fit)

        if self.output is None:
            raise typer.BadParameter("Output file not defined.")

        ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
        output_dict = Ncm.ObjDictStr.new()
        output_dict.add("model-set", self.mset)
        ser.dict_str_to_yaml_file(output_dict, self.output.absolute().as_posix())

        self.close_logging()


@dataclasses.dataclass(kw_only=True)
class CheckM2lnL(LoadExperiment):
    """Recompute -2ln(L) for every row of a catalog and compare against the
    stored value.

    Loads the experiment file's likelihood (the CURRENT code) and, for each
    catalog row, sets its free parameters into the experiment's model-set and
    recomputes -2ln(L), reporting how it compares to the value the catalog
    already has stored. Useful for validating catalogs produced by an older
    version of the code against the current one.

    Unlike the other `catalog` commands, this one genuinely needs an
    experiment file: it is the only one that re-evaluates the likelihood
    rather than just reading what the catalog already stored.
    """

    mcmc_file: Annotated[
        Path,
        typer.Argument(help="Path to the MCMC catalog file."),
    ]

    burnin: Annotated[
        int,
        typer.Option(
            help="Number of iterations (ensemble steps) to discard as burnin.",
            min=0,
        ),
    ] = 0

    tail: Annotated[
        Optional[int],
        typer.Option(
            help=(
                "Keep only the last N iterations (ensemble steps) instead of "
                "burning in from the start. Give at most one of --burnin/--tail."
            ),
            min=0,
        ),
    ] = None

    tolerance: Annotated[
        float,
        typer.Option(
            help="Absolute |delta(-2lnL)| above which a row is flagged as "
            "mismatched.",
        ),
    ] = 1.0e-6

    max_report: Annotated[
        int,
        typer.Option(
            help="Maximum number of mismatching rows to print individually.",
        ),
    ] = 20

    stride: Annotated[
        int,
        typer.Option(
            min=1,
            help="Only check every Nth row (1 = check every row). Each row "
            "requires a full likelihood recomputation, which can be slow for "
            "large catalogs -- use this for a quick sanity check instead of "
            "an exhaustive pass.",
        ),
    ] = 1

    max_rows: Annotated[
        Optional[int],
        typer.Option(
            help="Stop after checking this many rows (after --stride). If not "
            "given, all selected rows are checked.",
        ),
    ] = None

    def __post_init__(self) -> None:
        """Recompute -2ln(L) for every catalog row and report mismatches."""
        super().__post_init__()

        loaded: LoadedCatalog = load_catalog(self.mcmc_file, self.burnin, self.tail)

        if self.mset.fparams_len() != loaded.fparams_len:
            raise typer.BadParameter(
                f"Experiment {self.experiment} has {self.mset.fparams_len()} free "
                f"parameters but catalog {self.mcmc_file} has {loaded.fparams_len} "
                f"-- they are not compatible."
            )

        found, m2lnl_col = loaded.mcat.col_by_name(Ncm.MSET_CATALOG_M2LNL_COLNAME)
        if not found:
            raise typer.BadParameter(
                f"Catalog {self.mcmc_file} has no "
                f"{Ncm.MSET_CATALOG_M2LNL_COLNAME} column."
            )

        n_total = loaded.mcat.len()
        selected = list(range(0, n_total, self.stride))
        if self.max_rows is not None:
            selected = selected[: self.max_rows]
        n_rows = len(selected)

        self.console.print(
            f"# Checking {n_rows} of {n_total} row(s) of {self.mcmc_file} "
            f"against {self.experiment} (stride={self.stride})."
        )

        stored = np.empty(n_rows)
        recomputed = np.empty(n_rows)

        for k, i in track(
            list(enumerate(selected)),
            description="Recomputing -2ln(L)...",
            console=self.console,
        ):
            row = np.array(loaded.mcat.peek_row(i).dup_array(), dtype=np.float64)
            fparams = row[loaded.nadd_vals :]
            self.mset.fparams_set_array(fparams)
            stored[k] = row[m2lnl_col]
            recomputed[k] = self.likelihood.m2lnL_val(self.mset)

        diffs = recomputed - stored
        abs_diffs = np.abs(diffs)
        n_mismatch = int(np.sum(abs_diffs > self.tolerance))

        self.console.print(f"# Checked {n_rows} row(s).")
        self.console.print(
            f"# max|delta(-2lnL)| = {abs_diffs.max():.6g}, "
            f"mean|delta(-2lnL)| = {abs_diffs.mean():.6g}, "
            f"mismatched rows (tolerance {self.tolerance:.3g}): "
            f"{n_mismatch}/{n_rows}."
        )

        if n_mismatch > 0:
            worst = np.argsort(-abs_diffs)
            table = Table(title="Worst mismatches")
            table.add_column("row", justify="right")
            table.add_column("stored -2lnL", justify="right")
            table.add_column("recomputed -2lnL", justify="right")
            table.add_column("delta", justify="right")

            n_shown = 0
            for k in worst:
                if abs_diffs[k] <= self.tolerance:
                    break
                if n_shown >= self.max_report:
                    break
                table.add_row(
                    str(selected[k]),
                    f"{stored[k]:.6f}",
                    f"{recomputed[k]:.6f}",
                    f"{diffs[k]:.6g}",
                )
                n_shown += 1

            self.console.print(table)

        self.close_logging()


@dataclasses.dataclass(kw_only=True)
class DumpMset(AppLogging):
    """Dump the model-set stored in an MCMC catalog file as YAML.

    Catalog files store their model-set internally (in the primary FITS
    HDU for new files, or in a legacy `.mset` sidecar file for old ones).
    This command extracts it as a standalone YAML file, useful for
    inspection or as a starting point for a new fit.
    """

    mcmc_file: Annotated[
        Path,
        typer.Argument(help="Path to the MCMC catalog file."),
    ]

    output: Annotated[
        Optional[Path],
        typer.Option(
            "--output",
            "-o",
            help="Path to the output YAML file. If not given, prints to the console.",
        ),
    ] = None

    def __post_init__(self) -> None:
        """Load the catalog and dump its model-set as YAML."""
        super().__post_init__()

        if not self.mcmc_file.exists():
            raise RuntimeError(f"MCMC file {self.mcmc_file} not found.")

        mcat: Ncm.MSetCatalog = Ncm.MSetCatalog.new_from_file_ro(
            self.mcmc_file.absolute().as_posix(), 0
        )
        mset: Ncm.MSet = mcat.peek_mset()
        assert isinstance(mset, Ncm.MSet)

        ser = Ncm.Serialize.new(Ncm.SerializeOpt.NONE)
        if self.output is not None:
            ser.to_yaml_file(mset, self.output.absolute().as_posix())
            self.console.print(f"Model-set written to {self.output}.")
        else:
            self.console.print(ser.to_yaml(mset))

        self.close_logging()

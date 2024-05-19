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
import typer
from rich.table import Table
from rich.text import Text
import numpy as np
import matplotlib.pyplot as plt
import getdist
import getdist.plots

from .. import Ncm
from ..interpolation.stats_dist import (
    create_stats_dist,
    CrossValidationMethod,
    InterpolationKernel,
    InterpolationMethod,
)
from ..plotting.tools import set_rc_params_article, confidence_ellipse
from ..plotting.getdist import mcat_to_mcsamples
from .loading import LoadCatalog


@dataclasses.dataclass
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

            self.end_experiment()
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

        self.end_experiment()


@dataclasses.dataclass
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


@dataclasses.dataclass
class PlotCorner(LoadCatalog):
    """Plots the corner plot of the catalog."""

    plot_name: Annotated[
        Optional[str],
        typer.Option(
            help="Name of the plot file.",
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
            help="Mark the best-fit parameters.",
        ),
    ] = False

    title_limit: Annotated[
        int,
        typer.Option(
            help="Include the n-sigma limit in the title.",
            min=0,
        ),
    ] = 0

    def __post_init__(self) -> None:
        """Corner plot of the catalog."""
        super().__post_init__()

        mcat = self.mcat
        if self.plot_name is None:
            self.plot_name = str(self.mcmc_file)
        mcsample, _, _ = mcat_to_mcsamples(mcat, self.plot_name, indices=self.indices)

        set_rc_params_article(ncol=2)
        g = getdist.plots.get_subplot_plotter(
            width_inch=plt.rcParams["figure.figsize"][0]
        )
        g.settings.linewidth = 0.01
        bf = None
        if self.mark_bestfit:
            bf = np.array(mcat.get_bestfit_row().dup_array())[1:]
        g.triangle_plot(
            [mcsample], shaded=True, markers=bf, title_limit=self.title_limit
        )

        if self.output is not None:
            filename = self.output.with_suffix(".corner.pdf").absolute().as_posix()
            plt.savefig(filename, bbox_inches="tight")

        plt.show()

        self.end_experiment()

#
# wl_catalog.py
#
# Fri Nov 13 2026
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

"""NumCosmo APP commands that write weak lensing catalogs."""

from __future__ import annotations

import dataclasses
import time
from pathlib import Path
from typing import Annotated

import typer
from rich.console import Console

from numcosmo_py import Nc, Ncm
from numcosmo_py.experiments.cluster_wl import (
    GalaxyPopGen,
    GalaxyShapeFactorGenSpectral,
    IntegMethod,
    IntegMethodOptions,
    ShapeFactorGen,
    WLCatalogID,
    check_shape_pop_compat,
    load_cluster_wl,
)

from .generate import (
    ClusterWL,
    load_wl_obs,
    parse_pop_dist,
    parse_shape_factor,
)


def shape_interp_tables_path(output: Path) -> Path:
    """The interpolants file of the prepared catalog OUTPUT: its .gvar suffix
    replaced by .shape_interp.gvar."""
    return output.with_name(f"{output.stem}.shape_interp.gvar")


@dataclasses.dataclass(kw_only=True)
class TablesWLCatalog:
    """Write a catalog with tables of its galaxies.

    Reads a NcGalaxyWLObs (--catalog or --data-file), builds the tables chosen, at
    least one, in parallel with OMP_NUM_THREADS threads, and writes the catalog with
    them to OUTPUT. The rows are unchanged. The two tables are independent: a catalog
    may hold either or both, and one call builds both from one reading of the catalog.

    --cheb-grid builds the CHEB_WEIGHTS redshift grid of every galaxy for the cluster
    redshift of the catalog metadata or --cluster-z (NcGalaxyWLObs:cheb-grid), stored
    in OUTPUT. It depends only on the p(z) of the galaxies, their integration limits
    and the cluster redshift. A dataset under ``--integ-method cheb-weights`` takes it
    from the catalog instead of building it, and every copy made by a sampler shares
    it; a dataset with another cluster redshift aborts.

    --shape-interp builds the interpolant of ln P(eps_obs | g) of every galaxy for the
    population --pop-dist with the NcGalaxyShapeFactorSpectral of --shape-factor,
    writes the interpolants to the file OUTPUT with the suffix ``.shape_interp.gvar``
    in place of ``.gvar``, and stores in OUTPUT the record of what built them and the
    name of that file. ``numcosmo generate cluster-wl-load --data-file OUTPUT`` builds
    an experiment that uses them, with the shape factor and the population of the
    record, and a fit never builds them. The catalog and the experiment store the bare
    name of the interpolants file, so the experiment must be written in its directory,
    and the two files move together. Without --free-pop the population parameters are
    fixed and the interpolants hold for their values only. With --free-pop every
    population parameter is one more variable of the interpolants, over its bounds: an
    experiment built from OUTPUT may fit it (--parameter-list) or hold it fixed at any
    value in its range.

    Without --shape-interp, a catalog that refers to an interpolants file keeps it, and
    OUTPUT must be written in the directory of that file.
    """

    output: Annotated[
        Path, typer.Argument(help="Path of the catalog to write (.gvar).")
    ]

    cheb_grid: Annotated[
        bool,
        typer.Option(
            "--cheb-grid",
            help="Build the CHEB_WEIGHTS redshift grid of the galaxies.",
            rich_help_panel="Tables",
        ),
    ] = False

    shape_interp: Annotated[
        bool,
        typer.Option(
            "--shape-interp",
            help=("Build the shape interpolants of the galaxies; needs --pop-dist."),
            rich_help_panel="Tables",
        ),
    ] = False

    pop_dist: Annotated[
        str | None,
        typer.Option(
            help=GalaxyPopGen.get_help_text() + " Only with --shape-interp.",
            metavar=GalaxyPopGen.get_help_metavar(),
            rich_help_panel="Galaxy intrinsic-ellipticity population distribution",
        ),
    ] = None

    shape_factor: Annotated[
        str | None,
        typer.Option(
            help=(
                ShapeFactorGen.get_help_text()
                + " Must be spectral, only with --shape-interp. Default: spectral "
                "with the catalog's ellipticity convention and frame."
            ),
            metavar=ShapeFactorGen.get_help_metavar(),
            rich_help_panel="Galaxy shape-marginalization factor",
        ),
    ] = None

    catalog: Annotated[
        WLCatalogID | None,
        typer.Option(
            help=(
                "Curated catalog of the NumCosmo data file release. Mutually "
                "exclusive with --data-file; exactly one of the two is required."
            ),
        ),
    ] = None

    data_file: Annotated[
        Path | None,
        typer.Option(
            help=(
                "Path to a NcGalaxyWLObs (.gvar) file. Mutually exclusive with "
                "--catalog; exactly one of the two is required."
            ),
            exists=True,
            readable=True,
        ),
    ] = None

    cluster_ra: Annotated[
        float | None,
        typer.Option(help="Cluster right ascension. Falls back to catalog metadata."),
    ] = None

    cluster_dec: Annotated[
        float | None,
        typer.Option(help="Cluster declination. Falls back to catalog metadata."),
    ] = None

    cluster_z: Annotated[
        float | None,
        typer.Option(help="Cluster redshift. Falls back to catalog metadata."),
    ] = None

    cluster_c: Annotated[
        float | None,
        typer.Option(help="Cluster concentration. Falls back to catalog metadata."),
    ] = None

    free_pop: Annotated[
        bool,
        typer.Option(
            help=(
                "Make every parameter of the population a variable of the "
                "interpolants, over its bounds, so that a fit may vary it. Only "
                "with --shape-interp."
            ),
        ),
    ] = False

    def _check_options(self) -> None:
        """Reject an invalid output, no table, or options of a table not chosen."""
        if self.output.suffix != ".gvar":
            raise typer.BadParameter(f"Invalid output suffix: {self.output.suffix}")
        if not (self.cheb_grid or self.shape_interp):
            raise typer.BadParameter(
                "choose at least one table: --cheb-grid, --shape-interp"
            )
        if self.shape_interp and self.pop_dist is None:
            raise typer.BadParameter("--shape-interp needs --pop-dist")
        if not self.shape_interp and (
            self.pop_dist is not None or self.shape_factor is not None or self.free_pop
        ):
            raise typer.BadParameter(
                "--pop-dist, --shape-factor and --free-pop apply only to --shape-interp"
            )

    def _generators(self, obs: Nc.GalaxyWLObs):
        """The shape factor and population generators of the dataset."""
        if not self.shape_interp:
            # The grid does not depend on the shape factor or the population: the
            # cheapest pair builds the dataset.
            _shape_factor, shape_gen = parse_shape_factor(
                f"var_add "
                f"ellip_conv={obs.get_ellip_conv().value_nick} "
                f"ellip_coord={obs.get_coord().value_nick}"
            )
            _pop_dist, pop_gen = parse_pop_dist("gauss")
            return shape_gen, pop_gen

        shape_factor = self.shape_factor
        if shape_factor is None:
            shape_factor = (
                f"spectral "
                f"ellip_conv={obs.get_ellip_conv().value_nick} "
                f"ellip_coord={obs.get_coord().value_nick}"
            )
        _shape_factor, shape_gen = parse_shape_factor(shape_factor)
        assert self.pop_dist is not None
        _pop_dist, pop_gen = parse_pop_dist(self.pop_dist)

        if not isinstance(shape_gen, GalaxyShapeFactorGenSpectral):
            raise typer.BadParameter("the shape factor must be spectral")
        # The build does not know the priors of the fits that will use the catalog: a
        # table restricted to the shears reachable at the build's centre and mass would
        # leave the shears a fit with a wider prior reaches to the integrals. The tables
        # cover interp-g-max unless restrict_range is given explicitly.
        # (A generator builds its shape factor once, so the option is parsed again.)
        if "restrict_range" not in shape_gen.model_fields_set:
            _shape_factor, shape_gen = parse_shape_factor(
                f"{shape_factor} restrict_range=false"
            )
        check_shape_pop_compat(shape_gen, pop_gen)
        return shape_gen, pop_gen

    def __post_init__(self) -> None:
        """Build the tables chosen and write the catalog.

        Raises:
            typer.BadParameter: Invalid output suffix, no table chosen, options of a
                table not chosen, a shape factor other than spectral, or an output
                outside the directory of the interpolants file the catalog keeps.
        """
        Ncm.cfg_init()
        self._check_options()

        obs = load_wl_obs(self.catalog, self.data_file)

        if self.shape_interp:
            # A catalog prepared before carries interpolants; the new ones replace them.
            obs.set_property("shape-interp", None)
        else:
            shape_interp = obs.get_property("shape-interp")
            if shape_interp is not None and shape_interp.peek_filename() is not None:
                tables = Path(shape_interp.peek_filename())
                if tables.parent != self.output.absolute().parent:
                    raise typer.BadParameter(
                        f"the catalog refers to the interpolants file {tables}; "
                        f"write the output in {tables.parent}"
                    )

        shape_gen, pop_gen = self._generators(obs)
        exp = load_cluster_wl(
            obs=obs,
            cluster_ra=self.cluster_ra,
            cluster_dec=self.cluster_dec,
            cluster_z=self.cluster_z,
            cluster_mass=ClusterWL.cluster_mass,
            cluster_mass_min=ClusterWL.cluster_mass_min,
            cluster_mass_max=ClusterWL.cluster_mass_max,
            cluster_c=self.cluster_c,
            r_min=ClusterWL.r_min,
            r_max=ClusterWL.r_max,
            shape_gen=shape_gen,
            pop_gen=pop_gen,
            integ_options=IntegMethodOptions(integ_method=IntegMethod.CHEB_WEIGHTS),
            summary=False,
        )
        mset = exp.peek("model-set")
        likelihood = exp.peek("likelihood")
        assert isinstance(mset, Ncm.MSet)
        assert isinstance(likelihood, Ncm.Likelihood)
        data = likelihood.peek_dataset().get_data(0)
        assert isinstance(data, Nc.DataClusterWLFactor)

        built: list[str] = []
        if self.shape_interp:
            built.append(self._write_shape_interp(mset, data, obs, pop_gen))
        if self.cheb_grid:
            t_start = time.perf_counter()
            data.write_cheb_grid(mset, obs)
            grid = obs.peek_cheb_grid()
            assert grid is not None
            built.append(
                f"Chebyshev grid of {grid.len()} galaxies for z_cl = "
                f"{grid.get_z_cl()}, built in {time.perf_counter() - t_start:.1f} s"
            )

        ser = Ncm.Serialize.new(Ncm.SerializeOpt.KEEP_SHARED)
        ser.to_binfile(obs, self.output.absolute().as_posix())

        Console().print(f"Wrote {self.output}: " + "; ".join(built) + ".")

    def _write_shape_interp(
        self,
        mset: Ncm.MSet,
        data: Nc.DataClusterWLFactor,
        obs: Nc.GalaxyWLObs,
        pop_gen,
    ) -> str:
        """Build the interpolants into their file; return the summary of the build."""
        pop = mset.peek(Nc.GalaxyShapePop.id())
        assert isinstance(pop, Nc.GalaxyShapePop)
        free = ""
        if self.free_pop:
            for k in range(pop.len()):
                pop.param_set_ftype(k, Ncm.ParamType.FREE)
            free = (
                " with "
                + ", ".join(pop.param_name(k) for k in range(pop.len()))
                + " free"
            )

        # Each interpolant goes to a temporary records file as it is built, and the
        # interpolants file is assembled from them: they are never all in memory.
        # The catalog then refers to that file.
        tables = shape_interp_tables_path(self.output)
        t_start = time.perf_counter()
        data.write_shape_interp_file(mset, obs, tables.absolute().as_posix())
        t_build = time.perf_counter() - t_start

        return (
            f"{tables}, interpolants of {obs.len()} galaxies for {pop_gen!r}{free}, "
            f"built in {t_build:.1f} s"
        )

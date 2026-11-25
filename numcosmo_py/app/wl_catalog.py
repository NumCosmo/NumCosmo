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
    DEFAULT_INTEG_OPTIONS,
    GalaxyPopGen,
    GalaxyShapeFactorGenSpectral,
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


@dataclasses.dataclass(kw_only=True)
class ShapeInterpWLCatalog:
    """Write a catalog with the interpolants of ln P(eps_obs | g) of its galaxies.

    Reads a NcGalaxyWLObs (--catalog or --data-file), builds the interpolant of
    every galaxy for the population --pop-dist with the
    NcGalaxyShapeFactorSpectral of --shape-factor, in parallel with
    OMP_NUM_THREADS threads, and writes the catalog with the interpolants and
    the record of what built them to OUTPUT. The rows are unchanged. The
    interpolants are built only here: ``numcosmo generate cluster-wl-load
    --data-file OUTPUT`` builds an experiment that uses them, with the shape
    factor and the population of the record, and a fit never builds them.
    The experiment stores the absolute path of OUTPUT, not the interpolants,
    and reads them from OUTPUT when it is loaded: OUTPUT must stay at that path.
    Without --free-pop the population parameters are fixed and the interpolants
    hold for their values only. With --free-pop every population parameter is
    one more variable of the interpolants, over its bounds: an experiment built
    from OUTPUT may fit it (--parameter-list) or hold it fixed at any value in
    its range.
    """

    output: Annotated[
        Path, typer.Argument(help="Path of the catalog to write (.gvar).")
    ]

    pop_dist: Annotated[
        str,
        typer.Option(
            help=GalaxyPopGen.get_help_text(),
            metavar=GalaxyPopGen.get_help_metavar(),
            rich_help_panel="Galaxy intrinsic-ellipticity population distribution",
        ),
    ]

    shape_factor: Annotated[
        str | None,
        typer.Option(
            help=(
                ShapeFactorGen.get_help_text()
                + " Must be spectral. Default: spectral with the catalog's "
                "ellipticity convention and frame."
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
                "interpolants, over its bounds, so that a fit may vary it."
            ),
        ),
    ] = False

    def __post_init__(self) -> None:
        """Build the interpolants and write the catalog.

        Raises:
            typer.BadParameter: Invalid output suffix, or a shape factor other
                than spectral.
        """
        Ncm.cfg_init()

        if self.output.suffix != ".gvar":
            raise typer.BadParameter(f"Invalid output suffix: {self.output.suffix}")

        obs = load_wl_obs(self.catalog, self.data_file)
        # A catalog prepared before carries a record; the new one replaces it.
        obs.props.shape_interp_source = Ncm.ObjDictStr.new()

        if self.shape_factor is None:
            self.shape_factor = (
                f"spectral "
                f"ellip_conv={obs.get_ellip_conv().value_nick} "
                f"ellip_coord={obs.get_coord().value_nick}"
            )
        _shape_factor, shape_gen = parse_shape_factor(self.shape_factor)
        _pop_dist, pop_gen = parse_pop_dist(self.pop_dist)

        if not isinstance(shape_gen, GalaxyShapeFactorGenSpectral):
            raise typer.BadParameter("the shape factor must be spectral")
        check_shape_pop_compat(shape_gen, pop_gen)

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
            integ_options=DEFAULT_INTEG_OPTIONS,
            summary=False,
        )
        mset = exp.peek("model-set")
        likelihood = exp.peek("likelihood")
        assert isinstance(mset, Ncm.MSet)
        assert isinstance(likelihood, Ncm.Likelihood)
        data = likelihood.peek_dataset().get_data(0)
        assert isinstance(data, Nc.DataClusterWLFactor)

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
        # catalog is assembled from them: they are never all in memory.
        t_start = time.perf_counter()
        data.write_shape_interp_file(mset, obs, self.output.absolute().as_posix())
        t_build = time.perf_counter() - t_start

        Console().print(
            f"Wrote {self.output}: interpolants of {obs.len()} galaxies for "
            f"{pop_gen!r}{free}, built in {t_build:.1f} s."
        )

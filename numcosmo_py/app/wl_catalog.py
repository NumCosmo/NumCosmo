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
    NcGalaxyShapeFactorSpectral of --shape-factor (use_interp=true), in
    parallel with OMP_NUM_THREADS threads, and writes the catalog with the
    interpolants and the record of what built them to OUTPUT. The rows are
    unchanged. ``numcosmo generate cluster-wl-load --data-file OUTPUT`` then
    builds an experiment that loads the interpolants instead of building them,
    with the shape factor and the population of the record. The population
    parameters are fixed: the interpolants hold for their values only.
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
                + " Must be spectral with use_interp=true. Default: spectral "
                "use_interp=true with the catalog's ellipticity convention and "
                "frame."
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

    cluster_mass_max: Annotated[
        float,
        typer.Option(
            help=(
                "Upper bound of the cluster mass of the fits the catalog is for: the "
                "shear range of restrict_range=true is bounded at this mass."
            ),
            show_default=True,
        ),
    ] = ClusterWL.cluster_mass_max

    centre_box_deg: Annotated[
        float,
        typer.Option(
            help=(
                "Half-width in degrees of the box of the halo center (ra and dec) of the "
                "fits the catalog is for. Zero: the center is fixed. Used by "
                "restrict_range=true."
            ),
            min=0.0,
            show_default=True,
        ),
    ] = 0.0

    cluster_c_min: Annotated[
        float | None,
        typer.Option(
            help=(
                "Lower bound of the concentration of the fits the catalog is for, with "
                "--cluster-c-max: the shear range of restrict_range=true is the largest "
                "over the concentrations of the box. Default: the concentration is fixed."
            ),
        ),
    ] = None

    cluster_c_max: Annotated[
        float | None,
        typer.Option(help="Upper bound of the concentration, with --cluster-c-min."),
    ] = None

    pop_box: Annotated[
        list[str] | None,
        typer.Option(
            help=(
                "NAME=LO:HI, repeatable: frees the population parameter NAME over "
                "[LO, HI]. With pop_order in the spectral shape factor, the catalog "
                "holds the interpolants of the nodes of the expansion over that box, "
                "and fits built from it have those parameters free."
            ),
        ),
    ] = None

    def __post_init__(self) -> None:
        """Build the interpolants and write the catalog.

        Raises:
            typer.BadParameter: Invalid output suffix, or a shape factor other
                than spectral with use_interp=true.
        """
        Ncm.cfg_init()

        if self.output.suffix != ".gvar":
            raise typer.BadParameter(f"Invalid output suffix: {self.output.suffix}")

        obs = load_wl_obs(self.catalog, self.data_file)
        # A catalog prepared before carries a record; the new one replaces it.
        obs.props.shape_interp_source = Ncm.ObjDictStr.new()

        if self.shape_factor is None:
            self.shape_factor = (
                f"spectral use_interp=true "
                f"ellip_conv={obs.get_ellip_conv().value_nick} "
                f"ellip_coord={obs.get_coord().value_nick}"
            )
        _shape_factor, shape_gen = parse_shape_factor(self.shape_factor)
        _pop_dist, pop_gen = parse_pop_dist(self.pop_dist)

        if not (
            isinstance(shape_gen, GalaxyShapeFactorGenSpectral) and shape_gen.use_interp
        ):
            raise typer.BadParameter(
                "the shape factor must be spectral with use_interp=true"
            )
        check_shape_pop_compat(shape_gen, pop_gen)

        exp = load_cluster_wl(
            obs=obs,
            cluster_ra=self.cluster_ra,
            cluster_dec=self.cluster_dec,
            cluster_z=self.cluster_z,
            cluster_mass=ClusterWL.cluster_mass,
            cluster_mass_min=ClusterWL.cluster_mass_min,
            cluster_mass_max=self.cluster_mass_max,
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
        self._set_box(mset)

        t_start = time.perf_counter()
        data.write_shape_interp(mset, obs)
        t_build = time.perf_counter() - t_start

        ser = Ncm.Serialize.new(Ncm.SerializeOpt.KEEP_SHARED)
        ser.to_binfile(obs, self.output.absolute().as_posix())

        Console().print(
            f"Wrote {self.output}: interpolants of {obs.len()} galaxies for "
            f"{pop_gen!r}, built in {t_build:.1f} s."
        )

    def _set_box(self, mset: Ncm.MSet) -> None:
        """Free the halo center, the concentration and the population over the box.

        The box of the center and of the concentration feeds the bound of
        restrict_range=true; the box of the population is that of the expansion of
        pop_order, and the record of the catalog keeps it.

        Raises:
            typer.BadParameter: Only one of --cluster-c-min and --cluster-c-max, or
                a --pop-box that is not NAME=LO:HI for a parameter of the
                population.
        """
        hp = mset.peek(Nc.HaloPosition.id())
        hms = mset.peek(Nc.HaloMassSummary.id())
        pop = mset.peek(Nc.GalaxyShapePop.id())
        assert isinstance(pop, Nc.GalaxyShapePop)
        for item in self.pop_box or []:
            try:
                name, bounds = item.split("=", 1)
                lo, hi = (float(v) for v in bounds.split(":", 1))
            except ValueError as exc:
                raise typer.BadParameter(
                    f"--pop-box {item!r}: expected NAME=LO:HI"
                ) from exc
            found, idx = pop.param_index_from_name(name)
            if not found:
                raise typer.BadParameter(
                    f"--pop-box: {type(pop).__name__} has no parameter {name!r}"
                )
            if not lo < hi:
                raise typer.BadParameter(f"--pop-box {item!r}: LO must be below HI")
            pop.param_set_desc(name, {"lower-bound": lo, "upper-bound": hi})
            pop.param_set_ftype(idx, Ncm.ParamType.FREE)
            pop[name] = 0.5 * (lo + hi)
        assert isinstance(hp, Nc.HaloPosition)
        assert isinstance(hms, Nc.HaloMassSummary)
        if self.centre_box_deg > 0.0:
            for name in ("ra", "dec"):
                v = hp[name]
                hp.param_set_desc(
                    name,
                    {
                        "lower-bound": v - self.centre_box_deg,
                        "upper-bound": v + self.centre_box_deg,
                    },
                )
                _found, idx = hp.param_index_from_name(name)
                hp.param_set_ftype(idx, Ncm.ParamType.FREE)
        if (self.cluster_c_min is None) != (self.cluster_c_max is None):
            raise typer.BadParameter(
                "--cluster-c-min and --cluster-c-max must be given together"
            )
        if self.cluster_c_min is not None and self.cluster_c_max is not None:
            hms.param_set_desc(
                "cDelta",
                {"lower-bound": self.cluster_c_min, "upper-bound": self.cluster_c_max},
            )
            _found, idx = hms.param_index_from_name("cDelta")
            hms.param_set_ftype(idx, Ncm.ParamType.FREE)

#
# cluster_wl.py
#
# Tue Oct 22 12:50:10 2024
# Copyright  2024  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# cluster_wl.py
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

"""Factory functions for LSST Weak Lensing cluster forecast.

This module provides factory functions to create the experiment dictionary
for the LSS Weak Lensing cluster forecast.
"""

from enum import Enum
from dataclasses import dataclass
import numpy as np
from rich.console import Console
from rich.table import Table


from numcosmo_py import Ncm, Nc


class HaloProfileType(str, Enum):
    """Halo density profile types."""

    NFW = "nfw"
    EINASTO = "einasto"
    HERNQUIST = "hernquist"


class GalaxyZDist(str, Enum):
    """Galaxy redshift source distribution types."""

    SPEC = "spec"
    GAUSS = "gauss"


class GalaxySDShapeDist(str, Enum):
    """Galaxy shape source distribution types."""

    GAUSS = "gauss"


@dataclass(frozen=True, kw_only=True)
class HaloPositionData:
    """Cluster position parameters."""

    ra: float
    dec: float
    z: float


@dataclass(frozen=True, kw_only=True)
class GalaxyDistributionData:
    """Galaxy distribution parameters.

    The galaxy distribution is defined by the minimum and maximum values of right
    ascension and declination, and the galaxy density. The galaxy density is given in
    deg^-2.
    """

    ra_min: float
    ra_max: float
    dec_min: float
    dec_max: float
    z_min: float
    z_max: float
    density: float = 18 * 60 * 60  # 18 arcmin^-2 in deg^-2


class ClusterModel:
    """Cluster model parameters."""

    def __init__(
        self,
        *,
        mass_def: Nc.HaloMassSummaryMassDef = Nc.HaloMassSummaryMassDef.CRITICAL,
        delta: float = 200.0,
        profile_type: HaloProfileType = HaloProfileType.NFW,
        position: HaloPositionData = HaloPositionData(ra=12.34, dec=-55.123, z=0.2),
        cluster_c: float = 4.0,
        cluster_mass: float = 1.0e14,
        dist: None | Nc.Distance = None,
    ) -> None:
        """Initialize the cluster model."""
        if dist is None:
            dist = Nc.Distance.new(5.0)

        self.mass_def = mass_def
        self.mass_delta = delta
        self.profile_type = profile_type

        halo_mass_summary = Nc.HaloCMParam.new(mass_def, delta)
        match profile_type:
            case HaloProfileType.NFW:
                self.density_profile = Nc.HaloDensityProfileNFW.new(halo_mass_summary)
            case HaloProfileType.EINASTO:
                self.density_profile = Nc.HaloDensityProfileEinasto.new(
                    halo_mass_summary
                )
            case HaloProfileType.HERNQUIST:
                self.density_profile = Nc.HaloDensityProfileHernquist.new(
                    halo_mass_summary
                )
            case _:
                raise ValueError(f"Invalid halo profile type: {profile_type}")

        self.halo_mass_summary = halo_mass_summary
        self.surface_mass_density = Nc.WLSurfaceMassDensity.new(dist)
        self.halo_position = Nc.HaloPosition.new(dist)

        self.halo_position["ra"] = position.ra
        self.halo_position["dec"] = position.dec
        self.halo_position["z"] = position.z

        self.halo_mass_summary["cDelta"] = cluster_c
        self.halo_mass_summary["log10MDelta"] = np.log10(cluster_mass)

    @property
    def position_data(self) -> HaloPositionData:
        """Return the cluster position."""
        return HaloPositionData(
            ra=self.halo_position["ra"],
            dec=self.halo_position["dec"],
            z=self.halo_position["z"],
        )

    @property
    def mass(self) -> float:
        """Return the cluster mass."""
        return 10.0 ** self.halo_mass_summary["log10MDelta"]

    @property
    def concentration(self) -> float:
        """Return the cluster concentration."""
        return self.halo_mass_summary["cDelta"]

    def prepare(self, cosmo: Nc.HICosmo) -> None:
        """Prepare the cluster model."""
        self.surface_mass_density.prepare(cosmo)
        self.halo_position.prepare(cosmo)


class GalaxyDistributionModel:
    """Galaxy distribution model parameters."""

    def __init__(
        self,
        galaxies: GalaxyDistributionData,
        z_dist: GalaxyZDist = GalaxyZDist.GAUSS,
        sigma_z: float = 0.03,
        shape_type: GalaxySDShapeDist = GalaxySDShapeDist.GAUSS,
        galaxy_shape_e_rms: float = 2.0e-1,
        galaxy_shape_e_sigma: float = 1.0e-4,
    ) -> None:
        """Initialize the galaxy distribution model."""
        self.galaxy_position = Nc.GalaxySDPositionFlat.new(
            galaxies.ra_min, galaxies.ra_max, galaxies.dec_min, galaxies.dec_max
        )
        self.galaxy_redshift_true = Nc.GalaxySDTrueRedshiftLSSTSRD.new()
        # self.galaxy_redshift_true.set_property(
        #    "lim", Ncm.DTuple2.new(galaxies.z_min, galaxies.z_max)
        # )
        dist = self.galaxy_redshift_true.dist(1.0e-5, 1.0e-15)
        frac = dist.eval_pdf(galaxies.z_max) - dist.eval_pdf(galaxies.z_min)
        self.sky_area = (galaxies.ra_max - galaxies.ra_min) * (
            galaxies.dec_max - galaxies.dec_min
        )
        self.n_galaxies = int(galaxies.density * frac * self.sky_area)
        print(
            f"Number of galaxies: {self.n_galaxies} = "
            f"{galaxies.density} * {frac} * {self.sky_area}"
        )

        match z_dist:
            case GalaxyZDist.SPEC:
                self.galaxy_redshift = Nc.GalaxySDObsRedshiftSpec.new(
                    self.galaxy_redshift_true
                )

                def gen_z_spec(mset, z_data, rng):
                    self.galaxy_redshift.gen(mset, z_data, rng)

                self.gen_z = gen_z_spec

            case GalaxyZDist.GAUSS:
                self.galaxy_redshift = Nc.GalaxySDObsRedshiftGauss.new(
                    self.galaxy_redshift_true
                )

                def gen_z_gauss(mset, z_data, rng):
                    self.galaxy_redshift.gen(mset, z_data, sigma_z, rng)

                self.gen_z = gen_z_gauss
            case _:
                raise ValueError(f"Invalid galaxy redshift distribution: {z_dist}")

        match shape_type:
            case GalaxySDShapeDist.GAUSS:
                self.galaxy_shape = Nc.GalaxySDShapeGauss.new()
                self.galaxy_shape["e-rms"] = galaxy_shape_e_rms
                self.e_sigma = galaxy_shape_e_sigma
            case _:
                raise ValueError(f"Invalid galaxy shape distribution: {shape_type}")

    def generate_data(
        self, cosmo: Nc.HICosmo, cluster: ClusterModel, rng: Ncm.RNG
    ) -> tuple[Nc.DataClusterWL, Ncm.MSet]:
        """Generate the galaxy data."""
        cluster.prepare(cosmo)
        cluster_data = Nc.DataClusterWL.new()
        mset = Ncm.MSet.new_array(
            [
                cosmo,
                cluster.density_profile,
                cluster.surface_mass_density,
                cluster.halo_position,
                self.galaxy_redshift,
                self.galaxy_position,
                self.galaxy_shape,
            ]
        )
        z_data = Nc.GalaxySDObsRedshiftData.new(self.galaxy_redshift)
        p_data = Nc.GalaxySDPositionData.new(self.galaxy_position, z_data)
        s_data = Nc.GalaxySDShapeData.new(self.galaxy_shape, p_data)

        obs = Nc.GalaxyWLObs.new(
            Nc.GalaxyWLObsCoord.EUCLIDEAN,
            self.n_galaxies,
            list(s_data.required_columns()),
        )

        for i in range(self.n_galaxies):
            self.gen_z(mset, z_data, rng)
            self.galaxy_position.gen(mset, p_data, rng)
            self.galaxy_shape.gen(
                mset,
                s_data,
                self.e_sigma,
                self.e_sigma,
                Nc.GalaxyWLObsCoord.EUCLIDEAN,
                rng,
            )
            s_data.write_row(obs, i)

        cluster_data.set_obs(obs)

        return cluster_data, mset


def create_cosmo() -> Nc.HICosmo:
    """Create a cosmology for the cluster model."""
    cosmo = Nc.HICosmoDEXcdm()

    cosmo.params_set_default_ftype()
    cosmo.omega_x2omega_k()
    cosmo["H0"] = 67.81
    cosmo["Omegab"] = 0.0486
    cosmo["Omegac"] = 0.2612
    cosmo["w"] = -1.0
    cosmo["Omegak"] = 0.00

    prim = Nc.HIPrimPowerLaw.new()
    prim["ln10e10ASA"] = 3.02745
    prim["n_SA"] = 0.9660

    reion = Nc.HIReionCamb.new()

    cosmo.param_set_desc("H0", {"fit": False})
    cosmo.param_set_desc("Omegac", {"fit": False})
    cosmo.param_set_desc("Omegab", {"fit": False})
    cosmo.param_set_desc("w", {"fit": False})
    cosmo.param_set_desc("Omegak", {"fit": False})
    prim.param_set_desc("ln10e10ASA", {"fit": False})
    prim.param_set_desc("n_SA", {"fit": False})
    reion.param_set_desc("z_re", {"fit": False})

    cosmo.add_submodel(prim)
    cosmo.add_submodel(reion)

    return cosmo


def generate_lsst_cluster_wl(
    *,
    cluster_ra: float,
    cluster_dec: float,
    cluster_z: float,
    cluster_mass: float,
    cluster_mass_min: float,
    cluster_mass_max: float,
    cluster_c: float,
    ra_min: float,
    ra_max: float,
    dec_min: float,
    dec_max: float,
    z_min: float,
    z_max: float,
    z_dist: GalaxyZDist,
    sigma_z: float,
    shape_dist: GalaxySDShapeDist,
    galaxy_shape_e_rms: float,
    galaxy_shape_e_sigma: float,
    seed: None | int,
    summary: bool,
) -> Ncm.ObjDictStr:
    """Generate J-Pas forecast 2024 experiment dictionary."""
    halo_position = HaloPositionData(ra=cluster_ra, dec=cluster_dec, z=cluster_z)
    cluster = ClusterModel(
        position=halo_position,
        cluster_mass=cluster_mass,
        cluster_c=cluster_c,
    )
    galaxy_distribution = GalaxyDistributionData(
        ra_min=ra_min,
        ra_max=ra_max,
        dec_min=dec_min,
        dec_max=dec_max,
        z_min=z_min,
        z_max=z_max,
    )
    galaxy_model = GalaxyDistributionModel(
        galaxies=galaxy_distribution,
        z_dist=z_dist,
        sigma_z=sigma_z,
        shape_type=shape_dist,
        galaxy_shape_e_rms=galaxy_shape_e_rms,
        galaxy_shape_e_sigma=galaxy_shape_e_sigma,
    )
    if seed is not None:
        rng = Ncm.RNG.seeded_new("mt19937", seed)
    else:
        rng = Ncm.RNG.new("mt19937")
        rng.set_random_seed(True)

    cluster.halo_position.param_set_desc(
        "ra", {"lower-bound": ra_min, "upper-bound": ra_max}
    )
    cluster.halo_position.param_set_desc(
        "dec", {"lower-bound": dec_min, "upper-bound": dec_max}
    )
    cluster.halo_mass_summary.param_set_desc(
        "log10MDelta",
        {
            "lower-bound": float(np.log10(cluster_mass_min)),
            "upper-bound": float(np.log10(cluster_mass_max)),
        },
    )

    cluster_data, mset = galaxy_model.generate_data(create_cosmo(), cluster, rng)

    dset = Ncm.Dataset.new_array([cluster_data])
    likelihood = Ncm.Likelihood.new(dset)

    # Save experiment
    experiment = Ncm.ObjDictStr()

    experiment.set("likelihood", likelihood)
    experiment.set("model-set", mset)

    if summary:
        # Here we use typer and rich to print the summary
        # We print the Halo parameters
        # Then the Galaxy parameters
        # The number of galaxies
        console = Console()

        table = Table(title="Cluster Parameters")
        table.add_column("Parameter")
        table.add_column("Value")
        table.add_row("RA", f"{halo_position.ra}")
        table.add_row("DEC", f"{halo_position.dec}")
        table.add_row("z", f"{halo_position.z}")
        table.add_row("Mass Profile", f"{cluster.profile_type}")
        table.add_row("Mass Definition", f"{cluster.mass_def.value_name}")
        table.add_row("Mass Delta", f"{cluster.mass_delta}")
        table.add_row("Mass", f"{cluster.mass:.2e}")
        table.add_row("Concentration", f"{cluster.concentration}")
        console.print(table)

        table = Table(title="Galaxy Sample Parameters")
        table.add_column("Parameter")
        table.add_column("Value")
        table.add_row("RA min", f"{galaxy_distribution.ra_min}")
        table.add_row("RA max", f"{galaxy_distribution.ra_max}")
        table.add_row("DEC min", f"{galaxy_distribution.dec_min}")
        table.add_row("DEC max", f"{galaxy_distribution.dec_max}")
        table.add_row("z min", f"{galaxy_distribution.z_min}")
        table.add_row("z max", f"{galaxy_distribution.z_max}")
        table.add_row("Density", f"{galaxy_distribution.density}")
        table.add_row("z Distribution", f"{z_dist}")
        if z_dist != GalaxyZDist.SPEC:
            table.add_row("z Sigma", f"{sigma_z}")
        table.add_row("Shape Distribution", f"{shape_dist}")
        table.add_row("Shape e_rms", f"{galaxy_shape_e_rms}")
        table.add_row("Shape e_sigma", f"{galaxy_shape_e_sigma}")

        console.print(table)

        console.print(f"Number of galaxies: {galaxy_model.n_galaxies}")

    return experiment

"""Module to match objects in the sky halo-halo, cluster-halo, cluster-cluster."""

from typing import Any, TypedDict
from pathlib import Path
import tqdm
import numpy as np

from astropy.io import fits
from astropy.table import Table
from numcosmo_py import Ncm, Nc
from numcosmo_py.helper import npa_to_seq

Ncm.cfg_init()


class Coordinates(TypedDict, total=False):
    """Coordinates mapping.

    Dictionary to map the coordinates to the name of the columns in the FITS file.
    """

    RA: str
    DEC: str
    z: str


def _check_coordinates(coordinates: Coordinates) -> None:
    """Check if the coordinates are provided."""
    if "RA" not in coordinates or "DEC" not in coordinates:
        raise ValueError("RA and DEC coordinates must be provided.")


def _load_fits_data(catalog: Path) -> fits.FITS_rec:
    """Load FITS data from the provided catalog path.

    The function returns the first FITS table found in the provided catalog.
    """
    hdul1 = fits.open(catalog.as_posix())
    hdu1_data: fits.FITS_rec | None = None
    for hdu in hdul1:
        if isinstance(hdu, fits.TableHDU):
            hdu1_data = hdu.data
            break
    if hdu1_data is None:
        raise ValueError("No FITS table found in the provided catalog.")

    return hdu1_data


class NcSkyMatching:
    """Class to match objects in the sky halo-halo, cluster-halo, cluster-cluster."""

    def __init__(
        self,
        query_catalog_path: Path,
        query_coordinates: Coordinates,
        match_catalog_path: Path,
        match_coordinates: Coordinates,
        query_properties: list[str] | None = None,
        match_properties: list[str] | None = None,
    ):
        """Initialize the class."""
        self.query_catalog_path = query_catalog_path
        self.match_catalog_path = match_catalog_path

        _check_coordinates(query_coordinates)
        _check_coordinates(match_coordinates)

        self.query_coordinates = query_coordinates
        self.match_coordinates = match_coordinates

        self.query_properties = query_properties
        self.match_properties = match_properties

    def ra_dec_to_theta_phi(self, ra, dec):
        """Convert RA and DEC to theta and phi."""
        ra_rad = np.radians(ra)
        dec_rad = np.radians(dec)

        theta = np.pi / 2 - dec_rad
        phi = ra_rad

        return theta, phi

    def match_3d(
        self,
        cosmo: Nc.HICosmo,
        matching_distance: float,
        n_nearest_neighbours: int,
        verbose: bool = True,
    ) -> Table:
        """Match objects in the sky.

        The function matches objects in the sky using the provided matching distance and
        number of nearest neighbours. It considers the tridimensional space to match the
        objects computing the cosmological distances between the objects.
        """
        if ("z" not in self.match_coordinates) or ("z" not in self.query_coordinates):
            raise ValueError(
                "To perform a 3D matching, "
                "the redshift must be provided for both catalogs."
            )

        matched: dict[Any, Any] = {
            "ID": [],
            "RA": [],
            "DEC": [],
            "z": [],
            "ID_matched": [],
            "RA_matched": [],
            "DEC_matched": [],
            "z_matched": [],
            "distances Mpc": [],
        }
        if self.query_properties is not None:
            for prop in self.query_properties:
                matched[prop] = []

        if self.match_properties is not None:
            for prop in self.match_properties:
                matched[prop] = []

        query_data = _load_fits_data(self.query_catalog_path)
        match_data = _load_fits_data(self.match_catalog_path)

        # Print columns for each file
        theta_q, phi_q = self.ra_dec_to_theta_phi(
            query_data[self.query_coordinates["RA"]],
            query_data[self.query_coordinates["DEC"]],
        )
        z_q = query_data[self.query_coordinates["z"]]
        theta_m, phi_m = self.ra_dec_to_theta_phi(
            match_data[self.match_coordinates["RA"]],
            match_data[self.match_coordinates["DEC"]],
        )
        z_m = match_data[self.match_coordinates["z"]]

        snn = Ncm.SphereNN()
        dist = Nc.Distance.new(3.0)
        dist.prepare(cosmo)

        r_m = np.array(
            [dist.comoving(cosmo, z) for z in z_m],
            dtype=np.float64,
        )
        snn.insert_array(npa_to_seq(r_m), theta_m, phi_m)
        snn.rebuild()

        loop_arg = enumerate(zip(theta_q, phi_q, z_q))
        for i, (theta, phi, z) in (
            tqdm.tqdm(loop_arg, total=len(z_q)) if verbose else loop_arg
        ):
            r = dist.comoving(cosmo, z)
            distances_list, indices = snn.knn_search_distances(
                r, theta, phi, n_nearest_neighbours
            )
            distances = np.sqrt(distances_list) * cosmo.RH_Mpc()

            matched["RA"].append(match_data[self.match_coordinates["RA"]][i])
            matched["DEC"].append(match_data[self.match_coordinates["DEC"]][i])
            matched["z"].append(match_data[self.match_coordinates["z"]][i])
            matched["ID"].append(i)

            if self.match_properties is not None:
                for prop in self.match_properties:
                    matched[prop].append(match_data[prop][i])

            # Select only the halos that are within the matching distance
            matching_dist_indices = distances < matching_distance
            distances_matched = distances[matching_dist_indices]
            ID_matched = indices[matching_dist_indices]

            # Get the matched halos properties
            RA_matched = self.query_coordinates["RA"][ID_matched]
            DEC_matched = self.query_coordinates["DEC"][ID_matched]
            z_matched = self.query_coordinates["z"][ID_matched]
            if self.query_properties is not None:
                for prop in self.query_properties:
                    matched[prop].append(query_data[prop][ID_matched])

            matched["ID_matched"].append(ID_matched)
            matched["distances Mpc"].append(distances_matched)
            matched["RA_matched"].append(RA_matched)
            matched["DEC_matched"].append(DEC_matched)
            matched["z_matched"].append(z_matched)

        return Table(matched)

    def match_2d(
        self,
        matching_distance: float,
        n_nearest_neighbours: int,
        verbose: bool = True,
    ) -> Table:
        """Match objects in the sky.

        The function matches objects in the sky using the provided matching distance.
        """
        matched: dict[Any, Any] = {
            "ID": [],
            "RA": [],
            "DEC": [],
            "ID_matched": [],
            "RA_matched": [],
            "DEC_matched": [],
            "distances": [],
        }
        if self.query_properties is not None:
            for prop in self.query_properties:
                matched[prop] = []

        if self.match_properties is not None:
            for prop in self.match_properties:
                matched[prop] = []

        query_data = _load_fits_data(self.query_catalog_path)
        match_data = _load_fits_data(self.match_catalog_path)
        if (matching_distance < 0) or (matching_distance > np.pi):
            raise ValueError("The matching distance must be between 0 and pi.")

        # Print columns for each file
        theta_q, phi_q = self.ra_dec_to_theta_phi(
            query_data[self.query_coordinates["RA"]],
            query_data[self.query_coordinates["DEC"]],
        )
        theta_m, phi_m = self.ra_dec_to_theta_phi(
            match_data[self.match_coordinates["RA"]],
            match_data[self.match_coordinates["DEC"]],
        )
        r_m = np.ones_like(theta_m)

        snn = Ncm.SphereNN()
        snn.insert_array(r_m, phi_m, theta_m)
        snn.rebuild()

        loop_arg = enumerate(zip(theta_q, phi_q))
        for i, (theta, phi) in (
            tqdm.tqdm(loop_arg, total=len(theta_q)) if verbose else loop_arg
        ):
            distances_list, indices = snn.knn_search_distances(
                1.0, theta, phi, n_nearest_neighbours
            )

            # Below we convert the euclidean distances between two points in the sphere
            # with radius 1 to the angular distances between the two points.
            distances = 2 * np.arcsin(np.sqrt(distances_list) / 2)

            matched["RA"].append(match_data[self.match_coordinates["RA"]][i])
            matched["DEC"].append(match_data[self.match_coordinates["DEC"]][i])
            matched["ID"].append(i)

            if self.match_properties is not None:
                for prop in self.match_properties:
                    matched[prop].append(match_data[prop][i])

            matching_distances_indices = distances < matching_distance
            indices = indices[matching_distances_indices]
            distances = distances[matching_distances_indices]

            matched["ID_matched"].append(indices)
            matched["distances"].append(distances)
            matched["RA_matched"].append(
                match_data[self.match_coordinates["RA"]][indices]
            )
            matched["DEC_matched"].append(
                match_data[self.match_coordinates["DEC"]][indices]
            )

            if self.query_properties is not None:
                for prop in self.query_properties:
                    matched[prop].append(query_data[prop][indices])

        return Table(matched)

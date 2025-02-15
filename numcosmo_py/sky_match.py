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


# ------------------------------------------------------------------------------------#
# TypedDict: Coordinates: Class
# -------------------------------------------------------------------------------------------------#
class Coordinates(TypedDict, total=False):
    """Coordinates mapping.

    Dictionary to map the coordinates to the name of the columns in the FITS file.
    """

    RA: str
    DEC: str
    z: str


# -----Check if the coordinates are provided-------------------------------------------------------#
# _check_coordinates()
# fits.FITS_rec: table:       fits table to be used with the coordinates;
# Coordinates:   coordinates: coordinates we need to check to do the matching;
# Returns:
# Check if the coordinates are provided
# -------------------------------------------------------------------------------------------------#
def _check_coordinates(table: fits.FITS_rec, coordinates: Coordinates) -> None:
    """Check if the coordinates are provided."""
    if ("RA" not in coordinates) or ("DEC" not in coordinates):
        raise ValueError("RA and DEC coordinates must be provided.")

    if (coordinates["RA"] not in table.names) or (
        coordinates["DEC"] not in table.names
    ):
        raise ValueError(
            f"RA and DEC coordinates mapped by {coordinates} "
            f"not found in the provided catalog {table.names}."
        )

    if "z" in coordinates and coordinates["z"] not in table.names:
        raise ValueError(
            f"Redshift coordinate mapped by {coordinates} "
            f"not found in the provided catalog {table.names}."
        )


# -----Load FITS data from the provided catalog path-----------------------------------------------#
# _load_fits_data()
# Path: catalog:      Path of the catalog to load;
# Returns: hdu1_data: data from  the fits catalog
# -------------------------------------------------------------------------------------------------#
def _load_fits_data(catalog: Path) -> fits.FITS_rec:
    """Load FITS data from the provided catalog path.

    The function returns the first FITS table found in the provided catalog.
    """
    hdul1 = fits.open(catalog.as_posix())
    hdu1_data: fits.FITS_rec | None = None
    for hdu in hdul1:
        if isinstance(hdu, (fits.TableHDU, fits.BinTableHDU)):
            hdu1_data = hdu.data
            break
    if hdu1_data is None:
        raise ValueError("No FITS table found in the provided catalog.")

    return hdu1_data


# -------------------------------------------------------------------------------------------------#
# SkyMatch: Class
# Path:        query_catalog_path: Path of the catalog to be matched;
# Coordinates: query_coordinates:  Coordinates of the catalog to be matched;
# Path:        match_catalog_path: Path of the catalog we are searching the match;
# Coordinates: match_coordinates:  Coordinates of the catalog we are searching the match;
# list[str]:   query_properties:   Properties we want to keep from the catalog to be matched;
# list[str]:   match_properties:   Properties we want to keep from the catalog we are searching the match;
# -------------------------------------------------------------------------------------------------#
class SkyMatch:
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

        self.query_data = _load_fits_data(self.query_catalog_path)
        self.match_data = _load_fits_data(self.match_catalog_path)

        _check_coordinates(self.query_data, query_coordinates)
        _check_coordinates(self.match_data, match_coordinates)

        self.query_coordinates = query_coordinates
        self.match_coordinates = match_coordinates

        self.query_properties = query_properties
        self.match_properties = match_properties

    # ----- Convert RA and DEC to theta and phi -------------------------------------------------------#
    # ra_dec_to_theta_phi()
    # float: ra:  right ascencion angle in degrees;
    # float: dec: declination angle in degrees
    # Returns: float: theta: theta in radiands, phi: phi in radiands
    # -------------------------------------------------------------------------------------------------#

    def ra_dec_to_theta_phi(self, ra, dec):
        """Convert RA and DEC to theta and phi."""
        ra_rad = np.radians(ra)
        dec_rad = np.radians(dec)

        theta = np.pi / 2 - dec_rad
        phi = ra_rad

        return theta, phi

    # ----- Performs the 3d match between the two catalogs -----------------------------------------------------------------------------------------------#
    # match_3d()
    # NcHicosmo: cosmo:                cosmology object to convert z in distance
    # float:     matching_distance:    maximum distance to consider a match
    # int:       n_nearest_neighbours: number of nearest neighbours
    # bool:      verbose:              show the progress of the match
    # str:       selection_criteria:   selection criteria to choose the best canditate can be more_massive , distances, redshift_proximity
    # Returns: astropy_table: matched: table with all candidates of matched objects, best_matched: table with the  best candidate of matched objects
    # ----------------------------------------------------------------------------------------------------------------------------------------------------#
    def match_3d(
        self,
        cosmo: Nc.HICosmo,
        matching_distance: float,
        n_nearest_neighbours: int,
        verbose: bool = True,
        selection_criteria: str = "distances",
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
            "distances": [],
        }
        best_matched: dict[Any, Any] = {
            "ID": [],
            "RA": [],
            "DEC": [],
            "z": [],
            "ID_matched": [],
            "RA_matched": [],
            "DEC_matched": [],
            "z_matched": [],
            "distances": [],
        }
        # String necessary to count the number of objects that find at least one match
        matches = 0
        # Creating the columns for the final table
        if self.query_properties is not None:
            for prop in self.query_properties:
                matched[self.query_properties.get(prop)] = []
                best_matched[self.query_properties.get(prop)] = []

        if self.match_properties is not None:
            for prop in self.match_properties:
                matched[self.match_properties.get(prop)] = []
                best_matched[self.match_properties.get(prop)] = []

        theta_q, phi_q = self.ra_dec_to_theta_phi(
            self.query_data[self.query_coordinates["RA"]],
            self.query_data[self.query_coordinates["DEC"]],
        )
        z_q = self.query_data[self.query_coordinates["z"]]
        theta_m, phi_m = self.ra_dec_to_theta_phi(
            self.match_data[self.match_coordinates["RA"]],
            self.match_data[self.match_coordinates["DEC"]],
        )
        z_m = self.match_data[self.match_coordinates["z"]]

        # Preparing the k-nearest-neighbours and cosmology objects for the matching
        snn = Ncm.SphereNN()
        dist = Nc.Distance.new(3.0)
        dist.prepare(cosmo)
        RH_Mpc = cosmo.RH_Mpc()

        r_m = np.array(
            [dist.angular_diameter(cosmo, z) for z in z_m],
            dtype=np.float64,
        )

        snn.insert_array(r_m, theta_m, phi_m)
        snn.rebuild()

        loop_arg = enumerate(zip(theta_q, phi_q, z_q))
        for i, (theta, phi, z) in (
            tqdm.tqdm(loop_arg, total=len(theta_q)) if verbose else loop_arg
        ):
            r = dist.angular_diameter(cosmo, z)

            distances_list, indices_list = snn.knn_search_distances(
                r, theta, phi, n_nearest_neighbours
            )

            # Adding the coordinates and propertties of the object searching for matches
            indices = np.array(indices_list)
            matched["RA"].append(self.query_data[self.query_coordinates["RA"]][i])
            matched["DEC"].append(self.query_data[self.query_coordinates["DEC"]][i])
            matched["z"].append(self.query_data[self.query_coordinates["z"]][i])
            matched["ID"].append(i)

            best_matched["RA"].append(self.query_data[self.query_coordinates["RA"]][i])
            best_matched["DEC"].append(
                self.query_data[self.query_coordinates["DEC"]][i]
            )
            best_matched["z"].append(self.query_data[self.query_coordinates["z"]][i])
            best_matched["ID"].append(i)

            if self.query_properties is not None:
                for prop in self.query_properties:
                    matched[self.query_properties.get(prop)].append(
                        self.query_data[prop][i]
                    )
                    best_matched[self.query_properties.get(prop)].append(
                        self.query_data[prop][i]
                    )

            # Checking the k-nearest-neighbours that satisfies the distance condition

            distances = np.sqrt(distances_list) * RH_Mpc
            matching_distances_indices = distances <= matching_distance

            indices = indices[matching_distances_indices]
            distances = distances[matching_distances_indices]

            # Adding the matched coordinates and properties for the multiple matched table

            matched["ID_matched"].append(indices)
            matched["distances"].append(distances)
            matched["RA_matched"].append(
                self.match_data[self.match_coordinates["RA"]][indices]
            )
            matched["DEC_matched"].append(
                self.match_data[self.match_coordinates["DEC"]][indices]
            )
            matched["z_matched"].append(
                self.match_data[self.match_coordinates["z"]][indices]
            )

            if self.match_properties is not None:
                for prop in self.match_properties:
                    matched[self.match_properties.get(prop)].append(
                        self.match_data[prop][indices]
                    )

            # Selecting the best candidate based on the choosen criteria

            if len(indices) > 0:
                matches += 1
                match selection_criteria:
                    case "distances":
                        min_distances = min(distances)
                        index = np.where(distances == min_distances)[0][0]
                    case "redshift_proximity":
                        redshift_proximity = abs(
                            matched["z"][i]
                            - self.match_data[self.match_coordinates["z"]][indices]
                        )
                        min_redshift_proximity = min(redshift_proximity)
                        index = np.where(redshift_proximity == min_redshift_proximity)[
                            0
                        ][0]
                    case "more_massive":

                        if (
                            "mass1" not in self.match_properties.values()
                            and "mass2" not in self.match_properties.values()
                        ):
                            raise Exception(
                                "match_properties must have a mass1 or mass2 column to use this criteria"
                            )
                        else:
                            masses = list(
                                filter(
                                    lambda key: self.match_properties[key] == "mass1"
                                    or self.match_properties[key] == "mass2",
                                    self.match_properties,
                                )
                            )[0]

                            max_mass = max(self.match_data[masses][indices])
                            index = np.where(
                                self.match_data[masses][indices] == max_mass
                            )[0][0]

                    case _:
                        raise Exception(
                            "selection_criteria must be eihter distances , redshift_proximity or more_massive. Was given %s"
                            % (selection_criteria)
                        )

                # Adding the matched coordinates and properties for the best matched table

                best_matched["ID_matched"].append(indices[index])
                best_matched["distances"].append(distances[index])
                best_matched["RA_matched"].append(
                    self.match_data[self.match_coordinates["RA"]][indices[index]]
                )
                best_matched["DEC_matched"].append(
                    self.match_data[self.match_coordinates["DEC"]][indices[index]]
                )
                best_matched["z_matched"].append(
                    self.match_data[self.match_coordinates["z"]][indices[index]]
                )

                if self.match_properties is not None:
                    for prop in self.match_properties:
                        best_matched[self.match_properties.get(prop)].append(
                            self.match_data[prop][indices[index]]
                        )

            else:
                best_matched["ID_matched"].append(indices)
                best_matched["distances"].append(distances)
                best_matched["RA_matched"].append(
                    self.match_data[self.match_coordinates["RA"]][indices]
                )
                best_matched["DEC_matched"].append(
                    self.match_data[self.match_coordinates["DEC"]][indices]
                )
                best_matched["z_matched"].append(
                    self.match_data[self.match_coordinates["z"]][indices]
                )
                if self.match_properties is not None:
                    for prop in self.match_properties:
                        best_matched[self.match_properties.get(prop)].append(
                            self.match_data[prop][indices]
                        )

        # Summary of the number of matches found

        print(
            """
        Number of objects matching: %s
        Number of objects in the target catalog: %s
        Number of unmatched objects: %s
        Number of matched objects: %s
        """
            % (
                len(best_matched["ID"]),
                len(r_m),
                len(best_matched["ID"]) - matches,
                matches,
            )
        )

        return Table(matched), Table(best_matched)

    # ----- Performs the 2d match between the two catalogs -----------------------------------------------------------------------------------------------#
    # match_2d()
    # NcHicosmo: cosmo:                cosmology object to convert z in distance
    # float:     matching_distance:    maximum distance to consider a match
    # int:       n_nearest_neighbours: number of nearest neighbours
    # bool:      verbose:              show the progress of the match
    # bool:      match_dist_3d:        If converts from angular separation to physical distance
    # float:     delta_z:              Uncertainty associated with the measurament z to match
    # float:     n_delta_z:            Number of sigmas we expect to cut the objects in z space
    # str:       selection_criteria:   selection criteria to choose the best canditate can be more_massive , distances, redshift_proximity
    # bool:      use_zerr_match:       If use the z_err of the match catalog instead of delta_z
    # bool:      use_zerr_query:       If use the z_err of the query catalog instead of delta_z
    # str:       which_radius:         The radius to be used to convert from angular separation to physical distance can be query_radius , match_radius ,
    #                                  min_radius , max_radius
    # Returns: astropy_table: matched: table with all candidates of matched objects, best_matched: table with the  best candidate of matched objects
    # ----------------------------------------------------------------------------------------------------------------------------------------------------#

    def match_2d(
        self,
        cosmo: Nc.HICosmo,
        matching_distance: float,
        n_nearest_neighbours: int,
        verbose: bool = True,
        match_dist_3d: bool = False,
        delta_z: float = 0.1,
        n_delta_z: int = 1,
        selection_criteria: str = "distances",
        use_zerr_match: bool = False,
        use_zerr_query: bool = False,
        which_radius: str = "query_radius",
    ) -> Table:
        """Match objects in the sky.

        The function matches objects in the sky using the provided matching distance.
        """
        if ("z" not in self.match_coordinates) or ("z" not in self.query_coordinates):
            raise ValueError(
                "To perform a  matching, "
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
            "distances": [],
        }
        best_matched: dict[Any, Any] = {
            "ID": [],
            "RA": [],
            "DEC": [],
            "z": [],
            "ID_matched": [],
            "RA_matched": [],
            "DEC_matched": [],
            "z_matched": [],
            "distances": [],
        }
        # String necessary to count the number of objects that find at least one match

        matches = 0

        # Creating the columns for the final table
        if self.query_properties is not None:
            for prop in self.query_properties:
                matched[self.query_properties.get(prop)] = []
                best_matched[self.query_properties.get(prop)] = []

        if self.match_properties is not None:
            for prop in self.match_properties:
                matched[self.match_properties.get(prop)] = []
                best_matched[self.match_properties.get(prop)] = []

        if (matching_distance < 0) or (matching_distance > np.pi):
            raise ValueError("The matching distance must be between 0 and pi.")

        # Preparing the k-nearest-neighbours and cosmology objects for the matching

        theta_q, phi_q = self.ra_dec_to_theta_phi(
            self.query_data[self.query_coordinates["RA"]],
            self.query_data[self.query_coordinates["DEC"]],
        )
        z_q = self.query_data[self.query_coordinates["z"]]
        theta_m, phi_m = self.ra_dec_to_theta_phi(
            self.match_data[self.match_coordinates["RA"]],
            self.match_data[self.match_coordinates["DEC"]],
        )

        r_m_ones = np.ones_like(theta_m)

        snn = Ncm.SphereNN()
        dist = Nc.Distance.new(3.0)
        dist.prepare(cosmo)
        RH_Mpc = cosmo.RH_Mpc()

        snn.insert_array(r_m_ones, theta_m, phi_m)
        snn.rebuild()

        loop_arg = enumerate(zip(theta_q, phi_q, z_q))
        for i, (theta, phi, z) in (
            tqdm.tqdm(loop_arg, total=len(theta_q)) if verbose else loop_arg
        ):

            distances_list, indices_list = snn.knn_search_distances(
                1.0, theta, phi, n_nearest_neighbours
            )

            # Adding the coordinates and propertties of the object searching for matches

            indices = np.array(indices_list)
            matched["RA"].append(self.query_data[self.query_coordinates["RA"]][i])
            matched["DEC"].append(self.query_data[self.query_coordinates["DEC"]][i])
            matched["z"].append(self.query_data[self.query_coordinates["z"]][i])
            matched["ID"].append(i)

            best_matched["RA"].append(self.query_data[self.query_coordinates["RA"]][i])
            best_matched["DEC"].append(
                self.query_data[self.query_coordinates["DEC"]][i]
            )
            best_matched["z"].append(self.query_data[self.query_coordinates["z"]][i])
            best_matched["ID"].append(i)

            if self.query_properties is not None:
                for prop in self.query_properties:
                    matched[self.query_properties.get(prop)].append(
                        self.query_data[prop][i]
                    )
                    best_matched[self.query_properties.get(prop)].append(
                        self.query_data[prop][i]
                    )

            # Checking the k-nearest-neighbours that satisfies the distance condition based on the choosen radius

            if match_dist_3d:
                match which_radius:
                    case "query_radius":
                        r = dist.angular_diameter(cosmo, z) * RH_Mpc
                        distances = r * np.sqrt(distances_list)
                    case "match_radius":
                        distances = np.array(
                            [
                                RH_Mpc
                                * dist.angular_diameter(
                                    cosmo,
                                    self.match_data[self.match_coordinates["z"]][
                                        indices[i]
                                    ],
                                )
                                * np.sqrt(distances_list[i])
                                for i in range(len(indices))
                            ],
                            dtype=np.float64,
                        )
                    case "max_radius":
                        r = dist.angular_diameter(cosmo, z) * RH_Mpc
                        distances = np.array(
                            [
                                max(
                                    RH_Mpc
                                    * dist.angular_diameter(
                                        cosmo,
                                        self.match_data[self.match_coordinates["z"]][
                                            indices[i]
                                        ],
                                    )
                                    * np.sqrt(distances_list[i]),
                                    r * np.sqrt(distances_list[i]),
                                )
                                for i in range(len(indices))
                            ],
                            dtype=np.float64,
                        )
                    case "min_radius":
                        r = dist.angular_diameter(cosmo, z) * RH_Mpc
                        distances = np.array(
                            [
                                min(
                                    RH_Mpc
                                    * dist.angular_diameter(
                                        cosmo,
                                        self.match_data[self.match_coordinates["z"]][
                                            indices[i]
                                        ],
                                    )
                                    * np.sqrt(distances_list[i]),
                                    r * np.sqrt(distances_list[i]),
                                )
                                for i in range(len(indices))
                            ],
                            dtype=np.float64,
                        )
                    case _:
                        raise Exception(
                            "which_radius must be eihter query_radius , match_radius , max_radius or min_radius. Was given %s"
                            % (which_radius)
                        )

            else:
                distances = 2 * np.arcsin(np.sqrt(distances_list) / 2)

            matching_distances_indices = distances <= matching_distance

            # Creating the redshift separation filter to cut the candidates

            z_match_min, z_match_max = self.z_interval(
                use_zerr_match,
                self.match_data[self.match_coordinates["z"]][indices],
                delta_z,
                n_delta_z,
                z_err=(
                    self.match_data[
                        list(
                            filter(
                                lambda key: self.match_properties[key] == "z_err1"
                                or self.match_properties[key] == "z_err2",
                                self.match_properties,
                            )
                        )[0]
                    ][indices]
                    if use_zerr_match == True
                    else 0.0
                ),
            )
            z_query_min, z_query_max = self.z_interval(
                use_zerr_query,
                matched["z"][i],
                delta_z,
                n_delta_z,
                z_err=(
                    matched[
                        self.query_properties[
                            list(
                                filter(
                                    lambda key: self.query_properties[key] == "z_err1"
                                    or self.query_properties[key] == "z_err2",
                                    self.query_properties,
                                )
                            )[0]
                        ]
                    ][i]
                    if use_zerr_query == True
                    else 0.0
                ),
            )

            z_matching_min = z_match_max >= z_query_min
            z_matching_max = z_match_min <= z_query_max

            indices = indices[
                z_matching_min & z_matching_max & matching_distances_indices
            ]
            distances = distances[
                z_matching_min & z_matching_max & matching_distances_indices
            ]

            # Adding the matched coordinates and properties for the multiple matched table

            matched["ID_matched"].append(indices)
            matched["distances"].append(distances)
            matched["RA_matched"].append(
                self.match_data[self.match_coordinates["RA"]][indices]
            )
            matched["DEC_matched"].append(
                self.match_data[self.match_coordinates["DEC"]][indices]
            )
            matched["z_matched"].append(
                self.match_data[self.match_coordinates["z"]][indices]
            )

            if self.match_properties is not None:
                for prop in self.match_properties:
                    matched[self.match_properties.get(prop)].append(
                        self.match_data[prop][indices]
                    )

            # Selecting the best candidate based on the choosen criteria

            if len(indices) > 0:
                matches += 1
                match selection_criteria:
                    case "distances":
                        min_distances = min(distances)
                        index = np.where(distances == min_distances)[0][0]
                    case "redshift_proximity":
                        redshift_proximity = abs(
                            matched["z"][i]
                            - self.match_data[self.match_coordinates["z"]][indices]
                        )
                        min_redshift_proximity = min(redshift_proximity)
                        index = np.where(redshift_proximity == min_redshift_proximity)[
                            0
                        ][0]
                    case "more_massive":

                        if (
                            "mass1" not in self.match_properties.values()
                            and "mass2" not in self.match_properties.values()
                        ):
                            raise Exception(
                                "match_properties must have a mass1 or mass2 column to use this criteria"
                            )
                        else:
                            masses = list(
                                filter(
                                    lambda key: self.match_properties[key] == "mass1"
                                    or self.match_properties[key] == "mass2",
                                    self.match_properties,
                                )
                            )[0]

                            max_mass = max(self.match_data[masses][indices])
                            index = np.where(
                                self.match_data[masses][indices] == max_mass
                            )[0][0]

                    case _:
                        raise Exception(
                            "selection_criteria must be eihter distances , redshift_proximity or more_massive. Was given %s"
                            % (selection_criteria)
                        )

                # Adding the matched coordinates and properties for the best matched table

                best_matched["ID_matched"].append(indices[index])
                best_matched["distances"].append(distances[index])
                best_matched["RA_matched"].append(
                    self.match_data[self.match_coordinates["RA"]][indices[index]]
                )
                best_matched["DEC_matched"].append(
                    self.match_data[self.match_coordinates["DEC"]][indices[index]]
                )
                best_matched["z_matched"].append(
                    self.match_data[self.match_coordinates["z"]][indices[index]]
                )

                if self.match_properties is not None:
                    for prop in self.match_properties:
                        best_matched[self.match_properties.get(prop)].append(
                            self.match_data[prop][indices[index]]
                        )

            else:
                best_matched["ID_matched"].append(indices)
                best_matched["distances"].append(distances)
                best_matched["RA_matched"].append(
                    self.match_data[self.match_coordinates["RA"]][indices]
                )
                best_matched["DEC_matched"].append(
                    self.match_data[self.match_coordinates["DEC"]][indices]
                )
                best_matched["z_matched"].append(
                    self.match_data[self.match_coordinates["z"]][indices]
                )
                if self.match_properties is not None:
                    for prop in self.match_properties:
                        best_matched[self.match_properties.get(prop)].append(
                            self.match_data[prop][indices]
                        )

        # Summary of the number of matches found

        print(
            """
        Number of objects matching: %s
        Number of objects in the target catalog: %s
        Number of unmatched objects: %s
        Number of matched objects: %s
        """
            % (
                len(best_matched["ID"]),
                len(r_m_ones),
                len(best_matched["ID"]) - matches,
                matches,
            )
        )

        return Table(matched), Table(best_matched)

    # ----- Performs the 2d match between the two catalogs in both directions ------------------------------------------------------------------------------#
    # cross_match_2d()
    # NcHicosmo: cosmo:                 cosmology object to convert z in distance
    # float:     matching_distance1:    maximum distance to consider a match in  catalog 1
    # float:     matching_distance2:    maximum distance to consider a match in  catalog 2
    # int:       n_nearest_neighbours1: number of nearest neighbours in catalog 1
    # int:       n_nearest_neighbours2: number of nearest neighbours in catalog 2
    # bool:      verbose:               show the progress of the match
    # bool:      match_dist_3d:         If converts from angular separation to physical distance
    # float:     delta_z1:              Uncertainty associated with the measurament z to match in  catalog 1
    # float:     delta_z2:              Uncertainty associated with the measurament z to match in  catalog 2
    # float:     n_delta_z1:            Number of sigmas we expect to cut the objects in z space in  catalog 1
    # float:     n_delta_z2:            Number of sigmas we expect to cut the objects in z space in  catalog 2
    # str:       selection_criteria1:   selection criteria to choose the best canditate in  catalog 1 can be more_massive , distances, redshift_proximity
    # str:       selection_criteria2:   selection criteria to choose the best canditate in  catalog 2 can be more_massive , distances, redshift_proximity
    # bool:      use_zerr1:             If use the z_err of  catalog 1 instead of delta_z
    # bool:      use_zerr2:             If use the z_err of  catalog 2 instead of delta_z
    # str:       which_radius1:         The radius to be used to convert from angular separation to physical distance  in  catalog 1 can be query_radius ,
    #                                    match_radius , min_radius , max_radius
    # str:       which_radius2:         The radius to be used to convert from angular separation to physical distance  in  catalog 2 can be query_radius ,
    #                                    match_radius , min_radius , max_radius
    # Returns: astropy_table: mult_cat1: table with all candidates of matched objects of catalog1, mult_cat2: table with all candidates of matched objects of catalog2, best_cat1: table with the  best candidate of matched objects of catalog 1, best_cat2: table with the  best candidate of matched objects of catalog 2, cross: table with the objects that are best candidates in both directions
    # ----------------------------------------------------------------------------------------------------------------------------------------------------#
    def cross_match_2d(
        self,
        cosmo: Nc.HICosmo,
        matching_distance1: float,
        matching_distance2: float,
        n_nearest_neighbours1: int,
        n_nearest_neighbours2: int,
        verbose: bool = True,
        match_dist_3d: bool = False,
        delta_z1: float = 0.1,
        delta_z2: float = 0.1,
        n_delta_z1: int = 1,
        n_delta_z2: int = 1,
        selection_criteria1: str = "distances",
        selection_criteria2: str = "distances",
        use_zerr1: bool = False,
        use_zerr2: bool = False,
        which_radius1: str = "query_radius",
        which_radius2: str = "query_radius",
    ) -> Table:
        """Match objects in the sky.

        The function matches objects in the sky using the provided matching distance in both directions using a 2d match.
        """
        cat2 = SkyMatch(
            query_catalog_path=self.match_catalog_path,
            query_coordinates=self.match_coordinates,
            match_catalog_path=self.query_catalog_path,
            match_coordinates=self.query_coordinates,
            query_properties=self.match_properties,
            match_properties=self.query_properties,
        )

        mult_cat1, best_cat1 = self.match_2d(
            cosmo=cosmo,
            matching_distance=matching_distance1,
            n_nearest_neighbours=n_nearest_neighbours1,
            match_dist_3d=match_dist_3d,
            delta_z=delta_z1,
            n_delta_z=n_delta_z1,
            selection_criteria=selection_criteria1,
            use_zerr_query=use_zerr1,
            use_zerr_match=use_zerr2,
            which_radius=which_radius1,
        )

        mult_cat2, best_cat2 = cat2.match_2d(
            cosmo=cosmo,
            matching_distance=matching_distance2,
            n_nearest_neighbours=n_nearest_neighbours2,
            match_dist_3d=match_dist_3d,
            delta_z=delta_z2,
            n_delta_z=n_delta_z2,
            selection_criteria=selection_criteria2,
            use_zerr_query=use_zerr2,
            use_zerr_match=use_zerr1,
            which_radius=which_radius2,
        )

        cross = Table(
            names=(best_cat1.columns),
            dtype=tuple(
                [best_cat1.columns[i].dtype for i in range(len(best_cat1.columns))]
            ),
        )

        for i in tqdm.tqdm(range(len(best_cat1["ID"]))):
            if type(best_cat1["ID_matched"][i]) != np.ndarray:
                if (
                    best_cat1["ID"][i]
                    == best_cat2[best_cat2["ID"] == best_cat1["ID_matched"][i]][
                        "ID_matched"
                    ]
                ):
                    cross.add_row(best_cat1[i])

        print(
            """
        Number of cross matched objects: %s
        """
            % (len(cross["ID"]))
        )

        return mult_cat1, mult_cat2, best_cat1, best_cat2, cross

    # ----- Performs the 3d match between the two catalogs in both directions ------------------------------------------------------------------------------#
    # cross_match_3d()
    # NcHicosmo: cosmo:                 cosmology object to convert z in distance
    # float:     matching_distance1:    maximum distance to consider a match in  catalog 1
    # float:     matching_distance2:    maximum distance to consider a match in  catalog 2
    # int:       n_nearest_neighbours1: number of nearest neighbours in catalog 1
    # int:       n_nearest_neighbours2: number of nearest neighbours in catalog 2
    # bool:      verbose:               show the progress of the match
    # str:       selection_criteria1:   selection criteria to choose the best canditate in  catalog 1 can be more_massive , distances, redshift_proximity
    # str:       selection_criteria2:   selection criteria to choose the best canditate in  catalog 2 can be more_massive , distances, redshift_proximity
    # Returns: astropy_table: mult_cat1: table with all candidates of matched objects of catalog1, mult_cat2: table with all candidates of matched objects of catalog2, best_cat1: table with the  best candidate of matched objects of catalog 1, best_cat2: table with the  best candidate of matched objects of catalog 2, cross: table with the objects that are best candidates in both directions
    # ----------------------------------------------------------------------------------------------------------------------------------------------------#

    def cross_match_3d(
        self,
        cosmo: Nc.HICosmo,
        matching_distance1: float,
        matching_distance2: float,
        n_nearest_neighbours1: int,
        n_nearest_neighbours2: int,
        verbose: bool = True,
        selection_criteria1: str = "distances",
        selection_criteria2: str = "distances",
    ) -> Table:
        """Match objects in the sky.

        The function matches objects in the sky using the provided matching distance in both directions using a 3d match.
        """
        cat2 = SkyMatch(
            query_catalog_path=self.match_catalog_path,
            query_coordinates=self.match_coordinates,
            match_catalog_path=self.query_catalog_path,
            match_coordinates=self.query_coordinates,
            query_properties=self.match_properties,
            match_properties=self.query_properties,
        )

        mult_cat1, best_cat1 = self.match_3d(
            cosmo=cosmo,
            matching_distance=matching_distance1,
            n_nearest_neighbours=n_nearest_neighbours1,
            selection_criteria=selection_criteria1,
        )

        mult_cat2, best_cat2 = cat2.match_3d(
            cosmo=cosmo,
            matching_distance=matching_distance2,
            n_nearest_neighbours=n_nearest_neighbours2,
            selection_criteria=selection_criteria2,
        )

        cross = Table(
            names=(best_cat1.columns),
            dtype=tuple(
                [best_cat1.columns[i].dtype for i in range(len(best_cat1.columns))]
            ),
        )

        for i in tqdm.tqdm(range(len(best_cat1["ID"]))):
            if type(best_cat1["ID_matched"][i]) != np.ndarray:
                if (
                    best_cat1["ID"][i]
                    == best_cat2[best_cat2["ID"] == best_cat1["ID_matched"][i]][
                        "ID_matched"
                    ]
                ):
                    cross.add_row(best_cat1[i])

        print(
            """
        Number of cross matched objects: %s
        """
            % (len(cross["ID"]))
        )

        return mult_cat1, mult_cat2, best_cat1, best_cat2, cross

    # ----- Calculate the interval in redshift to consider the neighbour as possible candidate -----------------=------------------------------------------#
    # z_interval()
    # bool:  use_zerr:  If use the z_err of the match catalog instead of delta_z
    # float: z:         redshift of the object
    # float: delta_z:   Uncertainty associated with the redshift of the object
    # float: n_delta_z: Number of sigmas we expect to cut the objects in z space
    # float: z_err:     redshift error of object

    # Returns: float: z_min: lower bound of the interval , z_max: upper bound of the interval
    # ----------------------------------------------------------------------------------------------------------------------------------------------------#

    def z_interval(
        self,
        use_zerr: bool = False,
        z: float = 0.0,
        delta_z: float = 0.0,
        n_delta_z: int = 0.0,
        z_err: float = 0.0,
    ) -> float:
        """Calculate the interval in redshift to consider the neighbour as possible candidate."""
        if use_zerr == True:
            z_min = z - n_delta_z * z_err
            z_max = z + n_delta_z * z_err

        else:
            z_min = z - delta_z * n_delta_z * (1.0 + z)
            z_max = z + delta_z * n_delta_z * (1.0 + z)

        return z_min, z_max

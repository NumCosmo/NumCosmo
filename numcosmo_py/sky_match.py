"""Matching objects in the sky.

Module to match objects in the sky halo-halo, cluster-halo, cluster-cluster.
"""

from __future__ import annotations
import dataclasses
from typing import TypedDict, cast
from pathlib import Path
from enum import Enum, auto
import numpy as np
import numpy.typing as npt

from astropy.io import fits
from astropy.table import Table
from numcosmo_py import Ncm, Nc
from numcosmo_py.helper import npa_to_seq

Ncm.cfg_init()


class Coordinates(TypedDict, total=False):
    """Coordinates mapping.

    Dictionary to map the coordinates to the name of the columns in the FITS file.

    :param RA: Right Ascension.
    :param DEC: Declination.
    :param z: Redshift.
    """

    RA: str
    DEC: str
    z: str


class SelectionCriteria(str, Enum):
    """Selection criteria for the best candidate."""

    @staticmethod
    def _generate_next_value_(name, _start, _count, _last_values):
        return name.lower()

    DISTANCES = auto()
    REDSHIFT_PROXIMITY = auto()
    MORE_MASSIVE = auto()


class DistanceMethod(str, Enum):
    """Distance method to use in the matching.

    :param ANGULAR_SEPARATION: Angular separation between the objects.
    :param QUERY_RADIUS: 3d distance rescaled using the query object radii.
    :param MATCH_RADIUS: 3d distance rescaled using the match object radii.
    :param MIN_RADIUS: 3d distance rescaled using the minimum of the query and match
        object radii.
    :param MAX_RADIUS: 3d distance rescaled using the maximum of the query and match
        object radii.
    """

    @staticmethod
    def _generate_next_value_(name, _start, _count, _last_values):
        return name.lower()

    ANGULAR_SEPARATION = auto()
    QUERY_RADIUS = auto()
    MATCH_RADIUS = auto()
    MIN_RADIUS = auto()
    MAX_RADIUS = auto()


@dataclasses.dataclass(frozen=True)
class Mask:
    """Sky match mask.

    :param mask: Boolean mask of filtered matches.
    """

    mask: npt.NDArray[np.bool_]

    def __post_init__(self):
        """Check the mask shape."""
        assert len(self.mask.shape) == 2

    def __and__(self, other: Mask) -> Mask:
        """Logical and between two masks."""
        assert self.mask.shape == other.mask.shape
        return Mask(mask=self.mask & other.mask)

    @property
    def shape(self) -> tuple[int, ...]:
        """Shape of the mask."""
        return self.mask.shape

    @property
    def array(self) -> np.ndarray[tuple[int, int], np.dtype[np.bool_]]:
        """Return the mask array."""
        return cast(np.ndarray[tuple[int, int], np.dtype[np.bool_]], self.mask)


@dataclasses.dataclass(frozen=True)
class BestCandidates:
    """Sky match best candidates.

    :param query_filter: Boolean array to filter the query for which we have a best
        candidate.
    :param indices: Array with the best candidates indices.
    """

    query_filter: npt.NDArray[np.bool_]
    indices: npt.NDArray[np.int64]

    def __post_init__(self):
        """Check the shapes."""
        assert len(self.query_filter.shape) == 1
        assert len(self.indices.shape) == 1
        assert np.sum(self.query_filter) == len(self.indices)

    def get_cross_match_indices(self, inverse_best: BestCandidates) -> dict[int, int]:
        """Get the cross match indices."""
        query_indices = np.arange(len(self.query_filter))[self.query_filter]
        match_indices = np.arange(len(inverse_best.query_filter))[
            inverse_best.query_filter
        ]
        query_to_match = dict(zip(query_indices, self.indices))
        match_to_query = dict(zip(match_indices, inverse_best.indices))
        cross = {
            query_i: match_i
            for query_i, match_i in query_to_match.items()
            if match_i in match_to_query and match_to_query[match_i] == query_i
        }
        return cross


def _check_coordinates(table: fits.FITS_rec, coordinates: Coordinates) -> None:
    """Check if the coordinates are provided.

    :param table: fits table to be used with the coordinates.
    :param coordinates: coordinates we need to check to do the matching.

    The function checks if the coordinates are provided in the table.
    """
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


def _load_fits_data(catalog: Path) -> fits.FITS_rec:
    """Load FITS data from the provided catalog path.

    The function returns the first FITS table found in the provided catalog.

    :param Path catalog: Path of the catalog to load.
    :return: data from the fits catalog.
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


class SkyMatchResult:
    """Class to store the results of the sky match."""

    def __init__(
        self,
        sky_match: SkyMatch,
        nearest_neighbours_indices: np.ndarray[tuple[int, int], np.dtype[np.int64]],
        nearest_neighbours_distances: np.ndarray[tuple[int, int], np.dtype[np.float64]],
    ):
        """Initialize the SkyMatchResult class.

        :param sky_match: SkyMatch object used to perform the match.
        :param nearest_neighbours: Array of indices of the nearest neighbours.
        """
        self.sky_match = sky_match
        self.nearest_neighbours_indices = nearest_neighbours_indices
        self.nearest_neighbours_distances = nearest_neighbours_distances

        assert (
            self.nearest_neighbours_indices.shape
            == self.nearest_neighbours_distances.shape
        )

        assert self.nearest_neighbours_indices.shape[0] == len(self.sky_match.query_ra)

    def full_mask(self) -> Mask:
        """Return a mask with all the elements set to True."""
        return Mask(
            np.ones_like(self.nearest_neighbours_distances, dtype=bool),
        )

    def filter_mask_by_distance(
        self,
        max_distance: float,
        mask: Mask | None = None,
    ) -> Mask:
        """Filter the mask by distance.

        :param max_distance: Maximum distance to consider a match.
        """
        if mask is None:
            mask = self.full_mask()
        return mask & Mask(self.nearest_neighbours_distances < max_distance)

    def filter_mask_by_redshift_proximity(
        self,
        sigma_z: float = 0.1,
        n_sigma: float = 1.0,
        mask: Mask | None = None,
        *,
        query_sigma_z_column: str | None = None,
        match_sigma_z_column: str | None = None,
    ) -> Mask:
        """Filter the mask by redshift proximity.

        :param delta_z: Maximum delta_z to consider a match.
        """
        if (
            query_sigma_z_column is not None
            and query_sigma_z_column not in self.sky_match.query_data
        ):
            raise ValueError(f"Column {query_sigma_z_column} not found in query data.")

        if (
            match_sigma_z_column is not None
            and match_sigma_z_column not in self.sky_match.match_data
        ):
            raise ValueError(f"Column {match_sigma_z_column} not found in match data.")

        if mask is None:
            mask = self.full_mask()

        if query_sigma_z_column is None:
            query_delta_z = n_sigma * sigma_z * (1.0 + self.sky_match.query_z)
        else:
            query_delta_z = n_sigma * self.sky_match.query_data[query_sigma_z_column]

        if match_sigma_z_column is None:
            match_delta_z = n_sigma * sigma_z * (1.0 + self.sky_match.match_z)
        else:
            match_delta_z = n_sigma * self.sky_match.match_data[match_sigma_z_column]

        match_z_upper = self.sky_match.match_z + match_delta_z
        match_z_lower = self.sky_match.match_z - match_delta_z

        query_z_upper = self.sky_match.query_z + query_delta_z
        query_z_lower = self.sky_match.query_z - query_delta_z

        z_mask_array = np.array(
            [
                (match_z_upper[indx] >= z_l) & (match_z_lower[indx] <= z_u)
                for indx, z_l, z_u in zip(
                    self.nearest_neighbours_indices, query_z_lower, query_z_upper
                )
            ],
            dtype=bool,
        )

        return mask & Mask(z_mask_array)

    def select_best(
        self,
        selection_criteria: SelectionCriteria = SelectionCriteria.DISTANCES,
        more_massive_column: str | None = None,
        mask: Mask | None = None,
    ) -> BestCandidates:
        """Select the best matched objects.

        :param selection_criteria: Selection criteria to use.
        :param mask: Mask to use to select the best matched objects.
        """
        if selection_criteria == SelectionCriteria.MORE_MASSIVE:
            if more_massive_column is None or (
                more_massive_column not in self.sky_match.match_data
            ):
                raise ValueError(
                    f"A more_massive_column ({more_massive_column}) must "
                    f"be provided and present in the match data."
                )

        if mask is None:
            mask = self.full_mask()
        assert mask.shape == self.nearest_neighbours_distances.shape

        # Filter query objects with at least one match in the mask
        query_filter = np.sum(mask.array, axis=1).astype(bool)
        mask_array = mask.array[query_filter]
        assert mask_array is not None

        filtered_indices = [
            nni[m]
            for nni, m in zip(self.nearest_neighbours_indices[query_filter], mask_array)
        ]
        match selection_criteria:
            case SelectionCriteria.DISTANCES:
                filtered_distances = self.nearest_neighbours_distances[query_filter]
                best_candidates_indices = [
                    nni[np.argmin(distances[m])]
                    for nni, distances, m in zip(
                        filtered_indices, filtered_distances, mask_array
                    )
                ]
            case SelectionCriteria.REDSHIFT_PROXIMITY:
                match_z = self.sky_match.match_z
                query_z = self.sky_match.query_z
                filtered_z = query_z[query_filter]
                best_candidates_indices = [
                    nni[np.argmin(np.abs(match_z[nni] - z))]
                    for nni, z in zip(filtered_indices, filtered_z)
                ]
            case SelectionCriteria.MORE_MASSIVE:
                match_more_massive = self.sky_match.match_data[more_massive_column]
                best_candidates_indices = [
                    nni[np.argmax(match_more_massive[nni])] for nni in filtered_indices
                ]
            case _:
                raise ValueError(
                    "The selection_criteria must be a SelectionCriteria enum value."
                )
        return BestCandidates(
            query_filter=query_filter,
            indices=np.array(best_candidates_indices, dtype=np.int64),
        )

    def _get_by_indices(
        self,
        x: npt.NDArray,
        indices: np.ndarray[tuple[int, int], np.dtype[np.int64]],
        mask: Mask,
    ) -> list[npt.NDArray]:
        assert len(x.shape) == 1
        assert mask.shape == indices.shape
        return [x[i[m]] for i, m in zip(indices, mask.array)]

    def to_table_complete(
        self,
        mask: Mask | None = None,
        *,
        query_properties: dict[str, str] | None = None,
        match_properties: dict[str, str] | None = None,
    ) -> Table:
        """Convert the match result to a complete table.

        The function returns a table with all the properties of the query catalog
        and the properties of the match catalog for the best matched objects.
        """
        table = Table()
        table["ID"] = np.arange(len(self.sky_match.query_ra))
        table["RA"] = self.sky_match.query_ra
        table["DEC"] = self.sky_match.query_dec
        table["z"] = self.sky_match.query_z

        if mask is None:
            mask = self.full_mask()

        if query_properties is not None:
            assert isinstance(query_properties, dict)
            for key, value in query_properties.items():
                table[value] = self.sky_match.query_data[key]

        if match_properties is not None:
            assert isinstance(match_properties, dict)
            for key, value in match_properties.items():
                table[value] = self._get_by_indices(
                    self.sky_match.match_data[key],
                    self.nearest_neighbours_indices,
                    mask,
                )

        table["ID_matched"] = self._get_by_indices(
            np.arange(len(self.sky_match.match_ra)),
            self.nearest_neighbours_indices,
            mask,
        )
        table["distances"] = [
            d[m] for d, m in zip(self.nearest_neighbours_distances, mask.array)
        ]
        table["RA_matched"] = self._get_by_indices(
            self.sky_match.match_ra, self.nearest_neighbours_indices, mask
        )
        table["DEC_matched"] = self._get_by_indices(
            self.sky_match.match_dec, self.nearest_neighbours_indices, mask
        )
        table["z_matched"] = self._get_by_indices(
            self.sky_match.match_z, self.nearest_neighbours_indices, mask
        )
        return table

    def to_table_best(
        self,
        best: BestCandidates,
        *,
        query_properties: dict[str, str] | None = None,
        match_properties: dict[str, str] | None = None,
    ) -> Table:
        """Convert the match result to a table with only the best matched objects."""
        table = Table()

        assert len(best.query_filter) == len(self.sky_match.query_ra)
        query_filter = best.query_filter

        table["ID"] = np.arange(len(self.sky_match.query_ra))[query_filter]
        table["RA"] = self.sky_match.query_ra[query_filter]
        table["DEC"] = self.sky_match.query_dec[query_filter]
        table["z"] = self.sky_match.query_z[query_filter]

        table["ID_matched"] = best.indices
        table["RA_matched"] = self.sky_match.match_ra[best.indices]
        table["DEC_matched"] = self.sky_match.match_dec[best.indices]
        table["z_matched"] = self.sky_match.match_z[best.indices]

        if query_properties is not None:
            for key, value in query_properties.items():
                table[value] = self.sky_match.query_data[key][query_filter]

        if match_properties is not None:
            for key, value in match_properties.items():
                table[value] = self.sky_match.match_data[key][best.indices]

        return table


class SkyMatch:
    """Class to match objects in the sky halo-halo, cluster-halo, cluster-cluster."""

    def __init__(
        self,
        query_data: fits.FITS_rec,
        query_coordinates: Coordinates,
        match_data: fits.FITS_rec,
        match_coordinates: Coordinates,
    ) -> None:
        """Create a new SkyMatch object from FITS data."""
        self.query_data = query_data
        self.match_data = match_data
        _check_coordinates(self.query_data, query_coordinates)
        _check_coordinates(self.match_data, match_coordinates)
        self.query_coordinates = query_coordinates
        self.match_coordinates = match_coordinates

    @classmethod
    def new_from_fits(
        cls,
        query_catalog_path: Path,
        query_coordinates: Coordinates,
        match_catalog_path: Path,
        match_coordinates: Coordinates,
    ) -> SkyMatch:
        """Initialize the class.

        :param query_catalog_path: Path of the catalog to be matched.
        :param query_coordinates: Coordinates of the catalog to be matched.
        :param match_catalog_path: Path of the catalog we are searching the match.
        :param match_coordinates: Coordinates of the catalog we are searching the match.
        """
        query_data = _load_fits_data(query_catalog_path)
        match_data = _load_fits_data(match_catalog_path)
        return cls(query_data, query_coordinates, match_data, match_coordinates)

    def invert_query_match(self) -> SkyMatch:
        """Invert the query and match catalogs."""
        return SkyMatch(
            self.match_data,
            self.match_coordinates,
            self.query_data,
            self.query_coordinates,
        )

    def ra_dec_to_theta_phi(self, ra, dec):
        """Convert RA and DEC to theta and phi.

        :param float ra: right ascencion angle in degrees.
        :param float dec: declination angle in degrees.
        :return: float theta: theta in radiands, float phi: phi in radiands
        """
        ra_rad = np.radians(ra)
        dec_rad = np.radians(dec)

        theta = np.pi / 2 - dec_rad
        phi = ra_rad

        return theta, phi

    @property
    def query_ra(self) -> np.ndarray:
        """Return the RA coordinates of the query catalog."""
        return self.query_data[self.query_coordinates["RA"]]

    @property
    def query_dec(self) -> np.ndarray:
        """Return the DEC coordinates of the query catalog."""
        return self.query_data[self.query_coordinates["DEC"]]

    @property
    def query_z(self) -> np.ndarray:
        """Return the z coordinates of the query catalog."""
        return self.query_data[self.query_coordinates["z"]]

    @property
    def match_ra(self) -> np.ndarray:
        """Return the RA coordinates of the match catalog."""
        return self.match_data[self.match_coordinates["RA"]]

    @property
    def match_dec(self) -> np.ndarray:
        """Return the DEC coordinates of the match catalog."""
        return self.match_data[self.match_coordinates["DEC"]]

    @property
    def match_z(self) -> np.ndarray:
        """Return the z coordinates of the match catalog."""
        return self.match_data[self.match_coordinates["z"]]

    def match_3d(self, cosmo: Nc.HICosmo, n_nearest_neighbours: int) -> SkyMatchResult:
        """Match objects in the sky.

        The function matches objects in the sky using the provided matching distance and
        number of nearest neighbours. It considers the tridimensional space to match the
        objects computing the cosmological distances between the objects.

        :param cosmo: Cosmology object to convert z in distance.
        :param n_nearest_neighbours: Number of nearest neighbours.
        :return: A SkyMatchResult object.
        """
        if ("z" not in self.match_coordinates) or ("z" not in self.query_coordinates):
            raise ValueError(
                "To perform a 3D matching, "
                "the redshift must be provided for both catalogs."
            )

        query_theta, query_phi = self.ra_dec_to_theta_phi(self.query_ra, self.query_dec)
        query_z = self.query_z

        match_ra = self.match_ra
        match_dec = self.match_dec
        match_z = self.match_z
        match_theta, match_phi = self.ra_dec_to_theta_phi(match_ra, match_dec)

        # Preparing the k-nearest-neighbours and cosmology objects for the matching
        snn = Ncm.SphereNN()
        dist = Nc.Distance.new(9.0)
        dist.prepare(cosmo)
        RH_Mpc = cosmo.RH_Mpc()

        match_r = dist.angular_diameter_vector(cosmo, npa_to_seq(match_z))
        query_r = dist.angular_diameter_vector(cosmo, npa_to_seq(query_z))

        snn.insert_array(match_r, match_theta, match_phi)
        snn.rebuild()

        distances_list, indices_list = snn.knn_search_distances_batch(
            query_r, query_theta, query_phi, n_nearest_neighbours
        )
        distances = (np.sqrt(distances_list) * RH_Mpc).reshape(-1, n_nearest_neighbours)
        indices = np.array(indices_list, dtype=int).reshape(-1, n_nearest_neighbours)

        return SkyMatchResult(self, indices, distances)

    def match_2d(
        self,
        cosmo: Nc.HICosmo,
        n_nearest_neighbours: int,
        distance_method: DistanceMethod = DistanceMethod.QUERY_RADIUS,
    ) -> SkyMatchResult:
        """Match objects in the sky.

        The function matches objects in the sky using the provided matching distance.

        :param cosmo: Cosmology object to convert z in distance.
        :param matching_distance: Maximum distance to consider a match.
        :param n_nearest_neighbours: Number of nearest neighbours.
        :param distance_method: Method to compute the distance between the objects.
        :return: astropy_table: matched: table with all candidates of matched objects,
            best_matched: table with the  best candidate of matched objects
        """
        if ("z" not in self.match_coordinates) or ("z" not in self.query_coordinates):
            raise ValueError(
                "To perform a  matching, "
                "the redshift must be provided for both catalogs."
            )

        # Preparing the k-nearest-neighbours and cosmology objects for the matching

        query_theta, query_phi = self.ra_dec_to_theta_phi(self.query_ra, self.query_dec)
        match_theta, match_phi = self.ra_dec_to_theta_phi(self.match_ra, self.match_dec)

        r_m_ones = np.ones_like(match_theta)
        r_q_ones = np.ones_like(query_theta)

        snn = Ncm.SphereNN()
        dist = Nc.Distance.new(3.0)
        dist.prepare(cosmo)
        RH_Mpc = cosmo.RH_Mpc()

        snn.insert_array(r_m_ones, match_theta, match_phi)
        snn.rebuild()

        distances_list, indices_list = snn.knn_search_distances_batch(
            r_q_ones, query_theta, query_phi, n_nearest_neighbours
        )
        indices = np.array(indices_list, dtype=int).reshape(-1, n_nearest_neighbours)

        match distance_method:
            case DistanceMethod.ANGULAR_SEPARATION:
                distances = (2 * np.arcsin(np.sqrt(distances_list) / 2)).reshape(
                    -1, n_nearest_neighbours
                )
            case DistanceMethod.QUERY_RADIUS:
                distances = (np.sqrt(distances_list) * RH_Mpc).reshape(
                    -1, n_nearest_neighbours
                )
                query_r = np.array(
                    dist.angular_diameter_vector(cosmo, npa_to_seq(self.query_z))
                )
                # We multiply the distances by the angular diameter distance of the
                # query object to rescale the distances to physical distances.
                distances = cast(
                    np.ndarray[tuple[int, int], np.dtype[np.float64]],
                    query_r.reshape(-1, 1) * distances,
                )
            case DistanceMethod.MATCH_RADIUS:
                distances = (np.sqrt(distances_list) * RH_Mpc).reshape(
                    -1, n_nearest_neighbours
                )
                match_r = np.array(
                    dist.angular_diameter_vector(cosmo, npa_to_seq(self.match_z))
                )
                # We multiply the distances by the angular diameter distance of the
                # match object to rescale the distances to physical distances.
                distances = cast(
                    np.ndarray[tuple[int, int], np.dtype[np.float64]],
                    match_r[indices] * distances,
                )
            case DistanceMethod.MIN_RADIUS:
                distances = (np.sqrt(distances_list) * RH_Mpc).reshape(
                    -1, n_nearest_neighbours
                )
                query_r = np.array(
                    dist.angular_diameter_vector(cosmo, npa_to_seq(self.query_z))
                )
                query_distances = query_r.reshape(-1, 1) * distances
                match_r = np.array(
                    dist.angular_diameter_vector(cosmo, npa_to_seq(self.match_z))
                )
                match_distances = match_r[indices] * distances
                distances = cast(
                    np.ndarray[tuple[int, int], np.dtype[np.float64]],
                    np.minimum(query_distances, match_distances),
                )
            case DistanceMethod.MAX_RADIUS:
                distances = (np.sqrt(distances_list) * RH_Mpc).reshape(
                    -1, n_nearest_neighbours
                )
                query_r = np.array(
                    dist.angular_diameter_vector(cosmo, npa_to_seq(self.query_z))
                )
                query_distances = query_r.reshape(-1, 1) * distances
                match_r = np.array(
                    dist.angular_diameter_vector(cosmo, npa_to_seq(self.match_z))
                )
                match_distances = match_r[indices] * distances
                distances = cast(
                    np.ndarray[tuple[int, int], np.dtype[np.float64]],
                    np.maximum(query_distances, match_distances),
                )
            case _:
                raise ValueError(
                    "The distance method must be a DistanceMethod enum value."
                )

        return SkyMatchResult(self, indices, distances)

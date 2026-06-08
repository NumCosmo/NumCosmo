"""Matching objects in the sky.

Module to match objects in the sky halo-halo, cluster-halo, cluster-cluster.
"""

from __future__ import annotations
import dataclasses
from typing import TypedDict, cast
from pathlib import Path
from enum import StrEnum, auto
from typing_extensions import assert_never
import numpy as np
import numpy.typing as npt

from astropy.table import Table, Column, join
import networkx as nx

from numcosmo_py import Ncm, Nc

Ncm.cfg_init()


class Coordinates(TypedDict, total=False):
    """Coordinates mapping.

    Dictionary to map the coordinates to the name of the columns in the
    astropy.table.Table.

    :param RA: Right Ascension.
    :param DEC: Declination.
    :param z: Redshift.
    """

    RA: str
    DEC: str
    z: str


class IDs(TypedDict, total=False):
    """ID mapping.

    Dictionary to map the IDs to the names of the columns in the
    astropy.table.Table.

    :param ID: Object ID.
    :param MemberID: Member ID.
    :param pmem: Probability of membership.
    """

    ID: str
    MemberID: str
    pmem: str


class SelectionCriteria(StrEnum):
    """Selection criteria for the best candidate."""

    DISTANCES = auto()
    REDSHIFT_PROXIMITY = auto()
    MORE_MASSIVE = auto()


class SharedFractionMethod(StrEnum):
    """Method to evaluate shared fraction in the matching.

    :param NO_PMEM: Shared fraction computed as the number of shared members
        divided by the number of members of the object.
    :param QUERY_PMEM: Shared fraction computed as the sum of probability membership of
        the shared members divided by the probability of membership of the query
        object.
    :param MATCH_PMEM: Shared fraction computed as the sum of probability membership of
        the shared members divided by the probability of membership of the match
        object.
    :param PMEM: Shared fraction computed as the sum of probability membership of the
        shared members divided by the sum of probability membership of the query and
        match objects.
    """

    @staticmethod
    def _generate_next_value_(name, _start, _count, _last_values):
        return name.lower()

    NO_PMEM = auto()
    MATCH_PMEM = auto()
    QUERY_PMEM = auto()
    PMEM = auto()


class DistanceMethod(StrEnum):
    """Distance method to use in the matching.

    :param ANGULAR_SEPARATION: Angular separation between the objects.
    :param QUERY_RADIUS: 3d distance rescaled using the query object radii.
    :param MATCH_RADIUS: 3d distance rescaled using the match object radii.
    :param MIN_RADIUS: 3d distance rescaled using the minimum of the query and match
        object radii.
    :param MAX_RADIUS: 3d distance rescaled using the maximum of the query and match
        object radii.
    """

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

    def __invert__(self) -> Mask:
        """Logical not of the mask."""
        return Mask(mask=~self.mask)

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

    @property
    def query_indices(self) -> npt.NDArray[np.int64]:
        """Query indices."""
        return np.arange(len(self.query_filter))[self.query_filter]

    @property
    def query_match_dict(self) -> dict[int, int]:
        """Query match dictionary."""
        return dict(zip(self.query_indices, self.indices))

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


def _check_coordinates(table: Table, coordinates: Coordinates) -> None:
    """Check if the coordinates are provided.

    :param table: fits table to be used with the coordinates.
    :param coordinates: coordinates we need to check to do the matching.

    The function checks if the coordinates are provided in the table.
    """
    if ("RA" not in coordinates) or ("DEC" not in coordinates):
        raise ValueError("RA and DEC coordinates must be provided.")

    if (coordinates["RA"] not in table.columns) or (
        coordinates["DEC"] not in table.columns
    ):
        raise ValueError(
            f"RA and DEC coordinates mapped by {coordinates} "
            f"not found in the provided catalog {table.columns}."
        )

    if "z" in coordinates and coordinates["z"] not in table.columns:
        raise ValueError(
            f"Redshift coordinate mapped by {coordinates} "
            f"not found in the provided catalog {table.columns}."
        )


def _check_ID(table: Table, ids: IDs | None) -> None:
    """Check if the IDs are provided.

    :param table: fits table to be used with the IDs.
    :param ids: IDs we need to check to do the matching.

    The function checks if the IDs are provided in the table.
    """
    if ids is None or ("ID" not in ids) or ("MemberID" not in ids):
        raise ValueError("ID and MemberID must be provided.")

    if (ids["ID"] not in table.columns) or (ids["MemberID"] not in table.columns):
        raise ValueError(
            f"ID and MemberID coordinates mapped by {ids} "
            f"not found in the provided catalog {table.columns}."
        )


def _load_fits_data(catalog: Path) -> Table:
    """Load FITS data from the provided catalog path.

    The function returns the first FITS table found in the provided catalog.

    :param Path catalog: Path of the catalog to load.
    :return: data from the fits catalog.
    """
    return Table.read(catalog.as_posix(), format="fits")


class SkyMatchResult:
    """Result of a distance based sky match (:meth:`SkyMatch.match_2d` / ``match_3d``).

    Holds, for each of the ``n_query`` query objects, the row indices of its
    ``k`` nearest neighbours in the match catalog together with the corresponding
    distances. Both arrays are dense and have shape ``(n_query, k)``.
    """

    def __init__(
        self,
        sky_match: SkyMatch,
        nearest_neighbours_indices: np.ndarray[tuple[int, int], np.dtype[np.int64]],
        nearest_neighbours_distances: np.ndarray[tuple[int, int], np.dtype[np.float64]],
    ):
        """Initialize the SkyMatchResult class.

        :param sky_match: SkyMatch object used to perform the match.
        :param nearest_neighbours_indices: Array of indices of the nearest neighbours.
        :param nearest_neighbours_distances: Array of distances to the nearest
            neighbours.
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
        if query_sigma_z_column is not None and (
            query_sigma_z_column not in self.sky_match.query_data.columns
        ):
            raise ValueError(f"Column {query_sigma_z_column} not found in query data.")

        if match_sigma_z_column is not None and (
            match_sigma_z_column not in self.sky_match.match_data.columns
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
                more_massive_column not in self.sky_match.match_data.columns
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
            case _ as unreachable:  # pragma: no cover
                assert_never(unreachable)

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
            for key, value in query_properties.items():
                table[value] = self.sky_match.query_data[key]

        if match_properties is not None:
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


class SkyMatchIDResult:
    """Result of a membership based sky match (:meth:`SkyMatch.match_ID`).

    Unlike :class:`SkyMatchResult`, the candidate lists are ragged: each query
    object may share members with an arbitrary number of match objects (possibly
    none). For every query row we therefore store a variable length array of the
    matched object IDs (as they appear in the match catalog) and the corresponding
    linking coefficients.
    """

    def __init__(
        self,
        sky_match: SkyMatch,
        candidate_ids: npt.NDArray[np.object_],
        candidate_coefficients: npt.NDArray[np.object_],
    ):
        """Initialize the SkyMatchIDResult class.

        :param sky_match: SkyMatch object used to perform the match.
        :param candidate_ids: Object array, one entry per query row, each holding the
            match-catalog object IDs sharing members with that query object.
        :param candidate_coefficients: Object array, one entry per query row, each
            holding the linking coefficients aligned with ``candidate_ids``.
        """
        self.sky_match = sky_match
        self.candidate_ids = candidate_ids
        self.candidate_coefficients = candidate_coefficients

        assert candidate_ids.shape == candidate_coefficients.shape
        assert candidate_ids.shape[0] == len(self.sky_match.query_ra)

    def _match_id_to_row(self) -> dict[int, int]:
        """Map match-catalog object IDs to their row index in the match catalog."""
        return {mid: i for i, mid in enumerate(self.sky_match.match_id)}

    def select_best(self) -> BestCandidates:
        """Resolve a one-to-one matching maximizing the total linking coefficient.

        The candidate graph (query objects on one side, match objects on the other,
        edges weighted by the linking coefficient) is solved with a maximum weight
        matching. The returned :class:`BestCandidates` holds, for every query object
        that keeps a match, the row index of its match in the match catalog.
        """
        match_id_to_row = self._match_id_to_row()

        graph = nx.Graph()
        for query_row, (ids, coeffs) in enumerate(
            zip(self.candidate_ids, self.candidate_coefficients)
        ):
            for match_id, coeff in zip(ids, coeffs):
                graph.add_edge(
                    ("query", query_row),
                    ("match", match_id_to_row[match_id]),
                    weight=float(coeff),
                )

        matching = nx.max_weight_matching(graph, maxcardinality=False)

        query_to_match_row: dict[int, int] = {}
        for node_a, node_b in matching:
            query_node, match_node = (
                (node_a, node_b) if node_a[0] == "query" else (node_b, node_a)
            )
            query_to_match_row[query_node[1]] = match_node[1]

        n_query = len(self.sky_match.query_ra)
        query_filter = np.zeros(n_query, dtype=bool)
        indices = []
        for query_row in sorted(query_to_match_row):
            query_filter[query_row] = True
            indices.append(query_to_match_row[query_row])

        return BestCandidates(
            query_filter=query_filter,
            indices=np.array(indices, dtype=np.int64),
        )

    def _coefficient_for(self, query_row: int, match_id) -> float:
        """Return the linking coefficient stored for a given query/match pair."""
        ids = list(self.candidate_ids[query_row])
        return float(self.candidate_coefficients[query_row][ids.index(match_id)])

    def to_table_complete(
        self,
        *,
        query_properties: dict[str, str] | None = None,
        match_properties: dict[str, str] | None = None,
    ) -> Table:
        """Convert the match result to a complete table.

        The table has one row per query object, with the full (ragged) list of
        candidate matches and their linking coefficients.
        """
        match_id_to_row = self._match_id_to_row()
        candidate_rows = [
            np.array([match_id_to_row[mid] for mid in ids], dtype=np.int64)
            for ids in self.candidate_ids
        ]

        table = Table()
        table["query_id"] = np.asarray(self.sky_match.query_id)
        table["RA"] = self.sky_match.query_ra
        table["DEC"] = self.sky_match.query_dec
        table["z"] = self.sky_match.query_z

        if query_properties is not None:
            for key, value in query_properties.items():
                table[value] = self.sky_match.query_data[key]

        table["match_id"] = [np.asarray(ids) for ids in self.candidate_ids]
        table["linking_coefficient"] = [
            np.asarray(coeffs, dtype=float) for coeffs in self.candidate_coefficients
        ]

        match_ra = np.asarray(self.sky_match.match_ra)
        match_dec = np.asarray(self.sky_match.match_dec)
        match_z = np.asarray(self.sky_match.match_z)
        table["RA_matched"] = [match_ra[rows] for rows in candidate_rows]
        table["DEC_matched"] = [match_dec[rows] for rows in candidate_rows]
        table["z_matched"] = [match_z[rows] for rows in candidate_rows]

        if match_properties is not None:
            for key, value in match_properties.items():
                column = np.asarray(self.sky_match.match_data[key])
                table[value] = [column[rows] for rows in candidate_rows]

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

        table["query_id"] = np.asarray(self.sky_match.query_id)[query_filter]
        table["RA"] = np.asarray(self.sky_match.query_ra)[query_filter]
        table["DEC"] = np.asarray(self.sky_match.query_dec)[query_filter]
        table["z"] = np.asarray(self.sky_match.query_z)[query_filter]

        match_id = np.asarray(self.sky_match.match_id)
        table["match_id"] = match_id[best.indices]
        table["linking_coefficient"] = [
            self._coefficient_for(query_row, match_id[match_row])
            for query_row, match_row in best.query_match_dict.items()
        ]
        table["RA_matched"] = np.asarray(self.sky_match.match_ra)[best.indices]
        table["DEC_matched"] = np.asarray(self.sky_match.match_dec)[best.indices]
        table["z_matched"] = np.asarray(self.sky_match.match_z)[best.indices]

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
        query_data: Table,
        query_coordinates: Coordinates,
        match_data: Table,
        match_coordinates: Coordinates,
        query_member_data: Table | None = None,
        query_ids: IDs | None = None,
        match_member_data: Table | None = None,
        match_ids: IDs | None = None,
    ) -> None:
        """Create a new SkyMatch object from an astropy.table.Table."""
        self.query_data = query_data
        self.match_data = match_data
        self.query_member_data = query_member_data
        self.match_member_data = match_member_data
        self.query_ids = query_ids
        self.match_ids = match_ids

        _check_coordinates(self.query_data, query_coordinates)
        _check_coordinates(self.match_data, match_coordinates)
        self.query_coordinates = query_coordinates
        self.match_coordinates = match_coordinates

        if query_member_data is not None and match_member_data is not None:
            _check_ID(query_member_data, query_ids)
            _check_ID(match_member_data, match_ids)
            assert query_ids is not None and match_ids is not None

            query_consistent = len(
                np.unique(query_member_data[query_ids["ID"]])
            ) == len(np.unique(self.query_data[query_ids["ID"]]))
            match_consistent = len(
                np.unique(match_member_data[match_ids["ID"]])
            ) == len(np.unique(self.match_data[match_ids["ID"]]))
            if not (query_consistent and match_consistent):
                raise ValueError(
                    "The number of unique 'ID' in object and member "
                    "catalogs must be the same."
                )

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
            self.match_member_data,
            self.match_ids,
            self.query_member_data,
            self.query_ids,
        )

    @staticmethod
    def ra_dec_to_theta_phi(
        ra: npt.ArrayLike, dec: npt.ArrayLike
    ) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
        """Convert RA and DEC to theta and phi.

        :param ra: right ascension angle in degrees.
        :param dec: declination angle in degrees.
        :return: theta in radians, phi in radians.
        """
        ra_rad = np.radians(ra)
        dec_rad = np.radians(dec)

        theta = np.pi / 2 - dec_rad
        phi = ra_rad

        return theta, phi

    @staticmethod
    def object_count(id_column: npt.ArrayLike) -> Table:
        """Return the number of occurrences for each unique object ID.

        :param id_column: The array or column containing the object IDs to be counted.
        :return: A table containing two columns:
        - ``id``: The unique object IDs found in the input.
        - ``nmem``: The count of occurrences for each ID.
        """
        ids, counts = np.unique(id_column, return_counts=True)

        return Table({"id": ids, "nmem": counts})

    @property
    def query_ra(self) -> Column:
        """Return the RA coordinates of the query catalog."""
        return self.query_data[self.query_coordinates["RA"]]

    @property
    def query_dec(self) -> Column:
        """Return the DEC coordinates of the query catalog."""
        return self.query_data[self.query_coordinates["DEC"]]

    @property
    def query_z(self) -> Column:
        """Return the z coordinates of the query catalog."""
        return self.query_data[self.query_coordinates["z"]]

    @property
    def match_ra(self) -> Column:
        """Return the RA coordinates of the match catalog."""
        return self.match_data[self.match_coordinates["RA"]]

    @property
    def match_dec(self) -> Column:
        """Return the DEC coordinates of the match catalog."""
        return self.match_data[self.match_coordinates["DEC"]]

    @property
    def match_z(self) -> Column:
        """Return the z coordinates of the match catalog."""
        return self.match_data[self.match_coordinates["z"]]

    @property
    def match_id(self) -> Column:
        """Return the ID column of the match catalog."""
        assert self.match_ids is not None
        return self.match_data[self.match_ids["ID"]]

    @property
    def query_id(self) -> Column:
        """Return the ID column of the query catalog."""
        assert self.query_ids is not None
        return self.query_data[self.query_ids["ID"]]

    @property
    def query_pmem(self) -> Column:
        """Return the pmem values of the query catalog."""
        assert self.query_ids is not None
        return self.query_data[self.query_ids["pmem"]]

    @property
    def match_pmem(self) -> Column:
        """Return the pmem values of the match catalog."""
        assert self.match_ids is not None
        return self.match_data[self.match_ids["pmem"]]

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

        # Decide which distance is most appropriate to 3d matching
        match_r = dist.angular_diameter_array(cosmo, match_z)
        query_r = dist.angular_diameter_array(cosmo, query_z)

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
                "To perform a matching, "
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
                query_r = np.array(dist.angular_diameter_array(cosmo, self.query_z))
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
                match_r = np.array(dist.angular_diameter_array(cosmo, self.match_z))
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
                query_r = np.array(dist.angular_diameter_array(cosmo, self.query_z))
                query_distances = query_r.reshape(-1, 1) * distances
                match_r = np.array(dist.angular_diameter_array(cosmo, self.match_z))
                match_distances = match_r[indices] * distances
                distances = cast(
                    np.ndarray[tuple[int, int], np.dtype[np.float64]],
                    np.minimum(query_distances, match_distances),
                )
            case DistanceMethod.MAX_RADIUS:
                distances = (np.sqrt(distances_list) * RH_Mpc).reshape(
                    -1, n_nearest_neighbours
                )
                query_r = np.array(dist.angular_diameter_array(cosmo, self.query_z))
                query_distances = query_r.reshape(-1, 1) * distances
                match_r = np.array(dist.angular_diameter_array(cosmo, self.match_z))
                match_distances = match_r[indices] * distances
                distances = cast(
                    np.ndarray[tuple[int, int], np.dtype[np.float64]],
                    np.maximum(query_distances, match_distances),
                )
            case _ as unreachable:  # pragma: no cover
                assert_never(unreachable)

        return SkyMatchResult(self, indices, distances)

    def match_ID(
        self,
        use_shared_fraction: bool = False,
        shared_fraction_method: SharedFractionMethod = SharedFractionMethod.NO_PMEM,
    ) -> SkyMatchIDResult:
        """Match objects in the sky using shared members.

        The function matches objects by their members: two objects are candidate
        matches when they share at least one member. The strength of each candidate
        match is summarized by a linking coefficient.

        :param use_shared_fraction: If True, compute the linking coefficient from the
            shared member fractions (see ``shared_fraction_method``); otherwise use the
            Jaccard index of the member sets.
        :param shared_fraction_method: Method to compute the shared fraction of members.
        :return: A SkyMatchIDResult object.
        """
        assert self.query_member_data is not None
        assert self.match_member_data is not None
        assert self.query_ids is not None
        assert self.match_ids is not None

        # Preparing the member catalogs for the matching

        query_table = self.query_member_data.copy()
        query_table.rename_column(self.query_ids["ID"], "query_id")
        query_table.rename_column(self.query_ids["MemberID"], "MemberID")

        match_table = self.match_member_data.copy()
        match_table.rename_column(self.match_ids["ID"], "match_id")
        match_table.rename_column(self.match_ids["MemberID"], "MemberID")

        if "pmem" in query_table.colnames:
            query_table.rename_column("pmem", "pmem_query")
        if "pmem" in match_table.colnames:
            match_table.rename_column("pmem", "pmem_match")

        # Matching the catalogs by MemberID

        matched_catalog = join(
            query_table, match_table, keys="MemberID", join_type="inner"
        )

        matched_catalog_grouped = matched_catalog.group_by(["query_id", "match_id"])

        all_combinations = (
            matched_catalog_grouped.groups.keys
        )  # Table with the multiple matched query_id-match_id pairs

        # Counting 'number of members' for each object in the catalogs

        nmem_query = self.object_count(query_table["query_id"])
        nmem_query.rename_columns(["id", "nmem"], ["query_id", "nmem_query"])

        nmem_match = self.object_count(match_table["match_id"])
        nmem_match.rename_columns(["id", "nmem"], ["match_id", "nmem_match"])

        # Adding the nmem_query and nmem_match columns to all_combinations table

        all_combinations = join(
            all_combinations, nmem_query, keys="query_id", join_type="left"
        )
        all_combinations = join(
            all_combinations, nmem_match, keys="match_id", join_type="left"
        )

        # Counting number of shared members

        indices = matched_catalog_grouped.groups.indices
        counts = np.diff(indices)
        shared_count = matched_catalog_grouped.groups.keys
        shared_count["shared_count"] = counts

        all_combinations = join(
            all_combinations,
            shared_count,
            keys=["query_id", "match_id"],
            join_type="left",
        )

        # Shared fraction

        fraction_query = None
        fraction_match = None

        match shared_fraction_method:
            case SharedFractionMethod.NO_PMEM:
                fraction_query = (
                    all_combinations["shared_count"] / all_combinations["nmem_query"]
                )
                fraction_match = (
                    all_combinations["shared_count"] / all_combinations["nmem_match"]
                )

            case SharedFractionMethod.QUERY_PMEM:
                if "pmem_query" not in matched_catalog.colnames:
                    raise ValueError(
                        "To perform a matching, pmem column must "
                        "be provided for the query catalog."
                    )

                all_combinations["sum_shared_pmem"] = matched_catalog_grouped[
                    "query_id", "match_id", "pmem_query"
                ].groups.aggregate(np.sum)["pmem_query"]

                query_group = query_table.group_by(["query_id"])

                total_pmem = query_group["query_id", "pmem_query"].groups.aggregate(
                    np.sum
                )
                total_pmem.rename_column("pmem_query", "sum_total_pmem")

                all_combinations = join(all_combinations, total_pmem, keys="query_id")

                fraction_query = (
                    all_combinations["sum_shared_pmem"]
                    / all_combinations["sum_total_pmem"]
                )
                fraction_match = (
                    all_combinations["shared_count"] / all_combinations["nmem_match"]
                )

            case SharedFractionMethod.MATCH_PMEM:
                if "pmem_match" not in matched_catalog.colnames:
                    raise ValueError(
                        "To perform a matching, pmem column must "
                        "be provided for the match catalog."
                    )

                all_combinations["sum_shared_pmem"] = matched_catalog_grouped[
                    "query_id", "match_id", "pmem_match"
                ].groups.aggregate(np.sum)["pmem_match"]

                match_group = match_table.group_by(["match_id"])
                total_pmem = match_group["match_id", "pmem_match"].groups.aggregate(
                    np.sum
                )
                total_pmem.rename_column("pmem_match", "sum_total_pmem")

                all_combinations = join(all_combinations, total_pmem, keys="match_id")

                fraction_query = (
                    all_combinations["shared_count"] / all_combinations["nmem_query"]
                )
                fraction_match = (
                    all_combinations["sum_shared_pmem"]
                    / all_combinations["sum_total_pmem"]
                )

            case SharedFractionMethod.PMEM:
                if (
                    "pmem_query" not in matched_catalog.colnames
                    or "pmem_match" not in matched_catalog.colnames
                ):
                    raise ValueError(
                        "To perform a matching, pmem column must "
                        "be provided for both catalogs."
                    )

                all_combinations["sum_shared_pmem_query"] = matched_catalog_grouped[
                    "query_id", "match_id", "pmem_query"
                ].groups.aggregate(np.sum)["pmem_query"]
                all_combinations["sum_shared_pmem_match"] = matched_catalog_grouped[
                    "query_id", "match_id", "pmem_match"
                ].groups.aggregate(np.sum)["pmem_match"]

                query_group = query_table.group_by(["query_id"])
                sum_total_pmem_query = query_group[
                    "query_id", "pmem_query"
                ].groups.aggregate(np.sum)
                sum_total_pmem_query.rename_column("pmem_query", "sum_total_pmem_query")

                match_group = match_table.group_by(["match_id"])
                sum_total_pmem_match = match_group[
                    "match_id", "pmem_match"
                ].groups.aggregate(np.sum)
                sum_total_pmem_match.rename_column("pmem_match", "sum_total_pmem_match")

                all_combinations = join(
                    all_combinations, sum_total_pmem_query, keys="query_id"
                )
                all_combinations = join(
                    all_combinations, sum_total_pmem_match, keys="match_id"
                )

                fraction_query = (
                    all_combinations["sum_shared_pmem_query"]
                    / all_combinations["sum_total_pmem_query"]
                )
                fraction_match = (
                    all_combinations["sum_shared_pmem_match"]
                    / all_combinations["sum_total_pmem_match"]
                )

            case _ as unreachable:
                assert_never(unreachable)

        # Linking coefficient

        if use_shared_fraction:
            all_combinations["linking_coefficient"] = (
                fraction_query * (fraction_query + fraction_match) / 2
            )
        else:
            all_combinations["linking_coefficient"] = all_combinations[
                "shared_count"
            ] / (
                all_combinations["nmem_query"]
                + all_combinations["nmem_match"]
                - all_combinations["shared_count"]
            )

        # Structured result: gather candidate matches per query object

        grouped_comb = all_combinations.group_by("query_id")

        comb_dict = {}
        for group in grouped_comb.groups:
            qid = group["query_id"][0]
            comb_dict[qid] = (
                list(group["match_id"]),
                list(group["linking_coefficient"]),
            )

        candidate_ids = np.empty(len(self.query_data), dtype=object)
        candidate_coefficients = np.empty(len(self.query_data), dtype=object)
        for i, qid in enumerate(self.query_data[self.query_ids["ID"]]):
            ids, coeffs = comb_dict.get(qid, ([], []))
            candidate_ids[i] = np.array(ids)
            candidate_coefficients[i] = np.array(coeffs, dtype=float)

        return SkyMatchIDResult(self, candidate_ids, candidate_coefficients)

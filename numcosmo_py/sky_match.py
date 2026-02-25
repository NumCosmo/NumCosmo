"""Matching objects in the sky.

Module to match objects in the sky halo-halo, cluster-halo, cluster-cluster.
"""

from __future__ import annotations
import dataclasses
from typing import TypedDict, cast
from pathlib import Path
from enum import Enum, auto
from typing_extensions import assert_never
import numpy as np
import numpy.typing as npt

from astropy.table import Table, join
from numcosmo_py import Ncm, Nc
from numcosmo_py.helper import npa_to_seq

import networkx as nx
import random

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
    """

    ID: str
    MemberID: str
    pmem: str
    
class MatchingType(str,Enum):
    """Matching type.

    :param DISTANCE: Match by distance.
    :param ID: Match by ID.
    """

    @staticmethod
    def _generate_next_value_(name, _start, _count, _last_values):
        return name.lower()

    DISTANCE = auto()
    ID = auto()
    


class SelectionCriteria(str, Enum):
    """Selection criteria for the best candidate."""

    @staticmethod
    def _generate_next_value_(name, _start, _count, _last_values):
        return name.lower()

    DISTANCES = auto()
    REDSHIFT_PROXIMITY = auto()
    MORE_MASSIVE = auto()

class SharedFractionMethod(str, Enum):
    """Method to evaluate shared fraction in the matching.

    :param NO_PMEM: Shared fraction computed as the number of shared members
        divided by the number of members of the object.
    :param QUERY_PMEM: Shared fraction computed as the sum of probability membership of the shared members
        divided by the probability of membership of the query object.
    :param MATCH_PMEM: Shared fraction computed as the sum of probability membership of the shared members
        divided by the probability of membership of the match object.
    :param PMEM: Shared fraction computed as the sum of probability membership of the shared members
        divided by the sum of probability membership of the query and match objects.
    """

    @staticmethod
    def _generate_next_value_(name, _start, _count, _last_values):
        return name.lower()

    NO_PMEM = auto()
    MATCH_PMEM = auto()
    QUERY_PMEM = auto()
    PMEM = auto()

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

def _check_ID(table: Table, ids: IDs) -> None:
    """Check if the IDs are provided.

    :param table: fits table to be used with the IDs.
    :param ids: IDs we need to check to do the matching.

    The function checks if the IDs are provided in the table.
    """
    if ("ID" not in ids) or ("MemberID" not in ids):
        raise ValueError("ID and MemberID must be provided.")

    if (ids["ID"] not in table.columns) or (
        ids["MemberID"] not in table.columns
    ):
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
    """Class to store the results of the sky match."""

    def __init__(
        self,
        sky_match: SkyMatch,
        matching_type: MatchingType,
        nearest_neighbours_indices: np.ndarray[tuple[int, int], np.dtype[np.int64]],
        nearest_neighbours_distances: np.ndarray[tuple[int, int], np.dtype[np.float64]],
    ):
        """Initialize the SkyMatchResult class.

        :param sky_match: SkyMatch object used to perform the match.
        :param matching_type: Type of matching performed.
        :param nearest_neighbours: Array of indices of the nearest neighbours.
        """
        self.sky_match = sky_match
        self.matching_type = matching_type
        self.nearest_neighbours_indices = nearest_neighbours_indices
        self.nearest_neighbours_distances = nearest_neighbours_distances

        assert (
            self.nearest_neighbours_indices.shape
            == self.nearest_neighbours_distances.shape
        )

        if self.matching_type == MatchingType.DISTANCE:
            assert self.nearest_neighbours_indices.shape[0] == len(self.sky_match.query_ra)

    def full_mask(self) -> Mask:
        """Return a mask with all the elements set to True."""
        
        if self.matching_type == MatchingType.ID:
            raise ValueError("Full mask is not available for ID matching.")
        
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
        if self.matching_type == MatchingType.ID:
            raise ValueError("Filtering by distance is not available for ID matching.")
        
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
        if self.matching_type == MatchingType.ID:
            raise ValueError("Filtering by redshift proximity is not available for ID matching.")

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

    
    def _select_best_distance(
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
    
    
    def _select_best_id(
        self,
        ) -> BestCandidates:
       
        matched_query = np.repeat(self.sky_match.query_id, [len(i) for i in self.nearest_neighbours_indices])
        matched_match = np.concatenate(self.nearest_neighbours_indices)
        coefficients = np.concatenate(self.nearest_neighbours_distances)
                
        matches = Table(
            [matched_query, matched_match, coefficients],
            names=('query_id', 'match_id', 'linking_coeff')
        )
        
        matches = matches[matches['match_id'] != -1]
       
        G = nx.Graph()
        for id1, id2, w in matches:
            G.add_edge(("query_id",id1), ("match_id",id2), weight=w)
        
        random.seed(42)
        total_weight = 0.0
        global_matching = set()
        
        for i, nodes in enumerate(nx.connected_components(G)):
            
            subG = G.subgraph(nodes)
            m = nx.max_weight_matching(subG, maxcardinality=False) 
            w = sum(subG[u][v]['weight'] for u, v in m)

            total_weight += w
            global_matching |= m 

              
        query_id_unique = []
        match_id_unique = []
        linking_coef = []
        
        for u, v in global_matching:
            
            # garantir lado correto
            if u[0] == "query_id":
                node1, node2 = u, v
            else:
                node1, node2 = v, u
        
            id1 = node1[1]
            id2 = node2[1]
            weight = G[node1][node2]['weight']
        
            query_id_unique.append(id1)
            match_id_unique.append(id2)
            linking_coef.append(weight)
        
        best = Table()
        best['query_id'] = query_id_unique
        best['match_id'] = match_id_unique
        best['linking_coeff'] = linking_coef   


        mapping = dict(zip(best['query_id'], best['match_id']))

        combinations = Table(
            [self.sky_match.query_id, self.nearest_neighbours_indices],
            names=('query_id', 'matches')
        )

        combinations['matches'] = [ mapping.get(id_val, False) for id_val in combinations['query_id'] ]

        query_filter = (combinations['matches'] != False)

        valid_matches = list(combinations['matches'][query_filter])

        match_to_index = {val: i for i, val in enumerate(self.sky_match.match_id)}

        best_candidates_indices = np.array([
            match_to_index[m] for m in valid_matches if m in match_to_index
        ], dtype=np.int64)
        
            
        return BestCandidates(
            query_filter=query_filter,
            indices=np.array(best_candidates_indices, dtype=np.int64),
        )
    

    def select_best(
        self,
        selection_criteria: SelectionCriteria | None = None,
        more_massive_column: str | None = None,
        mask: Mask | None = None,
    ) -> BestCandidates:
        
        if self.matching_type == MatchingType.ID:
            return self._select_best_id()
        
        if self.matching_type == MatchingType.DISTANCE:
            return self._select_best_distance(selection_criteria, more_massive_column, mask)


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

    def to_table_best_id(
        self,
        query_properties: dict[str, str] | None = None,
        match_properties: dict[str, str] | None = None,
        ) -> Table:
        """Convert the match result to a table with only the best matched objects."""

        best = self.select_best()     
        
        query_table = Table()
        
        query_table["ID"] = self.sky_match.query_id
        query_table["RA"] =  self.sky_match.query_ra
        query_table["DEC"] = self.sky_match.query_dec
        query_table["z"] = self.sky_match.query_z

        match_table = Table()

        query_table["ID_matched"] = self.sky_match.query_id
        query_table["RA_matched"] =  self.sky_match.query_ra
        query_table["DEC_matched"] = self.sky_match.query_dec
        query_table["z_matched"] = self.sky_match.query_z

        
        table = Table()

        table["ID"] = self.sky_match.query_id
        table["RA"] =  self.sky_match.query_ra
        table["DEC"] = self.sky_match.query_dec
        table["z"] = self.sky_match.query_z

        table["ID_matched"] = best['match_id']
        table["RA_matched"] = best["query_id"]
        table["DEC_matched"] = best["query_id"]
        table["z_matched"] = best["query_id"]

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
        match_member_data: Table | None = None,
        query_ids: IDs | None = None,
        match_ids: IDs | None = None,
    ) -> None:
        """Create a new SkyMatch object from an astropy.table.Table."""
        self.query_data = query_data
        self.match_data = match_data
        self.query_member_data = query_member_data
        self.match_member_data = match_member_data
        # if query_ids or match_id is not None:
        #     _check_ID(self.query_data, query_ids)
        #     _check_ID(self.match_data, match_ids)
        if query_member_data is not None and match_member_data is not None:
            _check_ID(self.query_member_data, query_ids)
            _check_ID(self.match_member_data, match_ids)
        self.query_ids = query_ids
        self.match_ids = match_ids
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

    def object_count(self, id_column: np.ndarray) -> Table:
        """Return the number of occurrences for each unique object ID.

        :param id_column: The array or column containing the object IDs to be counted.
        :return: A table containing two columns:
        - ``id``: The unique object IDs found in the input.
        - ``nmem``: The count of occurrences for each ID.
        """

        ids, counts = np.unique(id_column, return_counts=True)

        return Table({'id': ids, 'nmem': counts})

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
    
    @property
    def match_id(self) -> np.ndarray:
        """Return the ID coordinates of the match catalog."""
        return self.match_data[self.match_ids["ID"]]

    @property
    def query_id(self) -> np.ndarray:
        """Return the ID coordinates of the query catalog."""
        return self.query_data[self.query_ids["ID"]]
    
    @property
    def query_pmem(self) -> np.ndarray:
        """Return the PMEM values of the query catalog."""
        return self.query_data[self.query_ids["PMEM"]]
    
    @property
    def match_pmem(self) -> np.ndarray:
        """Return the PMEM values of the match catalog."""
        return self.match_data[self.match_ids["PMEM"]]
    



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

        match_r = dist.angular_diameter_array(cosmo, npa_to_seq(match_z))
        query_r = dist.angular_diameter_array(cosmo, npa_to_seq(query_z))

        snn.insert_array(match_r, match_theta, match_phi)
        snn.rebuild()

        distances_list, indices_list = snn.knn_search_distances_batch(
            query_r, query_theta, query_phi, n_nearest_neighbours
        )
        distances = (np.sqrt(distances_list) * RH_Mpc).reshape(-1, n_nearest_neighbours)
        indices = np.array(indices_list, dtype=int).reshape(-1, n_nearest_neighbours)

        return SkyMatchResult(self, MatchingType.DISTANCE, indices, distances)

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
                query_r = np.array(
                    dist.angular_diameter_array(cosmo, npa_to_seq(self.query_z))
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
                    dist.angular_diameter_array(cosmo, npa_to_seq(self.match_z))
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
                    dist.angular_diameter_array(cosmo, npa_to_seq(self.query_z))
                )
                query_distances = query_r.reshape(-1, 1) * distances
                match_r = np.array(
                    dist.angular_diameter_array(cosmo, npa_to_seq(self.match_z))
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
                    dist.angular_diameter_array(cosmo, npa_to_seq(self.query_z))
                )
                query_distances = query_r.reshape(-1, 1) * distances
                match_r = np.array(
                    dist.angular_diameter_array(cosmo, npa_to_seq(self.match_z))
                )
                match_distances = match_r[indices] * distances
                distances = cast(
                    np.ndarray[tuple[int, int], np.dtype[np.float64]],
                    np.maximum(query_distances, match_distances),
                )
            case _ as unreachable:  # pragma: no cover
                assert_never(unreachable)

        return SkyMatchResult(self, MatchingType.DISTANCE, indices, distances)
    

    def match_ID(
        self,
        shared_fraction_method: SharedFractionMethod = SharedFractionMethod.NO_PMEM,
    ) -> SkyMatchResult:
        """Match objects in the sky.

        The function matches objects in the sky using the provided ids.

        :param shared_fraction_method: Method to compute the shared fraction of members.
        :return: astropy_table: matched: table with all candidates of matched objects,
        :best_matched: table with the best candidate of matched objects
        """
         
        # Preparing the member catalogs for the matching

        query_table = self.query_member_data.copy()        
        query_table.rename_column(self.query_ids['ID'], 'query_id')
        query_table.rename_column(self.query_ids['MemberID'], 'MemberID')
        
        match_table = self.match_member_data.copy()
        match_table.rename_column(self.match_ids['ID'], 'match_id')
        match_table.rename_column(self.match_ids['MemberID'], 'MemberID')

        if 'pmem' in query_table.colnames:
            query_table.rename_column('pmem', 'pmem_query')
        if 'pmem' in match_table.colnames:
            match_table.rename_column('pmem', 'pmem_match')

       
        # Matching the catalogs by MemberID

        matched_catalog = join(query_table, match_table, keys='MemberID', join_type='inner')
        
        matched_catalog_grouped = matched_catalog.group_by(['query_id', 'match_id'])
        
        all_combinations = matched_catalog_grouped.groups.keys # Table with the multiple matched query_id-match_id pairs

        
        # Evaluating 'number of members' for each object in the catalogs
        
        nmem_query = self.object_count(query_table['query_id'])
        nmem_query.rename_columns(['id', 'nmem'], ['query_id', 'nmem_query'])

        nmem_match = self.object_count(match_table['match_id'])
        nmem_match.rename_columns(['id', 'nmem'], ['match_id', 'nmem_match'])          
        
        
        # Adding the nmem_query and nmem_match columns to all_combinations table 

        all_combinations = join(all_combinations, nmem_query, keys='query_id', join_type='left')
        all_combinations = join(all_combinations, nmem_match, keys='match_id', join_type='left')

        
        # Evaluating number of shared members
        
        matched_pairs = np.column_stack((matched_catalog['query_id'], matched_catalog['match_id']))
        
        unique_pairs, counts = np.unique(matched_pairs, axis=0, return_counts=True)
        
        all_combinations['shared_count'] = counts

        
        # Shared fraction
        
        fraction_query = None
        fraction_match = None

        match shared_fraction_method:
            case SharedFractionMethod.NO_PMEM:
                fraction_query = all_combinations['shared_count'] / all_combinations['nmem_query']
                fraction_match = all_combinations['shared_count'] / all_combinations['nmem_match']

            case SharedFractionMethod.QUERY_PMEM:
                if "pmem_query" not in matched_catalog.colnames:
                    raise ValueError("To perform a matching, pmem column must be provided for the query catalog.")
                
                all_combinations['sum_shared_pmem'] = matched_catalog_grouped.groups.aggregate(np.sum)['pmem_query']

                query_group = query_table.group_by(['query_id'])

                total_pmem = query_group.groups.aggregate(np.sum)['query_id', 'pmem_query']
                total_pmem.rename_column('pmem_query', 'sum_total_pmem')
                
                all_combinations = join(all_combinations, total_pmem, keys='query_id')
                
                fraction_query = all_combinations['sum_shared_pmem'] / all_combinations['sum_total_pmem']
                fraction_match = all_combinations['shared_count'] / all_combinations['nmem_match']
                
            case SharedFractionMethod.MATCH_PMEM:
                if "pmem_match" not in matched_catalog.colnames:
                    raise ValueError("To perform a matching, pmem column must be provided for the match catalog.")

                all_combinations['sum_shared_pmem'] = matched_catalog_grouped.groups.aggregate(np.sum)['pmem_match']

                match_group = match_table.group_by(['match_id'])
                total_pmem = match_group.groups.aggregate(np.sum)['match_id', 'pmem_match']
                total_pmem.rename_column('pmem_match', 'sum_total_pmem')
                
                all_combinations = join(all_combinations, total_pmem, keys='match_id')

                fraction_query = all_combinations['shared_count'] / all_combinations['nmem_query']
                fraction_match = all_combinations['sum_shared_pmem'] / all_combinations['sum_total_pmem']

            case SharedFractionMethod.PMEM:
                if "pmem_query" not in matched_catalog.colnames or "pmem_match" not in matched_catalog.colnames:
                    raise ValueError("To perform a matching, pmem column must be provided for both catalogs.")
                
                all_combinations['sum_shared_pmem_query'] = matched_catalog_grouped.groups.aggregate(np.sum)['pmem_query']
                all_combinations['sum_shared_pmem_match'] = matched_catalog_grouped.groups.aggregate(np.sum)['pmem_match']
                
                query_group = query_table.group_by(['query_id'])
                sum_total_pmem_query = query_group.groups.aggregate(np.sum)['query_id', 'pmem_query']
                sum_total_pmem_query.rename_column('pmem_query', 'sum_total_pmem_query')
                
                match_group = match_table.group_by(['match_id'])
                sum_total_pmem_match = match_group.groups.aggregate(np.sum)['match_id', 'pmem_match']
                sum_total_pmem_match.rename_column('pmem_match', 'sum_total_pmem_match')
                
                all_combinations = join(all_combinations, sum_total_pmem_query, keys='query_id')
                all_combinations = join(all_combinations, sum_total_pmem_match, keys='match_id')
                
                fraction_query = all_combinations['sum_shared_pmem_query'] / all_combinations['sum_total_pmem_query']
                fraction_match = all_combinations['sum_shared_pmem_match'] / all_combinations['sum_total_pmem_match']

            case _ as unreachable:
                assert_never(unreachable)

        
        # Linking coefficient
        all_combinations["linking_coefficient"] = fraction_query * (fraction_query + fraction_match) / 2

        
        # Making a dictionary with the struture: 
        # {query_id: list(matched objects), list(linking coefficients)}
        
        grouped_comb = all_combinations.group_by('query_id')
        
        comb_dict = {}
        for group in grouped_comb.groups:
            qid = group['query_id'][0]
            comb_dict[qid] = (list(group['match_id']), list(group['linking_coefficient']))

        ids = self.query_data[self.query_ids['ID']]

        match_list, coeff_list = zip(*[comb_dict.get(qid, ([-1], [0])) for qid in ids])
        
        match_list = np.array(match_list, dtype=object)
        coeff_list = np.array(coeff_list, dtype=object)
      
        return SkyMatchResult(self, MatchingType.ID, match_list, coeff_list)
                
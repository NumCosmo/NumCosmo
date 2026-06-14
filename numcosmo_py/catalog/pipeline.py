"""End-to-end mock catalog pipeline built on the C generators.

This module assembles the existing NumCosmo C generators into a single,
linked mock: :class:`~numcosmo_py.Nc.HaloCatalogGenerator` draws the cluster
population (with sky positions and radius enabled) and
:class:`~numcosmo_py.Nc.HaloCatalogMemberGenerator` populates each cluster with
galaxy members through a halo occupation distribution. The orchestration layer
is intentionally thin: the heavy lifting lives in the C objects, and this only
wires them together and converts the results to astropy tables.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from astropy.table import Table

from numcosmo_py import Nc, Ncm

from .table import catalog_to_table

__all__ = ["MockCatalogs", "MockPipeline"]


@dataclass
class MockCatalogs:
    """A cluster catalog and its linked galaxy-member catalog.

    The two tables are linked by id: ``members["parent_id"]`` equals
    ``clusters["cluster_id"]`` (the cluster row index), so the members of a
    cluster are ``members[members["parent_id"] == clusters["cluster_id"][i]]``.
    """

    clusters: Table
    members: Table


class MockPipeline:
    """Chains the C catalog generators into a linked cluster + member mock.

    Holds the generation collaborators by composition: a cluster generator built
    from a :class:`~numcosmo_py.Nc.ClusterAbundance` over a
    :class:`~numcosmo_py.Ncm.SkyFootprint`, and a member generator driven by a
    :class:`~numcosmo_py.Nc.GalaxyHOD`. A single :meth:`generate` call produces
    both catalogs as astropy tables.
    """

    def __init__(
        self,
        abundance: Nc.ClusterAbundance,
        hod: Nc.GalaxyHOD,
        footprint: Ncm.SkyFootprint,
        *,
        distance: Nc.Distance | None = None,
    ) -> None:
        """Create a pipeline from its generation collaborators.

        :param abundance: the cluster abundance model driving the counts.
        :param hod: the halo occupation distribution placing member galaxies.
        :param footprint: the sky footprint over which positions are sampled.
        :param distance: optional distance object reused for the member angular
            placement; one is built on demand when omitted.
        """
        self.cluster_generator = Nc.HaloCatalogGenerator.new(abundance)
        self.cluster_generator.set_footprint(footprint)
        self.cluster_generator.set_with_radius(True)

        self.member_generator = Nc.HaloCatalogMemberGenerator.new(hod)
        if distance is not None:
            self.member_generator.set_distance(distance)

    def generate(self, mset: Ncm.MSet, rng: Ncm.RNG) -> MockCatalogs:
        """Generate a linked cluster and member mock.

        :param mset: model set holding the cosmology and the cluster
            redshift/mass proxies.
        :param rng: the random number generator.
        :return: the cluster and member catalogs as astropy tables.
        """
        clusters_hcat = self.cluster_generator.generate(mset, rng)
        members_hcat = self.member_generator.generate(clusters_hcat, mset, rng)
        # The cluster catalog always carries the required columns, so member
        # generation cannot fail here.
        assert members_hcat is not None

        clusters = catalog_to_table(clusters_hcat)
        # The member parent ids are the cluster row indices; expose them as an
        # explicit join key on the cluster side.
        clusters["cluster_id"] = np.arange(len(clusters))

        members = catalog_to_table(members_hcat)

        return MockCatalogs(clusters=clusters, members=members)

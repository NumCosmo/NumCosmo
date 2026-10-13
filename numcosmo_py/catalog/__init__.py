"""Catalog tools for NumCosmo.

This subpackage groups catalog-level utilities: matching objects in the sky
(:mod:`~numcosmo_py.catalog.sky_match`) and generating mock catalogs of halos,
clusters and galaxy members (:mod:`~numcosmo_py.catalog.mock`).
"""

from .confusion import (
    CatalogType,
    calculate_catalog_metrics,
    calculate_split_metrics,
    get_ratios,
)
from .mock import (
    CompletenessModel,
    ConstantCompleteness,
    ConstantPurity,
    MockGenerator,
    PurityModel,
    identity_scaling_relation,
)
from .pipeline import (
    MockCatalogs,
    MockPipeline,
)
from .sky_match import (
    BestCandidates,
    Coordinates,
    DistanceMethod,
    IDs,
    Mask,
    SelectionCriteria,
    SharedFractionMethod,
    SkyMatch,
    SkyMatchIDResult,
    SkyMatchResult,
)
from .table import (
    catalog_from_table,
    catalog_to_table,
)

__all__ = [
    "BestCandidates",
    "CatalogType",
    "CompletenessModel",
    "ConstantCompleteness",
    "ConstantPurity",
    "Coordinates",
    "DistanceMethod",
    "IDs",
    "Mask",
    "MockCatalogs",
    "MockGenerator",
    "MockPipeline",
    "PurityModel",
    "SelectionCriteria",
    "SharedFractionMethod",
    "SkyMatch",
    "SkyMatchIDResult",
    "SkyMatchResult",
    "calculate_catalog_metrics",
    "calculate_split_metrics",
    "catalog_from_table",
    "catalog_to_table",
    "get_ratios",
    "identity_scaling_relation",
]

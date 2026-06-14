"""Catalog tools for NumCosmo.

This subpackage groups catalog-level utilities: matching objects in the sky
(:mod:`~numcosmo_py.catalog.sky_match`) and generating mock catalogs of halos,
clusters and galaxy members (:mod:`~numcosmo_py.catalog.mock`).
"""

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
from .mock import (
    CompletenessModel,
    ConstantCompleteness,
    ConstantPurity,
    MockGenerator,
    PurityModel,
    identity_scaling_relation,
)
from .confusion import (
    CatalogType,
    calculate_catalog_metrics,
    calculate_split_metrics,
    get_ratios,
)
from .table import (
    catalog_from_table,
    catalog_to_table,
)
from .pipeline import (
    MockCatalogs,
    MockPipeline,
)

__all__ = [
    "BestCandidates",
    "Coordinates",
    "DistanceMethod",
    "IDs",
    "Mask",
    "SelectionCriteria",
    "SharedFractionMethod",
    "SkyMatch",
    "SkyMatchIDResult",
    "SkyMatchResult",
    "CompletenessModel",
    "ConstantCompleteness",
    "ConstantPurity",
    "MockGenerator",
    "PurityModel",
    "identity_scaling_relation",
    "CatalogType",
    "calculate_catalog_metrics",
    "calculate_split_metrics",
    "get_ratios",
    "catalog_from_table",
    "catalog_to_table",
    "MockCatalogs",
    "MockPipeline",
]

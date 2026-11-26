"""Columns of :class:`~gi.repository.NumCosmoMath.Catalog` as numpy arrays.

:class:`NcmCatalog` stores each column by its type
(:class:`~gi.repository.NumCosmoMath.CatalogColType`): the DOUBLE columns in the matrix
of ``peek_data()``, the INT columns as 64-bit integers and the BOOL columns as booleans.
These functions convert catalogs to and from :class:`astropy.table.Table`, reading and
writing each column through the storage of its type
(:func:`numcosmo_py.helper.catalog_column`, :func:`numcosmo_py.helper.fill_catalog`)
and inferring the types from the table column dtypes on the way in.
"""

from __future__ import annotations

import numpy as np
from astropy.table import Table

from numcosmo_py import Ncm
from numcosmo_py.helper import catalog_column, col_type_of, fill_catalog


def catalog_to_table(catalog: Ncm.Catalog) -> Table:
    """Convert an :class:`NcmCatalog` to an astropy :class:`~astropy.table.Table`.

    :param catalog: the catalog to convert.
    :return: a table whose columns carry the dtype implied by each column's
        :class:`~gi.repository.NumCosmoMath.CatalogColType`.
    """
    table = Table()

    for col in catalog.peek_columns():
        table[col] = catalog_column(catalog, col)

    return table


def catalog_from_table(table: Table) -> Ncm.Catalog:
    """Build an :class:`NcmCatalog` from an astropy :class:`~astropy.table.Table`.

    Column types are inferred from the table column dtypes: boolean columns map to
    :attr:`CatalogColType.BOOL`, integer columns to :attr:`CatalogColType.INT`, and
    everything else to :attr:`CatalogColType.DOUBLE`.

    :param table: the table to convert.
    :return: a new catalog holding the table data.
    """
    columns = list(table.colnames)
    arrays = {col: np.asarray(table[col]) for col in columns}
    col_types = [col_type_of(arrays[col]) for col in columns]

    catalog = Ncm.Catalog.new_full(len(table), columns, col_types)
    fill_catalog(catalog, arrays)

    return catalog

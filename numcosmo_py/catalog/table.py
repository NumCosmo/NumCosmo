"""Conversions between :class:`~gi.repository.NumCosmoMath.Catalog` and astropy.

:class:`NcmCatalog` stores every column in a single double matrix and tags each
column with a logical type (:class:`~gi.repository.NumCosmoMath.CatalogColType`).
These helpers move data to and from :class:`astropy.table.Table`, using the type
tags to restore the intended ``int``/``bool``/``float`` dtypes on the way out and
to infer them from the table column dtypes on the way in.
"""

from __future__ import annotations

import numpy as np
from astropy.table import Table

from numcosmo_py import Ncm


def catalog_to_table(catalog: Ncm.Catalog) -> Table:
    """Convert an :class:`NcmCatalog` to an astropy :class:`~astropy.table.Table`.

    :param catalog: the catalog to convert.
    :return: a table whose columns carry the dtype implied by each column's
        :class:`~gi.repository.NumCosmoMath.CatalogColType`.
    """
    columns = list(catalog.peek_columns())
    data = catalog.peek_data().to_numpy()
    table = Table()

    for j, col in enumerate(columns):
        col_type = catalog.get_col_type(col)
        column = data[:, j]

        if col_type == Ncm.CatalogColType.INT:
            table[col] = column.astype(np.int64)
        elif col_type == Ncm.CatalogColType.BOOL:
            table[col] = column.astype(bool)
        else:
            table[col] = column.astype(np.float64)

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
    nrows = len(table)
    col_types = []
    data = np.empty((nrows, len(columns)), dtype=np.float64)

    for j, col in enumerate(columns):
        column = np.asarray(table[col])

        if np.issubdtype(column.dtype, np.bool_):
            col_types.append(Ncm.CatalogColType.BOOL)
        elif np.issubdtype(column.dtype, np.integer):
            col_types.append(Ncm.CatalogColType.INT)
        else:
            col_types.append(Ncm.CatalogColType.DOUBLE)

        data[:, j] = column.astype(np.float64)

    catalog = Ncm.Catalog.new_full(nrows, columns, col_types)
    catalog.peek_data().set_from_array(data.reshape(-1).tolist())

    return catalog

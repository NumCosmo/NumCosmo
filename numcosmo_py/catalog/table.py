"""Columns of :class:`~gi.repository.NumCosmoMath.Catalog` as numpy arrays.

:class:`NcmCatalog` stores each column by its type
(:class:`~gi.repository.NumCosmoMath.CatalogColType`): the DOUBLE columns in the matrix
of ``peek_data()``, the INT columns as 64-bit integers and the BOOL columns as booleans.
These helpers read and write columns through the storage of their type, and convert
catalogs to and from :class:`astropy.table.Table`, inferring the types from the table
column dtypes on the way in.
"""

from __future__ import annotations

from collections.abc import Mapping

import numpy as np
from astropy.table import Table

from numcosmo_py import Ncm


def catalog_column(catalog: Ncm.Catalog, col: str) -> np.ndarray:
    """Return the column ``col`` of ``catalog`` as an array of the dtype of its type.

    :param catalog: the catalog.
    :param col: the column name.
    :return: ``float64`` for DOUBLE, ``int64`` for INT and ``bool`` for BOOL columns.
    """
    col_type = catalog.get_col_type(col)

    if col_type == Ncm.CatalogColType.INT:
        return np.array(catalog.get_int_column(col), dtype=np.int64)
    if col_type == Ncm.CatalogColType.BOOL:
        return np.array(catalog.get_bool_column(col), dtype=bool)

    _, j = catalog.get_data_index(col)
    return catalog.peek_data().to_numpy()[:, j].copy()


def fill_catalog(catalog: Ncm.Catalog, columns: Mapping[str, np.ndarray]) -> None:
    """Set every column of ``catalog`` from ``columns``, each through its type.

    :param catalog: the catalog.
    :param columns: one array per column of ``catalog``, keyed by column name, with
        one value per row.
    """
    names = list(catalog.peek_columns())
    nrows = catalog.len()
    double_names = [
        col for col in names if catalog.get_col_type(col) == Ncm.CatalogColType.DOUBLE
    ]

    if double_names and nrows > 0:
        data = np.column_stack(
            [np.asarray(columns[col], dtype=np.float64) for col in double_names]
        )
        catalog.peek_data().set_from_array(data.reshape(-1).tolist())
        catalog.mark_outdated()

    for col in names:
        col_type = catalog.get_col_type(col)
        if col_type == Ncm.CatalogColType.INT:
            catalog.set_int_column(
                col, np.asarray(columns[col], dtype=np.int64).tolist()
            )
        elif col_type == Ncm.CatalogColType.BOOL:
            catalog.set_bool_column(col, np.asarray(columns[col], dtype=bool).tolist())


def col_type_of(column: np.ndarray) -> Ncm.CatalogColType:
    """Return the catalog column type of a numpy column dtype.

    :param column: the column.
    :return: BOOL for booleans, INT for integers and DOUBLE otherwise.
    """
    if np.issubdtype(column.dtype, np.bool_):
        return Ncm.CatalogColType.BOOL
    if np.issubdtype(column.dtype, np.integer):
        return Ncm.CatalogColType.INT
    return Ncm.CatalogColType.DOUBLE


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

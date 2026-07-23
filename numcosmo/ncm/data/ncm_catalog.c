/***************************************************************************
 *            ncm_catalog.c
 *
 *  Sat Jun 13 14:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_catalog.c
 * Copyright (C) 2026 Sandro Dias Pinto Vitenti <vitenti@uel.br>
 *
 * numcosmo is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * numcosmo is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * NcmCatalog:
 *
 * Lightweight named-column catalog.
 *
 * A generic, lightweight container storing a #NcmMatrix of double values together
 * with a name for each column. It provides column-name based access to the data and
 * is meant to be used as a derivable base class for catalog-like objects (for
 * example galaxy and halo catalogs), as well as on its own as a simple table that
 * converts cheaply to/from external table formats.
 *
 * Each column also carries a logical type (#NcmCatalogColType): all values are
 * physically stored as doubles, but a column may be tagged as integer or boolean
 * so the intended dtype can be restored when converting to/from external table
 * formats (integers are exact below 2^53, booleans are stored as 0.0/1.0). String
 * and other heterogeneous columns are intentionally not handled here and are left
 * to specialized subclasses or future additive extensions.
 *
 * This is a pure data container with no cosmological content. It is unrelated to
 * #NcmMSetCatalog, which stores Monte Carlo chains.
 *
 * A catalog may also carry an optional #NcmVarDict of free-form scalar metadata
 * (see ncm_catalog_set_meta() / ncm_catalog_peek_meta()), e.g. provenance
 * information about the whole catalog that does not belong in any single row
 * (a survey field's known cluster redshift and position, say). It is %NULL
 * unless explicitly set.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"
#include "ncm/data/ncm_catalog.h"
#include "ncm/algebra/ncm_matrix.h"
#include "ncm/core/ncm_cfg.h"

/* *INDENT-OFF* */
G_DEFINE_QUARK (ncm-catalog-error, ncm_catalog_error) 
/* *INDENT-ON* */

typedef struct _NcmCatalogPrivate
{
  NcmMatrix *data;
  GStrv columns;
  GHashTable *columns_hash;
  GArray *col_types;
  guint ncols;
  guint len;
  NcmVarDict *meta;
} NcmCatalogPrivate;

enum
{
  PROP_0,
  PROP_DATA,
  PROP_COLUMNS,
  PROP_COL_TYPES,
  PROP_LEN,
  PROP_META,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcmCatalog, ncm_catalog, G_TYPE_OBJECT)

static void
ncm_catalog_init (NcmCatalog *catalog)
{
  NcmCatalogPrivate * const self = ncm_catalog_get_instance_private (catalog);

  self->data         = NULL;
  self->columns      = NULL;
  self->columns_hash = g_hash_table_new_full (g_str_hash, g_str_equal, g_free, g_free);
  self->col_types    = NULL;
  self->ncols        = 0;
  self->len          = 0;
  self->meta         = NULL;
}

static void _ncm_catalog_set_data (NcmCatalogPrivate *self, NcmMatrix *data);

static void
_ncm_catalog_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmCatalog *catalog            = NCM_CATALOG (object);
  NcmCatalogPrivate * const self = ncm_catalog_get_instance_private (catalog);

  switch (prop_id)
  {
    case PROP_DATA:
      _ncm_catalog_set_data (self, g_value_get_object (value));
      break;
    case PROP_COLUMNS:
      g_assert_null (self->columns);
      self->columns = g_value_dup_boxed (value);
      break;
    case PROP_COL_TYPES:
    {
      GVariant *var = g_value_get_variant (value);

      g_clear_pointer (&self->col_types, g_array_unref);

      if (var != NULL)
      {
        self->col_types = g_array_new (FALSE, TRUE, sizeof (guint32));
        ncm_cfg_array_set_variant (self->col_types, var);
      }

      break;
    }
    case PROP_LEN:
      self->len = g_value_get_uint (value);
      break;
    case PROP_META:
      ncm_catalog_set_meta (catalog, g_value_get_boxed (value));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_catalog_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmCatalog *catalog            = NCM_CATALOG (object);
  NcmCatalogPrivate * const self = ncm_catalog_get_instance_private (catalog);

  switch (prop_id)
  {
    case PROP_DATA:
      g_value_set_object (value, self->data);
      break;
    case PROP_COLUMNS:
      g_value_set_boxed (value, self->columns);
      break;
    case PROP_COL_TYPES:
      g_value_take_variant (value, ncm_cfg_array_to_variant (self->col_types, G_VARIANT_TYPE ("u")));
      break;
    case PROP_LEN:
      g_value_set_uint (value, ncm_catalog_len (catalog));
      break;
    case PROP_META:
      g_value_set_boxed (value, self->meta);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_catalog_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (ncm_catalog_parent_class)->constructed (object);
  {
    NcmCatalog *catalog            = NCM_CATALOG (object);
    NcmCatalogPrivate * const self = ncm_catalog_get_instance_private (catalog);
    guint i;

    g_assert_nonnull (self->columns);
    self->ncols = g_strv_length (self->columns);

    for (i = 0; i < self->ncols; i++)
    {
      const gchar *col = self->columns[i];
      guint *index     = g_new0 (guint, 1);

      *index = i;

      g_assert (col != NULL);
      g_assert (!g_hash_table_contains (self->columns_hash, col));
      g_hash_table_insert (self->columns_hash, g_strdup (col), index);
    }

    if (self->col_types == NULL)
    {
      self->col_types = g_array_sized_new (FALSE, TRUE, sizeof (guint32), self->ncols);
      g_array_set_size (self->col_types, self->ncols);
      /* zero-initialized => NCM_CATALOG_COL_TYPE_DOUBLE for every column. */
    }
    else
    {
      g_assert_cmpuint (self->col_types->len, ==, self->ncols);

      for (i = 0; i < self->ncols; i++)
        g_assert_cmpuint (g_array_index (self->col_types, guint32, i), <, NCM_CATALOG_COL_TYPE_LEN);
    }

    if (self->data == NULL)
    {
      self->data = ncm_matrix_new (self->len, self->ncols);
      ncm_matrix_set_zero (self->data);
    }
  }
}

static void
_ncm_catalog_dispose (GObject *object)
{
  NcmCatalog *catalog            = NCM_CATALOG (object);
  NcmCatalogPrivate * const self = ncm_catalog_get_instance_private (catalog);

  ncm_matrix_clear (&self->data);
  g_clear_pointer (&self->columns_hash, g_hash_table_unref);
  g_clear_pointer (&self->columns, g_strfreev);
  g_clear_pointer (&self->col_types, g_array_unref);
  g_clear_pointer (&self->meta, ncm_var_dict_unref);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_catalog_parent_class)->dispose (object);
}

static void
_ncm_catalog_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_catalog_parent_class)->finalize (object);
}

static void
ncm_catalog_class_init (NcmCatalogClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = &_ncm_catalog_set_property;
  object_class->get_property = &_ncm_catalog_get_property;
  object_class->constructed  = &_ncm_catalog_constructed;
  object_class->dispose      = &_ncm_catalog_dispose;
  object_class->finalize     = &_ncm_catalog_finalize;

  /**
   * NcmCatalog:data:
   *
   * The catalog data matrix (rows are entries, columns are named fields).
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_DATA,
                                   g_param_spec_object ("data",
                                                        "Data",
                                                        "Catalog data matrix",
                                                        NCM_TYPE_MATRIX,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_STRINGS));

  /**
   * NcmCatalog:columns:
   *
   * The catalog column names.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_COLUMNS,
                                   g_param_spec_boxed ("columns",
                                                       "Columns",
                                                       "Catalog column names",
                                                       G_TYPE_STRV,
                                                       G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_STRINGS));

  /**
   * NcmCatalog:col-types:
   *
   * The per-column logical types (#NcmCatalogColType values), as an array of
   * unsigned integers parallel to #NcmCatalog:columns. When unset every column
   * defaults to %NCM_CATALOG_COL_TYPE_DOUBLE.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_COL_TYPES,
                                   g_param_spec_variant ("col-types",
                                                         "Column types",
                                                         "Per-column logical types",
                                                         G_VARIANT_TYPE ("au"), NULL,
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_STRINGS));

  /**
   * NcmCatalog:len:
   *
   * The number of rows (entries) in the catalog.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_LEN,
                                   g_param_spec_uint ("len",
                                                      "Length",
                                                      "Number of catalog rows",
                                                      0, G_MAXUINT, 0,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_STRINGS));

  /**
   * NcmCatalog:meta:
   *
   * Free-form scalar metadata about the catalog as a whole (e.g. provenance
   * information that does not belong in any single row). %NULL unless
   * explicitly set.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_META,
                                   g_param_spec_boxed ("meta",
                                                       "Metadata",
                                                       "Catalog-wide free-form metadata",
                                                       NCM_TYPE_VAR_DICT,
                                                       G_PARAM_READWRITE | G_PARAM_STATIC_STRINGS));
}

static void
_ncm_catalog_set_data (NcmCatalogPrivate *self, NcmMatrix *data)
{
  ncm_matrix_clear (&self->data);

  if (data != NULL)
  {
    self->data = ncm_matrix_ref (data);
    self->len  = ncm_matrix_nrows (self->data);
  }
}

/**
 * ncm_catalog_new:
 * @nrows: the number of rows
 * @col_names: (array zero-terminated=1): the column names
 *
 * Creates a new #NcmCatalog with @nrows rows and the columns named by
 * @col_names, with all values initialized to zero.
 *
 * Returns: (transfer full): a new #NcmCatalog.
 */
NcmCatalog *
ncm_catalog_new (guint nrows, GStrv col_names)
{
  return g_object_new (NCM_TYPE_CATALOG,
                       "columns", col_names,
                       "len", nrows,
                       NULL);
}

/**
 * ncm_catalog_new_full:
 * @nrows: the number of rows
 * @col_names: (array zero-terminated=1): the column names
 * @col_types: (array length=n_types): the per-column logical types
 * @n_types: the number of entries in @col_types
 *
 * Creates a new #NcmCatalog with @nrows rows and the columns named by
 * @col_names, with all values initialized to zero. Each column carries the
 * logical type given in @col_types (which must have one entry per column).
 *
 * Returns: (transfer full): a new #NcmCatalog.
 */
NcmCatalog *
ncm_catalog_new_full (guint nrows, GStrv col_names, const NcmCatalogColType *col_types, guint n_types)
{
  GArray *types = g_array_sized_new (FALSE, TRUE, sizeof (guint32), n_types);
  GVariant *var;
  NcmCatalog *catalog;
  guint i;

  for (i = 0; i < n_types; i++)
  {
    const guint32 t = col_types[i];

    g_array_append_val (types, t);
  }

  var = ncm_cfg_array_to_variant (types, G_VARIANT_TYPE ("u"));

  catalog = g_object_new (NCM_TYPE_CATALOG,
                          "columns", col_names,
                          "col-types", var,
                          "len", nrows,
                          NULL);

  g_variant_unref (var);
  g_array_unref (types);

  return catalog;
}

/**
 * ncm_catalog_ref:
 * @catalog: a #NcmCatalog
 *
 * Increases the reference count of @catalog by one.
 *
 * Returns: (transfer full): @catalog.
 */
NcmCatalog *
ncm_catalog_ref (NcmCatalog *catalog)
{
  return g_object_ref (catalog);
}

/**
 * ncm_catalog_free:
 * @catalog: a #NcmCatalog
 *
 * Decreases the reference count of @catalog by one.
 *
 */
void
ncm_catalog_free (NcmCatalog *catalog)
{
  g_object_unref (catalog);
}

/**
 * ncm_catalog_clear:
 * @catalog: a #NcmCatalog
 *
 * If *@catalog is different from %NULL, decreases its reference count and sets
 * *@catalog to %NULL.
 *
 */
void
ncm_catalog_clear (NcmCatalog **catalog)
{
  g_clear_object (catalog);
}

/**
 * ncm_catalog_get_index:
 * @catalog: a #NcmCatalog
 * @col: a column name
 * @i: (out): a pointer to store the column index
 *
 * Gets the column index from the column name.
 *
 * Returns: %TRUE if the column was found, %FALSE otherwise.
 */
gboolean
ncm_catalog_get_index (NcmCatalog *catalog, const gchar *col, guint *i)
{
  NcmCatalogPrivate * const self = ncm_catalog_get_instance_private (catalog);
  guint *j;

  j = g_hash_table_lookup (self->columns_hash, col);

  if (j == NULL)
    return FALSE;

  *i = *j;

  return TRUE;
}

/**
 * ncm_catalog_has_column:
 * @catalog: a #NcmCatalog
 * @col: a column name
 *
 * Checks whether @catalog has a column named @col.
 *
 * Returns: %TRUE if the column exists, %FALSE otherwise.
 */
gboolean
ncm_catalog_has_column (NcmCatalog *catalog, const gchar *col)
{
  NcmCatalogPrivate * const self = ncm_catalog_get_instance_private (catalog);

  return g_hash_table_contains (self->columns_hash, col);
}

/**
 * ncm_catalog_get_col_type:
 * @catalog: a #NcmCatalog
 * @col: a column name
 * @error: a #GError for error reporting
 *
 * Gets the logical type of the column named @col.
 *
 * Returns: the #NcmCatalogColType of the column.
 */
NcmCatalogColType
ncm_catalog_get_col_type (NcmCatalog *catalog, const gchar *col, GError **error)
{
  NcmCatalogPrivate * const self = ncm_catalog_get_instance_private (catalog);
  guint *j;

  j = g_hash_table_lookup (self->columns_hash, col);

  if (j == NULL)
  {
    ncm_util_set_or_call_error (error, NCM_CATALOG_ERROR, NCM_CATALOG_ERROR_COLUMN_NOT_FOUND,
                                "Column '%s' not found.", col);

    return NCM_CATALOG_COL_TYPE_INVALID;
  }

  return g_array_index (self->col_types, guint32, *j);
}

/**
 * ncm_catalog_set:
 * @catalog: a #NcmCatalog
 * @col: a column name
 * @i: a row index
 * @val: the value to set
 * @error: a #GError for error reporting
 *
 * Sets the value at row @i and column @col.
 *
 */
void
ncm_catalog_set (NcmCatalog *catalog, const gchar *col, const guint i, gdouble val, GError **error)
{
  NcmCatalogPrivate * const self = ncm_catalog_get_instance_private (catalog);
  guint *j;

  j = g_hash_table_lookup (self->columns_hash, col);

  if (j == NULL)
  {
    ncm_util_set_or_call_error (error, NCM_CATALOG_ERROR, NCM_CATALOG_ERROR_COLUMN_NOT_FOUND,
                                "Column '%s' not found.", col);

    return;
  }

  ncm_matrix_set (self->data, i, *j, val);
}

/**
 * ncm_catalog_get:
 * @catalog: a #NcmCatalog
 * @col: a column name
 * @i: a row index
 * @error: a #GError for error reporting
 *
 * Gets the value at row @i and column @col.
 *
 * Returns: the value at row @i and column @col.
 */
gdouble
ncm_catalog_get (NcmCatalog *catalog, const gchar *col, const guint i, GError **error)
{
  NcmCatalogPrivate * const self = ncm_catalog_get_instance_private (catalog);
  guint *j;

  j = g_hash_table_lookup (self->columns_hash, col);

  if (j == NULL)
  {
    ncm_util_set_or_call_error (error, NCM_CATALOG_ERROR, NCM_CATALOG_ERROR_COLUMN_NOT_FOUND,
                                "Column '%s' not found.", col);

    return GSL_NAN;
  }

  return ncm_matrix_get (self->data, i, *j);
}

/**
 * ncm_catalog_set_int:
 * @catalog: a #NcmCatalog
 * @col: a column name
 * @i: a row index
 * @val: the value to set
 * @error: a #GError for error reporting
 *
 * Sets the integer value at row @i and column @col. The value is stored in the
 * double backing matrix (lossless for magnitudes below 2^53).
 *
 */
void
ncm_catalog_set_int (NcmCatalog *catalog, const gchar *col, const guint i, gint64 val, GError **error)
{
  ncm_catalog_set (catalog, col, i, (gdouble) val, error);
}

/**
 * ncm_catalog_get_int:
 * @catalog: a #NcmCatalog
 * @col: a column name
 * @i: a row index
 * @error: a #GError for error reporting
 *
 * Gets the integer value at row @i and column @col.
 *
 * Returns: the value at row @i and column @col as a #gint64.
 */
gint64
ncm_catalog_get_int (NcmCatalog *catalog, const gchar *col, const guint i, GError **error)
{
  return (gint64) ncm_catalog_get (catalog, col, i, error);
}

/**
 * ncm_catalog_set_bool:
 * @catalog: a #NcmCatalog
 * @col: a column name
 * @i: a row index
 * @val: the value to set
 * @error: a #GError for error reporting
 *
 * Sets the boolean value at row @i and column @col, stored as 0.0/1.0 in the
 * double backing matrix.
 *
 */
void
ncm_catalog_set_bool (NcmCatalog *catalog, const gchar *col, const guint i, gboolean val, GError **error)
{
  ncm_catalog_set (catalog, col, i, val ? 1.0 : 0.0, error);
}

/**
 * ncm_catalog_get_bool:
 * @catalog: a #NcmCatalog
 * @col: a column name
 * @i: a row index
 * @error: a #GError for error reporting
 *
 * Gets the boolean value at row @i and column @col.
 *
 * Returns: the value at row @i and column @col as a #gboolean.
 */
gboolean
ncm_catalog_get_bool (NcmCatalog *catalog, const gchar *col, const guint i, GError **error)
{
  return ncm_catalog_get (catalog, col, i, error) != 0.0;
}

/**
 * ncm_catalog_peek_columns:
 * @catalog: a #NcmCatalog
 *
 * Gets the catalog column names.
 *
 * Returns: (transfer none) (array zero-terminated=1): the column names.
 */
GStrv
ncm_catalog_peek_columns (NcmCatalog *catalog)
{
  NcmCatalogPrivate * const self = ncm_catalog_get_instance_private (catalog);

  return self->columns;
}

/**
 * ncm_catalog_peek_data:
 * @catalog: a #NcmCatalog
 *
 * Gets the catalog data matrix.
 *
 * Returns: (transfer none): the data #NcmMatrix.
 */
NcmMatrix *
ncm_catalog_peek_data (NcmCatalog *catalog)
{
  NcmCatalogPrivate * const self = ncm_catalog_get_instance_private (catalog);

  return self->data;
}

/**
 * ncm_catalog_len:
 * @catalog: a #NcmCatalog
 *
 * Gets the number of rows (entries) in the catalog.
 *
 * Returns: the number of rows.
 */
guint
ncm_catalog_len (NcmCatalog *catalog)
{
  NcmCatalogPrivate * const self = ncm_catalog_get_instance_private (catalog);

  return ncm_matrix_nrows (self->data);
}

/**
 * ncm_catalog_ncols:
 * @catalog: a #NcmCatalog
 *
 * Gets the number of columns in the catalog.
 *
 * Returns: the number of columns.
 */
guint
ncm_catalog_ncols (NcmCatalog *catalog)
{
  NcmCatalogPrivate * const self = ncm_catalog_get_instance_private (catalog);

  return self->ncols;
}

/**
 * ncm_catalog_set_meta:
 * @catalog: a #NcmCatalog
 * @meta: (nullable): a #NcmVarDict
 *
 * Sets the catalog-wide metadata to @meta, replacing any previous value.
 * Passing %NULL clears it.
 *
 */
void
ncm_catalog_set_meta (NcmCatalog *catalog, NcmVarDict *meta)
{
  NcmCatalogPrivate * const self = ncm_catalog_get_instance_private (catalog);

  g_clear_pointer (&self->meta, ncm_var_dict_unref);

  if (meta != NULL)
    self->meta = ncm_var_dict_ref (meta);
}

/**
 * ncm_catalog_peek_meta:
 * @catalog: a #NcmCatalog
 *
 * Gets the catalog-wide metadata.
 *
 * Returns: (transfer none) (nullable): the #NcmVarDict, or %NULL if unset.
 */
NcmVarDict *
ncm_catalog_peek_meta (NcmCatalog *catalog)
{
  NcmCatalogPrivate * const self = ncm_catalog_get_instance_private (catalog);

  return self->meta;
}


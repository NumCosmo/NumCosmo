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
 * A generic, lightweight container of named columns, each of a type
 * (#NcmCatalogColType). It provides column-name based access to the data and is meant
 * to be used as a derivable base class for catalog-like objects (for example galaxy
 * and halo catalogs), as well as on its own as a simple table that converts cheaply
 * to/from external table formats.
 *
 * The DOUBLE columns form the #NcmMatrix of ncm_catalog_peek_data(); every INT, BOOL and
 * DOUBLE_SEPARATE column has an array of its own, of 64-bit integers, #gboolean and
 * doubles. The accessors of a type act only on columns of that type
 * (ncm_catalog_set() and ncm_catalog_get() on DOUBLE and DOUBLE_SEPARATE, the `_int`
 * functions on INT, the `_bool` functions on BOOL); a column of another type is an error
 * (%NCM_CATALOG_ERROR_WRONG_TYPE).
 *
 * ncm_catalog_shallow_dup() makes a catalog that references the storage of another, the
 * matrix and every column array. ncm_catalog_replace_by_copy() gives one column outside the
 * matrix storage of its own, a copy of its array, and ncm_catalog_replace_data_by_copy() does
 * the same for the matrix, which holds every DOUBLE column. A catalog that writes only some
 * columns of a shared catalog thus copies only those. String and other heterogeneous columns are
 * intentionally not handled here and are left to specialized subclasses or future
 * additive extensions.
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
 * Every change of the data through the setters of the catalog or through the data
 * properties increments the key returned by ncm_catalog_get_key(), as a parameter change
 * does to the pkey of a #NcmModel. A change written directly into the matrix of
 * ncm_catalog_peek_data() must be followed by ncm_catalog_mark_outdated().
 *
 * Data derived from the rows and attached to the catalog, as by a subclass, saves the key of
 * the rows it describes and checks it with ncm_catalog_key_is_current(). The key is
 * meaningful only within one process, so such data set through a property, as by
 * deserialization, saves %NCM_CATALOG_KEY_PENDING instead: the order in which deserialization
 * sets the properties, the data matrix among them, is not fixed, and the first check takes
 * the key of the rows as they are then. Serialized derived data must be empty when it is not
 * current, so a file never carries data of other rows.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"
#include "ncm/data/ncm_catalog.h"
#include "ncm/algebra/ncm_matrix.h"
#include "ncm/core/ncm_cfg.h"

#include <string.h>

/* *INDENT-OFF* */
G_DEFINE_QUARK (ncm-catalog-error, ncm_catalog_error) 
/* *INDENT-ON* */

typedef struct _NcmCatalogPrivate
{
  NcmMatrix *data;      /* the DOUBLE columns, nrows x n_double */
  GPtrArray *col_data;  /* per column: NULL for DOUBLE, else its own array of len values */
  GVariant *col_data_v; /* #NcmCatalog:col-data set before construction */
  GStrv columns;
  GHashTable *columns_hash;
  GArray *col_types;
  GArray *col_pos; /* column of each DOUBLE column in the matrix */
  guint ncols;
  guint n_double;
  guint len;
  NcmVarDict *meta;
  guint64 key;
  gboolean constructed;
} NcmCatalogPrivate;

enum
{
  PROP_0,
  PROP_DATA,
  PROP_COL_DATA,
  PROP_COLUMNS,
  PROP_COL_TYPES,
  PROP_LEN,
  PROP_META,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcmCatalog, ncm_catalog, G_TYPE_OBJECT)

/* col_data holds NULL for the DOUBLE columns */
static void
_ncm_catalog_col_free (gpointer a)
{
  if (a != NULL)
    g_array_unref (a);
}

static void
ncm_catalog_init (NcmCatalog *catalog)
{
  NcmCatalogPrivate * const self = ncm_catalog_get_instance_private (catalog);

  self->data         = NULL;
  self->col_data     = g_ptr_array_new_with_free_func (&_ncm_catalog_col_free);
  self->col_data_v   = NULL;
  self->columns      = NULL;
  self->columns_hash = g_hash_table_new_full (g_str_hash, g_str_equal, g_free, g_free);
  self->col_types    = NULL;
  self->col_pos      = g_array_new (FALSE, TRUE, sizeof (guint));
  self->ncols        = 0;
  self->n_double     = 0;
  self->len          = 0;
  self->key          = 1;
  self->meta         = NULL;
  self->constructed  = FALSE;
}

static void _ncm_catalog_set_data (NcmCatalog *catalog, NcmMatrix *data);
static void _ncm_catalog_set_col_data (NcmCatalog *catalog, GVariant *var);
static GVariant *_ncm_catalog_col_data_to_variant (NcmCatalog *catalog);
static NcmCatalog *_ncm_catalog_shallow_dup (NcmCatalog *catalog);

static void
_ncm_catalog_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmCatalog *catalog            = NCM_CATALOG (object);
  NcmCatalogPrivate * const self = ncm_catalog_get_instance_private (catalog);

  switch (prop_id)
  {
    case PROP_DATA:
      _ncm_catalog_set_data (catalog, g_value_get_object (value));
      break;
    case PROP_COL_DATA:
      _ncm_catalog_set_col_data (catalog, g_value_get_variant (value));
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
    case PROP_COL_DATA:
      g_value_take_variant (value, _ncm_catalog_col_data_to_variant (catalog));
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

    if (self->columns == NULL)
      g_error ("_ncm_catalog_constructed: a %s needs its column names.", G_OBJECT_TYPE_NAME (object));

    self->ncols = g_strv_length (self->columns);

    for (i = 0; i < self->ncols; i++)
    {
      const gchar *col = self->columns[i];
      guint *index     = g_new0 (guint, 1);

      *index = i;

      if (g_hash_table_contains (self->columns_hash, col))
        g_error ("_ncm_catalog_constructed: column `%s' appears more than once.", col);

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
      if (self->col_types->len != self->ncols)
        g_error ("_ncm_catalog_constructed: %u column types for %u columns.", self->col_types->len, self->ncols);

      for (i = 0; i < self->ncols; i++)
      {
        const guint32 t = g_array_index (self->col_types, guint32, i);

        if (t >= NCM_CATALOG_COL_TYPE_LEN)
          g_error ("_ncm_catalog_constructed: column `%s' has the invalid type %u.", self->columns[i], t);
      }
    }

    g_array_set_size (self->col_pos, self->ncols);
    g_ptr_array_set_size (self->col_data, self->ncols);

    for (i = 0; i < self->ncols; i++)
    {
      const guint32 type = g_array_index (self->col_types, guint32, i);

      if (type == NCM_CATALOG_COL_TYPE_DOUBLE)
      {
        g_array_index (self->col_pos, guint, i) = self->n_double++;
      }
      else
      {
        const guint size = (type == NCM_CATALOG_COL_TYPE_INT) ? sizeof (gint64) :
                           (type == NCM_CATALOG_COL_TYPE_BOOL) ? sizeof (gboolean) : sizeof (gdouble);
        GArray *a = g_array_sized_new (FALSE, TRUE, size, self->len);

        g_array_set_size (a, self->len);
        g_ptr_array_index (self->col_data, i) = a;
      }
    }

    if (self->data == NULL)
    {
      self->data = ncm_matrix_new (self->len, self->n_double);
      ncm_matrix_set_zero (self->data);
    }
    else if (ncm_matrix_ncols (self->data) != self->n_double)
    {
      g_error ("_ncm_catalog_constructed: the catalog has %u double columns, but the matrix has %u.",
               self->n_double, ncm_matrix_ncols (self->data));
    }

    self->constructed = TRUE;

    if (self->col_data_v != NULL)
    {
      GVariant *var = self->col_data_v;

      self->col_data_v = NULL;
      _ncm_catalog_set_col_data (catalog, var);
      g_variant_unref (var);
    }
  }
}

static void
_ncm_catalog_dispose (GObject *object)
{
  NcmCatalog *catalog            = NCM_CATALOG (object);
  NcmCatalogPrivate * const self = ncm_catalog_get_instance_private (catalog);

  ncm_matrix_clear (&self->data);
  g_clear_pointer (&self->col_data, g_ptr_array_unref);
  g_clear_pointer (&self->col_data_v, g_variant_unref);
  g_clear_pointer (&self->columns_hash, g_hash_table_unref);
  g_clear_pointer (&self->columns, g_strfreev);
  g_clear_pointer (&self->col_types, g_array_unref);
  g_clear_pointer (&self->col_pos, g_array_unref);
  g_clear_pointer (&self->meta, ncm_var_dict_unref);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_catalog_parent_class)->dispose (object);
}

static void
ncm_catalog_class_init (NcmCatalogClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = &_ncm_catalog_set_property;
  object_class->get_property = &_ncm_catalog_get_property;
  object_class->constructed  = &_ncm_catalog_constructed;
  object_class->dispose      = &_ncm_catalog_dispose;
  klass->shallow_dup         = &_ncm_catalog_shallow_dup;

  /**
   * NcmCatalog:data:
   *
   * The matrix of the %NCM_CATALOG_COL_TYPE_DOUBLE columns, in catalog order (rows are
   * entries). Setting it replaces those columns; it must have one column per DOUBLE
   * column. In a catalog with columns outside the matrix it must keep the number of rows;
   * in one without, it replaces the number of rows.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_DATA,
                                   g_param_spec_object ("data",
                                                        "Data",
                                                        "Catalog double columns",
                                                        NCM_TYPE_MATRIX,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_STRINGS));

  /**
   * NcmCatalog:col-data:
   *
   * The columns that are not DOUBLE, in catalog order, one variant each: an array of
   * 64-bit integers for an INT column, of booleans for a BOOL column and of doubles for a
   * DOUBLE_SEPARATE column, with one value per row.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_COL_DATA,
                                   g_param_spec_variant ("col-data",
                                                         "Column data",
                                                         "Catalog columns outside the double matrix",
                                                         G_VARIANT_TYPE ("av"), NULL,
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
_ncm_catalog_set_data (NcmCatalog *catalog, NcmMatrix *data)
{
  NcmCatalogPrivate * const self = ncm_catalog_get_instance_private (catalog);

  if (data == NULL)
  {
    /* Unset during construction: constructed allocates the zero matrix. */
    if (self->data != NULL)
      g_error ("_ncm_catalog_set_data: the data matrix of a catalog cannot be unset.");

    return;
  }

  /* Before construction the column types are not known yet, so constructed checks the size */
  if (self->constructed)
  {
    if (ncm_matrix_ncols (data) != self->n_double)
      g_error ("_ncm_catalog_set_data: the catalog has %u double columns, but the matrix has %u.",
               self->n_double, ncm_matrix_ncols (data));

    if ((self->n_double < self->ncols) && (ncm_matrix_nrows (data) != self->len))
      g_error ("_ncm_catalog_set_data: the catalog has %u rows and columns outside the matrix, "
               "but the matrix has %u rows.", self->len, ncm_matrix_nrows (data));
  }

  ncm_matrix_ref (data);
  ncm_matrix_clear (&self->data);
  self->data = data;
  self->len  = ncm_matrix_nrows (self->data);
  self->key++;
}

/*
 * Sets the columns outside the matrix from @var, an array with one variant per such column
 * in catalog order. Before construction the column types are not known yet, so @var is
 * kept and constructed sets it.
 */
static void
_ncm_catalog_set_col_data (NcmCatalog *catalog, GVariant *var)
{
  NcmCatalogPrivate * const self = ncm_catalog_get_instance_private (catalog);
  guint i, k = 0;

  if (var == NULL)
    return;

  if (!self->constructed)
  {
    g_clear_pointer (&self->col_data_v, g_variant_unref);
    self->col_data_v = g_variant_ref (var);

    return;
  }

  for (i = 0; i < self->ncols; i++)
  {
    const guint32 type = g_array_index (self->col_types, guint32, i);
    GArray *a          = g_ptr_array_index (self->col_data, i);
    GVariant *cvar;

    if (a == NULL)
      continue;

    if (k >= g_variant_n_children (var))
      g_error ("_ncm_catalog_set_col_data: %zu arrays for the columns outside the matrix, but column `%s' is one.",
               g_variant_n_children (var), self->columns[i]);

    cvar = g_variant_get_child_value (var, k++);
    {
      GVariant *col = g_variant_get_variant (cvar);

      if (type == NCM_CATALOG_COL_TYPE_BOOL)
      {
        gsize n             = 0;
        const guchar *bytes = g_variant_get_fixed_array (col, &n, sizeof (guchar));
        gsize r;

        /* GVariant stores a boolean in one byte; in memory each value is a gboolean */
        g_array_set_size (a, n);

        for (r = 0; r < n; r++)
          g_array_index (a, gboolean, r) = (bytes[r] != 0);
      }
      else
      {
        ncm_cfg_array_set_variant (a, col);
      }

      if (a->len != self->len)
        g_error ("_ncm_catalog_set_col_data: column `%s' has %u values for %u rows.", self->columns[i], a->len, self->len);

      g_variant_unref (col);
    }
    g_variant_unref (cvar);
  }

  if (k != g_variant_n_children (var))
    g_error ("_ncm_catalog_set_col_data: %zu arrays for %u columns outside the matrix.", g_variant_n_children (var), k);

  self->key++;
}

static GVariant *
_ncm_catalog_col_data_to_variant (NcmCatalog *catalog)
{
  NcmCatalogPrivate * const self = ncm_catalog_get_instance_private (catalog);
  GVariantBuilder builder;
  guint i;

  g_variant_builder_init (&builder, G_VARIANT_TYPE ("av"));

  for (i = 0; i < self->ncols; i++)
  {
    const guint32 type = g_array_index (self->col_types, guint32, i);
    GArray *a          = g_ptr_array_index (self->col_data, i);
    GVariant *col;

    if (a == NULL)
      continue;

    if (type == NCM_CATALOG_COL_TYPE_BOOL)
    {
      guchar *bytes = g_new (guchar, MAX (a->len, 1));
      guint r;

      for (r = 0; r < a->len; r++)
        bytes[r] = g_array_index (a, gboolean, r) ? 1 : 0;

      col = g_variant_ref_sink (g_variant_new_fixed_array (G_VARIANT_TYPE_BOOLEAN, bytes, a->len, sizeof (guchar)));
      g_free (bytes);
    }
    else
    {
      col = ncm_cfg_array_to_variant (a, (type == NCM_CATALOG_COL_TYPE_INT) ? G_VARIANT_TYPE_INT64 : G_VARIANT_TYPE_DOUBLE);
    }

    g_variant_builder_add (&builder, "v", col);
    g_variant_unref (col);
  }

  return g_variant_ref_sink (g_variant_builder_end (&builder));
}

/* A catalog of the same type, columns and rows that references the matrix and the column
 * arrays of @catalog */
static NcmCatalog *
_ncm_catalog_shallow_dup (NcmCatalog *catalog)
{
  NcmCatalogPrivate * const self = ncm_catalog_get_instance_private (catalog);
  GVariant *types                = ncm_cfg_array_to_variant (self->col_types, G_VARIANT_TYPE ("u"));
  NcmCatalog *dup                = g_object_new (G_OBJECT_TYPE (catalog),
                                                 "columns", self->columns,
                                                 "col-types", types,
                                                 "len", self->len,
                                                 "data", self->data,
                                                 "meta", self->meta,
                                                 NULL);
  NcmCatalogPrivate * const dself = ncm_catalog_get_instance_private (dup);
  guint i;

  for (i = 0; i < self->ncols; i++)
  {
    GArray *a = g_ptr_array_index (self->col_data, i);

    if (a != NULL)
    {
      g_array_unref (g_ptr_array_index (dself->col_data, i));
      g_ptr_array_index (dself->col_data, i) = g_array_ref (a);
    }
  }

  g_variant_unref (types);

  return dup;
}

/**
 * ncm_catalog_shallow_dup: (virtual shallow_dup)
 * @catalog: a #NcmCatalog
 *
 * Creates a catalog of the same type, columns and rows that references the storage of
 * @catalog: the matrix of the DOUBLE columns and the array of every other column. A
 * subclass also gives the copy the rest of its state. A value written into a column of
 * either catalog is seen by both until ncm_catalog_replace_by_copy() or, for a DOUBLE
 * column, ncm_catalog_replace_data_by_copy() gives that column storage of its own.
 *
 * Returns: (transfer full): the new #NcmCatalog.
 */
NcmCatalog *
ncm_catalog_shallow_dup (NcmCatalog *catalog)
{
  return NCM_CATALOG_GET_CLASS (catalog)->shallow_dup (catalog);
}

/**
 * ncm_catalog_replace_by_copy:
 * @catalog: a #NcmCatalog
 * @col: a column name
 * @error: a #GError for error reporting
 *
 * Replaces the array of the column @col, a column outside the matrix, by a copy with the same
 * values, so a later write to @catalog no longer reaches a catalog that shares the array
 * (ncm_catalog_shallow_dup()). An unknown column sets @error, and so does a DOUBLE column,
 * whose storage is the matrix: ncm_catalog_replace_data_by_copy() copies it.
 */
void
ncm_catalog_replace_by_copy (NcmCatalog *catalog, const gchar *col, GError **error)
{
  NcmCatalogPrivate * const self = ncm_catalog_get_instance_private (catalog);
  guint *jp                      = g_hash_table_lookup (self->columns_hash, col);
  GArray *a, *copy;

  if (jp == NULL)
  {
    ncm_util_set_or_call_error (error, NCM_CATALOG_ERROR, NCM_CATALOG_ERROR_COLUMN_NOT_FOUND,
                                "Column '%s' not found.", col);

    return;
  }

  a = g_ptr_array_index (self->col_data, *jp);

  if (a == NULL)
  {
    ncm_util_set_or_call_error (error, NCM_CATALOG_ERROR, NCM_CATALOG_ERROR_WRONG_TYPE,
                                "Column '%s' is a DOUBLE column, stored in the matrix: "
                                "use ncm_catalog_replace_data_by_copy().", col);

    return;
  }

  copy = g_array_sized_new (FALSE, FALSE, g_array_get_element_size (a), a->len);
  g_array_append_vals (copy, a->data, a->len);
  g_array_unref (a);
  g_ptr_array_index (self->col_data, *jp) = copy;
}

/**
 * ncm_catalog_replace_data_by_copy:
 * @catalog: a #NcmCatalog
 *
 * Replaces the matrix of the DOUBLE columns (ncm_catalog_peek_data()) by a copy with the same
 * values, so a later write to a DOUBLE column of @catalog no longer reaches a catalog that
 * shares the matrix (ncm_catalog_shallow_dup()). Each call copies the matrix.
 */
void
ncm_catalog_replace_data_by_copy (NcmCatalog *catalog)
{
  NcmCatalogPrivate * const self = ncm_catalog_get_instance_private (catalog);
  NcmMatrix *data                = ncm_matrix_dup (self->data);

  ncm_matrix_free (self->data);
  self->data = data;
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
 * Decreases the reference count of @catalog by one. If the reference count
 * reaches zero, @catalog is freed.
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

static const gchar *
_ncm_catalog_col_type_name (const NcmCatalogColType type)
{
  switch (type)
  {
    case NCM_CATALOG_COL_TYPE_DOUBLE:
      return "DOUBLE";

    case NCM_CATALOG_COL_TYPE_INT:
      return "INT";

    case NCM_CATALOG_COL_TYPE_BOOL:
      return "BOOL";

    case NCM_CATALOG_COL_TYPE_DOUBLE_SEPARATE:
      return "DOUBLE_SEPARATE";

    default:                   /* LCOV_EXCL_LINE */
      g_assert_not_reached (); /* LCOV_EXCL_LINE */

      return NULL; /* LCOV_EXCL_LINE */
  }
}

/*
 * Finds the column @col and checks that its type is @type, DOUBLE standing for both DOUBLE
 * and DOUBLE_SEPARATE; sets @pos to the index of the column. On failure sets @error (or
 * aborts when @error is %NULL) and returns %FALSE.
 */
static gboolean
_ncm_catalog_locate_column (NcmCatalog *catalog, const gchar *col, const NcmCatalogColType type, guint *pos, GError **error)
{
  NcmCatalogPrivate * const self = ncm_catalog_get_instance_private (catalog);
  guint *jp                      = g_hash_table_lookup (self->columns_hash, col);
  NcmCatalogColType col_type;

  if (jp == NULL)
  {
    ncm_util_set_or_call_error (error, NCM_CATALOG_ERROR, NCM_CATALOG_ERROR_COLUMN_NOT_FOUND,
                                "Column '%s' not found.", col);

    return FALSE;
  }

  col_type = g_array_index (self->col_types, guint32, *jp);

  if ((col_type != type) &&
      !((type == NCM_CATALOG_COL_TYPE_DOUBLE) && (col_type == NCM_CATALOG_COL_TYPE_DOUBLE_SEPARATE)))
  {
    ncm_util_set_or_call_error (error, NCM_CATALOG_ERROR, NCM_CATALOG_ERROR_WRONG_TYPE,
                                "Column '%s' is %s, not %s.", col,
                                _ncm_catalog_col_type_name (col_type), _ncm_catalog_col_type_name (type));

    return FALSE;
  }

  *pos = *jp;

  return TRUE;
}

/*
 * As _ncm_catalog_locate_column(), and checks the row @i.
 */
static gboolean
_ncm_catalog_locate (NcmCatalog *catalog, const gchar *col, const guint i, const NcmCatalogColType type, guint *pos, GError **error)
{
  NcmCatalogPrivate * const self = ncm_catalog_get_instance_private (catalog);

  if (!_ncm_catalog_locate_column (catalog, col, type, pos, error))
    return FALSE;

  if (i >= self->len)
  {
    ncm_util_set_or_call_error (error, NCM_CATALOG_ERROR, NCM_CATALOG_ERROR_ROW_OUT_OF_RANGE,
                                "Row %u of column '%s' requested, but the catalog has %u rows.",
                                i, col, self->len);

    return FALSE;
  }

  return TRUE;
}

/**
 * ncm_catalog_get_data_index:
 * @catalog: a #NcmCatalog
 * @col: a column name
 * @j: (out): a pointer to store the matrix column
 *
 * Gets the column of ncm_catalog_peek_data() that holds the %NCM_CATALOG_COL_TYPE_DOUBLE
 * column @col.
 *
 * Returns: %TRUE if @col is a DOUBLE column, %FALSE if it is missing or of another type.
 */
gboolean
ncm_catalog_get_data_index (NcmCatalog *catalog, const gchar *col, guint *j)
{
  NcmCatalogPrivate * const self = ncm_catalog_get_instance_private (catalog);
  guint *jp                      = g_hash_table_lookup (self->columns_hash, col);

  if ((jp == NULL) || (g_array_index (self->col_types, guint32, *jp) != NCM_CATALOG_COL_TYPE_DOUBLE))
    return FALSE;

  *j = g_array_index (self->col_pos, guint, *jp);

  return TRUE;
}

/**
 * ncm_catalog_set:
 * @catalog: a #NcmCatalog
 * @col: a %NCM_CATALOG_COL_TYPE_DOUBLE or %NCM_CATALOG_COL_TYPE_DOUBLE_SEPARATE column name
 * @i: a row index
 * @val: the value to set
 * @error: a #GError for error reporting
 *
 * Sets the value at row @i of the real-valued column @col. An unknown column, a column of
 * another type or a row index not smaller than ncm_catalog_len() sets @error.
 *
 */
void
ncm_catalog_set (NcmCatalog *catalog, const gchar *col, const guint i, gdouble val, GError **error)
{
  NcmCatalogPrivate * const self = ncm_catalog_get_instance_private (catalog);
  guint j;

  if (_ncm_catalog_locate (catalog, col, i, NCM_CATALOG_COL_TYPE_DOUBLE, &j, error))
  {
    GArray *a = g_ptr_array_index (self->col_data, j);

    if (a == NULL)
      ncm_matrix_set (self->data, i, g_array_index (self->col_pos, guint, j), val);
    else
      g_array_index (a, gdouble, i) = val;

    self->key++;
  }
}

/**
 * ncm_catalog_get:
 * @catalog: a #NcmCatalog
 * @col: a %NCM_CATALOG_COL_TYPE_DOUBLE or %NCM_CATALOG_COL_TYPE_DOUBLE_SEPARATE column name
 * @i: a row index
 * @error: a #GError for error reporting
 *
 * Gets the value at row @i of the real-valued column @col. An unknown column, a column of
 * another type or a row index not smaller than ncm_catalog_len() sets @error and
 * returns NaN.
 *
 * Returns: the value at row @i and column @col.
 */
gdouble
ncm_catalog_get (NcmCatalog *catalog, const gchar *col, const guint i, GError **error)
{
  NcmCatalogPrivate * const self = ncm_catalog_get_instance_private (catalog);
  guint j;

  if (!_ncm_catalog_locate (catalog, col, i, NCM_CATALOG_COL_TYPE_DOUBLE, &j, error))
    return GSL_NAN;

  if (g_ptr_array_index (self->col_data, j) == NULL)
    return ncm_matrix_get (self->data, i, g_array_index (self->col_pos, guint, j));

  return g_array_index ((GArray *) g_ptr_array_index (self->col_data, j), gdouble, i);
}

/**
 * ncm_catalog_set_int:
 * @catalog: a #NcmCatalog
 * @col: a %NCM_CATALOG_COL_TYPE_INT column name
 * @i: a row index
 * @val: the value to set
 * @error: a #GError for error reporting
 *
 * Sets the value at row @i of the INT column @col, stored as a 64-bit integer. An
 * unknown column, a column of another type or a row index not smaller than
 * ncm_catalog_len() sets @error.
 *
 */
void
ncm_catalog_set_int (NcmCatalog *catalog, const gchar *col, const guint i, gint64 val, GError **error)
{
  NcmCatalogPrivate * const self = ncm_catalog_get_instance_private (catalog);
  guint j;

  if (_ncm_catalog_locate (catalog, col, i, NCM_CATALOG_COL_TYPE_INT, &j, error))
  {
    g_array_index ((GArray *) g_ptr_array_index (self->col_data, j), gint64, i) = val;
    self->key++;
  }
}

/**
 * ncm_catalog_get_int:
 * @catalog: a #NcmCatalog
 * @col: a %NCM_CATALOG_COL_TYPE_INT column name
 * @i: a row index
 * @error: a #GError for error reporting
 *
 * Gets the value at row @i of the INT column @col. An unknown column, a column of
 * another type or a row index not smaller than ncm_catalog_len() sets @error and
 * returns 0.
 *
 * Returns: the value at row @i and column @col.
 */
gint64
ncm_catalog_get_int (NcmCatalog *catalog, const gchar *col, const guint i, GError **error)
{
  NcmCatalogPrivate * const self = ncm_catalog_get_instance_private (catalog);
  guint j;

  if (!_ncm_catalog_locate (catalog, col, i, NCM_CATALOG_COL_TYPE_INT, &j, error))
    return 0;

  return g_array_index ((GArray *) g_ptr_array_index (self->col_data, j), gint64, i);
}

/**
 * ncm_catalog_set_bool:
 * @catalog: a #NcmCatalog
 * @col: a %NCM_CATALOG_COL_TYPE_BOOL column name
 * @i: a row index
 * @val: the value to set
 * @error: a #GError for error reporting
 *
 * Sets the value at row @i of the BOOL column @col. An unknown column, a column of
 * another type or a row index not smaller than ncm_catalog_len() sets @error.
 *
 */
void
ncm_catalog_set_bool (NcmCatalog *catalog, const gchar *col, const guint i, gboolean val, GError **error)
{
  NcmCatalogPrivate * const self = ncm_catalog_get_instance_private (catalog);
  guint j;

  if (_ncm_catalog_locate (catalog, col, i, NCM_CATALOG_COL_TYPE_BOOL, &j, error))
  {
    g_array_index ((GArray *) g_ptr_array_index (self->col_data, j), gboolean, i) = (val != FALSE);
    self->key++;
  }
}

/**
 * ncm_catalog_get_bool:
 * @catalog: a #NcmCatalog
 * @col: a %NCM_CATALOG_COL_TYPE_BOOL column name
 * @i: a row index
 * @error: a #GError for error reporting
 *
 * Gets the value at row @i of the BOOL column @col. An unknown column, a column of
 * another type or a row index not smaller than ncm_catalog_len() sets @error and
 * returns %FALSE.
 *
 * Returns: the value at row @i and column @col.
 */
gboolean
ncm_catalog_get_bool (NcmCatalog *catalog, const gchar *col, const guint i, GError **error)
{
  NcmCatalogPrivate * const self = ncm_catalog_get_instance_private (catalog);
  guint j;

  if (!_ncm_catalog_locate (catalog, col, i, NCM_CATALOG_COL_TYPE_BOOL, &j, error))
    return FALSE;

  return g_array_index ((GArray *) g_ptr_array_index (self->col_data, j), gboolean, i);
}

/**
 * ncm_catalog_set_int_column:
 * @catalog: a #NcmCatalog
 * @col: a %NCM_CATALOG_COL_TYPE_INT column name
 * @vals: (element-type gint64): the values, one per row
 * @error: a #GError for error reporting
 *
 * Sets every row of the INT column @col from @vals. An unknown column, a column of
 * another type or a length of @vals other than ncm_catalog_len() sets @error.
 *
 */
void
ncm_catalog_set_int_column (NcmCatalog *catalog, const gchar *col, GArray *vals, GError **error)
{
  NcmCatalogPrivate * const self = ncm_catalog_get_instance_private (catalog);
  guint j;

  if (!_ncm_catalog_locate_column (catalog, col, NCM_CATALOG_COL_TYPE_INT, &j, error))
    return;

  if (vals->len != self->len)
  {
    ncm_util_set_or_call_error (error, NCM_CATALOG_ERROR, NCM_CATALOG_ERROR_ROW_OUT_OF_RANGE,
                                "ncm_catalog_set_int_column: %u values for the %u rows of column '%s'.",
                                vals->len, self->len, col);

    return;
  }

  memcpy (((GArray *) g_ptr_array_index (self->col_data, j))->data, vals->data, sizeof (gint64) * self->len);

  self->key++;
}

/**
 * ncm_catalog_get_int_column:
 * @catalog: a #NcmCatalog
 * @col: a %NCM_CATALOG_COL_TYPE_INT column name
 * @vals: (out callee-allocates) (transfer full) (element-type gint64): the values
 * @error: a #GError for error reporting
 *
 * Gets every row of the INT column @col. A %NULL *@vals is created; otherwise it is
 * resized. An unknown column or a column of another type sets @error.
 *
 */
void
ncm_catalog_get_int_column (NcmCatalog *catalog, const gchar *col, GArray **vals, GError **error)
{
  NcmCatalogPrivate * const self = ncm_catalog_get_instance_private (catalog);
  guint j;

  if (!_ncm_catalog_locate_column (catalog, col, NCM_CATALOG_COL_TYPE_INT, &j, error))
    return;

  if (*vals == NULL)
    *vals = g_array_sized_new (FALSE, FALSE, sizeof (gint64), self->len);

  g_array_set_size (*vals, self->len);

  memcpy ((*vals)->data, ((GArray *) g_ptr_array_index (self->col_data, j))->data, sizeof (gint64) * self->len);
}

/**
 * ncm_catalog_set_bool_column:
 * @catalog: a #NcmCatalog
 * @col: a %NCM_CATALOG_COL_TYPE_BOOL column name
 * @vals: (element-type gboolean): the values, one per row
 * @error: a #GError for error reporting
 *
 * Sets every row of the BOOL column @col from @vals. An unknown column, a column of
 * another type or a length of @vals other than ncm_catalog_len() sets @error.
 *
 */
void
ncm_catalog_set_bool_column (NcmCatalog *catalog, const gchar *col, GArray *vals, GError **error)
{
  NcmCatalogPrivate * const self = ncm_catalog_get_instance_private (catalog);
  guint j, i;

  if (!_ncm_catalog_locate_column (catalog, col, NCM_CATALOG_COL_TYPE_BOOL, &j, error))
    return;

  if (vals->len != self->len)
  {
    ncm_util_set_or_call_error (error, NCM_CATALOG_ERROR, NCM_CATALOG_ERROR_ROW_OUT_OF_RANGE,
                                "ncm_catalog_set_bool_column: %u values for the %u rows of column '%s'.",
                                vals->len, self->len, col);

    return;
  }

  for (i = 0; i < self->len; i++)
    g_array_index ((GArray *) g_ptr_array_index (self->col_data, j), gboolean, i) = (g_array_index (vals, gboolean, i) != FALSE);

  self->key++;
}

/**
 * ncm_catalog_get_bool_column:
 * @catalog: a #NcmCatalog
 * @col: a %NCM_CATALOG_COL_TYPE_BOOL column name
 * @vals: (out callee-allocates) (transfer full) (element-type gboolean): the values
 * @error: a #GError for error reporting
 *
 * Gets every row of the BOOL column @col. A %NULL *@vals is created; otherwise it is
 * resized. An unknown column or a column of another type sets @error.
 *
 */
void
ncm_catalog_get_bool_column (NcmCatalog *catalog, const gchar *col, GArray **vals, GError **error)
{
  NcmCatalogPrivate * const self = ncm_catalog_get_instance_private (catalog);
  guint j;

  if (!_ncm_catalog_locate_column (catalog, col, NCM_CATALOG_COL_TYPE_BOOL, &j, error))
    return;

  if (*vals == NULL)
    *vals = g_array_sized_new (FALSE, FALSE, sizeof (gboolean), self->len);

  g_array_set_size (*vals, self->len);

  memcpy ((*vals)->data, ((GArray *) g_ptr_array_index (self->col_data, j))->data, sizeof (gboolean) * self->len);
}

/**
 * ncm_catalog_set_double_column:
 * @catalog: a #NcmCatalog
 * @col: a %NCM_CATALOG_COL_TYPE_DOUBLE or %NCM_CATALOG_COL_TYPE_DOUBLE_SEPARATE column name
 * @vals: (element-type gdouble): the values, one per row
 * @error: a #GError for error reporting
 *
 * Sets every row of the real-valued column @col from @vals. An unknown column, a column of
 * another type or a length of @vals other than ncm_catalog_len() sets @error.
 *
 */
void
ncm_catalog_set_double_column (NcmCatalog *catalog, const gchar *col, GArray *vals, GError **error)
{
  NcmCatalogPrivate * const self = ncm_catalog_get_instance_private (catalog);
  guint j, i;

  if (!_ncm_catalog_locate_column (catalog, col, NCM_CATALOG_COL_TYPE_DOUBLE, &j, error))
    return;

  if (vals->len != self->len)
  {
    ncm_util_set_or_call_error (error, NCM_CATALOG_ERROR, NCM_CATALOG_ERROR_ROW_OUT_OF_RANGE,
                                "ncm_catalog_set_double_column: %u values for the %u rows of column '%s'.",
                                vals->len, self->len, col);

    return;
  }

  if (g_ptr_array_index (self->col_data, j) == NULL)
    for (i = 0; i < self->len; i++)
      ncm_matrix_set (self->data, i, g_array_index (self->col_pos, guint, j), g_array_index (vals, gdouble, i));

  else
    memcpy (((GArray *) g_ptr_array_index (self->col_data, j))->data, vals->data, sizeof (gdouble) * self->len);

  self->key++;
}

/**
 * ncm_catalog_get_double_column:
 * @catalog: a #NcmCatalog
 * @col: a %NCM_CATALOG_COL_TYPE_DOUBLE or %NCM_CATALOG_COL_TYPE_DOUBLE_SEPARATE column name
 * @vals: (out callee-allocates) (transfer full) (element-type gdouble): the values
 * @error: a #GError for error reporting
 *
 * Gets every row of the real-valued column @col. A %NULL *@vals is created; otherwise it
 * is resized. An unknown column or a column of another type sets @error.
 *
 */
void
ncm_catalog_get_double_column (NcmCatalog *catalog, const gchar *col, GArray **vals, GError **error)
{
  NcmCatalogPrivate * const self = ncm_catalog_get_instance_private (catalog);
  guint j, i;

  if (!_ncm_catalog_locate_column (catalog, col, NCM_CATALOG_COL_TYPE_DOUBLE, &j, error))
    return;

  if (*vals == NULL)
    *vals = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), self->len);

  g_array_set_size (*vals, self->len);

  if (g_ptr_array_index (self->col_data, j) == NULL)
    for (i = 0; i < self->len; i++)
      g_array_index (*vals, gdouble, i) = ncm_matrix_get (self->data, i, g_array_index (self->col_pos, guint, j));

  else
    memcpy ((*vals)->data, ((GArray *) g_ptr_array_index (self->col_data, j))->data, sizeof (gdouble) * self->len);
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
 * Gets the matrix of the %NCM_CATALOG_COL_TYPE_DOUBLE columns, in catalog order;
 * ncm_catalog_get_data_index() gives the matrix column of a DOUBLE column. A change
 * written into it must be followed by ncm_catalog_mark_outdated().
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
 * ncm_catalog_get_key:
 * @catalog: a #NcmCatalog
 *
 * Gets the key of the data of @catalog, which every change of the data increments (see
 * ncm_catalog_mark_outdated()). The key is meaningful only within one process: a
 * deserialized catalog has keys of its own.
 *
 * Returns: the current key of the data.
 */
guint64
ncm_catalog_get_key (NcmCatalog *catalog)
{
  NcmCatalogPrivate * const self = ncm_catalog_get_instance_private (catalog);

  return self->key;
}

/**
 * ncm_catalog_mark_outdated:
 * @catalog: a #NcmCatalog
 *
 * Increments the key of the data of @catalog. Call it after writing directly into the
 * matrix of ncm_catalog_peek_data(); the setters of the catalog call it themselves.
 */
void
ncm_catalog_mark_outdated (NcmCatalog *catalog)
{
  NcmCatalogPrivate * const self = ncm_catalog_get_instance_private (catalog);

  self->key++;
}

/**
 * ncm_catalog_key_is_current:
 * @catalog: a #NcmCatalog
 * @saved_key: (inout): the key saved with data derived from the rows
 *
 * Checks whether data saved with @saved_key describes the current rows of @catalog. A
 * @saved_key equal to %NCM_CATALOG_KEY_PENDING takes the current key first.
 *
 * Returns: whether @saved_key equals ncm_catalog_get_key()
 */
gboolean
ncm_catalog_key_is_current (NcmCatalog *catalog, guint64 *saved_key)
{
  NcmCatalogPrivate * const self = ncm_catalog_get_instance_private (catalog);

  if (*saved_key == NCM_CATALOG_KEY_PENDING)
    *saved_key = self->key;

  return *saved_key == self->key;
}

/**
 * ncm_catalog_get_data_hash:
 * @catalog: a #NcmCatalog
 *
 * Computes the SHA-256 of the column names, the column types and the values of @catalog,
 * in this order: each name followed by a zero byte, the types as 32-bit integers, the
 * DOUBLE columns of each row as doubles, row by row, then every other column in catalog
 * order (INT as 64-bit integers, BOOL as #gboolean, DOUBLE_SEPARATE as doubles), in the
 * byte order of the machine. Two catalogs
 * with the same columns and the same values give the same hash on machines of the same
 * byte order; #NcmCatalog:meta does not enter it.
 *
 * Returns: (transfer full): the hash as a string of 64 hexadecimal digits.
 */
gchar *
ncm_catalog_get_data_hash (NcmCatalog *catalog)
{
  NcmCatalogPrivate * const self = ncm_catalog_get_instance_private (catalog);
  GChecksum *sum                 = g_checksum_new (G_CHECKSUM_SHA256);
  gchar *hash;
  guint i;

  for (i = 0; i < self->ncols; i++)
    g_checksum_update (sum, (const guchar *) self->columns[i], strlen (self->columns[i]) + 1);

  g_checksum_update (sum, (const guchar *) self->col_types->data, sizeof (guint32) * self->ncols);

  if (self->n_double > 0)
    for (i = 0; i < self->len; i++)
      g_checksum_update (sum, (const guchar *) ncm_matrix_ptr (self->data, i, 0), sizeof (gdouble) * self->n_double);


  for (i = 0; i < self->ncols; i++)
  {
    GArray *a = g_ptr_array_index (self->col_data, i);

    if (a != NULL)
      g_checksum_update (sum, (const guchar *) a->data, g_array_get_element_size (a) * a->len);
  }

  hash = g_strdup (g_checksum_get_string (sum));
  g_checksum_free (sum);

  return hash;
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

  return self->len;
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


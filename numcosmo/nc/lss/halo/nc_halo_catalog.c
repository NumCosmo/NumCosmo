/***************************************************************************
 *            nc_halo_catalog.c
 *
 *  Sat Jun 13 16:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_halo_catalog.c
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
 * NcHaloCatalog:
 *
 * Catalog of halos, clusters or member objects.
 *
 * A thin specialization of #NcmCatalog that tags its rows with a logical
 * #NcHaloCatalogKind (halo, cluster or member) and knows which of its columns
 * hold the row id and the parent id. The data itself is stored in the generic
 * #NcmCatalog backend; this class only adds the schema metadata and convenience
 * accessors to follow the id/parent-id linkage between related catalogs (for
 * example clusters pointing back to their halos, or members to their host).
 *
 * It does not generate or sample anything: it is a container assembled by
 * higher-level code from the existing sampling machinery.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"
#include "nc/lss/halo/nc_halo_catalog.h"
#include "ncm/data/ncm_catalog.h"
#include "ncm/core/ncm_cfg.h"
#include "nc_enum_types.h"

/* *INDENT-OFF* */
G_DEFINE_QUARK (nc-halo-catalog-error, nc_halo_catalog_error)
/* *INDENT-ON* */

typedef struct _NcHaloCatalogPrivate
{
  NcHaloCatalogKind kind;
  gchar *id_col;
  gchar *parent_id_col;
} NcHaloCatalogPrivate;

enum
{
  PROP_0,
  PROP_KIND,
  PROP_ID_COL,
  PROP_PARENT_ID_COL,
};

struct _NcHaloCatalog
{
  NcmCatalog parent_instance;
};

G_DEFINE_TYPE_WITH_PRIVATE (NcHaloCatalog, nc_halo_catalog, NCM_TYPE_CATALOG)

static void
nc_halo_catalog_init (NcHaloCatalog *hcat)
{
  NcHaloCatalogPrivate * const self = nc_halo_catalog_get_instance_private (hcat);

  self->kind          = NC_HALO_CATALOG_KIND_HALO;
  self->id_col        = NULL;
  self->parent_id_col = NULL;
}

static void
_nc_halo_catalog_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcHaloCatalog *hcat               = NC_HALO_CATALOG (object);
  NcHaloCatalogPrivate * const self = nc_halo_catalog_get_instance_private (hcat);

  switch (prop_id)
  {
    case PROP_KIND:
      self->kind = g_value_get_enum (value);
      break;
    case PROP_ID_COL:
      g_free (self->id_col);
      self->id_col = g_value_dup_string (value);
      break;
    case PROP_PARENT_ID_COL:
      g_free (self->parent_id_col);
      self->parent_id_col = g_value_dup_string (value);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_halo_catalog_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcHaloCatalog *hcat               = NC_HALO_CATALOG (object);
  NcHaloCatalogPrivate * const self = nc_halo_catalog_get_instance_private (hcat);

  switch (prop_id)
  {
    case PROP_KIND:
      g_value_set_enum (value, self->kind);
      break;
    case PROP_ID_COL:
      g_value_set_string (value, self->id_col);
      break;
    case PROP_PARENT_ID_COL:
      g_value_set_string (value, self->parent_id_col);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_halo_catalog_finalize (GObject *object)
{
  NcHaloCatalog *hcat               = NC_HALO_CATALOG (object);
  NcHaloCatalogPrivate * const self = nc_halo_catalog_get_instance_private (hcat);

  g_clear_pointer (&self->id_col, g_free);
  g_clear_pointer (&self->parent_id_col, g_free);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_halo_catalog_parent_class)->finalize (object);
}

static void
nc_halo_catalog_class_init (NcHaloCatalogClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = &_nc_halo_catalog_set_property;
  object_class->get_property = &_nc_halo_catalog_get_property;
  object_class->finalize     = &_nc_halo_catalog_finalize;

  /**
   * NcHaloCatalog:kind:
   *
   * The kind of objects stored in the catalog.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_KIND,
                                   g_param_spec_enum ("kind",
                                                      "Kind",
                                                      "Kind of objects stored in the catalog",
                                                      NC_TYPE_HALO_CATALOG_KIND,
                                                      NC_HALO_CATALOG_KIND_HALO,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_STRINGS));

  /**
   * NcHaloCatalog:id-col:
   *
   * Name of the column holding each row's own id, or %NULL when the catalog has
   * no id column.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_ID_COL,
                                   g_param_spec_string ("id-col",
                                                        "Id column",
                                                        "Name of the column holding each row id",
                                                        NULL,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_STRINGS));

  /**
   * NcHaloCatalog:parent-id-col:
   *
   * Name of the column holding each row's parent id, or %NULL when the catalog
   * has no parent linkage (for example a top-level halo catalog).
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_PARENT_ID_COL,
                                   g_param_spec_string ("parent-id-col",
                                                        "Parent id column",
                                                        "Name of the column holding each row parent id",
                                                        NULL,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_STRINGS));
}

/**
 * nc_halo_catalog_new:
 * @kind: the #NcHaloCatalogKind of the catalog
 * @id_col: (nullable): the name of the id column, or %NULL
 * @parent_id_col: (nullable): the name of the parent id column, or %NULL
 * @nrows: the number of rows
 * @col_names: (array zero-terminated=1): the column names
 * @col_types: (array length=n_types) (nullable): the per-column logical types, or %NULL for all-double
 * @n_types: the number of entries in @col_types
 *
 * Creates a new #NcHaloCatalog with @nrows rows and the columns named by
 * @col_names, all values initialized to zero. @id_col and @parent_id_col, when
 * given, name the columns used to follow the linkage between catalogs.
 *
 * Returns: (transfer full): a new #NcHaloCatalog.
 */
NcHaloCatalog *
nc_halo_catalog_new (NcHaloCatalogKind kind, const gchar *id_col, const gchar *parent_id_col, guint nrows, GStrv col_names, const NcmCatalogColType *col_types, guint n_types)
{
  NcHaloCatalog *hcat;
  GVariant *var = NULL;

  if (col_types != NULL)
  {
    GArray *types = g_array_sized_new (FALSE, TRUE, sizeof (guint32), n_types);
    guint i;

    for (i = 0; i < n_types; i++)
    {
      const guint32 t = col_types[i];

      g_array_append_val (types, t);
    }

    var = ncm_cfg_array_to_variant (types, G_VARIANT_TYPE ("u"));
    g_array_unref (types);
  }

  hcat = g_object_new (NC_TYPE_HALO_CATALOG,
                       "kind", kind,
                       "id-col", id_col,
                       "parent-id-col", parent_id_col,
                       "columns", col_names,
                       "col-types", var,
                       "len", nrows,
                       NULL);

  g_clear_pointer (&var, g_variant_unref);

  return hcat;
}

/**
 * nc_halo_catalog_ref:
 * @hcat: a #NcHaloCatalog
 *
 * Increases the reference count of @hcat by one.
 *
 * Returns: (transfer full): @hcat.
 */
NcHaloCatalog *
nc_halo_catalog_ref (NcHaloCatalog *hcat)
{
  return g_object_ref (hcat);
}

/**
 * nc_halo_catalog_free:
 * @hcat: a #NcHaloCatalog
 *
 * Decreases the reference count of @hcat by one.
 *
 */
void
nc_halo_catalog_free (NcHaloCatalog *hcat)
{
  g_object_unref (hcat);
}

/**
 * nc_halo_catalog_clear:
 * @hcat: a #NcHaloCatalog
 *
 * If *@hcat is different from %NULL, decreases its reference count and sets
 * *@hcat to %NULL.
 *
 */
void
nc_halo_catalog_clear (NcHaloCatalog **hcat)
{
  g_clear_object (hcat);
}

/**
 * nc_halo_catalog_get_kind:
 * @hcat: a #NcHaloCatalog
 *
 * Gets the kind of objects stored in @hcat.
 *
 * Returns: the #NcHaloCatalogKind of the catalog.
 */
NcHaloCatalogKind
nc_halo_catalog_get_kind (NcHaloCatalog *hcat)
{
  NcHaloCatalogPrivate * const self = nc_halo_catalog_get_instance_private (hcat);

  return self->kind;
}

/**
 * nc_halo_catalog_peek_id_col:
 * @hcat: a #NcHaloCatalog
 *
 * Gets the name of the id column.
 *
 * Returns: (transfer none) (nullable): the id column name, or %NULL.
 */
const gchar *
nc_halo_catalog_peek_id_col (NcHaloCatalog *hcat)
{
  NcHaloCatalogPrivate * const self = nc_halo_catalog_get_instance_private (hcat);

  return self->id_col;
}

/**
 * nc_halo_catalog_peek_parent_id_col:
 * @hcat: a #NcHaloCatalog
 *
 * Gets the name of the parent id column.
 *
 * Returns: (transfer none) (nullable): the parent id column name, or %NULL.
 */
const gchar *
nc_halo_catalog_peek_parent_id_col (NcHaloCatalog *hcat)
{
  NcHaloCatalogPrivate * const self = nc_halo_catalog_get_instance_private (hcat);

  return self->parent_id_col;
}

/**
 * nc_halo_catalog_get_id:
 * @hcat: a #NcHaloCatalog
 * @i: a row index
 * @error: a #GError for error reporting
 *
 * Gets the id of row @i, read from the configured id column.
 *
 * Returns: the id of row @i, or -1 on error.
 */
gint64
nc_halo_catalog_get_id (NcHaloCatalog *hcat, const guint i, GError **error)
{
  NcHaloCatalogPrivate * const self = nc_halo_catalog_get_instance_private (hcat);

  g_return_val_if_fail (error == NULL || *error == NULL, -1);

  if (self->id_col == NULL)
  {
    g_set_error_literal (error, NC_HALO_CATALOG_ERROR, NC_HALO_CATALOG_ERROR_NO_LINKAGE,
                         "No id column configured.");

    return -1;
  }

  return ncm_catalog_get_int (NCM_CATALOG (hcat), self->id_col, i, error);
}

/**
 * nc_halo_catalog_get_parent_id:
 * @hcat: a #NcHaloCatalog
 * @i: a row index
 * @error: a #GError for error reporting
 *
 * Gets the parent id of row @i, read from the configured parent id column.
 *
 * Returns: the parent id of row @i, or -1 on error.
 */
gint64
nc_halo_catalog_get_parent_id (NcHaloCatalog *hcat, const guint i, GError **error)
{
  NcHaloCatalogPrivate * const self = nc_halo_catalog_get_instance_private (hcat);

  g_return_val_if_fail (error == NULL || *error == NULL, -1);

  if (self->parent_id_col == NULL)
  {
    g_set_error_literal (error, NC_HALO_CATALOG_ERROR, NC_HALO_CATALOG_ERROR_NO_LINKAGE,
                         "No parent id column configured.");

    return -1;
  }

  return ncm_catalog_get_int (NCM_CATALOG (hcat), self->parent_id_col, i, error);
}

/**
 * nc_halo_catalog_find_children:
 * @hcat: a #NcHaloCatalog
 * @parent_id: the parent id to match
 * @error: a #GError for error reporting
 *
 * Finds the rows whose parent id equals @parent_id, using the configured parent
 * id column.
 *
 * Returns: (transfer full) (element-type guint): an array of matching row
 *     indices, or %NULL on error.
 */
GArray *
nc_halo_catalog_find_children (NcHaloCatalog *hcat, const gint64 parent_id, GError **error)
{
  NcHaloCatalogPrivate * const self = nc_halo_catalog_get_instance_private (hcat);
  GArray *indices;
  guint len;
  guint i;

  g_return_val_if_fail (error == NULL || *error == NULL, NULL);

  if (self->parent_id_col == NULL)
  {
    g_set_error_literal (error, NC_HALO_CATALOG_ERROR, NC_HALO_CATALOG_ERROR_NO_LINKAGE,
                         "No parent id column configured.");

    return NULL;
  }

  len     = ncm_catalog_len (NCM_CATALOG (hcat));
  indices = g_array_new (FALSE, FALSE, sizeof (guint));

  for (i = 0; i < len; i++)
  {
    GError *local_error = NULL;
    const gint64 pid    = ncm_catalog_get_int (NCM_CATALOG (hcat), self->parent_id_col, i, &local_error);

    if (local_error != NULL)
    {
      g_propagate_error (error, local_error);
      g_array_unref (indices);

      return NULL;
    }

    if (pid == parent_id)
      g_array_append_val (indices, i);
  }

  return indices;
}

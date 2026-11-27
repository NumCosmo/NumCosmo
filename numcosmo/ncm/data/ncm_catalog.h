/***************************************************************************
 *            ncm_catalog.h
 *
 *  Sat Jun 13 14:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_catalog.h
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

#ifndef _NCM_CATALOG_H
#define _NCM_CATALOG_H

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/ncm/algebra/ncm_matrix.h>
#include <numcosmo/ncm/core/ncm_obj_array.h>

G_BEGIN_DECLS

#define NCM_TYPE_CATALOG (ncm_catalog_get_type ())

G_DECLARE_DERIVABLE_TYPE (NcmCatalog, ncm_catalog, NCM, CATALOG, GObject)

struct _NcmCatalogClass
{
  /*< private >*/
  GObjectClass parent_class;

  NcmCatalog *(*shallow_dup) (NcmCatalog *catalog);

  /* Padding to allow 18 virtual functions without breaking ABI. */
  gpointer padding[17];
};

/**
 * NcmCatalogColType:
 * @NCM_CATALOG_COL_TYPE_INVALID: invalid column type.
 * @NCM_CATALOG_COL_TYPE_DOUBLE: a real-valued column, stored in the double matrix.
 * @NCM_CATALOG_COL_TYPE_INT: an integer-valued column, stored in its own array of 64-bit
 *   integers.
 * @NCM_CATALOG_COL_TYPE_BOOL: a boolean-valued column, stored in its own array of #gboolean.
 * @NCM_CATALOG_COL_TYPE_DOUBLE_SEPARATE: a real-valued column, stored in its own array,
 *   outside the double matrix.
 *
 * Type of a catalog column. The DOUBLE columns form one matrix; every other column has an
 * array of its own.
 *
 */
typedef enum _NcmCatalogColType /*< prefix=NCM_CATALOG_COL_TYPE >*/
{
  NCM_CATALOG_COL_TYPE_INVALID = -1,
  NCM_CATALOG_COL_TYPE_DOUBLE  = 0,
  NCM_CATALOG_COL_TYPE_INT,
  NCM_CATALOG_COL_TYPE_BOOL,
  NCM_CATALOG_COL_TYPE_DOUBLE_SEPARATE,
  /* < private > */
  NCM_CATALOG_COL_TYPE_LEN, /*< skip >*/
} NcmCatalogColType;

/* Error domain for the NcmCatalog class */

/**
 * NcmCatalogError:
 * @NCM_CATALOG_ERROR_COLUMN_NOT_FOUND: The column was not found.
 * @NCM_CATALOG_ERROR_ROW_OUT_OF_RANGE: The row index is not smaller than the number of rows.
 * @NCM_CATALOG_ERROR_WRONG_TYPE: The column is not of the type of the accessor.
 *
 * Error codes returned by the #NcmCatalog class.
 *
 */
/**
 * NCM_CATALOG_KEY_PENDING:
 *
 * A saved key that takes the current key of the catalog at the first
 * ncm_catalog_key_is_current(); never a key of the catalog, which starts at one.
 */
#define NCM_CATALOG_KEY_PENDING (0)

typedef enum _NcmCatalogError /*< enum,prefix=NCM_CATALOG_ERROR >*/
{
  NCM_CATALOG_ERROR_COLUMN_NOT_FOUND,
  NCM_CATALOG_ERROR_ROW_OUT_OF_RANGE,
  NCM_CATALOG_ERROR_WRONG_TYPE,
} NcmCatalogError;

GQuark ncm_catalog_error_quark (void);

#define NCM_CATALOG_ERROR (ncm_catalog_error_quark ())


NcmCatalog *ncm_catalog_new (guint nrows, GStrv col_names);
NcmCatalog *ncm_catalog_new_full (guint nrows, GStrv col_names, const NcmCatalogColType *col_types, guint n_types);
NcmCatalog *ncm_catalog_ref (NcmCatalog *catalog);
NcmCatalog *ncm_catalog_shallow_dup (NcmCatalog *catalog);
void ncm_catalog_replace_by_copy (NcmCatalog *catalog, const gchar *col, GError **error);
void ncm_catalog_replace_data_by_copy (NcmCatalog *catalog);

void ncm_catalog_free (NcmCatalog *catalog);
void ncm_catalog_clear (NcmCatalog **catalog);

gboolean ncm_catalog_get_index (NcmCatalog *catalog, const gchar *col, guint *i);
gboolean ncm_catalog_get_data_index (NcmCatalog *catalog, const gchar *col, guint *j);
gboolean ncm_catalog_has_column (NcmCatalog *catalog, const gchar *col);

NcmCatalogColType ncm_catalog_get_col_type (NcmCatalog *catalog, const gchar *col, GError **error);

void ncm_catalog_set (NcmCatalog *catalog, const gchar *col, const guint i, gdouble val, GError **error);
gdouble ncm_catalog_get (NcmCatalog *catalog, const gchar *col, const guint i, GError **error);

void ncm_catalog_set_int (NcmCatalog *catalog, const gchar *col, const guint i, gint64 val, GError **error);
gint64 ncm_catalog_get_int (NcmCatalog *catalog, const gchar *col, const guint i, GError **error);

void ncm_catalog_set_bool (NcmCatalog *catalog, const gchar *col, const guint i, gboolean val, GError **error);
gboolean ncm_catalog_get_bool (NcmCatalog *catalog, const gchar *col, const guint i, GError **error);

void ncm_catalog_set_int_column (NcmCatalog *catalog, const gchar *col, GArray *vals, GError **error);
void ncm_catalog_get_int_column (NcmCatalog *catalog, const gchar *col, GArray **vals, GError **error);
void ncm_catalog_set_bool_column (NcmCatalog *catalog, const gchar *col, GArray *vals, GError **error);
void ncm_catalog_get_bool_column (NcmCatalog *catalog, const gchar *col, GArray **vals, GError **error);
void ncm_catalog_set_double_column (NcmCatalog *catalog, const gchar *col, GArray *vals, GError **error);
void ncm_catalog_get_double_column (NcmCatalog *catalog, const gchar *col, GArray **vals, GError **error);

GStrv ncm_catalog_peek_columns (NcmCatalog *catalog);
NcmMatrix *ncm_catalog_peek_data (NcmCatalog *catalog);
guint64 ncm_catalog_get_key (NcmCatalog *catalog);
void ncm_catalog_mark_outdated (NcmCatalog *catalog);
gchar *ncm_catalog_get_data_hash (NcmCatalog *catalog);
gboolean ncm_catalog_key_is_current (NcmCatalog *catalog, guint64 *saved_key);

guint ncm_catalog_len (NcmCatalog *catalog);
guint ncm_catalog_ncols (NcmCatalog *catalog);

void ncm_catalog_set_meta (NcmCatalog *catalog, NcmVarDict *meta);
NcmVarDict *ncm_catalog_peek_meta (NcmCatalog *catalog);

G_END_DECLS

#endif /* _NCM_CATALOG_H */


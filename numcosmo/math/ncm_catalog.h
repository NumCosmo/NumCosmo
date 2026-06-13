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
#include <numcosmo/math/ncm_matrix.h>

G_BEGIN_DECLS

#define NCM_TYPE_CATALOG (ncm_catalog_get_type ())

G_DECLARE_DERIVABLE_TYPE (NcmCatalog, ncm_catalog, NCM, CATALOG, GObject)

struct _NcmCatalogClass
{
  /*< private >*/
  GObjectClass parent_class;

  /* Padding to allow 18 virtual functions without breaking ABI. */
  gpointer padding[18];
};

/**
 * NcmCatalogColType:
 * @NCM_CATALOG_COL_TYPE_INVALID: invalid column type.
 * @NCM_CATALOG_COL_TYPE_DOUBLE: a real-valued column.
 * @NCM_CATALOG_COL_TYPE_INT: an integer-valued column (stored losslessly in the double matrix).
 * @NCM_CATALOG_COL_TYPE_BOOL: a boolean-valued column (stored as 0.0/1.0 in the double matrix).
 *
 * Logical type of a catalog column. All values are physically stored in the
 * double backing matrix; the type is schema metadata used to restore the
 * intended dtype when converting to/from external table formats.
 *
 */
typedef enum _NcmCatalogColType
{
  NCM_CATALOG_COL_TYPE_INVALID = -1,
  NCM_CATALOG_COL_TYPE_DOUBLE  = 0,
  NCM_CATALOG_COL_TYPE_INT,
  NCM_CATALOG_COL_TYPE_BOOL,
  /* < private > */
  NCM_CATALOG_COL_TYPE_LEN, /*< skip >*/
} NcmCatalogColType;

/* Error domain for the NcmCatalog class */

/**
 * NcmCatalogError:
 * @NCM_CATALOG_ERROR_COLUMN_NOT_FOUND: The column was not found.
 *
 * Error codes returned by the #NcmCatalog class.
 *
 */
typedef enum _NcmCatalogError
{
  NCM_CATALOG_ERROR_COLUMN_NOT_FOUND,
} NcmCatalogError;

GQuark ncm_catalog_error_quark (void);

#define NCM_CATALOG_ERROR (ncm_catalog_error_quark ())


NcmCatalog *ncm_catalog_new (guint nrows, GStrv col_names);
NcmCatalog *ncm_catalog_new_full (guint nrows, GStrv col_names, const NcmCatalogColType *col_types, guint n_types);
NcmCatalog *ncm_catalog_ref (NcmCatalog *catalog);

void ncm_catalog_free (NcmCatalog *catalog);
void ncm_catalog_clear (NcmCatalog **catalog);

gboolean ncm_catalog_get_index (NcmCatalog *catalog, const gchar *col, guint *i);
gboolean ncm_catalog_has_column (NcmCatalog *catalog, const gchar *col);

NcmCatalogColType ncm_catalog_get_col_type (NcmCatalog *catalog, const gchar *col, GError **error);

void ncm_catalog_set (NcmCatalog *catalog, const gchar *col, const guint i, gdouble val, GError **error);
gdouble ncm_catalog_get (NcmCatalog *catalog, const gchar *col, const guint i, GError **error);

void ncm_catalog_set_int (NcmCatalog *catalog, const gchar *col, const guint i, gint64 val, GError **error);
gint64 ncm_catalog_get_int (NcmCatalog *catalog, const gchar *col, const guint i, GError **error);

void ncm_catalog_set_bool (NcmCatalog *catalog, const gchar *col, const guint i, gboolean val, GError **error);
gboolean ncm_catalog_get_bool (NcmCatalog *catalog, const gchar *col, const guint i, GError **error);

GStrv ncm_catalog_peek_columns (NcmCatalog *catalog);
NcmMatrix *ncm_catalog_peek_data (NcmCatalog *catalog);

guint ncm_catalog_len (NcmCatalog *catalog);
guint ncm_catalog_ncols (NcmCatalog *catalog);

G_END_DECLS

#endif /* _NCM_CATALOG_H */


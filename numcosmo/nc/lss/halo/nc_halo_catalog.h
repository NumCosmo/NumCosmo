/***************************************************************************
 *            nc_halo_catalog.h
 *
 *  Sat Jun 13 16:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_halo_catalog.h
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

#ifndef _NC_HALO_CATALOG_H
#define _NC_HALO_CATALOG_H

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/ncm/data/ncm_catalog.h>

G_BEGIN_DECLS

#define NC_TYPE_HALO_CATALOG (nc_halo_catalog_get_type ())

G_DECLARE_FINAL_TYPE (NcHaloCatalog, nc_halo_catalog, NC, HALO_CATALOG, NcmCatalog)

/**
 * NcHaloCatalogKind:
 * @NC_HALO_CATALOG_KIND_HALO: a catalog of halos (top-level objects).
 * @NC_HALO_CATALOG_KIND_CLUSTER: a catalog of clusters (detections).
 * @NC_HALO_CATALOG_KIND_MEMBER: a catalog of member objects (e.g. galaxies).
 *
 * The kind of objects stored in a #NcHaloCatalog. It is schema metadata only:
 * all kinds share the same generic column storage; the linkage between catalogs
 * is expressed through the id and parent-id columns.
 *
 */
typedef enum _NcHaloCatalogKind
{
  NC_HALO_CATALOG_KIND_HALO = 0,
  NC_HALO_CATALOG_KIND_CLUSTER,
  NC_HALO_CATALOG_KIND_MEMBER,
  /* < private > */
  NC_HALO_CATALOG_KIND_LEN, /*< skip >*/
} NcHaloCatalogKind;

/* Error domain for the NcHaloCatalog class */

/**
 * NcHaloCatalogError:
 * @NC_HALO_CATALOG_ERROR_NO_LINKAGE: the requested linkage column is not configured.
 *
 * Error codes returned by the #NcHaloCatalog class.
 *
 */
typedef enum _NcHaloCatalogError /*< enum,prefix=NC_HALO_CATALOG_ERROR >*/
{
  NC_HALO_CATALOG_ERROR_NO_LINKAGE,
} NcHaloCatalogError;

GQuark nc_halo_catalog_error_quark (void);

#define NC_HALO_CATALOG_ERROR (nc_halo_catalog_error_quark ())

NcHaloCatalog *nc_halo_catalog_new (NcHaloCatalogKind kind, const gchar *id_col, const gchar *parent_id_col, guint nrows, GStrv col_names, const NcmCatalogColType *col_types, guint n_types);
NcHaloCatalog *nc_halo_catalog_ref (NcHaloCatalog *hcat);

void nc_halo_catalog_free (NcHaloCatalog *hcat);
void nc_halo_catalog_clear (NcHaloCatalog **hcat);

NcHaloCatalogKind nc_halo_catalog_get_kind (NcHaloCatalog *hcat);
const gchar *nc_halo_catalog_peek_id_col (NcHaloCatalog *hcat);
const gchar *nc_halo_catalog_peek_parent_id_col (NcHaloCatalog *hcat);

gint64 nc_halo_catalog_get_id (NcHaloCatalog *hcat, const guint i, GError **error);
gint64 nc_halo_catalog_get_parent_id (NcHaloCatalog *hcat, const guint i, GError **error);

GArray *nc_halo_catalog_find_children (NcHaloCatalog *hcat, const gint64 parent_id, GError **error);

G_END_DECLS

#endif /* _NC_HALO_CATALOG_H */

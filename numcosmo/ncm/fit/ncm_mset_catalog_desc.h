/***************************************************************************
 *            ncm_mset_catalog_desc.h
 *
 *  Thu Oct 29 12:00:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_mset_catalog_desc.h
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

#ifndef _NCM_MSET_CATALOG_DESC_H_
#define _NCM_MSET_CATALOG_DESC_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/ncm/model/ncm_mset.h>

G_BEGIN_DECLS

#define NCM_TYPE_MSET_CATALOG_DESC (ncm_mset_catalog_desc_get_type ())

G_DECLARE_FINAL_TYPE (NcmMSetCatalogDesc, ncm_mset_catalog_desc, NCM, MSET_CATALOG_DESC, GObject)

/**
 * NcmMSetCatalogDescError:
 * @NCM_MSET_CATALOG_DESC_ERROR_INVALID_MSET: the #NcmMSet is missing or its free-parameter
 *   map is not prepared
 * @NCM_MSET_CATALOG_DESC_ERROR_INVALID_ADD_VALS: the additional-value names and symbols
 *   differ in number
 * @NCM_MSET_CATALOG_DESC_ERROR_INVALID_M2LNP_VAR: the $-2\ln(p)$ column is not one of the
 *   additional values
 * @NCM_MSET_CATALOG_DESC_ERROR_INCOMPATIBLE: two descriptions are not compatible
 * @NCM_MSET_CATALOG_DESC_ERROR_ORDER: two descriptions list the free parameters in different
 *   orders
 *
 * Error codes returned by the #NcmMSetCatalogDesc class.
 *
 */
typedef enum _NcmMSetCatalogDescError /*< prefix=NCM_MSET_CATALOG_DESC_ERROR >*/
{
  NCM_MSET_CATALOG_DESC_ERROR_INVALID_MSET,
  NCM_MSET_CATALOG_DESC_ERROR_INVALID_ADD_VALS,
  NCM_MSET_CATALOG_DESC_ERROR_INVALID_M2LNP_VAR,
  NCM_MSET_CATALOG_DESC_ERROR_INCOMPATIBLE,
  NCM_MSET_CATALOG_DESC_ERROR_ORDER,
} NcmMSetCatalogDescError;

GQuark ncm_mset_catalog_desc_error_quark (void);

#define NCM_MSET_CATALOG_DESC_ERROR (ncm_mset_catalog_desc_error_quark ())

NcmMSetCatalogDesc *ncm_mset_catalog_desc_new (NcmMSet *mset, guint nchains, gboolean weighted, const gchar * const *names, const gchar * const *symbols, gint m2lnp_var, const gchar *run_type, GError **error);
NcmMSetCatalogDesc *ncm_mset_catalog_desc_ref (NcmMSetCatalogDesc *desc);
void ncm_mset_catalog_desc_free (NcmMSetCatalogDesc *desc);
void ncm_mset_catalog_desc_clear (NcmMSetCatalogDesc **desc);

NcmMSet *ncm_mset_catalog_desc_peek_mset (NcmMSetCatalogDesc *desc);
guint ncm_mset_catalog_desc_nchains (NcmMSetCatalogDesc *desc);
gboolean ncm_mset_catalog_desc_weighted (NcmMSetCatalogDesc *desc);
gint ncm_mset_catalog_desc_get_m2lnp_var (NcmMSetCatalogDesc *desc);
const gchar *ncm_mset_catalog_desc_get_run_type (NcmMSetCatalogDesc *desc);

guint ncm_mset_catalog_desc_nadd_vals (NcmMSetCatalogDesc *desc);
guint ncm_mset_catalog_desc_ncols (NcmMSetCatalogDesc *desc);
const gchar *ncm_mset_catalog_desc_col_name (NcmMSetCatalogDesc *desc, guint i);
const gchar *ncm_mset_catalog_desc_col_full_name (NcmMSetCatalogDesc *desc, guint i);
const gchar *ncm_mset_catalog_desc_col_symb (NcmMSetCatalogDesc *desc, guint i);
gboolean ncm_mset_catalog_desc_col_by_name (NcmMSetCatalogDesc *desc, const gchar *name, guint *col_index);

gboolean ncm_mset_catalog_desc_check_compatible (NcmMSetCatalogDesc *desc, NcmMSetCatalogDesc *other, GError **error);
gboolean ncm_mset_catalog_desc_check_same_order (NcmMSetCatalogDesc *desc, NcmMSetCatalogDesc *other, GError **error);

#define NCM_MSET_CATALOG_DESC_WEIGHT_COLNAME "NcmMSetCatalog:Row-weights"
#define NCM_MSET_CATALOG_DESC_WEIGHT_SYMBOL "W"

G_END_DECLS

#endif /* _NCM_MSET_CATALOG_DESC_H_ */


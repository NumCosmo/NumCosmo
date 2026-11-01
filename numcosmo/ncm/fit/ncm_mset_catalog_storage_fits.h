/***************************************************************************
 *            ncm_mset_catalog_storage_fits.h
 *
 *  Sat Oct 31 14:00:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_mset_catalog_storage_fits.h
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

#ifndef _NCM_MSET_CATALOG_STORAGE_FITS_H_
#define _NCM_MSET_CATALOG_STORAGE_FITS_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/ncm/fit/ncm_mset_catalog_storage.h>

G_BEGIN_DECLS

#define NCM_TYPE_MSET_CATALOG_STORAGE_FITS (ncm_mset_catalog_storage_fits_get_type ())

G_DECLARE_FINAL_TYPE (NcmMSetCatalogStorageFits, ncm_mset_catalog_storage_fits, NCM, MSET_CATALOG_STORAGE_FITS, NcmMSetCatalogStorage)

void ncm_mset_catalog_storage_fits_free (NcmMSetCatalogStorageFits *storage_fits);
void ncm_mset_catalog_storage_fits_clear (NcmMSetCatalogStorageFits **storage_fits);
NcmMSetCatalogStorageFits *ncm_mset_catalog_storage_fits_new (const gchar *filename, gboolean read_only, GError **error);
NcmMSetCatalogStorageFits *ncm_mset_catalog_storage_fits_ref (NcmMSetCatalogStorageFits *storage_fits);

const gchar *ncm_mset_catalog_storage_fits_peek_filename (NcmMSetCatalogStorageFits *storage_fits);
gboolean ncm_mset_catalog_storage_fits_get_read_only (NcmMSetCatalogStorageFits *storage_fits);
gint ncm_mset_catalog_storage_fits_get_first_id (NcmMSetCatalogStorageFits *storage_fits, GError **error);

#define NCM_MSET_CATALOG_STORAGE_FITS_NSEGS_LABEL "NSEGS"
#define NCM_MSET_CATALOG_STORAGE_FITS_SEG_START_LABEL "SEGS"
#define NCM_MSET_CATALOG_STORAGE_FITS_SEG_SAMPLER_LABEL "SEGN"
#define NCM_MSET_CATALOG_STORAGE_FITS_SEG_OPTIONS_LABEL "SEGO"

G_END_DECLS

#endif /* _NCM_MSET_CATALOG_STORAGE_FITS_H_ */


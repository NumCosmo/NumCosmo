/***************************************************************************
 *            ncm_mset_catalog_storage_memory.h
 *
 *  Thu Oct 29 12:00:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_mset_catalog_storage_memory.h
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

#ifndef _NCM_MSET_CATALOG_STORAGE_MEMORY_H_
#define _NCM_MSET_CATALOG_STORAGE_MEMORY_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/ncm/fit/ncm_mset_catalog_storage.h>

G_BEGIN_DECLS

#define NCM_TYPE_MSET_CATALOG_STORAGE_MEMORY (ncm_mset_catalog_storage_memory_get_type ())

G_DECLARE_FINAL_TYPE (NcmMSetCatalogStorageMemory, ncm_mset_catalog_storage_memory, NCM, MSET_CATALOG_STORAGE_MEMORY, NcmMSetCatalogStorage)

NcmMSetCatalogStorageMemory *ncm_mset_catalog_storage_memory_new (void);
NcmMSetCatalogStorageMemory *ncm_mset_catalog_storage_memory_ref (NcmMSetCatalogStorageMemory *storage_mem);
void ncm_mset_catalog_storage_memory_free (NcmMSetCatalogStorageMemory *storage_mem);
void ncm_mset_catalog_storage_memory_clear (NcmMSetCatalogStorageMemory **storage_mem);

G_END_DECLS

#endif /* _NCM_MSET_CATALOG_STORAGE_MEMORY_H_ */

/***************************************************************************
 *            nc_halo_catalog_member_generator.h
 *
 *  Sun Jun 14 14:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_halo_catalog_member_generator.h
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

#ifndef _NC_HALO_CATALOG_MEMBER_GENERATOR_H
#define _NC_HALO_CATALOG_MEMBER_GENERATOR_H

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_mset.h>
#include <numcosmo/math/ncm_rng.h>
#include <numcosmo/nc/background/nc_distance.h>
#include <numcosmo/nc/lss/galaxy/nc_galaxy_hod.h>
#include <numcosmo/nc/lss/halo/nc_halo_catalog.h>

G_BEGIN_DECLS

#define NC_TYPE_HALO_CATALOG_MEMBER_GENERATOR (nc_halo_catalog_member_generator_get_type ())

G_DECLARE_FINAL_TYPE (NcHaloCatalogMemberGenerator, nc_halo_catalog_member_generator, NC, HALO_CATALOG_MEMBER_GENERATOR, GObject)

/* Error domain for the NcHaloCatalogMemberGenerator class */

/**
 * NcHaloCatalogMemberGeneratorError:
 * @NC_HALO_CATALOG_MEMBER_GENERATOR_ERROR_MISSING_COLUMN: the host catalog lacks a required column.
 *
 * Error codes returned by the #NcHaloCatalogMemberGenerator class.
 *
 */
typedef enum _NcHaloCatalogMemberGeneratorError /*< enum,prefix=NC_HALO_CATALOG_MEMBER_GENERATOR_ERROR >*/
{
  NC_HALO_CATALOG_MEMBER_GENERATOR_ERROR_MISSING_COLUMN,
} NcHaloCatalogMemberGeneratorError;

GQuark nc_halo_catalog_member_generator_error_quark (void);

#define NC_HALO_CATALOG_MEMBER_GENERATOR_ERROR (nc_halo_catalog_member_generator_error_quark ())

NcHaloCatalogMemberGenerator *nc_halo_catalog_member_generator_new (NcGalaxyHOD *hod);
NcHaloCatalogMemberGenerator *nc_halo_catalog_member_generator_ref (NcHaloCatalogMemberGenerator *memgen);

void nc_halo_catalog_member_generator_free (NcHaloCatalogMemberGenerator *memgen);
void nc_halo_catalog_member_generator_clear (NcHaloCatalogMemberGenerator **memgen);

NcGalaxyHOD *nc_halo_catalog_member_generator_peek_hod (NcHaloCatalogMemberGenerator *memgen);

void nc_halo_catalog_member_generator_set_distance (NcHaloCatalogMemberGenerator *memgen, NcDistance *dist);
NcDistance *nc_halo_catalog_member_generator_peek_distance (NcHaloCatalogMemberGenerator *memgen);

NcHaloCatalog *nc_halo_catalog_member_generator_generate (NcHaloCatalogMemberGenerator *memgen, NcHaloCatalog *host, NcmMSet *mset, NcmRNG *rng, GError **error);

G_END_DECLS

#endif /* _NC_HALO_CATALOG_MEMBER_GENERATOR_H */

/***************************************************************************
 *            nc_halo_catalog_generator.h
 *
 *  Sat Jun 13 17:30 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_halo_catalog_generator.h
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

#ifndef _NC_HALO_CATALOG_GENERATOR_H
#define _NC_HALO_CATALOG_GENERATOR_H

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_mset.h>
#include <numcosmo/math/ncm_rng.h>
#include <numcosmo/math/ncm_sky_footprint.h>
#include <numcosmo/lss/nc_cluster_abundance.h>
#include <numcosmo/lss/nc_halo_catalog.h>

G_BEGIN_DECLS

#define NC_TYPE_HALO_CATALOG_GENERATOR (nc_halo_catalog_generator_get_type ())

G_DECLARE_FINAL_TYPE (NcHaloCatalogGenerator, nc_halo_catalog_generator, NC, HALO_CATALOG_GENERATOR, GObject)

NcHaloCatalogGenerator *nc_halo_catalog_generator_new (NcClusterAbundance *cad);
NcHaloCatalogGenerator *nc_halo_catalog_generator_ref (NcHaloCatalogGenerator *gen);

void nc_halo_catalog_generator_free (NcHaloCatalogGenerator *gen);
void nc_halo_catalog_generator_clear (NcHaloCatalogGenerator **gen);

NcClusterAbundance *nc_halo_catalog_generator_peek_abundance (NcHaloCatalogGenerator *gen);

void nc_halo_catalog_generator_set_footprint (NcHaloCatalogGenerator *gen, NcmSkyFootprint *footprint);
NcmSkyFootprint *nc_halo_catalog_generator_peek_footprint (NcHaloCatalogGenerator *gen);

NcHaloCatalog *nc_halo_catalog_generator_generate (NcHaloCatalogGenerator *gen, NcmMSet *mset, NcmRNG *rng);

G_END_DECLS

#endif /* _NC_HALO_CATALOG_GENERATOR_H */

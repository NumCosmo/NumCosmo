/***************************************************************************
 *            ncm_sky_footprint.h
 *
 *  Sat Jun 13 19:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_sky_footprint.h
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

#ifndef _NCM_SKY_FOOTPRINT_H
#define _NCM_SKY_FOOTPRINT_H

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/ncm/core/ncm_rng.h>

G_BEGIN_DECLS

#define NCM_TYPE_SKY_FOOTPRINT (ncm_sky_footprint_get_type ())

G_DECLARE_DERIVABLE_TYPE (NcmSkyFootprint, ncm_sky_footprint, NCM, SKY_FOOTPRINT, GObject)

struct _NcmSkyFootprintClass
{
  /*< private >*/
  GObjectClass parent_class;

  void (*gen_ra_dec) (NcmSkyFootprint *footprint, NcmRNG *rng, gdouble *ra, gdouble *dec);
  gboolean (*contains) (NcmSkyFootprint *footprint, const gdouble ra, const gdouble dec);
  gdouble (*get_area) (NcmSkyFootprint *footprint);
  gdouble (*density) (NcmSkyFootprint *footprint, const gdouble ra, const gdouble dec);
  gdouble (*ln_density) (NcmSkyFootprint *footprint, const gdouble ra, const gdouble dec);

  /* Padding to allow 18 virtual functions without breaking ABI. */
  gpointer padding[13];
};

NcmSkyFootprint *ncm_sky_footprint_ref (NcmSkyFootprint *footprint);

void ncm_sky_footprint_free (NcmSkyFootprint *footprint);
void ncm_sky_footprint_clear (NcmSkyFootprint **footprint);

void ncm_sky_footprint_gen_ra_dec (NcmSkyFootprint *footprint, NcmRNG *rng, gdouble *ra, gdouble *dec);
gboolean ncm_sky_footprint_contains (NcmSkyFootprint *footprint, const gdouble ra, const gdouble dec);
gdouble ncm_sky_footprint_get_area (NcmSkyFootprint *footprint);
gdouble ncm_sky_footprint_density (NcmSkyFootprint *footprint, const gdouble ra, const gdouble dec);
gdouble ncm_sky_footprint_ln_density (NcmSkyFootprint *footprint, const gdouble ra, const gdouble dec);

G_END_DECLS

#endif /* _NCM_SKY_FOOTPRINT_H */


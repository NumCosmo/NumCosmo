/***************************************************************************
 *            ncm_sky_footprint_rectangular.h
 *
 *  Sat Jun 13 19:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_sky_footprint_rectangular.h
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

#ifndef _NCM_SKY_FOOTPRINT_RECTANGULAR_H
#define _NCM_SKY_FOOTPRINT_RECTANGULAR_H

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_sky_footprint.h>

G_BEGIN_DECLS

#define NCM_TYPE_SKY_FOOTPRINT_RECTANGULAR (ncm_sky_footprint_rectangular_get_type ())

G_DECLARE_FINAL_TYPE (NcmSkyFootprintRectangular, ncm_sky_footprint_rectangular, NCM, SKY_FOOTPRINT_RECTANGULAR, NcmSkyFootprint)

NcmSkyFootprintRectangular *ncm_sky_footprint_rectangular_new (const gdouble ra_min, const gdouble ra_max, const gdouble dec_min, const gdouble dec_max);
NcmSkyFootprintRectangular *ncm_sky_footprint_rectangular_ref (NcmSkyFootprintRectangular *rect);

void ncm_sky_footprint_rectangular_free (NcmSkyFootprintRectangular *rect);
void ncm_sky_footprint_rectangular_clear (NcmSkyFootprintRectangular **rect);

void ncm_sky_footprint_rectangular_set_ra_lim (NcmSkyFootprintRectangular *rect, const gdouble ra_min, const gdouble ra_max);
void ncm_sky_footprint_rectangular_get_ra_lim (NcmSkyFootprintRectangular *rect, gdouble *ra_min, gdouble *ra_max);

void ncm_sky_footprint_rectangular_set_dec_lim (NcmSkyFootprintRectangular *rect, const gdouble dec_min, const gdouble dec_max);
void ncm_sky_footprint_rectangular_get_dec_lim (NcmSkyFootprintRectangular *rect, gdouble *dec_min, gdouble *dec_max);

G_END_DECLS

#endif /* _NCM_SKY_FOOTPRINT_RECTANGULAR_H */

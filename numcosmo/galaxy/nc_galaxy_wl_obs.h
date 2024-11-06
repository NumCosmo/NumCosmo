/***************************************************************************
 *           nc_galaxy_wl_obs.h
 *
 *  Tue Jul 16 06:25:17 2024
 *  Copyright  2024 Caio Lima de Oliveira
 *  <caiooliveiracode@pm.me>
 *  Copyright  2024  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_galaxy_wl_obs.h
 * Copyright (C) 2024 Caio Lima de Oliveira <caiooliveiracode@pm.me>
 * Copyright (C) 2024 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

#ifndef _NC_GALAXY_WL_OBS_H
#define _NC_GALAXY_WL_OBS_H

#include <glib.h>
#include <glib-object.h>
#include <stdarg.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/nc_enum_types.h>
#include <numcosmo/math/ncm_matrix.h>
#include <numcosmo/math/ncm_obj_array.h>
#include <numcosmo/math/ncm_spline.h>
#include <numcosmo/math/ncm_mset.h>

G_BEGIN_DECLS

#define NC_TYPE_GALAXY_WL_OBS (nc_galaxy_wl_obs_get_type ())

G_DECLARE_FINAL_TYPE (NcGalaxyWLObs, nc_galaxy_wl_obs, NC, GALAXY_WL_OBS, GObject)

typedef struct _NcGalaxyWLObsPrivate NcGalaxyWLObsPrivate;

typedef enum _NcGalaxyWLObsCoord
{
  NC_GALAXY_WL_OBS_COORD_CELESTIAL,
  NC_GALAXY_WL_OBS_COORD_EUCLIDEAN,
} NcGalaxyWLObsCoord;

NcGalaxyWLObs *nc_galaxy_wl_obs_new (NcGalaxyWLObsCoord coord, guint nrows, GStrv col_names);
NcGalaxyWLObs *nc_galaxy_wl_obs_ref (NcGalaxyWLObs *obs);

void nc_galaxy_wl_obs_free (NcGalaxyWLObs *obs);
void nc_galaxy_wl_obs_clear (NcGalaxyWLObs **obs);

gboolean nc_galaxy_wl_obs_get_index (NcGalaxyWLObs *obs, const gchar *col, guint *i);

void nc_galaxy_wl_obs_set (NcGalaxyWLObs *obs, const gchar *col, const guint i, gdouble val);
void nc_galaxy_wl_obs_set_pz (NcGalaxyWLObs *obs, const guint i, NcmSpline *pz);

gdouble nc_galaxy_wl_obs_get (NcGalaxyWLObs *obs, const gchar *col, const guint i);
NcmSpline *nc_galaxy_wl_obs_peek_pz (NcGalaxyWLObs *obs, const guint i);

GStrv nc_galaxy_wl_obs_peek_columns (NcGalaxyWLObs *obs);

void nc_galaxy_wl_obs_set_coord (NcGalaxyWLObs *obs, NcGalaxyWLObsCoord coord);
NcGalaxyWLObsCoord nc_galaxy_wl_obs_get_coord (NcGalaxyWLObs *obs);

guint nc_galaxy_wl_obs_len (NcGalaxyWLObs *obs);

G_END_DECLS

#endif /* _NC_GALAXY_WL_OBS_H */


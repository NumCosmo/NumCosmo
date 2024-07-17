/***************************************************************************
 *           nc_galaxy_wl_obs.h
 *
 *  Tue Jul 16 06:25:17 2024
 *  Copyright  2024 Caio Lima de Oliveira
 *  <caiooliveiracode@pm.me>
 ****************************************************************************/
/*
 * nc_galaxy_wl_obs.h
 * Copyright (C) 2024 Caio Lima de Oliveira <caiooliveiracode@pm.me>
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

#ifndef _nc_galaxy_wl_obs_H
#define _nc_galaxy_wl_obs_H

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_matrix.h>
#include <numcosmo/nc_enum_types.h>

G_BEGIN_DECLS

#define NC_TYPE_GALAXY_WL_OBS (nc_galaxy_wl_obs_get_type ())

G_DECLARE_FINAL_TYPE (NcGalaxyWLObs, nc_galaxy_wl_obs, NC, GALAXY_WL_OBS, GObject)

typedef struct _NcGalaxyWLObsPrivate NcGalaxyWLObsPrivate;

typedef enum _NcGalaxyWLObsCoord
{
  NC_GALAXY_WL_OBS_COORD_SKY,
  NC_GALAXY_WL_OBS_COORD_PIXEL,
} NcGalaxyWLObsCoord;


struct _NcGalaxyWLObs
{
  /*< private >*/
  GObject parent_instance;
  NcGalaxyWLObsPrivate *priv;
  NcmMatrix *data;
  NcGalaxyWLObsCoord coord;
  gdouble len;
};

NcGalaxyWLObs *nc_galaxy_wl_obs_new (const guint nrows, NcGalaxyWLObsCoord coord);

void nc_galaxy_wl_obs_set (NcGalaxyWLObs *obs, const guint i, const guint j, gdouble val);
gdouble nc_galaxy_wl_obs_get (NcGalaxyWLObs *obs, const guint i, const guint j);

NcmMatrix *nc_galaxy_wl_obs_get_data (NcGalaxyWLObs *obs);
void nc_galaxy_wl_obs_set_data (NcGalaxyWLObs *obs, NcmMatrix *data);

NcGalaxyWLObsCoord nc_galaxy_wl_obs_get_coord (NcGalaxyWLObs *obs);
void nc_galaxy_wl_obs_set_coord (NcGalaxyWLObs *obs, NcGalaxyWLObsCoord coord);

gdouble nc_galaxy_wl_obs_len (NcGalaxyWLObs *obs);

static void nc_galaxy_wl_obs_dispose (GObject *object);
static void nc_galaxy_wl_obs_finalize (GObject *object);
void nc_galaxy_wl_obs_clear (NcGalaxyWLObs **obs);

G_END_DECLS

#endif /* _nc_galaxy_wl_obs_H */


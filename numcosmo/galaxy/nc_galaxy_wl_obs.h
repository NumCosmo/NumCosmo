/***************************************************************************
 *           nc_galaxy_wl_obs.h
 *
 *  Tue Jul 16 06:25:17 2024
 *  Copyright  2024 Caio Lima de Oliveira, Sandro Dias Pinto Vitenti
 *  <caiooliveiracode@pm.me>, <vitenti@uel.br>
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

typedef struct _NcHICosmo NcHICosmo;
typedef struct _NcHaloDensityProfile NcHaloDensityProfile;
typedef struct _NcHaloPosition NcHaloPosition;
typedef struct _NcWLSurfaceMassDensity NcWLSurfaceMassDensity;
typedef struct _NcGalaxySDPosition NcGalaxySDPosition;
typedef struct _NcGalaxySDObsRedshift NcGalaxySDObsRedshift;
typedef struct _NcGalaxySDShape NcGalaxySDShape;


G_BEGIN_DECLS

#define NC_TYPE_GALAXY_WL_OBS (nc_galaxy_wl_obs_get_type ())

G_DECLARE_FINAL_TYPE (NcGalaxyWLObs, nc_galaxy_wl_obs, NC, GALAXY_WL_OBS, GObject)

typedef struct _NcGalaxyWLObsPrivate NcGalaxyWLObsPrivate;

typedef enum _NcGalaxyWLObsCoord
{
  NC_GALAXY_WL_OBS_COORD_CELESTIAL,
  NC_GALAXY_WL_OBS_COORD_EUCLIDEAN,
} NcGalaxyWLObsCoord;

typedef struct _NcGalaxyWLObsModels NcGalaxyWLObsModels;

/**
 * NcGalaxyWLObsModels:
 * @cosmo: a #NcHICosmo object.
 * @density_profile: a #NcHaloDensityProfile object.
 * @surface_mass_density: a #NcWLSurfaceMassDensity object.
 * @galaxy_position: a #NcGalaxySDPosition object.
 * @galaxy_redshift: a #NcGalaxySDObsRedshift object.
 * @galaxy_shape: a #NcGalaxySDShape object.
 *
 * A structure to store the models used to analyze the weak lensing galaxy samples.
 * This is a simple structure to store the models, it will not handle the memory
 * management of the models. The user must control the reference count of the models.
 *
 */
struct _NcGalaxyWLObsModels
{
  NcHICosmo *cosmo;
  NcHaloDensityProfile *density_profile;
  NcHaloPosition *halo_position;
  NcWLSurfaceMassDensity *surface_mass_density;
  NcGalaxySDPosition *galaxy_position;
  NcGalaxySDObsRedshift *galaxy_redshift;
  NcGalaxySDShape *galaxy_shape;
};

#define NC_TYPE_GALAXY_WL_OBS_MODELS (nc_galaxy_wl_obs_models_get_type ())
GType nc_galaxy_wl_obs_models_get_type (void) G_GNUC_CONST;

NcGalaxyWLObsModels *nc_galaxy_wl_obs_models_new ();
NcGalaxyWLObsModels *nc_galaxy_wl_obs_models_dup (const NcGalaxyWLObsModels *models);
void nc_galaxy_wl_obs_models_free (NcGalaxyWLObsModels *models);


NcGalaxyWLObs *nc_galaxy_wl_obs_new (NcGalaxyWLObsCoord coord, guint nrows, GStrv col_names);
NcGalaxyWLObs *nc_galaxy_wl_obs_ref (NcGalaxyWLObs *obs);

void nc_galaxy_wl_obs_free (NcGalaxyWLObs *obs);
void nc_galaxy_wl_obs_clear (NcGalaxyWLObs **obs);

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


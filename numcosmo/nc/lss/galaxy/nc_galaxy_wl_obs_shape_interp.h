/***************************************************************************
 *            nc_galaxy_wl_obs_shape_interp.h
 *
 *  Wed November 25 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_galaxy_wl_obs_shape_interp.h
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

#ifndef _NC_GALAXY_WL_OBS_SHAPE_INTERP_H_
#define _NC_GALAXY_WL_OBS_SHAPE_INTERP_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/ncm/core/ncm_obj_array.h>
#include <numcosmo/ncm/algebra/ncm_spectral_tensor.h>

G_BEGIN_DECLS

#define NC_TYPE_GALAXY_WL_OBS_SHAPE_INTERP (nc_galaxy_wl_obs_shape_interp_get_type ())

G_DECLARE_FINAL_TYPE (NcGalaxyWLObsShapeInterp, nc_galaxy_wl_obs_shape_interp, NC, GALAXY_WL_OBS_SHAPE_INTERP, GObject)

/**
 * NcGalaxyWLObsShapeInterpError:
 * @NC_GALAXY_WL_OBS_SHAPE_INTERP_ERROR_FILE: A file cannot be created, read, written, renamed or removed, or its name is empty.
 * @NC_GALAXY_WL_OBS_SHAPE_INTERP_ERROR_FORMAT: A file does not hold what it should: it is empty, holds another type of object, refers to another file, or its records are truncated.
 * @NC_GALAXY_WL_OBS_SHAPE_INTERP_ERROR_DATA_HASH: The file holds the interpolants of another data hash.
 * @NC_GALAXY_WL_OBS_SHAPE_INTERP_ERROR_STATE: The call is not valid in the state of the object.
 *
 * Error codes returned by the #NcGalaxyWLObsShapeInterp class.
 *
 */
typedef enum _NcGalaxyWLObsShapeInterpError /*< enum,prefix=NC_GALAXY_WL_OBS_SHAPE_INTERP_ERROR >*/
{
  NC_GALAXY_WL_OBS_SHAPE_INTERP_ERROR_FILE,
  NC_GALAXY_WL_OBS_SHAPE_INTERP_ERROR_FORMAT,
  NC_GALAXY_WL_OBS_SHAPE_INTERP_ERROR_DATA_HASH,
  NC_GALAXY_WL_OBS_SHAPE_INTERP_ERROR_STATE,
} NcGalaxyWLObsShapeInterpError;

GQuark nc_galaxy_wl_obs_shape_interp_error_quark (void);

#define NC_GALAXY_WL_OBS_SHAPE_INTERP_ERROR (nc_galaxy_wl_obs_shape_interp_error_quark ())

NcGalaxyWLObsShapeInterp *nc_galaxy_wl_obs_shape_interp_new (void);
NcGalaxyWLObsShapeInterp *nc_galaxy_wl_obs_shape_interp_new_from_file (const gchar *filename, const gchar *data_hash, NcmObjDictStr *source);
NcGalaxyWLObsShapeInterp *nc_galaxy_wl_obs_shape_interp_create_empty (const gchar *filename, GError **error);
NcGalaxyWLObsShapeInterp *nc_galaxy_wl_obs_shape_interp_ref (NcGalaxyWLObsShapeInterp *si);
void nc_galaxy_wl_obs_shape_interp_free (NcGalaxyWLObsShapeInterp *si);
void nc_galaxy_wl_obs_shape_interp_clear (NcGalaxyWLObsShapeInterp **si);

void nc_galaxy_wl_obs_shape_interp_set_interp (NcGalaxyWLObsShapeInterp *si, const guint i, NcmSpectralTensor *st, GError **error);
NcmSpectralTensor *nc_galaxy_wl_obs_shape_interp_peek_interp (NcGalaxyWLObsShapeInterp *si, const guint i, GError **error);
void nc_galaxy_wl_obs_shape_interp_set_source (NcGalaxyWLObsShapeInterp *si, NcmObjDictStr *source, GError **error);
NcmObjDictStr *nc_galaxy_wl_obs_shape_interp_peek_source (NcGalaxyWLObsShapeInterp *si);
void nc_galaxy_wl_obs_shape_interp_set_data_hash (NcGalaxyWLObsShapeInterp *si, const gchar *data_hash, GError **error);
const gchar *nc_galaxy_wl_obs_shape_interp_peek_data_hash (NcGalaxyWLObsShapeInterp *si);
const gchar *nc_galaxy_wl_obs_shape_interp_peek_filename (NcGalaxyWLObsShapeInterp *si);

void nc_galaxy_wl_obs_shape_interp_close (NcGalaxyWLObsShapeInterp *si, GError **error);

G_END_DECLS

#endif /* _NC_GALAXY_WL_OBS_SHAPE_INTERP_H_ */


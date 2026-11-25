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

NcGalaxyWLObsShapeInterp *nc_galaxy_wl_obs_shape_interp_new (void);
NcGalaxyWLObsShapeInterp *nc_galaxy_wl_obs_shape_interp_new_from_file (const gchar *filename, const gchar *data_hash, NcmObjDictStr *source);
NcGalaxyWLObsShapeInterp *nc_galaxy_wl_obs_shape_interp_ref (NcGalaxyWLObsShapeInterp *si);
void nc_galaxy_wl_obs_shape_interp_free (NcGalaxyWLObsShapeInterp *si);
void nc_galaxy_wl_obs_shape_interp_clear (NcGalaxyWLObsShapeInterp **si);

void nc_galaxy_wl_obs_shape_interp_set_interp (NcGalaxyWLObsShapeInterp *si, const guint i, NcmSpectralTensor *st);
NcmSpectralTensor *nc_galaxy_wl_obs_shape_interp_peek_interp (NcGalaxyWLObsShapeInterp *si, const guint i);
void nc_galaxy_wl_obs_shape_interp_set_source (NcGalaxyWLObsShapeInterp *si, NcmObjDictStr *source);
NcmObjDictStr *nc_galaxy_wl_obs_shape_interp_peek_source (NcGalaxyWLObsShapeInterp *si);
void nc_galaxy_wl_obs_shape_interp_set_data_hash (NcGalaxyWLObsShapeInterp *si, const gchar *data_hash);
const gchar *nc_galaxy_wl_obs_shape_interp_peek_data_hash (NcGalaxyWLObsShapeInterp *si);
const gchar *nc_galaxy_wl_obs_shape_interp_peek_filename (NcGalaxyWLObsShapeInterp *si);

G_END_DECLS

#endif /* _NC_GALAXY_WL_OBS_SHAPE_INTERP_H_ */


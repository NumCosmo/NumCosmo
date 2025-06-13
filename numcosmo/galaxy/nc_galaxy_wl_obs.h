/***************************************************************************
 *           nc_galaxy_wl_obs.h
 *
 *  Tue Jul 16 06:25:17 2024
 *  Copyright  2024 Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 *  Copyright  2024  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_galaxy_wl_obs.h
 * Copyright (C) 2024 Caio Lima de Oliveira <caiolimadeoliveira@pm.me>
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

/**
 * NcGalaxyWLObsCoord:
 * @NC_GALAXY_WL_OBS_COORD_CELESTIAL: Celestial coordinates.
 * @NC_GALAXY_WL_OBS_COORD_EUCLIDEAN: Euclidean coordinates.
 *
 * Coordinate system for the galaxy weak lensing ellipticity data.
 *
 */
typedef enum _NcGalaxyWLObsCoord
{
  NC_GALAXY_WL_OBS_COORD_CELESTIAL,
  NC_GALAXY_WL_OBS_COORD_EUCLIDEAN,
} NcGalaxyWLObsCoord;


/**
 * NcGalaxyWLObsEllipConv:
 * @NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE: Normalization by the quadrupole trace.
 * @NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET: Normalization by the quadrupole trace and determinant.
 *
 * Method used to compute the ellipticity from the quadrupole moment matrix.
 * The quadrupole matrix can be normalized either by its trace or by a combination
 * of trace and determinant. In both cases, the shape can be described by an ellipse
 * parameterized as
 * $$
 * \left(\frac{x\cos\theta + y\sin\theta}{a}\right)^2 +
 * \left(\frac{y\cos\theta - x\sin\theta}{b}\right)^2 = 1,
 * $$
 * where $a$ and $b$ are the semi-major and semi-minor axes, and $\theta$ is the
 * orientation angle.
 *
 * For the trace normalization
 * (@NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE), the ellipticity is commonly defined as
 * $$
 * \chi = \frac{a^2 - b^2}{a^2 + b^2} e^{2i\theta},
 * $$
 * sometimes referred to as the distortion.
 *
 * For the normalization involving both trace and determinant
 * (@NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET), the ellipticity is typically defined as
 * $$
 * \epsilon = \frac{a - b}{a + b} e^{2i\theta},
 * $$
 * which corresponds to the complex ellipticity based on axis ratio.
 *
 */
typedef enum _NcGalaxyWLObsEllipConv
{
  NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE,
  NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET,
  /*< private >*/
  NC_GALAXY_WL_OBS_ELLIP_CONV_LEN
} NcGalaxyWLObsEllipConv;

NcGalaxyWLObs *nc_galaxy_wl_obs_new (NcGalaxyWLObsEllipConv ellip_conv, NcGalaxyWLObsCoord coord, guint nrows, GStrv col_names);
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

NcGalaxyWLObsEllipConv nc_galaxy_wl_obs_get_ellip_conv (NcGalaxyWLObs *obs);

guint nc_galaxy_wl_obs_len (NcGalaxyWLObs *obs);

G_END_DECLS

#endif /* _NC_GALAXY_WL_OBS_H */


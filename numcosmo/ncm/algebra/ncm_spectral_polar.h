/***************************************************************************
 *            ncm_spectral_polar.h
 *
 *  Wed Nov 11 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2026 <vitenti@uel.br>
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

#ifndef _NCM_SPECTRAL_POLAR_H_
#define _NCM_SPECTRAL_POLAR_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/ncm/algebra/ncm_matrix.h>
#include <numcosmo/ncm/algebra/ncm_spectral_polar_panels.h>

G_BEGIN_DECLS

/**
 * NcmSpectralPolarF:
 * @user_data: user data
 * @r: the radius
 * @n: number of angles
 * @theta: (array length=n): the angles
 * @cos_theta: (array length=n): the cosines of @theta
 * @sin_theta: (array length=n): the sines of @theta
 * @fval: (array length=n): $G (r, \theta_j)$, filled in place
 *
 * Function of a radius and an angle, evaluated at one radius and @n angles.
 */
typedef void (*NcmSpectralPolarF) (gpointer user_data, gdouble r, guint n, const gdouble *theta, const gdouble *cos_theta, const gdouble *sin_theta, gdouble *fval);

#define NCM_TYPE_SPECTRAL_POLAR (ncm_spectral_polar_get_type ())
#define NCM_SPECTRAL_POLAR_ERROR (ncm_spectral_polar_error_quark ())

G_DECLARE_FINAL_TYPE (NcmSpectralPolar, ncm_spectral_polar, NCM, SPECTRAL_POLAR, GObject)

/**
 * NcmSpectralPolarError:
 * @NCM_SPECTRAL_POLAR_ERROR_MAX_DEPTH: a panel reached #NcmSpectralPolar:max-depth bisections without passing the test
 *
 * Error codes of ncm_spectral_polar_integrate() and ncm_spectral_polar_build().
 */
typedef enum _NcmSpectralPolarError /*< prefix=NCM_SPECTRAL_POLAR_ERROR >*/
{
  NCM_SPECTRAL_POLAR_ERROR_MAX_DEPTH,
} NcmSpectralPolarError;

GQuark ncm_spectral_polar_error_quark (void);

NcmSpectralPolar *ncm_spectral_polar_new (void);
NcmSpectralPolar *ncm_spectral_polar_ref (NcmSpectralPolar *sp);

void ncm_spectral_polar_free (NcmSpectralPolar *sp);
void ncm_spectral_polar_clear (NcmSpectralPolar **sp);

void ncm_spectral_polar_set_level_min (NcmSpectralPolar *sp, guint level_min);
void ncm_spectral_polar_set_level_cap (NcmSpectralPolar *sp, guint level_cap);
void ncm_spectral_polar_set_max_depth (NcmSpectralPolar *sp, guint max_depth);
void ncm_spectral_polar_set_edge_ratio (NcmSpectralPolar *sp, gdouble edge_ratio);

guint ncm_spectral_polar_get_level_min (NcmSpectralPolar *sp);
guint ncm_spectral_polar_get_level_cap (NcmSpectralPolar *sp);
guint ncm_spectral_polar_get_max_depth (NcmSpectralPolar *sp);
gdouble ncm_spectral_polar_get_edge_ratio (NcmSpectralPolar *sp);

gdouble ncm_spectral_polar_integrate (NcmSpectralPolar *sp, NcmSpectralPolarF G, gpointer user_data, gdouble R, gdouble a, gdouble b, GArray *r_edges, GArray *theta_edges, gdouble reltol, gdouble abstol, guint64 *n_eval, NcmMatrix **panels, GError **error);
void ncm_spectral_polar_build (NcmSpectralPolar *sp, NcmSpectralPolarF F, gpointer user_data, GArray *r_edges, GArray *theta_edges, gdouble reltol, gdouble abstol, NcmSpectralPolarPanels *spp, guint64 *n_eval, GError **error);
void ncm_spectral_polar_peak_edges (NcmSpectralPolar *sp, gdouble R, gdouble r_peak, gdouble theta_peak, gdouble L_r, gdouble L_theta, GArray *r_edges, GArray *theta_edges);

G_END_DECLS

#endif /* _NCM_SPECTRAL_POLAR_H_ */


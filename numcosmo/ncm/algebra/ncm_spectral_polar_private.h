/***************************************************************************
 *            ncm_spectral_polar_private.h
 *
 *  Tue Nov 17 2026
 *  Copyright  2026  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) 2026 Caio Lima de Oliveira <caiolimadeoliveira@pm.me>
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

#ifndef _NCM_SPECTRAL_POLAR_PRIVATE_H_
#define _NCM_SPECTRAL_POLAR_PRIVATE_H_

#include <glib.h>
#include "ncm/algebra/ncm_spectral_polar.h"

G_BEGIN_DECLS

/*
 * NcmSpectralPolarRadialF:
 * @user_data: user data
 * @n: number of radii
 * @r: the radii r_i
 * @rho: set to the factor rho (r_i)
 *
 * A factor rho of the integrand that depends only on the radius, at n radii at once. The
 * integrator multiplies its weights by rho, so the callback only sets its values.
 */
typedef void (*NcmSpectralPolarRadialF) (gpointer user_data, guint n, const gdouble *r, gdouble *rho);

/*
 * NcmSpectralPolarMultiRadialF:
 * @user_data: user data
 * @n: number of radii
 * @r: the radii
 * @n_w: number of weights
 * @w: the n_w x n factors, set
 *
 * The radial factors rho_k of the weights of _ncm_spectral_polar_integrate_radial_multi():
 * sets w[k n + i] to rho_k (r_i).
 */
typedef void (*NcmSpectralPolarMultiRadialF) (gpointer user_data, guint n, const gdouble *r, guint n_w, gdouble *w);

typedef struct _NcmSpectralPolarMulti NcmSpectralPolarMulti;

/*
 * NcmSpectralPolarMultiF:
 * @user_data: user data
 * @r: the radius
 * @n: number of angles
 * @theta: the angles
 * @cos_theta: their cosines
 * @sin_theta: their sines
 * @n_out: number of functions
 * @fval: the n_out x n values, set
 *
 * Several functions at one radius and n angles, for _ncm_spectral_polar_refit_multi():
 * sets fval[q n + j] to F_q (r, theta_j).
 */
typedef void (*NcmSpectralPolarMultiF) (gpointer user_data, gdouble r, guint n, const gdouble *theta, const gdouble *cos_theta, const gdouble *sin_theta, guint n_out, gdouble *fval);

gdouble _ncm_spectral_polar_integrate_radial (NcmSpectralPolar *sp, NcmSpectralPolarF G, NcmSpectralPolarRadialF radial, gpointer user_data, gdouble R, gdouble a, gdouble b, GArray *r_edges, GArray *theta_edges, gdouble reltol, gdouble abstol, guint64 *n_eval, NcmMatrix **panels, GError **error);
NcmSpectralPolarMulti *_ncm_spectral_polar_multi_new (guint level_cap, guint n_w, const gdouble *a, const gdouble *b);
void _ncm_spectral_polar_multi_free (NcmSpectralPolarMulti *multi);
gdouble _ncm_spectral_polar_integrate_radial_multi (NcmSpectralPolar *sp, NcmSpectralPolarF G, NcmSpectralPolarRadialF radial, gpointer user_data, gdouble R, gdouble a, gdouble b, GArray *r_edges, GArray *theta_edges, gdouble reltol, gdouble abstol, NcmSpectralPolarMulti *multi, NcmSpectralPolarMultiRadialF mradial, gdouble *out, guint64 *n_eval, GError **error);
void _ncm_spectral_polar_refit_multi (NcmSpectralPolar *sp, NcmSpectralPolarPanels *src, NcmSpectralPolarMultiF F, gpointer user_data, guint n_out, NcmSpectralPolarPanels **out, guint64 *n_eval);

G_END_DECLS

#endif /* _NCM_SPECTRAL_POLAR_PRIVATE_H_ */


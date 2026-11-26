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
 * NcmSpectralPolarRadialSeries:
 * @ra: lower radial end of the panel
 * @rb: upper radial end of the panel
 * @ta: lower angular end of the panel
 * @tb: upper angular end of the panel
 * @n_coeffs: number of coefficients of the series
 * @offset: index of the first coefficient in the array of coefficients
 *
 * The Chebyshev series in r of the angular integral of G on one accepted panel: with
 * r = (ra + rb) / 2 + x (rb - ra) / 2, the integral of G (r, theta) over [ta, tb] is the sum
 * of c^k T_k (x), k < n_coeffs, with c^k at offset + k of the array of coefficients.
 */
typedef struct _NcmSpectralPolarRadialSeries
{
  gdouble ra;
  gdouble rb;
  gdouble ta;
  gdouble tb;
  guint n_coeffs;
  gsize offset;
} NcmSpectralPolarRadialSeries;

gdouble _ncm_spectral_polar_integrate_radial (NcmSpectralPolar *sp, NcmSpectralPolarF G, NcmSpectralPolarRadialF radial, gpointer user_data, gdouble R, gdouble a, gdouble b, GArray *r_edges, GArray *theta_edges, gdouble reltol, gdouble abstol, guint64 *n_eval, NcmMatrix **panels, GError **error);
gdouble _ncm_spectral_polar_integrate_radial_series (NcmSpectralPolar *sp, NcmSpectralPolarF G, NcmSpectralPolarRadialF radial, gpointer user_data, gdouble R, gdouble a, gdouble b, GArray *r_edges, GArray *theta_edges, gdouble reltol, gdouble abstol, guint64 *n_eval, GArray *series, GArray *coeffs, GError **error);

/*
 * NcmSpectralPolarFactorF:
 * @idx: index of the factor
 * @n: number of radii
 * @r: the radii r_i
 * @f: set to the factor at the radii
 * @user_data: user data
 *
 * One of the factors u_j or v_l of the weights of a NcmSpectralPolarWeightGrid, at n radii
 * at once, without the power law of its end.
 */
typedef void (*NcmSpectralPolarFactorF) (guint idx, guint n, const gdouble *r, gdouble *f, gpointer user_data);

/*
 * NcmSpectralPolarWeightGrid:
 *
 * The radial weights w_jl (r) = r^a_j (R - r)^b_l u_j (r) v_l (r), j < Ka and l < Kb, and the
 * tables to integrate them against the radial series of
 * _ncm_spectral_polar_integrate_radial_series(). One object serves one thread at a time.
 */
typedef struct _NcmSpectralPolarWeightGrid NcmSpectralPolarWeightGrid;

NcmSpectralPolarWeightGrid *_ncm_spectral_polar_weight_grid_new (NcmSpectralPolar *sp, gdouble R, guint Ka, const gdouble *a, guint Kb, const gdouble *b);
void _ncm_spectral_polar_weight_grid_free (NcmSpectralPolarWeightGrid *wg);
void _ncm_spectral_polar_weight_grid_integrate (NcmSpectralPolarWeightGrid *wg, GArray *series, GArray *coeffs, NcmSpectralPolarFactorF fa, NcmSpectralPolarFactorF fb, gpointer user_data, gdouble *res);

G_END_DECLS

#endif /* _NCM_SPECTRAL_POLAR_PRIVATE_H_ */


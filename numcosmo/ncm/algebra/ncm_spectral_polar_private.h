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

gdouble _ncm_spectral_polar_integrate_radial (NcmSpectralPolar *sp, NcmSpectralPolarF G, NcmSpectralPolarRadialF radial, gpointer user_data, gdouble R, gdouble a, gdouble b, GArray *r_edges, GArray *theta_edges, gdouble reltol, gdouble abstol, guint64 *n_eval, NcmMatrix **panels, GError **error);

G_END_DECLS

#endif /* _NCM_SPECTRAL_POLAR_PRIVATE_H_ */


/***************************************************************************
 *            nc_galaxy_shape_factor_moment_series_private.h
 *
 *  Thu Sep 3 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 *  Copyright  2026  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * nc_galaxy_shape_factor_moment_series_private.h
 * Copyright (C) 2026 Sandro Dias Pinto Vitenti <vitenti@uel.br>
 * Copyright (C) 2026 Caio Lima de Oliveira <caiolimadeoliveira@pm.me>
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
 * with this program. If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef _NC_GALAXY_SHAPE_FACTOR_MOMENT_SERIES_PRIVATE_H_
#define _NC_GALAXY_SHAPE_FACTOR_MOMENT_SERIES_PRIVATE_H_

#include <glib.h>
#include "nc/lss/galaxy/nc_galaxy_shape_factor_moment_series.h"
#include "nc/lss/galaxy/nc_galaxy_wl_obs.h"

G_BEGIN_DECLS

/*
 * Shared with NcGalaxyShapeFactorTiltedSeries, which needs the same
 * population-independent shear-map series algebra to build its own target
 * moment series (TILT_SERIES.md sec. 3 / tilt_series.tex sec. 5): the
 * target t(g) = (mu, C_t+mu^2, C_x) this class already computes per galaxy
 * is exactly what the tilt's moment conditions are matched against, so the
 * O(n^3) table build (population-independent, a function only of
 * (ellip-conv, trunc-order)) is shared rather than duplicated. Not part of
 * the public API -- mirrors nc_galaxy_shape_pop_gauss_private.h's own
 * precedent of sharing plain functions between related, non-inheriting
 * classes.
 *
 * Poly2 represents a formal 2D polynomial in (chi_I, chibar_I): c[a*sz+b]
 * is the coefficient of chi_I^a * chibar_I^b. See
 * nc_galaxy_shape_factor_moment_series.c's own top-of-file note for the
 * three-term recursion this builds on and the TRACE_DET map-error history.
 */
typedef struct _NcGalaxyShapeFactorMomentSeriesPoly2
{
  gdouble *c; /* c[a*sz+b]: coefficient of chi_I^a * chibar_I^b */
  guint sz;
} Poly2;

Poly2 _nc_galaxy_shape_factor_moment_series_poly2_new (guint sz);
void _nc_galaxy_shape_factor_moment_series_poly2_clear (Poly2 *p);
gdouble _nc_galaxy_shape_factor_moment_series_poly2_get (const Poly2 *p, guint a, guint b);

/*
 * c[j], j=0..n: the g^j coefficient (as a Poly2 in chi_I, chibar_I) of the
 * forward shear map S(g,chi_I) in the tangential gauge (real g), one
 * recursion per NcGalaxyWLObsEllipConv. @c must point to n+1
 * uninitialized elements (g_new(Poly2, n+1), not g_new0 -- each c[j] is
 * poly2_new()'d here); sz must be >= n+3. The caller frees each c[j] with
 * _nc_galaxy_shape_factor_moment_series_poly2_clear() and then @c itself.
 */
void _nc_galaxy_shape_factor_moment_series_build_c_trace (Poly2 *c, guint n, guint sz);
void _nc_galaxy_shape_factor_moment_series_build_c_trace_det (Poly2 *c, guint n, guint sz);

/*
 * Builds the mean/covariance g-series coefficient tables for the given
 * (ellip-conv, trunc-order) -- see nc_galaxy_shape_factor_moment_series.c's
 * own _nc_galaxy_shape_factor_moment_series_build_tables() for the layout
 * of @tab_m_out/@tab_v_out/@tab_w_out (n_m/n_v rows of n_moments columns
 * each, row-major). Out-params rather than writing into
 * NcGalaxyShapeFactorMomentSeriesPrivate directly, so a caller outside this
 * class can reuse the same build; all three *_out arrays are freshly
 * g_new()'d and owned by the caller (g_free()).
 */
void _nc_galaxy_shape_factor_moment_series_build_tables (guint trunc_order, NcGalaxyWLObsEllipConv ellip_conv,
                                                         guint *n_m_out, guint *n_v_out, guint *n_moments_out,
                                                         gdouble **tab_m_out, gdouble **tab_v_out, gdouble **tab_w_out);

G_END_DECLS

#endif /* _NC_GALAXY_SHAPE_FACTOR_MOMENT_SERIES_PRIVATE_H_ */

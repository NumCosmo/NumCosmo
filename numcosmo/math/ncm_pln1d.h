/***************************************************************************
 *            ncm_pln1d.h
 *
 *  Fri Nov 28 11:24:36 2025
 *  Copyright  2025  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_pln1d.h
 * Copyright (C) 2025 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

#ifndef _NCM_PLN1D_H_
#define _NCM_PLN1D_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_vector.h>

G_BEGIN_DECLS

#define NCM_TYPE_PLN1D (ncm_pln1d_get_type ())

/*
 * NcmPLN1D:
 *
 * A simple Poisson–Lognormal 1D integrator using
 * mode finding (GSL Lambert-W), shifted Gauss–Hermite,
 * or Laplace fallback.
 *
 * Fields are private.
 */
G_DECLARE_FINAL_TYPE (NcmPLN1D, ncm_pln1d, NCM, PLN1D, GObject)

/* Constructors and references */
NcmPLN1D *ncm_pln1d_new (guint gh_order);
NcmPLN1D *ncm_pln1d_ref (NcmPLN1D *pln1d);
void ncm_pln1d_free (NcmPLN1D *pln1d);
void ncm_pln1d_clear (NcmPLN1D **pln1d);

/* Configuration */
void ncm_pln1d_set_order (NcmPLN1D *pln1d, guint gh_order);
guint ncm_pln1d_get_order (NcmPLN1D *pln1d);

gdouble ncm_pln1d_mode (gdouble R, gdouble mu, gdouble sigma);
gdouble ncm_pln1d_eval_p (NcmPLN1D *pln, gdouble R, gdouble mu, gdouble sigma);
gdouble ncm_pln1d_eval_lnp (NcmPLN1D *pln, gdouble R, gdouble mu, gdouble sigma);
gdouble ncm_pln1d_eval_range_sum (NcmPLN1D *pln, guint R_min, guint R_max, gdouble mu, gdouble sigma);
gdouble ncm_pln1d_eval_range_sum_lnp (NcmPLN1D *pln, guint R_min, guint R_max, gdouble mu, gdouble sigma);

G_END_DECLS

#endif /* _NCM_PLN1D_H_ */


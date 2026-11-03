/***************************************************************************
 *            nc_xcor_cov_gauss.h
 *
 *  Tue November 03 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
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

#ifndef _NC_XCOR_COV_GAUSS_H_
#define _NC_XCOR_COV_GAUSS_H_

#include <glib-object.h>
#include <glib.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/ncm/algebra/ncm_matrix.h>

G_BEGIN_DECLS

#define NC_TYPE_XCOR_COV_GAUSS (nc_xcor_cov_gauss_get_type ())

G_DECLARE_FINAL_TYPE (NcXcorCovGauss, nc_xcor_cov_gauss, NC, XCOR_COV_GAUSS, GObject)

NcXcorCovGauss *nc_xcor_cov_gauss_new (guint n_tracers, guint ell_min, guint ell_max, gdouble f_sky);
NcXcorCovGauss *nc_xcor_cov_gauss_ref (NcXcorCovGauss *cov);
void nc_xcor_cov_gauss_free (NcXcorCovGauss *cov);
void nc_xcor_cov_gauss_clear (NcXcorCovGauss **cov);

guint nc_xcor_cov_gauss_get_n_tracers (NcXcorCovGauss *cov);
guint nc_xcor_cov_gauss_get_ell_min (NcXcorCovGauss *cov);
guint nc_xcor_cov_gauss_get_ell_max (NcXcorCovGauss *cov);
guint nc_xcor_cov_gauss_get_n_ell (NcXcorCovGauss *cov);

void nc_xcor_cov_gauss_set_f_sky (NcXcorCovGauss *cov, gdouble f_sky);
gdouble nc_xcor_cov_gauss_get_f_sky (NcXcorCovGauss *cov);

guint nc_xcor_cov_gauss_get_n_pairs (NcXcorCovGauss *cov);
guint nc_xcor_cov_gauss_get_pair_index (NcXcorCovGauss *cov, guint a, guint b);

void nc_xcor_cov_gauss_add_spectrum (NcXcorCovGauss *cov, guint a, guint b, NcmMatrix *window);
guint nc_xcor_cov_gauss_get_n_spectra (NcXcorCovGauss *cov);
guint nc_xcor_cov_gauss_get_n_data (NcXcorCovGauss *cov);

void nc_xcor_cov_gauss_compute (NcXcorCovGauss *cov, NcmMatrix *cl, NcmMatrix *noise);
NcmMatrix *nc_xcor_cov_gauss_peek_matrix (NcXcorCovGauss *cov);

G_END_DECLS

#endif /* _NC_XCOR_COV_GAUSS_H_ */


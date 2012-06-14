/***************************************************************************
 *            util.h
 *
 *  Mon Jul 16 18:02:22 2007
 *  Copyright  2007  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <sandro@isoftware.com.br>
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

#ifndef _NCM_H
#define _NCM_H

#include <glib.h>
#include <gmp.h>
#include <mpfr.h>

G_BEGIN_DECLS

gulong ncm_random_seed (void);
gsl_rng *ncm_get_rng (void);

void ncm_finite_diff_calc_J (NcmModel *model, NcData *data, NcmMatrix *jac);

gdouble ncm_get_bao_Omega_m (NcDistance *dist, NcmModel *model);

gdouble *ncm_smoothd (gdouble *in, size_t N, size_t points, size_t pass);
gboolean ncm_get_uniform_sample (NcmMSet *mset, NcmMSetFunc *func, gdouble x0, gdouble x1, NcmVector *sample);
#ifdef NUMCOSMO_HAVE_FFTW3
gboolean ncm_get_smoothed_uniform_sample (NcmMSet *mset, NcmMSetFunc *func, gdouble x0, gdouble x1, gdouble delta, NcmVector *sample);
#endif /* NUMCOSMO_HAVE_FFTW3 */

gdouble ncm_topology_comoving_a0_lss (guint n, gdouble alpha);
gdouble ncm_topology_sigma_comoving_a0_lss (guint n, gdouble alpha, gdouble sigma_alpha);

void ncm_rational_coarce_double (gdouble x, mpq_t q);

gdouble ncm_sphPlm_x (gint l, gint m, gint order);
gdouble ncm_sphPlm_test_theta (gdouble theta, gint lmax, gint *lmin_data);

gsize ncm_mpfr_out_raw (FILE *stream, mpfr_t op);
gsize ncm_mpfr_inp_raw (mpfr_t rop, FILE *stream);
gsize ncm_mpq_out_raw (FILE *f, mpq_t q);
gsize ncm_mpq_inp_raw (mpq_t q, FILE *f);

void ncm_mpz_inits (mpz_t z, ...);
void ncm_mpz_clears (mpz_t z, ...);

gdouble ncm_sum (gdouble *d, gulong n);
gdouble ncm_numdiff_1 (gsl_function *F, const gdouble x, const gdouble ho, gdouble *err);
gdouble ncm_numdiff_2 (gsl_function *F, gdouble *ofx, const gdouble x, const gdouble ho, gdouble *err);
gdouble ncm_numdiff_2_err (gsl_function *F, gdouble *ofx, const gdouble x, const gdouble ho, gdouble err, gdouble *ferr);
gdouble ncm_sqrt1px_m1 (gdouble x);

gdouble ncm_userdef0_f (NcmModel *model, gdouble z);
gdouble ncm_userdef1_f (NcmModel *model, gdouble z);
gdouble ncm_userdef2_f (NcmModel *model, gdouble z);
gdouble ncm_userdef3_f (NcmModel *model, gdouble z);
gdouble ncm_userdef4_f (NcmModel *model, gdouble z);
gboolean ncm_userdef0_df (NcmModel *model, gdouble z, NcmVector *grad);
gboolean ncm_userdef1_df (NcmModel *model, gdouble z, NcmVector *grad);
gboolean ncm_userdef2_df (NcmModel *model, gdouble z, NcmVector *grad);
gboolean ncm_userdef3_df (NcmModel *model, gdouble z, NcmVector *grad);
gboolean ncm_userdef4_df (NcmModel *model, gdouble z, NcmVector *grad);
gboolean ncm_model_userdef0_df (NcmModel *model, gdouble z, NcmVector *grad);
gboolean ncm_model_userdef1_df (NcmModel *model, gdouble z, NcmVector *grad);
gboolean ncm_model_userdef2_df (NcmModel *model, gdouble z, NcmVector *grad);
gboolean ncm_model_userdef3_df (NcmModel *model, gdouble z, NcmVector *grad);
gboolean ncm_model_userdef4_df (NcmModel *model, gdouble z, NcmVector *grad);

G_END_DECLS

#endif /* _UTIL_H */

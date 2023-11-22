/***************************************************************************
 *            ncm_stats_vec.h
 *
 *  Fri August 02 13:41:15 2013
 *  Copyright  2013  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_stats_vec.h
 * Copyright (C) 2013 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

#ifndef _NCM_STATS_VEC_H_
#define _NCM_STATS_VEC_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_vector.h>
#include <numcosmo/math/ncm_matrix.h>

#ifndef NUMCOSMO_GIR_SCAN
#include <math.h>
#include <complex.h>
#endif /* NUMCOSMO_GIR_SCAN */

G_BEGIN_DECLS

#define NCM_TYPE_STATS_VEC (ncm_stats_vec_get_type ())

G_DECLARE_FINAL_TYPE (NcmStatsVec, ncm_stats_vec, NCM, STATS_VEC, GObject)

/**
 * NcmStatsVecType:
 * @NCM_STATS_VEC_MEAN: Calculates mean only.
 * @NCM_STATS_VEC_VAR: Calculates mean and variance.
 * @NCM_STATS_VEC_COV: Calculates mean, variance and covariance.
 *
 * FIXME
 *
 */
typedef enum
{
  NCM_STATS_VEC_MEAN = 0,
  NCM_STATS_VEC_VAR,
  NCM_STATS_VEC_COV,
  /* < private > */
  NCM_STATS_VEC_TYPES_LEN, /*< skip >*/
} NcmStatsVecType;

typedef void (*NcmStatsVecUpdateFunc) (NcmStatsVec *svec, const gdouble w, NcmVector *x);

/**
 * NcmStatsVecARType:
 * @NCM_STATS_VEC_AR_NONE: Calculates using the required order.
 * @NCM_STATS_VEC_AR_FPE: Uses the FPE criterium to choose the ar order.
 * @NCM_STATS_VEC_AR_AIC: Uses the AIC criterium to choose the ar order.
 * @NCM_STATS_VEC_AR_AICC: Uses the AICc criterium to choose the ar order.
 *
 * FIXME
 *
 */
typedef enum
{
  NCM_STATS_VEC_AR_NONE = 0,
  NCM_STATS_VEC_AR_FPE,
  NCM_STATS_VEC_AR_AIC,
  NCM_STATS_VEC_AR_AICC,
  /* < private > */
  NCM_STATS_VEC_AR_LEN, /*< skip >*/
} NcmStatsVecARType;

GType ncm_stats_vec_get_type (void) G_GNUC_CONST;

NcmStatsVec *ncm_stats_vec_new (guint len, NcmStatsVecType t, gboolean save_x);
NcmStatsVec *ncm_stats_vec_ref (NcmStatsVec *svec);
void ncm_stats_vec_free (NcmStatsVec *svec);
void ncm_stats_vec_clear (NcmStatsVec **svec);

void ncm_stats_vec_reset (NcmStatsVec *svec, gboolean rm_saved);
void ncm_stats_vec_update_weight (NcmStatsVec *svec, gdouble w);

void ncm_stats_vec_append_weight (NcmStatsVec *svec, NcmVector *x, gdouble w, gboolean dup);
void ncm_stats_vec_prepend_weight (NcmStatsVec *svec, NcmVector *x, gdouble w, gboolean dup);

void ncm_stats_vec_append (NcmStatsVec *svec, NcmVector *x, gboolean dup);
void ncm_stats_vec_prepend (NcmStatsVec *svec, NcmVector *x, gboolean dup);
void ncm_stats_vec_append_data (NcmStatsVec *svec, GPtrArray *data, gboolean dup);
void ncm_stats_vec_prepend_data (NcmStatsVec *svec, GPtrArray *data, gboolean dup);

void ncm_stats_vec_enable_quantile (NcmStatsVec *svec, gdouble p);
void ncm_stats_vec_disable_quantile (NcmStatsVec *svec);
gdouble ncm_stats_vec_get_quantile (NcmStatsVec *svec, guint i);
gdouble ncm_stats_vec_get_quantile_spread (NcmStatsVec *svec, guint i);

NcmVector *ncm_stats_vec_get_autocorr (NcmStatsVec *svec, guint p);
NcmVector *ncm_stats_vec_get_subsample_autocorr (NcmStatsVec *svec, guint p, guint subsample);
gdouble ncm_stats_vec_get_autocorr_tau (NcmStatsVec *svec, guint p, const guint max_lag);
gdouble ncm_stats_vec_get_subsample_autocorr_tau (NcmStatsVec *svec, guint p, guint subsample, const guint max_lag);

gboolean ncm_stats_vec_fit_ar_model (NcmStatsVec *svec, guint p, const guint order, NcmStatsVecARType ar_crit, NcmVector **rho, NcmVector **pacf, gdouble *ivar, guint *c_order);
gdouble ncm_stats_vec_ar_ess (NcmStatsVec *svec, guint p, NcmStatsVecARType ar_crit, gdouble *spec0, guint *c_order);
gdouble ncm_stats_vec_estimate_const_break (NcmStatsVec *svec, guint p);

NcmVector *ncm_stats_vec_max_ess_time (NcmStatsVec *svec, const guint ntests, gint *bindex, guint *wp, guint *wp_order, gdouble *wp_ess);
NcmVector *ncm_stats_vec_heidel_diag (NcmStatsVec *svec, const guint ntests, const gdouble pvalue, gint *bindex, guint *wp, guint *wp_order, gdouble *wp_pvalue);
NcmVector *ncm_stats_vec_visual_heidel_diag (NcmStatsVec *svec, const guint p, const guint fi, gdouble *mean, gdouble *var);

GPtrArray *ncm_stats_vec_dup_saved_x (NcmStatsVec *svec);

NcmMatrix *ncm_stats_vec_compute_cov_robust_diag (NcmStatsVec *svec);
NcmMatrix *ncm_stats_vec_compute_cov_robust_ogk (NcmStatsVec *svec);

NcmVector *ncm_stats_vec_peek_x (NcmStatsVec *svec);
void ncm_stats_vec_set (NcmStatsVec *svec, guint i, gdouble x_i);
gdouble ncm_stats_vec_get (NcmStatsVec *svec, guint i);
void ncm_stats_vec_update (NcmStatsVec *svec);
guint ncm_stats_vec_len (NcmStatsVec *svec);
gdouble ncm_stats_vec_get_mean (NcmStatsVec *svec, guint i);
gdouble ncm_stats_vec_get_var (NcmStatsVec *svec, guint i);
gdouble ncm_stats_vec_get_sd (NcmStatsVec *svec, guint i);
gdouble ncm_stats_vec_get_cov (NcmStatsVec *svec, guint i, guint j);
gdouble ncm_stats_vec_get_cor (NcmStatsVec *svec, guint i, guint j);
gdouble ncm_stats_vec_get_weight (NcmStatsVec *svec);
void ncm_stats_vec_get_mean_vector (NcmStatsVec *svec, NcmVector *mean, guint offset);
NcmVector *ncm_stats_vec_peek_mean (NcmStatsVec *svec);
void ncm_stats_vec_get_cov_matrix (NcmStatsVec *svec, NcmMatrix *m, guint offset);
NcmMatrix *ncm_stats_vec_peek_cov_matrix (NcmStatsVec *svec, guint offset);
guint ncm_stats_vec_nrows (NcmStatsVec *svec);
guint ncm_stats_vec_nitens (NcmStatsVec *svec);
NcmVector *ncm_stats_vec_peek_row (NcmStatsVec *svec, guint i);
gdouble ncm_stats_vec_get_param_at (NcmStatsVec *svec, guint i, guint p);

#define NCM_STATS_VEC_HEIDEL_PVAL_COR(pvalue, n) (1.0 - pow (1.0 - (pvalue), 1.0 / ((gdouble) (n))))

G_END_DECLS

#endif /* _NCM_STATS_VEC_H_ */


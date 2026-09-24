/***************************************************************************
 *            ncm_stats_dist.h
 *
 *  Wed November 07 16:02:25 2018
 *  Copyright  2018  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_stats_dist.h
 * Copyright (C) 2018 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

#ifndef _NCM_STATS_DIST_H_
#define _NCM_STATS_DIST_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/ncm/core/ncm_rng.h>
#include <numcosmo/ncm/algebra/ncm_vector.h>
#include <numcosmo/ncm/algebra/ncm_matrix.h>
#include <numcosmo/ncm/stats/ncm_stats_vec.h>
#include <numcosmo/ncm/stats/ncm_stats_dist_kernel.h>

G_BEGIN_DECLS

#define NCM_TYPE_STATS_DIST (ncm_stats_dist_get_type ())

G_DECLARE_DERIVABLE_TYPE (NcmStatsDist, ncm_stats_dist, NCM, STATS_DIST, GObject)

struct _NcmStatsDistClass
{
  /*< private >*/
  GObjectClass parent_class;

  void (*set_dim) (NcmStatsDist *sd, const guint dim);
  gdouble (*bandwidth) (NcmStatsDist *sd);
  void (*prepare_shapes) (NcmStatsDist *sd, GPtrArray *sample_array);
  void (*prepare_kernels) (NcmStatsDist *sd);
  void (*compute_IM) (NcmStatsDist *sd, NcmMatrix *IM);
  gdouble (*amise) (NcmStatsDist *sd);
  NcmMatrix *(*peek_cov_decomp) (NcmStatsDist *sd, guint i);
  NcmMatrix *(*peek_full_cov_decomp) (NcmStatsDist *sd);
  NcmMatrix *(*peek_full_cov) (NcmStatsDist *sd);
  gdouble (*get_lnnorm) (NcmStatsDist *sd, guint i);
  gdouble (*eval_weights) (NcmStatsDist *sd, NcmVector *weights, NcmVector *x);
  gdouble (*eval_weights_m2lnp) (NcmStatsDist *sd, NcmVector *weights, NcmVector *x);
  void (*eval_weights_m2lnp_vec) (NcmStatsDist *sd, NcmVector *weights, GPtrArray *x_a, NcmVector *m2lnp);
  void (*reset) (NcmStatsDist *sd);

  /* Padding to allow 18 virtual functions without breaking ABI. */
  gpointer padding[6];
};

/**
 * NcmStatsDistCV:
 * @NCM_STATS_DIST_CV_NONE: no cross validation, the rule-of-thumb bandwidth
 * @NCM_STATS_DIST_CV_SPLIT_M2LNP: sample split, bandwidth by the out-of-sample density
 * @NCM_STATS_DIST_CV_LOO: leave-one-out, bandwidth by the AMISE estimate
 * @NCM_STATS_DIST_CV_SPLIT_ACCEPT: sample split, bandwidth by the out-of-sample acceptance estimate
 * @NCM_STATS_DIST_CV_LOO_M2LNP: leave-one-out, bandwidth by the leave-one-out density
 *
 * Cross-validation method to be applied.
 *
 */
typedef enum _NcmStatsDistCV /*< prefix=NCM_STATS_DIST_CV >*/
{
  NCM_STATS_DIST_CV_NONE,
  NCM_STATS_DIST_CV_SPLIT_M2LNP,
  NCM_STATS_DIST_CV_LOO,
  NCM_STATS_DIST_CV_SPLIT_ACCEPT,
  NCM_STATS_DIST_CV_LOO_M2LNP,
  /* < private > */
  NCM_STATS_DIST_CV_LEN, /*< skip >*/
} NcmStatsDistCV;

NcmStatsDist *ncm_stats_dist_ref (NcmStatsDist *sd);
void ncm_stats_dist_free (NcmStatsDist *sd);
void ncm_stats_dist_clear (NcmStatsDist **sd);

void ncm_stats_dist_set_kernel (NcmStatsDist *sd, NcmStatsDistKernel *sdk);
NcmStatsDistKernel *ncm_stats_dist_peek_kernel (NcmStatsDist *sd);
NcmStatsDistKernel *ncm_stats_dist_get_kernel (NcmStatsDist *sd);

guint ncm_stats_dist_get_dim (NcmStatsDist *sd);
guint ncm_stats_dist_get_sample_size (NcmStatsDist *sd);
guint ncm_stats_dist_get_n_kernels (NcmStatsDist *sd);
gdouble ncm_stats_dist_get_href (NcmStatsDist *sd);

void ncm_stats_dist_set_over_smooth (NcmStatsDist *sd, const gdouble over_smooth);
gdouble ncm_stats_dist_get_over_smooth (NcmStatsDist *sd);

void ncm_stats_dist_set_split_frac (NcmStatsDist *sd, const gdouble split_frac);
gdouble ncm_stats_dist_get_split_frac (NcmStatsDist *sd);


void ncm_stats_dist_set_center_shrink (NcmStatsDist *sd, const gboolean center_shrink);
gboolean ncm_stats_dist_get_center_shrink (NcmStatsDist *sd);
void ncm_stats_dist_set_auto_kernel (NcmStatsDist *sd, gboolean auto_kernel);
gboolean ncm_stats_dist_get_auto_kernel (NcmStatsDist *sd);
gdouble ncm_stats_dist_get_center_shrink_factor (NcmStatsDist *sd);
NcmMatrix *ncm_stats_dist_peek_center_shrink_matrix (NcmStatsDist *sd);
void ncm_stats_dist_set_uniform_weights (NcmStatsDist *sd, const gboolean uniform_weights);
gboolean ncm_stats_dist_get_uniform_weights (NcmStatsDist *sd);
void ncm_stats_dist_set_defensive_frac (NcmStatsDist *sd, const gdouble frac);
gdouble ncm_stats_dist_get_defensive_frac (NcmStatsDist *sd);
void ncm_stats_dist_set_defensive_scale (NcmStatsDist *sd, const gdouble scale);
gdouble ncm_stats_dist_get_defensive_scale (NcmStatsDist *sd);
void ncm_stats_dist_set_defensive_nu (NcmStatsDist *sd, const gdouble nu);
gdouble ncm_stats_dist_get_defensive_nu (NcmStatsDist *sd);

void ncm_stats_dist_set_print_fit (NcmStatsDist *sd, const gboolean print_fit);
gboolean ncm_stats_dist_get_print_fit (NcmStatsDist *sd);

void ncm_stats_dist_set_cv_type (NcmStatsDist *sd, const NcmStatsDistCV cv_type);
NcmStatsDistCV ncm_stats_dist_get_cv_type (NcmStatsDist *sd);

void ncm_stats_dist_set_use_threads (NcmStatsDist *sd, const gboolean use_threads);
gboolean ncm_stats_dist_get_use_threads (NcmStatsDist *sd);

void ncm_stats_dist_prepare_shapes (NcmStatsDist *sd, GPtrArray *sample_array);
void ncm_stats_dist_prepare (NcmStatsDist *sd, NcmVector *m2lnL);

gdouble ncm_stats_dist_eval (NcmStatsDist *sd, NcmVector *x);
gdouble ncm_stats_dist_eval_m2lnp (NcmStatsDist *sd, NcmVector *x);
void ncm_stats_dist_eval_m2lnp_vec (NcmStatsDist *sd, GPtrArray *x_a, NcmVector *m2lnp);

guint ncm_stats_dist_kernel_choose (NcmStatsDist *sd, NcmRNG *rng);
void ncm_stats_dist_sample (NcmStatsDist *sd, NcmVector *x, NcmRNG *rng);

gdouble ncm_stats_dist_get_rnorm (NcmStatsDist *sd);

void ncm_stats_dist_add_obs (NcmStatsDist *sd, NcmVector *y);

GPtrArray *ncm_stats_dist_peek_sample_array (NcmStatsDist *sd);
GPtrArray *ncm_stats_dist_peek_center_array (NcmStatsDist *sd);
NcmMatrix *ncm_stats_dist_peek_cov_decomp (NcmStatsDist *sd, guint i);
NcmMatrix *ncm_stats_dist_peek_full_cov_decomp (NcmStatsDist *sd);
NcmMatrix *ncm_stats_dist_peek_full_cov (NcmStatsDist *sd);
gdouble ncm_stats_dist_get_lnnorm (NcmStatsDist *sd, guint i);
NcmVector *ncm_stats_dist_peek_weights (NcmStatsDist *sd);

void ncm_stats_dist_get_Ki (NcmStatsDist *sd, const guint i, NcmVector **y_i, NcmMatrix **cov_i, gdouble *n_i, gdouble *w_i);

void ncm_stats_dist_reset (NcmStatsDist *sd);

G_END_DECLS

#endif /* _NCM_STATS_DIST_H_ */


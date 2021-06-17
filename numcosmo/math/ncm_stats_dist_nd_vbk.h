/* -*- Mode: C; indent-tabs-mode: nil; c-basic-offset: 2; tab-width: 2 -*-  */
/***************************************************************************
 *            ncm_stats_dist_nd_vbk.h
 *
 *  Wed November 07 16:02:25 2018
 *  Copyright  2018  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * ncm_stats_dist_nd_vbk.h
 * Copyright (C) 2018 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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

#ifndef _NCM_STATS_DIST_ND_VBK_H_
#define _NCM_STATS_DIST_ND_VBK_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_rng.h>
#include <numcosmo/math/ncm_vector.h>
#include <numcosmo/math/ncm_matrix.h>
#include <numcosmo/math/ncm_stats_vec.h>

G_BEGIN_DECLS

#define NCM_TYPE_STATS_DIST_ND_VBK             (ncm_stats_dist_nd_vbk_get_type ())
#define NCM_STATS_DIST_ND_VBK(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_STATS_DIST_ND_VBK, NcmStatsDistNdVBK))
#define NCM_STATS_DIST_ND_VBK_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_STATS_DIST_ND_VBK, NcmStatsDistNdVBKClass))
#define NCM_IS_STATS_DIST_ND_VBK(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_STATS_DIST_ND_VBK))
#define NCM_IS_STATS_DIST_ND_VBK_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_STATS_DIST_ND_VBK))
#define NCM_STATS_DIST_ND_VBK_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_STATS_DIST_ND_VBK, NcmStatsDistNdVBKClass))

typedef struct _NcmStatsDistNdVBKClass NcmStatsDistNdVBKClass;
typedef struct _NcmStatsDistNdVBK NcmStatsDistNdVBK;
typedef struct _NcmStatsDistNdVBKPrivate NcmStatsDistNdVBKPrivate;

struct _NcmStatsDistNdVBKClass
{
  /*< private >*/ 
  GObjectClass parent_class;
  void (*set_dim) (NcmStatsDistNdVBK *dnd, const guint dim);
  gdouble (*get_rot_bandwidth) (NcmStatsDistNdVBK *dnd, const guint d, const gdouble n);
  gdouble (*get_kernel_lnnorm) (NcmStatsDistNdVBK *dnd, NcmMatrix *cov_decomp, const guint d, const gdouble n, const NcmVector *href);
  void (*prepare_kernel_args) (NcmStatsDistNdVBK *dnd, NcmStatsVec *sample);
  void (*prepare_IM) (NcmStatsDistNdVBK *dnd, GPtrArray *Us, const gint d, const gint n, const NcmVector *href, NcmMatrix *IM, GPtrArray *sample_array, GArray *norm);
  void (*prepare) (NcmStatsDistNdVBK *dnd);
  void (*prepare_interp) (NcmStatsDistNdVBK *dnd, NcmVector *m2lnp);
  gdouble (*eval) (NcmStatsDistNdVBK *dnd, NcmVector *weights, NcmVector *y, GPtrArray *sample_array, const gint d, const gint n, const NcmVector *href, GPtrArray *cov_array, GArray *norm_array);
  gdouble (*eval_m2lnp) (NcmStatsDistNdVBK *dnd, NcmVector *weights, NcmVector *invUy, GPtrArray *invUsample, const gint d, const gint n, const NcmVector *href, GPtrArray *cov_array);
  void (*kernel_sample) (NcmStatsDistNdVBK *dnd, NcmMatrix *cov_decomp, const guint d, NcmVector *y, NcmVector *mu, const NcmVector *href, NcmRNG *rng);
  gdouble (*kernel_eval_m2lnp) (NcmStatsDistNdVBK *dnd, NcmMatrix *cov_decomp, const guint d, NcmVector *x, NcmVector *y, NcmVector *v, const NcmVector *href);
  void (*reset) (NcmStatsDistNdVBK *dnd);
};

struct _NcmStatsDistNdVBK
{
  /*< private >*/
  GObject parent_instance;
  NcmStatsDistNdVBKPrivate *priv;
};

/**
 * NcmStatsDistNdVBKCV:
 * @NCM_STATS_DIST_ND_VBK_CV_NONE: No cross validation
 * @NCM_STATS_DIST_ND_CV_VBK_SPLIT: Sample split cross validation
 * @NCM_STATS_DIST_ND_CV_VBK_SPLIT_FITD: Sample split cross validation fitting all diagonal elements
 *
 * Cross-validation method to be applied.
 *
 */
typedef enum _NcmStatsDistNdVBKCV
{
  NCM_STATS_DIST_ND_VBK_CV_NONE,
  NCM_STATS_DIST_ND_VBK_CV_SPLIT,
  NCM_STATS_DIST_ND_VBK_CV_SPLIT_FITD,
  /* < private > */
  NCM_STATS_DIST_ND_VBK_CV_LEN, /*< skip >*/

} NcmStatsDistNdVBKCV;

GType ncm_stats_dist_nd_vbk_get_type (void) G_GNUC_CONST;

NcmStatsDistNdVBK *ncm_stats_dist_nd_vbk_ref (NcmStatsDistNdVBK *dnd);
void ncm_stats_dist_nd_vbk_free (NcmStatsDistNdVBK *dnd);
void ncm_stats_dist_nd_vbk_clear (NcmStatsDistNdVBK **dnd);

guint ncm_stats_dist_nd_vbk_get_dim (NcmStatsDistNdVBK *dnd);

gdouble ncm_stats_dist_nd_vbk_get_rot_bandwidth (NcmStatsDistNdVBK *dnd, const guint d, const gdouble n);
gdouble ncm_stats_dist_nd_vbk_get_kernel_lnnorm (NcmStatsDistNdVBK *dnd, NcmMatrix *cov_decomp, const guint d, const gdouble n, const NcmVector *href);

void ncm_stats_dist_nd_vbk_set_over_smooth (NcmStatsDistNdVBK *dnd, const gdouble over_smooth);
gdouble ncm_stats_dist_nd_vbk_get_over_smooth (NcmStatsDistNdVBK *dnd);

void ncm_stats_dist_nd_vbk_set_split_frac (NcmStatsDistNdVBK *dnd, const gdouble split_frac);
gdouble ncm_stats_dist_nd_vbk_get_split_frac (NcmStatsDistNdVBK *dnd);

void ncm_stats_dist_nd_vbk_set_nearPD_maxiter (NcmStatsDistNdVBK *dnd, const guint maxiter);
guint ncm_stats_dist_nd_vbk_get_nearPD_maxiter (NcmStatsDistNdVBK *dnd);

void ncm_stats_dist_nd_vbk_set_cv_type (NcmStatsDistNdVBK *dnd, const NcmStatsDistNdVBKCV cv_type);
NcmStatsDistNdVBKCV ncm_stats_dist_nd_vbk_get_cv_type (NcmStatsDistNdVBK *dnd);

void ncm_stats_dist_nd_vbk_prepare (NcmStatsDistNdVBK *dnd);
void ncm_stats_dist_nd_vbk_prepare_interp (NcmStatsDistNdVBK *dnd, NcmVector *m2lnp);
gdouble ncm_stats_dist_nd_vbk_eval (NcmStatsDistNdVBK *dnd, NcmVector *x);
gdouble ncm_stats_dist_nd_vbk_eval_m2lnp (NcmStatsDistNdVBK *dnd, NcmVector *x);
void ncm_stats_dist_nd_vbk_sample (NcmStatsDistNdVBK *dnd, NcmVector *x, NcmRNG *rng);

gdouble ncm_stats_dist_nd_vbk_get_rnorm (NcmStatsDistNdVBK *dnd);

void ncm_stats_dist_nd_vbk_kernel_sample (NcmStatsDistNdVBK *dnd, NcmMatrix *cov_decomp, NcmVector *x, NcmVector *mu, const NcmVector *href, NcmRNG *rng);
gdouble ncm_stats_dist_nd_vbk_kernel_eval_m2lnp (NcmStatsDistNdVBK *dnd, NcmVector *x, NcmVector *y, const NcmVector *href);

void ncm_stats_dist_nd_vbk_add_obs_weight (NcmStatsDistNdVBK *dndg, NcmVector *y, const gdouble w);
void ncm_stats_dist_nd_vbk_add_obs (NcmStatsDistNdVBK *dndg, NcmVector *y);

void ncm_stats_dist_nd_vbk_reset (NcmStatsDistNdVBK *dnd);

G_END_DECLS

#endif /* _NCM_STATS_DIST_ND_H_ */

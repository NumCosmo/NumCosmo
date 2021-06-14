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
#define NCM_STATS_DIST_ND_VBK(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_STATS_DIST_ND_VBK, NcmStatsDistNdVbk))
#define NCM_STATS_DIST_ND_VBK_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_STATS_DIST_ND_VBK, NcmStatsDistNdVbkClass))
#define NCM_IS_STATS_DIST_ND_VBK(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_STATS_DIST_ND_VBK))
#define NCM_IS_STATS_DIST_ND_VBK_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_STATS_DIST_ND_VBK))
#define NCM_STATS_DIST_ND_VBK_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_STATS_DIST_ND_VBK, NcmStatsDistNdVbkClass))

typedef struct _NcmStatsDistNdVbkClass NcmStatsDistNdVbkClass;
typedef struct _NcmStatsDistNdVbk NcmStatsDistNdVbk;
typedef struct _NcmStatsDistNdVbkPrivate NcmStatsDistNdVbkPrivate;

struct _NcmStatsDistNdVbkClass
{
  /*< private >*/ 
  GObjectClass parent_class;
  void (*set_dim) (NcmStatsDistNdVbk *dnd, const guint dim);
  gdouble (*get_rot_bandwidth) (NcmStatsDistNdVbk *dnd, const guint d, const gdouble n);
  gdouble (*get_kernel_lnnorm) (NcmStatsDistNdVbk *dnd, NcmMatrix *cov_decomp, const guint d, const gdouble n, const NcmVector *href);
  void (*prepare_kernel_args) (NcmStatsDistNdVbk *dnd, NcmStatsVec *sample);
  void (*prepare_IM) (NcmStatsDistNdVbk *dnd, GPtrArray *invUsample, const gint d, const gint n, const NcmVector *href, NcmMatrix *IM);
  void (*prepare) (NcmStatsDistNdVbk *dnd);
  void (*prepare_interp) (NcmStatsDistNdVbk *dnd, NcmVector *m2lnp);
  gdouble (*eval) (NcmStatsDistNdVbk *dnd, NcmVector *weights, NcmVector *invUy, GPtrArray *invUsample, const gint d, const gint n, const NcmVector *href, GPtrArray *cov_array);
  gdouble (*eval_m2lnp) (NcmStatsDistNdVbk *dnd, NcmVector *weights, NcmVector *invUy, GPtrArray *invUsample, const gint d, const gint n, const NcmVector *href, GPtrArray *cov_array);
  void (*kernel_sample) (NcmStatsDistNdVbk *dnd, NcmMatrix *cov_decomp, const guint d, NcmVector *y, NcmVector *mu, const NcmVector *href, NcmRNG *rng);
  gdouble (*kernel_eval_m2lnp) (NcmStatsDistNdVbk *dnd, NcmMatrix *cov_decomp, const guint d, NcmVector *x, NcmVector *y, NcmVector *v, const NcmVector *href);
  void (*reset) (NcmStatsDistNdVbk *dnd);
};

struct _NcmStatsDistNdVbk
{
  /*< private >*/
  GObject parent_instance;
  NcmStatsDistNdVbkPrivate *priv;
};

/**
 * NcmStatsDistNdVbkCV:
 * @NCM_STATS_DIST_ND_VBK_CV_NONE: No cross validation
 * @NCM_STATS_DIST_ND_CV_VBK_SPLIT: Sample split cross validation
 * @NCM_STATS_DIST_ND_CV_VBK_SPLIT_FITD: Sample split cross validation fitting all diagonal elements
 *
 * Cross-validation method to be applied.
 *
 */
typedef enum _NcmStatsDistNdVbkCV
{
  NCM_STATS_DIST_ND_VBK_CV_NONE,
  NCM_STATS_DIST_ND_VBK_CV_SPLIT,
  NCM_STATS_DIST_ND_VBK_CV_SPLIT_FITD,
  /* < private > */
  NCM_STATS_DIST_ND_VBK_CV_LEN, /*< skip >*/

} NcmStatsDistNdVbkCV;

GType ncm_stats_dist_nd_vbk_get_type (void) G_GNUC_CONST;

NcmStatsDistNdVbk *ncm_stats_dist_nd_vbk_ref (NcmStatsDistNdVbk *dnd);
void ncm_stats_dist_nd_vbk_free (NcmStatsDistNdVbk *dnd);
void ncm_stats_dist_nd_vbk_clear (NcmStatsDistNdVbk **dnd);

guint ncm_stats_dist_nd_vbk_get_dim (NcmStatsDistNdVbk *dnd);

gdouble ncm_stats_dist_nd_vbk_get_rot_bandwidth (NcmStatsDistNdVbk *dnd, const guint d, const gdouble n);
gdouble ncm_stats_dist_nd_vbk_get_kernel_lnnorm (NcmStatsDistNdVbk *dnd, NcmMatrix *cov_decomp, const guint d, const gdouble n, const NcmVector *href);

void ncm_stats_dist_nd_vbk_set_over_smooth (NcmStatsDistNdVbk *dnd, const gdouble over_smooth);
gdouble ncm_stats_dist_nd_vbk_get_over_smooth (NcmStatsDistNdVbk *dnd);

void ncm_stats_dist_nd_vbk_set_split_frac (NcmStatsDistNdVbk *dnd, const gdouble split_frac);
gdouble ncm_stats_dist_nd_vbk_get_split_frac (NcmStatsDistNdVbk *dnd);

void ncm_stats_dist_nd_vbk_set_nearPD_maxiter (NcmStatsDistNdVbk *dnd, const guint maxiter);
guint ncm_stats_dist_nd_vbk_get_nearPD_maxiter (NcmStatsDistNdVbk *dnd);

void ncm_stats_dist_nd_vbk_set_cv_type (NcmStatsDistNdVbk *dnd, const NcmStatsDistNdVbkCV cv_type);
NcmStatsDistNdVbkCV ncm_stats_dist_nd_vbk_get_cv_type (NcmStatsDistNdVbk *dnd);

void ncm_stats_dist_nd_vbk_prepare (NcmStatsDistNdVbk *dnd);
void ncm_stats_dist_nd_vbk_prepare_interp (NcmStatsDistNdVbk *dnd, NcmVector *m2lnp);
gdouble ncm_stats_dist_nd_vbk_eval (NcmStatsDistNdVbk *dnd, NcmVector *x);
gdouble ncm_stats_dist_nd_vbk_eval_m2lnp (NcmStatsDistNdVbk *dnd, NcmVector *x);
void ncm_stats_dist_nd_vbk_sample (NcmStatsDistNdVbk *dnd, NcmVector *x, NcmRNG *rng);

gdouble ncm_stats_dist_nd_vbk_get_rnorm (NcmStatsDistNdVbk *dnd);

void ncm_stats_dist_nd_vbk_kernel_sample (NcmStatsDistNdVbk *dnd, NcmVector *x, NcmVector *mu, const NcmVector *href, NcmRNG *rng);
gdouble ncm_stats_dist_nd_vbk_kernel_eval_m2lnp (NcmStatsDistNdVbk *dnd, NcmVector *x, NcmVector *y, const NcmVector *href);

void ncm_stats_dist_nd_vbk_add_obs_weight (NcmStatsDistNdVbk *dndg, NcmVector *y, const gdouble w);
void ncm_stats_dist_nd_vbk_add_obs (NcmStatsDistNdVbk *dndg, NcmVector *y);

void ncm_stats_dist_nd_vbk_reset (NcmStatsDistNdVbk *dnd);

G_END_DECLS

#endif /* _NCM_STATS_DIST_ND_H_ */

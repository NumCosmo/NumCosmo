/* -*- Mode: C; indent-tabs-mode: nil; c-basic-offset: 2; tab-width: 2 -*-  */
/***************************************************************************
 *            ncm_stats_dist_nd.h
 *
 *  Wed November 07 16:02:25 2018
 *  Copyright  2018  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * ncm_stats_dist_nd.h
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

#ifndef _NCM_STATS_DIST_ND_H_
#define _NCM_STATS_DIST_ND_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_rng.h>
#include <numcosmo/math/ncm_vector.h>
#include <numcosmo/math/ncm_matrix.h>
#include <numcosmo/math/ncm_stats_vec.h>

G_BEGIN_DECLS

#define NCM_TYPE_STATS_DIST_ND             (ncm_stats_dist_nd_get_type ())
#define NCM_STATS_DIST_ND(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_STATS_DIST_ND, NcmStatsDistNd))
#define NCM_STATS_DIST_ND_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_STATS_DIST_ND, NcmStatsDistNdClass))
#define NCM_IS_STATS_DIST_ND(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_STATS_DIST_ND))
#define NCM_IS_STATS_DIST_ND_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_STATS_DIST_ND))
#define NCM_STATS_DIST_ND_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_STATS_DIST_ND, NcmStatsDistNdClass))

typedef struct _NcmStatsDistNdClass NcmStatsDistNdClass;
typedef struct _NcmStatsDistNd NcmStatsDistNd;
typedef struct _NcmStatsDistNdPrivate NcmStatsDistNdPrivate;

struct _NcmStatsDistNdClass
{
  /*< private >*/ 
  GObjectClass parent_class;
  void (*set_dim) (NcmStatsDistNd *dnd, const guint dim);
  gdouble (*get_rot_bandwidth) (NcmStatsDistNd *dnd, const guint d, const gdouble n);
  gdouble (*get_kernel_lnnorm) (NcmStatsDistNd *dnd, NcmMatrix *cov_decomp, const guint d, const gdouble n, const gdouble href);
  void (*prepare_kernel_args) (NcmStatsDistNd *dnd, NcmStatsVec *sample);
  void (*prepare_IM) (NcmStatsDistNd *dnd, GPtrArray *invUsample, const gint d, const gint n, const gdouble href, const gdouble href2, NcmMatrix *IM);
  void (*prepare) (NcmStatsDistNd *dnd);
  void (*prepare_interp) (NcmStatsDistNd *dnd, NcmVector *m2lnp);
  gdouble (*eval) (NcmStatsDistNd *dnd, NcmVector *weights, NcmVector *invUy, GPtrArray *invUsample, const gint d, const gint n, const gdouble href, const gdouble href2);
  void (*kernel_sample) (NcmStatsDistNd *dnd, NcmMatrix *cov_decomp, const guint d, NcmVector *y, NcmVector *mu, gdouble href, NcmRNG *rng);
  gdouble (*kernel_eval_m2lnp) (NcmStatsDistNd *dnd, NcmMatrix *cov_decomp, const guint d, NcmVector *x, NcmVector *y, NcmVector *v, const gdouble href, const gdouble href2);
  void (*reset) (NcmStatsDistNd *dnd);
};

struct _NcmStatsDistNd
{
  /*< private >*/
  GObject parent_instance;
  NcmStatsDistNdPrivate *priv;
};

/**
 * NcmStatsDistNdCV:
 * @NCM_STATS_DIST_ND_CV_NONE: No cross validation
 * @NCM_STATS_DIST_ND_CV_SPLIT: Sample split cross validation
 *
 * Cross-validation method to be applied.
 *
 */
typedef enum _NcmStatsDistNdCV
{
  NCM_STATS_DIST_ND_CV_NONE,
  NCM_STATS_DIST_ND_CV_SPLIT,
  /* < private > */
  NCM_STATS_DIST_ND_CV_LEN, /*< skip >*/

} NcmStatsDistNdCV;

GType ncm_stats_dist_nd_get_type (void) G_GNUC_CONST;

NcmStatsDistNd *ncm_stats_dist_nd_ref (NcmStatsDistNd *dnd);
void ncm_stats_dist_nd_free (NcmStatsDistNd *dnd);
void ncm_stats_dist_nd_clear (NcmStatsDistNd **dnd);

guint ncm_stats_dist_nd_get_dim (NcmStatsDistNd *dnd);

gdouble ncm_stats_dist_nd_get_rot_bandwidth (NcmStatsDistNd *dnd, const guint d, const gdouble n);
gdouble ncm_stats_dist_nd_get_kernel_lnnorm (NcmStatsDistNd *dnd, NcmMatrix *cov_decomp, const guint d, const gdouble n, const gdouble href);

void ncm_stats_dist_nd_set_over_smooth (NcmStatsDistNd *dnd, const gdouble over_smooth);
gdouble ncm_stats_dist_nd_get_over_smooth (NcmStatsDistNd *dnd);

void ncm_stats_dist_nd_set_split_frac (NcmStatsDistNd *dnd, const gdouble split_frac);
gdouble ncm_stats_dist_nd_get_split_frac (NcmStatsDistNd *dnd);

void ncm_stats_dist_nd_set_nearPD_maxiter (NcmStatsDistNd *dnd, const guint maxiter);
guint ncm_stats_dist_nd_get_nearPD_maxiter (NcmStatsDistNd *dnd);

void ncm_stats_dist_nd_set_cv_type (NcmStatsDistNd *dnd, const NcmStatsDistNdCV cv_type);
NcmStatsDistNdCV ncm_stats_dist_nd_get_cv_type (NcmStatsDistNd *dnd);

void ncm_stats_dist_nd_prepare (NcmStatsDistNd *dnd);
void ncm_stats_dist_nd_prepare_interp (NcmStatsDistNd *dnd, NcmVector *m2lnp);
gdouble ncm_stats_dist_nd_eval (NcmStatsDistNd *dnd, NcmVector *x);
gdouble ncm_stats_dist_nd_eval_m2lnp (NcmStatsDistNd *dnd, NcmVector *x);
void ncm_stats_dist_nd_sample (NcmStatsDistNd *dnd, NcmVector *x, NcmRNG *rng);

gdouble ncm_stats_dist_nd_get_rnorm (NcmStatsDistNd *dnd);

void ncm_stats_dist_nd_kernel_sample (NcmStatsDistNd *dnd, NcmVector *x, NcmVector *mu, const gdouble href, NcmRNG *rng);
gdouble ncm_stats_dist_nd_kernel_eval_m2lnp (NcmStatsDistNd *dnd, NcmVector *x, NcmVector *y, const gdouble href);

void ncm_stats_dist_nd_add_obs_weight (NcmStatsDistNd *dndg, NcmVector *y, const gdouble w);
void ncm_stats_dist_nd_add_obs (NcmStatsDistNd *dndg, NcmVector *y);

void ncm_stats_dist_nd_reset (NcmStatsDistNd *dnd);

G_END_DECLS

#endif /* _NCM_STATS_DIST_ND_H_ */

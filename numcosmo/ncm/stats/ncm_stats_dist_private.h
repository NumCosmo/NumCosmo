/***************************************************************************
 *            ncm_stats_dist_private.h
 *
 *  Thu July 22 15:12:38 2021
 *  Copyright  2021  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_stats_dist_private.h
 * Copyright (C) 2021 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

#ifndef _NCM_STATS_DIST_PRIVATE_H_
#define _NCM_STATS_DIST_PRIVATE_H_

#include <glib.h>
#include "ncm/stats/ncm_stats_dist.h"
#include "ncm/algebra/ncm_nnls.h"
#include "ncm/core/ncm_memory_pool.h"
#include "ncm/stats/ncm_stats_vec.h"

G_BEGIN_DECLS

/*
 * Centre shrinkage (West's kernel shrinkage): the transform A with
 * A (C + kappa h^2 Sigma-bar) A^T = C, for the sample covariance C = U_C^T U_C and the
 * mean kernel covariance Sigma-bar, so that the mixture's covariance matches the
 * sample's. Two stages: the basis W, W^{-1} and the eigenvalues s once per prepare
 * (shape stage), then A, a = det(A)^{1/d} and Ahat = A / a once per bandwidth (kernel
 * stage). Centres are m + A (x_i - m); kernel factors follow Ahat.
 */
typedef struct _NcmStatsDistShrink
{
  gboolean on;
  NcmMatrix *A;
  NcmMatrix *Ahat;
  gdouble scale;
  gboolean is_isotropic;
  NcmVector *eigval;
  NcmVector *diag;
  NcmMatrix *W;
  NcmMatrix *Winv;
  NcmMatrix *tmp;
  NcmLapackWS *ws;
} NcmStatsDistShrink;

typedef struct _NcmStatsDistPrivate
{
  /*< private >*/
  NcmStatsDistKernel *kernel;
  GPtrArray *sample_array;
  NcmVector *weights;
  NcmVector *wcum;
  gboolean wcum_ready;
  gboolean print_fit;
  gdouble over_smooth;
  NcmStatsDistCV cv_type;
  gboolean use_threads;
  gdouble split_frac;
  gdouble m2lnL_min;
  gdouble href;
  gdouble rnorm;
  gboolean auto_kernel;
  GPtrArray *center_array;
  NcmVector *sample_mean;
  NcmMatrix *sample_decomp;
  NcmMatrix *kernel_cov;
  NcmStatsDistShrink shrink;
  NcmMatrix *refactor_M;
  NcmMatrix *refactor_B;
  gdouble defensive_frac;
  gdouble defensive_scale;
  gdouble defensive_nu;
  NcmStatsDistKernel *defensive_kernel;
  NcmMatrix *defensive_decomp;
  gdouble defensive_lnnorm;
  NcmMemoryPool *mp_dx;
  guint n_obs;
  guint n_kernels;
  guint d;
  NcmNNLS *nnls;
  NcmMatrix *IM;
  NcmVector *target;
  NcmVector *ones;
  NcmVector *cv_m2lnp;
  GPtrArray *cv_x;
  NcmVector *m2lnL;
  NcmVector *cv_w;
  gboolean uniform_weights;
  gboolean fit_weights;
  gboolean cut_weights;
  GArray *kernel_order;
  GArray *kernel_density;
  NcmVector *amise_x1;
  NcmVector *amise_x2;
  NcmStatsVec *amise_stats;
  NcmRNG *amise_rng;
} NcmStatsDistPrivate;

/* Center shrinkage protocol between NcmStatsDist and its subclasses. */
void _ncm_stats_dist_refactor_decomp (NcmStatsDist *sd, NcmMatrix *U0, NcmMatrix *U);
void _ncm_stats_dist_cholesky (NcmMatrix *decomp, const NcmMatrix *cov, const guint maxiter, const gchar *what);
gdouble _ncm_stats_dist_amise (NcmStatsDist *sd);

G_END_DECLS

#endif /* _NCM_STATS_DIST_PRIVATE_H_ */


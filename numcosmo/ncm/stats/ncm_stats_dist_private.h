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

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_multimin.h>
#endif /* NUMCOSMO_GIR_SCAN */

G_BEGIN_DECLS

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
  gdouble shrink;
  gdouble min_m2lnp;
  gdouble max_m2lnp;
  gdouble href;
  gdouble rnorm;
  gboolean center_shrink;
  gboolean auto_kernel;
  GPtrArray *center_array;
  NcmVector *center_mean;
  NcmMatrix *center_C_decomp;
  NcmMatrix *center_mean_cov;
  NcmMatrix *center_A;
  NcmMatrix *center_Ahat;
  gboolean center_Ahat_identity;
  gdouble center_a;
  gdouble defensive_frac;
  gdouble defensive_scale;
  gdouble defensive_nu;
  NcmStatsDistKernel *defensive_kernel;
  NcmMatrix *defensive_decomp;
  gdouble defensive_lnnorm;
  guint n_obs;
  guint n_kernels;
  guint alloc_n_obs;
  guint alloc_n_kernels;
  gboolean alloc_subs;
  guint d;
  GArray *sampling;
  NcmNNLS *nnls;
  NcmMatrix *IM;
  NcmMatrix *sub_IM;
  NcmVector *sub_x;
  NcmVector *f;
  NcmVector *f1;
  NcmVector *cv_m2lnp;
  gdouble *levmar_workz;
  guint levmar_n;
  gsl_multimin_fminimizer *fmin;
  GArray *m2lnp_sort;
  GArray *m2lnp;
  NcmRNG *rng;
} NcmStatsDistPrivate;

/* Center shrinkage protocol between NcmStatsDist and its subclasses. */
void _ncm_stats_dist_center_matrices (NcmStatsDist *sd, NcmMatrix **C_decomp, NcmMatrix **mean_cov);
void _ncm_stats_dist_refactor_decomp (NcmStatsDist *sd, NcmMatrix *U0, NcmMatrix *U);
gboolean _ncm_stats_dist_center_transform_is_identity (NcmStatsDist *sd);
void _ncm_stats_dist_zero_strict_lower (NcmMatrix *U);

G_END_DECLS

#endif /* _NCM_STATS_DIST_PRIVATE_H_ */


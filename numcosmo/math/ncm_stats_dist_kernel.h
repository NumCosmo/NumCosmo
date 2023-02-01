/* -*- Mode: C; indent-tabs-mode: nil; c-basic-offset: 2; tab-width: 2 -*-  */

/***************************************************************************
 *            ncm_stats_dist_kernel.h
 *
 *  Wed November 07 17:41:38 2018
 *  Copyright  2018  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_stats_dist_kernel.h
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

#ifndef _NCM_STATS_DIST_KERNEL_H_
#define _NCM_STATS_DIST_KERNEL_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_vector.h>
#include <numcosmo/math/ncm_matrix.h>

G_BEGIN_DECLS

#define NCM_TYPE_STATS_DIST_KERNEL             (ncm_stats_dist_kernel_get_type ())
#define NCM_STATS_DIST_KERNEL(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_STATS_DIST_KERNEL, NcmStatsDistKernel))
#define NCM_STATS_DIST_KERNEL_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_STATS_DIST_KERNEL, NcmStatsDistKernelClass))
#define NCM_IS_STATS_DIST_KERNEL(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_STATS_DIST_KERNEL))
#define NCM_IS_STATS_DIST_KERNEL_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_STATS_DIST_KERNEL))
#define NCM_STATS_DIST_KERNEL_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_STATS_DIST_KERNEL, NcmStatsDistKernelClass))

typedef struct _NcmStatsDistKernelClass NcmStatsDistKernelClass;
typedef struct _NcmStatsDistKernel NcmStatsDistKernel;
typedef struct _NcmStatsDistKernelPrivate NcmStatsDistKernelPrivate;


struct _NcmStatsDistKernelClass
{
  GObjectClass parent_class;
  
  void (*set_dim) (NcmStatsDistKernel *sdk, const guint dim);
  guint (*get_dim) (NcmStatsDistKernel *sdk);
  gdouble (*get_rot_bandwidth) (NcmStatsDistKernel *sdk, const gdouble n);
  gdouble (*get_lnnorm) (NcmStatsDistKernel *sdk, NcmMatrix *cov_decomp);
  gdouble (*eval_unnorm) (NcmStatsDistKernel *sdk, const gdouble chi2);
  void (*eval_unnorm_vec) (NcmStatsDistKernel *sdk, NcmVector *chi2, NcmVector *Ku);
  void (*eval_sum0_gamma_lambda) (NcmStatsDistKernel *sdk, NcmVector *chi2, NcmVector *weights, NcmVector *lnnorms, gdouble *gamma, gdouble *lambda);
  void (*eval_sum1_gamma_lambda) (NcmStatsDistKernel *sdk, NcmVector *chi2, NcmVector *weights, gdouble lnnorm, gdouble *gamma, gdouble *lambda);
  void (*sample) (NcmStatsDistKernel *sdk, NcmMatrix *cov_decomp, const gdouble href, NcmVector *mu, NcmVector *y, NcmRNG *rng);
};

struct _NcmStatsDistKernel
{
  GObject parent_instance;
  
  NcmStatsDistKernelPrivate *priv;
};

GType ncm_stats_dist_kernel_get_type (void) G_GNUC_CONST;

NcmStatsDistKernel *ncm_stats_dist_kernel_ref (NcmStatsDistKernel *sdk);
void ncm_stats_dist_kernel_free (NcmStatsDistKernel *sdk);
void ncm_stats_dist_kernel_clear (NcmStatsDistKernel **sdk);

guint ncm_stats_dist_kernel_get_dim (NcmStatsDistKernel *sdk);

gdouble ncm_stats_dist_kernel_get_rot_bandwidth (NcmStatsDistKernel *sdk, const gdouble n);
gdouble ncm_stats_dist_kernel_get_lnnorm (NcmStatsDistKernel *sdk, NcmMatrix *cov_decomp);

gdouble ncm_stats_dist_kernel_eval_unnorm (NcmStatsDistKernel *sdk, const gdouble chi2);
void ncm_stats_dist_kernel_eval_unnorm_vec (NcmStatsDistKernel *sdk, NcmVector *chi2, NcmVector *Ku);
void ncm_stats_dist_kernel_eval_sum0_gamma_lambda (NcmStatsDistKernel *sdk, NcmVector *chi2, NcmVector *weights, NcmVector *lnnorms, gdouble *gamma, gdouble *lambda);
void ncm_stats_dist_kernel_eval_sum1_gamma_lambda (NcmStatsDistKernel *sdk, NcmVector *chi2, NcmVector *weights, gdouble lnnorm, gdouble *gamma, gdouble *lambda);

void ncm_stats_dist_kernel_sample (NcmStatsDistKernel *sdk, NcmMatrix *cov_decomp, const gdouble href, NcmVector *mu, NcmVector *y, NcmRNG *rng);

G_END_DECLS

#endif /* _NCM_STATS_DIST_KERNEL_H_ */


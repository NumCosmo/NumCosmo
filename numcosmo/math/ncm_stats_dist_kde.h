/* -*- Mode: C; indent-tabs-mode: nil; c-basic-offset: 2; tab-width: 2 -*-  */

/***************************************************************************
 *            ncm_stats_dist_kde.h
 *
 *  Wed November 07 16:02:25 2018
 *  Copyright  2018  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_stats_dist_kde.h
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

#ifndef _NCM_STATS_DIST_KDE_H_
#define _NCM_STATS_DIST_KDE_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_rng.h>
#include <numcosmo/math/ncm_vector.h>
#include <numcosmo/math/ncm_matrix.h>
#include <numcosmo/math/ncm_stats_vec.h>
#include <numcosmo/math/ncm_stats_dist.h>
#include <numcosmo/math/ncm_stats_dist_kernel.h>

G_BEGIN_DECLS

#define NCM_TYPE_STATS_DIST_KDE (ncm_stats_dist_kde_get_type ())

G_DECLARE_DERIVABLE_TYPE (NcmStatsDistKDE, ncm_stats_dist_kde, NCM, STATS_DIST_KDE, NcmStatsDist)

struct _NcmStatsDistKDEClass
{
  /*< private >*/
  NcmStatsDistClass parent_class;

  /* Padding to allow 18 virtual functions without breaking ABI. */
  gpointer padding[7];
};

/**
 * NcmStatsDistKDECovType:
 * @NCM_STATS_DIST_KDE_COV_TYPE_SAMPLE: Use sample covariance.
 * @NCM_STATS_DIST_KDE_COV_TYPE_FIXED: Use a fixed covariance matrix.
 * @NCM_STATS_DIST_KDE_COV_TYPE_ROBUST_DIAG: Use an 1D robust estimator to build a diagonal covariance.
 * @NCM_STATS_DIST_KDE_COV_TYPE_ROBUST: Use the OGK method to build a covariance.
 *
 * Selects the covariance type to use in the kernel interpolation.
 *
 */
typedef enum _NcmStatsDistKDECovType
{
  NCM_STATS_DIST_KDE_COV_TYPE_SAMPLE,
  NCM_STATS_DIST_KDE_COV_TYPE_FIXED,
  NCM_STATS_DIST_KDE_COV_TYPE_ROBUST_DIAG,
  NCM_STATS_DIST_KDE_COV_TYPE_ROBUST,
  /* < private > */
  NCM_STATS_DIST_KDE_COV_TYPE_LEN, /*< skip >*/
} NcmStatsDistKDECovType;

NcmStatsDistKDE *ncm_stats_dist_kde_new (NcmStatsDistKernel *sdk, NcmStatsDistCV CV_type);
NcmStatsDistKDE *ncm_stats_dist_kde_ref (NcmStatsDistKDE *sdkde);
void ncm_stats_dist_kde_free (NcmStatsDistKDE *sdkde);
void ncm_stats_dist_kde_clear (NcmStatsDistKDE **sdkde);

void ncm_stats_dist_kde_set_nearPD_maxiter (NcmStatsDistKDE *sdkde, const guint maxiter);
guint ncm_stats_dist_kde_get_nearPD_maxiter (NcmStatsDistKDE *sdkde);

void ncm_stats_dist_kde_set_cov_type (NcmStatsDistKDE *sdkde, NcmStatsDistKDECovType cov_type);
NcmStatsDistKDECovType ncm_stats_dist_kde_get_cov_type (NcmStatsDistKDE *sdkde);

void ncm_stats_dist_kde_set_cov_fixed (NcmStatsDistKDE *sdkde, NcmMatrix *cov_fixed);
NcmMatrix *ncm_stats_dist_kde_peek_cov_fixed (NcmStatsDistKDE *sdkde);

G_END_DECLS

#endif /* _NCM_STATS_DIST_KDE_H_ */


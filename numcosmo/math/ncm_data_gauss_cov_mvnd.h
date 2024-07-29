/* -*- Mode: C; indent-tabs-mode: nil; c-basic-offset: 2; tab-width: 2 -*-  */

/***************************************************************************
 *            ncm_data_gauss_cov_mvnd.h
 *
 *  Sun February 04 15:07:40 2018
 *  Copyright  2018  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_data_gauss_cov_mvnd.h
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

#ifndef _NCM_DATA_GAUSS_COV_MVND_H_
#define _NCM_DATA_GAUSS_COV_MVND_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_data_gauss_cov.h>
#include <numcosmo/math/ncm_stats_vec.h>

G_BEGIN_DECLS

#define NCM_TYPE_DATA_GAUSS_COV_MVND (ncm_data_gauss_cov_mvnd_get_type ())
G_DECLARE_FINAL_TYPE (NcmDataGaussCovMVND, ncm_data_gauss_cov_mvnd, NCM, DATA_GAUSS_COV_MVND, NcmDataGaussCov)

typedef gboolean (*NcmDataGaussCovMVNDBound) (gpointer obj, NcmVector *y);

NcmDataGaussCovMVND *ncm_data_gauss_cov_mvnd_new (const guint dim);
NcmDataGaussCovMVND *ncm_data_gauss_cov_mvnd_new_full (const guint dim, const gdouble sigma_min, const gdouble sigma_max, const gdouble cor_level, const gdouble mean_min, const gdouble mean_max, NcmRNG *rng);
NcmDataGaussCovMVND *ncm_data_gauss_cov_mvnd_ref (NcmDataGaussCovMVND *data_mvnd);
void ncm_data_gauss_cov_mvnd_free (NcmDataGaussCovMVND *data_mvnd);
void ncm_data_gauss_cov_mvnd_clear (NcmDataGaussCovMVND **data_mvnd);

void ncm_data_gauss_cov_mvnd_gen_cov_mean (NcmDataGaussCovMVND *data_mvnd, const gdouble sigma_min, const gdouble sigma_max, const gdouble cor_level, const gdouble mean_min, const gdouble mean_max, NcmRNG *rng);
void ncm_data_gauss_cov_mvnd_set_cov_mean (NcmDataGaussCovMVND *data_mvnd, NcmVector *mean, NcmMatrix *cov);
NcmVector *ncm_data_gauss_cov_mvnd_peek_mean (NcmDataGaussCovMVND *data_mvnd);

NcmVector *ncm_data_gauss_cov_mvnd_gen (NcmDataGaussCovMVND *data_mvnd, NcmMSet *mset, gpointer obj, NcmDataGaussCovMVNDBound bound, NcmRNG *rng, gulong *N);
gdouble ncm_data_gauss_cov_mvnd_est_ratio (NcmDataGaussCovMVND *data_mvnd, NcmMSet *mset, gpointer obj, NcmDataGaussCovMVNDBound bound, gulong *N, gulong *Nin, const gdouble reltol, NcmRNG *rng);

void ncm_data_gauss_cov_mvnd_log_info (NcmDataGaussCovMVND *data_mvnd);

NcmStatsVec *ncm_data_gauss_cov_mvnd_stats_vec (NcmDataGaussCovMVND *data_mvnd, NcmMSet *mset, const guint n, const glong maxiter, NcmVector *lower, NcmVector *upper, gboolean save_realizations, NcmRNG *rng);

G_END_DECLS

#endif /* _NCM_DATA_GAUSS_COV_MVND_H_ */


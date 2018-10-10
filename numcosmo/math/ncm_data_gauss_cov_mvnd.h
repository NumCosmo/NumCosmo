/* -*- Mode: C; indent-tabs-mode: nil; c-basic-offset: 2; tab-width: 2 -*-  */
/***************************************************************************
 *            ncm_data_gauss_cov_mvnd.h
 *
 *  Sun February 04 15:07:40 2018
 *  Copyright  2018  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * ncm_data_gauss_cov_mvnd.h
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

#ifndef _NCM_DATA_GAUSS_COV_MVND_H_
#define _NCM_DATA_GAUSS_COV_MVND_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_data_gauss_cov.h>

G_BEGIN_DECLS

#define NCM_TYPE_DATA_GAUSS_COV_MVND             (ncm_data_gauss_cov_mvnd_get_type ())
#define NCM_DATA_GAUSS_COV_MVND(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_DATA_GAUSS_COV_MVND, NcmDataGaussCovMVND))
#define NCM_DATA_GAUSS_COV_MVND_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_DATA_GAUSS_COV_MVND, NcmDataGaussCovMVNDClass))
#define NCM_IS_DATA_GAUSS_COV_MVND(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_DATA_GAUSS_COV_MVND))
#define NCM_IS_DATA_GAUSS_COV_MVND_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_DATA_GAUSS_COV_MVND))
#define NCM_DATA_GAUSS_COV_MVND_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_DATA_GAUSS_COV_MVND, NcmDataGaussCovMVNDClass))

typedef struct _NcmDataGaussCovMVNDClass NcmDataGaussCovMVNDClass;
typedef struct _NcmDataGaussCovMVND NcmDataGaussCovMVND;
typedef struct _NcmDataGaussCovMVNDPrivate NcmDataGaussCovMVNDPrivate;

struct _NcmDataGaussCovMVNDClass
{
  /*< private >*/
  NcmDataGaussCovClass parent_class;
};

struct _NcmDataGaussCovMVND
{
  /*< private >*/
  NcmDataGaussCov parent_instance;
  NcmDataGaussCovMVNDPrivate *priv;
};

typedef gboolean (*NcmDataGaussCovMVNDBound) (gpointer obj, NcmVector *y);

GType ncm_data_gauss_cov_mvnd_get_type (void) G_GNUC_CONST;

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

G_END_DECLS

#endif /* _NCM_DATA_GAUSS_COV_MVND_H_ */

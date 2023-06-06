/***************************************************************************
 *            ncm_data_gauss_cov.h
 *
 *  Tue November 20 18:45:55 2012
 *  Copyright  2012  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) 2012 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

#ifndef _NCM_DATA_GAUSS_COV_H_
#define _NCM_DATA_GAUSS_COV_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_data.h>
#include <numcosmo/math/ncm_bootstrap.h>

G_BEGIN_DECLS

#define NCM_TYPE_DATA_GAUSS_COV (ncm_data_gauss_cov_get_type ())

G_DECLARE_DERIVABLE_TYPE (NcmDataGaussCov, ncm_data_gauss_cov, NCM, DATA_GAUSS_COV, NcmData)

struct _NcmDataGaussCovClass
{
  /*< private >*/
  NcmDataClass parent_class;

  void (*mean_func) (NcmDataGaussCov *gauss, NcmMSet *mset, NcmVector *vp);
  gboolean (*cov_func) (NcmDataGaussCov *gauss, NcmMSet *mset, NcmMatrix *cov);
  void (*lnNorma2) (NcmDataGaussCov *gauss, NcmMSet *mset, gdouble *m2lnL);
  void (*lnNorma2_bs) (NcmDataGaussCov *gauss, NcmMSet *mset, NcmBootstrap *bstrap, gdouble *m2lnL);
  void (*set_size) (NcmDataGaussCov *gauss, guint np);
  guint (*get_size) (NcmDataGaussCov *gauss);

  /* Padding to allow 18 virtual functions without breaking ABI. */
  gpointer padding[12];
};

void ncm_data_gauss_cov_set_size (NcmDataGaussCov *gauss, guint np);
guint ncm_data_gauss_cov_get_size (NcmDataGaussCov *gauss);

void ncm_data_gauss_cov_use_norma (NcmDataGaussCov *gauss, gboolean use_norma);

void ncm_data_gauss_cov_replace_mean (NcmDataGaussCov *gauss, NcmVector *mean);
NcmVector *ncm_data_gauss_cov_peek_mean (NcmDataGaussCov *gauss);
NcmMatrix *ncm_data_gauss_cov_peek_cov (NcmDataGaussCov *gauss);

gdouble ncm_data_gauss_cov_get_log_norma (NcmDataGaussCov *gauss, NcmMSet *mset);

G_END_DECLS

#endif /* _NCM_DATA_GAUSS_COV_H_ */


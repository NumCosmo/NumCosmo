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

#define NCM_TYPE_DATA_GAUSS_COV             (ncm_data_gauss_cov_get_type ())
#define NCM_DATA_GAUSS_COV(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_DATA_GAUSS_COV, NcmDataGaussCov))
#define NCM_DATA_GAUSS_COV_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_DATA_GAUSS_COV, NcmDataGaussCovClass))
#define NCM_IS_DATA_GAUSS_COV(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_DATA_GAUSS_COV))
#define NCM_IS_DATA_GAUSS_COV_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_DATA_GAUSS_COV))
#define NCM_DATA_GAUSS_COV_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_DATA_GAUSS_COV, NcmDataGaussCovClass))

typedef struct _NcmDataGaussCovClass NcmDataGaussCovClass;
typedef struct _NcmDataGaussCov NcmDataGaussCov;

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
};

struct _NcmDataGaussCov
{
  /*< private >*/
  NcmData parent_instance;
  guint np;
  NcmVector *y;
  NcmVector *v;
  NcmMatrix *cov;
  NcmMatrix *LLT;
  gboolean prepared_LLT;
  gboolean use_norma;
};

GType ncm_data_gauss_cov_get_type (void) G_GNUC_CONST;

void ncm_data_gauss_cov_set_size (NcmDataGaussCov *gauss, guint np);
guint ncm_data_gauss_cov_get_size (NcmDataGaussCov *gauss);

void ncm_data_gauss_cov_use_norma (NcmDataGaussCov *gauss, gboolean use_norma);

G_END_DECLS

#endif /* _NCM_DATA_GAUSS_COV_H_ */


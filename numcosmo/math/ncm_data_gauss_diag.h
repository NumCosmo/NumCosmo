/***************************************************************************
 *            ncm_data_gauss_diag.h
 *
 *  Fri Mar 19 14:57:54 2010
 *  Copyright  2010  Sandro Dias Pinto Vitenti
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

#ifndef _NCM_DATA_GAUSS_DIAG_H_
#define _NCM_DATA_GAUSS_DIAG_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_data.h>

G_BEGIN_DECLS

#define NCM_TYPE_DATA_GAUSS_DIAG             (ncm_data_gauss_diag_get_type ())
#define NCM_DATA_GAUSS_DIAG(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_DATA_GAUSS_DIAG, NcmDataGaussDiag))
#define NCM_DATA_GAUSS_DIAG_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_DATA_GAUSS_DIAG, NcmDataGaussDiagClass))
#define NCM_IS_DATA_GAUSS_DIAG(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_DATA_GAUSS_DIAG))
#define NCM_IS_DATA_GAUSS_DIAG_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_DATA_GAUSS_DIAG))
#define NCM_DATA_GAUSS_DIAG_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_DATA_GAUSS_DIAG, NcmDataGaussDiagClass))

typedef struct _NcmDataGaussDiagClass NcmDataGaussDiagClass;
typedef struct _NcmDataGaussDiag NcmDataGaussDiag;

struct _NcmDataGaussDiagClass
{
  /*< private >*/
  NcmDataClass parent_class;

  void (*mean_func) (NcmDataGaussDiag *diag, NcmMSet *mset, NcmVector *vp);
  gboolean (*sigma_func) (NcmDataGaussDiag *diag, NcmMSet *mset, NcmVector *var);
  void (*set_size) (NcmDataGaussDiag *diag, guint np);
  guint (*get_size) (NcmDataGaussDiag *diag);
};

struct _NcmDataGaussDiag
{
  /*< private >*/
  NcmData parent_instance;
  guint np;
  NcmVector *y;
  NcmVector *v;
  NcmVector *sigma;
  NcmVector *weight;
  gdouble wt;
  gboolean prepared_w;
  gboolean wmean;
};

GType ncm_data_gauss_diag_get_type (void) G_GNUC_CONST;

void ncm_data_gauss_diag_set_size (NcmDataGaussDiag *diag, guint np);
guint ncm_data_gauss_diag_get_size (NcmDataGaussDiag *diag);

G_END_DECLS

#endif /* _NCM_DATA_GAUSS_DIAG_H_ */


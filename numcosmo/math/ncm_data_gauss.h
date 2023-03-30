/***************************************************************************
 *            ncm_data_gauss.h
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

#ifndef _NCM_DATA_GAUSS_H_
#define _NCM_DATA_GAUSS_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_data.h>

G_BEGIN_DECLS

#define NCM_TYPE_DATA_GAUSS             (ncm_data_gauss_get_type ())

G_DECLARE_DERIVABLE_TYPE (NcmDataGauss, ncm_data_gauss, NCM, DATA_GAUSS, NcmData)

struct _NcmDataGaussClass
{
  /*< private >*/
  NcmDataClass parent_class;

  void (*mean_func) (NcmDataGauss *gauss, NcmMSet *mset, NcmVector *vp);
  gboolean (*inv_cov_func) (NcmDataGauss *gauss, NcmMSet *mset, NcmMatrix *inv_cov);
  void (*set_size) (NcmDataGauss *gauss, guint np);
  guint (*get_size) (NcmDataGauss *gauss);

  /* Padding to allow 18 virtual functions without breaking ABI. */
  gpointer padding[14];
};

void ncm_data_gauss_set_size (NcmDataGauss *gauss, guint np);
guint ncm_data_gauss_get_size (NcmDataGauss *gauss);

NcmMatrix *ncm_data_gauss_peek_inv_cov (NcmDataGauss *gauss);
NcmVector *ncm_data_gauss_peek_mean (NcmDataGauss *gauss);

G_END_DECLS

#endif /* _NCM_DATA_GAUSS_H_ */


/***************************************************************************
 *            nc_reduced_shear_calib.h
 *
 *  Tue December 04 23:40:36 2018
 *  Copyright  2018  Mariana Penna Lima
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * nc_reduced_shear_calib.h
 * Copyright (C) 2018 Mariana Penna Lima <pennalima@gmail.com>
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

#ifndef _NC_REDUCED_SHEAR_CALIB_H_
#define _NC_REDUCED_SHEAR_CALIB_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/ncm/model/ncm_mset.h>
#include <numcosmo/ncm/model/ncm_model.h>

G_BEGIN_DECLS

#define NC_TYPE_REDUCED_SHEAR_CALIB             (nc_reduced_shear_calib_get_type ())

G_DECLARE_DERIVABLE_TYPE (NcReducedShearCalib, nc_reduced_shear_calib, NC, REDUCED_SHEAR_CALIB, NcmModel)

struct _NcReducedShearCalibClass
{
  /*< private >*/
  NcmModelClass parent_class;
  gdouble (*eval) (NcReducedShearCalib *rs_calib, const gdouble g_th, const gdouble psf_size, const gdouble gal_size);

  /* Padding to allow adding up to 17 more virtual functions without breaking ABI. */
  gpointer padding[17];
};


NCM_MSET_MODEL_DECLARE_ID (nc_reduced_shear_calib);

NcReducedShearCalib *nc_reduced_shear_calib_ref (NcReducedShearCalib *rs_calib);
void nc_reduced_shear_calib_free (NcReducedShearCalib *rs_calib);
void nc_reduced_shear_calib_clear (NcReducedShearCalib **rs_calib);

gdouble nc_reduced_shear_calib_eval (NcReducedShearCalib *rs_calib, const gdouble g_th, const gdouble psf_size, const gdouble gal_size);

G_END_DECLS

#endif /* _NC_REDUCED_SHEAR_CALIB_H_ */


/* -*- Mode: C; indent-tabs-mode: nil; c-basic-offset: 2; tab-width: 2 -*-  */
/***************************************************************************
 *            nc_reduced_shear_calib_wtg.h
 *
 *  Tue December 04 23:47:36 2018
 *  Copyright  2018  Mariana Penna Lima
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * nc_reduced_shear_calib_wtg.h
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

#ifndef _NC_REDUCED_SHEAR_CALIB_WTG_H_
#define _NC_REDUCED_SHEAR_CALIB_WTG_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_mset.h>
#include <numcosmo/math/ncm_model.h>
#include <numcosmo/lss/nc_reduced_shear_calib.h>

G_BEGIN_DECLS

#define NC_TYPE_REDUCED_SHEAR_CALIB_WTG             (nc_reduced_shear_calib_wtg_get_type ())
#define NC_REDUCED_SHEAR_CALIB_WTG(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_REDUCED_SHEAR_CALIB_WTG, NcReducedShearCalibWtg))
#define NC_REDUCED_SHEAR_CALIB_WTG_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_REDUCED_SHEAR_CALIB_WTG, NcReducedShearCalibWtgClass))
#define NC_IS_REDUCED_SHEAR_CALIB_WTG(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_REDUCED_SHEAR_CALIB_WTG))
#define NC_IS_REDUCED_SHEAR_CALIB_WTG_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_REDUCED_SHEAR_CALIB_WTG))
#define NC_REDUCED_SHEAR_CALIB_WTG_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_REDUCED_SHEAR_CALIB_WTG, NcReducedShearCalibWtgClass))

typedef struct _NcReducedShearCalibWtgClass NcReducedShearCalibWtgClass;
typedef struct _NcReducedShearCalibWtg NcReducedShearCalibWtg;
typedef struct _NcReducedShearCalibWtgPrivate NcReducedShearCalibWtgPrivate;

struct _NcReducedShearCalibWtgClass
{
  /*< private >*/
  NcReducedShearCalibClass parent_class;
};

struct _NcReducedShearCalibWtg
{
  /*< private >*/
  NcReducedShearCalib parent_instance;
  NcReducedShearCalibWtgPrivate *priv;
};

/**
 * NcReducedShearCalibWtgSParams:
 * @NC_REDUCED_SHEAR_CALIB_WTG_MSLOPE: FIXME
 * @NC_REDUCED_SHEAR_CALIB_WTG_MB: FIXME
 * @NC_REDUCED_SHEAR_CALIB_WTG_C: FIXME
 * @NC_REDUCED_SHEAR_CALIB_WTG_SIZE_RATIO: ratio between galaxy and psf sizes 
 *
 * WtG calibration parameters.
 * 
 */
typedef enum _NcReducedShearCalibWtgSParams
{
  NC_REDUCED_SHEAR_CALIB_WTG_MSLOPE, 
  NC_REDUCED_SHEAR_CALIB_WTG_MB, 
  NC_REDUCED_SHEAR_CALIB_WTG_C, 
  NC_REDUCED_SHEAR_CALIB_WTG_SIZE_RATIO, 
  /* < private > */
  NNC_REDUCED_SHEAR_CALIB_WTG_SPARAM_LEN, /*< skip >*/
} NcReducedShearCalibWtgSParams;

GType nc_reduced_shear_calib_wtg_get_type (void) G_GNUC_CONST;

NcReducedShearCalibWtg *nc_reduced_shear_calib_wtg_new (void);
NcReducedShearCalibWtg *nc_reduced_shear_calib_wtg_ref (NcReducedShearCalibWtg *rs_wtg);
void nc_reduced_shear_calib_wtg_free (NcReducedShearCalibWtg *rs_wtg);
void nc_reduced_shear_calib_wtg_clear (NcReducedShearCalibWtg **rs_wtg);

G_END_DECLS

#endif /* _NC_REDUCED_SHEAR_CALIB_WTG_H_ */

/***************************************************************************
 *            nc_powspec_ml_fix_spline.h
 *
 *  Sun Jun 18 10:26:38 2017
 *  Copyright  2017  Mariana Penna Lima
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * nc_powspec_ml_fix_spline.h
 * Copyright (C) 2017 Mariana Penna Lima <pennalima@gmail.com>
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

#ifndef _NC_POWSPEC_ML_FIX_SPLINE_H_
#define _NC_POWSPEC_ML_FIX_SPLINE_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_spline.h>
#include <numcosmo/lss/nc_growth_func.h>
#include <numcosmo/nc_powspec_ml.h>

G_BEGIN_DECLS

#define NC_TYPE_POWSPEC_ML_FIX_SPLINE             (nc_powspec_ml_fix_spline_get_type ())
#define NC_POWSPEC_ML_FIX_SPLINE(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_POWSPEC_ML_FIX_SPLINE, NcPowspecMLFixSpline))
#define NC_POWSPEC_ML_FIX_SPLINE_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_POWSPEC_ML_FIX_SPLINE, NcPowspecMLFixSplineClass))
#define NC_IS_POWSPEC_ML_FIX_SPLINE(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_POWSPEC_ML_FIX_SPLINE))
#define NC_IS_POWSPEC_ML_FIX_SPLINE_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_POWSPEC_ML_FIX_SPLINE))
#define NC_POWSPEC_ML_FIX_SPLINE_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_POWSPEC_ML_FIX_SPLINE, NcPowspecMLFixSplineClass))

typedef struct _NcPowspecMLFixSplineClass NcPowspecMLFixSplineClass;
typedef struct _NcPowspecMLFixSpline NcPowspecMLFixSpline;

struct _NcPowspecMLFixSplineClass
{
  /*< private > */
  NcPowspecMLClass parent_class;
};

struct _NcPowspecMLFixSpline
{
  /*< private > */
  NcPowspecML parent_instance;
  NcmSerialize *ser;
  NcmSpline *Pk;
  NcGrowthFunc *gf;
  gchar *filename;
};

GType nc_powspec_ml_fix_spline_get_type (void) G_GNUC_CONST;

NcPowspecMLFixSpline *nc_powspec_ml_fix_spline_new (const gchar *filename);

void nc_powspec_ml_fix_spline_set_file (NcPowspecMLFixSpline *ps_fs, const gchar *filename);

G_END_DECLS

#endif /* _NC_POWSPEC_ML_FIX_SPLINE_H_ */


/***************************************************************************
 *            nc_powspec_ml_spline.h
 *
 *  Sun Jun 18 10:26:38 2017
 *  Copyright  2017  Mariana Penna Lima
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * nc_powspec_ml_spline.h
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

#ifndef _NC_POWSPEC_ML_SPLINE_H_
#define _NC_POWSPEC_ML_SPLINE_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_spline.h>
#include <numcosmo/lss/nc_growth_func.h>
#include <numcosmo/nc_powspec_ml.h>

G_BEGIN_DECLS

#define NC_TYPE_POWSPEC_ML_SPLINE (nc_powspec_ml_spline_get_type ())

G_DECLARE_FINAL_TYPE (NcPowspecMLSpline, nc_powspec_ml_spline, NC, POWSPEC_ML_SPLINE, NcPowspecML)

NcPowspecMLSpline *nc_powspec_ml_spline_new (const gchar *filename);

void nc_powspec_ml_spline_set_spline (NcPowspecMLSpline *ps_fs, NcmSpline *Pk);
NcmSpline *nc_powspec_ml_spline_peek_spline (NcPowspecMLSpline *ps_fs);

void nc_powspec_ml_spline_load_from_file (NcPowspecMLSpline *ps_fs, const gchar *filename);

G_END_DECLS

#endif /* _NC_POWSPEC_ML_SPLINE_H_ */


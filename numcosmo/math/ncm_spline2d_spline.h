/***************************************************************************
 *            ncm_spline2d_spline.h
 *
 *  Sun Aug  1 17:17:20 2010
 *  Copyright  2010  Mariana Penna Lima & Sandro Dias Pinto Vitenti
 *  <pennalima@gmail.com>, <vitenti@uel.br>
 ****************************************************************************/

/*
 * numcosmo
 * Copyright (C) Mariana Penna Lima & Sandro Dias Pinto Vitenti 2012 <pennalima@gmail.com>, <vitenti@uel.br>
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

#ifndef _NCM_SPLINE2D_SPLINE_H_
#define _NCM_SPLINE2D_SPLINE_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_spline2d.h>

G_BEGIN_DECLS

#define NCM_TYPE_SPLINE2D_SPLINE (ncm_spline2d_spline_get_type ())

G_DECLARE_FINAL_TYPE (NcmSpline2dSpline, ncm_spline2d_spline, NCM, SPLINE2D_SPLINE, NcmSpline2d)

NcmSpline2d *ncm_spline2d_spline_new (NcmSpline * s);

G_END_DECLS

#endif /* _NCM_SPLINE2D_SPLINE_H_ */


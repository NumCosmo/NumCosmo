/***************************************************************************
 *            ncm_spline_cubic_d2.h
 *
 *  Wed Aug 13 21:13:59 2008
 *  Copyright  2008  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/

/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <vitenti@uel.br>
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

#ifndef _NCM_SPLINE_CUBIC_D2_H_
#define _NCM_SPLINE_CUBIC_D2_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_spline_cubic.h>

G_BEGIN_DECLS

#define NCM_TYPE_SPLINE_CUBIC_D2 (ncm_spline_cubic_d2_get_type ())

G_DECLARE_FINAL_TYPE (NcmSplineCubicD2, ncm_spline_cubic_d2, NCM, SPLINE_CUBIC_D2, NcmSplineCubic)

NcmSplineCubicD2 *ncm_spline_cubic_d2_new (NcmVector * xv, NcmVector * yv, NcmVector * d2yv, gboolean init);

G_END_DECLS

#endif /* _NCM_SPLINE_CUBIC_D2_H_ */


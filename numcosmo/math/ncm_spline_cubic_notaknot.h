/***************************************************************************
 *            ncm_spline_cubic_notaknot.h
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

#ifndef _NCM_SPLINE_CUBIC_NOTAKNOT_H_
#define _NCM_SPLINE_CUBIC_NOTAKNOT_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_spline_cubic.h>

G_BEGIN_DECLS

#define NCM_TYPE_SPLINE_CUBIC_NOTAKNOT (ncm_spline_cubic_notaknot_get_type ())

G_DECLARE_FINAL_TYPE (NcmSplineCubicNotaknot, ncm_spline_cubic_notaknot, NCM, SPLINE_CUBIC_NOTAKNOT, NcmSplineCubic)

NcmSpline *ncm_spline_cubic_notaknot_new (void);
NcmSpline *ncm_spline_cubic_notaknot_new_full (NcmVector *xv, NcmVector *yv, gboolean init);

G_END_DECLS

#endif /* _NCM_SPLINE_CUBIC_NOTAKNOT_H_ */


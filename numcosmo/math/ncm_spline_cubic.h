/***************************************************************************
 *            ncm_spline_cubic.h
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

#ifndef _NCM_SPLINE_CUBIC_H_
#define _NCM_SPLINE_CUBIC_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_spline.h>

G_BEGIN_DECLS

#define NCM_TYPE_SPLINE_CUBIC (ncm_spline_cubic_get_type ())

G_DECLARE_DERIVABLE_TYPE (NcmSplineCubic, ncm_spline_cubic, NCM, SPLINE_CUBIC, NcmSpline)

struct _NcmSplineCubicClass
{
  /*< private >*/
  NcmSplineClass parent_class;

  /* Padding to allow 18 virtual functions without breaking ABI. */
  gpointer padding[18];
};

NcmVector *ncm_spline_cubic_peek_b_vec (const NcmSplineCubic *s);
NcmVector *ncm_spline_cubic_peek_c_vec (const NcmSplineCubic *s);
NcmVector *ncm_spline_cubic_peek_d_vec (const NcmSplineCubic *s);

NcmVector *ncm_spline_cubic_peek_diag_vec (const NcmSplineCubic *s);
NcmVector *ncm_spline_cubic_peek_offdiag_vec (const NcmSplineCubic *s);

G_END_DECLS

#endif /* _NCM_SPLINE_CUBIC_H_ */


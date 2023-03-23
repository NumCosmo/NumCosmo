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

#define NCM_TYPE_SPLINE_CUBIC_D2             (ncm_spline_cubic_d2_get_type ())
#define NCM_SPLINE_CUBIC_D2(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_SPLINE_CUBIC_D2, NcmSplineCubicD2))
#define NCM_SPLINE_CUBIC_D2_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_SPLINE_CUBIC_D2, NcmSplineCubicD2Class))
#define NCM_IS_SPLINE_CUBIC_D2(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_SPLINE_CUBIC_D2))
#define NCM_IS_SPLINE_CUBIC_D2_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_SPLINE_CUBIC_D2))
#define NCM_SPLINE_CUBIC_D2_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_SPLINE_CUBIC_D2, NcmSplineCubicD2Class))

typedef struct _NcmSplineCubicD2Class NcmSplineCubicD2Class;
typedef struct _NcmSplineCubicD2 NcmSplineCubicD2;

struct _NcmSplineCubicD2Class
{
  /*< private >*/
  NcmSplineCubicClass parent_class;
};

struct _NcmSplineCubicD2
{
  /*< private >*/
  NcmSplineCubic parent_instance;
  NcmVector *d2;
  gsize len;
};

GType ncm_spline_cubic_d2_get_type (void) G_GNUC_CONST;

NcmSplineCubicD2 *ncm_spline_cubic_d2_new (NcmVector *xv, NcmVector *yv, NcmVector *d2yv, gboolean init);

G_END_DECLS

#endif /* _NCM_SPLINE_CUBIC_D2_H_ */


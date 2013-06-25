/***************************************************************************
 *            ncm_spline2d_spline.h
 *
 *  Sun Aug  1 17:17:20 2010
 *  Copyright  2010  Mariana Penna Lima & Sandro Dias Pinto Vitenti
 *  <pennalima@gmail.com>, <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Mariana Penna Lima & Sandro Dias Pinto Vitenti 2012 <pennalima@gmail.com>, <sandro@isoftware.com.br>
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

#define NCM_TYPE_SPLINE2D_SPLINE             (ncm_spline2d_spline_get_type ())
#define NCM_SPLINE2D_SPLINE(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_SPLINE2D_SPLINE, NcmSpline2dSpline))
#define NCM_SPLINE2D_SPLINE_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_SPLINE2D_SPLINE, NcmSpline2dSplineClass))
#define NCM_IS_SPLINE2D_SPLINE(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_SPLINE2D_SPLINE))
#define NCM_IS_SPLINE2D_SPLINE_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_SPLINE2D_SPLINE))
#define NCM_SPLINE2D_SPLINE_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_SPLINE2D_SPLINE, NcmSpline2dSplineClass))

typedef struct _NcmSpline2dSplineClass NcmSpline2dSplineClass;
typedef struct _NcmSpline2dSpline NcmSpline2dSpline;

struct _NcmSpline2dSplineClass
{
	/*< private >*/
	NcmSpline2dClass parent_class;
};

struct _NcmSpline2dSpline
{
	/*< private >*/
	NcmSpline2d parent_instance;
	gboolean first_prepare;
	gboolean first_prepare_integ;
	gdouble last_x;
	gdouble last_xl;
	gdouble last_xu;
	gdouble last_yl;
	gdouble last_yu;
	NcmVector *vertv;
	NcmVector *vertintv;
	NcmSpline **s_hor;
	NcmSpline *s_ver;
	NcmSpline *s_ver_integ;
  guint s_hor_len;
};

GType ncm_spline2d_spline_get_type (void) G_GNUC_CONST;

NcmSpline2d *ncm_spline2d_spline_new (NcmSpline *s);

G_END_DECLS

#endif /* _NCM_SPLINE2D_SPLINE_H_ */

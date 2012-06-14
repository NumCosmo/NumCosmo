/***************************************************************************
 *            ncm_spline_cubic.h
 *
 *  Wed Aug 13 21:13:59 2008
 *  Copyright  2008  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <sandro@isoftware.com.br>
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

#include <glib-object.h>

G_BEGIN_DECLS

#define NCM_TYPE_SPLINE_CUBIC             (ncm_spline_cubic_get_type ())
#define NCM_SPLINE_CUBIC(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_SPLINE_CUBIC, NcmSplineCubic))
#define NCM_SPLINE_CUBIC_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_SPLINE_CUBIC, NcmSplineCubicClass))
#define NCM_IS_SPLINE_CUBIC(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_SPLINE_CUBIC))
#define NCM_IS_SPLINE_CUBIC_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_SPLINE_CUBIC))
#define NCM_SPLINE_CUBIC_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_SPLINE_CUBIC, NcmSplineCubicClass))

typedef struct _NcmSplineCubicClass NcmSplineCubicClass;
typedef struct _NcmSplineCubic NcmSplineCubic;

struct _NcmSplineCubicClass
{
  /*< private >*/
	NcmSplineClass parent_class;
};

struct _NcmSplineCubic
{
  /*< private >*/
  NcmSpline parent_instance;
  NcmVector *b;
  NcmVector *c;
  NcmVector *d;
  NcmVector *g;
  NcmVector *diag;
  NcmVector *offdiag;
  gboolean init;
  gsize len;
};

GType ncm_spline_cubic_get_type (void) G_GNUC_CONST;

G_END_DECLS

#endif /* _NCM_SPLINE_CUBIC_H_ */

/***************************************************************************
 *            ncm_spline_cubic_notaknot.h
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

#ifndef _NCM_SPLINE_CUBIC_NOTAKNOT_H_
#define _NCM_SPLINE_CUBIC_NOTAKNOT_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_spline_cubic.h>

G_BEGIN_DECLS

#define NCM_TYPE_SPLINE_CUBIC_NOTAKNOT             (ncm_spline_cubic_notaknot_get_type ())
#define NCM_SPLINE_CUBIC_NOTAKNOT(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_SPLINE_CUBIC_NOTAKNOT, NcmSplineCubicNotaknot))
#define NCM_SPLINE_CUBIC_NOTAKNOT_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_SPLINE_CUBIC_NOTAKNOT, NcmSplineCubicNotaknotClass))
#define NCM_IS_SPLINE_CUBIC_NOTAKNOT(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_SPLINE_CUBIC_NOTAKNOT))
#define NCM_IS_SPLINE_CUBIC_NOTAKNOT_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_SPLINE_CUBIC_NOTAKNOT))
#define NCM_SPLINE_CUBIC_NOTAKNOT_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_SPLINE_CUBIC_NOTAKNOT, NcmSplineCubicNotaknotClass))

typedef struct _NcmSplineCubicNotaknotClass NcmSplineCubicNotaknotClass;
typedef struct _NcmSplineCubicNotaknot NcmSplineCubicNotaknot;

struct _NcmSplineCubicNotaknotClass
{
  /*< private >*/
  NcmSplineCubicClass parent_class;
};

struct _NcmSplineCubicNotaknot
{
  /*< private >*/
  NcmSplineCubic parent_instance;
};

GType ncm_spline_cubic_notaknot_get_type (void) G_GNUC_CONST;

NcmSpline *ncm_spline_cubic_notaknot_new (void);
NcmSpline *ncm_spline_cubic_notaknot_new_full (NcmVector *xv, NcmVector *yv, gboolean init);

G_END_DECLS

#endif /* _NCM_SPLINE_CUBIC_NOTAKNOT_H_ */

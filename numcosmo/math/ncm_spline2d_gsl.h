/***************************************************************************
 *            ncm_spline2d_gsl.h
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

#ifndef _NCM_SPLINE2D_GSL_H_
#define _NCM_SPLINE2D_GSL_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_spline2d.h>

G_BEGIN_DECLS

#define NCM_TYPE_SPLINE2D_GSL             (ncm_spline2d_gsl_get_type ())
#define NCM_SPLINE2D_GSL(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_SPLINE2D_GSL, NcmSpline2dGsl))
#define NCM_SPLINE2D_GSL_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_SPLINE2D_GSL, NcmSpline2dGslClass))
#define NCM_IS_SPLINE2D_GSL(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_SPLINE2D_GSL))
#define NCM_IS_SPLINE2D_GSL_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_SPLINE2D_GSL))
#define NCM_SPLINE2D_GSL_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_SPLINE2D_GSL, NcmSpline2dGslClass))

typedef struct _NcmSpline2dGslClass NcmSpline2dGslClass;
typedef struct _NcmSpline2dGsl NcmSpline2dGsl;

struct _NcmSpline2dGsl
{
  /*< private >*/
  NcmSpline2d parent_instance;
  NcmMatrix *zdiff;
  NcmVector *vertv;
  NcmVector *vertintv;
  NcmSpline **s_hor;
  NcmSpline **s_dzdy;
  NcmSpline *s_ver;
  NcmSpline *s_ver_integ;
  guint s_hor_len;
};

struct _NcmSpline2dGslClass
{
  /*< private >*/
  NcmSpline2dClass parent_class;
};

GType ncm_spline2d_gsl_get_type (void) G_GNUC_CONST;

NcmSpline2d *ncm_spline2d_gsl_new (NcmSpline *s);
NcmSpline2d *ncm_spline2d_gsl_natural_new (void);

G_END_DECLS

#endif /* _NCM_SPLINE2D_GSL_H_ */

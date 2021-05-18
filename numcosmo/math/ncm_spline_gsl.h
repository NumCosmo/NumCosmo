/***************************************************************************
 *            ncm_spline_gsl.h
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

#ifndef _NCM_SPLINE_GSL_H_
#define _NCM_SPLINE_GSL_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_spline.h>

G_BEGIN_DECLS

#define NCM_TYPE_SPLINE_GSL             (ncm_spline_gsl_get_type ())
#define NCM_SPLINE_GSL(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_SPLINE_GSL, NcmSplineGsl))
#define NCM_SPLINE_GSL_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_SPLINE_GSL, NcmSplineGslClass))
#define NCM_IS_SPLINE_GSL(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_SPLINE_GSL))
#define NCM_IS_SPLINE_GSL_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_SPLINE_GSL))
#define NCM_SPLINE_GSL_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_SPLINE_GSL, NcmSplineGslClass))

typedef struct _NcmSplineGslClass NcmSplineGslClass;
typedef struct _NcmSplineGsl NcmSplineGsl;

struct _NcmSplineGslClass
{
  /*< private >*/
  NcmSplineClass parent_class;
};

/**
 * NcmSplineGslType:
 * @NCM_SPLINE_GSL_LINEAR: Uses [gsl_interp_linear](https://www.gnu.org/software/gsl/doc/html/interp.html#c.gsl_interp_linear) interpolation method. 
 * @NCM_SPLINE_GSL_POLYNOMIAL: Uses [gsl_interp_polynomial](https://www.gnu.org/software/gsl/doc/html/interp.html#c.gsl_interp_polynomial) interpolation method. 
 * @NCM_SPLINE_GSL_CSPLINE: Uses [gsl_interp_cspline](https://www.gnu.org/software/gsl/doc/html/interp.html#c.gsl_interp_cspline) interpolation method. 
 * @NCM_SPLINE_GSL_CSPLINE_PERIODIC: Uses [gsl_interp_cspline_periodic](https://www.gnu.org/software/gsl/doc/html/interp.html#c.gsl_interp_cspline_periodic) interpolation method.
 * @NCM_SPLINE_GSL_AKIMA: Uses [gsl_interp_akima](https://www.gnu.org/software/gsl/doc/html/interp.html#c.gsl_interp_akima) interpolation method. 
 * @NCM_SPLINE_GSL_AKIMA_PERIODIC: Uses [gsl_interp_akima_periodic](https://www.gnu.org/software/gsl/doc/html/interp.html#c.gsl_interp_akima_periodic) interpolation method. 
 *
 * Enumeration to choose which [GSL interpolation method](https://www.gnu.org/software/gsl/doc/html/interp.html#interpolation) as backend to be used by the object. It can be used with the function ncm_spline_gsl_new_by_id() when a new #NcmSplineGsl is created.  
 */
typedef enum _NcmSplineGslType
{
  NCM_SPLINE_GSL_LINEAR = 0,
  NCM_SPLINE_GSL_POLYNOMIAL,
  NCM_SPLINE_GSL_CSPLINE,
  NCM_SPLINE_GSL_CSPLINE_PERIODIC,
  NCM_SPLINE_GSL_AKIMA,
  NCM_SPLINE_GSL_AKIMA_PERIODIC,
  /* < private > */
  NCM_SPLINE_GSL_TYPES_LEN, /*< skip >*/
} NcmSplineGslType;

struct _NcmSplineGsl
{
  /*< private >*/
  NcmSpline parent_instance;
  gsl_interp *interp;
  NcmSplineGslType type_id;
  gchar *inst_name;
  const gsl_interp_type *type;
};

GType ncm_spline_gsl_get_type (void) G_GNUC_CONST;

NcmSpline *ncm_spline_gsl_new (const gsl_interp_type *type);
NcmSpline *ncm_spline_gsl_new_full (const gsl_interp_type *type, NcmVector *xv, NcmVector *yv, gboolean init);
NcmSpline *ncm_spline_gsl_new_by_id (NcmSplineGslType type_id);
NcmSpline *ncm_spline_gsl_new_full_by_id (NcmSplineGslType type_id, NcmVector *xv, NcmVector *yv, gboolean init);
void ncm_spline_gsl_set_type (NcmSplineGsl *sg, const gsl_interp_type *type);
void ncm_spline_gsl_set_type_by_id (NcmSplineGsl *sg, NcmSplineGslType type_id);
void ncm_spline_gsl_set_type_by_name (NcmSplineGsl *sg, const gchar *type_name);

G_END_DECLS

#endif /* _NCM_SPLINE_GSL_H_ */

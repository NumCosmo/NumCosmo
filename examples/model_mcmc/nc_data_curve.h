/***************************************************************************
 *            nc_data_curve.h
 *
 *  Wed July 12 13:46:10 2017
 *  Copyright  2017  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_data_curve.h
 * Copyright (C) 2017 Sandro Dias Pinto Vitenti <vitenti@uel.br>
 *
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

#ifndef _NC_DATA_CURVE_H_
#define _NC_DATA_CURVE_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>

/*
 * Here we are writing an object external to NumCosmo so we must use
 * the full header numcosmo/numcosmo.h. Otherwise we could include
 * only the actual header we need, i.e.,
 * #include <numcosmo/math/ncm_data_gauss_diag.h>
 * 
 */
#include <numcosmo/numcosmo.h>

G_BEGIN_DECLS

/* These are the basic macros useful for the GObject framework */
#define NC_TYPE_DATA_CURVE             (nc_data_curve_get_type ())
#define NC_DATA_CURVE(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_DATA_CURVE, NcDataCurve))
#define NC_DATA_CURVE_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_DATA_CURVE, NcDataCurveClass))
#define NC_IS_DATA_CURVE(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_DATA_CURVE))
#define NC_IS_DATA_CURVE_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_DATA_CURVE))
#define NC_DATA_CURVE_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_DATA_CURVE, NcDataCurveClass))

/* 
 * This is the class struct, it will contains everything that is
 * that is defined across instances.
 * 
 */
typedef struct _NcDataCurveClass NcDataCurveClass;

/* 
 * This is the instance struct, it will contains everything pertaining 
 * to a given instance.
 * 
 */
typedef struct _NcDataCurve NcDataCurve;

/*
 * The first item is just the parent object structure.
 * 
 */
typedef struct _NcDataCurvePrivate NcDataCurvePrivate;

/*
 * The first item is just the parent class structure.
 * 
 */
struct _NcDataCurveClass
{
  /*< private >*/
  NcmDataGaussDiagClass parent_class;
};

/*
 * No elements needed apart from our private struct.
 * The first item is just the parent object structure.
 * 
 */
struct _NcDataCurve
{
  /*< private >*/
  NcmDataGaussDiag parent_instance;
  NcDataCurvePrivate *priv;
};

GType nc_data_curve_get_type (void) G_GNUC_CONST;

/* 
 * The default methods that all objects should implement.
 * 
 */
NcDataCurve *nc_data_curve_new (const guint np);
NcDataCurve *nc_data_curve_ref (NcDataCurve *data_curve);
void nc_data_curve_free (NcDataCurve *data_curve);
void nc_data_curve_clear (NcDataCurve **data_curve);

/* 
 * Here we create a simple method that can be used to initialize 
 * our data object.
 * 
 */
void nc_data_curve_init_data (NcDataCurve *data_curve, NcmVector *xv, NcmVector *yv, NcmVector *sigmav);

G_END_DECLS

#endif /* _NC_DATA_CURVE_H_ */


/***************************************************************************
 *            nc_curve_linear.h
 *
 *  Wed July 12 10:03:02 2017
 *  Copyright  2017  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * nc_curve_linear.h
 * Copyright (C) 2017 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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

#ifndef _NC_CURVE_LINEAR_H_
#define _NC_CURVE_LINEAR_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>

/*
 * We need the parents head to define its child.
 * 
 */
#include "nc_curve.h"

G_BEGIN_DECLS

/* These are the basic macros useful for the GObject framework */
#define NC_TYPE_CURVE_LINEAR             (nc_curve_linear_get_type ())
#define NC_CURVE_LINEAR(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_CURVE_LINEAR, NcCurveLinear))
#define NC_CURVE_LINEAR_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_CURVE_LINEAR, NcCurveLinearClass))
#define NC_IS_CURVE_LINEAR(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_CURVE_LINEAR))
#define NC_IS_CURVE_LINEAR_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_CURVE_LINEAR))
#define NC_CURVE_LINEAR_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_CURVE_LINEAR, NcCurveLinearClass))

/* 
 * This is the class struct, it will contains everything that is
 * that is defined across instances.
 * 
 */
typedef struct _NcCurveLinearClass NcCurveLinearClass;

/* 
 * This is the instance struct, it will contains everything pertaining 
 * to a given instance.
 * 
 */
typedef struct _NcCurveLinear NcCurveLinear;

/*
 * The first item is just the parent object structure.
 * 
 */
struct _NcCurveLinearClass
{
  /*< private >*/
  NcCurveClass parent_class;
};

/*
 * Since this object has parameters we need to define its
 * parameters enumerator.
 * 
 */
/**
 * NcCurveLinearParams:
 * @NC_CURVE_LINEAR_A: $a$
 * @NC_CURVE_LINEAR_B: $b$
 *
 * The two necessary parameters for our 
 * linear model $f(x) = a + b x$.
 * 
 */
typedef enum _NcCurveLinearParams
{
  NC_CURVE_LINEAR_A = 0,
  NC_CURVE_LINEAR_B,          /*< private >*/
  NC_CURVE_LINEAR_SPARAM_LEN, /*< skip >*/
} NcCurveLinearParams;

/*
 * The first item is just the parent object structure.
 * 
 */
struct _NcCurveLinear
{
  /*< private >*/
  NcCurve parent_instance;
};

GType nc_curve_linear_get_type (void) G_GNUC_CONST;

/* 
 * The default methods that all objects should implement.
 * 
 */
NcCurveLinear *nc_curve_linear_new (const gdouble xl, const gdouble xu);
NcCurveLinear *nc_curve_linear_ref (NcCurveLinear *curve_linear);
void nc_curve_linear_free (NcCurveLinear *curve_linear);
void nc_curve_linear_clear (NcCurveLinear **curve_linear);

/*
 * In this simple case we don't need any additional methods.
 * 
 */

G_END_DECLS

#endif /* _NC_CURVE_LINEAR_H_ */

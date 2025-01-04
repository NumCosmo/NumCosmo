/***************************************************************************
 *            nc_curve.h
 *
 *  Tue July 11 16:30:36 2017
 *  Copyright  2017  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_curve.h
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

#ifndef _NC_CURVE_H_
#define _NC_CURVE_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>

/*
 * Here we are writing an object external to NumCosmo so we must use
 * the full header numcosmo/numcosmo.h. Otherwise we could include
 * only the actual header we need, i.e.,
 * #include <numcosmo/math/ncm_model.h>
 *
 */
#include <numcosmo/numcosmo.h>

G_BEGIN_DECLS

/* These are the basic macros useful for the GObject framework */
#define NC_TYPE_CURVE             (nc_curve_get_type ())
#define NC_CURVE(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_CURVE, NcCurve))
#define NC_CURVE_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_CURVE, NcCurveClass))
#define NC_IS_CURVE(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_CURVE))
#define NC_IS_CURVE_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_CURVE))
#define NC_CURVE_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_CURVE, NcCurveClass))

/*
 * This is the class struct, it will contains everything that is
 * that is defined across instances.
 *
 */
typedef struct _NcCurveClass NcCurveClass;

/*
 * This is the instance struct, it will contains everything pertaining
 * to a given instance.
 *
 */
typedef struct _NcCurve NcCurve;

/*
 * This is the instance private struct, it will contains everything
 * pertaining that should never be seen outside of the compilation unit.
 *
 */
typedef struct _NcCurvePrivate NcCurvePrivate;

/*
 * Here we define a type to describe the virtual function which represets
 * the actual f(x).
 *
 */
typedef gdouble (*NcCurveF) (NcCurve *curve, const gdouble x);

/*
 * The class struct contains just one element `f' that must be assigned
 * by implementations. The first item is just the parent class structure.
 *
 */
struct _NcCurveClass
{
  /*< private >*/
  NcmModelClass parent_class;
  NcCurveF f;
};

/*
 * No elements needed apart from our private struct.
 * The first item is just the parent object structure.
 *
 */
struct _NcCurve
{
  /*< private >*/
  NcmModel parent_instance;
  NcCurvePrivate *priv;
};

/**
 * NcCurveImpl:
 * @NC_CURVE_F: The curve function $f(x)$.
 *
 * These flags are used to control which functions each child
 * implement. In this case the object is very simple and has
 * only one function that may be implement. In more complex cases
 * a child could implement only partially the abstract model
 * and it would the enumerator to inform which functions are
 * actually implement.
 *
 */
typedef enum _NcCurveImpl
{
  NC_CURVE_IMPL_f = 0,
} NcCurveImpl;

#define NC_CURVE_IMPL_ALL NCM_MODEL_CLASS_IMPL_ALL

GType nc_curve_get_type (void) G_GNUC_CONST;

/*
 * Since this is a abstract model we need to define its ID
 * in order to allow it to be found inside of a #NcmMSet
 * object.
 *
 */
NCM_MSET_MODEL_DECLARE_ID (nc_curve);

/*
 * This function must be used by implementations to set the virtual function f.
 *
 */
void nc_curve_set_f_impl (NcCurveClass *curve_class, NcCurveF f);

/*
 * The default methods that all objects should implement.
 * We also chose to implement a constructor nc_curve_new_from_name()
 * to allow creating any child directly.
 *
 */
NcCurve *nc_curve_new_from_name (const gchar *curve_name);
NcCurve *nc_curve_ref (NcCurve *curve);
void nc_curve_free (NcCurve *curve);
void nc_curve_clear (NcCurve **curve);

/*
 * Property accessors prototypes.
 *
 */
void nc_curve_set_xl (NcCurve *curve, const gdouble xl);
void nc_curve_set_xu (NcCurve *curve, const gdouble xu);
gdouble nc_curve_get_xl (NcCurve *curve);
gdouble nc_curve_get_xu (NcCurve *curve);

/*
 * The only quantity calculated by this model is the value
 * of $f(x)$, which is represented by the following virtual
 * method.
 *
 */
gdouble nc_curve_f (NcCurve *curve, const gdouble x);

G_END_DECLS

#endif /* _NC_CURVE_H_ */


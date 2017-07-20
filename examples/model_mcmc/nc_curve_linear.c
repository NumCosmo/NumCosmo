/***************************************************************************
 *            nc_curve_linear.c
 *
 *  Wed July 12 10:10:02 2017
 *  Copyright  2017  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * nc_curve_linear.c
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

/**
 * SECTION:nc_curve_linear
 * @title: NcCurveLinear
 * @short_description: Linear curve implementation of #NcCurve
 * 
 * This is a very simple implementation of #NcCurve where we 
 * will model $f(x) = a + b x$.
 * 
 * (This is a gtk-doc comment, which always start with two **)
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc_curve_linear.h"

/*
 * Properties enumerator, note that we added a last item PROP_SIZE,
 * this is useful to obtain the actual number of properties in the
 * object without hard coding it. PROP_0 is a special property of
 * the GObject system we should always be there.
 * 
 */
enum
{
  PROP_0,
  PROP_SIZE,
};

/*
 * Here we define the object GType using the macro G_DEFINE_TYPE.
 * This macro basically defines the function nc_curve_linear_get_type (void) and
 * everything necessary to define an object in the GLib type system. The
 * since it is not an ABSTRACT object, it creates a GType that can be instantiated.
 * 
 * The last argument provides the GType of the parent object.
 * 
 */
G_DEFINE_TYPE (NcCurveLinear, nc_curve_linear, NC_TYPE_CURVE);

/*
 * Nothing to do here.
 * 
 */
static void
nc_curve_linear_init (NcCurveLinear *curve_linear)
{
}

/*
 * Here we must de-allocate any memory allocated *inside* gobject framework,
 * i.e., we must unref any outside object contained in our object.
 * 
 * Nothing to do!
 * 
 */
static void
_nc_curve_linear_dispose (GObject *object)
{

  /* 
   * The following comment is always included to remark that at this point the parent
   * method must be called, chaining down the function call, i.e., first we finalize
   * the properties of the child, if any, then the parent and parent's parent, etc.
   * 
   */
  /* Chain up : end */
  G_OBJECT_CLASS (nc_curve_linear_parent_class)->dispose (object);
}

/*
 * Here we must de-allocate any memory allocated *outside* gobject framework.
 * Nothing to do!
 * 
 */
static void
_nc_curve_linear_finalize (GObject *object)
{

  /* 
   * The following comment is always included to remark that at this point the parent
   * method must be called, chaining down the function call, i.e., first we finalize
   * the properties of the child, if any, then the parent and parent's parent, etc.
   * 
   */
  /* Chain up : end */
  G_OBJECT_CLASS (nc_curve_linear_parent_class)->finalize (object);
}

/*
 * Note that we don't have additional properties and we don't need 
 * _nc_curve_linear_set_property/_nc_curve_linear_get_property
 * functions.
 * 
 */

/*
 * Prototype of our implementation of $f(x)$.
 * 
 */
static gdouble _nc_curve_linear_f (NcCurve *curve, const gdouble x);

/*
 * At _class_init we will define all properties and parameters we should
 * also include a default implementation for our virtual function `f'.
 * 
 */
static void
nc_curve_linear_class_init (NcCurveLinearClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcmModelClass *model_class = NCM_MODEL_CLASS (klass);
  NcCurveClass *curve_class  = NC_CURVE_CLASS (klass);

  /*
   * The usual GObject hooks assigned to GObjectClass.
   * 
   */  
  object_class->dispose     = &_nc_curve_linear_dispose;
  object_class->finalize    = &_nc_curve_linear_finalize;

  /*
   * First we set the model nick `CurveLinear-f' and its name.
   * 
   */
  ncm_model_class_set_name_nick (model_class, "CurveLinear-f", "NcCurveLinear");

  /*
   * Now we inform that we have NC_CURVE_LINEAR_SPARAM_LEN == 2 scalar parameters 
   * and 0 vector parameters. PROP_SIZE here is 1 and informs that we don't have any
   * additional property.
   * 
   */
  ncm_model_class_add_params (model_class, NC_CURVE_LINEAR_SPARAM_LEN, 0, PROP_SIZE);

  /*
   * Here we set the first parameter NC_CURVE_LINEAR_A, its symbol `$a$' coincide with 
   * its name `a', the allowed interval is $[-10, 10]$, a rough scale of variation is $1$.
   * 
   * Here we choose the absolute scale equal to zero. A absolute tolerance different 
   * from zero is only useful if you know a priori the order of magnitude of the
   * parameter.
   * 
   * Finally we choose the default value equal to $1$ and the default state of the parameter
   * is FIXED.
   * 
   */
  ncm_model_class_set_sparam (model_class, NC_CURVE_LINEAR_A, "a", "a",
                               -10.0, 10.0, 1.0,
                               0.0, 1.0,
                               NCM_PARAM_TYPE_FIXED);

  /*
   * Here we set the second parameter NC_CURVE_LINEAR_B, its symbol `$b$' coincide with 
   * its name `b', the allowed interval is $[1, 2]$, a rough scale of variation is $0.1$.
   * 
   * Here we choose the absolute scale equal to zero. A absolute tolerance different 
   * from zero is only useful if you know a priori the order of magnitude of the
   * parameter.
   * 
   * Finally we choose the default value equal to $1/2$ and the default state of the parameter
   * is FIXED.
   * 
   */
  ncm_model_class_set_sparam (model_class, NC_CURVE_LINEAR_B, "b", "b",
                               0.0, 2.0, 0.1,
                               0.0, 0.5,
                               NCM_PARAM_TYPE_FIXED);

  /* 
   * Check for errors in parameters initialization.
   * 
   */
  ncm_model_class_check_params_info (model_class);

  nc_curve_set_f_impl (curve_class, _nc_curve_linear_f);
}

/*
 * Here we introduce three helper macros to easily 
 * get the parameters.
 * 
 */
#define VECTOR (NCM_MODEL (curve_linear)->params)
#define A      (ncm_vector_get (VECTOR, NC_CURVE_LINEAR_A))
#define B      (ncm_vector_get (VECTOR, NC_CURVE_LINEAR_B))

/*
 * The actual implementation of $f(x) = a + b x$.
 * 
 */
static gdouble 
_nc_curve_linear_f (NcCurve *curve, const gdouble x)
{
  /*
   * Some notes:
   * 
   * - This is the implementation of NcCurveLinear, so in this case
   *   curve is actually a NcCurveLinear object, so we can recast
   *   it using NC_CURVE_LINEAR.
   * 
   * - We first assign the parameters to constant double variables,
   *   this can make the code more readable and do not make it slower,
   *   almost any compiler with optimization flags turned on will just 
   *   skip these steps.
   * 
   */
  NcCurveLinear *curve_linear = NC_CURVE_LINEAR (curve); 
  const gdouble a = A;
  const gdouble b = B;
  const gdouble f = a + b * x;

  return f;
}

/*
 * This defines the constructor for #NcCurveLinear 
 * 
 */ 
/**
 * nc_curve_linear_new:
 * @xl: lower bound $x_l$
 * @xu: upper bound $x_u$
 * 
 * Creates a new #NcCurveLinear.
 * 
 * Returns: (transfer full): a new #NcCurveLinear
 */ 
NcCurveLinear *
nc_curve_linear_new (const gdouble xl, const gdouble xu)
{
  /*
   * Here we use the actual GObject constructor to instantiate
   * an object of GType NC_TYPE_CURVE_LINEAR. Note that it 
   * inherits the #NcCurve properties xl and xu, and they must
   * be assigned during construction.
   * 
   */
  NcCurveLinear *curve_linear = g_object_new (NC_TYPE_CURVE_LINEAR,
                                              "xl", xl,
                                              "xu", xu,
                                              NULL);
  return curve_linear;
}

/*
 * The reference increasing function
 * 
 */
/**
 * nc_curve_linear_ref:
 * @curve_linear: a #NcCurveLinear
 * 
 * Increase reference count by one.
 * 
 * Returns: (transfer full): @curve_linear.
 */
NcCurveLinear *
nc_curve_linear_ref (NcCurveLinear *curve_linear)
{
  return g_object_ref (curve_linear);
}

/*
 * The reference decreasing function
 * 
 */
/**
 * nc_curve_linear_free:
 * @curve_linear: a #NcCurveLinear
 * 
 * Decrease reference count by one.
 * 
 */
void 
nc_curve_linear_free (NcCurveLinear *curve_linear)
{
  g_object_unref (curve_linear);
}

/*
 * This function decreases the reference count
 * by one only if *curve_linear != NULL, and in that case
 * it sets *curve_linear to NULL after decreasing the 
 * reference count. It is useful to use these 
 * functions in the dispose hooks.
 * 
 * 
 */
/**
 * nc_curve_linear_clear:
 * @curve_linear: a #NcCurveLinear
 * 
 * Decrease reference count by one if *@curve_linear != NULL
 * and sets @curve_linear to NULL. 
 * 
 */
void 
nc_curve_linear_clear (NcCurveLinear **curve_linear)
{
  g_clear_object (curve_linear);
}

/***************************************************************************
 *            nc_curve.c
 *
 *  Tue July 11 16:30:26 2017
 *  Copyright  2017  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_curve.c
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

/**
 * SECTION:nc_curve
 * @title: NcCurve
 * @short_description: Abstract class for curves!
 *
 * NcCurve is the abstract class designed to include the functions
 * that any simple curve should implement, see NcCurveImpl.
 * Its parent_class is NcmModel.
 *
 * (This is a gtk-doc comment, which always start with two **)
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc_curve.h"

/*
 * All implementations of NcCurve will required the function interval.
 * We could have let the interval itself as model parameters in order
 * to allow it to be fit too, this approach will be presented in a more
 * advanced example.
 *
 */
struct _NcCurvePrivate
{
  gdouble xl;
  gdouble xu;
};

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
  PROP_XL,
  PROP_XU,
  PROP_SIZE,
};

/*
 * Here we define the object GType using the macro G_DEFINE_ABSTRACT_TYPE.
 * This macro basically defines the function nc_curve_get_type (void) and
 * everything necessary to define an object in the GLib type system. The
 * ABSTRACT version of the macro creates a GType that cannot be instantiated
 * this means that to use this object we *must* define a child.
 *
 * The last argument provides the GType of the parent object.
 *
 */
G_DEFINE_ABSTRACT_TYPE (NcCurve, nc_curve, NCM_TYPE_MODEL);

static void
nc_curve_init (NcCurve *curve)
{
  /*
   * The first step is the creation of the private structure, all allocation and
   * de-allocation is automatically performed by the GObject framework.
   */
  curve->priv = G_TYPE_INSTANCE_GET_PRIVATE (curve, NC_TYPE_CURVE, NcCurvePrivate);

  /*
   * Here we initialize all structure member to null/zero.
   *
   */
  curve->priv->xl = 0.0;
  curve->priv->xu = 0.0;
}

/*
 * Here we call all property accessors in order to set property values.
 *
 */
static void
_nc_curve_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcCurve *curve = NC_CURVE (object);

  g_return_if_fail (NC_IS_CURVE (object));

  switch (prop_id)
  {
    case PROP_XL:
      nc_curve_set_xl (curve, g_value_get_double (value));
      break;
    case PROP_XU:
      nc_curve_set_xu (curve, g_value_get_double (value));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

/*
 * Here we call all property accessors in order to get property values.
 *
 */
static void
_nc_curve_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcCurve *curve = NC_CURVE (object);

  g_return_if_fail (NC_IS_CURVE (object));

  switch (prop_id)
  {
    case PROP_XL:
      g_value_set_double (value, nc_curve_get_xl (curve));
      break;
    case PROP_XU:
      g_value_set_double (value, nc_curve_get_xu (curve));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

/*
 * Here we must de-allocate any memory allocated *inside* gobject framework,
 * i.e., we must unref any outside object contained in our object.
 *
 * Nothing to do!
 *
 */
static void
_nc_curve_dispose (GObject *object)
{
  /*
   * The following comment is always included to remark that at this point the parent
   * method must be called, chaining down the function call, i.e., first we finalize
   * the properties of the child, if any, then the parent and parent's parent, etc.
   *
   */
  /* Chain up : end */
  G_OBJECT_CLASS (nc_curve_parent_class)->dispose (object);
}

/*
 * Here we must de-allocate any memory allocated *outside* gobject framework.
 * Nothing to do!
 *
 */
static void
_nc_curve_finalize (GObject *object)
{
  /*
   * The following comment is always included to remark that at this point the parent
   * method must be called, chaining down the function call, i.e., first we finalize
   * the properties of the child, if any, then the parent and parent's parent, etc.
   *
   */
  /* Chain up : end */
  G_OBJECT_CLASS (nc_curve_parent_class)->finalize (object);
}

/*
 * Registry the ID in NumCosmo model system.
 *
 */
NCM_MSET_MODEL_REGISTER_ID (nc_curve, NC_TYPE_CURVE);

/*
 * Prototype of the default implementation.
 *
 */
static gdouble _nc_curve_f (NcCurve *curve, const gdouble x);

/*
 * At _class_init we will define all properties and parameters we should
 * also include a default implementation for our virtual function `f'.
 *
 */
static void
nc_curve_class_init (NcCurveClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);
  NcmModelClass *model_class = NCM_MODEL_CLASS (klass);

  /*
   * Tells the type system that we have a private struct
   *
   */
  g_type_class_add_private (klass, sizeof (NcCurvePrivate));

  /*
   * The NcmModel class takes cares of the parameters, thus,
   * the set/get functions above must be set in NcmModelClass
   * structure.
   *
   */
  model_class->set_property = &_nc_curve_set_property;
  model_class->get_property = &_nc_curve_get_property;

  /*
   * The other functions are the usual GObject hooks and
   * must be assigned to GObjectClass.
   *
   */
  object_class->dispose  = &_nc_curve_dispose;
  object_class->finalize = &_nc_curve_finalize;

  /*
   * First we set the model nick `Curve-f' and its name
   *
   */
  ncm_model_class_set_name_nick (model_class, "Curve-f", "NcCurve");

  /*
   * Now we inform that we have 0 scalar parameters and 0 vector parameters.
   * Note that here we should include parameters that *all* implementations
   * would share, here we have none. Finallt, with the last argument we assert
   * that we have PROP_SIZE == 2 additional properties.
   *
   */
  ncm_model_class_add_params (model_class, 0, 0, PROP_SIZE);

  /*
   * Next step is to register the object in NumCosmo's model set object class.
   * We don't include a long description (first NULL below). The last arguments
   * FALSE and NCM_MSET_MODEL_MAIN determines that this object in non-stackable
   * and it is a MAIN model, both concepts will be discussed in an advanced example.
   *
   */
  ncm_mset_model_register_id (model_class,
                              "NcCurve",
                              "Curve model.",
                              NULL,
                              FALSE,
                              NCM_MSET_MODEL_MAIN);


  /*
   * The property PROP_XL is allowed in range [-G_MAXDOUBLE, G_MAXDOUBLE]
   * and its default value is 0.0. This property can be set only during
   * the object construction.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_XL,
                                   g_param_spec_double ("xl",
                                                        NULL,
                                                        "x lower bound",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /*
   * The property PROP_XU is allowed in range [-G_MAXDOUBLE, G_MAXDOUBLE]
   * and its default value is 1.0. This property can be set only during
   * the object construction.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_XU,
                                   g_param_spec_double ("xu",
                                                        NULL,
                                                        "x upper bound",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 1.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /*
   * Check for errors in parameters initialization.
   *
   */
  ncm_model_class_check_params_info (model_class);

  curve_class = &_nc_curve_f;
}

/*
 * Default implementation. It raises an error if the method is called
 * and the child didn't implement it.
 *
 */
static gdouble
_nc_curve_f (NcCurve *curve, const gdouble x)
{
  g_error ("nc_curve_f: model `%s' does not implement this function.",
           G_OBJECT_TYPE_NAME (curve));

  return 0.0;
}

/*
 * This creates a function to be used by the children to
 * implement this method.
 *
 */
/**
 * nc_curve_set_f_impl: (skip)
 * @curve_class: a #NcCurveClass
 * @f: function $f(x)$
 *
 * Sets the implementation of the curve $f(x)$.
 *
 */
NCM_MODEL_SET_IMPL_FUNC (NC_CURVE, NcCurve, nc_curve, NcCurveF, f)

/*
 * This defines a generic constructor for the subclasses
 *
 */
/**
 * nc_curve_new_from_name:
 * @curve_name: #NcCurve child type name
 *
 * Creates a new #NcCurve of the type described by @curve_name.
 *
 * Returns: (transfer full): a new #NcCurve
 */
NcCurve *
nc_curve_new_from_name (const gchar *curve_name)
{
  /*
   * We use the serialization object to transform a string into an instance.
   * This also allows us to chose the parameters through the string.
   *
   */
  GObject *obj = ncm_serialize_global_from_string (curve_name);

  /*
   * This gets the GType of this new instance.
   *
   */
  GType curve_type = G_OBJECT_TYPE (obj);

  /*
   * Check if the string represent an actual child of #NcCurve.
   *
   */
  if (!g_type_is_a (curve_type, NC_TYPE_CLUSTER_MASS))
    g_error ("nc_curve_new_from_name: NcCurve `%s' do not descend from `%s'.",
             curve_name, g_type_name (NC_TYPE_CURVE));

  /*
   * Returns the correct cast.
   *
   */
  return NC_CURVE (obj);
}

/*
 * The reference increasing function
 *
 */

/**
 * nc_curve_ref:
 * @curve: a #NcCurve
 *
 * Increase reference count by one.
 *
 * Returns: (transfer full): @curve.
 */
NcCurve *
nc_curve_ref (NcCurve *curve)
{
  return g_object_ref (curve);
}

/*
 * The reference decreasing function
 *
 */

/**
 * nc_curve_free:
 * @curve: a #NcCurve
 *
 * Decrease reference count by one.
 *
 */
void
nc_curve_free (NcCurve *curve)
{
  g_object_unref (curve);
}

/*
 * This function decreases the reference count
 * by one only if *curve != NULL, and in that case
 * it sets *curve to NULL after decreasing the
 * reference count. It is useful to use these
 * functions in the dispose hooks.
 *
 *
 */

/**
 * nc_curve_clear:
 * @curve: a #NcCurve
 *
 * Decrease reference count by one if *@curve != NULL
 * and sets @curve to NULL.
 *
 */
void
nc_curve_clear (NcCurve **curve)
{
  g_clear_object (curve);
}

/*
 * Below we implement the accessor functions.
 *
 */

/**
 * nc_curve_set_xl:
 * @curve: a #NcCurve
 * @xl: new $x$ lower bound
 *
 * Sets $x$ lower bound to @xl.
 *
 */
void
nc_curve_set_xl (NcCurve *curve, const gdouble xl)
{
  /*
   * Checking if the object is still in the unintitalized state,
   * if that's the case just assign the value to the private
   * struct.
   *
   */
  if ((curve->priv->xl == 0.0) && (curve->priv->xl == curve->priv->xu))
  {
    curve->priv->xl = xl;
  }
  else
  {
    /*
     * Otherwise assert that xl will be less than xu.
     *
     */
    g_assert_cmpfloat (xl, <, curve->priv->xu);
    curve->priv->xl = xl;
  }
}

/**
 * nc_curve_set_xu:
 * @curve: a #NcCurve
 * @xu: new $x$ upper bound
 *
 * Sets $x$ upper bound to @xu.
 *
 */
void
nc_curve_set_xu (NcCurve *curve, const gdouble xu)
{
  /*
   * Checking if the object is still in the unintitalized state,
   * if that's the case just assign the value to the private
   * struct.
   *
   */
  if ((curve->priv->xl == 0.0) && (curve->priv->xl == curve->priv->xu))
  {
    curve->priv->xu = xu;
  }
  else
  {
    /*
     * Otherwise assert that xl will be less than xu.
     *
     */
    g_assert_cmpfloat (curve->priv->xl, <, xu);
    curve->priv->xu = xu;
  }
}

/**
 * nc_curve_get_xl:
 * @curve: a #NcCurve
 *
 * Sets $x$ lower bound.
 *
 * Returns: $x_l$.
 */
gdouble
nc_curve_get_xl (NcCurve *curve)
{
  return curve->priv->xl;
}

/**
 * nc_curve_get_xu:
 * @curve: a #NcCurve
 *
 * Sets $x$ lower bound.
 *
 * Returns: $x_l$.
 */
gdouble
nc_curve_get_xu (NcCurve *curve)
{
  return curve->priv->xu;
}

/*
 * Finally, we implement the generic caller for the
 * virtual function `f'.
 *
 */

/**
 * nc_curve_f: (virtual f)
 * @curve: a #NcCurve
 * @x: $x$
 *
 * Computes $f(x)$.
 *
 * Returns: the value of $f(x)$.
 */
gdouble
nc_curve_f (NcCurve *curve, const gdouble x)
{
  /*
   * This function call by pointer guarantees that
   * the correct virtual function will be called.
   *
   */
  return NC_CURVE_GET_CLASS (curve)->f (curve, x);
}


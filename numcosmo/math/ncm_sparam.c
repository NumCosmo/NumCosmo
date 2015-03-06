/***************************************************************************
 *            ncm_sparam.c
 *
 *  Fri February 24 20:14:00 2012
 *  Copyright  2012  Sandro Dias Pinto Vitenti
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

/**
 * SECTION:ncm_sparam
 * @title: NcmSParam
 * @short_description: Properties of a scalar parameter.
 * 
 * This object comprises the necessary properties to define a scalar parameter.
 * It is used by #NcmModel to store the description of the scalar model parameters.
 * 
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_sparam.h"
#include "math/ncm_cfg.h"
#include "nc_hicosmo.h"
#include "ncm_enum_types.h"

enum
{
  PROP_0,
  PROP_NAME,
  PROP_SYMBOL,
  PROP_TYPE,
  PROP_LOWER_BOUND,
  PROP_UPPER_BOUND,
  PROP_SCALE,
  PROP_ABSOLUTE_TOLERANCE,
  PROP_DEFAULT_VALUE,
  PROP_FIT_TYPE
};

G_DEFINE_TYPE (NcmSParam, ncm_sparam, G_TYPE_OBJECT);

/**
 * ncm_sparam_new:
 * @name: #NcmSParam:name.
 * @symbol: #NcmSParam:symbol.
 * @lower_bound: value of #NcmSParam:lower-bound.
 * @upper_bound: value of #NcmSParam:upper-bound.
 * @scale: value of #NcmSParam:scale.
 * @abstol: value of #NcmSParam:absolute-tolerance.
 * @default_val: value of #NcmSParam:default-value.
 * @ftype: a #NcmParamType.
 *
 * This function allocates memory for a new #NcmSParam object and sets its properties to the values from
 * the input arguments.
 *
 * The @name parameter is restricted to the interval [@lower_bound, @upper_bound].
 * @scale is an initial step for the statistical algorithms.
 * @abstol is the absolute error tolerance of the parameter.
 * @ftype indicates if the parameter will be fitted or not.
 *
 * Returns: A new #NcmSParam.
 */
NcmSParam *
ncm_sparam_new (const gchar *name, const gchar *symbol, gdouble lower_bound, gdouble upper_bound, gdouble scale, gdouble abstol, gdouble default_val, NcmParamType ftype)
{
  NcmSParam *sparam = g_object_new (NCM_TYPE_SPARAM,
                                    "name",               name,
                                    "symbol",             symbol,
                                    "lower-bound",        lower_bound,
                                    "upper-bound",        upper_bound,
                                    "scale",              scale,
                                    "absolute-tolerance", abstol,
                                    "default-value",      default_val,
                                    "fit-type",           ftype,
                                    NULL);
  return sparam;
}

/**
 * ncm_sparam_copy:
 * @sparam: a #NcmSParam.
 *
 * Duplicates the #NcmSParam object setting the same values of the original propertities.
 *
 * Returns: (transfer full): A new #NcmSParam.
   */
NcmSParam *
ncm_sparam_copy (NcmSParam *sparam)
{
  return ncm_sparam_new (sparam->name, sparam->symbol, sparam->lower_bound,
                         sparam->upper_bound, sparam->scale, sparam->abstol,
                         sparam->default_val, sparam->ftype);
}

/**
 * ncm_sparam_free:
 * @sparam: a #NcmSParam.
 *
 * Atomically decrements the reference count of @sparam by one. If the reference count drops to 0,
 * all memory allocated by @sparam is released.
 *
 */
void
ncm_sparam_free (NcmSParam *sparam)
{
  g_object_unref (sparam);
}

/**
 * ncm_sparam_clear:
 * @sparam: a #NcmSParam.
 *
 * Atomically decrements the reference count of @sparam by one. If the reference count drops to 0,
 * all memory allocated by @sparam is released. Set the pointer to NULL.
 *
 */
void
ncm_sparam_clear (NcmSParam **sparam)
{
  g_clear_object (sparam);
}

/**
 * ncm_sparam_ref:
 * @sparam: a #NcmSParam.
 *
 * Atomically increase the reference count of @sparam by one.
 *
 * Returns: (transfer full): @sparam
   */
NcmSParam *
ncm_sparam_ref (NcmSParam *sparam)
{
  return g_object_ref (sparam);
}

/**
 * ncm_sparam_set_lower_bound:
 * @sparam: a #NcmSParam.
 * @lb: value of #NcmSParam:lower-bound.
 *
 * Sets the value @lb to the #NcmSParam:lower-bound property.
 *
 */
void
ncm_sparam_set_lower_bound (NcmSParam *sparam, const gdouble lb)
{
  g_assert (lb < sparam->upper_bound);
  sparam->lower_bound = lb;
}

/**
 * ncm_sparam_get_lower_bound:
 * @sparam: a #NcmSParam.
 *
 * Returns: the value of #NcmSParam:lower-bound property.
 */
gdouble
ncm_sparam_get_lower_bound (const NcmSParam *sparam)
{
  return sparam->lower_bound;
}

/**
 * ncm_sparam_set_upper_bound:
 * @sparam: a #NcmSParam.
 * @ub: value of #NcmSParam:upper-bound.
 *
 * Sets the value @ub to the #NcmSParam:upper-bound property.
 *
 */
void
ncm_sparam_set_upper_bound (NcmSParam *sparam, const gdouble ub)
{
  g_assert (ub > sparam->lower_bound);
  sparam->upper_bound = ub;
}

/**
 * ncm_sparam_get_upper_bound:
 * @sparam: a #NcmSParam.
 *
 * Returns: The value of #NcmSParam:upper-bound property.
 */
gdouble
ncm_sparam_get_upper_bound (const NcmSParam *sparam)
{
  return sparam->upper_bound;
}

/**
 * ncm_sparam_set_scale:
 * @sparam: a #NcmSParam.
 * @scale: value of #NcmSParam:scale.
 *
 * Sets the value @scale to the #NcmSParam:scale property.
 *
 */
void
ncm_sparam_set_scale (NcmSParam *sparam, const gdouble scale)
{
  g_assert (scale > 0);
  sparam->scale = scale;
}


/**
 * ncm_sparam_get_scale:
 * @sparam: a #NcmSParam.
 *
 * Returns: The value of #NcmSParam:scale property.
 */
gdouble
ncm_sparam_get_scale (const NcmSParam *sparam)
{
  return sparam->scale;
}

/**
 * ncm_sparam_set_absolute_tolerance:
 * @sparam: a #NcmSParam.
 * @abstol: value of #NcmSParam:absolute-tolerance.
 *
 * Sets the value @abstol to the #NcmSParam:absolute-tolerance property.
 *
 */
void
ncm_sparam_set_absolute_tolerance (NcmSParam *sparam, gdouble abstol)
{
  g_assert (abstol >= 0);
  sparam->abstol = abstol;
}

/**
 * ncm_sparam_get_absolute_tolerance:
 * @sparam: a #NcmSParam.
 *
 * Returns: the value of #NcmSParam:absolute_tolerance property.
 */
gdouble
ncm_sparam_get_absolute_tolerance (const NcmSParam *sparam)
{
  return sparam->abstol;
}

/**
 * ncm_sparam_set_default_value:
 * @sparam: a #NcmSParam.
 * @default_val: value of #NcmSParam:default-value.
 *
 * Sets the value @default_val to the #NcmSParam:default-value property.
 *
 */
void
ncm_sparam_set_default_value (NcmSParam *sparam, gdouble default_val)
{
  sparam->default_val = default_val;
}

/**
 * ncm_sparam_get_default_value:
 * @sparam: a #NcmSParam.
 *
 * Returns: the value of #NcmSParam:default-value property.
 */
gdouble
ncm_sparam_get_default_value (const NcmSParam *sparam)
{
  return sparam->default_val;
}

/**
 * ncm_sparam_set_fit_type:
 * @sparam: a #NcmSParam.
 * @ftype: a #NcmParamType.
 *
 * Sets the value @ftype to the #NcmSParam:fit-type property.
 *
 */
void
ncm_sparam_set_fit_type (NcmSParam *sparam, NcmParamType ftype)
{
  sparam->ftype = ftype;
}

/**
 * ncm_sparam_get_fit_type:
 * @sparam: a #NcmSParam.
 *
 * Returns: the #NcmParamType value of #NcmSParam:fit-type property.
 */
NcmParamType
ncm_sparam_get_fit_type (const NcmSParam *sparam)
{
  return sparam->ftype;
}

/**
 * ncm_sparam_take_name:
 * @sparam: a #NcmSParam.
 * @name: a string
 *
 * Take @name as the name string.
 * The caller doesn't have to free it any more.
 *
 */
void
ncm_sparam_take_name (NcmSParam *sparam, gchar *name)
{
  g_free (sparam->name);
  sparam->name = name;
}

/**
 * ncm_sparam_take_symbol:
 * @sparam: a #NcmSParam.
 * @symbol: a string
 *
 * Take @symbol as the symbol string.
 * The caller doesn't have to free it any more.
 *
 */
void
ncm_sparam_take_symbol (NcmSParam *sparam, gchar *symbol)
{
  g_free (sparam->symbol);
  sparam->symbol = symbol;
}

/**
 * ncm_sparam_name:
 * @sparam: a #NcmSParam.
 *
 * Returns: the internal name string. The caller must not free it.
 */
const gchar *
ncm_sparam_name (const NcmSParam *sparam)
{
  return sparam->name;
}

/**
 * ncm_sparam_symbol:
 * @sparam: a #NcmSParam.
 *
 * Returns: the internal symbol string. The caller must not free it.
 */
const gchar *
ncm_sparam_symbol (const NcmSParam *sparam)
{
  return sparam->symbol;
}

static void
ncm_sparam_init (NcmSParam *sparam)
{
  sparam->name = NULL;
  sparam->symbol = NULL;
  sparam->lower_bound = -G_MAXDOUBLE;
  sparam->upper_bound =  G_MAXDOUBLE;
  sparam->scale = 0.0;
  sparam->abstol = NC_HICOSMO_DEFAULT_PARAMS_ABSTOL;
  sparam->default_val = 0.0;
  sparam->ftype = NCM_PARAM_TYPE_FREE;
}

static void
_ncm_sparam_finalize (GObject *object)
{
  NcmSParam *sparam = NCM_SPARAM (object);
  if (sparam->name)
    g_free (sparam->name);
  if (sparam->symbol)
    g_free (sparam->symbol);

  G_OBJECT_CLASS (ncm_sparam_parent_class)->finalize (object);
}

static void
_ncm_sparam_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmSParam *sparam = NCM_SPARAM (object);
  g_return_if_fail (NCM_IS_SPARAM (object));

  switch (prop_id)
  {
    case PROP_NAME:
      sparam->name = g_value_dup_string (value);
      break;
    case PROP_SYMBOL:
      sparam->symbol = g_value_dup_string (value);
      break;
    case PROP_LOWER_BOUND:
      ncm_sparam_set_lower_bound (sparam, g_value_get_double (value));
      break;
    case PROP_UPPER_BOUND:
      ncm_sparam_set_upper_bound (sparam, g_value_get_double (value));
      break;
    case PROP_SCALE:
      ncm_sparam_set_scale (sparam, g_value_get_double (value));
      break;
    case PROP_ABSOLUTE_TOLERANCE:
      ncm_sparam_set_absolute_tolerance (sparam, g_value_get_double (value));
      break;
    case PROP_DEFAULT_VALUE:
      ncm_sparam_set_default_value (sparam, g_value_get_double (value));
      break;
    case PROP_FIT_TYPE:
      ncm_sparam_set_fit_type (sparam, g_value_get_enum (value));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_sparam_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmSParam *sparam = NCM_SPARAM (object);
  g_return_if_fail (NCM_IS_SPARAM (object));

  switch (prop_id)
  {
    case PROP_NAME:
      g_value_set_string (value, sparam->name);
      break;
    case PROP_SYMBOL:
      g_value_set_string (value, sparam->symbol);
      break;
    case PROP_LOWER_BOUND:
      g_value_set_double (value, ncm_sparam_get_lower_bound (sparam));
      break;
    case PROP_UPPER_BOUND:
      g_value_set_double (value, ncm_sparam_get_upper_bound (sparam));
      break;
    case PROP_SCALE:
      g_value_set_double (value, ncm_sparam_get_scale (sparam));
      break;
    case PROP_ABSOLUTE_TOLERANCE:
      g_value_set_double (value, ncm_sparam_get_absolute_tolerance (sparam));
      break;
    case PROP_DEFAULT_VALUE:
      g_value_set_double (value, ncm_sparam_get_default_value (sparam));
      break;
    case PROP_FIT_TYPE:
      g_value_set_enum (value, ncm_sparam_get_fit_type (sparam));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
ncm_sparam_class_init (NcmSParamClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = _ncm_sparam_set_property;
  object_class->get_property = _ncm_sparam_get_property;
  object_class->finalize = _ncm_sparam_finalize;

  /**
   * NcmSParam:name:
   *
   * The parameter' s name must be a string written using only ASCII and -.
   */
  g_object_class_install_property (object_class,
                                   PROP_NAME,
                                   g_param_spec_string ("name",
                                                        NULL,
                                                        "Name (only ASCII plus -)",
                                                        NULL,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmSParam:symbol:
   *
   * Parameter's name written in a usual form (including latex).
     */
  g_object_class_install_property (object_class,
                                   PROP_SYMBOL,
                                   g_param_spec_string ("symbol",
                                                        NULL,
                                                        "Symbol (latex)",
                                                        NULL,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmSParam:lower-bound:
   *
   * Lower parameter threshold whose value is restricted to [-G_MAXDOUBLE, G_MAXDOUBLE].
     */
  g_object_class_install_property (object_class,
                                   PROP_LOWER_BOUND,
                                   g_param_spec_double ("lower-bound",
                                                        NULL,
                                                        "Lower bound",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmSParam:upper-bound:
   *
   * Upper parameter threshold whose value is restricted to [-G_MAXDOUBLE, G_MAXDOUBLE].
     */
  g_object_class_install_property (object_class,
                                   PROP_UPPER_BOUND,
                                   g_param_spec_double ("upper-bound",
                                                        NULL,
                                                        "Upper bound",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmSParam:scale:
   *
   * Scale, whose value is restricted to [0, G_MAXDOUBLE], is the step used by #NcmFit to increment the value of the parameter.
     */
  g_object_class_install_property (object_class,
                                   PROP_SCALE,
                                   g_param_spec_double ("scale",
                                                        NULL,
                                                        "Scale in which the model varies",
                                                        0.0, G_MAXDOUBLE, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmSParam:absolute-tolerance:
   *
   * Absolute tolerance, whose value is restricted to [0, G_MAXDOUBLE], is the size of the error used by #NcmFit.
     */
  g_object_class_install_property (object_class,
                                   PROP_ABSOLUTE_TOLERANCE,
                                   g_param_spec_double ("absolute-tolerance",
                                                        NULL,
                                                        "Absolute tolerance",
                                                        0.0, G_MAXDOUBLE, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  /**
   * NcmSParam:default-value:
   *
   * Parameter's default value.
   */
  g_object_class_install_property (object_class,
                                   PROP_DEFAULT_VALUE,
                                   g_param_spec_double ("default-value",
                                                        NULL,
                                                        "Default value",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmSParam:fit-type:
   *
   * Parameter's fit type: FIXED or FREE.
   */
  g_object_class_install_property (object_class,
                                   PROP_FIT_TYPE,
                                   g_param_spec_enum ("fit-type",
                                                      NULL,
                                                      "Fit Type",
                                                      NCM_TYPE_PARAM_TYPE, NCM_PARAM_TYPE_FREE,
                                                      G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

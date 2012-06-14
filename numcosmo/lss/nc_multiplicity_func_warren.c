/***************************************************************************
 *            nc_multiplicity_func_warren.c
 *
 *  Mon Jun 28 15:09:13 2010
 *  Copyright  2010  Mariana Penna Lima
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Mariana Penna Lima 2012 <pennalima@gmail.com>
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
 * SECTION:nc_multiplicity_func_warren
 * @title: Warren Multiplicity Function
 * @short_description: Dark Matter Halo FIXME
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include <numcosmo/numcosmo.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_integration.h>

G_DEFINE_TYPE (NcMultiplicityFuncWarren, nc_multiplicity_func_warren, NC_TYPE_MULTIPLICITY_FUNC);

enum
{
  PROP_0,
  PROP_A,
  PROP_A1,
  PROP_B,
  PROP_C
};

/**
 * nc_multiplicity_func_warren_new:
 * @A: FIXME
 * @a: FIXME
 * @b: FIXME
 * @c: FIXME 
 *   
 * FIXME
 *
 * Returns: A new #NcMultiplicityFunc.
 */
NcMultiplicityFunc *
nc_multiplicity_func_warren_new (gdouble A, gdouble a, gdouble b, gdouble c)
{
  return g_object_new (NC_TYPE_MULTIPLICITY_FUNC_WARREN,
                       "A", A,
                       "a", a, 
                       "b", b,
                       "c", c,
                       NULL);
}

static gdouble
_nc_multiplicity_func_warren_eval (NcMultiplicityFunc *mulf, NcHICosmo *model, gdouble sigma, gdouble z)   /* $f(\sigma)$ Warren: astro-ph/0506395 */
{
  NcMultiplicityFuncWarren *mulf_warren = NC_MULTIPLICITY_FUNC_WARREN (mulf);
  gdouble f_Warren = mulf_warren->A * (pow(sigma, - mulf_warren->a) + mulf_warren->b) * exp(-(mulf_warren->c)/ (sigma * sigma) );

  return f_Warren;
}

/**
 * nc_multiplicity_func_warren_set_A:
 * @mulf_warren: a #NcMultiplicityFuncWarren.
 * @A: value of #NcMultiplicityFuncWarren:A.
 *
 * Sets the value @A to the #NcMultiplicityFuncWarren:A property.
 *
 */
void
nc_multiplicity_func_warren_set_A (NcMultiplicityFuncWarren *mulf_warren, gdouble A)
{
  g_assert (A >= 0);
  mulf_warren->A = A;
}

/**
 * nc_multiplicity_func_warren_get_A:
 * @mulf_warren: a #NcMultiplicityFuncWarren.
 *
 * Returns: the value of #NcMultiplicityFuncWarren:A property.
 */
gdouble
nc_multiplicity_func_warren_get_A (const NcMultiplicityFuncWarren *mulf_warren)
{
  return mulf_warren->A;
}

/**
 * nc_multiplicity_func_warren_set_a:
 * @mulf_warren: a #NcMultiplicityFuncWarren.
 * @a: value of #NcMultiplicityFuncWarren:a.
 *
 * Sets the value @a to the #NcMultiplicityFuncWarren:a property.
 *
 */
void
nc_multiplicity_func_warren_set_a (NcMultiplicityFuncWarren *mulf_warren, gdouble a)
{
  g_assert (a >= 0);
  mulf_warren->a = a;
}

/**
 * nc_multiplicity_func_warren_get_a:
 * @mulf_warren: a #NcMultiplicityFuncWarren.
 *
 * Returns: the value of #NcMultiplicityFuncWarren:a property.
 */
gdouble
nc_multiplicity_func_warren_get_a (const NcMultiplicityFuncWarren *mulf_warren)
{
  return mulf_warren->a;
}

/**
 * nc_multiplicity_func_warren_set_b:
 * @mulf_warren: a #NcMultiplicityFuncWarren.
 * @b: value of #NcMultiplicityFuncWarren:b.
 *
 * Sets the value @b to the #NcMultiplicityFuncWarren:b property.
 *
 */
void
nc_multiplicity_func_warren_set_b (NcMultiplicityFuncWarren *mulf_warren, gdouble b)
{
  g_assert (b >= 0);
  mulf_warren->b = b;
}

/**
 * nc_multiplicity_func_warren_get_b:
 * @mulf_warren: a #NcMultiplicityFuncWarren.
 *
 * Returns: the value of #NcMultiplicityFuncWarren:b property.
 */
gdouble
nc_multiplicity_func_warren_get_b (const NcMultiplicityFuncWarren *mulf_warren)
{
  return mulf_warren->b;
}

/**
 * nc_multiplicity_func_warren_set_c:
 * @mulf_warren: a #NcMultiplicityFuncWarren.
 * @c: value of #NcMultiplicityFuncWarren:c.
 *
 * Sets the value @c to the #NcMultiplicityFuncWarren:c property.
 *
 */
void
nc_multiplicity_func_warren_set_c (NcMultiplicityFuncWarren *mulf_warren, gdouble c)
{
  g_assert (c >= 0);
  mulf_warren->c = c;
}

/**
 * nc_multiplicity_func_warren_get_c:
 * @mulf_warren: a #NcMultiplicityFuncWarren.
 *
 * Returns: the value of #NcMultiplicityFuncWarren:c property.
 */
gdouble
nc_multiplicity_func_warren_get_c (const NcMultiplicityFuncWarren *mulf_warren)
{
  return mulf_warren->c;
}

// _NC_MULTIPLICITY_FUNCTION_WARREN_DATASET_0506395 = {0.7234, 1.625, 0.2538, 1.1982};

static void
nc_multiplicity_func_warren_init (NcMultiplicityFuncWarren *mulf_warren)
{
  /* TODO: Add initialization code here */
  mulf_warren->A = 0.7234;
  mulf_warren->a = 1.625;
  mulf_warren->b = 0.2538;
  mulf_warren->c = 1.1982;
}

static void
_nc_multiplicity_func_warren_finalize (GObject *object)
{
  /* TODO: Add deinitalization code here */

  G_OBJECT_CLASS (nc_multiplicity_func_warren_parent_class)->finalize (object);
}

static void
_nc_multiplicity_func_warren_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcMultiplicityFuncWarren *mulf_warren = NC_MULTIPLICITY_FUNC_WARREN (object);
  g_return_if_fail (NC_IS_MULTIPLICITY_FUNC_WARREN (object));

  switch (prop_id)
  {
	case PROP_A:
	  mulf_warren->A = g_value_get_double (value);
	  break;
	case PROP_A1:
	  mulf_warren->a = g_value_get_double (value);
	  break;
	case PROP_B:
	  mulf_warren->b = g_value_get_double (value);
	  break;
	case PROP_C:
	  mulf_warren->c = g_value_get_double (value);
	  break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_multiplicity_func_warren_get_property (GObject * object, guint prop_id, GValue * value, GParamSpec * pspec)
{
  NcMultiplicityFuncWarren *mulf_warren = NC_MULTIPLICITY_FUNC_WARREN (object);
  g_return_if_fail (NC_IS_MULTIPLICITY_FUNC_WARREN (object));

  switch (prop_id)
  {
	case PROP_A:
	  g_value_set_double (value, mulf_warren->A);
	  break;
	case PROP_A1:
	  g_value_set_double (value, mulf_warren->a);
	  break;
	case PROP_B:
	  g_value_set_double (value, mulf_warren->b);
	  break;
    case PROP_C:
      g_value_set_double (value, mulf_warren->c);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_multiplicity_func_warren_class_init (NcMultiplicityFuncWarrenClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcMultiplicityFuncClass* parent_class = NC_MULTIPLICITY_FUNC_CLASS (klass);

  parent_class->eval = &_nc_multiplicity_func_warren_eval;

  object_class->finalize = _nc_multiplicity_func_warren_finalize;
  object_class->set_property = _nc_multiplicity_func_warren_set_property;
  object_class->get_property = _nc_multiplicity_func_warren_get_property;

  /**
   * NcMultiplicityFuncWarren:A:
   *
   * FIXME Set correct values (limits)
   */
  g_object_class_install_property (object_class,
                                   PROP_A,
                                   g_param_spec_double ("A",
                                                        NULL,
                                                        "A",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 0.7234,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY |G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  /**
   * NcMultiplicityFuncWarren:a:
   *
   * FIXME Set correct values (limits)
   */
  g_object_class_install_property (object_class,
                                   PROP_A1,
                                   g_param_spec_double ("a",
                                                        NULL,
                                                        "a",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 1.625,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY |G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  /**
   * NcMultiplicityFuncWarren:b:
   *
   * FIXME Set correct values (limits)
   */
  g_object_class_install_property (object_class,
                                   PROP_B,
                                   g_param_spec_double ("b",
                                                        NULL,
                                                        "b",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 0.2538,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY |G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  /**
   * NcMultiplicityFuncWarren:c:
   *
   * FIXME Set correct values (limits)
   */
  g_object_class_install_property (object_class,
                                   PROP_C,
                                   g_param_spec_double ("c",
                                                        NULL,
                                                        "c",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 1.1982,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY |G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}


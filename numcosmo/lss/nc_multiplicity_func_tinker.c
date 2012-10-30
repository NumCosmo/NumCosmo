/***************************************************************************
 *            nc_multiplicity_func_tinker.c
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
 * SECTION:nc_multiplicity_func_tinker
 * @title: Tinker Multiplicity Function
 * @short_description: Dark Matter Halo FIXME
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "lss/nc_multiplicity_func_tinker.h"

G_DEFINE_TYPE (NcMultiplicityFuncTinker, nc_multiplicity_func_tinker, NC_TYPE_MULTIPLICITY_FUNC);

enum
{
  PROP_0,
  PROP_A0,
  PROP_A1,
  PROP_B0,
  PROP_C,
  PROP_DELTA
};

/**
 * nc_multiplicity_func_tinker_new:
 * @A0: FIXME
 * @a0: FIXME 
 * @b0: FIXME
 * @c: FIXME
 * @Delta: FIXME 
 *   
 * FIXME
 *
 * Returns: A new #NcMultiplicityFunc.
 */
NcMultiplicityFunc *
nc_multiplicity_func_tinker_new (gdouble A0, gdouble a0, gdouble b0, gdouble c, gdouble Delta)
{
  return g_object_new (NC_TYPE_MULTIPLICITY_FUNC_TINKER,
                       "A0", A0,
                       "a0", a0,
                       "b0", b0, 
                       "c", c,
                       "Delta", Delta,
                       NULL);
}

static gdouble
_nc_multiplicity_func_tinker_eval (NcMultiplicityFunc *mulf, NcHICosmo *model, gdouble sigma, gdouble z)   /* $f(\sigma)$ Tinker: astro-ph/0803.2706 */
{
  NcMultiplicityFuncTinker *mulf_tinker = NC_MULTIPLICITY_FUNC_TINKER (mulf);
  const gdouble A = mulf_tinker->A0 * pow(1.0 + z, -0.14);
  const gdouble a = mulf_tinker->a0 * pow(1.0 + z, -0.06);
  const gdouble log10alpha = -pow(0.75 / log10 (mulf_tinker->Delta / 75.0), 1.2);
  const gdouble alpha = pow(10.0, log10alpha);
  const gdouble b = mulf_tinker->b0 * pow(1.0 + z, -alpha);
  const gdouble f_Tinker = A * (pow(sigma/b, -a) + 1.0) * exp(-(mulf_tinker->c) / (sigma * sigma));

  //	printf ("%.15g %.15g %.15g %.15g | %.15g %.15g %.15g %.15g\n", mulf_tinker->A0, mulf_tinker->a0, mulf_tinker->b0, mulf_tinker->c,
  //	        A, a, b, mulf_tinker->c);
  return f_Tinker;
}

/**
 * nc_multiplicity_func_tinker_set_A0:
 * @mulf_tinker: a #NcMultiplicityFuncTinker.
 * @A0: value of #NcMultiplicityFuncTinker:A0.
 *
 * Sets the value @A0 to the #NcMultiplicityFuncTinker:A0 property.
 *
 */
void
nc_multiplicity_func_tinker_set_A0 (NcMultiplicityFuncTinker *mulf_tinker, gdouble A0)
{
  g_assert (A0 >= 0);
  mulf_tinker->A0 = A0;
}

/**
 * nc_multiplicity_func_tinker_get_A0:
 * @mulf_tinker: a #NcMultiplicityFuncTinker.
 *
 * Returns: the value of #NcMultiplicityFuncTinker:A0 property.
 */
gdouble
nc_multiplicity_func_tinker_get_A0 (const NcMultiplicityFuncTinker *mulf_tinker)
{
  return mulf_tinker->A0;
}

/**
 * nc_multiplicity_func_tinker_set_a0:
 * @mulf_tinker: a #NcMultiplicityFuncTinker.
 * @a0: value of #NcMultiplicityFuncTinker:a0.
 *
 * Sets the value @a0 to the #NcMultiplicityFuncTinker:a0 property.
 *
 */
void
nc_multiplicity_func_tinker_set_a0 (NcMultiplicityFuncTinker *mulf_tinker, gdouble a0)
{
  g_assert (a0 >= 0);
  mulf_tinker->a0 = a0;
}

/**
 * nc_multiplicity_func_tinker_get_a0:
 * @mulf_tinker: a #NcMultiplicityFuncTinker.
 *
 * Returns: the value of #NcMultiplicityFuncTinker:a0 property.
 */
gdouble
nc_multiplicity_func_tinker_get_a0 (const NcMultiplicityFuncTinker *mulf_tinker)
{
  return mulf_tinker->a0;
}

/**
 * nc_multiplicity_func_tinker_set_b0:
 * @mulf_tinker: a #NcMultiplicityFuncTinker.
 * @b0: value of #NcMultiplicityFuncTinker:b0.
 *
 * Sets the value @b0 to the #NcMultiplicityFuncTinker:b0 property.
 *
 */
void
nc_multiplicity_func_tinker_set_b0 (NcMultiplicityFuncTinker *mulf_tinker, gdouble b0)
{
  g_assert (b0 >= 0);
  mulf_tinker->b0 = b0;
}

/**
 * nc_multiplicity_func_tinker_get_b0:
 * @mulf_tinker: a #NcMultiplicityFuncTinker.
 *
 * Returns: the value of #NcMultiplicityFuncTinker:b0 property.
 */
gdouble
nc_multiplicity_func_tinker_get_b0 (const NcMultiplicityFuncTinker *mulf_tinker)
{
  return mulf_tinker->b0;
}

/**
 * nc_multiplicity_func_tinker_set_c:
 * @mulf_tinker: a #NcMultiplicityFuncTinker.
 * @c: value of #NcMultiplicityFuncTinker:c.
 *
 * Sets the value @c to the #NcMultiplicityFuncTinker:c property.
 *
 */
void
nc_multiplicity_func_tinker_set_c (NcMultiplicityFuncTinker *mulf_tinker, gdouble c)
{
  g_assert (c >= 0);
  mulf_tinker->c = c;
}

/**
 * nc_multiplicity_func_tinker_get_c:
 * @mulf_tinker: a #NcMultiplicityFuncTinker.
 *
 * Returns: the value of #NcMultiplicityFuncTinker:c property.
 */
gdouble
nc_multiplicity_func_tinker_get_c (const NcMultiplicityFuncTinker *mulf_tinker)
{
  return mulf_tinker->c;
}

/**
 * nc_multiplicity_func_tinker_set_Delta:
 * @mulf_tinker: a #NcMultiplicityFuncTinker.
 * @Delta: value of #NcMultiplicityFuncTinker:Delta.
 *
 * Sets the value @Delta to the #NcMultiplicityFuncTinker:Delta property.
 *
 */
void
nc_multiplicity_func_tinker_set_Delta (NcMultiplicityFuncTinker *mulf_tinker, gdouble Delta)
{
  g_assert (Delta >= 0);
  mulf_tinker->Delta = Delta;
}

/**
 * nc_multiplicity_func_tinker_get_Delta:
 * @mulf_tinker: a #NcMultiplicityFuncTinker.
 *
 * Returns: the value of #NcMultiplicityFuncTinker:Delta property.
 */
gdouble
nc_multiplicity_func_tinker_get_Delta (const NcMultiplicityFuncTinker *mulf_tinker)
{
  return mulf_tinker->Delta;
}

// _NC_MULTIPLICITY_FUNCTION_TINKER_DATASET_0803_2706_DELTA200 = {0.186, 1.47, 2.57, 1.19, 200.0};
// _NC_MULTIPLICITY_FUNCTION_TINKER_DATASET_0803_2706_DELTA800 = {0.248, 1.87, 1.59, 1.58, 800.0};
// _NC_MULTIPLICITY_FUNCTION_TINKER_DATASET_0803_2706_DELTA3200 = {0.26, 2.66, 1.41, 2.44, 3200.0};

static void
nc_multiplicity_func_tinker_init (NcMultiplicityFuncTinker *mulf_tinker)
{
  /* TODO: Add initialization code here */
  mulf_tinker->A0 = 0.186;
  mulf_tinker->a0 = 1.47;
  mulf_tinker->b0 = 2.57;
  mulf_tinker->c = 1.19;
  mulf_tinker->Delta = 200.0;
}

static void
_nc_multiplicity_func_tinker_finalize (GObject *object)
{
  /* TODO: Add deinitalization code here */

  G_OBJECT_CLASS (nc_multiplicity_func_tinker_parent_class)->finalize (object);
}

static void
_nc_multiplicity_func_tinker_set_property (GObject * object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcMultiplicityFuncTinker *mulf_tinker = NC_MULTIPLICITY_FUNC_TINKER (object);
  g_return_if_fail (NC_IS_MULTIPLICITY_FUNC_TINKER (object));

  switch (prop_id)
  {
	case PROP_A0:
	  mulf_tinker->A0 = g_value_get_double (value);
	  break;
	case PROP_A1:
	  mulf_tinker->a0 = g_value_get_double (value);
	  break;
	case PROP_B0:
	  mulf_tinker->b0 = g_value_get_double (value);
	  break;
	case PROP_C:
	  mulf_tinker->c = g_value_get_double (value);
	  break;
	case PROP_DELTA:
	  mulf_tinker->Delta = g_value_get_double (value);
	  break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_multiplicity_func_tinker_get_property (GObject * object, guint prop_id, GValue * value, GParamSpec * pspec)
{
  NcMultiplicityFuncTinker *mulf_tinker = NC_MULTIPLICITY_FUNC_TINKER (object);
  g_return_if_fail (NC_IS_MULTIPLICITY_FUNC_TINKER (object));

  switch (prop_id)
  {
	case PROP_A0:
	  g_value_set_double (value, mulf_tinker->A0);
	  break;
	case PROP_A1:
	  g_value_set_double (value, mulf_tinker->a0);
	  break;
	case PROP_B0:
	  g_value_set_double (value, mulf_tinker->b0);
	  break;
	case PROP_C:
	  g_value_set_double (value, mulf_tinker->c);
	  break;
    case PROP_DELTA:
      g_value_set_double (value, mulf_tinker->Delta);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_multiplicity_func_tinker_class_init (NcMultiplicityFuncTinkerClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcMultiplicityFuncClass* parent_class = NC_MULTIPLICITY_FUNC_CLASS (klass);

  parent_class->eval = &_nc_multiplicity_func_tinker_eval;

  object_class->finalize = _nc_multiplicity_func_tinker_finalize;
  object_class->set_property = _nc_multiplicity_func_tinker_set_property;
  object_class->get_property = _nc_multiplicity_func_tinker_get_property;

  /**
   * NcMultiplicityFuncTinker:A0:
   *
   * FIXME Set correct values (limits)
   */
  g_object_class_install_property (object_class,
                                   PROP_A0,
                                   g_param_spec_double ("A0",
                                                        NULL,
                                                        "A0",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 0.186,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY |G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  /**
   * NcMultiplicityFuncTinker:a0:
   *
   * FIXME Set correct values (limits)
   */
  g_object_class_install_property (object_class,
                                   PROP_A1,
                                   g_param_spec_double ("a0",
                                                        NULL,
                                                        "a0",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 1.47,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY |G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcMultiplicityFuncTinker:b0:
   *
   * FIXME Set correct values (limits)
   */
  g_object_class_install_property (object_class,
                                   PROP_B0,
                                   g_param_spec_double ("b0",
                                                        NULL,
                                                        "b0",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 2.57,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY |G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  /**
   * NcMultiplicityFuncTinker:c:
   *
   * FIXME Set correct values (limits)
   */
  g_object_class_install_property (object_class,
                                   PROP_C,
                                   g_param_spec_double ("c",
                                                        NULL,
                                                        "c",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 1.19,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY |G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  /**
   * NcMultiplicityFuncTinker:Delta:
   *
   * FIXME Set correct values (limits)
   */
  g_object_class_install_property (object_class,
                                   PROP_DELTA,
                                   g_param_spec_double ("Delta",
                                                        NULL,
                                                        "Delta",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 200.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY |G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}


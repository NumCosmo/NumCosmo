/***************************************************************************
 *            nc_halo_bias_st_spher.c
 *
 *  Tue June 28 15:41:57 2011
 *  Copyright  2011  Mariana Penna Lima
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
 * SECTION:nc_halo_bias_st_spher
 * @title: NcHaloBiasSTSpher
 * @short_description: Sheth-Tormen spherical halo bias function type.
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "lss/nc_halo_bias_st_spher.h"
#include "math/ncm_cfg.h"
#include "math/ncm_util.h"

G_DEFINE_TYPE (NcHaloBiasSTSpher, nc_halo_bias_st_spher, NC_TYPE_HALO_BIAS);

enum
{
  PROP_0,
  PROP_DELTA_C,
  PROP_A,
  PROP_P
};

/**
 * nc_halo_bias_st_spher_new:
 * @delta_c: FIXME
 * @a: FIXME
 * @p: FIXME
 *
 * FIXME
 *
 * Returns: A new #NcHaloBias.
 */
NcHaloBias *
nc_halo_bias_st_spher_new (gdouble delta_c, gdouble a, gdouble p)
{
  return g_object_new (NC_TYPE_HALO_BIAS_ST_SPHER,
                       "critical-delta", delta_c,
                       "a", a,
                       "p", p,
                       NULL);
}

static gdouble
_nc_halo_bias_st_spher_eval (NcHaloBias *biasf, gdouble sigma, gdouble z)
{
  NcHaloBiasSTSpher *bias_st_spher = NC_HALO_BIAS_ST_SPHER (biasf);
  const gdouble a = bias_st_spher->a;
  const gdouble p = bias_st_spher->p;
  gdouble x = bias_st_spher->delta_c / sigma;
  gdouble x2 = x * x;
  gdouble b_ST_spher = 1.0  + ((a * x2 - 1.0) + (2.0 * p) / (1.0 + pow (a * x2, p))) / bias_st_spher->delta_c;

  NCM_UNUSED (z);
//  printf ("a = %.5g, p=%.5g, delta_c= %.5g\n", a, p, bias_st_spher->delta_c);

  return b_ST_spher;
}

/**
 * nc_halo_bias_st_spher_set_delta_c:
 * @biasf_st_spher: a #NcHaloBiasSTSpher.
 * @delta_c: value of #NcHaloBiasSTSpher:critical-delta.
 *
 * Sets the value @delta_c to the #NcHaloBiasSTSpher:critical-delta property.
 *
 */
void
nc_halo_bias_st_spher_set_delta_c (NcHaloBiasSTSpher *biasf_st_spher, gdouble delta_c)
{
  g_assert (delta_c >= 0);
  biasf_st_spher->delta_c = delta_c;
}

/**
 * nc_halo_bias_st_spher_get_delta_c:
 * @biasf_st_spher: a #NcHaloBiasSTSpher.
 *
 * Returns: the value of #NcHaloBiasSTSpher:critical_delta property.
 */
gdouble
nc_halo_bias_st_spher_get_delta_c (const NcHaloBiasSTSpher *biasf_st_spher)
{
  return biasf_st_spher->delta_c;
}

/**
 * nc_halo_bias_st_spher_set_a:
 * @biasf_st_spher: a #NcHaloBiasSTSpher.
 * @a: value of #NcHaloBiasSTSpher:a.
 *
 * Sets the value @a to the #NcHaloBiasSTSpher:a property.
 *
 */
void
nc_halo_bias_st_spher_set_a (NcHaloBiasSTSpher *biasf_st_spher, gdouble a)
{
  g_assert (a >= 0);
  biasf_st_spher->a = a;
}

/**
 * nc_halo_bias_st_spher_get_a:
 * @biasf_st_spher: a #NcHaloBiasSTSpher.
 *
 * Returns: the value of #NcHaloBiasSTSpher:a property.
 */
gdouble
nc_halo_bias_st_spher_get_a (const NcHaloBiasSTSpher *biasf_st_spher)
{
  return biasf_st_spher->a;
}

/**
 * nc_halo_bias_st_spher_set_p:
 * @biasf_st_spher: a #NcHaloBiasSTSpher.
 * @p: value of #NcHaloBiasSTSpher:p.
 *
 * Sets the value @p to the #NcHaloBiasSTSpher:p property.
 *
 */
void
nc_halo_bias_st_spher_set_p (NcHaloBiasSTSpher *biasf_st_spher, gdouble p)
{
  g_assert (p >= 0);
  biasf_st_spher->p = p;
}

/**
 * nc_halo_bias_st_spher_get_p:
 * @biasf_st_spher: a #NcHaloBiasSTSpher.
 *
 * Returns: the value of #NcHaloBiasSTSpher:p property.
 */
gdouble
nc_halo_bias_st_spher_get_p (const NcHaloBiasSTSpher *biasf_st_spher)
{
  return biasf_st_spher->p;
}

// _NC_BIAS_FUNCTION_ST_SPHER_DATASET_9901122 = {1.686, 0.75, 0.3};

static void
nc_halo_bias_st_spher_init (NcHaloBiasSTSpher *biasf_st_spher)
{
  /* TODO: Add initialization code here */
  biasf_st_spher->delta_c = 1.686;
  biasf_st_spher->a = 0.75;
  biasf_st_spher->p = 0.3;
}

static void
_nc_halo_bias_st_spher_finalize (GObject *object)
{
  /* TODO: Add deinitalization code here */

  G_OBJECT_CLASS (nc_halo_bias_st_spher_parent_class)->finalize (object);
}

static void
_nc_halo_bias_st_spher_set_property (GObject * object, guint prop_id, const GValue * value, GParamSpec * pspec)
{
  NcHaloBiasSTSpher *biasf_st_spher = NC_HALO_BIAS_ST_SPHER (object);
  g_return_if_fail (NC_IS_HALO_BIAS_ST_SPHER (object));

  switch (prop_id)
  {
    case PROP_DELTA_C:
      biasf_st_spher->delta_c = g_value_get_double (value);
      break;
	case PROP_A:
      biasf_st_spher->a = g_value_get_double (value);
      break;
	case PROP_P:
      biasf_st_spher->p = g_value_get_double (value);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_halo_bias_st_spher_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcHaloBiasSTSpher *biasf_st_spher = NC_HALO_BIAS_ST_SPHER (object);
  g_return_if_fail (NC_IS_HALO_BIAS_ST_SPHER (object));

  switch (prop_id)
  {
    case PROP_DELTA_C:
      g_value_set_double (value, biasf_st_spher->delta_c);
      break;
	case PROP_A:
      g_value_set_double (value, biasf_st_spher->a);
      break;
	case PROP_P:
      g_value_set_double (value, biasf_st_spher->p);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_halo_bias_st_spher_class_init (NcHaloBiasSTSpherClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcHaloBiasClass* parent_class = NC_HALO_BIAS_CLASS (klass);

  parent_class->eval = &_nc_halo_bias_st_spher_eval;

  object_class->finalize = _nc_halo_bias_st_spher_finalize;
  object_class->set_property = _nc_halo_bias_st_spher_set_property;
  object_class->get_property = _nc_halo_bias_st_spher_get_property;

  /**
   * NcHaloBiasSTSpher:critical_delta:
   *
   * FIXME (check limits values)
   */
  g_object_class_install_property (object_class,
                                   PROP_DELTA_C,
                                   g_param_spec_double ("critical-delta",
                                                        NULL,
                                                        "Critical delta",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 1.686,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY |G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  /**
   * NcHaloBiasSTSpher:a:
   *
   * FIXME (check limits values)
   */
  g_object_class_install_property (object_class,
                                   PROP_A,
                                   g_param_spec_double ("a",
                                                        NULL,
                                                        "a",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 0.75,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY |G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  /**
   * NcHaloBiasSTSpher:p:
   *
   * FIXME (check limits values)
   */
  g_object_class_install_property (object_class,
                                   PROP_P,
                                   g_param_spec_double ("p",
                                                        NULL,
                                                        "p",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 0.75,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY |G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}


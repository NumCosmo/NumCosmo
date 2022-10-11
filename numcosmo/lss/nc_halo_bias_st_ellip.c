/***************************************************************************
 *            nc_halo_bias_st_ellip.c
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
 * SECTION:nc_halo_bias_st_ellip
 * @title: NcHaloBiasSTEllip
 * @short_description: Sheth-Tormen elliptical halo bias function type.
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "lss/nc_halo_bias_st_ellip.h"
#include "math/ncm_cfg.h"
#include "math/ncm_util.h"

G_DEFINE_TYPE (NcHaloBiasSTEllip, nc_halo_bias_st_ellip, NC_TYPE_HALO_BIAS);

enum
{
  PROP_0,
  PROP_DELTA_C,
  PROP_A,
  PROP_B,
  PROP_C
};

/**
 * nc_halo_bias_st_ellip_new:
 * @delta_c: FIXME
 * @a: FIXME
 * @b: FIXME
 * @c: FIXME
 *
 * FIXME
 *
 * Returns: A new #NcHaloBias.
 */
NcHaloBias *
nc_halo_bias_st_ellip_new (gdouble delta_c, gdouble a, gdouble b, gdouble c)
{
  return g_object_new (NC_TYPE_HALO_BIAS_ST_ELLIP,
                       "critical-delta", delta_c,
                       "a", a,
                       "b", b,
                       "c", c,
                       NULL);
}

static gdouble
_nc_halo_bias_st_ellip_eval (NcHaloBias *biasf, gdouble sigma, gdouble z)
{
  NcHaloBiasSTEllip *bias_st_ellip = NC_HALO_BIAS_ST_ELLIP (biasf);
  const gdouble a = bias_st_ellip->a;
  const gdouble b = bias_st_ellip->b;
  const gdouble c = bias_st_ellip->c;
  gdouble x = bias_st_ellip->delta_c / sigma;
  gdouble x2 = x * x;
  gdouble ax2_c = pow(a * x2, c);
  gdouble b_ST_ellip = 1.0  + (a * x2 + b * pow(a * x2, (1.0 - c)) - ax2_c / (sqrt(a) * (ax2_c + b * (1.0 - c) * (1.0 - c/2.0)))) / bias_st_ellip->delta_c;

  NCM_UNUSED (z);
  
  return b_ST_ellip;
}

/**
 * nc_halo_bias_st_ellip_set_delta_c:
 * @biasf_st_ellip: a #NcHaloBiasSTEllip.
 * @delta_c: value of #NcHaloBiasSTEllip:critical-delta.
 *
 * Sets the value @delta_c to the #NcHaloBiasSTEllip:critical-delta property.
 *
 */
void
nc_halo_bias_st_ellip_set_delta_c (NcHaloBiasSTEllip *biasf_st_ellip, gdouble delta_c)
{
  g_assert (delta_c >= 0);
  biasf_st_ellip->delta_c = delta_c;
}

/**
 * nc_halo_bias_st_ellip_get_delta_c:
 * @biasf_st_ellip: a #NcHaloBiasSTEllip.
 *
 * Returns: the value of #NcHaloBiasSTEllip:critical_delta property.
 */
gdouble
nc_halo_bias_st_ellip_get_delta_c (const NcHaloBiasSTEllip *biasf_st_ellip)
{
  return biasf_st_ellip->delta_c;
}

/**
 * nc_halo_bias_st_ellip_set_a:
 * @biasf_st_ellip: a #NcHaloBiasSTEllip.
 * @a: value of #NcHaloBiasSTEllip:a.
 *
 * Sets the value @a to the #NcHaloBiasSTEllip:a property.
 *
 */
void
nc_halo_bias_st_ellip_set_a (NcHaloBiasSTEllip *biasf_st_ellip, gdouble a)
{
  g_assert (a >= 0);
  biasf_st_ellip->a = a;
}

/**
 * nc_halo_bias_st_ellip_get_a:
 * @biasf_st_ellip: a #NcHaloBiasSTEllip.
 *
 * Returns: the value of #NcHaloBiasSTEllip:a property.
 */
gdouble
nc_halo_bias_st_ellip_get_a (const NcHaloBiasSTEllip *biasf_st_ellip)
{
  return biasf_st_ellip->a;
}

/**
 * nc_halo_bias_st_ellip_set_b:
 * @biasf_st_ellip: a #NcHaloBiasSTEllip.
 * @b: value of #NcHaloBiasSTEllip:b.
 *
 * Sets the value @b to the #NcHaloBiasSTEllip:b property.
 *
 */
void
nc_halo_bias_st_ellip_set_b (NcHaloBiasSTEllip *biasf_st_ellip, gdouble b)
{
  g_assert (b >= 0);
  biasf_st_ellip->b = b;
}

/**
 * nc_halo_bias_st_ellip_get_b:
 * @biasf_st_ellip: a #NcHaloBiasSTEllip.
 *
 * Returns: the value of #NcHaloBiasSTEllip:b property.
 */
gdouble
nc_halo_bias_st_ellip_get_b (const NcHaloBiasSTEllip *biasf_st_ellip)
{
  return biasf_st_ellip->b;
}

/**
 * nc_halo_bias_st_ellip_set_c:
 * @biasf_st_ellip: a #NcHaloBiasSTEllip.
 * @c: value of #NcHaloBiasSTEllip:c.
 *
 * Sets the value @c to the #NcHaloBiasSTEllip:c property.
 *
 */
void
nc_halo_bias_st_ellip_set_c (NcHaloBiasSTEllip *biasf_st_ellip, gdouble c)
{
  g_assert (c >= 0);
  biasf_st_ellip->c = c;
}

/**
 * nc_halo_bias_st_ellip_get_c:
 * @biasf_st_ellip: a #NcHaloBiasSTEllip.
 *
 * Returns: the value of #NcHaloBiasSTEllip:c property.
 */
gdouble
nc_halo_bias_st_ellip_get_c (const NcHaloBiasSTEllip *biasf_st_ellip)
{
  return biasf_st_ellip->c;
}

// _NC_BIAS_FUNCTION_ST_ELLIP_DATASET_1001_3162 = {1.686, 0.707, 0.5, 0.6};

static void
nc_halo_bias_st_ellip_init (NcHaloBiasSTEllip *biasf_st_ellip)
{
  /* TODO: Add initialization code here */
  biasf_st_ellip->delta_c = 1.686;
  biasf_st_ellip->a = 0.707;
  biasf_st_ellip->b = 0.5;
  biasf_st_ellip->c = 0.6;
}

static void
_nc_halo_bias_st_ellip_finalize (GObject *object)
{
  /* TODO: Add deinitalization code here */

  G_OBJECT_CLASS (nc_halo_bias_st_ellip_parent_class)->finalize (object);
}

static void
_nc_halo_bias_st_ellip_set_property (GObject * object, guint prop_id, const GValue * value, GParamSpec * pspec)
{
  NcHaloBiasSTEllip *biasf_st_ellip = NC_HALO_BIAS_ST_ELLIP (object);
  g_return_if_fail (NC_IS_HALO_BIAS_ST_ELLIP (object));

  switch (prop_id)
  {
    case PROP_DELTA_C:
      biasf_st_ellip->delta_c = g_value_get_double (value);
      break;
	case PROP_A:
      biasf_st_ellip->a = g_value_get_double (value);
      break;
	case PROP_B:
      biasf_st_ellip->b = g_value_get_double (value);
      break;
	case PROP_C:
      biasf_st_ellip->c = g_value_get_double (value);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_halo_bias_st_ellip_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcHaloBiasSTEllip *biasf_st_ellip = NC_HALO_BIAS_ST_ELLIP (object);
  g_return_if_fail (NC_IS_HALO_BIAS_ST_ELLIP (object));

  switch (prop_id)
  {
    case PROP_DELTA_C:
      g_value_set_double (value, biasf_st_ellip->delta_c);
      break;
	case PROP_A:
      g_value_set_double (value, biasf_st_ellip->a);
      break;
	case PROP_B:
      g_value_set_double (value, biasf_st_ellip->b);
      break;
	case PROP_C:
      g_value_set_double (value, biasf_st_ellip->c);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_halo_bias_st_ellip_class_init (NcHaloBiasSTEllipClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcHaloBiasClass* parent_class = NC_HALO_BIAS_CLASS (klass);

  parent_class->eval = &_nc_halo_bias_st_ellip_eval;

  object_class->finalize = _nc_halo_bias_st_ellip_finalize;
  object_class->set_property = _nc_halo_bias_st_ellip_set_property;
  object_class->get_property = _nc_halo_bias_st_ellip_get_property;

  /**
   * NcHaloBiasSTEllip:critical_delta:
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
   * NcHaloBiasSTEllip:a:
   *
   * FIXME (check limits values)
   */
  g_object_class_install_property (object_class,
                                   PROP_A,
                                   g_param_spec_double ("a",
                                                        NULL,
                                                        "a",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 0.707,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY |G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  /**
   * NcHaloBiasSTEllip:b:
   *
   * FIXME (check limits values)
   */
  g_object_class_install_property (object_class,
                                   PROP_B,
                                   g_param_spec_double ("b",
                                                        NULL,
                                                        "b",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 0.5,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY |G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  /**
   * NcHaloBiasSTEllip:c:
   *
   * FIXME (check limits values)
   */
  g_object_class_install_property (object_class,
                                   PROP_C,
                                   g_param_spec_double ("c",
                                                        NULL,
                                                        "c",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 0.6,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY |G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}


/***************************************************************************
 *            nc_halo_bias_tinker.c
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
 * SECTION:nc_halo_bias_tinker
 * @title: NcHaloBiasTinker
 * @short_description: Tinker halo bias function type.
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "lss/nc_halo_bias_tinker.h"
#include "math/ncm_cfg.h"
#include "math/ncm_util.h"

G_DEFINE_TYPE (NcHaloBiasTinker, nc_halo_bias_tinker, NC_TYPE_HALO_BIAS);

enum
{
  PROP_0,
  PROP_DELTA_C,
  PROP_B0,
  PROP_B1,
  PROP_C,
  PROP_SIZE
};

/**
 * nc_halo_bias_tinker_new:
 * @delta_c: FIXME
 * @B: FIXME
 * @b: FIXME
 * @c: FIXME
 *
 * FIXME
 *
 * Returns: A new #NcHaloBias.
 */
NcHaloBias *
nc_halo_bias_tinker_new (gdouble delta_c, gdouble B, gdouble b, gdouble c)
{
  return g_object_new (NC_TYPE_HALO_BIAS_TINKER,
                       "critical-delta", delta_c,
                       "B", B,
                       "b", b,
                       "c", c,
                       NULL);
}

static gdouble
_nc_halo_bias_tinker_eval (NcHaloBias *biasf, NcHICosmo *cosmo ,gdouble sigma, gdouble z)
{
  NcHaloBiasTinker *bias_tinker = NC_HALO_BIAS_TINKER (biasf);
  NcMultiplicityFunc *mulf = nc_halo_mass_function_peek_multiplicity_function(biasf->mfp);

  const gdouble Delta = nc_multiplicity_func_get_matter_Delta (mulf, cosmo, z);
  const gdouble y = log10(Delta);
  const gdouble u = exp(- pow ( 4.0 / y, 4.0));
  const gdouble A = 1.0 + 0.24 * y * u;
  const gdouble a = 0.44 * y - 0.88;
  const gdouble B = bias_tinker->B;
  const gdouble b = bias_tinker->b;
  const gdouble C = 0.019 + 0.107 * y + 0.19 * u;
  const gdouble c = bias_tinker->c;
  gdouble x = bias_tinker->delta_c / sigma;
  gdouble b_Tinker = 1.0  - A * pow(x, a) / (pow(x, a) + pow(bias_tinker->delta_c, a)) + B * pow(x, b) + C * pow(x, c);
//  printf ("A = %.5g, a=%.5g, B=%.5g, b=%.5g, C=%.5g, c=%.5g, delta_c= %.5g Delta=%.5g\n", A, a, B, b, C, c, bias_tinker->delta_c, log10(x));

  return b_Tinker;
}

/**
 * nc_halo_bias_tinker_set_delta_c:
 * @biasf_tinker: a #NcHaloBiasTinker.
 * @delta_c: value of #NcHaloBiasTinker:critical-delta.
 *
 * Sets the value @delta_c to the #NcHaloBiasTinker:critical-delta property.
 *
 */
void
nc_halo_bias_tinker_set_delta_c (NcHaloBiasTinker *biasf_tinker, gdouble delta_c)
{
  g_assert (delta_c >= 0);
  biasf_tinker->delta_c = delta_c;
}

/**
 * nc_halo_bias_tinker_get_delta_c:
 * @biasf_tinker: a #NcHaloBiasTinker.
 *
 * Returns: the value of #NcHaloBiasTinker:critical_delta property.
 */
gdouble
nc_halo_bias_tinker_get_delta_c (const NcHaloBiasTinker *biasf_tinker)
{
  return biasf_tinker->delta_c;
}

/**
 * nc_halo_bias_tinker_set_B:
 * @biasf_tinker: a #NcHaloBiasTinker.
 * @B: value of #NcHaloBiasTinker:B.
 *
 * Sets the value @B to the #NcHaloBiasTinker:B property.
 *
 */
void
nc_halo_bias_tinker_set_B (NcHaloBiasTinker *biasf_tinker, gdouble B)
{
  g_assert (B >= 0);
  biasf_tinker->B = B;
}

/**
 * nc_halo_bias_tinker_get_B:
 * @biasf_tinker: a #NcHaloBiasTinker.
 *
 * Returns: the value of #NcHaloBiasTinker:B property.
 */
gdouble
nc_halo_bias_tinker_get_B (const NcHaloBiasTinker *biasf_tinker)
{
  return biasf_tinker->B;
}

/**
 * nc_halo_bias_tinker_set_b:
 * @biasf_tinker: a #NcHaloBiasTinker.
 * @b: value of #NcHaloBiasTinker:b.
 *
 * Sets the value @b to the #NcHaloBiasTinker:b property.
 *
 */
void
nc_halo_bias_tinker_set_b (NcHaloBiasTinker *biasf_tinker, gdouble b)
{
  g_assert (b >= 0);
  biasf_tinker->b = b;
}

/**
 * nc_halo_bias_tinker_get_b:
 * @biasf_tinker: a #NcHaloBiasTinker.
 *
 * Returns: the value of #NcHaloBiasTinker:b property.
 */
gdouble
nc_halo_bias_tinker_get_b (const NcHaloBiasTinker *biasf_tinker)
{
  return biasf_tinker->b;
}

/**
 * nc_halo_bias_tinker_set_c:
 * @biasf_tinker: a #NcHaloBiasTinker.
 * @c: value of #NcHaloBiasTinker:c.
 *
 * Sets the value @c to the #NcHaloBiasTinker:c property.
 *
 */
void
nc_halo_bias_tinker_set_c (NcHaloBiasTinker *biasf_tinker, gdouble c)
{
  g_assert (c >= 0);
  biasf_tinker->c = c;
}

/**
 * nc_halo_bias_tinker_get_c:
 * @biasf_tinker: a #NcHaloBiasTinker.
 *
 * Returns: the value of #NcHaloBiasTinker:c property.
 */
gdouble
nc_halo_bias_tinker_get_c (const NcHaloBiasTinker *biasf_tinker)
{
  return biasf_tinker->c;
}


// _NC_BIAS_FUNCTION_TINKER_DATASET_1001_3162_DELTA = {1.686, 0.183, 1.5, 2.4, 200.0};

static void
nc_halo_bias_tinker_init (NcHaloBiasTinker *biasf_tinker)
{
  /* TODO: Add initialization code here */
  biasf_tinker->delta_c = 1.686;
  biasf_tinker->B = 0.183;
  biasf_tinker->b = 1.5;
  biasf_tinker->c = 2.4;
}

static void
_nc_halo_bias_tinker_finalize (GObject *object)
{
  /* TODO: Add deinitalization code here */

  G_OBJECT_CLASS (nc_halo_bias_tinker_parent_class)->finalize (object);
}

static void
_nc_halo_bias_tinker_set_property (GObject * object, guint prop_id, const GValue * value, GParamSpec * pspec)
{
  NcHaloBiasTinker *biasf_tinker = NC_HALO_BIAS_TINKER (object);
  g_return_if_fail (NC_IS_HALO_BIAS_TINKER (object));

  switch (prop_id)
  {
    case PROP_DELTA_C:
      biasf_tinker->delta_c = g_value_get_double (value);
      break;
	case PROP_B0:
      biasf_tinker->B = g_value_get_double (value);
      break;
	case PROP_B1:
      biasf_tinker->b = g_value_get_double (value);
      break;
	case PROP_C:
      biasf_tinker->c = g_value_get_double (value);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_halo_bias_tinker_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcHaloBiasTinker *biasf_tinker = NC_HALO_BIAS_TINKER (object);
  g_return_if_fail (NC_IS_HALO_BIAS_TINKER (object));

  switch (prop_id)
  {
    case PROP_DELTA_C:
      g_value_set_double (value, biasf_tinker->delta_c);
      break;
	case PROP_B0:
      g_value_set_double (value, biasf_tinker->B);
      break;
	case PROP_B1:
      g_value_set_double (value, biasf_tinker->b);
      break;
	case PROP_C:
      g_value_set_double (value, biasf_tinker->c);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_halo_bias_tinker_class_init (NcHaloBiasTinkerClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcHaloBiasClass* parent_class = NC_HALO_BIAS_CLASS (klass);

  parent_class->eval = &_nc_halo_bias_tinker_eval;

  object_class->finalize = _nc_halo_bias_tinker_finalize;
  object_class->set_property = _nc_halo_bias_tinker_set_property;
  object_class->get_property = _nc_halo_bias_tinker_get_property;

  /**
   * NcHaloBiasTinker:critical_delta:
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
   * NcHaloBiasTinker:B:
   *
   * FIXME (check limits values)
   */
  g_object_class_install_property (object_class,
                                   PROP_B0,
                                   g_param_spec_double ("B",
                                                        NULL,
                                                        "B",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 0.183,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY |G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  /**
   * NcHaloBiasTinker:b:
   *
   * FIXME (check limits values)
   */
  g_object_class_install_property (object_class,
                                   PROP_B1,
                                   g_param_spec_double ("b",
                                                        NULL,
                                                        "b",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 1.5,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY |G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  /**
   * NcHaloBiasTinker:c:
   *
   * FIXME (check limits values)
   */
  g_object_class_install_property (object_class,
                                   PROP_C,
                                   g_param_spec_double ("c",
                                                        NULL,
                                                        "c",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 2.4,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY |G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

}


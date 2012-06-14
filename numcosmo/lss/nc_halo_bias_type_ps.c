/***************************************************************************
 *            nc_halo_bias_type_ps.c
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
 * SECTION:nc_halo_bias_type_ps
 * @title: PS Halo Bias Function Type
 * @short_description: Press-Schechter FIXME
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
#include <glib.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_const_mksa.h>
#include <glib.h>

G_DEFINE_TYPE (NcHaloBiasTypePS, nc_halo_bias_type_ps, NC_TYPE_HALO_BIAS_TYPE);

enum
{
  PROP_0,
  PROP_DELTA_C
};

/**
 * nc_halo_bias_type_ps_new:
 * @delta_c: FIXME
 *
 * FIXME
 *
 * Returns: A new #NcHaloBiasType.
 */
NcHaloBiasType *
nc_halo_bias_type_ps_new (gdouble delta_c)
{
  return g_object_new (NC_TYPE_HALO_BIAS_TYPE_PS,
                       "critical-delta", delta_c,
                       NULL);
}

static gdouble
_nc_halo_bias_type_ps_eval (NcHaloBiasType*biasf, gdouble sigma, gdouble z)
{
  NcHaloBiasTypePS *bias_ps = NC_HALO_BIAS_TYPE_PS (biasf);
  gdouble x = bias_ps->delta_c / sigma;        /* \delta_c \sigma^{-1} */
  gdouble x2 = x * x;
  gdouble b_PS = 1.0 + (x2 - 1.0) / bias_ps->delta_c;

  return b_PS;
}

/**
 * nc_halo_bias_type_ps_set_delta_c:
 * @biasf_ps: a #NcHaloBiasTypePS.
 * @delta_c: value of #NcHaloBiasTypePS:critical-delta.
 *
 * Sets the value @delta_c to the #NcHaloBiasTypePS:critical-delta property.
 *
 */
void
nc_halo_bias_type_ps_set_delta_c (NcHaloBiasTypePS *biasf_ps, gdouble delta_c)
{
  g_assert (delta_c >= 0);
  biasf_ps->delta_c = delta_c;
}

/**
 * nc_halo_bias_type_ps_get_delta_c:
 * @biasf_ps: a #NcHaloBiasTypePS.
 *
 * Returns: the value of #NcHaloBiasTypePS:critical_delta property.
 */
gdouble
nc_halo_bias_type_ps_get_delta_c (const NcHaloBiasTypePS *biasf_ps)
{
  return biasf_ps->delta_c;
}

static void
nc_halo_bias_type_ps_init (NcHaloBiasTypePS *biasf_ps)
{
  /* TODO: Add initialization code here */
  biasf_ps->delta_c = 1.686;
}

static void
_nc_halo_bias_type_ps_finalize (GObject *object)
{
  /* TODO: Add deinitalization code here */

  G_OBJECT_CLASS (nc_halo_bias_type_ps_parent_class)->finalize (object);
}

static void
_nc_halo_bias_type_ps_set_property (GObject * object, guint prop_id, const GValue * value, GParamSpec * pspec)
{
  NcHaloBiasTypePS *biasf_ps = NC_HALO_BIAS_TYPE_PS (object);
  g_return_if_fail (NC_IS_HALO_BIAS_TYPE_PS (object));

  switch (prop_id)
  {
    case PROP_DELTA_C:
      biasf_ps->delta_c = g_value_get_double (value);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_halo_bias_type_ps_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcHaloBiasTypePS *biasf_ps = NC_HALO_BIAS_TYPE_PS (object);
  g_return_if_fail (NC_IS_HALO_BIAS_TYPE_PS (object));

  switch (prop_id)
  {
    case PROP_DELTA_C:
      g_value_set_double (value, biasf_ps->delta_c);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_halo_bias_type_ps_class_init (NcHaloBiasTypePSClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcHaloBiasTypeClass* parent_class = NC_HALO_BIAS_TYPE_CLASS (klass);

  parent_class->eval = &_nc_halo_bias_type_ps_eval;

  object_class->finalize = _nc_halo_bias_type_ps_finalize;
  object_class->set_property = _nc_halo_bias_type_ps_set_property;
  object_class->get_property = _nc_halo_bias_type_ps_get_property;

  /**
   * NcHaloBiasTypePS:critical_delta:
   *
   * FIXME
   */
  g_object_class_install_property (object_class,
                                   PROP_DELTA_C,
                                   g_param_spec_double ("critical-delta",
                                                        NULL,
                                                        "Critical delta",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 1.686,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY |G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}


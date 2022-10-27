/***************************************************************************
 *            nc_halo_bias_ps.c
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
 * SECTION:nc_halo_bias_ps
 * @title: NcHaloBiasPS
 * @short_description: Press-Schechter halo bias function type.
 *
 * FIXME
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "lss/nc_halo_bias_ps.h"
#include "math/ncm_cfg.h"
#include "math/ncm_util.h"

G_DEFINE_TYPE (NcHaloBiasPS, nc_halo_bias_ps, NC_TYPE_HALO_BIAS);

enum
{
  PROP_0,
  PROP_DELTA_C
};

static void
nc_halo_bias_ps_init (NcHaloBiasPS *biasf_ps)
{
  biasf_ps->delta_c = 0.0;
}

static void
_nc_halo_bias_ps_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_halo_bias_ps_parent_class)->finalize (object);
}

static void
_nc_halo_bias_ps_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcHaloBiasPS *biasf_ps = NC_HALO_BIAS_PS (object);

  g_return_if_fail (NC_IS_HALO_BIAS_PS (object));

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
_nc_halo_bias_ps_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcHaloBiasPS *biasf_ps = NC_HALO_BIAS_PS (object);

  g_return_if_fail (NC_IS_HALO_BIAS_PS (object));

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

static gdouble _nc_halo_bias_ps_eval (NcHaloBias *biasf, NcHICosmo *cosmo, gdouble sigma, gdouble z);

static void
nc_halo_bias_ps_class_init (NcHaloBiasPSClass *klass)
{
  GObjectClass *object_class    = G_OBJECT_CLASS (klass);
  NcHaloBiasClass *parent_class = NC_HALO_BIAS_CLASS (klass);

  object_class->finalize     = _nc_halo_bias_ps_finalize;
  object_class->set_property = _nc_halo_bias_ps_set_property;
  object_class->get_property = _nc_halo_bias_ps_get_property;

  /**
   * NcHaloBiasPS:critical_delta:
   *
   * FIXME
   */
  g_object_class_install_property (object_class,
                                   PROP_DELTA_C,
                                   g_param_spec_double ("critical-delta",
                                                        NULL,
                                                        "Critical delta",
                                                        0.0, G_MAXDOUBLE, 1.686,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  parent_class->eval = &_nc_halo_bias_ps_eval;
}

static gdouble
_nc_halo_bias_ps_eval (NcHaloBias *biasf, NcHICosmo *cosmo, gdouble sigma, gdouble z)
{
  NcHaloBiasPS *bias_ps = NC_HALO_BIAS_PS (biasf);
  gdouble x             = bias_ps->delta_c / sigma; /* \delta_c \sigma^{-1} */
  gdouble x2            = x * x;
  gdouble b_PS          = 1.0 + (x2 - 1.0) / bias_ps->delta_c;

  NCM_UNUSED (z);

  return b_PS;
}

/**
 * nc_halo_bias_ps_new:  (constructor)
 * @mfp: a #NcHaloMassFunction
 *
 * FIXME
 *
 * Returns: A new #NcHaloBias.
 */
NcHaloBiasPS *
nc_halo_bias_ps_new (NcHaloMassFunction *mfp)
{
  return g_object_new (NC_TYPE_HALO_BIAS_PS,
                       "mass-function", mfp,
                       NULL);
}

/**
 * nc_halo_bias_ps_new_full: (constructor)
 * @mfp: a #NcHaloMassFunction
 * @delta_c: FIXME
 *
 * FIXME
 *
 * Returns: A new #NcHaloBias.
 */
NcHaloBiasPS *
nc_halo_bias_ps_new_full (NcHaloMassFunction *mfp, gdouble delta_c)
{
  return g_object_new (NC_TYPE_HALO_BIAS_PS,
                       "mass-function", mfp,
                       "critical-delta", delta_c,
                       NULL);
}

/**
 * nc_halo_bias_ps_set_delta_c:
 * @biasf_ps: a #NcHaloBiasPS.
 * @delta_c: value of #NcHaloBiasPS:critical-delta.
 *
 * Sets the value @delta_c to the #NcHaloBiasPS:critical-delta property.
 *
 */
void
nc_halo_bias_ps_set_delta_c (NcHaloBiasPS *biasf_ps, gdouble delta_c)
{
  g_assert (delta_c >= 0);
  biasf_ps->delta_c = delta_c;
}

/**
 * nc_halo_bias_ps_get_delta_c:
 * @biasf_ps: a #NcHaloBiasPS.
 *
 * Returns: the value of #NcHaloBiasPS:critical_delta property.
 */
gdouble
nc_halo_bias_ps_get_delta_c (const NcHaloBiasPS *biasf_ps)
{
  return biasf_ps->delta_c;
}


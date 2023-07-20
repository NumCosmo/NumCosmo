/***************************************************************************
 *            ncm_integralnd.c
 *
 *  Thu July 20 08:39:30 2023
 *  Copyright  2023 Eduardo José Barroso
 *  <eduardo.jsbarroso@uel.br>
 ****************************************************************************/
/*
 * ncm_integralnd.c
 * Copyright (C) 2023 Eduardo José Barroso <eduardo.jsbarroso@uel.br>
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
 * SECTION:ncm_integralnd
 * @title: NcmIntegralnd
 * @short_description: One dimensional integration object.
 * @stability: Stable
 * @include: numcosmo/math/ncm_integralnd.h
 *
 * FIXME
 * 
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_integralnd.h"
#include "math/ncm_c.h"
#include "math/ncm_cfg.h"

#include "misc/cubature.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#endif /* NUMCOSMO_GIR_SCAN */

struct _NcmIntegralndPrivate
{
  guint partition;
  gdouble reltol;
  gdouble abstol;
  guint rule;
  gsl_integration_workspace *ws;
  gsl_integration_cquad_workspace *cquad_ws;
};

G_DEFINE_ABSTRACT_TYPE_WITH_PRIVATE (NcmIntegralnd, ncm_integralnd, G_TYPE_OBJECT);

enum
{
  PROP_0,
  PROP_INTEGRAND,
  PROP_PARTITION,
  PROP_RULE,
  PROP_RELTOL,
  PROP_ABSTOL,
  PROP_SIZE,
};

static void
ncm_integralnd_init (NcmIntegralnd *intnd)
{
  NcmIntegralndPrivate * const self = intnd->priv = ncm_integralnd_get_instance_private (intnd);
  self->partition = 0;
  self->rule      = 0;
  self->reltol    = 0.0;
  self->abstol    = 0.0;
  self->ws        = NULL;
  self->cquad_ws  = NULL;
}

static void
ncm_integralnd_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmIntegralnd *intnd = NCM_INTEGRALND (object);
  g_return_if_fail (NCM_IS_INTEGRALND (object));

  switch (prop_id)
  {
    case PROP_PARTITION:
      ncm_integralnd_set_partition (intnd, g_value_get_uint (value));
      break;
    case PROP_RULE:
      ncm_integralnd_set_rule (intnd, g_value_get_uint (value));
      break;
    case PROP_RELTOL:
      ncm_integralnd_set_reltol (intnd, g_value_get_double (value));
      break;
    case PROP_ABSTOL:
      ncm_integralnd_set_abstol (intnd, g_value_get_double (value));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
ncm_integralnd_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmIntegralnd *intnd = NCM_INTEGRALND (object);
  g_return_if_fail (NCM_IS_INTEGRALND (object));

  switch (prop_id)
  {
    case PROP_PARTITION:
      g_value_set_uint (value, ncm_integralnd_get_partition (intnd));
      break;
    case PROP_RULE:
      g_value_set_uint (value, ncm_integralnd_get_rule (intnd));
      break;
    case PROP_RELTOL:
      g_value_set_double (value, ncm_integralnd_get_reltol (intnd));
      break;
    case PROP_ABSTOL:
      g_value_set_double (value, ncm_integralnd_get_abstol (intnd));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
ncm_integralnd_finalize (GObject *object)
{
  NcmIntegralnd *intnd = NCM_INTEGRALND (object);
  NcmIntegralndPrivate * const self = intnd->priv;

  g_clear_pointer (&self->ws,       gsl_integration_workspace_free);
  g_clear_pointer (&self->cquad_ws, gsl_integration_cquad_workspace_free);
  
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_integralnd_parent_class)->finalize (object);
}

static void
ncm_integralnd_class_init (NcmIntegralndClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = &ncm_integralnd_set_property;
  object_class->get_property = &ncm_integralnd_get_property;
  object_class->finalize     = &ncm_integralnd_finalize;

  g_object_class_install_property (object_class,
                                   PROP_PARTITION,
                                   g_param_spec_uint ("partition",
                                                      NULL,
                                                      "Integral maximum partititon",
                                                      10, G_MAXUINT32, NCM_INTEGRALND_DEFAULT_PARTITION,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT| G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));  
  g_object_class_install_property (object_class,
                                   PROP_RULE,
                                   g_param_spec_uint ("rule",
                                                      NULL,
                                                      "Integration rule",
                                                      1, 6, NCM_INTEGRALND_DEFAULT_ALG,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT| G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));  
  g_object_class_install_property (object_class,
                                   PROP_RELTOL,
                                   g_param_spec_double ("reltol",
                                                        NULL,
                                                        "Integral relative tolerance",
                                                        0.0, 1.0, NCM_INTEGRALND_DEFAULT_RELTOL,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT| G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));  
  g_object_class_install_property (object_class,
                                   PROP_ABSTOL,
                                   g_param_spec_double ("abstol",
                                                        NULL,
                                                        "Integral absolute tolerance",
                                                        0.0, G_MAXDOUBLE, NCM_INTEGRALND_DEFAULT_ABSTOL,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT| G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

/**
 * ncm_integralnd_ref:
 * @intnd: a #NcmIntegralnd
 * 
 * Increases the reference count of @intnd by one.
 * 
 * Returns: (transfer full): @intnd.
 */
NcmIntegralnd *
ncm_integralnd_ref (NcmIntegralnd *intnd)
{
  return g_object_ref (intnd);
}

/**
 * ncm_integralnd_free:
 * @intnd: a #NcmIntegralnd
 * 
 * Decreases the reference count of @intnd by one.
 * 
 */
void 
ncm_integralnd_free (NcmIntegralnd *intnd)
{
  g_object_unref (intnd);
}

/**
 * ncm_integralnd_clear:
 * @intnd: a #NcmIntegralnd
 * 
 * If *@intnd is different from NULL, decreases the reference 
 * count of *@intnd by one and sets *@intnd to NULL.
 * 
 */
void 
ncm_integralnd_clear (NcmIntegralnd **intnd)
{
  g_clear_object (intnd);
}

/**
 * ncm_integralnd_set_partition:
 * @intnd: a #NcmIntegralnd
 * @partition: max number of subintervals
 * 
 * Sets the max number of subintervals to @partition.
 * 
 */
void 
ncm_integralnd_set_partition (NcmIntegralnd *intnd, guint partition)
{
  NcmIntegralndPrivate * const self = intnd->priv;

  g_assert_cmpuint (partition, >=, 10);
  if (self->partition != partition)
  {
    g_clear_pointer (&self->ws, gsl_integration_workspace_free);

    self->ws        = gsl_integration_workspace_alloc (partition);
    self->partition = partition;
  }
}

/**
 * ncm_integralnd_set_rule:
 * @intnd: a #NcmIntegralnd
 * @rule: Gauss-Kronrod rule
 * 
 * Sets the Gauss-Kronrod @rule to use.
 * 
 */
void 
ncm_integralnd_set_rule (NcmIntegralnd *intnd, guint rule)
{
  NcmIntegralndPrivate * const self = intnd->priv;

  self->rule = rule;
}

/**
 * ncm_integralnd_set_reltol:
 * @intnd: a #NcmIntegralnd
 * @reltol: relative tolerance
 * 
 * Sets the relative tolerance @reltol to use.
 * 
 */
void 
ncm_integralnd_set_reltol (NcmIntegralnd *intnd, gdouble reltol)
{
  NcmIntegralndPrivate * const self = intnd->priv;

  self->reltol = reltol;
}

/**
 * ncm_integralnd_set_abstol:
 * @intnd: a #NcmIntegralnd
 * @abstol: absolute tolerance
 * 
 * Sets the absolute tolerance @reltol to use.
 * 
 */
void 
ncm_integralnd_set_abstol (NcmIntegralnd *intnd, gdouble abstol)
{
  NcmIntegralndPrivate * const self = intnd->priv;

  self->abstol = abstol;
}

/**
 * ncm_integralnd_get_partition:
 * @intnd: a #NcmIntegralnd
 * 
 * Returns: the maximum number of subdivisions used.
 */
guint 
ncm_integralnd_get_partition (NcmIntegralnd *intnd)
{
  NcmIntegralndPrivate * const self = intnd->priv;

  return self->partition;
}

/**
 * ncm_integralnd_get_rule:
 * @intnd: a #NcmIntegralnd
 * 
 * Returns: the Gauss-Kronrod rule used.
 */
guint 
ncm_integralnd_get_rule (NcmIntegralnd *intnd)
{
  NcmIntegralndPrivate * const self = intnd->priv;

  return self->rule;
}

/**
 * ncm_integralnd_get_reltol:
 * @intnd: a #NcmIntegralnd
 * 
 * Returns: the relative tolerance used.
 */
gdouble 
ncm_integralnd_get_reltol (NcmIntegralnd *intnd)
{
  NcmIntegralndPrivate * const self = intnd->priv;

  return self->reltol;
}

/**
 * ncm_integralnd_get_abstol:
 * @intnd: a #NcmIntegralnd
 * 
 * Returns: the absolute tolerance used.
 */
gdouble 
ncm_integralnd_get_abstol (NcmIntegralnd *intnd)
{
  NcmIntegralndPrivate * const self = intnd->priv;

  return self->abstol;
}

/**
 * ncm_integralnd_eval_h:
 * @intnd: a #NcmIntegralnd
 * @xi: inferior integration limit $x_i$
 * @xf: superior integration limit $x_f$
 * @err: (out): the error in the integration
 * 
 * Evaluated the integral $I_F(x_i, x_f) = \int_{x_i}^{x_f}F(x)\mathrm{d}x$.
 * 
 * Returns: the value of the integral $I_F(x_i, x_f)$.
 */
gdouble 
ncm_integralnd_eval_h (NcmIntegralnd *intnd, const gdouble xi, const gdouble xf, gdouble *err)
{
  NcmIntegralndPrivate * const self = intnd->priv;
  gdouble result = 0.0;
  gsl_function F;
  gint ret;
  
  F.function = &ncm_integralnd_integrand;
  F.params   = &intnd;
  result = 1.0;
  
  return result;
}


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
 * This object is used to perform n-dimensional integration of a function
 * using different methods.
 *
 * The integration can be performed using the cubature library. The cubature
 * library is a library for adaptive multidimensional integration.
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
  gdouble reltol;
  gdouble abstol;
};

G_DEFINE_ABSTRACT_TYPE_WITH_PRIVATE (NcmIntegralnd, ncm_integralnd, G_TYPE_OBJECT);

enum
{
  PROP_0,
  PROP_RELTOL,
  PROP_ABSTOL,
  PROP_SIZE,
};

static void
ncm_integralnd_init (NcmIntegralnd *intnd)
{
  NcmIntegralndPrivate * const self = intnd->priv = ncm_integralnd_get_instance_private (intnd);

  self->reltol = 0.0;
  self->abstol = 0.0;
}

static void
ncm_integralnd_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmIntegralnd *intnd = NCM_INTEGRALND (object);

  g_return_if_fail (NCM_IS_INTEGRALND (object));

  switch (prop_id)
  {
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
  NcmIntegralnd *intnd              = NCM_INTEGRALND (object);
  NcmIntegralndPrivate * const self = intnd->priv;

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_integralnd_parent_class)->finalize (object);
}

static void
ncm_integralnd_class_init (NcmIntegralndClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = &ncm_integralnd_set_property;
  object_class->get_property = &ncm_integralnd_get_property;
  object_class->finalize     = &ncm_integralnd_finalize;

  g_object_class_install_property (object_class,
                                   PROP_RELTOL,
                                   g_param_spec_double ("reltol",
                                                        NULL,
                                                        "Integral relative tolerance",
                                                        0.0, 1.0, NCM_INTEGRALND_DEFAULT_RELTOL,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_ABSTOL,
                                   g_param_spec_double ("abstol",
                                                        NULL,
                                                        "Integral absolute tolerance",
                                                        0.0, G_MAXDOUBLE, NCM_INTEGRALND_DEFAULT_ABSTOL,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
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
 * ncm_integralnd_eval:
 * @intnd: a #NcmIntegralnd
 * @xi: (array) (element-type double): inferior integration limit $x_i$
 * @xf: (array) (element-type double): superior integration limit $x_f$
 * @err: (array) (element-type double) (out callee-allocates): the error in the integration
 *
 * Evaluated the integral $I_F(x_i, x_f) = \int_{x_i}^{x_f}F(x)\mathrm{d}x$.
 *
 * Returns: (transfer full) (array) (element-type double): the value of the integral $I_F(x_i, x_f)$.
 */
gdouble *
ncm_integralnd_eval (NcmIntegralnd *intnd, const gdouble *xi, const gdouble *xf, gdouble **err)
{
  NcmIntegralndPrivate * const self = intnd->priv;
  gdouble *result                   = g_new (gdouble, 1);

  printf ("xi = %g, xf = %g\n", xi[0], xf[0]);
  fflush (stdout);

  *err = g_new (gdouble, 1);

  return result;
}


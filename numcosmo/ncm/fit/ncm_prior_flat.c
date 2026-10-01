/***************************************************************************
 *            ncm_prior_flat.c
 *
 *  Wed August 03 16:58:19 2016
 *  Copyright  2016  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_prior_flat.c
 * Copyright (C) 2016 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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
 * NcmPriorFlat:
 *
 * Base class for flat prior distributions.
 *
 * This object subclasses #NcmPrior and is the base class of the flat priors used by
 * #NcmLikelihood, on a parameter (#NcmPriorFlatParam) or on a derived quantity
 * (#NcmPriorFlatFunc). The subclass provides the quantity $x$; the prior returns the
 * least-squares form
 * $$
 * f = \exp\left[\frac{h_0}{s}\left(x_0 - x\right) + \frac{h_0}{2}\right] +
 *     \exp\left[\frac{h_0}{s}\left(x - x_1\right) + \frac{h_0}{2}\right],
 * \qquad -2\ln P(x) = f^2,
 * $$
 * where $x_0$ and $x_1$ are the lower and upper limits, $s$ is the width of the walls
 * and $h_0$ their height (#NcmPriorFlat:h0, default 20). Near $x_0$,
 * $-2\ln P \simeq \exp[(2 h_0/s)(x_0 - x + s/2)]$: it is $e^{h_0}$ at $x_0$, one at
 * $x_0 + s/2$ and $e^{-h_0}$ at $x_0 + s$, and the same holds mirrored around $x_1$.
 * Inside $[x_0 + s, x_1 - s]$ the prior adds less than $e^{-h_0}$ to $-2\ln L$.
 *
 * The prior is not normalized. It replaces a hard cut, which an analysis sensitive to
 * discontinuities cannot use, by walls that are smooth but steep.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/fit/ncm_prior_flat.h"

enum
{
  PROP_0,
  PROP_X_LOW,
  PROP_X_UPP,
  PROP_S,
  PROP_H0,
  PROP_VARIABLE,
};


typedef struct _NcmPriorFlatPrivate
{
  /*< private >*/
  gdouble x_low;
  gdouble x_upp;
  gdouble s;
  gdouble var;
  gdouble h0;
} NcmPriorFlatPrivate;


G_DEFINE_TYPE_WITH_PRIVATE (NcmPriorFlat, ncm_prior_flat, NCM_TYPE_PRIOR)

static void
ncm_prior_flat_init (NcmPriorFlat *pf)
{
  NcmPriorFlatPrivate * const self = ncm_prior_flat_get_instance_private (pf);

  self->x_low = 0.0;
  self->x_upp = 0.0;
  self->s     = 0.0;
  self->h0    = 0.0;
}

static void
_ncm_prior_flat_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmPriorFlat *pf = NCM_PRIOR_FLAT (object);

  g_return_if_fail (NCM_IS_PRIOR_FLAT (object));

  switch (prop_id)
  {
    case PROP_X_LOW:
      ncm_prior_flat_set_x_low (pf, g_value_get_double (value));
      break;
    case PROP_X_UPP:
      ncm_prior_flat_set_x_upp (pf, g_value_get_double (value));
      break;
    case PROP_S:
      ncm_prior_flat_set_scale (pf, g_value_get_double (value));
      break;
    case PROP_VARIABLE:
      ncm_prior_flat_set_var (pf, g_value_get_double (value));
      break;
    case PROP_H0:
      ncm_prior_flat_set_h0 (pf, g_value_get_double (value));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_prior_flat_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmPriorFlat *pf = NCM_PRIOR_FLAT (object);

  g_return_if_fail (NCM_IS_PRIOR_FLAT (object));

  switch (prop_id)
  {
    case PROP_X_LOW:
      g_value_set_double (value, ncm_prior_flat_get_x_low (pf));
      break;
    case PROP_X_UPP:
      g_value_set_double (value, ncm_prior_flat_get_x_upp (pf));
      break;
    case PROP_S:
      g_value_set_double (value, ncm_prior_flat_get_scale (pf));
      break;
    case PROP_VARIABLE:
      g_value_set_double (value, ncm_prior_flat_get_var (pf));
      break;
    case PROP_H0:
      g_value_set_double (value, ncm_prior_flat_get_h0 (pf));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void _ncm_prior_flat_eval (NcmMSetFunc *func, NcmMSet *mset, const gdouble *x, gdouble *res);

static gdouble
_ncm_prior_flat_mean (NcmPriorFlat *pf, NcmMSet *mset)
{
  g_error ("method mean not implemented by %s.", G_OBJECT_TYPE_NAME (pf));

  return 0.0;
}

static void
ncm_prior_flat_class_init (NcmPriorFlatClass *klass)
{
  GObjectClass *object_class        = G_OBJECT_CLASS (klass);
  NcmMSetFuncClass *mset_func_class = NCM_MSET_FUNC_CLASS (klass);

  object_class->set_property = &_ncm_prior_flat_set_property;
  object_class->get_property = &_ncm_prior_flat_get_property;

  /**
   * NcmPriorFlat:x-low:
   *
   * The lower limit $x_0$. Default: 0.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_X_LOW,
                                   g_param_spec_double ("x-low",
                                                        NULL,
                                                        "Lower limit",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmPriorFlat:x-upp:
   *
   * The upper limit $x_1$. Default: 1.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_X_UPP,
                                   g_param_spec_double ("x-upp",
                                                        NULL,
                                                        "Upper limit",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 1.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmPriorFlat:scale:
   *
   * The width $s$ of the walls. Default: $10^{-10}$.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_S,
                                   g_param_spec_double ("scale",
                                                        NULL,
                                                        "Width of the walls",
                                                        G_MINDOUBLE, G_MAXDOUBLE, 1.0e-10,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmPriorFlat:h0:
   *
   * The height $h_0$ of the walls: $-2\ln P = e^{h_0}$ at the limits. Default: 20.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_H0,
                                   g_param_spec_double ("h0",
                                                        NULL,
                                                        "Height of the walls",
                                                        1.0, G_MAXDOUBLE, 20.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmPriorFlat:variable:
   *
   * The argument passed to the mean function of #NcmPriorFlatFunc; the other
   * subclasses do not read it. Default: 0.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_VARIABLE,
                                   g_param_spec_double ("variable",
                                                        NULL,
                                                        "Argument of the mean function",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  NCM_PRIOR_CLASS (klass)->is_m2lnL = FALSE;
  mset_func_class->eval             = &_ncm_prior_flat_eval;
  klass->mean                       = &_ncm_prior_flat_mean;
}

static void
_ncm_prior_flat_eval (NcmMSetFunc *func, NcmMSet *mset, const gdouble *x, gdouble *res)
{
  NcmPriorFlat *pf                 = NCM_PRIOR_FLAT (func);
  NcmPriorFlatPrivate * const self = ncm_prior_flat_get_instance_private (pf);
  const gdouble mean               = NCM_PRIOR_FLAT_GET_CLASS (pf)->mean (pf, mset);

  /* f, with f^2 = -2 ln P, see the class description. */
  res[0] = exp (self->h0 / self->s * (self->x_low - mean) + 0.5 * self->h0) +
           exp (self->h0 / self->s * (mean - self->x_upp) + 0.5 * self->h0);
}

/**
 * ncm_prior_flat_ref:
 * @pf: a #NcmPriorFlat
 *
 * Increases the reference count of @pf atomically.
 *
 * Returns: (transfer full): @pf.
 */
NcmPriorFlat *
ncm_prior_flat_ref (NcmPriorFlat *pf)
{
  return g_object_ref (pf);
}

/**
 * ncm_prior_flat_free:
 * @pf: a #NcmPriorFlat
 *
 * Decreases the reference count of @pf atomically.
 *
 */
void
ncm_prior_flat_free (NcmPriorFlat *pf)
{
  g_object_unref (pf);
}

/**
 * ncm_prior_flat_clear:
 * @pf: a #NcmPriorFlat
 *
 * Decreases the reference count of *@pf and sets *@pf to NULL.
 *
 */
void
ncm_prior_flat_clear (NcmPriorFlat **pf)
{
  g_clear_object (pf);
}

/**
 * ncm_prior_flat_set_x_low:
 * @pf: a #NcmPriorFlat
 * @x_low: lower limit
 *
 * Sets the lower limit of @pf.
 *
 */
void
ncm_prior_flat_set_x_low (NcmPriorFlat *pf, const gdouble x_low)
{
  NcmPriorFlatPrivate * const self = ncm_prior_flat_get_instance_private (pf);

  self->x_low = x_low;
}

/**
 * ncm_prior_flat_set_x_upp:
 * @pf: a #NcmPriorFlat
 * @x_upp: upper limit
 *
 * Sets the upper limit of @pf.
 *
 */
void
ncm_prior_flat_set_x_upp (NcmPriorFlat *pf, const gdouble x_upp)
{
  NcmPriorFlatPrivate * const self = ncm_prior_flat_get_instance_private (pf);

  self->x_upp = x_upp;
}

/**
 * ncm_prior_flat_set_scale:
 * @pf: a #NcmPriorFlat
 * @scale: width of the walls
 *
 * Sets #NcmPriorFlat:scale.
 *
 */
void
ncm_prior_flat_set_scale (NcmPriorFlat *pf, const gdouble scale)
{
  NcmPriorFlatPrivate * const self = ncm_prior_flat_get_instance_private (pf);

  self->s = scale;
}

/**
 * ncm_prior_flat_set_var:
 * @pf: a #NcmPriorFlat
 * @var: argument of the mean function
 *
 * Sets #NcmPriorFlat:variable.
 *
 */
void
ncm_prior_flat_set_var (NcmPriorFlat *pf, const gdouble var)
{
  NcmPriorFlatPrivate * const self = ncm_prior_flat_get_instance_private (pf);

  self->var = var;
}

/**
 * ncm_prior_flat_set_h0:
 * @pf: a #NcmPriorFlat
 * @h0: height of the walls
 *
 * Sets #NcmPriorFlat:h0.
 *
 */
void
ncm_prior_flat_set_h0 (NcmPriorFlat *pf, const gdouble h0)
{
  NcmPriorFlatPrivate * const self = ncm_prior_flat_get_instance_private (pf);

  self->h0 = h0;
}

/**
 * ncm_prior_flat_get_x_low:
 * @pf: a #NcmPriorFlat
 *
 * Returns: the lower limit of @pf.
 */
gdouble
ncm_prior_flat_get_x_low (NcmPriorFlat *pf)
{
  NcmPriorFlatPrivate * const self = ncm_prior_flat_get_instance_private (pf);

  return self->x_low;
}

/**
 * ncm_prior_flat_get_x_upp:
 * @pf: a #NcmPriorFlat
 *
 * Returns: the upper limit of @pf.
 */
gdouble
ncm_prior_flat_get_x_upp (NcmPriorFlat *pf)
{
  NcmPriorFlatPrivate * const self = ncm_prior_flat_get_instance_private (pf);

  return self->x_upp;
}

/**
 * ncm_prior_flat_get_scale:
 * @pf: a #NcmPriorFlat
 *
 * Returns: #NcmPriorFlat:scale.
 */
gdouble
ncm_prior_flat_get_scale (NcmPriorFlat *pf)
{
  NcmPriorFlatPrivate * const self = ncm_prior_flat_get_instance_private (pf);

  return self->s;
}

/**
 * ncm_prior_flat_get_var:
 * @pf: a #NcmPriorFlat
 *
 * Returns: #NcmPriorFlat:variable.
 */
gdouble
ncm_prior_flat_get_var (NcmPriorFlat *pf)
{
  NcmPriorFlatPrivate * const self = ncm_prior_flat_get_instance_private (pf);

  return self->var;
}

/**
 * ncm_prior_flat_get_h0:
 * @pf: a #NcmPriorFlat
 *
 * Returns: #NcmPriorFlat:h0.
 */
gdouble
ncm_prior_flat_get_h0 (NcmPriorFlat *pf)
{
  NcmPriorFlatPrivate * const self = ncm_prior_flat_get_instance_private (pf);

  return self->h0;
}


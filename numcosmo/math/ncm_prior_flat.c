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
 * This object subclasses #NcmPrior and defines a base class for flat priors used by
 * NcmLikelihood. These objects describe flat prior distributions applicable to
 * parameters or any derived quantity.
 *
 * The prior is defined as:
 * $$
 * -2\ln P(x) = \exp\left(\frac{2 h_0}{s} \left(\left(x_0 - x\right) + \frac{s}{2.0}\right)\right) +
 *              \exp\left(\frac{2 h_0}{s} \left(\left(x - x_1\right) + \frac{s}{2.0}\right)\right),
 * $$
 * where $x_0$ and $x_1$ are the lower and upper limits, respectively, and $s$ is the
 * scale of the prior. The variable $h_0$ has a default value of $20.0$. Note that for
 * $x = x_0$, the first term is $e^{h_0}$ and grows exponentially with $x$ decreasing.
 * For $x = x_0 + s$, the first term is $e^{-h_0}$ and decreases exponentially with $x$
 * increasing. The same happens for the second term around $x_1$.
 *
 * The prior is not normalized. It is useful for defining a flat prior when the
 * analysis is sensitive to discontinuities in the prior.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_prior_flat.h"

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
  NcmPrior parent_instance;
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

static void
_ncm_prior_flat_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_prior_flat_parent_class)->finalize (object);
}

static void _ncm_prior_flat_eval (NcmMSetFunc *func, NcmMSet *mset, const gdouble *x, gdouble *res);

static void
ncm_prior_flat_class_init (NcmPriorFlatClass *klass)
{
  GObjectClass *object_class        = G_OBJECT_CLASS (klass);
  NcmMSetFuncClass *mset_func_class = NCM_MSET_FUNC_CLASS (klass);

  object_class->set_property = &_ncm_prior_flat_set_property;
  object_class->get_property = &_ncm_prior_flat_get_property;
  object_class->finalize     = &_ncm_prior_flat_finalize;

  g_object_class_install_property (object_class,
                                   PROP_X_LOW,
                                   g_param_spec_double ("x-low",
                                                        NULL,
                                                        "lower limit",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_X_UPP,
                                   g_param_spec_double ("x-upp",
                                                        NULL,
                                                        "upper limit",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 1.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_S,
                                   g_param_spec_double ("scale",
                                                        NULL,
                                                        "border scale",
                                                        G_MINDOUBLE, G_MAXDOUBLE, 1.0e-10,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_H0,
                                   g_param_spec_double ("h0",
                                                        NULL,
                                                        "Cut magnitude",
                                                        1.0, G_MAXDOUBLE, 20.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_VARIABLE,
                                   g_param_spec_double ("variable",
                                                        NULL,
                                                        "variable",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  NCM_PRIOR_CLASS (klass)->is_m2lnL = FALSE;
  mset_func_class->eval             = &_ncm_prior_flat_eval;
}

static void
_ncm_prior_flat_eval (NcmMSetFunc *func, NcmMSet *mset, const gdouble *x, gdouble *res)
{
  NcmPriorFlat *pf                 = NCM_PRIOR_FLAT (func);
  NcmPriorFlatPrivate * const self = ncm_prior_flat_get_instance_private (pf);
  const gdouble mean               = NCM_PRIOR_FLAT_GET_CLASS (pf)->mean (pf, mset);

  res[0] = 0.5 * (exp (-self->h0 / self->s * (mean - self->x_low) + self->h0) +
                  exp (+self->h0 / self->s * (mean - self->x_upp) + self->h0));
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
 * @scale: border scale
 *
 * Sets the border scale of @pf.
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
 * @var: variable
 *
 * Sets the variable of @pf.
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
 * @h0: Cut magnitude
 *
 * Sets the cut magnitude of @pf.
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
 * Returns: the border scale of @pf.
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
 * Returns: the variable of @pf.
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
 * Returns: the cut magnitude of @pf.
 */
gdouble
ncm_prior_flat_get_h0 (NcmPriorFlat *pf)
{
  NcmPriorFlatPrivate * const self = ncm_prior_flat_get_instance_private (pf);

  return self->h0;
}


/***************************************************************************
 *            ncm_prior_gauss.c
 *
 *  Wed August 03 10:19:46 2016
 *  Copyright  2016  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_prior_gauss.c
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
 * NcmPriorGauss:
 *
 * A gaussian prior for NcmLikelihood.
 *
 * This object is a subclass of #NcmPrior, serving as a base class for Gaussian priors
 * used by NcmLikelihood. These objects describe Gaussian prior distributions
 * applicable to parameters or any derived quantity.
 *
 * The Gaussian prior is defined as:
 * $$
 * -2\ln P(x) = \frac{\left(x - \mu\right)^2}{\sigma},
 * $$
 * where $\mu$ is the mean and $\sigma$ is the standard deviation. This Gaussian prior
 * places higher probability density around the mean and decreases exponentially as the
 * parameter deviates from the mean.
 *
 * The prior is not normalized. It is particularly useful for defining a Gaussian prior
 * when the analysis benefits from a smooth and symmetric distribution around the mean.
 * Additionally, it is compatible with least-squares based analysis.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_prior_gauss.h"

enum
{
  PROP_0,
  PROP_MU,
  PROP_SIGMA,
  PROP_VARIABLE,
};

typedef struct _NcmPriorGaussPrivate
{
  /*< private >*/
  NcmPrior parent_instance;
  gdouble mu;
  gdouble sigma;
  gdouble var;
} NcmPriorGaussPrivate;

G_DEFINE_TYPE_WITH_PRIVATE (NcmPriorGauss, ncm_prior_gauss, NCM_TYPE_PRIOR)

static void
ncm_prior_gauss_init (NcmPriorGauss *pg)
{
  NcmPriorGaussPrivate * const self = ncm_prior_gauss_get_instance_private (pg);

  self->mu    = 0.0;
  self->sigma = 0.0;
  self->var   = 0.0;
}

static void
_ncm_prior_gauss_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmPriorGauss *pg = NCM_PRIOR_GAUSS (object);

  g_return_if_fail (NCM_IS_PRIOR_GAUSS (object));

  switch (prop_id)
  {
    case PROP_MU:
      ncm_prior_gauss_set_mean (pg, g_value_get_double (value));
      break;
    case PROP_SIGMA:
      ncm_prior_gauss_set_sigma (pg, g_value_get_double (value));
      break;
    case PROP_VARIABLE:
      ncm_prior_gauss_set_var (pg, g_value_get_double (value));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_prior_gauss_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmPriorGauss *pg = NCM_PRIOR_GAUSS (object);

  g_return_if_fail (NCM_IS_PRIOR_GAUSS (object));

  switch (prop_id)
  {
    case PROP_MU:
      g_value_set_double (value, ncm_prior_gauss_get_mean (pg));
      break;
    case PROP_SIGMA:
      g_value_set_double (value, ncm_prior_gauss_get_sigma (pg));
      break;
    case PROP_VARIABLE:
      g_value_set_double (value, ncm_prior_gauss_get_var (pg));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_prior_gauss_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_prior_gauss_parent_class)->finalize (object);
}

static void _ncm_prior_gauss_eval (NcmMSetFunc *func, NcmMSet *mset, const gdouble *x, gdouble *res);

static void
ncm_prior_gauss_class_init (NcmPriorGaussClass *klass)
{
  GObjectClass *object_class        = G_OBJECT_CLASS (klass);
  NcmMSetFuncClass *mset_func_class = NCM_MSET_FUNC_CLASS (klass);

  object_class->set_property = &_ncm_prior_gauss_set_property;
  object_class->get_property = &_ncm_prior_gauss_get_property;
  object_class->finalize     = &_ncm_prior_gauss_finalize;

  g_object_class_install_property (object_class,
                                   PROP_MU,
                                   g_param_spec_double ("mu",
                                                        NULL,
                                                        "mean",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_SIGMA,
                                   g_param_spec_double ("sigma",
                                                        NULL,
                                                        "standard deviation",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 1.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_VARIABLE,
                                   g_param_spec_double ("variable",
                                                        NULL,
                                                        "variable",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  NCM_PRIOR_CLASS (klass)->is_m2lnL = FALSE;
  mset_func_class->eval             = &_ncm_prior_gauss_eval;
}

static void
_ncm_prior_gauss_eval (NcmMSetFunc *func, NcmMSet *mset, const gdouble *x, gdouble *res)
{
  NcmPriorGauss *pg                 = NCM_PRIOR_GAUSS (func);
  NcmPriorGaussPrivate * const self = ncm_prior_gauss_get_instance_private (pg);
  const gdouble mean                = NCM_PRIOR_GAUSS_GET_CLASS (pg)->mean (pg, mset);

  res[0] = (mean - self->mu) / self->sigma;
}

/**
 * ncm_prior_gauss_ref:
 * @pg: a #NcmPriorGauss
 *
 * Increases the reference count of @pg atomically.
 *
 * Returns: (transfer full): @pg.
 */
NcmPriorGauss *
ncm_prior_gauss_ref (NcmPriorGauss *pg)
{
  return g_object_ref (pg);
}

/**
 * ncm_prior_gauss_free:
 * @pg: a #NcmPriorGauss
 *
 * Decreases the reference count of @pg atomically.
 *
 */
void
ncm_prior_gauss_free (NcmPriorGauss *pg)
{
  g_object_unref (pg);
}

/**
 * ncm_prior_gauss_clear:
 * @pg: a #NcmPriorGauss
 *
 * Decreases the reference count of *@pg and sets *@pg to NULL.
 *
 */
void
ncm_prior_gauss_clear (NcmPriorGauss **pg)
{
  g_clear_object (pg);
}

/**
 * ncm_prior_gauss_set_mean:
 * @pg: a #NcmPriorGauss
 * @mean: mean
 *
 * Sets the mean of @pg.
 *
 */
void
ncm_prior_gauss_set_mean (NcmPriorGauss *pg, const gdouble mean)
{
  NcmPriorGaussPrivate * const self = ncm_prior_gauss_get_instance_private (pg);

  self->mu = mean;
}

/**
 * ncm_prior_gauss_set_sigma:
 * @pg: a #NcmPriorGauss
 * @sigma: standard deviation
 *
 * Sets the standard deviation of @pg.
 *
 */
void
ncm_prior_gauss_set_sigma (NcmPriorGauss *pg, const gdouble sigma)
{
  NcmPriorGaussPrivate * const self = ncm_prior_gauss_get_instance_private (pg);

  g_assert_cmpfloat (sigma, >, 0.0);
  self->sigma = sigma;
}

/**
 * ncm_prior_gauss_set_var:
 * @pg: a #NcmPriorGauss
 * @var: variable
 *
 * Sets the variable of @pg.
 *
 */
void
ncm_prior_gauss_set_var (NcmPriorGauss *pg, const gdouble var)
{
  NcmPriorGaussPrivate * const self = ncm_prior_gauss_get_instance_private (pg);

  self->var = var;
}

/**
 * ncm_prior_gauss_get_mean:
 * @pg: a #NcmPriorGauss
 *
 * Returns: the mean of @pg.
 */
gdouble
ncm_prior_gauss_get_mean (NcmPriorGauss *pg)
{
  NcmPriorGaussPrivate * const self = ncm_prior_gauss_get_instance_private (pg);

  return self->mu;
}

/**
 * ncm_prior_gauss_get_sigma:
 * @pg: a #NcmPriorGauss
 *
 * Returns: the standard deviation of @pg.
 */
gdouble
ncm_prior_gauss_get_sigma (NcmPriorGauss *pg)
{
  NcmPriorGaussPrivate * const self = ncm_prior_gauss_get_instance_private (pg);

  return self->sigma;
}

/**
 * ncm_prior_gauss_get_var:
 * @pg: a #NcmPriorGauss
 *
 * Returns: the variable of @pg.
 */
gdouble
ncm_prior_gauss_get_var (NcmPriorGauss *pg)
{
  NcmPriorGaussPrivate * const self = ncm_prior_gauss_get_instance_private (pg);

  return self->var;
}


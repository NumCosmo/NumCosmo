/***************************************************************************
 *            ncm_mset_trans_kern.c
 *
 *  Fri August 29 18:57:07 2014
 *  Copyright  2014  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_mset_trans_kern.c
 * Copyright (C) 2014 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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
 * NcmMSetTransKern:
 *
 * Abstract proposal distribution over the free parameters of a #NcmMSet.
 *
 * A kernel draws a point $\theta^\star$ given the current point $\theta$
 * (ncm_mset_trans_kern_generate()) and evaluates the density $q(\theta^\star|\theta)$
 * (ncm_mset_trans_kern_pdf()). #NcmFitMCMC uses it for its proposals. With a fixed
 * center (ncm_mset_trans_kern_set_prior()) it serves as a prior sampler, e.g. for the
 * initial points of #NcmFitESMCMC (ncm_mset_trans_kern_prior_sample()).
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/fit/ncm_mset_trans_kern.h"

enum
{
  PROP_0,
  PROP_MSET,
  PROP_SIZE
};

typedef struct _NcmMSetTransKernPrivate
{
  /*< private >*/
  NcmMSet *mset;
  NcmVector *theta;
} NcmMSetTransKernPrivate;


G_DEFINE_ABSTRACT_TYPE_WITH_PRIVATE (NcmMSetTransKern, ncm_mset_trans_kern, G_TYPE_OBJECT)

static void
ncm_mset_trans_kern_init (NcmMSetTransKern *tkern)
{
  NcmMSetTransKernPrivate *self = ncm_mset_trans_kern_get_instance_private (tkern);

  self->mset  = NULL;
  self->theta = NULL;
}

static void
ncm_mset_trans_kern_dispose (GObject *object)
{
  NcmMSetTransKern *tkern       = NCM_MSET_TRANS_KERN (object);
  NcmMSetTransKernPrivate *self = ncm_mset_trans_kern_get_instance_private (tkern);

  ncm_mset_clear (&self->mset);
  ncm_vector_clear (&self->theta);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_mset_trans_kern_parent_class)->dispose (object);
}

static void
ncm_mset_trans_kern_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_mset_trans_kern_parent_class)->finalize (object);
}

static void
ncm_mset_trans_kern_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmMSetTransKern *tkern = NCM_MSET_TRANS_KERN (object);

  g_return_if_fail (NCM_IS_MSET_TRANS_KERN (object));

  switch (prop_id)
  {
    case PROP_MSET:
      ncm_mset_trans_kern_set_mset (tkern, g_value_get_object (value));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
ncm_mset_trans_kern_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmMSetTransKern *tkern       = NCM_MSET_TRANS_KERN (object);
  NcmMSetTransKernPrivate *self = ncm_mset_trans_kern_get_instance_private (tkern);

  g_return_if_fail (NCM_IS_MSET_TRANS_KERN (object));

  switch (prop_id)
  {
    case PROP_MSET:
      g_value_set_object (value, self->mset);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void _ncm_mset_trans_kern_reset (NcmMSetTransKern *tkern);
static void _ncm_mset_trans_kern_set_mset (NcmMSetTransKern *tkern, NcmMSet *mset);
static void _ncm_mset_trans_kern_generate (NcmMSetTransKern *tkern, NcmVector *theta, NcmVector *thetastar, NcmRNG *rng);
static gdouble _ncm_mset_trans_kern_pdf (NcmMSetTransKern *tkern, NcmVector *theta, NcmVector *thetastar);
static const gchar *_ncm_mset_trans_kern_get_name (NcmMSetTransKern *tkern);

static void
ncm_mset_trans_kern_class_init (NcmMSetTransKernClass *klass)
{
  GObjectClass *object_class         = G_OBJECT_CLASS (klass);
  NcmMSetTransKernClass *tkern_class = NCM_MSET_TRANS_KERN_CLASS (klass);

  object_class->set_property = ncm_mset_trans_kern_set_property;
  object_class->get_property = ncm_mset_trans_kern_get_property;

  object_class->dispose  = ncm_mset_trans_kern_dispose;
  object_class->finalize = ncm_mset_trans_kern_finalize;

  g_object_class_install_property (object_class,
                                   PROP_MSET,
                                   g_param_spec_object ("mset",
                                                        NULL,
                                                        "NcmMSet",
                                                        NCM_TYPE_MSET,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  tkern_class->set_mset = &_ncm_mset_trans_kern_set_mset;
  tkern_class->generate = &_ncm_mset_trans_kern_generate;
  tkern_class->pdf      = &_ncm_mset_trans_kern_pdf;
  tkern_class->reset    = &_ncm_mset_trans_kern_reset;
  tkern_class->get_name = &_ncm_mset_trans_kern_get_name;
}

static void
_ncm_mset_trans_kern_reset (NcmMSetTransKern *tkern)
{
}

/* LCOV_EXCL_START */

static void
_ncm_mset_trans_kern_set_mset (NcmMSetTransKern *tkern, NcmMSet *mset)
{
  g_error ("method set_mset not implemented by %s.", G_OBJECT_TYPE_NAME (tkern));
}

static void
_ncm_mset_trans_kern_generate (NcmMSetTransKern *tkern, NcmVector *theta, NcmVector *thetastar, NcmRNG *rng)
{
  g_error ("method generate not implemented by %s.", G_OBJECT_TYPE_NAME (tkern));
}

static gdouble
_ncm_mset_trans_kern_pdf (NcmMSetTransKern *tkern, NcmVector *theta, NcmVector *thetastar)
{
  g_error ("method pdf not implemented by %s.", G_OBJECT_TYPE_NAME (tkern));

  return 0.0;
}

static const gchar *
_ncm_mset_trans_kern_get_name (NcmMSetTransKern *tkern)
{
  g_error ("method get_name not implemented by %s.", G_OBJECT_TYPE_NAME (tkern));

  return NULL;
}

/* LCOV_EXCL_STOP */

/**
 * ncm_mset_trans_kern_ref:
 * @tkern: a #NcmMSetTransKern.
 *
 * Increases the reference count of @tkern.
 *
 * Returns: (transfer full): @tkern.
 */
NcmMSetTransKern *
ncm_mset_trans_kern_ref (NcmMSetTransKern *tkern)
{
  return g_object_ref (tkern);
}

/**
 * ncm_mset_trans_kern_free:
 * @tkern: a #NcmMSetTransKern.
 *
 * Decreases the reference count of @tkern.
 *
 */
void
ncm_mset_trans_kern_free (NcmMSetTransKern *tkern)
{
  g_object_unref (tkern);
}

/**
 * ncm_mset_trans_kern_clear:
 * @tkern: a #NcmMSetTransKern.
 *
 * If *@tkern is not %NULL, unrefs it and sets *@tkern to %NULL.
 *
 */
void
ncm_mset_trans_kern_clear (NcmMSetTransKern **tkern)
{
  g_clear_object (tkern);
}

/**
 * ncm_mset_trans_kern_set_mset: (virtual set_mset)
 * @tkern: a #NcmMSetTransKern.
 * @mset: a #NcmMSet.
 *
 * Sets the @mset as the internal set #NcmMSet to be used by the transition kernel.
 *
 */
void
ncm_mset_trans_kern_set_mset (NcmMSetTransKern *tkern, NcmMSet *mset)
{
  NcmMSetTransKernPrivate *self = ncm_mset_trans_kern_get_instance_private (tkern);

  ncm_mset_clear (&self->mset);

  if (mset != NULL)
  {
    g_assert (ncm_mset_fparam_map_valid (mset));

    if (ncm_mset_fparam_len (mset) == 0)
      g_error ("ncm_mset_trans_kern_set_mset: invalid mset, no free parameters.");

    self->mset = ncm_mset_ref (mset);
    NCM_MSET_TRANS_KERN_GET_CLASS (tkern)->set_mset (tkern, mset);
  }
}

/**
 * ncm_mset_trans_kern_peek_mset:
 * @tkern: a #NcmMSetTransKern.
 *
 * Returns: (transfer none): the internal set #NcmMSet.
 */
NcmMSet *
ncm_mset_trans_kern_peek_mset (NcmMSetTransKern *tkern)
{
  NcmMSetTransKernPrivate *self = ncm_mset_trans_kern_get_instance_private (tkern);

  return self->mset;
}

/**
 * ncm_mset_trans_kern_set_prior:
 * @tkern: a #NcmMSetTransKern.
 * @theta: a #NcmVector of free-parameter values
 *
 * Makes @theta the fixed center of the kernel, so that ncm_mset_trans_kern_prior_sample()
 * and ncm_mset_trans_kern_prior_pdf() use it as a prior.
 *
 */
void
ncm_mset_trans_kern_set_prior (NcmMSetTransKern *tkern, NcmVector *theta)
{
  NcmMSetTransKernPrivate *self = ncm_mset_trans_kern_get_instance_private (tkern);

  ncm_vector_clear (&self->theta);
  self->theta = ncm_vector_ref (theta);
}

/**
 * ncm_mset_trans_kern_set_prior_from_mset:
 * @tkern: a #NcmMSetTransKern.
 *
 * As ncm_mset_trans_kern_set_prior() but uses the values present in the
 * internal set #NcmMSet.
 *
 */
void
ncm_mset_trans_kern_set_prior_from_mset (NcmMSetTransKern *tkern)
{
  NcmMSetTransKernPrivate *self = ncm_mset_trans_kern_get_instance_private (tkern);

  g_assert (self->mset != NULL);
  {
    guint fparams_len = ncm_mset_fparams_len (self->mset);
    NcmVector *theta  = ncm_vector_new (fparams_len);

    ncm_mset_fparams_get_vector (self->mset, theta);
    ncm_mset_trans_kern_set_prior (tkern, theta);
    ncm_vector_free (theta);
  }
}

/**
 * ncm_mset_trans_kern_generate: (virtual generate)
 * @tkern: a #NcmMSetTransKern.
 * @theta: current point.
 * @thetastar: try point.
 * @rng: a #NcmRNG.
 *
 * Generates a new point @thetastar from @theta using the transition kernel.
 *
 */
void
ncm_mset_trans_kern_generate (NcmMSetTransKern *tkern, NcmVector *theta, NcmVector *thetastar, NcmRNG *rng)
{
  NCM_MSET_TRANS_KERN_GET_CLASS (tkern)->generate (tkern, theta, thetastar, rng);
}

/**
 * ncm_mset_trans_kern_pdf: (virtual pdf)
 * @tkern: a #NcmMSetTransKern.
 * @theta: current point.
 * @thetastar: try point.
 *
 * Computes the density $q(\theta^\star|\theta)$ of the kernel.
 *
 * Returns: the density of @thetastar given @theta
 */
gdouble
ncm_mset_trans_kern_pdf (NcmMSetTransKern *tkern, NcmVector *theta, NcmVector *thetastar)
{
  return NCM_MSET_TRANS_KERN_GET_CLASS (tkern)->pdf (tkern, theta, thetastar);
}

#define NCM_MSET_TRANS_KERN_PRIOR_MAX_ITER 1000

/**
 * ncm_mset_trans_kern_prior_sample:
 * @tkern: a #NcmMSetTransKern.
 * @thetastar: try point.
 * @rng: a #NcmRNG.
 *
 * Sample from the transition kernel using it as a prior. To use as a prior one must
 * call one of the functions ncm_mset_trans_kern_set_prior_* first. Draws outside the
 * parameter bounds are drawn again; 1000 failed draws abort.
 *
 */
void
ncm_mset_trans_kern_prior_sample (NcmMSetTransKern *tkern, NcmVector *thetastar, NcmRNG *rng)
{
  NcmMSetTransKernPrivate *self = ncm_mset_trans_kern_get_instance_private (tkern);
  NcmMSet *mset                 = ncm_mset_trans_kern_peek_mset (tkern);
  guint iter, i;

  g_assert (self->theta != NULL);

  for (iter = 0; iter < NCM_MSET_TRANS_KERN_PRIOR_MAX_ITER; iter++)
  {
    ncm_mset_trans_kern_generate (tkern, self->theta, thetastar, rng);

    if (ncm_mset_fparam_valid_bounds (mset, thetastar))
      return;
  }

  for (i = 0; i < ncm_mset_fparam_len (mset); i++)
  {
    const gdouble lb  = ncm_mset_fparam_get_lower_bound (mset, i);
    const gdouble ub  = ncm_mset_fparam_get_upper_bound (mset, i);
    const gdouble val = ncm_vector_get (thetastar, i);

    if ((val < lb) || (val > ub))
      g_warning ("ncm_mset_trans_kern_prior_sample: parameter %u (%s) is out of bounds [%.16g, %.16g]: %.16g",
                 i, ncm_mset_fparam_name (mset, i), lb, ub, val);
  }

  g_error ("ncm_mset_trans_kern_prior_sample: failed to draw a sample within the bounds after %u draws.",
           NCM_MSET_TRANS_KERN_PRIOR_MAX_ITER);
}

/**
 * ncm_mset_trans_kern_prior_pdf:
 * @tkern: a #NcmMSetTransKern.
 * @thetastar: try point.
 *
 * Computes the density of the kernel at @thetastar about the center set by
 * ncm_mset_trans_kern_set_prior() or ncm_mset_trans_kern_set_prior_from_mset(),
 * which must be called first.
 *
 * Returns: the density at @thetastar
 */
gdouble
ncm_mset_trans_kern_prior_pdf (NcmMSetTransKern *tkern, NcmVector *thetastar)
{
  NcmMSetTransKernPrivate *self = ncm_mset_trans_kern_get_instance_private (tkern);

  g_assert (self->theta != NULL);

  return NCM_MSET_TRANS_KERN_GET_CLASS (tkern)->pdf (tkern, self->theta, thetastar);
}

/**
 * ncm_mset_trans_kern_reset: (virtual reset)
 * @tkern: a #NcmMSetTransKern.
 *
 * Resets the transition kernel.
 *
 */
void
ncm_mset_trans_kern_reset (NcmMSetTransKern *tkern)
{
  NCM_MSET_TRANS_KERN_GET_CLASS (tkern)->reset (tkern);
}

/**
 * ncm_mset_trans_kern_get_name: (virtual get_name)
 * @tkern: a #NcmMSetTransKern.
 *
 * Returns: the name of the sampler.
 *
 */
const gchar *
ncm_mset_trans_kern_get_name (NcmMSetTransKern *tkern)
{
  return NCM_MSET_TRANS_KERN_GET_CLASS (tkern)->get_name (tkern);
}


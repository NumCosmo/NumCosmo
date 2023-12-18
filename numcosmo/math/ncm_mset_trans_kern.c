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
 * SECTION:ncm_mset_trans_kern
 * @title: NcmMSetTransKern
 * @short_description: Abstract Class for a transition kernel and prior.
 *
 * This object defines the abstract class for a transition kernel and prior. It serves
 * as the base class for all transition kernels and priors, with two main purposes:
 *
 * - To define the interface for all transition kernels for use in the NcmFitMCMC
 *   object.
 * - To define the interface for all priors, generating random parameter vectors with
 *   multivariate parameters.
 *
 * Notably, it acts as a prior sampler for NcmFitESMCMC, generating the initial
 * population's first set of random parameter vectors.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_mset_trans_kern.h"

enum
{
  PROP_0,
  PROP_MSET,
  PROP_SIZE
};

typedef struct _NcmMSetTransKernPrivate
{
  /*< private >*/
  GObject parent_instance;
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
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
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
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void _ncm_mset_trans_kern_reset (NcmMSetTransKern *tkern);

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

  tkern_class->bernoulli_scheme = FALSE;
  tkern_class->set_mset         = NULL;
  tkern_class->generate         = NULL;
  tkern_class->reset            = &_ncm_mset_trans_kern_reset;
}

static void
_ncm_mset_trans_kern_reset (NcmMSetTransKern *tkern)
{
}

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
 * Increases the reference count of @tkern.
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
 * @theta: a #NcmMSet.
 *
 * Sets the @theta as the prior mean. This allows the transition kernel to
 * be used as a prior sampler.
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
 * Computes the value of the kernel at (@theta, @thetastar).
 *
 * Returns: the value of the kernel at (@theta, @thetastar).
 */
gdouble
ncm_mset_trans_kern_pdf (NcmMSetTransKern *tkern, NcmVector *theta, NcmVector *thetastar)
{
  return NCM_MSET_TRANS_KERN_GET_CLASS (tkern)->pdf (tkern, theta, thetastar);
}

/**
 * ncm_mset_trans_kern_prior_sample:
 * @tkern: a #NcmMSetTransKern.
 * @thetastar: try point.
 * @rng: a #NcmRNG.
 *
 * Sample from the transition kernel using it as a prior. To use as a prior one must
 * call one of the functions ncm_mset_trans_kern_set_prior_* first.
 *
 */
void
ncm_mset_trans_kern_prior_sample (NcmMSetTransKern *tkern, NcmVector *thetastar, NcmRNG *rng)
{
  NcmMSetTransKernPrivate *self = ncm_mset_trans_kern_get_instance_private (tkern);

  g_assert (self->theta != NULL);
  ncm_mset_trans_kern_generate (tkern, self->theta, thetastar, rng);
}

/**
 * ncm_mset_trans_kern_prior_pdf:
 * @tkern: a #NcmMSetTransKern.
 * @thetastar: try point.
 *
 * FIXME
 *
 * Returns: the value of the kernel at (@ktern->theta, @thetastar).
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


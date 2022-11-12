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
 * FIXME
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

G_DEFINE_ABSTRACT_TYPE (NcmMSetTransKern, ncm_mset_trans_kern, G_TYPE_OBJECT);

static void
ncm_mset_trans_kern_init (NcmMSetTransKern *tkern)
{
  tkern->mset  = NULL;
  tkern->theta = NULL;
}

static void
ncm_mset_trans_kern_dispose (GObject *object)
{
  NcmMSetTransKern *tkern = NCM_MSET_TRANS_KERN (object);

  ncm_mset_clear (&tkern->mset);
  ncm_vector_clear (&tkern->theta);
  
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
  NcmMSetTransKern *tkern = NCM_MSET_TRANS_KERN (object);
  g_return_if_fail (NCM_IS_MSET_TRANS_KERN (object));

  switch (prop_id)
  {
    case PROP_MSET:
      g_value_set_object (value, tkern->mset);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
ncm_mset_trans_kern_class_init (NcmMSetTransKernClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);
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
  tkern_class->set_mset = NULL;
  tkern_class->generate = NULL;
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
 * FIXME
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
 * FIXME
 * 
 */
void
ncm_mset_trans_kern_set_mset (NcmMSetTransKern *tkern, NcmMSet *mset)
{
  ncm_mset_clear (&tkern->mset);
  if (mset != NULL)
  {
    g_assert (mset->valid_map);
    if (ncm_mset_fparam_len (mset) == 0)
      g_error ("ncm_mset_trans_kern_set_mset: invalid mset, no free parameters.");
    
    tkern->mset = ncm_mset_ref (mset);
    NCM_MSET_TRANS_KERN_GET_CLASS (tkern)->set_mset (tkern, mset);
  }
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
  ncm_vector_clear (&tkern->theta);
  tkern->theta = ncm_vector_ref (theta);
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
  g_assert (tkern->mset != NULL);
  {
    guint fparams_len = ncm_mset_fparams_len (tkern->mset);
    NcmVector *theta = ncm_vector_new (fparams_len);
    ncm_mset_fparams_get_vector (tkern->mset, theta);
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
 * FIXME
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
 * FIXME
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
  g_assert (tkern->theta != NULL);
  ncm_mset_trans_kern_generate (tkern, tkern->theta, thetastar, rng);
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
  g_assert (tkern->theta != NULL);
  return NCM_MSET_TRANS_KERN_GET_CLASS (tkern)->pdf (tkern, tkern->theta, thetastar);
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

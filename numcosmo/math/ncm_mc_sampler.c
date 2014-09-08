/***************************************************************************
 *            ncm_mc_sampler.c
 *
 *  Fri August 29 18:57:07 2014
 *  Copyright  2014  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * ncm_mc_sampler.c
 * Copyright (C) 2014 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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
 * SECTION:ncm_mc_sampler
 * @title: Markov Chain Sampler
 * @short_description: Object implementing a generic chain sampler.
 *
 * FIXME
 * 
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_mc_sampler.h"

enum
{
   PROP_0,
   PROP_MSET,
   PROP_SIZE
};

G_DEFINE_ABSTRACT_TYPE (NcmMCSampler, ncm_mc_sampler, G_TYPE_OBJECT);

static void
ncm_mc_sampler_init (NcmMCSampler *ncm_mc_sampler)
{
}

static void
ncm_mc_sampler_dispose (GObject *object)
{
  NcmMCSampler *mcs = NCM_MC_SAMPLER (object);

  ncm_mset_clear (&mcs->mset);
  
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_mc_sampler_parent_class)->dispose (object);
}

static void
ncm_mc_sampler_finalize (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_mc_sampler_parent_class)->finalize (object);
}

static void
ncm_mc_sampler_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmMCSampler *mcs = NCM_MC_SAMPLER (object);
  g_return_if_fail (NCM_IS_MC_SAMPLER (object));

  switch (prop_id)
  {
    case PROP_MSET:
      ncm_mc_sampler_set_mset (mcs, g_value_get_object (value));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
ncm_mc_sampler_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmMCSampler *mcs = NCM_MC_SAMPLER (object);
  g_return_if_fail (NCM_IS_MC_SAMPLER (object));

  switch (prop_id)
  {
    case PROP_MSET:
      g_value_set_object (value, mcs->mset);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
ncm_mc_sampler_class_init (NcmMCSamplerClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);
  NcmMCSamplerClass *mcs_class = NCM_MC_SAMPLER_CLASS (klass);

  object_class->set_property = ncm_mc_sampler_set_property;
  object_class->get_property = ncm_mc_sampler_get_property;

  object_class->dispose  = ncm_mc_sampler_dispose;
  object_class->finalize = ncm_mc_sampler_finalize;

  g_object_class_install_property (object_class,
                                   PROP_MSET,
                                   g_param_spec_object ("mset",
                                                        NULL,
                                                        "NcmMSet",
                                                        NCM_TYPE_MSET,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  mcs_class->bernoulli_scheme = FALSE;
  mcs_class->set_mset = NULL;
  mcs_class->generate = NULL;
}

/**
 * ncm_mc_sampler_ref:
 * @mcs: a #NcmMCSampler.
 *
 * Increases the reference count of @mcs.
 *
 * Returns: (transfer full): @mcs.
 */
NcmMCSampler *
ncm_mc_sampler_ref (NcmMCSampler *mcs)
{
  return g_object_ref (mcs);
}

/**
 * ncm_mc_sampler_free:
 * @mcs: a #NcmMCSampler.
 *
 * Increases the reference count of @mcs.
 *
 */
void 
ncm_mc_sampler_free (NcmMCSampler *mcs)
{
  g_object_unref (mcs);
}

/**
 * ncm_mc_sampler_clear:
 * @mcs: a #NcmMCSampler.
 *
 * FIXME
 *
 */
void 
ncm_mc_sampler_clear (NcmMCSampler **mcs)
{
  g_clear_object (mcs);
}

/**
 * ncm_mc_sampler_set_mset: (virtual set_mset)
 * @mcs: a #NcmMCSampler.
 * @mset: a #NcmMSet.
 *
 * FIXME
 * 
 */
void
ncm_mc_sampler_set_mset (NcmMCSampler *mcs, NcmMSet *mset)
{
  ncm_mset_clear (&mcs->mset);
  if (mset != NULL)
  {
    g_assert (mset->valid_map);
    mcs->mset = ncm_mset_ref (mset);
    NCM_MC_SAMPLER_GET_CLASS (mcs)->set_mset (mcs, mset);
  }
}

/**
 * ncm_mc_sampler_generate: (virtual generate)
 * @mcs: a #NcmMCSampler.
 * @theta: current point.
 * @thetastar: try point.
 * @rng: a #NcmRNG.
 *
 * FIXME
 * 
 */
void
ncm_mc_sampler_generate (NcmMCSampler *mcs, NcmVector *theta, NcmVector *thetastar, NcmRNG *rng)
{
  NCM_MC_SAMPLER_GET_CLASS (mcs)->generate (mcs, theta, thetastar, rng);
}

/**
 * ncm_mc_sampler_get_name: (virtual get_name)
 * @mcs: a #NcmMCSampler.
 *
 * Returns: the name of the sampler. 
 * 
 */
const gchar *
ncm_mc_sampler_get_name (NcmMCSampler *mcs)
{
  return NCM_MC_SAMPLER_GET_CLASS (mcs)->get_name (mcs);
}

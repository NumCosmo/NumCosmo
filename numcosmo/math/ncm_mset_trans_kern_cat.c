/***************************************************************************
 *            ncm_mset_trans_kern_cat.c
 *
 *  Fri September 23 12:36:45 2016
 *  Copyright  2016  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * ncm_mset_trans_kern_cat.c
 * Copyright (C) 2016 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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
 * SECTION:ncm_mset_trans_kern_cat
 * @title: NcmMSetTransKernCat
 * @short_description: Catalog sampler.
 *
 * FIXME
 * 
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_mset_trans_kern_cat.h"
#include "math/ncm_c.h"
#include "math/ncm_stats_dist_nd_kde_gauss.h"
#include "ncm_enum_types.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_randist.h>
#endif /* NUMCOSMO_GIR_SCAN */

struct _NcmMSetTransKernCatPrivate
{
  NcmMSetCatalog *mcat;
  NcmMSetTransKernCatSampling stype;
  NcmStatsDistNdKDEGauss *rbf;
  gboolean rbf_prep;
};

enum
{
  PROP_0,
  PROP_MCAT,
  PROP_STYPE,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcmMSetTransKernCat, ncm_mset_trans_kern_cat, NCM_TYPE_MSET_TRANS_KERN);

static void
ncm_mset_trans_kern_cat_init (NcmMSetTransKernCat *tcat)
{
  NcmMSetTransKernCatPrivate * const self = tcat->priv = ncm_mset_trans_kern_cat_get_instance_private (tcat);
  
  self->mcat     = NULL;
  self->stype    = NCM_MSET_TRANS_KERN_CAT_SAMPLING_LEN;
  self->rbf      = NULL;
  self->rbf_prep = FALSE;
}

static void
_ncm_mset_trans_kern_cat_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmMSetTransKernCat *tcat = NCM_MSET_TRANS_KERN_CAT (object);
  NcmMSetTransKernCatPrivate * const self = tcat->priv;
  g_return_if_fail (NCM_IS_MSET_TRANS_KERN_CAT (object));

  switch (prop_id)
  {
    case PROP_MCAT:
      ncm_mset_catalog_clear (&self->mcat);
      self->mcat = g_value_dup_object (value);
      g_assert (self->mcat != NULL);
      break;
    case PROP_STYPE:
      ncm_mset_trans_kern_cat_set_sampling (tcat, g_value_get_enum (value));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_mset_trans_kern_cat_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmMSetTransKernCat *tcat = NCM_MSET_TRANS_KERN_CAT (object);
  NcmMSetTransKernCatPrivate * const self = tcat->priv;
  g_return_if_fail (NCM_IS_MSET_TRANS_KERN_CAT (object));

  switch (prop_id)
  {
    case PROP_MCAT:
      g_value_set_object (value, self->mcat);
      break;
    case PROP_STYPE:
      g_value_set_enum (value, ncm_mset_trans_kern_cat_get_sampling (tcat));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_mset_trans_kern_cat_dispose (GObject *object)
{
  NcmMSetTransKernCat *tcat = NCM_MSET_TRANS_KERN_CAT (object);
  NcmMSetTransKernCatPrivate * const self = tcat->priv;

  ncm_mset_catalog_clear (&self->mcat);
  ncm_stats_dist_nd_kde_gauss_clear (&self->rbf);
  
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_mset_trans_kern_cat_parent_class)->dispose (object);
}

static void
_ncm_mset_trans_kern_cat_finalize (GObject *object)
{
  
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_mset_trans_kern_cat_parent_class)->finalize (object);
}

static void _ncm_mset_trans_kern_cat_set_mset (NcmMSetTransKern *tkern, NcmMSet *mset);
static void _ncm_mset_trans_kern_cat_generate (NcmMSetTransKern *tkern, NcmVector *theta, NcmVector *thetastar, NcmRNG *rng);
static gdouble _ncm_mset_trans_kern_cat_pdf (NcmMSetTransKern *tkern, NcmVector *theta, NcmVector *thetastar);
static const gchar *_ncm_mset_trans_kern_cat_get_name (NcmMSetTransKern *tkern);

static void
ncm_mset_trans_kern_cat_class_init (NcmMSetTransKernCatClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcmMSetTransKernClass *tkern_class = NCM_MSET_TRANS_KERN_CLASS (klass);

  object_class->set_property = &_ncm_mset_trans_kern_cat_set_property;
  object_class->get_property = &_ncm_mset_trans_kern_cat_get_property;
  object_class->dispose      = &_ncm_mset_trans_kern_cat_dispose;
  object_class->finalize     = &_ncm_mset_trans_kern_cat_finalize;

  g_object_class_install_property (object_class,
                                   PROP_MCAT,
                                   g_param_spec_object ("catalog",
                                                        NULL,
                                                        "catalog",
                                                        NCM_TYPE_MSET_CATALOG,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_STYPE,
                                   g_param_spec_enum ("sampling-type",
                                                      NULL,
                                                      "Sampling method to use",
                                                      NCM_TYPE_MSET_TRANS_KERN_CAT_SAMPLING, NCM_MSET_TRANS_KERN_CAT_SAMPLING_RBF_INTERP,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  
  tkern_class->set_mset = &_ncm_mset_trans_kern_cat_set_mset;
  tkern_class->generate = &_ncm_mset_trans_kern_cat_generate;
  tkern_class->pdf      = &_ncm_mset_trans_kern_cat_pdf;
  tkern_class->get_name = &_ncm_mset_trans_kern_cat_get_name;
}

static void
_ncm_mset_trans_kern_cat_set_mset (NcmMSetTransKern *tkern, NcmMSet *mset)
{
  NcmMSetTransKernCat *tcat = NCM_MSET_TRANS_KERN_CAT (tkern);
  NcmMSetTransKernCatPrivate * const self = tcat->priv;
  NcmMSet *mcat_mset        = ncm_mset_catalog_peek_mset (self->mcat);
  NCM_UNUSED (mset);
  
  if (!ncm_mset_cmp (mset, mcat_mset, FALSE))
    g_error ("_ncm_mset_trans_kern_cat_set_mset: incompatible mset.");  
}

static void
_ncm_mset_trans_kern_cat_generate_choose (NcmMSetTransKern *tkern, NcmVector *theta, NcmVector *thetastar, NcmRNG *rng)
{
  NcmMSetTransKernCat *tcat = NCM_MSET_TRANS_KERN_CAT (tkern);
  NcmMSetTransKernCatPrivate * const self = tcat->priv;
  NcmMSet *mcat_mset        = ncm_mset_catalog_peek_mset (self->mcat);
  const guint cat_len       = ncm_mset_catalog_len (self->mcat);
  const guint theta_size    = ncm_mset_fparams_len (mcat_mset);

  g_assert_cmpuint (theta_size, ==, ncm_vector_len (theta));
  
  ncm_rng_lock (rng);
  while (TRUE)
  {
    gulong i = gsl_rng_uniform_int (rng->r, cat_len);
    NcmVector *row = ncm_mset_catalog_peek_row (self->mcat, i);

    ncm_vector_memcpy2 (thetastar, row, 0, ncm_mset_catalog_nadd_vals (self->mcat), theta_size);

    if (ncm_mset_fparam_valid_bounds (tkern->mset, thetastar))
      break;
  }
  ncm_rng_unlock (rng);
}

static void
_ncm_mset_trans_kern_cat_generate_rbf_interp (NcmMSetTransKern *tkern, NcmVector *theta, NcmVector *thetastar, NcmRNG *rng)
{
  NcmMSetTransKernCat *tcat = NCM_MSET_TRANS_KERN_CAT (tkern);
  NcmMSetTransKernCatPrivate * const self = tcat->priv;
  NcmMSet *mcat_mset        = ncm_mset_catalog_peek_mset (self->mcat);
  const guint cat_len       = ncm_mset_catalog_len (self->mcat);
  const guint theta_size    = ncm_mset_fparams_len (mcat_mset);

  g_assert_cmpuint (theta_size, ==, ncm_vector_len (theta));

  if (!self->rbf_prep)
  {
    const guint nchains = ncm_mset_catalog_nchains (self->mcat);
    const guint np      = nchains > 1 ? nchains : ((cat_len > 1000) ? (cat_len / 10) : cat_len);
    const guint n       = MIN (cat_len, np);
    const guint nadd    = ncm_mset_catalog_nadd_vals (self->mcat);
    const guint m2lnp_i = ncm_mset_catalog_get_m2lnp_var (self->mcat);
    NcmVector *m2lnp    = ncm_vector_new (np);
    NcmVector *last_row = NULL;
    gint i, j;

    ncm_stats_dist_nd_kde_gauss_clear (&self->rbf);

    self->rbf = ncm_stats_dist_nd_kde_gauss_new (theta_size, FALSE);

    j = 0;
    for (i = cat_len - n; i < cat_len; i++)
    {
      NcmVector *row_i = ncm_mset_catalog_peek_row (self->mcat, i);
      if ((last_row != NULL) && (ncm_vector_get (last_row, 0) == ncm_vector_get (row_i, 0)))
        continue;
      else
      {
        NcmVector *srow_i = ncm_vector_get_subvector (row_i, nadd, theta_size);

        ncm_stats_dist_nd_kde_gauss_add_obs (self->rbf, srow_i);
        ncm_vector_set (m2lnp, j, ncm_vector_get (row_i, m2lnp_i));
        j++;
        
        ncm_vector_free (srow_i);
      }
      last_row = row_i;
    }
    ncm_stats_dist_nd_prepare_interp (NCM_STATS_DIST_ND (self->rbf), m2lnp);
    ncm_vector_free (m2lnp);
  }

  ncm_rng_lock (rng);
  while (TRUE)
  {
    ncm_stats_dist_nd_sample (NCM_STATS_DIST_ND (self->rbf), thetastar, rng);

    if (ncm_mset_fparam_valid_bounds (tkern->mset, thetastar))
      break;
  }
  ncm_rng_unlock (rng);
}

static void
_ncm_mset_trans_kern_cat_generate_kde (NcmMSetTransKern *tkern, NcmVector *theta, NcmVector *thetastar, NcmRNG *rng)
{
  NcmMSetTransKernCat *tcat = NCM_MSET_TRANS_KERN_CAT (tkern);
  NcmMSetTransKernCatPrivate * const self = tcat->priv;
  NcmMSet *mcat_mset        = ncm_mset_catalog_peek_mset (self->mcat);
  const guint cat_len       = ncm_mset_catalog_len (self->mcat);
  const guint theta_size    = ncm_mset_fparams_len (mcat_mset);

  g_assert_cmpuint (theta_size, ==, ncm_vector_len (theta));

  if (!self->rbf_prep)
  {
    const guint nchains = ncm_mset_catalog_nchains (self->mcat);
    const guint np      = nchains > 1 ? nchains : ((cat_len > 1000) ? (cat_len / 10) : cat_len);
    const guint n       = MIN (cat_len, np);
    const guint nadd    = ncm_mset_catalog_nadd_vals (self->mcat);
    NcmVector *last_row = NULL;
    gint i;

    ncm_stats_dist_nd_kde_gauss_clear (&self->rbf);

    self->rbf = ncm_stats_dist_nd_kde_gauss_new (theta_size, FALSE);

    for (i = cat_len - n; i < cat_len; i++)
    {
      NcmVector *row_i = ncm_mset_catalog_peek_row (self->mcat, i);
      if ((last_row != NULL) && (ncm_vector_get (last_row, 0) == ncm_vector_get (row_i, 0)))
        continue;
      else
      {
        NcmVector *srow_i = ncm_vector_get_subvector (row_i, nadd, theta_size);

        ncm_stats_dist_nd_kde_gauss_add_obs (self->rbf, srow_i);
        
        ncm_vector_free (srow_i);
      }
      last_row = row_i;
    }
    ncm_stats_dist_nd_prepare (NCM_STATS_DIST_ND (self->rbf));
  }

  ncm_rng_lock (rng);
  while (TRUE)
  {
    ncm_stats_dist_nd_sample (NCM_STATS_DIST_ND (self->rbf), thetastar, rng);

    if (ncm_mset_fparam_valid_bounds (tkern->mset, thetastar))
      break;
  }
  ncm_rng_unlock (rng);
}

static void
_ncm_mset_trans_kern_cat_generate (NcmMSetTransKern *tkern, NcmVector *theta, NcmVector *thetastar, NcmRNG *rng)
{
  NcmMSetTransKernCat *tcat = NCM_MSET_TRANS_KERN_CAT (tkern);
  NcmMSetTransKernCatPrivate * const self = tcat->priv;
  switch (self->stype)
  {
    case NCM_MSET_TRANS_KERN_CAT_SAMPLING_CHOOSE:
      _ncm_mset_trans_kern_cat_generate_choose (tkern, theta, thetastar, rng);
      break;
    case NCM_MSET_TRANS_KERN_CAT_SAMPLING_RBF_INTERP:
      _ncm_mset_trans_kern_cat_generate_rbf_interp (tkern, theta, thetastar, rng);
      break;
    case NCM_MSET_TRANS_KERN_CAT_SAMPLING_KDE:
      _ncm_mset_trans_kern_cat_generate_kde (tkern, theta, thetastar, rng);
      break;
    default:
      g_assert_not_reached ();
      break;
  }
}

static gdouble
_ncm_mset_trans_kern_cat_pdf (NcmMSetTransKern *tkern, NcmVector *theta, NcmVector *thetastar)
{
  g_error ("NcmMSetTransKernCat: does not implement pdf.");
  return 0.0;
}

static const gchar *
_ncm_mset_trans_kern_cat_get_name (NcmMSetTransKern *tkern)
{
  return "Catalog Sampler";
}

/**
 * ncm_mset_trans_kern_cat_new:
 * @mcat: a #NcmMSetCatalog
 *
 * New NcmMSetTransKernCat from @mcat catalog .
 *
 * Returns: (transfer full): a new #NcmMSetTransKernCat.
 *
 */
NcmMSetTransKernCat *
ncm_mset_trans_kern_cat_new (NcmMSetCatalog *mcat)
{
  NcmMSetTransKernCat *tcat = g_object_new (NCM_TYPE_MSET_TRANS_KERN_CAT,
                                            "catalog", mcat,
                                            NULL);
  return tcat;
}

/**
 * ncm_mset_trans_kern_cat_set_sampling:
 * @tcat: a #NcmMSetTransKernCat
 * @sampling: a #NcmMSetTransKernCatSampling
 *
 * Sets the sampling type to @sampling.
 * 
 */
void 
ncm_mset_trans_kern_cat_set_sampling (NcmMSetTransKernCat *tcat, NcmMSetTransKernCatSampling sampling)
{
  NcmMSetTransKernCatPrivate * const self = tcat->priv;
  self->stype = sampling;
}

/**
 * ncm_mset_trans_kern_cat_get_sampling:
 * @tcat: a #NcmMSetTransKernCat
 *
 * Gets the sampling type.
 * 
 * Returns: the current sampling type #NcmMSetTransKernCat.
 */
NcmMSetTransKernCatSampling 
ncm_mset_trans_kern_cat_get_sampling (NcmMSetTransKernCat *tcat)
{
  NcmMSetTransKernCatPrivate * const self = tcat->priv;
  return self->stype;
}

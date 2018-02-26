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

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_randist.h>
#endif /* NUMCOSMO_GIR_SCAN */

enum
{
  PROP_0,
  PROP_MCAT,
};

G_DEFINE_TYPE (NcmMSetTransKernCat, ncm_mset_trans_kern_cat, NCM_TYPE_MSET_TRANS_KERN);

static void
ncm_mset_trans_kern_cat_init (NcmMSetTransKernCat *tcat)
{
  tcat->mcat = NULL;
}

static void
_ncm_mset_trans_kern_cat_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmMSetTransKernCat *tcat = NCM_MSET_TRANS_KERN_CAT (object);
  g_return_if_fail (NCM_IS_MSET_TRANS_KERN_CAT (object));

  switch (prop_id)
  {
    case PROP_MCAT:
      ncm_mset_catalog_clear (&tcat->mcat);
      tcat->mcat = g_value_dup_object (value);
      g_assert (tcat->mcat != NULL);
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
  g_return_if_fail (NCM_IS_MSET_TRANS_KERN_CAT (object));

  switch (prop_id)
  {
    case PROP_MCAT:
      g_value_set_object (value, tcat->mcat);
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

  ncm_mset_catalog_clear (&tcat->mcat);
  
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
  
  tkern_class->set_mset = &_ncm_mset_trans_kern_cat_set_mset;
  tkern_class->generate = &_ncm_mset_trans_kern_cat_generate;
  tkern_class->pdf      = &_ncm_mset_trans_kern_cat_pdf;
  tkern_class->get_name = &_ncm_mset_trans_kern_cat_get_name;
}

static void
_ncm_mset_trans_kern_cat_set_mset (NcmMSetTransKern *tkern, NcmMSet *mset)
{
  NcmMSetTransKernCat *tcat = NCM_MSET_TRANS_KERN_CAT (tkern);
  NcmMSet *mcat_mset        = ncm_mset_catalog_peek_mset (tcat->mcat);
  NCM_UNUSED (mset);
  
  if (!ncm_mset_cmp (mset, mcat_mset, FALSE))
    g_error ("_ncm_mset_trans_kern_cat_set_mset: incompatible mset.");  
}

static void
_ncm_mset_trans_kern_cat_generate (NcmMSetTransKern *tkern, NcmVector *theta, NcmVector *thetastar, NcmRNG *rng)
{
  NcmMSetTransKernCat *tcat = NCM_MSET_TRANS_KERN_CAT (tkern);
  NcmMSet *mcat_mset        = ncm_mset_catalog_peek_mset (tcat->mcat);
  const guint cat_len       = ncm_mset_catalog_len (tcat->mcat);
  const guint theta_size    = ncm_mset_fparams_len (mcat_mset);

  g_assert_cmpuint (theta_size, ==, ncm_vector_len (theta));
  
  ncm_rng_lock (rng);
  while (TRUE)
  {
    gulong i = gsl_rng_uniform_int (rng->r, cat_len);
    NcmVector *row = ncm_mset_catalog_peek_row (tcat->mcat, i);

    ncm_vector_memcpy2 (thetastar, row, 0, ncm_mset_catalog_nadd_vals (tcat->mcat), theta_size);

    if (ncm_mset_fparam_valid_bounds (tkern->mset, thetastar))
      break;
  }
  ncm_rng_unlock (rng);
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

/***************************************************************************
 *            ncm_mset_trans_kern_cat.c
 *
 *  Fri September 23 12:36:45 2016
 *  Copyright  2016  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_mset_trans_kern_cat.c
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
#include "math/ncm_stats_dist_kde.h"
#include "math/ncm_stats_dist_vkde.h"
#include "math/ncm_stats_dist_kernel_gauss.h"
#include "math/ncm_stats_dist_kernel_st.h"
#include "ncm_enum_types.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rstat.h>
#endif /* NUMCOSMO_GIR_SCAN */

struct _NcmMSetTransKernCatPrivate
{
  NcmMSetCatalog *mcat;
  NcmMSetTransKernCatSampling stype;
  NcmStatsDist *sd;
  gboolean sd_prep;
  gdouble m2lnL_reltol;
  GTree *m2lnL_tree;
  gdouble choose_nsigma;
  gdouble cur_choose_nsigma;
};

enum
{
  PROP_0,
  PROP_MCAT,
  PROP_SD,
  PROP_STYPE,
  PROP_M2LNL_RELTOL,
  PROP_CHOOSE_NSIGMA,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcmMSetTransKernCat, ncm_mset_trans_kern_cat, NCM_TYPE_MSET_TRANS_KERN);

static gint gdouble_compare (gconstpointer a, gconstpointer b, gpointer user_data);

static void
ncm_mset_trans_kern_cat_init (NcmMSetTransKernCat *tcat)
{
  NcmMSetTransKernCatPrivate * const self = tcat->priv = ncm_mset_trans_kern_cat_get_instance_private (tcat);

  self->mcat         = NULL;
  self->stype        = NCM_MSET_TRANS_KERN_CAT_SAMPLING_LEN;
  self->sd           = NULL;
  self->sd_prep      = FALSE;
  self->m2lnL_reltol = 0.0;
  self->m2lnL_tree   = g_tree_new_full (gdouble_compare,
                                        &self->m2lnL_reltol,
                                        g_free,
                                        NULL);
  self->choose_nsigma     = 0.0;
  self->cur_choose_nsigma = 0.0;
}

static void
_ncm_mset_trans_kern_cat_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmMSetTransKernCat *tcat               = NCM_MSET_TRANS_KERN_CAT (object);
  NcmMSetTransKernCatPrivate * const self = tcat->priv;

  g_return_if_fail (NCM_IS_MSET_TRANS_KERN_CAT (object));

  switch (prop_id)
  {
    case PROP_MCAT:
    {
      ncm_mset_catalog_clear (&self->mcat);
      self->mcat = g_value_dup_object (value);
      g_assert (self->mcat != NULL);
      break;
    }
    case PROP_SD:
      self->sd = g_value_dup_object (value);
      break;
    case PROP_STYPE:
      ncm_mset_trans_kern_cat_set_sampling (tcat, g_value_get_enum (value));
      break;
    case PROP_M2LNL_RELTOL:
      self->m2lnL_reltol = g_value_get_double (value);
      break;
    case PROP_CHOOSE_NSIGMA:
      self->choose_nsigma     = g_value_get_double (value);
      self->cur_choose_nsigma = self->choose_nsigma;
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_mset_trans_kern_cat_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmMSetTransKernCat *tcat               = NCM_MSET_TRANS_KERN_CAT (object);
  NcmMSetTransKernCatPrivate * const self = tcat->priv;

  g_return_if_fail (NCM_IS_MSET_TRANS_KERN_CAT (object));

  switch (prop_id)
  {
    case PROP_MCAT:
      g_value_set_object (value, self->mcat);
      break;
    case PROP_SD:
      g_value_set_object (value, self->sd);
      break;
    case PROP_STYPE:
      g_value_set_enum (value, ncm_mset_trans_kern_cat_get_sampling (tcat));
      break;
    case PROP_M2LNL_RELTOL:
      g_value_set_double (value, self->m2lnL_reltol);
      break;
    case PROP_CHOOSE_NSIGMA:
      g_value_set_double (value, self->choose_nsigma);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_mset_trans_kern_cat_dispose (GObject *object)
{
  NcmMSetTransKernCat *tcat               = NCM_MSET_TRANS_KERN_CAT (object);
  NcmMSetTransKernCatPrivate * const self = tcat->priv;

  ncm_mset_catalog_clear (&self->mcat);
  ncm_stats_dist_clear (&self->sd);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_mset_trans_kern_cat_parent_class)->dispose (object);
}

static void
_ncm_mset_trans_kern_cat_finalize (GObject *object)
{
  NcmMSetTransKernCat *tcat               = NCM_MSET_TRANS_KERN_CAT (object);
  NcmMSetTransKernCatPrivate * const self = tcat->priv;

  g_tree_destroy (self->m2lnL_tree);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_mset_trans_kern_cat_parent_class)->finalize (object);
}

static void _ncm_mset_trans_kern_cat_set_mset (NcmMSetTransKern *tkern, NcmMSet *mset);
static void _ncm_mset_trans_kern_cat_generate (NcmMSetTransKern *tkern, NcmVector *theta, NcmVector *thetastar, NcmRNG *rng);
static gdouble _ncm_mset_trans_kern_cat_pdf (NcmMSetTransKern *tkern, NcmVector *theta, NcmVector *thetastar);
void _ncm_mset_trans_kern_cat_reset (NcmMSetTransKern *tkern);
static const gchar *_ncm_mset_trans_kern_cat_get_name (NcmMSetTransKern *tkern);

static void
ncm_mset_trans_kern_cat_class_init (NcmMSetTransKernCatClass *klass)
{
  GObjectClass *object_class         = G_OBJECT_CLASS (klass);
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
                                   PROP_SD,
                                   g_param_spec_object ("stats-dist",
                                                        NULL,
                                                        "NcmStatsDist object",
                                                        NCM_TYPE_STATS_DIST,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_STYPE,
                                   g_param_spec_enum ("sampling-type",
                                                      NULL,
                                                      "Sampling method to use",
                                                      NCM_TYPE_MSET_TRANS_KERN_CAT_SAMPLING, NCM_MSET_TRANS_KERN_CAT_SAMPLING_RBF_INTERP,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_M2LNL_RELTOL,
                                   g_param_spec_double ("m2lnL-reltol",
                                                        NULL,
                                                        "Relative tolerance for m2lnL",
                                                        GSL_DBL_EPSILON, 1.0e-3, 1.0e-7,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_CHOOSE_NSIGMA,
                                   g_param_spec_double ("choose-nsigma",
                                                        NULL,
                                                        "Max number of sigmas to use when choosing rows",
                                                        1.0e-1, G_MAXDOUBLE, 3.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  tkern_class->set_mset = &_ncm_mset_trans_kern_cat_set_mset;
  tkern_class->generate = &_ncm_mset_trans_kern_cat_generate;
  tkern_class->pdf      = &_ncm_mset_trans_kern_cat_pdf;
  tkern_class->reset    = &_ncm_mset_trans_kern_cat_reset;
  tkern_class->get_name = &_ncm_mset_trans_kern_cat_get_name;
}

static gint
gdouble_compare (gconstpointer a, gconstpointer b, gpointer user_data)
{
  gdouble *reltol = (gdouble *) user_data;

  return gsl_fcmp (*((gdouble *) a), *((gdouble *) b), reltol[0]);
}

static void
_ncm_mset_trans_kern_cat_set_mset (NcmMSetTransKern *tkern, NcmMSet *mset)
{
  NcmMSetTransKernCat *tcat               = NCM_MSET_TRANS_KERN_CAT (tkern);
  NcmMSetTransKernCatPrivate * const self = tcat->priv;
  NcmMSet *mcat_mset                      = ncm_mset_catalog_peek_mset (self->mcat);

  NCM_UNUSED (mset);

  if (!ncm_mset_cmp (mset, mcat_mset, FALSE))
    g_error ("_ncm_mset_trans_kern_cat_set_mset: incompatible mset.");
}

static void
_ncm_mset_trans_kern_cat_generate_choose (NcmMSetTransKern *tkern, NcmVector *theta, NcmVector *thetastar, NcmRNG *rng)
{
  NcmMSetTransKernCat *tcat               = NCM_MSET_TRANS_KERN_CAT (tkern);
  NcmMSetTransKernCatPrivate * const self = tcat->priv;
  NcmMSet *mcat_mset                      = ncm_mset_catalog_peek_mset (self->mcat);
  const guint cat_len                     = ncm_mset_catalog_len (self->mcat);
  const guint theta_size                  = ncm_mset_fparams_len (mcat_mset);
  const gdouble m2lnL_bestfit             = ncm_mset_catalog_get_bestfit_m2lnL (self->mcat);
  const gdouble m2lnL_sigma               = sqrt (2.0 * theta_size);
  const guint m2lnL_index                 = ncm_mset_catalog_get_m2lnp_var (self->mcat);
  const guint max_iter                    = 100000;
  guint under                             = 0;
  guint iter                              = 0;

  g_assert_cmpuint (theta_size, ==, ncm_vector_len (theta));

  ncm_rng_lock (rng);

  while (++iter < max_iter)
  {
    gulong i        = gsl_rng_uniform_int (rng->r, cat_len);
    NcmVector *row  = ncm_mset_catalog_peek_row (self->mcat, i);
    gdouble m2lnL_i = ncm_vector_get (row, m2lnL_index);

    if (m2lnL_i > m2lnL_bestfit + self->cur_choose_nsigma * m2lnL_sigma)
    {
      under++;

      if (under > 100)
      {
        under                    = 0;
        self->cur_choose_nsigma += 1.0;
      }

      /*g_message ("# Too far % 22.15g <=> % 22.15g\n", m2lnL_i, m2lnL_bestfit + self->choose_nsigma * m2lnL_sigma); */
      continue;
    }

    under = 0;

    ncm_vector_memcpy2 (thetastar, row, 0, ncm_mset_catalog_nadd_vals (self->mcat), theta_size);

    if (ncm_mset_fparam_valid_bounds (tkern->mset, thetastar))
    {
      if (g_tree_lookup_extended (self->m2lnL_tree, &m2lnL_i, NULL, NULL))
      {
        continue;
      }
      else
      {
        gdouble *m2lnL_i_copy = g_new (gdouble, 1);

        m2lnL_i_copy[0] = m2lnL_i;

        g_tree_insert (self->m2lnL_tree, m2lnL_i_copy, NULL);
        break;
      }
    }
  }

  if (iter >= max_iter)
    g_error ("_ncm_mset_trans_kern_cat_generate_choose: max_iter reached.");

  ncm_rng_unlock (rng);
}

static void
_ncm_mset_trans_kern_cat_generate_rbf_interp (NcmMSetTransKern *tkern, NcmVector *theta, NcmVector *thetastar, NcmRNG *rng)
{
  NcmMSetTransKernCat *tcat               = NCM_MSET_TRANS_KERN_CAT (tkern);
  NcmMSetTransKernCatPrivate * const self = tcat->priv;
  NcmMSet *mcat_mset                      = ncm_mset_catalog_peek_mset (self->mcat);
  const guint cat_len                     = ncm_mset_catalog_len (self->mcat);
  const guint theta_size                  = ncm_mset_fparams_len (mcat_mset);

  g_assert_cmpuint (theta_size, ==, ncm_vector_len (theta));

  if (!self->sd_prep)
  {
    const guint nchains = ncm_mset_catalog_nchains (self->mcat);
    const guint np = nchains > 1 ? nchains : ((cat_len > 1000) ? (cat_len / 10) : cat_len);
    const guint n = MIN (cat_len, np);
    const guint nadd = ncm_mset_catalog_nadd_vals (self->mcat);
    const guint m2lnp_i = ncm_mset_catalog_get_m2lnp_var (self->mcat);
    NcmVector *m2lnp = ncm_vector_new (np);
    NcmVector *last_row = NULL;
    gint i, j;

    ncm_stats_dist_reset (self->sd);

    j = 0;

    for (i = cat_len - n; i < cat_len; i++)
    {
      NcmVector *row_i = ncm_mset_catalog_peek_row (self->mcat, i);

      if ((last_row != NULL) && (ncm_vector_get (last_row, 0) == ncm_vector_get (row_i, 0)))
      {
        continue;
      }
      else
      {
        NcmVector *srow_i = ncm_vector_get_subvector (row_i, nadd, theta_size);

        ncm_stats_dist_add_obs (self->sd, srow_i);
        ncm_vector_set (m2lnp, j, ncm_vector_get (row_i, m2lnp_i));
        j++;

        ncm_vector_free (srow_i);
      }

      last_row = row_i;
    }

    ncm_stats_dist_prepare_interp (self->sd, m2lnp);
    ncm_vector_free (m2lnp);

    self->sd_prep = TRUE;
  }

  ncm_rng_lock (rng);

  while (TRUE)
  {
    ncm_stats_dist_sample (self->sd, thetastar, rng);

    if (ncm_mset_fparam_valid_bounds (tkern->mset, thetastar))
      break;
  }

  ncm_rng_unlock (rng);
}

static void
_ncm_mset_trans_kern_cat_generate_kde (NcmMSetTransKern *tkern, NcmVector *theta, NcmVector *thetastar, NcmRNG *rng)
{
  NcmMSetTransKernCat *tcat               = NCM_MSET_TRANS_KERN_CAT (tkern);
  NcmMSetTransKernCatPrivate * const self = tcat->priv;
  NcmMSet *mcat_mset                      = ncm_mset_catalog_peek_mset (self->mcat);
  const guint cat_len                     = ncm_mset_catalog_len (self->mcat);
  const guint theta_size                  = ncm_mset_fparams_len (mcat_mset);

  g_assert_cmpuint (theta_size, ==, ncm_vector_len (theta));

  if (!self->sd_prep)
  {
    const guint nchains = ncm_mset_catalog_nchains (self->mcat);
    const guint np      = nchains > 1 ? nchains : ((cat_len > 1000) ? (cat_len / 10) : cat_len);
    const guint n       = MIN (cat_len, np);
    const guint nadd    = ncm_mset_catalog_nadd_vals (self->mcat);
    NcmVector *last_row = NULL;
    gint i;

    ncm_stats_dist_reset (self->sd);

    for (i = cat_len - n; i < cat_len; i++)
    {
      NcmVector *row_i = ncm_mset_catalog_peek_row (self->mcat, i);

      if ((last_row != NULL) && (ncm_vector_get (last_row, 0) == ncm_vector_get (row_i, 0)))
      {
        continue;
      }
      else
      {
        NcmVector *srow_i = ncm_vector_get_subvector (row_i, nadd, theta_size);

        ncm_stats_dist_add_obs (self->sd, srow_i);

        ncm_vector_free (srow_i);
      }

      last_row = row_i;
    }

    ncm_stats_dist_prepare (self->sd);
    self->sd_prep = TRUE;
  }

  ncm_rng_lock (rng);

  while (TRUE)
  {
    ncm_stats_dist_sample (self->sd, thetastar, rng);

    if (ncm_mset_fparam_valid_bounds (tkern->mset, thetastar))
      break;
  }

  ncm_rng_unlock (rng);
}

static void
_ncm_mset_trans_kern_cat_generate (NcmMSetTransKern *tkern, NcmVector *theta, NcmVector *thetastar, NcmRNG *rng)
{
  NcmMSetTransKernCat *tcat               = NCM_MSET_TRANS_KERN_CAT (tkern);
  NcmMSetTransKernCatPrivate * const self = tcat->priv;

  switch (self->stype)
  {
    case NCM_MSET_TRANS_KERN_CAT_SAMPLING_CHOOSE:
      _ncm_mset_trans_kern_cat_generate_choose (tkern, theta, thetastar, rng);
      break;
    case NCM_MSET_TRANS_KERN_CAT_SAMPLING_RBF_INTERP:

      if (self->sd == NULL)
        g_error ("NcmMSetTransKernCat: RBF interpolation requires a NcmStatsDist.");

      _ncm_mset_trans_kern_cat_generate_rbf_interp (tkern, theta, thetastar, rng);
      break;
    case NCM_MSET_TRANS_KERN_CAT_SAMPLING_KDE:

      if (self->sd == NULL)
        g_error ("NcmMSetTransKernCat: KDE requires a NcmStatsDist.");

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

void
_ncm_mset_trans_kern_cat_reset (NcmMSetTransKern *tkern)
{
  NcmMSetTransKernCat *tcat               = NCM_MSET_TRANS_KERN_CAT (tkern);
  NcmMSetTransKernCatPrivate * const self = tcat->priv;

  self->sd_prep = FALSE;

  if (self->sd != NULL)
    ncm_stats_dist_reset (self->sd);

  self->cur_choose_nsigma = self->choose_nsigma;
  g_tree_remove_all (self->m2lnL_tree);
}

static const gchar *
_ncm_mset_trans_kern_cat_get_name (NcmMSetTransKern *tkern)
{
  return "Catalog Sampler";
}

/**
 * ncm_mset_trans_kern_cat_new:
 * @mcat: a #NcmMSetCatalog
 * @sd: (allow-none): a #NcmStatsDist
 *
 * New NcmMSetTransKernCat from @mcat catalog with the interpolation
 * object @sd.
 *
 * Returns: (transfer full): a new #NcmMSetTransKernCat.
 *
 */
NcmMSetTransKernCat *
ncm_mset_trans_kern_cat_new (NcmMSetCatalog *mcat, NcmStatsDist *sd)
{
  NcmMSetTransKernCat *tcat = g_object_new (NCM_TYPE_MSET_TRANS_KERN_CAT,
                                            "catalog", mcat,
                                            "stats-dist", sd,
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

  self->stype   = sampling;
  self->sd_prep = FALSE;
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


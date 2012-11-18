/***************************************************************************
 *            ncm_likelihood.c
 *
 *  Wed May 30 15:36:41 2007
 *  Copyright  2007  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) 2012 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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
 * SECTION:ncm_likelihood
 * @title: Likelihood
 * @short_description: FIXME
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_likelihood.h"

enum
{
  PROP_0,
  PROP_DATASET,
  PROP_SIZE,
};

G_DEFINE_TYPE (NcmLikelihood, ncm_likelihood, G_TYPE_OBJECT);

static void
ncm_likelihood_init (NcmLikelihood *lh)
{
  lh->dset   = NULL;
  lh->priors = g_ptr_array_sized_new (10);
  g_ptr_array_set_free_func (lh->priors, (GDestroyNotify) &ncm_mset_func_free);
}

static void
ncm_likelihood_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmLikelihood *lh = NCM_LIKELIHOOD (object);
  g_return_if_fail (NCM_IS_LIKELIHOOD (object));

  switch (prop_id)
  {
    case PROP_DATASET:
      ncm_dataset_clear (&lh->dset);
      lh->dset = g_value_dup_object (value);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
ncm_likelihood_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmLikelihood *lh = NCM_LIKELIHOOD (object);
  g_return_if_fail (NCM_IS_LIKELIHOOD (object));

  switch (prop_id)
  {
    case PROP_DATASET:
      g_value_set_object (value, lh->dset);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
ncm_likelihood_dispose (GObject *object)
{
  NcmLikelihood *lh = NCM_LIKELIHOOD (object);

  if (lh->priors != NULL)
  {
    g_ptr_array_unref (lh->priors);
    lh->priors = NULL;
  }

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_likelihood_parent_class)->dispose (object);
}

static void
ncm_likelihood_finalize (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_likelihood_parent_class)->finalize (object);
}

static void
ncm_likelihood_class_init (NcmLikelihoodClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = &ncm_likelihood_set_property;
  object_class->get_property = &ncm_likelihood_get_property;
  object_class->dispose      = &ncm_likelihood_dispose;
  object_class->finalize     = &ncm_likelihood_finalize;

  g_object_class_install_property (object_class,
                                   PROP_DATASET,
                                   g_param_spec_object ("dataset",
                                                        NULL,
                                                        "Dataset object",
                                                        NCM_TYPE_DATASET,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  
}

/**
 * ncm_likelihood_new:
 * @dset: a #NcmDataset.
 *
 * FIXME
 *
 * Returns: FIXME
 */
NcmLikelihood *
ncm_likelihood_new (NcmDataset *dset)
{
  NcmLikelihood *lh = g_object_new (NCM_TYPE_LIKELIHOOD, 
                                    "dataset", dset,
                                    NULL);
  return lh;
}

/**
 * ncm_likelihood_ref:
 * @lh: FIXME
 *
 * FIXME
 *
 * Returns: (transfer full): FIXME
   */
NcmLikelihood *
ncm_likelihood_ref (NcmLikelihood *lh)
{
  return g_object_ref (lh);
}

/**
 * ncm_likelihood_dup:
 * @lh: a #NcmLikelihood.
 *
 * Duplicates the object and all of its content.
 *
 * Returns: (transfer full): A duplicate of @lh.
 */
NcmLikelihood *
ncm_likelihood_dup (NcmLikelihood *lh)
{
  NcmDataset *dset = ncm_dataset_dup (lh->dset);
  NcmLikelihood *lh_dup = ncm_likelihood_new (dset);
  gint i;
  
  ncm_dataset_free (dset);

  g_ptr_array_set_size (lh_dup->priors, 0);
  for (i = 0; i < lh->priors->len; i++)
  {
    NcmMSetFunc *func = NCM_MSET_FUNC (g_ptr_array_index (lh->priors, i));
    ncm_likelihood_priors_add (lh_dup, ncm_mset_func_ref (func));
  }

  return lh_dup;
}


/**
 * ncm_likelihood_copy:
 * @lh: a #NcmLikelihood.
 *
 * Duplicates the object and gets a reference for its content.
 *
 * Returns: (transfer full): A copy of @lh.
 */
NcmLikelihood *
ncm_likelihood_copy (NcmLikelihood *lh)
{
  NcmDataset *dset = ncm_dataset_copy (lh->dset);
  NcmLikelihood *lh_ref = ncm_likelihood_new (dset);
  gint i;
  
  ncm_dataset_free (dset);

  g_ptr_array_set_size (lh_ref->priors, 0);
  for (i = 0; i < lh->priors->len; i++)
  {
    NcmMSetFunc *func = NCM_MSET_FUNC (g_ptr_array_index (lh->priors, i));
    ncm_likelihood_priors_add (lh_ref, ncm_mset_func_ref (func));
  }

  return lh_ref;
}

/**
 * ncm_likelihood_free:
 * @lh: FIXME
 *
 * FIXME
 */
void
ncm_likelihood_free (NcmLikelihood *lh)
{
  g_object_unref (lh);
}

/**
 * ncm_likelihood_clear:
 * @lh: FIXME
 *
 * FIXME
 */
void
ncm_likelihood_clear (NcmLikelihood **lh)
{
  g_clear_object (lh);
}

/**
 * ncm_likelihood_priors_add:
 * @lh: FIXME
 * @prior: FIXME
 *
 * FIXME
 */
void
ncm_likelihood_priors_add (NcmLikelihood *lh, NcmMSetFunc *prior)
{
  g_ptr_array_add (lh->priors, ncm_mset_func_ref (prior));
}

/**
 * ncm_likelihood_priors_peek:
 * @lh: FIXME
 * @i: FIXME
 *
 * FIXME
 * 
 * Returns: (transfer none): FIXME
 */
NcmMSetFunc *
ncm_likelihood_priors_peek (NcmLikelihood *lh, guint i)
{
  return g_ptr_array_index (lh->priors, i);
}

/**
 * ncm_likelihood_priors_length:
 * @lh: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
guint
ncm_likelihood_priors_length (NcmLikelihood *lh)
{
  return lh->priors->len;
}

/**
 * ncm_likelihood_has_leastsquares_J:
 * @lh: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gboolean
ncm_likelihood_has_leastsquares_J (NcmLikelihood *lh)
{
  const gdouble has_priors = !(ncm_likelihood_priors_length (lh) == 0);
  const gdouble dataset_has_leastsquares_J = ncm_dataset_has_leastsquares_J (lh->dset);
  return dataset_has_leastsquares_J && !has_priors;
}

/**
 * ncm_likelihood_has_m2lnL_grad:
 * @lh: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gboolean
ncm_likelihood_has_m2lnL_grad (NcmLikelihood *lh)
{
  const gdouble has_priors = !(ncm_likelihood_priors_length (lh) == 0);
  const gdouble dataset_has_m2lnL_grad = ncm_dataset_has_m2lnL_grad (lh->dset);
  return dataset_has_m2lnL_grad && !has_priors;
}

/**
 * ncm_likelihood_priors_leastsquares_f:
 * @lh: a #NcmLikelihood.
 * @mset: a #NcmMSet.
 * @priors_f: a #NcmVector.
 *
 * FIXME
 */
void
ncm_likelihood_priors_leastsquares_f (NcmLikelihood *lh, NcmMSet *mset, NcmVector *priors_f)
{
  guint i = 0;

  for (i = 0; i < lh->priors->len; i++)
  {
    NcmMSetFunc *func = NCM_MSET_FUNC (g_ptr_array_index (lh->priors, i));
    const gdouble prior = ncm_mset_func_eval0 (func, mset);
    ncm_vector_set (priors_f, i, prior);
  }

  return;
}

/**
 * ncm_likelihood_leastsquares_f:
 * @lh: a #NcmLikelihood.
 * @mset: a #NcmMSet.
 * @f: a #NcmVector.
 *
 * FIXME
 * 
 */
void
ncm_likelihood_leastsquares_f (NcmLikelihood *lh, NcmMSet *mset, NcmVector *f)
{
  guint data_size = ncm_dataset_get_length (lh->dset);
  guint priors_size = ncm_likelihood_priors_length (lh);
  
  if (data_size)
  {
    NcmVector *data_f = ncm_vector_get_subvector (f, 0, data_size);
    ncm_dataset_leastsquares_f (lh->dset, mset, data_f);
    ncm_vector_free (data_f);
  }

  if (priors_size)
  {
    NcmVector *priors_f = ncm_vector_get_subvector (f, data_size, priors_size);
    ncm_likelihood_priors_leastsquares_f (lh, mset, priors_f);
    ncm_vector_free (priors_f);
  }
  
  return;
}

/**
 * ncm_likelihood_leastsquares_J:
 * @lh: a #NcmLikelihood.
 * @mset: a #NcmMSet.
 * @J: a #NcmMatrix.
 *
 * FIXME
 * 
 */
void
ncm_likelihood_leastsquares_J (NcmLikelihood *lh, NcmMSet *mset, NcmMatrix *J)
{
  g_assert (lh->priors->len == 0);
  ncm_dataset_leastsquares_J (lh->dset, mset, J);
}

/**
 * ncm_likelihood_leastsquares_f_J:
 * @lh: a #NcmLikelihood.
 * @mset: a #NcmMSet.
 * @f: a #NcmVector.
 * @J: a #NcmMatrix.
 *
 * FIXME
 */
void
ncm_likelihood_leastsquares_f_J (NcmLikelihood *lh, NcmMSet *mset, NcmVector *f, NcmMatrix *J)
{
  g_assert (lh->priors->len == 0);
  ncm_dataset_leastsquares_f_J (lh->dset, mset, f, J);
}

/**
 * ncm_likelihood_priors_m2lnL_val:
 * @lh: a #NcmLikelihood.
 * @mset: a #NcmMSet.
 * @priors_m2lnL: (out): FIXME
 * 
 * FIXME
 */
void
ncm_likelihood_priors_m2lnL_val (NcmLikelihood *lh, NcmMSet *mset, gdouble *priors_m2lnL)
{
  guint i = 0;

  for (i = 0; i < lh->priors->len; i++)
  {
    NcmMSetFunc *func = NCM_MSET_FUNC (g_ptr_array_index (lh->priors, i));
    const gdouble prior = ncm_mset_func_eval0 (func, mset);
    *priors_m2lnL += prior * prior;
  }

  return;
}

/**
 * ncm_likelihood_m2lnL_val:
 * @lh: a #NcmLikelihood.
 * @mset: a #NcmMSet.
 * @m2lnL: (out): FIXME
 * 
 * FIXME
 * 
 */
void
ncm_likelihood_m2lnL_val (NcmLikelihood *lh, NcmMSet *mset, gdouble *m2lnL)
{
  gdouble data_m2lnL;
  gdouble priors_m2lnL;
  
  ncm_dataset_m2lnL_val (lh->dset, mset, &data_m2lnL);
  ncm_likelihood_priors_m2lnL_val (lh, mset, &priors_m2lnL);

  *m2lnL = data_m2lnL + priors_m2lnL;
  return;
}

/**
 * ncm_likelihood_m2lnL_grad:
 * @lh: a #NcmLikelihood.
 * @mset: a #NcmMSet.
 * @grad: a #NcmVector.
 *
 * FIXME
 */
void
ncm_likelihood_m2lnL_grad (NcmLikelihood *lh, NcmMSet *mset, NcmVector *grad)
{
  g_assert (lh->priors->len == 0);
  ncm_dataset_m2lnL_grad (lh->dset, mset, grad);
}

/**
 * ncm_likelihood_m2lnL_val_grad:
 * @lh: a #NcmLikelihood.
 * @mset: a #NcmMSet.
 * @m2lnL: (out): FIXME
 * @grad: a #NcmVector.
 *
 * FIXME
 */
void
ncm_likelihood_m2lnL_val_grad (NcmLikelihood *lh, NcmMSet *mset, gdouble *m2lnL, NcmVector *grad)
{
  g_assert (lh->priors->len == 0);
  ncm_dataset_m2lnL_val_grad (lh->dset, mset, m2lnL, grad);
}

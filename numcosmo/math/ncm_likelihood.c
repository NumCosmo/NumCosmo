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
 * @title: NcmLikelihood
 * @short_description: Likelihood combining a NcmDataset and priors.
 * 
 * FIXME
 * 
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_likelihood.h"
#include "math/ncm_cfg.h"

enum
{
  PROP_0,
  PROP_DATASET,
  PROP_PRIORS_M2LNL,
  PROP_PRIORS_F,
  PROP_SIZE,
};

G_DEFINE_TYPE (NcmLikelihood, ncm_likelihood, G_TYPE_OBJECT);

static void
ncm_likelihood_init (NcmLikelihood *lh)
{
  lh->dset   = NULL;
  lh->priors_f = g_ptr_array_sized_new (10);
  lh->priors_m2lnL = g_ptr_array_sized_new (10);
  g_ptr_array_set_free_func (lh->priors_f, (GDestroyNotify) &ncm_mset_func_free);
  g_ptr_array_set_free_func (lh->priors_m2lnL, (GDestroyNotify) &ncm_mset_func_free);
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
    case PROP_PRIORS_M2LNL:
    {
      guint p = g_value_get_int (value);
      GPtrArray *priors_m2lnL = GINT_TO_POINTER (p);
      if (priors_m2lnL != lh->priors_m2lnL)
      {
        g_ptr_array_unref (lh->priors_m2lnL);
        lh->priors_m2lnL = g_ptr_array_ref (priors_m2lnL);
      }
      break;
    }
    case PROP_PRIORS_F:
    {
      guint p = g_value_get_int (value);
      GPtrArray *priors_f = GINT_TO_POINTER (p);
      if (priors_f != lh->priors_f)
      {
        g_ptr_array_unref (lh->priors_f);
        lh->priors_f = g_ptr_array_ref (priors_f);
      }
      break;
    }
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
    case PROP_PRIORS_M2LNL:
      g_value_set_int (value, GPOINTER_TO_INT (lh->priors_m2lnL));
      break;
    case PROP_PRIORS_F:
      g_value_set_int (value, GPOINTER_TO_INT (lh->priors_f));
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

  if (lh->priors_f != NULL)
  {
    g_ptr_array_unref (lh->priors_f);
    lh->priors_f = NULL;
  }
  if (lh->priors_m2lnL != NULL)
  {
    g_ptr_array_unref (lh->priors_m2lnL);
    lh->priors_m2lnL = NULL;
  }

  ncm_dataset_clear (&lh->dset);

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
  g_object_class_install_property (object_class,
                                   PROP_PRIORS_M2LNL,
                                   g_param_spec_int ("priors-m2lnL-ptr",
                                                     NULL,
                                                     "Priors m2lnL pointer",
                                                     G_MININT, G_MAXINT, 0,
                                                     G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_PRIORS_F,
                                   g_param_spec_int ("priors-f-ptr",
                                                     NULL,
                                                     "Priors f pointer",
                                                     G_MININT, G_MAXINT, 0,
                                                     G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
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
 * @ser: a #NcmSerialize.
 *
 * Duplicates the object and all of its content.
 *
 * Returns: (transfer full): A duplicate of @lh.
 */
NcmLikelihood *
ncm_likelihood_dup (NcmLikelihood *lh, NcmSerialize *ser)
{
  return NCM_LIKELIHOOD (ncm_serialize_dup_obj (ser, G_OBJECT (lh)));
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
  guint i;
  
  ncm_dataset_free (dset);

  g_ptr_array_set_size (lh_ref->priors_f, 0);
  g_ptr_array_set_size (lh_ref->priors_m2lnL, 0);

  for (i = 0; i < lh->priors_f->len; i++)
  {
    NcmMSetFunc *func = NCM_MSET_FUNC (g_ptr_array_index (lh->priors_f, i));
    ncm_likelihood_priors_add (lh_ref, ncm_mset_func_ref (func), FALSE);
  }
  for (i = 0; i < lh->priors_m2lnL->len; i++)
  {
    NcmMSetFunc *func = NCM_MSET_FUNC (g_ptr_array_index (lh->priors_m2lnL, i));
    ncm_likelihood_priors_add (lh_ref, ncm_mset_func_ref (func), TRUE);
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
 * @is_m2lnL: FIXME
 *
 * FIXME
 */
void
ncm_likelihood_priors_add (NcmLikelihood *lh, NcmMSetFunc *prior, gboolean is_m2lnL)
{
  if (is_m2lnL)
    g_ptr_array_add (lh->priors_m2lnL, ncm_mset_func_ref (prior));
  else
    g_ptr_array_add (lh->priors_f, ncm_mset_func_ref (prior));
}

/**
 * ncm_likelihood_priors_peek_f:
 * @lh: FIXME
 * @i: FIXME
 *
 * FIXME
 * 
 * Returns: (transfer none): FIXME
 */
NcmMSetFunc *
ncm_likelihood_priors_peek_f (NcmLikelihood *lh, guint i)
{
  return g_ptr_array_index (lh->priors_f, i);
}

/**
 * ncm_likelihood_priors_peek_m2lnL:
 * @lh: FIXME
 * @i: FIXME
 *
 * FIXME
 * 
 * Returns: (transfer none): FIXME
 */
NcmMSetFunc *
ncm_likelihood_priors_peek_m2lnL (NcmLikelihood *lh, guint i)
{
  return g_ptr_array_index (lh->priors_m2lnL, i);
}

/**
 * ncm_likelihood_priors_length_f:
 * @lh: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
guint
ncm_likelihood_priors_length_f (NcmLikelihood *lh)
{
  return lh->priors_f->len;
}

/**
 * ncm_likelihood_priors_length_m2lnL:
 * @lh: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
guint
ncm_likelihood_priors_length_m2lnL (NcmLikelihood *lh)
{
  return lh->priors_m2lnL->len;
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
  const gboolean has_priors = (ncm_likelihood_priors_length_f (lh) + 
    ncm_likelihood_priors_length_m2lnL (lh)) != 0;
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
  const gboolean has_priors = (ncm_likelihood_priors_length_f (lh) + 
    ncm_likelihood_priors_length_m2lnL (lh)) != 0;
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
  gboolean has_prior_m2lnL = ncm_likelihood_priors_length_m2lnL (lh) != 0;
  if (has_prior_m2lnL)
    g_error ("ncm_likelihood_priors_leastsquares_f: cannot calculate leastsquares f, the likelihood contains m2lnL priors.");
  
  for (i = 0; i < lh->priors_f->len; i++)
  {
    NcmMSetFunc *func = NCM_MSET_FUNC (g_ptr_array_index (lh->priors_f, i));
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
  guint data_size = ncm_dataset_get_n (lh->dset);
  guint priors_f_size = ncm_likelihood_priors_length_f (lh);
  gboolean has_prior_m2lnL = ncm_likelihood_priors_length_m2lnL (lh) != 0;
  if (has_prior_m2lnL)
    g_error ("ncm_likelihood_priors_leastsquares_f: cannot calculate leastsquares f, the likelihood contains m2lnL priors.");
  
  if (data_size)
  {
    NcmVector *data_f = ncm_vector_get_subvector (f, 0, data_size);
    ncm_dataset_leastsquares_f (lh->dset, mset, data_f);
    ncm_vector_free (data_f);
  }

  if (priors_f_size)
  {
    NcmVector *priors_f = ncm_vector_get_subvector (f, data_size, priors_f_size);
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
  g_assert (lh->priors_f->len == 0 && lh->priors_m2lnL->len == 0);
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
  g_assert (lh->priors_f->len == 0 && lh->priors_m2lnL->len == 0);
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

  *priors_m2lnL = 0.0;

  for (i = 0; i < lh->priors_f->len; i++)
  {
    NcmMSetFunc *func = NCM_MSET_FUNC (g_ptr_array_index (lh->priors_f, i));
    const gdouble prior = ncm_mset_func_eval0 (func, mset);
    *priors_m2lnL += prior * prior;
  }

  for (i = 0; i < lh->priors_m2lnL->len; i++)
  {
    NcmMSetFunc *func = NCM_MSET_FUNC (g_ptr_array_index (lh->priors_m2lnL, i));
    const gdouble prior = ncm_mset_func_eval0 (func, mset);
    *priors_m2lnL += prior;
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
  gdouble data_m2lnL = 0.0;
  gdouble priors_m2lnL = 0.0;
  
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
  g_assert (lh->priors_f->len == 0 && lh->priors_m2lnL->len == 0);
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
  g_assert (lh->priors_f->len == 0 && lh->priors_m2lnL->len == 0);
  ncm_dataset_m2lnL_val_grad (lh->dset, mset, m2lnL, grad);
}

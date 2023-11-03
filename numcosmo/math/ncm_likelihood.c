/***************************************************************************
 *            ncm_likelihood.c
 *
 *  Wed May 30 15:36:41 2007
 *  Copyright  2007  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) 2012 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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
#include "math/ncm_prior_gauss_param.h"
#include "math/ncm_prior_gauss_func.h"
#include "math/ncm_prior_flat_param.h"
#include "math/ncm_prior_flat_func.h"

enum
{
  PROP_0,
  PROP_DATASET,
  PROP_PRIORS_M2LNL,
  PROP_PRIORS_F,
  PROP_M2LNL_V,
  PROP_SIZE,
};

G_DEFINE_TYPE (NcmLikelihood, ncm_likelihood, G_TYPE_OBJECT);

static void
ncm_likelihood_init (NcmLikelihood *lh)
{
  lh->dset         = NULL;
  lh->m2lnL_v      = NULL;
  lh->priors_f     = ncm_obj_array_sized_new (10);
  lh->priors_m2lnL = ncm_obj_array_sized_new (10);
  lh->ls_f_data    = ncm_vector_new_data_static (GINT_TO_POINTER (1), 1, 1); /* Invalid vectors */
  lh->ls_f_priors  = ncm_vector_new_data_static (GINT_TO_POINTER (1), 1, 1); /* Invalid vectors */
}

static void
_ncm_likelihood_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
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
      NcmObjArray *priors_m2lnL_old = lh->priors_m2lnL;

      lh->priors_m2lnL = ncm_obj_array_ref (g_value_get_boxed (value));
      ncm_obj_array_unref (priors_m2lnL_old);
      break;
    }
    case PROP_PRIORS_F:
    {
      NcmObjArray *priors_f_old = lh->priors_f;

      lh->priors_f = ncm_obj_array_ref (g_value_get_boxed (value));
      ncm_obj_array_unref (priors_f_old);
      break;
    }
    case PROP_M2LNL_V:
    {
      NcmVector *m2lnL_v = g_value_get_object (value);

      if (m2lnL_v != NULL)
      {
        ncm_vector_clear (&lh->m2lnL_v);
        lh->m2lnL_v = ncm_vector_ref (m2lnL_v);
      }

      break;
    }
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_likelihood_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmLikelihood *lh = NCM_LIKELIHOOD (object);

  g_return_if_fail (NCM_IS_LIKELIHOOD (object));

  switch (prop_id)
  {
    case PROP_DATASET:
      g_value_set_object (value, lh->dset);
      break;
    case PROP_PRIORS_M2LNL:
      g_value_set_boxed (value, lh->priors_m2lnL);
      break;
    case PROP_PRIORS_F:
      g_value_set_boxed (value, lh->priors_f);
      break;
    case PROP_M2LNL_V:
      g_value_set_object (value, lh->m2lnL_v);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_likelihood_dispose (GObject *object)
{
  NcmLikelihood *lh = NCM_LIKELIHOOD (object);

  ncm_obj_array_clear (&lh->priors_f);
  ncm_obj_array_clear (&lh->priors_m2lnL);

  ncm_vector_clear (&lh->m2lnL_v);
  ncm_dataset_clear (&lh->dset);

  ncm_vector_clear (&lh->ls_f_data);
  ncm_vector_clear (&lh->ls_f_priors);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_likelihood_parent_class)->dispose (object);
}

static void
_ncm_likelihood_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_likelihood_parent_class)->finalize (object);
}

static void
ncm_likelihood_class_init (NcmLikelihoodClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = &_ncm_likelihood_set_property;
  object_class->get_property = &_ncm_likelihood_get_property;
  object_class->dispose      = &_ncm_likelihood_dispose;
  object_class->finalize     = &_ncm_likelihood_finalize;

  g_object_class_install_property (object_class,
                                   PROP_DATASET,
                                   g_param_spec_object ("dataset",
                                                        NULL,
                                                        "Dataset object",
                                                        NCM_TYPE_DATASET,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_PRIORS_M2LNL,
                                   g_param_spec_boxed ("priors-m2lnL",
                                                       NULL,
                                                       "Priors m2lnL array",
                                                       NCM_TYPE_OBJ_ARRAY,
                                                       G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_PRIORS_F,
                                   g_param_spec_boxed ("priors-f",
                                                       NULL,
                                                       "Priors f array",
                                                       NCM_TYPE_OBJ_ARRAY,
                                                       G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_M2LNL_V,
                                   g_param_spec_object ("m2lnL-v",
                                                        NULL,
                                                        "m2lnL vector",
                                                        NCM_TYPE_VECTOR,
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
 * @lh: a #NcmLikelihood
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
 * ncm_likelihood_free:
 * @lh: a #NcmLikelihood
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
 * @lh: a #NcmLikelihood
 *
 * FIXME
 *
 */
void
ncm_likelihood_clear (NcmLikelihood **lh)
{
  g_clear_object (lh);
}

/**
 * ncm_likelihood_priors_add:
 * @lh: a #NcmLikelihood
 * @prior: a #NcmPrior
 *
 * FIXME
 *
 */
void
ncm_likelihood_priors_add (NcmLikelihood *lh, NcmPrior *prior)
{
  if (ncm_prior_is_m2lnL (prior))
    g_ptr_array_add (lh->priors_m2lnL, ncm_prior_ref (prior));
  else
    g_ptr_array_add (lh->priors_f, ncm_prior_ref (prior));
}

/**
 * ncm_likelihood_priors_add_gauss_param:
 * @lh: a #NcmLikelihood
 * @mid: FIXME
 * @pid: FIXME
 * @mu: FIXME
 * @sigma: FIXME
 *
 * FIXME
 *
 */
void
ncm_likelihood_priors_add_gauss_param (NcmLikelihood *lh, NcmModelID mid, guint pid, gdouble mu, gdouble sigma)
{
  NcmPrior *prior = NCM_PRIOR (ncm_prior_gauss_param_new (mid, pid, mu, sigma));

  ncm_likelihood_priors_add (lh, prior);
  ncm_prior_clear (&prior);
}

/**
 * ncm_likelihood_priors_add_gauss_param_pindex:
 * @lh: a #NcmLikelihood
 * @pi: FIXME
 * @mu: FIXME
 * @sigma: FIXME
 *
 * FIXME
 *
 */
void
ncm_likelihood_priors_add_gauss_param_pindex (NcmLikelihood *lh, const NcmMSetPIndex *pi, gdouble mu, gdouble sigma)
{
  NcmPrior *prior = NCM_PRIOR (ncm_prior_gauss_param_new_pindex (pi, mu, sigma));

  ncm_likelihood_priors_add (lh, prior);
  ncm_prior_clear (&prior);
}

/**
 * ncm_likelihood_priors_add_gauss_param_name:
 * @lh: a #NcmLikelihood
 * @mset: FIXME
 * @name: FIXME
 * @mu: FIXME
 * @sigma: FIXME
 *
 * FIXME
 *
 */
void
ncm_likelihood_priors_add_gauss_param_name (NcmLikelihood *lh, NcmMSet *mset, const gchar *name, gdouble mu, gdouble sigma)
{
  NcmPrior *prior = NCM_PRIOR (ncm_prior_gauss_param_new_name (mset, name, mu, sigma));

  ncm_likelihood_priors_add (lh, prior);
  ncm_prior_clear (&prior);
}

/**
 * ncm_likelihood_priors_add_gauss_func:
 * @lh: a #NcmLikelihood
 * @mean_func: FIXME
 * @mu: FIXME
 * @sigma: FIXME
 * @var: FIXME
 *
 * FIXME
 *
 */
void
ncm_likelihood_priors_add_gauss_func (NcmLikelihood *lh, NcmMSetFunc *mean_func, gdouble mu, gdouble sigma, gdouble var)
{
  NcmPrior *prior = NCM_PRIOR (ncm_prior_gauss_func_new (mean_func, mu, sigma, var));

  ncm_likelihood_priors_add (lh, prior);
  ncm_prior_clear (&prior);
}

/**
 * ncm_likelihood_priors_add_flat_param:
 * @lh: a #NcmLikelihood
 * @mid: FIXME
 * @pid: FIXME
 * @x_low: FIXME
 * @x_upp: FIXME
 * @scale: FIXME
 *
 * FIXME
 *
 */
void
ncm_likelihood_priors_add_flat_param (NcmLikelihood *lh, NcmModelID mid, guint pid, gdouble x_low, gdouble x_upp, gdouble scale)
{
  NcmPrior *prior = NCM_PRIOR (ncm_prior_flat_param_new (mid, pid, x_low, x_upp, scale));

  ncm_likelihood_priors_add (lh, prior);
  ncm_prior_clear (&prior);
}

/**
 * ncm_likelihood_priors_add_flat_param_pindex:
 * @lh: a #NcmLikelihood
 * @pi: FIXME
 * @x_low: FIXME
 * @x_upp: FIXME
 * @scale: FIXME
 *
 * FIXME
 *
 */
void
ncm_likelihood_priors_add_flat_param_pindex (NcmLikelihood *lh, const NcmMSetPIndex *pi, gdouble x_low, gdouble x_upp, gdouble scale)
{
  NcmPrior *prior = NCM_PRIOR (ncm_prior_flat_param_new_pindex (pi, x_low, x_upp, scale));

  ncm_likelihood_priors_add (lh, prior);
  ncm_prior_clear (&prior);
}

/**
 * ncm_likelihood_priors_add_flat_param_name:
 * @lh: a #NcmLikelihood
 * @mset: FIXME
 * @name: FIXME
 * @x_low: FIXME
 * @x_upp: FIXME
 * @scale: FIXME
 *
 * FIXME
 *
 */
void
ncm_likelihood_priors_add_flat_param_name (NcmLikelihood *lh, NcmMSet *mset, const gchar *name, gdouble x_low, gdouble x_upp, gdouble scale)
{
  NcmPrior *prior = NCM_PRIOR (ncm_prior_flat_param_new_name (mset, name, x_low, x_upp, scale));

  ncm_likelihood_priors_add (lh, prior);
  ncm_prior_clear (&prior);
}

/**
 * ncm_likelihood_priors_add_flat_func:
 * @lh: a #NcmLikelihood
 * @mean_func: FIXME
 * @x_low: FIXME
 * @x_upp: FIXME
 * @scale: FIXME
 * @variable: FIXME
 *
 * FIXME
 *
 */
void
ncm_likelihood_priors_add_flat_func (NcmLikelihood *lh, NcmMSetFunc *mean_func, gdouble x_low, gdouble x_upp, gdouble scale, gdouble variable)
{
  NcmPrior *prior = NCM_PRIOR (ncm_prior_flat_func_new (mean_func, x_low, x_upp, scale, variable));

  ncm_likelihood_priors_add (lh, prior);
  ncm_prior_clear (&prior);
}

/**
 * ncm_likelihood_priors_peek_f:
 * @lh: a #NcmLikelihood
 * @i: FIXME
 *
 * FIXME
 *
 * Returns: (transfer none): FIXME
 */
NcmPrior *
ncm_likelihood_priors_peek_f (NcmLikelihood *lh, guint i)
{
  return g_ptr_array_index (lh->priors_f, i);
}

/**
 * ncm_likelihood_priors_peek_m2lnL:
 * @lh: a #NcmLikelihood
 * @i: FIXME
 *
 * FIXME
 *
 * Returns: (transfer none): FIXME
 */
NcmPrior *
ncm_likelihood_priors_peek_m2lnL (NcmLikelihood *lh, guint i)
{
  return g_ptr_array_index (lh->priors_m2lnL, i);
}

/**
 * ncm_likelihood_priors_length_f:
 * @lh: a #NcmLikelihood
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
 * @lh: a #NcmLikelihood
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
  guint i                  = 0;
  gboolean has_prior_m2lnL = ncm_likelihood_priors_length_m2lnL (lh) != 0;

  if (has_prior_m2lnL)
    g_error ("ncm_likelihood_priors_leastsquares_f: cannot calculate leastsquares f, the likelihood contains m2lnL priors.");

  for (i = 0; i < lh->priors_f->len; i++)
  {
    NcmPrior *prior         = NCM_PRIOR (ncm_obj_array_peek (lh->priors_f, i));
    const gdouble prior_val = ncm_mset_func_eval0 (NCM_MSET_FUNC (prior), mset);

    ncm_vector_set (priors_f, i, prior_val);
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
  guint data_size          = ncm_dataset_get_n (lh->dset);
  guint priors_f_size      = ncm_likelihood_priors_length_f (lh);
  gboolean has_prior_m2lnL = ncm_likelihood_priors_length_m2lnL (lh) != 0;

  if (has_prior_m2lnL)
    g_error ("ncm_likelihood_priors_leastsquares_f: cannot calculate leastsquares f, the likelihood contains m2lnL priors.");

  if (data_size)
  {
    ncm_vector_get_subvector2 (lh->ls_f_data, f, 0, data_size);

    ncm_dataset_leastsquares_f (lh->dset, mset, lh->ls_f_data);
  }

  if (priors_f_size)
  {
    ncm_vector_get_subvector2 (lh->ls_f_priors, f, data_size, priors_f_size);

    ncm_likelihood_priors_leastsquares_f (lh, mset, lh->ls_f_priors);
  }

  return;
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
    NcmPrior *prior         = NCM_PRIOR (ncm_obj_array_peek (lh->priors_f, i));
    const gdouble prior_val = ncm_mset_func_eval0 (NCM_MSET_FUNC (prior), mset);

    *priors_m2lnL += prior_val * prior_val;
  }

  for (i = 0; i < lh->priors_m2lnL->len; i++)
  {
    NcmPrior *prior         = NCM_PRIOR (ncm_obj_array_peek (lh->priors_m2lnL, i));
    const gdouble prior_val = ncm_mset_func_eval0 (NCM_MSET_FUNC (prior), mset);

    *priors_m2lnL += prior_val;
  }

  return;
}

/**
 * ncm_likelihood_priors_m2lnL_vec:
 * @lh: a #NcmLikelihood.
 * @mset: a #NcmMSet.
 * @priors_m2lnL_v: a #NcmVector
 *
 * FIXME
 */
void
ncm_likelihood_priors_m2lnL_vec (NcmLikelihood *lh, NcmMSet *mset, NcmVector *priors_m2lnL_v)
{
  guint i = 0, j = 0;

  g_assert_cmpuint (ncm_vector_len (priors_m2lnL_v), >=, lh->priors_f->len + lh->priors_m2lnL->len);

  for (i = 0; i < lh->priors_f->len; i++)
  {
    NcmPrior *prior         = NCM_PRIOR (g_ptr_array_index (lh->priors_f, i));
    const gdouble prior_val = ncm_mset_func_eval0 (NCM_MSET_FUNC (prior), mset);

    ncm_vector_set (priors_m2lnL_v, j, prior_val * prior_val);
    j++;
  }

  for (i = 0; i < lh->priors_m2lnL->len; i++)
  {
    NcmPrior *prior         = NCM_PRIOR (g_ptr_array_index (lh->priors_m2lnL, i));
    const gdouble prior_val = ncm_mset_func_eval0 (NCM_MSET_FUNC (prior), mset);

    ncm_vector_set (priors_m2lnL_v, j, prior_val);
    j++;
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
  const guint data_length  = ncm_dataset_get_length (lh->dset);
  const guint prior_length = lh->priors_f->len + lh->priors_m2lnL->len;
  const guint v_size       = prior_length + data_length;

  if (lh->m2lnL_v != NULL)
  {
    if (ncm_vector_len (lh->m2lnL_v) != v_size)
    {
      ncm_vector_clear (&lh->m2lnL_v);
      lh->m2lnL_v = ncm_vector_new (v_size);
    }
  }
  else
  {
    lh->m2lnL_v = ncm_vector_new (v_size);
  }

  if (prior_length > 0)
  {
    NcmVector *priors_m2lnL_v = ncm_vector_get_subvector (lh->m2lnL_v, data_length, prior_length);

    ncm_dataset_m2lnL_vec (lh->dset, mset, lh->m2lnL_v);
    ncm_likelihood_priors_m2lnL_vec (lh, mset, priors_m2lnL_v);
    ncm_vector_free (priors_m2lnL_v);
  }
  else
  {
    ncm_dataset_m2lnL_vec (lh->dset, mset, lh->m2lnL_v);
  }

  /*ncm_vector_log_vals (lh->m2lnL_v, "m2lnL: ", "% 22.15g", TRUE);*/
  *m2lnL = ncm_vector_sum_cpts (lh->m2lnL_v);

  return;
}


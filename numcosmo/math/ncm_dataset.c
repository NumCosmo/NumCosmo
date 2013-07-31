/***************************************************************************
 *            ncm_dataset.c
 *
 *  Tue May 29 19:28:48 2007
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
 * SECTION:ncm_dataset
 * @title: Data Set
 * @short_description: Object representing a set of NcmData objects
 *
 * FIXME
 * 
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_dataset.h"
#include "math/ncm_cfg.h"

G_DEFINE_TYPE (NcmDataset, ncm_dataset, G_TYPE_OBJECT);

#define _NCM_DATASET_INITIAL_ALLOC 10

static void
ncm_dataset_init (NcmDataset *dset)
{
  dset->bstype    = NCM_DATASET_BSTRAP_DISABLE;
  dset->data      = g_ptr_array_sized_new (_NCM_DATASET_INITIAL_ALLOC);
  dset->data_prob = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), _NCM_DATASET_INITIAL_ALLOC);
  dset->bstrap    = g_array_sized_new (FALSE, FALSE, sizeof (guint), _NCM_DATASET_INITIAL_ALLOC);
  g_ptr_array_set_free_func (dset->data, (GDestroyNotify) &ncm_data_free);
}

static void
ncm_dataset_dispose (GObject *object)
{
  NcmDataset *dset = NCM_DATASET (object);

  if (dset->data != NULL)
  {
    g_ptr_array_unref (dset->data);
    dset->data = NULL;
  }
  if (dset->data_prob != NULL)
  {
    g_array_unref (dset->data_prob);
    dset->data_prob = NULL;
  }
  if (dset->bstrap != NULL)
  {
    g_array_unref (dset->bstrap);
    dset->bstrap = NULL;
  }

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_dataset_parent_class)->dispose (object);
}

static void
ncm_dataset_finalize (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_dataset_parent_class)->finalize (object);
}

static void
ncm_dataset_class_init (NcmDatasetClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);

  object_class->dispose  = ncm_dataset_dispose;
  object_class->finalize = ncm_dataset_finalize;
}

/**
 * ncm_dataset_new:
 *
 * Creates a new empty #NcmDataset object. 
 *
 * Returns: a new #NcmDataset.
 */
NcmDataset *
ncm_dataset_new (void)
{
  NcmDataset *dset = g_object_new (NCM_TYPE_DATASET, NULL);
  return dset;
}

/**
 * ncm_dataset_ref:
 * @dset: pointer to type defined by #NcmDataset
 *
 * Increases the reference count of @dset by one.
 *
 * Returns: (transfer full): @dset.
 */
NcmDataset *
ncm_dataset_ref (NcmDataset *dset)
{
  return g_object_ref (dset);
}

/**
 * ncm_dataset_dup:
 * @dset: pointer to type defined by #NcmDataset
 *
 * Duplicates the object and all of its content.
 *
 * Returns: (transfer full): the duplicate of @dset.
 */
NcmDataset *
ncm_dataset_dup (NcmDataset *dset)
{
  NcmDataset *dset_dup = ncm_dataset_new ();
  gint i;

  for (i = 0; i < dset->data->len; i++)
  {
    NcmData *data = ncm_dataset_peek_data (dset, i);
    NcmData *data_dup = ncm_data_dup (data);
    g_ptr_array_add (dset_dup->data, data_dup);
  }

  return dset_dup;
}

/**
 * ncm_dataset_copy:
 * @dset: pointer to type defined by #NcmDataset
 *
 * Duplicates the object getting a reference of its content.
 *
 * Returns: (transfer full): the duplicate of @dset, new container.
 */
NcmDataset *
ncm_dataset_copy (NcmDataset *dset)
{
  NcmDataset *dset_dup = ncm_dataset_new ();
  gint i;

  for (i = 0; i < dset->data->len; i++)
  {
    NcmData *data = ncm_dataset_peek_data (dset, i);
    NcmData *data_ref = ncm_data_ref (data);
    g_ptr_array_add (dset_dup->data, data_ref);
  }

  return dset_dup;
}

static void
_ncm_dataset_update_bstrap (NcmDataset *dset)
{
  if (dset->bstype == NCM_DATASET_BSTRAP_TOTAL)
  {
    guint n = ncm_dataset_get_n (dset);
    guint i;
    
    g_array_set_size (dset->data_prob, dset->data->len);
    g_array_set_size (dset->bstrap, dset->data->len);

    for (i = 0; i < dset->data->len; i++)
    {
      NcmData *data = ncm_dataset_peek_data (dset, i);
      gdouble p_i = ncm_data_get_length (data) * 1.0 / n;
      g_array_index (dset->data_prob, gdouble, i) = p_i;
    }
  }
}

/**
 * ncm_dataset_append_data:
 * @dset: pointer to type defined by #NcmDataset
 * @data: #NcmData object to be appended to #NcmDataset
 *
 * Appends @data to @dset.
 *
 */
void
ncm_dataset_append_data (NcmDataset *dset, NcmData *data)
{
  gboolean enable = (dset->bstype != NCM_DATASET_BSTRAP_DISABLE) ? TRUE : FALSE;

  data = ncm_data_ref (data);
  g_ptr_array_add (dset->data, data);

  ncm_data_bootstrap_set (data, enable);

  _ncm_dataset_update_bstrap (dset);
}

/**
 * ncm_dataset_get_n:
 * @dset: pointer to type defined by #NcmDataset
 *
 * Calculates the total number of data set points.
 *
 * Returns: total number of data set points.
 */
guint
ncm_dataset_get_n (NcmDataset *dset)
{
  gint i;
  guint n = 0;

  for (i = 0; i < dset->data->len; i++)
  {
    NcmData *data = ncm_dataset_peek_data (dset, i);
    n += ncm_data_get_length (data);
  }
  
  return n;
}

/**
 * ncm_dataset_get_dof:
 * @dset: pointer to type defined by #NcmDataset
 *
 * Calculate the total degrees of freedom associated with all #NcmData
 * objects.
 *
 * Returns: summed degrees of freedom of all #NcmData in @dset. 
 */
guint
ncm_dataset_get_dof (NcmDataset *dset)
{
  gint i;
  guint dof = 0;

  for (i = 0; i < dset->data->len; i++)
  {
    NcmData *data = ncm_dataset_peek_data (dset, i);
    dof += ncm_data_get_dof (data);
  }
  
  return dof;
}


/**
 * ncm_dataset_all_init:
 * @dset: pointer to type defined by #NcmDataset
 *
 * Checks whenever all #NcmData in @dset are initiated.
 *
 * Returns: whenever @dset is initiated.
 */
gboolean
ncm_dataset_all_init (NcmDataset *dset)
{
  gint i;

  for (i = 0; i < dset->data->len; i++)
  {
    NcmData *data = ncm_dataset_peek_data (dset, i);
    if (!data->init)
      return FALSE;
  }
  return TRUE;
}

/**
 * ncm_dataset_get_length:
 * @dset: pointer to type defined by #NcmDataset
 *
 * Number of diferent #NcmData in @dset.
 *
 * Returns: number of #NcmData objects in the set
 */
guint
ncm_dataset_get_length (NcmDataset *dset)
{
  return dset->data->len;
}

/**
 * ncm_dataset_get_data:
 * @dset: pointer to type defined by #NcmDataset
 * @n: the #NcmData index.
 *
 * Gets the @n-th #NcmData in @dset and increses its reference count by one.
 *
 * Returns: (transfer full): the #NcmData associated with @n.
 */
NcmData *
ncm_dataset_get_data (NcmDataset *dset, guint n)
{
  return ncm_data_ref (ncm_dataset_peek_data (dset, n));
}

/**
 * ncm_dataset_peek_data:
 * @dset: pointer to type defined by #NcmDataset
 * @n: the #NcmData index.
 *
 * Gets the @n-th #NcmData in @dset.
 *
 * Returns: (transfer none): the #NcmData associated with @n.
 */
NcmData *
ncm_dataset_peek_data (NcmDataset *dset, guint n)
{
  g_assert (n < dset->data->len);
  return NCM_DATA (g_ptr_array_index (dset->data, n));
}

/**
 * ncm_dataset_free:
 * @dset: pointer to type defined by #NcmDataset
 *
 * Decreses the reference count of @dset by one.
 */
void
ncm_dataset_free (NcmDataset *dset)
{
  g_object_unref (dset);
}

/**
 * ncm_dataset_clear:
 * @dset: pointer to type defined by #NcmDataset
 *
 * Decreses the reference count of *@dset by one, and sets *@dset to NULL.
 */
void
ncm_dataset_clear (NcmDataset **dset)
{
  g_clear_object (dset);
}

/**
 * ncm_dataset_resample:
 * @dset: a #NcmDataset.
 * @mset: a #NcmMSet.
 *
 * Resamples every #NcmData in @dset with the models contained in @mset.
 *
 */
void
ncm_dataset_resample (NcmDataset *dset, NcmMSet *mset)
{
  gint i;

  for (i = 0; i < dset->data->len; i++)
  {
    NcmData *data = ncm_dataset_peek_data (dset, i);
    ncm_data_resample (data, mset);
  }
}

/**
 * ncm_dataset_bootstrap_set:
 * @dset: a #NcmDataset.
 * @bstype: a #NcmDatasetBStrapType.
 *
 * Disable or sets bootstrap method for @dset.
 *
 */
void
ncm_dataset_bootstrap_set (NcmDataset *dset, NcmDatasetBStrapType bstype)
{
  gint i;
  gboolean enable = (bstype != NCM_DATASET_BSTRAP_DISABLE) ? TRUE : FALSE;
  
  for (i = 0; i < dset->data->len; i++)
  {
    NcmData *data = ncm_dataset_peek_data (dset, i);
    ncm_data_bootstrap_set (data, enable);
  }

  _ncm_dataset_update_bstrap (dset);
}

/**
 * ncm_dataset_bootstrap_resample:
 * @dset: a #NcmDataset.
 *
 * Perform one bootstrap as in ncm_data_bootstrap_resample() in every #NcmData 
 * in @dset.
 *
 */
void
ncm_dataset_bootstrap_resample (NcmDataset *dset)
{
  gint i;
  switch (dset->bstype)
  {
    case NCM_DATASET_BSTRAP_PARTIAL:
    {
      for (i = 0; i < dset->data->len; i++)
      {
        NcmData *data = ncm_dataset_peek_data (dset, i);
        ncm_bootstrap_set_bsize (data->bstrap, data->bstrap->fsize);
        ncm_data_bootstrap_resample (data);
      }
      break;
    }
    case NCM_DATASET_BSTRAP_TOTAL:
    {
      NcmRNG *rng = ncm_rng_pool_get (NCM_BOOTSTRAP_RNG_NAME);
      guint n = ncm_dataset_get_n (dset);
      ncm_rng_lock (rng);
      gsl_ran_multinomial (rng->r, dset->data->len, n, 
                           (gdouble *)dset->data_prob->data, 
                           (guint *)dset->bstrap->data);
      ncm_rng_unlock (rng);
      ncm_rng_free (rng);
      
      for (i = 0; i < dset->data->len; i++)
      {
        NcmData *data = ncm_dataset_peek_data (dset, i);
        guint bsize = g_array_index (dset->bstrap, guint, i);
        ncm_bootstrap_set_bsize (data->bstrap, bsize);
        if (bsize > 0)
          ncm_data_bootstrap_resample (data);
      }
      break;
    }
    default:
      g_error ("ncm_dataset_bootstrap_resample: bootstrap is disabled.");
      break;
  }
}

/**
 * ncm_dataset_log_info:
 * @dset: a #NcmDataset
 *
 * Prints in the log the informations associated with every #NcmData in @dset.
 * 
 */
void
ncm_dataset_log_info (NcmDataset *dset)
{
  gint i;

  ncm_cfg_msg_sepa ();
  g_message ("# Data used:\n");
  for (i = 0; i < dset->data->len; i++)
  {
    NcmData *data = ncm_dataset_peek_data (dset, i);
    gchar *desc = ncm_data_get_desc (data);
    ncm_message_ww (desc, 
                    "#   - ", 
                    "#       ", 
                    80);
  }

  return;
}

/**
 * ncm_dataset_has_leastsquares_f:
 * @dset: a #NcmDataset
 *
 * FIXME
 * 
 * Returns: FIXME
 */
gboolean 
ncm_dataset_has_leastsquares_f (NcmDataset *dset)
{
  if (dset->data->len == 0)
    return FALSE;
  else
  {
    gint i;
    for (i = 0; i < dset->data->len; i++)
    {
      NcmData *data = NCM_DATA (g_ptr_array_index (dset->data, i));
      if (!NCM_DATA_GET_CLASS (data)->leastsquares_f)
        return FALSE;
    }
  }
  return TRUE;
}

/**
 * ncm_dataset_has_leastsquares_J:
 * @dset: a #NcmDataset
 *
 * FIXME
 * 
 * Returns: FIXME
 */
gboolean 
ncm_dataset_has_leastsquares_J (NcmDataset *dset)
{
  if (dset->data->len == 0)
    return FALSE;
  else
  {
    gint i;
    for (i = 0; i < dset->data->len; i++)
    {
      NcmData *data = NCM_DATA (g_ptr_array_index (dset->data, i));
      if (!NCM_DATA_GET_CLASS (data)->leastsquares_J)
        return FALSE;
    }
  }
  return TRUE;
}

/**
 * ncm_dataset_has_leastsquares_f_J:
 * @dset: a #NcmDataset
 *
 * FIXME
 * 
 * Returns: FIXME
 */
gboolean 
ncm_dataset_has_leastsquares_f_J (NcmDataset *dset)
{
  if (dset->data->len == 0)
    return FALSE;
  else
  {
    gint i;
    for (i = 0; i < dset->data->len; i++)
    {
      NcmData *data = NCM_DATA (g_ptr_array_index (dset->data, i));
      if (!NCM_DATA_GET_CLASS (data)->leastsquares_f_J)
        return FALSE;
    }
  }
  return TRUE;
}

/**
 * ncm_dataset_has_m2lnL_val:
 * @dset: a #NcmDataset
 *
 * FIXME
 * 
 * Returns: FIXME
 * 
 */
gboolean 
ncm_dataset_has_m2lnL_val (NcmDataset *dset)
{
  if (dset->data->len == 0)
    return FALSE;
  else
  {
    gint i;
    for (i = 0; i < dset->data->len; i++)
    {
      NcmData *data = NCM_DATA (g_ptr_array_index (dset->data, i));
      if (!NCM_DATA_GET_CLASS (data)->m2lnL_val)
        return FALSE;
    }
  }
  return TRUE;
}

/**
 * ncm_dataset_has_m2lnL_grad:
 * @dset: a #NcmDataset
 *
 * FIXME
 * 
 * Returns: FIXME
 * 
 */
gboolean 
ncm_dataset_has_m2lnL_grad (NcmDataset *dset)
{
  if (dset->data->len == 0)
    return FALSE;
  else
  {
    gint i;
    for (i = 0; i < dset->data->len; i++)
    {
      NcmData *data = NCM_DATA (g_ptr_array_index (dset->data, i));
      if (!NCM_DATA_GET_CLASS (data)->m2lnL_grad)
        return FALSE;
    }
  }
  return TRUE;
}

/**
 * ncm_dataset_has_m2lnL_val_grad:
 * @dset: a #NcmDataset
 *
 * FIXME
 * 
 * Returns: FIXME
 * 
 */
gboolean 
ncm_dataset_has_m2lnL_val_grad (NcmDataset *dset)
{
  if (dset->data->len == 0)
    return FALSE;
  else
  {
    gint i;
    for (i = 0; i < dset->data->len; i++)
    {
      NcmData *data = NCM_DATA (g_ptr_array_index (dset->data, i));
      if (!NCM_DATA_GET_CLASS (data)->m2lnL_val_grad)
        return FALSE;
    }
  }
  return TRUE;
}

/**
 * ncm_dataset_data_leastsquares_f:
 * @dset: a #NcmLikelihood.
 * @mset: a #NcmMSet.
 * @f: a #NcmVector.
 *
 * FIXME
 */
void
ncm_dataset_leastsquares_f (NcmDataset *dset, NcmMSet *mset, NcmVector *f)
{
  guint pos = 0, i;

  for (i = 0; i < dset->data->len; i++)
  {
    NcmData *data = NCM_DATA (g_ptr_array_index (dset->data, i));
    guint n = ncm_data_get_length (data);

    if (!NCM_DATA_GET_CLASS (data)->leastsquares_f)
      g_error ("ncm_dataset_leastsquares_f: %s dont implement leastsquares vector f", G_OBJECT_TYPE_NAME (data));
    else
    {
      NcmVector *f_i = ncm_vector_get_subvector (f, pos, n);
      ncm_data_prepare (data, mset);
      NCM_DATA_GET_CLASS (data)->leastsquares_f (data, mset, f_i);
      pos += n;
      ncm_vector_free (f_i);
    }
  }
  
  return;
}

/**
 * ncm_dataset_leastsquares_J:
 * @dset: a #NcmLikelihood.
 * @mset: a #NcmMSet.
 * @J: a #NcmMatrix.
 *
 * FIXME
 * 
 */
void
ncm_dataset_leastsquares_J (NcmDataset *dset, NcmMSet *mset, NcmMatrix *J)
{
  gint pos = 0, i;

  for (i = 0; i < dset->data->len; i++)
  {
    NcmData *data = NCM_DATA (g_ptr_array_index (dset->data, i));
    guint n = ncm_data_get_length (data);

    if (!NCM_DATA_GET_CLASS (data)->leastsquares_J)
      g_error ("ncm_dataset_leastsquares_J: %s dont implement leastsquares matrix J", G_OBJECT_TYPE_NAME (data));
    else
    {
      NcmMatrix *J_i = ncm_matrix_get_submatrix (J, pos, 0, n, NCM_MATRIX_NCOLS (J));
      ncm_data_prepare (data, mset);
      
      NCM_DATA_GET_CLASS (data)->leastsquares_J (data, mset, J_i);
      
      pos += n;
      ncm_matrix_free (J_i);
    }
  }
  
  return;
}

/**
 * ncm_dataset_leastsquares_f_J:
 * @dset: a #NcmLikelihood.
 * @mset: a #NcmMSet.
 * @f: a #NcmVector.
 * @J: a #NcmMatrix.
 *
 * FIXME
 */
void
ncm_dataset_leastsquares_f_J (NcmDataset *dset, NcmMSet *mset, NcmVector *f, NcmMatrix *J)
{
  gint pos = 0, i;

  for (i = 0; i < dset->data->len; i++)
  {
    NcmData *data = NCM_DATA (g_ptr_array_index (dset->data, i));
    guint n = ncm_data_get_length (data);
    NcmMatrix *J_i = ncm_matrix_get_submatrix (J, pos, 0, n, NCM_MATRIX_NCOLS (J));
    NcmVector *f_i = ncm_vector_get_subvector (f, pos, n);
    ncm_data_prepare (data, mset);

    if (NCM_DATA_GET_CLASS (data)->leastsquares_f_J != NULL)
      NCM_DATA_GET_CLASS (data)->leastsquares_f_J (data, mset, f_i, J_i);
    else if (NCM_DATA_GET_CLASS (data)->leastsquares_f != NULL && NCM_DATA_GET_CLASS (data)->leastsquares_J != NULL)
    {
      NCM_DATA_GET_CLASS (data)->leastsquares_f (data, mset, f_i);
      NCM_DATA_GET_CLASS (data)->leastsquares_J (data, mset, J_i);
    }
    else
      g_error ("ncm_dataset_leastsquares_f_J: %s dont implement leastsquares f J", G_OBJECT_TYPE_NAME (data));

    pos += n;
    ncm_matrix_free (J_i);
    ncm_vector_free (f_i);
  }
}

/**
 * ncm_dataset_m2lnL_val:
 * @dset: a #NcmLikelihood.
 * @mset: a #NcmMSet.
 * @m2lnL: (out): FIXME
 * 
 * FIXME
 */
void
ncm_dataset_m2lnL_val (NcmDataset *dset, NcmMSet *mset, gdouble *m2lnL)
{
  gint i;
  *m2lnL = 0.0;

  for (i = 0; i < dset->data->len; i++)
  {
    NcmData *data = NCM_DATA (g_ptr_array_index (dset->data, i));

    if (!NCM_DATA_GET_CLASS (data)->m2lnL_val)
      g_error ("ncm_dataset_m2lnL_val: %s dont implement m2lnL", G_OBJECT_TYPE_NAME (data));
    else
    {
      gdouble m2lnL_i;
      ncm_data_prepare (data, mset);
      NCM_DATA_GET_CLASS (data)->m2lnL_val (data, mset, &m2lnL_i);
      *m2lnL += m2lnL_i;
    }
  }

  return;
}

/**
 * ncm_dataset_m2lnL_grad:
 * @dset: a #NcmLikelihood.
 * @mset: a #NcmMSet.
 * @grad: a #NcmVector.
 *
 * FIXME
 */
void
ncm_dataset_m2lnL_grad (NcmDataset *dset, NcmMSet *mset, NcmVector *grad)
{
  gint i;
  guint free_params_len = ncm_mset_fparams_len (mset);
  NcmVector *grad_i = ncm_vector_new (free_params_len);

  ncm_vector_set_zero (grad);

  for (i = 0; i < dset->data->len; i++)
  {
    NcmData *data = NCM_DATA (g_ptr_array_index (dset->data, i));

    if (!NCM_DATA_GET_CLASS (data)->m2lnL_grad)
      g_error ("ncm_dataset_m2lnL_grad: %s dont implement m2lnL grad", G_OBJECT_TYPE_NAME (data));
    else
    {
      ncm_data_prepare (data, mset);
      NCM_DATA_GET_CLASS (data)->m2lnL_grad (data, mset, grad_i);
      ncm_vector_add (grad, grad_i);
    }
  }

  ncm_vector_free (grad_i);

  return;
}

/**
 * ncm_dataset_m2lnL_val_grad:
 * @dset: a #NcmLikelihood.
 * @mset: a #NcmMSet.
 * @m2lnL: (out): FIXME
 * @grad: a #NcmVector.
 *
 * FIXME
 */
void
ncm_dataset_m2lnL_val_grad (NcmDataset *dset, NcmMSet *mset, gdouble *m2lnL, NcmVector *grad)
{
  gint i;
  guint free_params_len = ncm_mset_fparams_len (mset);
  NcmVector *grad_i = ncm_vector_new (free_params_len);

  ncm_vector_set_zero (grad);
  *m2lnL = 0.0;

  for (i = 0; i < dset->data->len; i++)
  {
    NcmData *data = NCM_DATA (g_ptr_array_index (dset->data, i));
    gdouble m2lnL_i;

    ncm_data_prepare (data, mset);
    
    if (NCM_DATA_GET_CLASS (data)->m2lnL_val_grad != NULL)
      NCM_DATA_GET_CLASS (data)->m2lnL_val_grad (data, mset, &m2lnL_i, grad_i);
    else if (NCM_DATA_GET_CLASS (data)->m2lnL_val != NULL && NCM_DATA_GET_CLASS (data)->m2lnL_grad != NULL)
    {
      NCM_DATA_GET_CLASS (data)->m2lnL_val (data, mset, &m2lnL_i);
      NCM_DATA_GET_CLASS (data)->m2lnL_grad (data, mset, grad_i);
    }
    else
      g_error ("ncm_dataset_m2lnL_val_grad: %s dont implement m2lnL val grad", G_OBJECT_TYPE_NAME (data));

    *m2lnL += m2lnL_i;
    ncm_vector_add (grad, grad_i);
  }

  ncm_vector_free (grad_i);  
}

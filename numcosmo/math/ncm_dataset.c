/***************************************************************************
 *            ncm_dataset.c
 *
 *  Tue May 29 19:28:48 2007
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
 * SECTION:ncm_dataset
 * @title: NcmDataset
 * @short_description: A set of NcmData objects
 *
 * The purpose of this class is to define a collection of #NcmData objects. These
 * objects serve as containers for #NcmData intended for use within the NumCosmo
 * library. Each individual #NcmData object is responsible for defining a distinct
 * data likelihood function. It is essential to note that all #NcmData objects are
 * designed to be entirely statistically independent of each other.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_dataset.h"
#include "math/ncm_cfg.h"
#include "ncm_enum_types.h"

enum
{
  PROP_0,
  PROP_BSTYPE,
  PROP_OA,
  PROP_SIZE,
};

struct _NcmDataset
{
  /*< private >*/
  GObject parent_instance;
  NcmObjArray *oa;
  NcmDatasetBStrapType bstype;
  GArray *data_prob;
  GArray *bstrap;
};

G_DEFINE_TYPE (NcmDataset, ncm_dataset, G_TYPE_OBJECT);

#define _NCM_DATASET_INITIAL_ALLOC 10

static void
ncm_dataset_init (NcmDataset *dset)
{
  dset->bstype    = NCM_DATASET_BSTRAP_DISABLE;
  dset->oa        = ncm_obj_array_sized_new (_NCM_DATASET_INITIAL_ALLOC);
  dset->data_prob = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), _NCM_DATASET_INITIAL_ALLOC);
  dset->bstrap    = g_array_sized_new (FALSE, FALSE, sizeof (guint), _NCM_DATASET_INITIAL_ALLOC);
}

static void
ncm_dataset_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmDataset *dset = NCM_DATASET (object);

  g_return_if_fail (NCM_IS_DATASET (object));

  switch (prop_id)
  {
    case PROP_BSTYPE:
      ncm_dataset_bootstrap_set (dset, g_value_get_enum (value));
      break;
    case PROP_OA:
      ncm_dataset_set_data_array (dset, (NcmObjArray *) g_value_get_boxed (value));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
ncm_dataset_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmDataset *dset = NCM_DATASET (object);

  g_return_if_fail (NCM_IS_DATASET (object));

  switch (prop_id)
  {
    case PROP_BSTYPE:
      g_value_set_enum (value, dset->bstype);
      break;
    case PROP_OA:
      g_value_set_boxed (value, ncm_dataset_peek_data_array (dset));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
ncm_dataset_dispose (GObject *object)
{
  NcmDataset *dset = NCM_DATASET (object);

  ncm_obj_array_clear (&dset->oa);

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
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = &ncm_dataset_set_property;
  object_class->get_property = &ncm_dataset_get_property;
  object_class->dispose      = &ncm_dataset_dispose;
  object_class->finalize     = &ncm_dataset_finalize;

  /**
   * NcmData:bootstrap-type:
   *
   * Bootstrap method to be used.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_BSTYPE,
                                   g_param_spec_enum ("bootstrap-type",
                                                      NULL,
                                                      "Bootstrap type",
                                                      NCM_TYPE_DATASET_BSTRAP_TYPE, NCM_DATASET_BSTRAP_DISABLE,
                                                      G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmData:data-array:
   *
   * The #NcmData array.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_OA,
                                   g_param_spec_boxed ("data-array",
                                                       NULL,
                                                       "NcmData array",
                                                       NCM_TYPE_OBJ_ARRAY,
                                                       G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
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
 * ncm_dataset_new_list:
 * @data0: first #NcmData to be added.
 * @...: a NULL ended list of #NcmData
 *
 * Creates a new #NcmDataset object and adds a NULL ended list
 * of #NcmData.
 *
 * Returns: a new #NcmDataset.
 */
NcmDataset *
ncm_dataset_new_list (gpointer data0, ...)
{
  va_list ap;
  NcmDataset *dset = ncm_dataset_new ();

  if (data0 != NULL)
  {
    NcmData *data = NULL;

    va_start (ap, data0);

    ncm_dataset_append_data (dset, data0);

    while ((data = va_arg (ap, NcmData *)) != NULL)
      ncm_dataset_append_data (dset, data);

    va_end (ap);
  }

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

static void
_ncm_dataset_update_bstrap (NcmDataset *dset)
{
  if (dset->bstype == NCM_DATASET_BSTRAP_TOTAL)
  {
    guint n = ncm_dataset_get_n (dset);
    guint i;

    g_array_set_size (dset->data_prob, dset->oa->len);
    g_array_set_size (dset->bstrap, dset->oa->len);

    for (i = 0; i < dset->oa->len; i++)
    {
      NcmData *data = ncm_dataset_peek_data (dset, i);
      gdouble p_i   = ncm_data_get_length (data) * 1.0 / n;

      g_array_index (dset->data_prob, gdouble, i) = p_i;
    }
  }
}

/**
 * ncm_dataset_dup:
 * @dset: a #NcmDataset
 * @ser: a #NcmSerialize
 *
 * Duplicates the object and all of its content.
 *
 * Returns: (transfer full): the duplicate of @dset.
 */
NcmDataset *
ncm_dataset_dup (NcmDataset *dset, NcmSerialize *ser)
{
  return NCM_DATASET (ncm_serialize_dup_obj (ser, G_OBJECT (dset)));
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
  guint i;

  for (i = 0; i < dset->oa->len; i++)
  {
    NcmData *data = ncm_dataset_peek_data (dset, i);

    ncm_obj_array_add (dset_dup->oa, G_OBJECT (data));
  }

  return dset_dup;
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

  g_assert (NCM_IS_DATA (data));
  ncm_obj_array_add (dset->oa, G_OBJECT (data));

  if (enable)
    ncm_data_bootstrap_create (data);
  else
    ncm_data_bootstrap_remove (data);

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
  guint i;
  guint n = 0;

  for (i = 0; i < dset->oa->len; i++)
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
  guint i;
  guint dof = 0;

  for (i = 0; i < dset->oa->len; i++)
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
  guint i;

  for (i = 0; i < dset->oa->len; i++)
  {
    NcmData *data = ncm_dataset_peek_data (dset, i);

    if (!ncm_data_is_init (data))
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
  return dset->oa->len;
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
 * @dset: a #NcmDataset
 * @n: the #NcmData index.
 *
 * Gets the @n-th #NcmData in @dset.
 *
 * Returns: (transfer none): the #NcmData associated with @n.
 */
NcmData *
ncm_dataset_peek_data (NcmDataset *dset, guint n)
{
  g_assert_cmpuint (n, <, dset->oa->len);

  return NCM_DATA (ncm_obj_array_peek (dset->oa, n));
}

/**
 * ncm_dataset_get_ndata:
 * @dset: a #NcmDataset
 *
 * Gets number of #NcmData in @dset.
 *
 * Returns: number of #NcmData objects in @dset.
 */
guint
ncm_dataset_get_ndata (NcmDataset *dset)
{
  return dset->oa->len;
}

/**
 * ncm_dataset_set_data_array:
 * @dset: a #NcmDataset
 * @oa: a #NcmObjArray containing #NcmData objects.
 *
 * Sets the @dset with @oa.
 *
 */
void
ncm_dataset_set_data_array (NcmDataset *dset, NcmObjArray *oa)
{
  guint i;
  NcmObjArray *old_oa = dset->oa;

  dset->oa = ncm_obj_array_ref (oa);
  ncm_obj_array_unref (old_oa);

  for (i = 0; i < dset->oa->len; i++)
  {
    NcmData *data = ncm_dataset_peek_data (dset, i);

    if (dset->bstype == NCM_DATASET_BSTRAP_DISABLE)
      ncm_data_bootstrap_remove (data);
    else
      ncm_data_bootstrap_create (data);
  }

  _ncm_dataset_update_bstrap (dset);
}

/**
 * ncm_dataset_peek_data_array:
 * @dset: a #NcmDataset
 *
 * Gets the #NcmObjArray from @dset.
 *
 * Returns: (transfer none): the array of #NcmData.
 */
NcmObjArray *
ncm_dataset_peek_data_array (NcmDataset *dset)
{
  return dset->oa;
}

/**
 * ncm_dataset_get_data_array:
 * @dset: a #NcmDataset
 *
 * Gets the #NcmObjArray from @dset.
 *
 * Returns: (transfer full): the array of #NcmData.
 */
NcmObjArray *
ncm_dataset_get_data_array (NcmDataset *dset)
{
  return ncm_obj_array_ref (ncm_dataset_peek_data_array (dset));
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
 * @dset: a #NcmDataset
 * @mset: a #NcmMSet
 * @rng: a #NcmRNG
 *
 * Resamples every #NcmData in @dset with the models contained in @mset.
 *
 */
void
ncm_dataset_resample (NcmDataset *dset, NcmMSet *mset, NcmRNG *rng)
{
  guint i;

  for (i = 0; i < dset->oa->len; i++)
  {
    NcmData *data = ncm_dataset_peek_data (dset, i);

    ncm_data_resample (data, mset, rng);
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
  if (dset->bstype != bstype)
  {
    guint i;
    gboolean enable = (bstype != NCM_DATASET_BSTRAP_DISABLE) ? TRUE : FALSE;

    dset->bstype = bstype;

    for (i = 0; i < dset->oa->len; i++)
    {
      NcmData *data = ncm_dataset_peek_data (dset, i);

      if (enable)
        ncm_data_bootstrap_create (data);
      else
        ncm_data_bootstrap_remove (data);
    }

    _ncm_dataset_update_bstrap (dset);
  }
}

/**
 * ncm_dataset_bootstrap_resample:
 * @dset: a #NcmDataset.
 * @rng: a #NcmRNG.
 *
 * Perform one bootstrap as in ncm_data_bootstrap_resample() in every #NcmData
 * in @dset.
 *
 */
void
ncm_dataset_bootstrap_resample (NcmDataset *dset, NcmRNG *rng)
{
  guint i;

  switch (dset->bstype)
  {
    case NCM_DATASET_BSTRAP_PARTIAL:
    {
      for (i = 0; i < dset->oa->len; i++)
      {
        NcmData *data        = ncm_dataset_peek_data (dset, i);
        NcmBootstrap *bstrap = ncm_data_peek_bootstrap (data);
        const guint fsize    = ncm_bootstrap_get_fsize (bstrap);

        ncm_bootstrap_set_bsize (bstrap, fsize);
        ncm_data_bootstrap_resample (data, rng);
      }

      break;
    }
    case NCM_DATASET_BSTRAP_TOTAL:
    {
      guint n = ncm_dataset_get_n (dset);

      ncm_rng_lock (rng);
      gsl_ran_multinomial (rng->r, dset->oa->len, n,
                           (gdouble *) dset->data_prob->data,
                           (guint *) dset->bstrap->data);
      ncm_rng_unlock (rng);

      for (i = 0; i < dset->oa->len; i++)
      {
        NcmData *data        = ncm_dataset_peek_data (dset, i);
        NcmBootstrap *bstrap = ncm_data_peek_bootstrap (data);

        guint bsize = g_array_index (dset->bstrap, guint, i);

        ncm_bootstrap_set_bsize (bstrap, bsize);

        if (bsize > 0)
          ncm_data_bootstrap_resample (data, rng);
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
  guint i;

  ncm_cfg_msg_sepa ();
  g_message ("# Data used:\n");

  for (i = 0; i < dset->oa->len; i++)
  {
    NcmData *data     = ncm_dataset_peek_data (dset, i);
    const gchar *desc = ncm_data_peek_desc (data);

    ncm_message_ww (desc,
                    "#   - ",
                    "#       ",
                    80);
  }

  return;
}

/**
 * ncm_dataset_get_info:
 * @dset: a #NcmDataset
 *
 * Obtains the informations associated with every #NcmData in @dset.
 *
 * Returns: (transfer full): @dset description
 */
gchar *
ncm_dataset_get_info (NcmDataset *dset)
{
  guint i;
  GString *desc = g_string_new ("# Data used:\n");

  for (i = 0; i < dset->oa->len; i++)
  {
    NcmData *data       = ncm_dataset_peek_data (dset, i);
    const gchar *desc_i = ncm_data_peek_desc (data);
    gchar *desc_i_ww    = ncm_string_ww (desc_i,
                                         "#   - ",
                                         "#       ",
                                         80);

    g_string_append (desc, desc_i_ww);
    g_free (desc_i_ww);
  }

  return g_string_free (desc, FALSE);
}

/**
 * ncm_dataset_has_leastsquares_f:
 * @dset: a #NcmDataset
 *
 * Whether all the #NcmData in @dset have a ncm_data_leastsquares_f() method.
 *
 * Returns: %TRUE if all the #NcmData in @dset have a ncm_data_leastsquares_f() method.
 */
gboolean
ncm_dataset_has_leastsquares_f (NcmDataset *dset)
{
  if (dset->oa->len == 0)
  {
    return FALSE;
  }
  else
  {
    guint i;

    for (i = 0; i < dset->oa->len; i++)
    {
      NcmData *data = ncm_dataset_peek_data (dset, i);

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
 * Whether all the #NcmData in @dset have a ncm_data_leastsquares_J() method.
 *
 * Returns: %TRUE if all the #NcmData in @dset have a ncm_data_leastsquares_J() method.
 */
gboolean
ncm_dataset_has_leastsquares_J (NcmDataset *dset)
{
  if (dset->oa->len == 0)
  {
    return FALSE;
  }
  else
  {
    guint i;

    for (i = 0; i < dset->oa->len; i++)
    {
      NcmData *data = ncm_dataset_peek_data (dset, i);

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
 * Whether all the #NcmData in @dset have a ncm_data_leastsquares_f_J() method.
 *
 * Returns: %TRUE if all the #NcmData in @dset have a ncm_data_leastsquares_f_J() method.
 */
gboolean
ncm_dataset_has_leastsquares_f_J (NcmDataset *dset)
{
  if (dset->oa->len == 0)
  {
    return FALSE;
  }
  else
  {
    guint i;

    for (i = 0; i < dset->oa->len; i++)
    {
      NcmData *data = ncm_dataset_peek_data (dset, i);

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
 * Whether all the #NcmData in @dset have a ncm_data_m2lnL_val() method.
 *
 * Returns: %TRUE if all the #NcmData in @dset have a ncm_data_m2lnL_val() method.
 */
gboolean
ncm_dataset_has_m2lnL_val (NcmDataset *dset)
{
  if (dset->oa->len == 0)
  {
    return FALSE;
  }
  else
  {
    guint i;

    for (i = 0; i < dset->oa->len; i++)
    {
      NcmData *data = ncm_dataset_peek_data (dset, i);

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
 * Whether all the #NcmData in @dset have a ncm_data_m2lnL_grad() method.
 *
 * Returns: %TRUE if all the #NcmData in @dset have a ncm_data_m2lnL_grad() method.
 *
 */
gboolean
ncm_dataset_has_m2lnL_grad (NcmDataset *dset)
{
  if (dset->oa->len == 0)
  {
    return FALSE;
  }
  else
  {
    guint i;

    for (i = 0; i < dset->oa->len; i++)
    {
      NcmData *data = ncm_dataset_peek_data (dset, i);

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
 * Whether all the #NcmData in @dset have a ncm_data_m2lnL_val_grad() method.
 *
 * Returns: %TRUE if all the #NcmData in @dset have a ncm_data_m2lnL_val_grad() method.
 *
 */
gboolean
ncm_dataset_has_m2lnL_val_grad (NcmDataset *dset)
{
  if (dset->oa->len == 0)
  {
    return FALSE;
  }
  else
  {
    guint i;

    for (i = 0; i < dset->oa->len; i++)
    {
      NcmData *data = ncm_dataset_peek_data (dset, i);

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
 * Computes the leastsquares vector f for the data @data.
 * The vector @f must be allocated with the correct size.
 * The vector @f is filled with the values of the leastsquares vector f.
 *
 */
void
ncm_dataset_leastsquares_f (NcmDataset *dset, NcmMSet *mset, NcmVector *f)
{
  guint pos = 0, i;

  for (i = 0; i < dset->oa->len; i++)
  {
    NcmData *data = ncm_dataset_peek_data (dset, i);
    guint n       = ncm_data_get_length (data);

    if (!NCM_DATA_GET_CLASS (data)->leastsquares_f)
    {
      g_error ("ncm_dataset_leastsquares_f: %s dont implement leastsquares vector f", G_OBJECT_TYPE_NAME (data));
    }
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
 * Computes the leastsquares matrix J for the data @data.
 * The matrix @J must be allocated with the correct size.
 * The matrix @J is filled with the values of the leastsquares matrix J.
 *
 */
void
ncm_dataset_leastsquares_J (NcmDataset *dset, NcmMSet *mset, NcmMatrix *J)
{
  guint pos = 0, i;

  for (i = 0; i < dset->oa->len; i++)
  {
    NcmData *data = ncm_dataset_peek_data (dset, i);
    guint n       = ncm_data_get_length (data);

    if (!NCM_DATA_GET_CLASS (data)->leastsquares_J)
    {
      g_error ("ncm_dataset_leastsquares_J: %s dont implement leastsquares matrix J", G_OBJECT_TYPE_NAME (data));
    }
    else
    {
      NcmMatrix *J_i = ncm_matrix_get_submatrix (J, pos, 0, n, ncm_matrix_ncols (J));

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
 * Computes the leastsquares vector f and matrix J for the data @data.
 * The vector @f must be allocated with the correct size.
 * The matrix @J must be allocated with the correct size.
 * The vector @f is filled with the values of the leastsquares vector f.
 * The matrix @J is filled with the values of the leastsquares matrix J.
 *
 */
void
ncm_dataset_leastsquares_f_J (NcmDataset *dset, NcmMSet *mset, NcmVector *f, NcmMatrix *J)
{
  guint pos = 0, i;

  for (i = 0; i < dset->oa->len; i++)
  {
    NcmData *data  = ncm_dataset_peek_data (dset, i);
    guint n        = ncm_data_get_length (data);
    NcmMatrix *J_i = ncm_matrix_get_submatrix (J, pos, 0, n, ncm_matrix_ncols (J));
    NcmVector *f_i = ncm_vector_get_subvector (f, pos, n);

    ncm_data_prepare (data, mset);

    if (NCM_DATA_GET_CLASS (data)->leastsquares_f_J != NULL)
    {
      NCM_DATA_GET_CLASS (data)->leastsquares_f_J (data, mset, f_i, J_i);
    }
    else if ((NCM_DATA_GET_CLASS (data)->leastsquares_f != NULL) && (NCM_DATA_GET_CLASS (data)->leastsquares_J != NULL))
    {
      NCM_DATA_GET_CLASS (data)->leastsquares_f (data, mset, f_i);
      NCM_DATA_GET_CLASS (data)->leastsquares_J (data, mset, J_i);
    }
    else
    {
      g_error ("ncm_dataset_leastsquares_f_J: %s dont implement leastsquares f J", G_OBJECT_TYPE_NAME (data));
    }

    pos += n;
    ncm_matrix_free (J_i);
    ncm_vector_free (f_i);
  }
}

/**
 * ncm_dataset_m2lnL_val:
 * @dset: a #NcmLikelihood.
 * @mset: a #NcmMSet.
 * @m2lnL: (out): a pointer to a double.
 *
 * Computes the value of the m2lnL for the data @data.
 * The value of the m2lnL is stored in @m2lnL.
 *
 */
void
ncm_dataset_m2lnL_val (NcmDataset *dset, NcmMSet *mset, gdouble *m2lnL)
{
  guint i;

  *m2lnL = 0.0;

  for (i = 0; i < dset->oa->len; i++)
  {
    NcmData *data = ncm_dataset_peek_data (dset, i);

    if (!NCM_DATA_GET_CLASS (data)->m2lnL_val)
    {
      g_error ("ncm_dataset_m2lnL_val: %s dont implement m2lnL", G_OBJECT_TYPE_NAME (data));
    }
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
 * ncm_dataset_m2lnL_vec:
 * @dset: a #NcmLikelihood.
 * @mset: a #NcmMSet.
 * @m2lnL_v: a #NcmVector
 *
 * Computes the value of $-2\ln L$ for every element in the @dset.
 * The values of $-2\ln L$ are stored in the #NcmVector @m2lnL_v.
 *
 */
void
ncm_dataset_m2lnL_vec (NcmDataset *dset, NcmMSet *mset, NcmVector *m2lnL_v)
{
  guint i;

  g_assert_cmpuint (ncm_vector_len (m2lnL_v), >=, dset->oa->len);

  for (i = 0; i < dset->oa->len; i++)
  {
    NcmData *data = ncm_dataset_peek_data (dset, i);

    if (!NCM_DATA_GET_CLASS (data)->m2lnL_val)
    {
      g_error ("ncm_dataset_m2lnL_val: %s dont implement m2lnL", G_OBJECT_TYPE_NAME (data));
    }
    else
    {
      gdouble m2lnL_i;

      ncm_data_prepare (data, mset);
      NCM_DATA_GET_CLASS (data)->m2lnL_val (data, mset, &m2lnL_i);
      ncm_vector_set (m2lnL_v, i, m2lnL_i);
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
 * Computes the gradient of the m2lnL for the data @data.
 * The gradient of the m2lnL is stored in @grad.
 *
 */
void
ncm_dataset_m2lnL_grad (NcmDataset *dset, NcmMSet *mset, NcmVector *grad)
{
  const guint fparams_len = ncm_mset_fparams_len (mset);
  NcmVector *grad_i       = ncm_vector_new (fparams_len);
  guint i;

  ncm_vector_set_zero (grad);

  for (i = 0; i < dset->oa->len; i++)
  {
    NcmData *data = ncm_dataset_peek_data (dset, i);

    if (!NCM_DATA_GET_CLASS (data)->m2lnL_grad)
    {
      g_error ("ncm_dataset_m2lnL_grad: %s dont implement m2lnL grad", G_OBJECT_TYPE_NAME (data));
    }
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
 * @m2lnL: (out): a pointer to a double.
 * @grad: a #NcmVector.
 *
 * Computes the value of the m2lnL for the data @data.
 * The value of the m2lnL is stored in @m2lnL.
 * The gradient of the m2lnL is stored in @grad.
 * The #NcmVector @grad must be previously allocated.
 */
void
ncm_dataset_m2lnL_val_grad (NcmDataset *dset, NcmMSet *mset, gdouble *m2lnL, NcmVector *grad)
{
  const guint fparams_len = ncm_mset_fparams_len (mset);
  NcmVector *grad_i       = ncm_vector_new (fparams_len);
  guint i;

  ncm_vector_set_zero (grad);
  *m2lnL = 0.0;

  for (i = 0; i < dset->oa->len; i++)
  {
    NcmData *data = ncm_dataset_peek_data (dset, i);
    gdouble m2lnL_i;

    ncm_data_prepare (data, mset);

    if (NCM_DATA_GET_CLASS (data)->m2lnL_val_grad != NULL)
    {
      NCM_DATA_GET_CLASS (data)->m2lnL_val_grad (data, mset, &m2lnL_i, grad_i);
    }
    else if ((NCM_DATA_GET_CLASS (data)->m2lnL_val != NULL) && (NCM_DATA_GET_CLASS (data)->m2lnL_grad != NULL))
    {
      NCM_DATA_GET_CLASS (data)->m2lnL_val (data, mset, &m2lnL_i);
      NCM_DATA_GET_CLASS (data)->m2lnL_grad (data, mset, grad_i);
    }
    else
    {
      g_error ("ncm_dataset_m2lnL_val_grad: %s dont implement m2lnL val grad", G_OBJECT_TYPE_NAME (data));
    }

    *m2lnL += m2lnL_i;
    ncm_vector_add (grad, grad_i);
  }

  ncm_vector_free (grad_i);
}

/**
 * ncm_dataset_m2lnL_i_val:
 * @dset: a #NcmLikelihood
 * @mset: a #NcmMSet
 * @i: an integer
 * @m2lnL_i: (out): a pointer to a double
 *
 * Get the value of the @i-th data in the dataset.
 *
 */
void
ncm_dataset_m2lnL_i_val (NcmDataset *dset, NcmMSet *mset, guint i, gdouble *m2lnL_i)
{
  *m2lnL_i = 0.0;

  g_assert_cmpuint (i, <, dset->oa->len);
  {
    NcmData *data = ncm_dataset_peek_data (dset, i);

    if (!NCM_DATA_GET_CLASS (data)->m2lnL_val)
    {
      g_error ("ncm_dataset_m2lnL_val: %s dont implement m2lnL", G_OBJECT_TYPE_NAME (data));
    }
    else
    {
      ncm_data_prepare (data, mset);
      NCM_DATA_GET_CLASS (data)->m2lnL_val (data, mset, m2lnL_i);
    }
  }

  return;
}

/**
 * ncm_dataset_fisher_matrix:
 * @dset: a #NcmDataset
 * @mset: a #NcmMSet
 * @IM: (out) (transfer full): The fisher matrix
 *
 * Calculates the Fisher-information matrix @I adding
 * the individual ones from each #NcmData in @dset.
 * If the #NcmMatrix pointer in *@IM is NULL a new #NcmMatrix
 * will be allocated otherwise *@IM will be used.
 *
 */
void
ncm_dataset_fisher_matrix (NcmDataset *dset, NcmMSet *mset, NcmMatrix **IM)
{
  const guint fparams_len = ncm_mset_fparams_len (mset);
  NcmMatrix *IM0          = ncm_matrix_new (fparams_len, fparams_len);
  guint i;

  *IM = ncm_matrix_new (fparams_len, fparams_len);

  ncm_matrix_set_zero (*IM);

  for (i = 0; i < dset->oa->len; i++)
  {
    NcmData *data = ncm_dataset_peek_data (dset, i);

    ncm_data_fisher_matrix (data, mset, &IM0);

    ncm_matrix_add (*IM, IM0);
  }

  ncm_matrix_free (IM0);
}


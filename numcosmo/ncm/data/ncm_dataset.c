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
 * NcmDataset:
 *
 * A set of statistically independent #NcmData.
 *
 * A #NcmDataset holds the #NcmData of an analysis. They are statistically
 * independent: $-2\ln L$, the Fisher matrix and the bias vector of the set are the
 * sums of those of its members, and the least-squares and mean vectors are their
 * concatenation.
 *
 * Every evaluation prepares all the #NcmData before evaluating any of them. Members
 * that share a resource, as the CMB likelihoods share a Boltzmann solver, state what
 * they need from it in their prepare, so each is evaluated with the requirements of
 * the whole set.
 *
 * The set can also be bootstrapped, see #NcmDatasetBStrapType.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/data/ncm_dataset.h"
#include "ncm/core/ncm_cfg.h"
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
  NcmVector *ls_f;
};

G_DEFINE_TYPE (NcmDataset, ncm_dataset, G_TYPE_OBJECT)

#define _NCM_DATASET_INITIAL_ALLOC 10

static void
ncm_dataset_init (NcmDataset *dset)
{
  dset->bstype    = NCM_DATASET_BSTRAP_DISABLE;
  dset->oa        = ncm_obj_array_sized_new (_NCM_DATASET_INITIAL_ALLOC);
  dset->data_prob = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), _NCM_DATASET_INITIAL_ALLOC);
  dset->bstrap    = g_array_sized_new (FALSE, FALSE, sizeof (guint), _NCM_DATASET_INITIAL_ALLOC);
  dset->ls_f      = ncm_vector_new_data_static (GINT_TO_POINTER (1), 1, 1);
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
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
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
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
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

  ncm_vector_clear (&dset->ls_f);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_dataset_parent_class)->dispose (object);
}

static void
ncm_dataset_class_init (NcmDatasetClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = &ncm_dataset_set_property;
  object_class->get_property = &ncm_dataset_get_property;
  object_class->dispose      = &ncm_dataset_dispose;

  /**
   * NcmDataset:bootstrap-type:
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
   * NcmDataset:data-array:
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
 * Returns: (transfer full): a new #NcmDataset.
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
 * Creates a new #NcmDataset object and adds a %NULL-terminated list of #NcmData.
 *
 * Returns: (transfer full): a new #NcmDataset.
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
 * ncm_dataset_new_array:
 * @data_array: (array length=len) (element-type NcmData): array of #NcmData to be added
 * @len: length of @data_array
 *
 * Creates a new #NcmDataset object and adds @len #NcmData from @data_array.
 *
 * Returns: (transfer full): a new #NcmDataset.
 */
NcmDataset *
ncm_dataset_new_array (NcmData **data_array, guint len)
{
  NcmDataset *dset = ncm_dataset_new ();
  guint i;

  for (i = 0; i < len; i++)
  {
    ncm_dataset_append_data (dset, data_array[i]);
  }

  return dset;
}

/**
 * ncm_dataset_ref:
 * @dset: a #NcmDataset
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
 * @dset: a #NcmDataset
 *
 * Creates a new #NcmDataset holding the same #NcmData as @dset, which are shared,
 * not duplicated, and the same bootstrap type.
 *
 * Returns: (transfer full): the copy of @dset.
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

  dset_dup->bstype = dset->bstype;
  _ncm_dataset_update_bstrap (dset_dup);

  return dset_dup;
}

/**
 * ncm_dataset_append_data:
 * @dset: a #NcmDataset
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
 * @dset: a #NcmDataset
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
 * @dset: a #NcmDataset
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
 * @dset: a #NcmDataset
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
 * @dset: a #NcmDataset
 *
 * Number of different #NcmData in @dset.
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
 * @dset: a #NcmDataset
 * @n: the #NcmData index.
 *
 * Gets the @n-th #NcmData in @dset and increases its reference count by one.
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
 * @dset: a #NcmDataset
 *
 * Decreases the reference count of @dset by one. If the reference count reaches
 * zero, @dset is freed.
 */
void
ncm_dataset_free (NcmDataset *dset)
{
  g_object_unref (dset);
}

/**
 * ncm_dataset_clear:
 * @dset: a #NcmDataset
 *
 * If *@dset is not %NULL, decreases the reference count of *@dset by one and sets
 * *@dset to %NULL.
 */
void
ncm_dataset_clear (NcmDataset **dset)
{
  g_clear_object (dset);
}

static void _ncm_dataset_prepare_all (NcmDataset *dset, NcmMSet *mset);
static void _ncm_dataset_check_fisher_matrix (NcmMatrix **IM, const guint fparams_len);

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

  /* Same reason as in the evaluation paths: a block must not be resampled against a
   * shared resource that the later blocks have not yet placed their requirements on.
   */
  _ncm_dataset_prepare_all (dset, mset);

  for (i = 0; i < dset->oa->len; i++)
  {
    NcmData *data = ncm_dataset_peek_data (dset, i);

    ncm_data_resample (data, mset, rng);
  }
}

/**
 * ncm_dataset_register_shared:
 * @dset: a #NcmDataset
 * @ser: a #NcmSerialize
 *
 * Calls ncm_data_register_shared() on every #NcmData in @dset.
 *
 */
void
ncm_dataset_register_shared (NcmDataset *dset, NcmSerialize *ser)
{
  guint i;

  for (i = 0; i < dset->oa->len; i++)
  {
    NcmData *data = ncm_dataset_peek_data (dset, i);

    ncm_data_register_shared (data, ser);
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
      ncm_rng_multinomial (rng, dset->oa->len, n,
                           (gdouble *) dset->data_prob->data,
                           (guint *) dset->bstrap->data);
      ncm_rng_unlock (rng);

      for (i = 0; i < dset->oa->len; i++)
      {
        NcmData *data        = ncm_dataset_peek_data (dset, i);
        NcmBootstrap *bstrap = ncm_data_peek_bootstrap (data);
        guint bsize          = g_array_index (dset->bstrap, guint, i);

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
 * ncm_dataset_data_leastsquares_f:
 * @dset: a #NcmDataset
 * @mset: a #NcmMSet.
 * @f: a #NcmVector.
 *
 * Computes the least-squares vector of @dset, concatenating those of its #NcmData,
 * into @f, which must have ncm_dataset_get_n() components.
 *
 */

/*
 * Every data block must be prepared before any of them is evaluated. Blocks that share a
 * resource declare what they need from it in their own prepare(): the CMB likelihoods call
 * nc_hipert_boltzmann_require() there, each asking for its own spectra and its own lmax.
 * Preparing and evaluating one block at a time therefore scored the first block against a
 * solution configured with only that block's requirements, and the answer depended on how
 * many blocks had been prepared before it. On the native Planck TT set that shifted the
 * first evaluation of a fresh experiment by 0.0067 in -2lnL, after which every call agreed.
 */
static void
_ncm_dataset_prepare_all (NcmDataset *dset, NcmMSet *mset)
{
  guint i;

  for (i = 0; i < dset->oa->len; i++)
    ncm_data_prepare (ncm_dataset_peek_data (dset, i), mset);
}

void
ncm_dataset_leastsquares_f (NcmDataset *dset, NcmMSet *mset, NcmVector *f)
{
  guint pos = 0, i;

  _ncm_dataset_prepare_all (dset, mset);

  for (i = 0; i < dset->oa->len; i++)
  {
    NcmData *data = ncm_dataset_peek_data (dset, i);
    guint n       = ncm_data_get_length (data);

    if (!NCM_DATA_GET_CLASS (data)->leastsquares_f)
    {
      g_error ("ncm_dataset_leastsquares_f: data `%s' does not implement leastsquares_f.", ncm_data_peek_desc (data));
    }
    else
    {
      ncm_vector_get_subvector2 (dset->ls_f, f, pos, n);

      NCM_DATA_GET_CLASS (data)->leastsquares_f (data, mset, dset->ls_f);
      pos += n;
    }
  }

  return;
}

/**
 * ncm_dataset_m2lnL_val:
 * @dset: a #NcmDataset
 * @mset: a #NcmMSet.
 * @m2lnL: (out): a pointer to a double.
 *
 * Computes $-2\ln L$ of @dset, the sum of those of its #NcmData, and stores it in
 * @m2lnL. Every #NcmData is prepared before any is evaluated.
 *
 */
void
ncm_dataset_m2lnL_val (NcmDataset *dset, NcmMSet *mset, gdouble *m2lnL)
{
  guint i;

  *m2lnL = 0.0;

  _ncm_dataset_prepare_all (dset, mset);

  for (i = 0; i < dset->oa->len; i++)
  {
    NcmData *data = ncm_dataset_peek_data (dset, i);

    if (!NCM_DATA_GET_CLASS (data)->m2lnL_val)
    {
      g_error ("ncm_dataset_m2lnL_val: data `%s' does not implement m2lnL_val.", ncm_data_peek_desc (data));
    }
    else
    {
      gdouble m2lnL_i;

      NCM_DATA_GET_CLASS (data)->m2lnL_val (data, mset, &m2lnL_i);
      *m2lnL += m2lnL_i;
    }
  }

  return;
}

/**
 * ncm_dataset_m2lnL_vec:
 * @dset: a #NcmDataset
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

  _ncm_dataset_prepare_all (dset, mset);

  g_assert_cmpuint (ncm_vector_len (m2lnL_v), >=, dset->oa->len);

  for (i = 0; i < dset->oa->len; i++)
  {
    NcmData *data = ncm_dataset_peek_data (dset, i);

    if (!NCM_DATA_GET_CLASS (data)->m2lnL_val)
    {
      g_error ("ncm_dataset_m2lnL_val: data `%s' does not implement m2lnL_val.", ncm_data_peek_desc (data));
    }
    else
    {
      gdouble m2lnL_i;

      NCM_DATA_GET_CLASS (data)->m2lnL_val (data, mset, &m2lnL_i);
      ncm_vector_set (m2lnL_v, i, m2lnL_i);
    }
  }

  return;
}

/**
 * ncm_dataset_m2lnL_i_val:
 * @dset: a #NcmDataset
 * @mset: a #NcmMSet
 * @i: an integer
 * @m2lnL_i: (out): a pointer to a double
 *
 * Computes $-2\ln L$ of the @i-th #NcmData in @dset; every #NcmData is prepared
 * first.
 *
 */
void
ncm_dataset_m2lnL_i_val (NcmDataset *dset, NcmMSet *mset, guint i, gdouble *m2lnL_i)
{
  *m2lnL_i = 0.0;

  g_assert_cmpuint (i, <, dset->oa->len);

  /* Every block, not only the requested one: see _ncm_dataset_prepare_all. */
  _ncm_dataset_prepare_all (dset, mset);
  {
    NcmData *data = ncm_dataset_peek_data (dset, i);

    if (!NCM_DATA_GET_CLASS (data)->m2lnL_val)
      g_error ("ncm_dataset_m2lnL_val: data `%s' does not implement m2lnL_val.", ncm_data_peek_desc (data));
    else
      NCM_DATA_GET_CLASS (data)->m2lnL_val (data, mset, m2lnL_i);
  }

  return;
}

/**
 * ncm_dataset_has_mean_vector:
 * @dset: a #NcmDataset
 *
 * Whether all the #NcmData in @dset have a ncm_data_mean_vector() method.
 *
 * Returns: %TRUE if all the #NcmData in @dset have a ncm_data_mean_vector() method.
 */
gboolean
ncm_dataset_has_mean_vector (NcmDataset *dset)
{
  guint i;

  for (i = 0; i < dset->oa->len; i++)
  {
    NcmData *data = ncm_dataset_peek_data (dset, i);

    if (!ncm_data_has_mean_vector (data))
      return FALSE;
  }

  return TRUE;
}

/**
 * ncm_dataset_mean_vector:
 * @dset: a #NcmDataset
 * @mset: a #NcmMSet
 * @mu: a #NcmVector
 *
 * Calculates the mean vector @mu, of ncm_dataset_get_n() components, concatenating
 * the individual ones from each #NcmData in @dset. Every #NcmData is prepared before
 * any mean is evaluated, see ncm_dataset_m2lnL_val().
 *
 */
void
ncm_dataset_mean_vector (NcmDataset *dset, NcmMSet *mset, NcmVector *mu)
{
  const guint total_n = ncm_dataset_get_n (dset);
  guint pos           = 0;
  guint i;

  if (ncm_vector_len (mu) != total_n)
    g_error ("ncm_dataset_mean_vector: the dataset has %u points, but the vector has %u.",
             total_n, ncm_vector_len (mu));

  _ncm_dataset_prepare_all (dset, mset);

  for (i = 0; i < dset->oa->len; i++)
  {
    NcmData *data = ncm_dataset_peek_data (dset, i);
    guint n       = ncm_data_get_length (data);

    ncm_vector_get_subvector2 (dset->ls_f, mu, pos, n);

    ncm_data_mean_vector (data, mset, dset->ls_f);
    pos += n;
  }
}

/**
 * ncm_dataset_fisher_matrix:
 * @dset: a #NcmDataset
 * @mset: a #NcmMSet
 * @IM: (inout) (allow-none) (transfer full): the Fisher matrix
 *
 * Calculates the Fisher-information matrix of @dset, the sum of those of its
 * #NcmData, see ncm_data_fisher_matrix(). If *@IM is %NULL a new matrix is
 * allocated, otherwise *@IM must be a square matrix of the number of free parameters
 * and is overwritten. Without free parameters *@IM is freed and set to %NULL. Every
 * #NcmData is prepared before any is evaluated, see ncm_dataset_m2lnL_val().
 *
 */
void
ncm_dataset_fisher_matrix (NcmDataset *dset, NcmMSet *mset, NcmMatrix **IM)
{
  const guint fparams_len = ncm_mset_fparams_len (mset);

  if (fparams_len == 0)
  {
    ncm_matrix_clear (IM);
  }
  else
  {
    NcmMatrix *IM0 = NULL;
    guint i;

    _ncm_dataset_check_fisher_matrix (IM, fparams_len);
    ncm_matrix_set_zero (*IM);

    _ncm_dataset_prepare_all (dset, mset);

    for (i = 0; i < dset->oa->len; i++)
    {
      NcmData *data = ncm_dataset_peek_data (dset, i);

      ncm_data_fisher_matrix (data, mset, &IM0);
      ncm_matrix_add (*IM, IM0);
    }

    ncm_matrix_clear (&IM0);
  }
}

/**
 * ncm_dataset_fisher_matrix_bias:
 * @dset: a #NcmDataset
 * @mset: a #NcmMSet
 * @f_true: a #NcmVector
 * @IM: (inout) (allow-none) (transfer full): the Fisher matrix
 * @delta_theta: (inout) (allow-none) (transfer full): the parameter shift vector
 *
 * Calculates the Fisher-information matrix, as ncm_dataset_fisher_matrix(), and the
 * parameter shift @delta_theta obtained when the true mean is @f_true, of
 * ncm_dataset_get_n() components, adding those of each #NcmData in @dset, see
 * ncm_data_fisher_matrix_bias(). *@delta_theta is allocated when %NULL and
 * overwritten otherwise; without free parameters both are freed and set to %NULL.
 *
 */
void
ncm_dataset_fisher_matrix_bias (NcmDataset *dset, NcmMSet *mset, NcmVector *f_true, NcmMatrix **IM, NcmVector **delta_theta)
{
  const guint fparams_len = ncm_mset_fparams_len (mset);
  const guint total_n     = ncm_dataset_get_n (dset);

  if (ncm_vector_len (f_true) != total_n)
    g_error ("ncm_dataset_fisher_matrix_bias: the dataset has %u points, but f_true has %u.",
             total_n, ncm_vector_len (f_true));

  if (fparams_len == 0)
  {
    ncm_matrix_clear (IM);
    ncm_vector_clear (delta_theta);
  }
  else
  {
    NcmMatrix *IM0          = NULL;
    NcmVector *delta_theta0 = NULL;
    guint pos               = 0;
    guint i;

    _ncm_dataset_check_fisher_matrix (IM, fparams_len);

    if (*delta_theta == NULL)
      *delta_theta = ncm_vector_new (fparams_len);
    else if (ncm_vector_len (*delta_theta) != fparams_len)
      g_error ("ncm_dataset_fisher_matrix_bias: the bias vector has %u components, "
               "but there are %u free parameters.",
               ncm_vector_len (*delta_theta), fparams_len);

    ncm_matrix_set_zero (*IM);
    ncm_vector_set_zero (*delta_theta);

    _ncm_dataset_prepare_all (dset, mset);

    for (i = 0; i < dset->oa->len; i++)
    {
      NcmData *data = ncm_dataset_peek_data (dset, i);
      guint n       = ncm_data_get_length (data);

      ncm_vector_get_subvector2 (dset->ls_f, f_true, pos, n);

      ncm_data_fisher_matrix_bias (data, mset, dset->ls_f, &IM0, &delta_theta0);
      ncm_matrix_add (*IM, IM0);
      ncm_vector_add (*delta_theta, delta_theta0);
      pos += n;
    }

    ncm_matrix_clear (&IM0);
    ncm_vector_clear (&delta_theta0);
  }
}

/* Allocates *IM, or checks that the matrix passed in is fparams_len x fparams_len. */
static void
_ncm_dataset_check_fisher_matrix (NcmMatrix **IM, const guint fparams_len)
{
  if (*IM == NULL)
    *IM = ncm_matrix_new (fparams_len, fparams_len);
  else if ((ncm_matrix_nrows (*IM) != fparams_len) || (ncm_matrix_ncols (*IM) != fparams_len))
    g_error ("ncm_dataset_fisher_matrix: the Fisher matrix passed in is %u x %u, "
             "but there are %u free parameters.",
             ncm_matrix_nrows (*IM), ncm_matrix_ncols (*IM), fparams_len);
}


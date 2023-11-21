/***************************************************************************
 *            ncm_data.c
 *
 *  Sat Mar 29 15:51:46 2008
 *  Copyright  2008  Sandro Dias Pinto Vitenti
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
 * SECTION:ncm_data
 * @title: NcmData
 * @short_description: Abstract class for implementing data objects.
 *
 * The #NcmData object represent generic data. This is the root object used when
 * building a statistical analysis. Every implementation of #NcmData envolves
 * the methods described in #NcmDataClass.
 *
 * A #NcmData must implement, at least, the method ncm_data_m2lnL_val() or
 * ncm_data_leastsquares_f() to perform respectively likelihood or least
 * squares analysis.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_data.h"
#include "math/ncm_cfg.h"

enum
{
  PROP_0,
  PROP_NAME,
  PROP_DESC,
  PROP_LONG_DESC,
  PROP_INIT,
  PROP_BSTRAP,
  PROP_SIZE,
};

typedef struct _NcmDataPrivate
{
  gchar *desc;
  gchar *long_desc;
  gboolean init;
  gboolean begin;
  NcmBootstrap *bstrap;
  NcmDiff *diff;
} NcmDataPrivate;

G_DEFINE_ABSTRACT_TYPE_WITH_PRIVATE (NcmData, ncm_data, G_TYPE_OBJECT)

static void
ncm_data_init (NcmData *data)
{
  NcmDataPrivate * const self = ncm_data_get_instance_private (data);

  self->desc      = NULL;
  self->long_desc = NULL;
  self->init      = FALSE;
  self->begin     = FALSE;
  self->diff      = ncm_diff_new ();
}

static void
_ncm_data_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmData *data               = NCM_DATA (object);
  NcmDataPrivate * const self = ncm_data_get_instance_private (data);

  g_return_if_fail (NCM_IS_DATA (object));

  switch (prop_id)
  {
    case PROP_DESC:
      ncm_data_set_desc (data, g_value_get_string (value));
      break;
    case PROP_LONG_DESC:
      g_clear_pointer (&self->long_desc, g_free);
      self->long_desc = g_value_dup_string (value);
      break;
    case PROP_INIT:
      ncm_data_set_init (data, g_value_get_boolean (value));
      break;
    case PROP_BSTRAP:
    {
      NcmBootstrap *bstrap = g_value_get_object (value);

      ncm_data_bootstrap_set (data, bstrap);
      break;
    }
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_data_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmData *data               = NCM_DATA (object);
  NcmDataPrivate * const self = ncm_data_get_instance_private (data);
  NcmDataClass *data_class    = NCM_DATA_GET_CLASS (object);

  g_return_if_fail (NCM_IS_DATA (object));

  switch (prop_id)
  {
    case PROP_NAME:
      g_value_set_string (value, data_class->name);
      break;
    case PROP_DESC:
      g_value_set_string (value, ncm_data_peek_desc (data));
      break;
    case PROP_LONG_DESC:
      g_value_set_string (value, self->long_desc);
      break;
    case PROP_INIT:
      g_value_set_boolean (value, self->init);
      break;
    case PROP_BSTRAP:
    {
      g_value_set_object (value, self->bstrap);
      break;
    }
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_data_dispose (GObject *object)
{
  NcmData *data               = NCM_DATA (object);
  NcmDataPrivate * const self = ncm_data_get_instance_private (data);

  ncm_bootstrap_clear (&self->bstrap);
  ncm_diff_clear (&self->diff);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_data_parent_class)->dispose (object);
}

static void
_ncm_data_finalize (GObject *object)
{
  NcmData *data               = NCM_DATA (object);
  NcmDataPrivate * const self = ncm_data_get_instance_private (data);

  g_clear_pointer (&self->desc, g_free);
  g_clear_pointer (&self->long_desc, g_free);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_data_parent_class)->finalize (object);
}

static void _ncm_data_fisher_matrix (NcmData *data, NcmMSet *mset, NcmMatrix **IM);

static void
ncm_data_class_init (NcmDataClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);
  NcmDataClass *data_class   = NCM_DATA_CLASS (klass);

  object_class->set_property = &_ncm_data_set_property;
  object_class->get_property = &_ncm_data_get_property;
  object_class->dispose      = &_ncm_data_dispose;
  object_class->finalize     = &_ncm_data_finalize;

  /**
   * NcmData:name:
   *
   * Name of the data object.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_NAME,
                                   g_param_spec_string ("name",
                                                        NULL,
                                                        "Data type name",
                                                        NULL,
                                                        G_PARAM_READABLE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmData:desc:
   *
   * Description of the data object.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_DESC,
                                   g_param_spec_string ("desc",
                                                        NULL,
                                                        "Data description",
                                                        NULL,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmData:long-desc:
   *
   * Description of the data object.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_LONG_DESC,
                                   g_param_spec_string ("long-desc",
                                                        NULL,
                                                        "Data detailed description",
                                                        NULL,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmData:initialized:
   *
   * Whether the #NcmData is initialized.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_INIT,
                                   g_param_spec_boolean ("init",
                                                         NULL,
                                                         "Data initialized state",
                                                         FALSE,
                                                         G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmData:bootstrap:
   *
   * The #NcmData bootstrap object if any.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_BSTRAP,
                                   g_param_spec_object ("bootstrap",
                                                        NULL,
                                                        "Data bootstrap object",
                                                        NCM_TYPE_BOOTSTRAP,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  data_class->name           = NULL;
  data_class->get_length     = NULL;
  data_class->begin          = NULL;
  data_class->prepare        = NULL;
  data_class->resample       = NULL;
  data_class->leastsquares_f = NULL;
  data_class->m2lnL_val      = NULL;

  data_class->mean_vector = NULL;
  data_class->inv_cov_UH  = NULL;

  data_class->fisher_matrix = &_ncm_data_fisher_matrix;
}

typedef struct _NcmDataDiffArg
{
  NcmMSet *mset;
  NcmData *data;
} NcmDataDiffArg;

void
_ncm_data_diff_f (NcmVector *x, NcmVector *y, gpointer user_data)
{
  NcmDataDiffArg *arg = (NcmDataDiffArg *) user_data;

  ncm_mset_fparams_set_vector (arg->mset, x);
  ncm_data_mean_vector (arg->data, arg->mset, y);
}

static void
_ncm_data_fisher_matrix (NcmData *data, NcmMSet *mset, NcmMatrix **IM)
{
  NcmDataPrivate * const self = ncm_data_get_instance_private (data);
  const guint fparams_len     = ncm_mset_fparams_len (mset);
  NcmVector *x_v              = ncm_vector_new (fparams_len);
  NcmDataDiffArg arg          = {mset, data};
  const guint dim             = ncm_data_get_length (data);

  if (*IM == NULL)
  {
    *IM = ncm_matrix_new (fparams_len, fparams_len);
  }
  else
  {
    g_assert_cmpuint (ncm_matrix_ncols (*IM), ==, ncm_matrix_nrows (*IM));
    g_assert_cmpuint (ncm_matrix_ncols (*IM), ==, fparams_len);
  }

  ncm_mset_fparams_get_vector (mset, x_v);
  {
    GArray *x_a    = ncm_vector_dup_array (x_v);
    GArray *dmu_a  = ncm_diff_rf_d1_N_to_M (self->diff, x_a, dim, _ncm_data_diff_f, &arg, NULL);
    NcmMatrix *dmu = ncm_matrix_new_array (dmu_a, dim);

    ncm_data_inv_cov_UH (data, mset, dmu);

    ncm_matrix_dgemm (*IM, 'N', 'T', 1.0, dmu, dmu, 0.0);

    g_array_unref (dmu_a);
    g_array_unref (x_a);
    ncm_matrix_free (dmu);
  }
  ncm_mset_fparams_set_vector (mset, x_v);
  ncm_vector_free (x_v);
}

/**
 * ncm_data_ref:
 * @data: a #NcmData.
 *
 * Increase the reference count of @data.
 *
 * Returns: (transfer full): @data.
 */
NcmData *
ncm_data_ref (NcmData *data)
{
  return g_object_ref (data);
}

/**
 * ncm_data_free:
 * @data: a #NcmData.
 *
 * Decrease the reference count of @data.
 *
 */
void
ncm_data_free (NcmData *data)
{
  g_object_unref (data);
}

/**
 * ncm_data_clear:
 * @data: a #NcmData.
 *
 * Decrease the reference count of *@data and sets the pointer *@data to NULL.
 *
 */
void
ncm_data_clear (NcmData **data)
{
  g_clear_object (data);
}

/**
 * ncm_data_dup:
 * @data: a #NcmData.
 * @ser_obj: a #NcmSerialize.
 *
 * Duplicate the @data object.
 *
 * Returns: (transfer full): a duplicate of @data.
 */
NcmData *
ncm_data_dup (NcmData *data, NcmSerialize *ser_obj)
{
  return NCM_DATA (ncm_serialize_dup_obj (ser_obj, G_OBJECT (data)));
}

/**
 * ncm_data_new_from_file:
 * @filename: file containing a serialized #NcmData child.
 *
 * Creates a new #NcmData from @filename.
 *
 * Returns: (transfer full): the newly created #NcmData.
 */
NcmData *
ncm_data_new_from_file (const gchar *filename)
{
  NcmData *data = NCM_DATA (ncm_serialize_global_from_file (filename));

  g_assert (NCM_IS_DATA (data));

  return data;
}

/**
 * ncm_data_get_length: (virtual get_length)
 * @data: a #NcmData.
 *
 * Return a integer representing the number of data points.
 *
 * Returns: number of data points.
 */
guint
ncm_data_get_length (NcmData *data)
{
  g_assert (NCM_DATA_GET_CLASS (data)->get_length != NULL);

  return NCM_DATA_GET_CLASS (data)->get_length (data);
}

/**
 * ncm_data_get_dof: (virtual get_dof)
 * @data: a #NcmData.
 *
 * Calculates the degrees of freedom associated with the data.
 *
 * Returns: degrees of freedom of the data.
 */
guint
ncm_data_get_dof (NcmData *data)
{
  g_assert ((NCM_DATA_GET_CLASS (data)->get_dof != NULL) ||
            (NCM_DATA_GET_CLASS (data)->get_length != NULL));

  if (NCM_DATA_GET_CLASS (data)->get_dof != NULL)
    return NCM_DATA_GET_CLASS (data)->get_dof (data);
  else
    return NCM_DATA_GET_CLASS (data)->get_length (data);
}

/**
 * ncm_data_set_init:
 * @data: a #NcmData
 * @state: a boolean
 *
 * Sets the @data to initialized or not @state.
 *
 */
void
ncm_data_set_init (NcmData *data, gboolean state)
{
  NcmDataPrivate * const self = ncm_data_get_instance_private (data);

  if (self->init)
  {
    if (!state)
    {
      self->init  = FALSE;
      self->begin = FALSE;
    }
  }
  else
  {
    if (state)
    {
      self->init  = TRUE;
      self->begin = FALSE;
    }
  }
}

/**
 * ncm_data_is_init:
 * @data: a #NcmData
 *
 * Returns: whether the @data object is initialized.
 */
gboolean
ncm_data_is_init (NcmData *data)
{
  NcmDataPrivate * const self = ncm_data_get_instance_private (data);

  return self->init;
}

/**
 * ncm_data_set_desc:
 * @data: a #NcmData.
 * @desc: description.
 *
 * Sets the @data description. It gets a copy of desc.
 *
 */
void
ncm_data_set_desc (NcmData *data, const gchar *desc)
{
  NcmDataPrivate * const self = ncm_data_get_instance_private (data);

  g_clear_pointer (&self->desc, g_free);
  self->desc = g_strdup (desc);
}

/**
 * ncm_data_take_desc:
 * @data: a #NcmData.
 * @desc: description.
 *
 * Sets the @data description @desc without copying it, the @desc memory will
 * be freed (g_free()) when the object is freed.
 *
 */
void
ncm_data_take_desc (NcmData *data, gchar *desc)
{
  NcmDataPrivate * const self = ncm_data_get_instance_private (data);

  g_clear_pointer (&self->desc, g_free);
  self->desc = desc;
}

/**
 * ncm_data_peek_desc:
 * @data: a #NcmData.
 *
 * Gets @data description.
 *
 * Returns: (transfer none): internal @data description.
 */
const gchar *
ncm_data_peek_desc (NcmData *data)
{
  NcmDataPrivate * const self = ncm_data_get_instance_private (data);

  if (self->desc == NULL)
  {
    NcmDataClass *data_class = NCM_DATA_GET_CLASS (data);

    if (data_class->name == NULL)
      self->desc = g_strdup (G_OBJECT_TYPE_NAME (data));
    else
      self->desc = g_strdup (data_class->name);
  }

  return self->desc;
}

/**
 * ncm_data_get_desc:
 * @data: a #NcmData.
 *
 * Gets @data description.
 *
 * Returns: (transfer full): copy of the @data description.
 */
gchar *
ncm_data_get_desc (NcmData *data)
{
  return g_strdup (ncm_data_peek_desc (data));
}

static void
_ncm_data_prepare (NcmData *data, NcmMSet *mset)
{
  NcmDataPrivate * const self = ncm_data_get_instance_private (data);

  if ((NCM_DATA_GET_CLASS (data)->begin != NULL) && !self->begin)
  {
    NCM_DATA_GET_CLASS (data)->begin (data);
    self->begin = TRUE;
  }

  if (NCM_DATA_GET_CLASS (data)->prepare != NULL)
    NCM_DATA_GET_CLASS (data)->prepare (data, mset);
}

/**
 * ncm_data_prepare: (virtual prepare)
 * @data: a #NcmData.
 * @mset: a #NcmMSet.
 *
 * Prepare all models in @data necessary for the statistical calculations.
 *
 */
void
ncm_data_prepare (NcmData *data, NcmMSet *mset)
{
  NcmDataPrivate * const self = ncm_data_get_instance_private (data);

  g_assert (self->init);
  _ncm_data_prepare (data, mset);
}

/**
 * ncm_data_resample: (virtual resample)
 * @data: a #NcmData
 * @mset: a #NcmMSet
 * @rng: a #NcmRNG
 *
 * Resample data in @data from the models contained in @mset.
 *
 */
void
ncm_data_resample (NcmData *data, NcmMSet *mset, NcmRNG *rng)
{
  NcmDataPrivate * const self = ncm_data_get_instance_private (data);

  if (NCM_DATA_GET_CLASS (data)->resample == NULL)
    g_error ("ncm_data_resample: The data (%s) does not implement resample.",
             ncm_data_get_desc (data));

  self->begin = TRUE;
  _ncm_data_prepare (data, mset);

  NCM_DATA_GET_CLASS (data)->resample (data, mset, rng);
  self->begin = FALSE;

  if ((NCM_DATA_GET_CLASS (data)->begin != NULL) && !self->begin)
  {
    NCM_DATA_GET_CLASS (data)->begin (data);
    self->begin = TRUE;
  }

  ncm_data_set_init (data, TRUE);
}

/**
 * ncm_data_bootstrap_create:
 * @data: a #NcmData.
 *
 * Creates a bootstrap object inside of @data. Uses the default bsize == fsize.
 *
 */
void
ncm_data_bootstrap_create (NcmData *data)
{
  NcmDataPrivate * const self = ncm_data_get_instance_private (data);

  if (!NCM_DATA_GET_CLASS (data)->bootstrap)
    g_error ("ncm_data_bootstrap_create: The data (%s) does not implement bootstrap.",
             ncm_data_get_desc (data));

  g_assert (self->init);

  if (self->bstrap == NULL)
  {
    self->bstrap = ncm_bootstrap_sized_new (ncm_data_get_length (data));
  }
  else
  {
    ncm_bootstrap_set_fsize (self->bstrap, ncm_data_get_length (data));
    ncm_bootstrap_set_bsize (self->bstrap, ncm_data_get_length (data));
  }
}

/**
 * ncm_data_bootstrap_remove:
 * @data: a #NcmData.
 *
 * Removes a bootstrap object inside of @data if any.
 *
 */
void
ncm_data_bootstrap_remove (NcmData *data)
{
  NcmDataPrivate * const self = ncm_data_get_instance_private (data);

  ncm_bootstrap_clear (&self->bstrap);
}

/**
 * ncm_data_bootstrap_set:
 * @data: a #NcmData.
 * @bstrap: a #NcmBootstrap.
 *
 * Sets the @bstrap object in @data checking if they are compatible.
 *
 */
void
ncm_data_bootstrap_set (NcmData *data, NcmBootstrap *bstrap)
{
  NcmDataPrivate * const self = ncm_data_get_instance_private (data);

  if (!NCM_DATA_GET_CLASS (data)->bootstrap)
    g_error ("ncm_data_bootstrap_set: The data (%s) does not implement bootstrap.",
             ncm_data_get_desc (data));

  g_assert (self->init);
  g_assert (bstrap != NULL);

  g_assert_cmpuint (ncm_bootstrap_get_fsize (bstrap), ==, ncm_data_get_length (data));
  ncm_bootstrap_ref (bstrap);
  ncm_bootstrap_clear (&self->bstrap);
  self->bstrap = bstrap;
}

/**
 * ncm_data_bootstrap_resample:
 * @data: a #NcmData.
 * @rng: a #NcmRNG.
 *
 * Perform one bootstrap, i.e., resample the data with replacement.
 *
 */
void
ncm_data_bootstrap_resample (NcmData *data, NcmRNG *rng)
{
  NcmDataPrivate * const self = ncm_data_get_instance_private (data);

  if (!NCM_DATA_GET_CLASS (data)->bootstrap)
    g_error ("ncm_data_bootstrap_resample: The data (%s) does not implement bootstrap.",
             ncm_data_get_desc (data));

  if (self->bstrap == NULL)
    g_error ("ncm_data_bootstrap_resample: Bootstrap of %s is not enabled.",
             ncm_data_get_desc (data));

  ncm_bootstrap_resample (self->bstrap, rng);
}

/**
 * ncm_data_bootstrap_enabled:
 * @data: a #NcmData.
 *
 * Checks whether bootstrap is enabled in @data.
 *
 * Returns: if bootstrap is enabled in @data.
 */
gboolean
ncm_data_bootstrap_enabled (NcmData *data)
{
  NcmDataPrivate * const self = ncm_data_get_instance_private (data);

  if (NCM_DATA_GET_CLASS (data)->bootstrap && (self->bstrap != NULL))
    return TRUE;
  else
    return FALSE;
}

/**
 * ncm_data_peek_bootstrap:
 * @data: a #NcmData.
 *
 * Returns: (transfer none): the current #NcmBootstrap object or NULL.
 */
NcmBootstrap *
ncm_data_peek_bootstrap (NcmData *data)
{
  NcmDataPrivate * const self = ncm_data_get_instance_private (data);

  return self->bstrap;
}

/**
 * ncm_data_leastsquares_f: (virtual leastsquares_f)
 * @data: a #NcmData.
 * @mset: a #NcmMSet.
 * @f: a #NcmVector
 *
 * Calculates the least squares vector $\vec{f}$ using the models contained in
 * @mset and set the results in @f.
 *
 */
void
ncm_data_leastsquares_f (NcmData *data, NcmMSet *mset, NcmVector *f)
{
  ncm_data_prepare (data, mset);

  if (NCM_DATA_GET_CLASS (data)->leastsquares_f == NULL)
    g_error ("ncm_data_leastsquares_f: The data (%s) does not implement leastsquares_f.",
             ncm_data_get_desc (data));

  NCM_DATA_GET_CLASS (data)->leastsquares_f (data, mset, f);
}

/**
 * ncm_data_m2lnL_val: (virtual m2lnL_val)
 * @data: a #NcmData.
 * @mset: a #NcmMSet.
 * @m2lnL: (out): a #double
 *
 * Calculates the value of $-2\ln(L)$, where $L$ represents the likelihood of
 * the data given the models in @mset. The result is stored in @m2lnL.
 *
 */
void
ncm_data_m2lnL_val (NcmData *data, NcmMSet *mset, gdouble *m2lnL)
{
  ncm_data_prepare (data, mset);

  if (NCM_DATA_GET_CLASS (data)->m2lnL_val == NULL)
    g_error ("ncm_data_m2lnL_val: The data (%s) does not implement m2lnL_val.",
             ncm_data_get_desc (data));

  NCM_DATA_GET_CLASS (data)->m2lnL_val (data, mset, m2lnL);
}

/**
 * ncm_data_mean_vector: (virtual mean_vector)
 * @data: a #NcmData
 * @mset: a #NcmMSet
 * @mu: the mean output #NcmVector
 *
 * Calculates the Gaussian mean vector (for non-Gaussian distribution
 * it should calculate the Gaussian approximated mean of the actual
 * distribution).
 *
 */
void
ncm_data_mean_vector (NcmData *data, NcmMSet *mset, NcmVector *mu)
{
  ncm_data_prepare (data, mset);

  if (NCM_DATA_GET_CLASS (data)->mean_vector == NULL)
    g_error ("ncm_data_mean_vector: The data (%s) does not implement mean_vector.",
             ncm_data_get_desc (data));

  NCM_DATA_GET_CLASS (data)->mean_vector (data, mset, mu);
}

/**
 * ncm_data_inv_cov_UH: (virtual inv_cov_UH)
 * @data: a #NcmData
 * @mset: a #NcmMSet
 * @H: a #NcmMatrix
 *
 *
 * Given the Cholesky decomposition of the inverse covariance $C^{-1} = L\cdotU$
 * this function returns in-place the product $U\cdotH$.
 *
 */
void
ncm_data_inv_cov_UH (NcmData *data, NcmMSet *mset, NcmMatrix *H)
{
  ncm_data_prepare (data, mset);

  if (NCM_DATA_GET_CLASS (data)->inv_cov_UH == NULL)
    g_error ("ncm_data_inv_cov_UH: The data (%s) does not implement inv_cov_UH.",
             ncm_data_get_desc (data));

  NCM_DATA_GET_CLASS (data)->inv_cov_UH (data, mset, H);
}

/**
 * ncm_data_fisher_matrix: (virtual fisher_matrix)
 * @data: a #NcmData
 * @mset: a #NcmMSet
 * @IM: (out): The fisher matrix
 *
 * Calculates the Fisher-information matrix @I.
 *
 */
void
ncm_data_fisher_matrix (NcmData *data, NcmMSet *mset, NcmMatrix **IM)
{
  ncm_data_prepare (data, mset);

  if (NCM_DATA_GET_CLASS (data)->fisher_matrix == NULL)
    g_error ("ncm_data_fisher_matrix: The data (%s) does not implement fisher_matrix.",
             ncm_data_get_desc (data));

  NCM_DATA_GET_CLASS (data)->fisher_matrix (data, mset, IM);
}


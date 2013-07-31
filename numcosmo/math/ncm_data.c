/***************************************************************************
 *            ncm_data.c
 *
 *  Sat Mar 29 15:51:46 2008
 *  Copyright  2008  Sandro Dias Pinto Vitenti
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
 * SECTION:ncm_data
 * @title: Data Abstract Class
 * @short_description: Base class for implementing data objects
 *
 * The #NcmData object represent generic data. This is the root object used when
 * building a statistical analysis. Every implementation of #NcmData envolves 
 * the methods described in #NcmDataClass.
 * 
 * A #NcmData must implement, at least, the method #NcmDataClass.m2lnL_val or 
 * #NcmDataClass.leastsquares_f to perform respectively likelihood or least
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
  PROP_SIZE,
};

G_DEFINE_ABSTRACT_TYPE (NcmData, ncm_data, G_TYPE_OBJECT);

static void
ncm_data_init (NcmData *data)
{
  data->desc  = NULL;
  data->init  = FALSE;
  data->begin = FALSE;
}

static void
ncm_data_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  g_return_if_fail (NCM_IS_DATA (object));

  switch (prop_id)
  {
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
ncm_data_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmDataClass *data_class = NCM_DATA_GET_CLASS (object);

  g_return_if_fail (NCM_IS_DATA (object));

  switch (prop_id)
  {
    case PROP_NAME:
      g_value_set_string (value, data_class->name);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
ncm_data_finalize (GObject *object)
{
  NcmData *data = NCM_DATA (object);

  if (data->desc != NULL)
    g_free (data->desc);
  ncm_bootstrap_clear (&data->bstrap);
  
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_data_parent_class)->finalize (object);
}

static void _ncm_data_copyto (NcmData *data, NcmData *data_dest);

static void
ncm_data_class_init (NcmDataClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcmDataClass* data_class = NCM_DATA_CLASS (klass);

  object_class->finalize = ncm_data_finalize;
  object_class->set_property = ncm_data_set_property;
  object_class->get_property = ncm_data_get_property;

  /**
   * NcmData:name:
   *
   * Description of the data object.
   * 
   */  
  g_object_class_install_property (object_class,
                                   PROP_NAME,
                                   g_param_spec_string ("name",
                                                        NULL,
                                                        "Data type name",
                                                        NULL,
                                                        G_PARAM_READABLE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  
  data_class->name               = NULL;
  data_class->get_length         = NULL;
  data_class->copyto             = &_ncm_data_copyto;
  data_class->dup                = NULL;
  data_class->begin              = NULL;
  data_class->prepare            = NULL;
  data_class->resample           = NULL;
  data_class->leastsquares_f     = NULL;
  data_class->leastsquares_J     = NULL;
  data_class->leastsquares_f_J   = NULL;
  data_class->m2lnL_val          = NULL;
  data_class->m2lnL_grad         = NULL;
  data_class->m2lnL_val_grad     = NULL;
}

static void
_ncm_data_copyto (NcmData *data, NcmData *data_dest)
{
  if (data->desc != NULL)
    data_dest->desc         = g_strdup (data->desc);
  else
    data_dest->desc         = NULL;

  data_dest->init           = data->init;
  data_dest->begin          = data->begin;
  data_dest->bootstrap      = data->bootstrap;
  data_dest->bootstrap_init = data->bootstrap_init;

  if (data->bootstrap)
  {
    if (data_dest->bstrap == NULL)
      data_dest->bstrap = ncm_bootstrap_full_new (ncm_bootstrap_get_fsize (data->bstrap),
                                                  ncm_bootstrap_get_bsize (data->bstrap));
    else
    {
      ncm_bootstrap_set_fsize (data_dest->bstrap, ncm_bootstrap_get_fsize (data->bstrap));
      ncm_bootstrap_set_bsize (data_dest->bstrap, ncm_bootstrap_get_bsize (data->bstrap));
    }

    if (data->bootstrap_init)
    {
      NCM_GARRAY_MEMCPY (data_dest->bstrap->increasing_index, 
                         data->bstrap->increasing_index);
    }
  }
  else
    ncm_bootstrap_clear (&data_dest->bstrap);
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
 *
 * Duplicate the @data object.
 * 
 * Returns: (transfer full): a duplicate of @data.
 * 
 * Virtual: dup
 */
NcmData *
ncm_data_dup (NcmData *data)
{
  g_assert (NCM_DATA_GET_CLASS (data)->dup != NULL);
  return NCM_DATA_GET_CLASS (data)->dup (data);
}

/**
 * ncm_data_copyto:
 * @data: a #NcmData.
 * @data_dest: a #NcmData.
 * 
 * Copy the contents of @data to @dest_data. The objects must be compatibles.
 * 
 * Virtual: copyto
 */
void
ncm_data_copyto (NcmData *data, NcmData *data_dest)
{
  g_assert (NCM_DATA_GET_CLASS (data)->copyto != NULL);
  g_assert_cmpuint (G_OBJECT_TYPE (data), ==, G_OBJECT_TYPE (data_dest));

  if (&_ncm_data_copyto == NCM_DATA_GET_CLASS (data)->copyto)
    g_error ("ncm_data_copyto: object %s does not implement copyto.", G_OBJECT_TYPE_NAME (data));
  
  NCM_DATA_GET_CLASS (data)->copyto (data, data_dest);
}

/**
 * ncm_data_get_length:
 * @data: a #NcmData.
 *
 * Return a integer representing the number of data points.
 * 
 * Returns: number of data points.
 * 
 * Virtual: get_length
 */
guint
ncm_data_get_length (NcmData *data)
{
  g_assert (NCM_DATA_GET_CLASS (data)->get_length != NULL);
  return NCM_DATA_GET_CLASS (data)->get_length (data);
}

/**
 * ncm_data_get_dof:
 * @data: a #NcmData.
 *
 * Calculates the degrees of freedom associated with the data.
 * 
 * Returns: degrees of freedom of the data.
 * 
 * Virtual: get_dof
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
 * @data: a #NcmData.
 *
 * Sets the @data to an initiated state. 
 * 
 */
void
ncm_data_set_init (NcmData *data)
{
  data->init           = TRUE;
  data->begin          = FALSE;
  data->bootstrap      = FALSE;
  data->bootstrap_init = FALSE;
}

/**
 * ncm_data_prepare:
 * @data: a #NcmData.
 * @mset: a #NcmMSet.
 *
 * Prepare all models in @data necessary for the statistical calculations.
 * 
 * Virtual: prepare
 */
void
ncm_data_prepare (NcmData *data, NcmMSet *mset)
{
  g_assert (data->init);
  
  if (NCM_DATA_GET_CLASS (data)->begin != NULL && !data->begin)
  {
    NCM_DATA_GET_CLASS (data)->begin (data);
    data->begin = TRUE;
  }

  if (NCM_DATA_GET_CLASS (data)->prepare != NULL)
    NCM_DATA_GET_CLASS (data)->prepare (data, mset);
}

/**
 * ncm_data_resample:
 * @data: a #NcmData.
 * @mset: a #NcmMSet.
 *
 * Resample data in @data from the models contained in @mset.
 * 
 * Virtual: resample
 */
void
ncm_data_resample (NcmData *data, NcmMSet *mset)
{
  if (NCM_DATA_GET_CLASS (data)->resample == NULL)
    g_error ("ncm_data_resample: The data (%s) does not implement resample.\n", 
             NCM_DATA_GET_CLASS (data)->name);

  ncm_data_prepare (data, mset);
  
  NCM_DATA_GET_CLASS (data)->resample (data, mset);

  ncm_data_set_init (data);
}

/**
 * ncm_data_bootstrap_set:
 * @data: a #NcmData.
 * @enable: Whenever to enable or disable bootstrap.
 *
 * Enables or disable bootstrap analysis for the NcmData.
 * 
 */
void
ncm_data_bootstrap_set (NcmData *data, gboolean enable)
{
  if (!NCM_DATA_GET_CLASS (data)->bootstrap && enable)
    g_error ("ncm_data_bootstrap_set: The data (%s) does not implement bootstrap.\n", 
             NCM_DATA_GET_CLASS (data)->name);
  g_assert (data->init);
  data->bootstrap      = enable;
  data->bootstrap_init = FALSE;
  if (enable)
  {
    if (data->bstrap == NULL)
      data->bstrap = ncm_bootstrap_sized_new (ncm_data_get_length (data));
    else
      ncm_bootstrap_set_fsize (data->bstrap, ncm_data_get_length (data));
  }
  else
    ncm_bootstrap_clear (&data->bstrap);
}

/**
 * ncm_data_bootstrap_resample:
 * @data: a #NcmData.
 *
 * Perform one bootstrap, i.e., resample the data with replacement.
 * 
 */
void
ncm_data_bootstrap_resample (NcmData *data)
{
  if (!NCM_DATA_GET_CLASS (data)->bootstrap)
    g_error ("ncm_data_bootstrap_resample: The data (%s) does not implement bootstrap.\n", 
             NCM_DATA_GET_CLASS (data)->name);
  if (!data->bootstrap)
    g_error ("ncm_data_bootstrap_resample: Bootstrap of %s is not enabled.\n", 
             NCM_DATA_GET_CLASS (data)->name);

  data->bootstrap_init = TRUE;
  ncm_bootstrap_resample (data->bstrap);
}

/**
 * ncm_data_leastsquares_f:
 * @data: a #NcmData.
 * @mset: a #NcmMSet.
 * @f: a #NcmVector
 *
 * Calculates the least squares vector $\vec{f}$ using the models contained in
 * @mset and set the results in @f.
 * 
 * Virtual: leastsquares_f
 */
void 
ncm_data_leastsquares_f (NcmData *data, NcmMSet *mset, NcmVector *f)
{
  ncm_data_prepare (data, mset);

  if (NCM_DATA_GET_CLASS (data)->leastsquares_f == NULL)
    g_error ("ncm_data_leastsquares_f: The data (%s) does not implement leastsquares_f.\n", 
             NCM_DATA_GET_CLASS (data)->name);

  NCM_DATA_GET_CLASS (data)->leastsquares_f (data, mset, f);
}

/**
 * ncm_data_leastsquares_J:
 * @data: a #NcmData.
 * @mset: a #NcmMSet.
 * @J: a #NcmMatrix.
 *
 * Calculates the least squares jacobian matrix $$J_{ij} = \frac{df_i}{dx_j},$$
 * where $f_i$ is the component of the least squares vector $\vec{f}$ and $x_j$
 * is the j-th parameter.  
 * 
 * Virtual: leastsquares_J
 */
void 
ncm_data_leastsquares_J (NcmData *data, NcmMSet *mset, NcmMatrix *J)
{
  ncm_data_prepare (data, mset);

  if (NCM_DATA_GET_CLASS (data)->leastsquares_J == NULL)
    g_error ("ncm_data_leastsquares_J: The data (%s) does not implement leastsquares_J.\n", 
             NCM_DATA_GET_CLASS (data)->name);

  NCM_DATA_GET_CLASS (data)->leastsquares_J (data, mset, J);
}

/**
 * ncm_data_leastsquares_f_J:
 * @data: a #NcmData.
 * @mset: a #NcmMSet.
 * @f: a #NcmVector.
 * @J: a #NcmMatrix
 *
 * Calculates both least squares vector and matrix as in ncm_data_leastsquares_f()
 * and ncm_data_leastsquares_J().
 * 
 * Virtual: leastsquares_f_J
 */
void 
ncm_data_leastsquares_f_J (NcmData *data, NcmMSet *mset, NcmVector *f, NcmMatrix *J)
{
  ncm_data_prepare (data, mset);

  if (NCM_DATA_GET_CLASS (data)->leastsquares_f_J == NULL)
    g_error ("ncm_data_leastsquares_f_J: The data (%s) does not implement leastsquares_f_J.\n", 
             NCM_DATA_GET_CLASS (data)->name);

  NCM_DATA_GET_CLASS (data)->leastsquares_f_J (data, mset, f, J);
}

/**
 * ncm_data_m2lnL_val:
 * @data: a #NcmData.
 * @mset: a #NcmMSet.
 * @m2lnL: (out): a #double
 *
 * Calculates the value of $-2\ln(L)$, where $L$ represents the likelihood of
 * the data given the models in @mset. The result is stored in @m2lnL.
 * 
 * Virtual: m2lnL_val
 */
void 
ncm_data_m2lnL_val (NcmData *data, NcmMSet *mset, gdouble *m2lnL)
{
  ncm_data_prepare (data, mset);

  if (NCM_DATA_GET_CLASS (data)->m2lnL_val == NULL)
    g_error ("ncm_data_m2lnL_val: The data (%s) does not implement m2lnL_val.\n", 
             NCM_DATA_GET_CLASS (data)->name);

  NCM_DATA_GET_CLASS (data)->m2lnL_val (data, mset, m2lnL);
}

/**
 * ncm_data_m2lnL_grad:
 * @data: a #NcmData.
 * @mset: a #NcmMSet.
 * @grad: a #NcmVector.
 *
 * Calculates the gradient of $-2\ln(L)$, i.e., $$g_i = -2\frac{d\ln(L)}{dx_i}.$$ 
 * where $L$ represents the likelihood of the data given the models in @mset. 
 * The result is stored in @grad.
 * 
 * Virtual: m2lnL_grad
 */
void 
ncm_data_m2lnL_grad (NcmData *data, NcmMSet *mset, NcmVector *grad)
{
  ncm_data_prepare (data, mset);

  if (NCM_DATA_GET_CLASS (data)->m2lnL_grad == NULL)
    g_error ("ncm_data_m2lnL_grad: The data (%s) does not implement m2lnL_grad.\n", 
             NCM_DATA_GET_CLASS (data)->name);

  NCM_DATA_GET_CLASS (data)->m2lnL_grad (data, mset, grad);
}

/**
 * ncm_data_m2lnL_val_grad:
 * @data: a #NcmData.
 * @mset: a #NcmMSet.
 * @m2lnL: (out): a #double.
 * @grad: a #NcmVector.
 *
 * Calculates both the value and the gradient of $-2\ln(L)$ as in ncm_data_m2lnL_val() and
 * ncm_data_m2lnL_grad().
 * 
 * Virtual: m2lnL_val_grad
 */
void ncm_data_m2lnL_val_grad (NcmData *data, NcmMSet *mset, gdouble *m2lnL, NcmVector *grad)
{
  ncm_data_prepare (data, mset);

  if (NCM_DATA_GET_CLASS (data)->m2lnL_val_grad == NULL)
    g_error ("ncm_data_m2lnL_val_grad: The data (%s) does not implement m2lnL_val_grad.\n", 
             NCM_DATA_GET_CLASS (data)->name);

  NCM_DATA_GET_CLASS (data)->m2lnL_val_grad (data, mset, m2lnL, grad);
}

/**
 * ncm_data_get_desc:
 * @data: a #NcmData.
 *
 * Returns the description of the data in @data.
 * 
 * Returns: (transfer none): A short description of @data.
 */
gchar *
ncm_data_get_desc (NcmData *data)
{
  if (data->desc == NULL)
  {
    NcmDataClass *data_class = NCM_DATA_GET_CLASS (data);
    if (data_class->name == NULL)
      data->desc = g_strdup (G_OBJECT_TYPE_NAME (data));
    else
      data->desc = g_strdup (data_class->name);
  }
  return data->desc;
}


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
 * @title: Generic Data Object
 * @short_description: Object representing generic data
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_data.h"

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
  
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_data_parent_class)->finalize (object);
}

static void
ncm_data_class_init (NcmDataClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcmDataClass* data_class = NCM_DATA_CLASS (klass);

  object_class->finalize = ncm_data_finalize;
  object_class->set_property = ncm_data_set_property;
  object_class->get_property = ncm_data_get_property;

  g_object_class_install_property (object_class,
                                   PROP_NAME,
                                   g_param_spec_string ("name",
                                                        NULL,
                                                        "Data type name",
                                                        NULL,
                                                        G_PARAM_READABLE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  
  data_class->name             = NULL;
  data_class->get_length       = NULL;
  data_class->copyto           = NULL;
  data_class->dup              = NULL;
  data_class->begin            = NULL;
  data_class->prepare          = NULL;
  data_class->resample         = NULL;
  data_class->leastsquares_f   = NULL;
  data_class->leastsquares_J   = NULL;
  data_class->leastsquares_f_J = NULL;
  data_class->m2lnL_val        = NULL;
  data_class->m2lnL_grad       = NULL;
  data_class->m2lnL_val_grad   = NULL;
}

/**
 * ncm_data_ref:
 * @data: a #NcmData.
 *
 * FIXME
 * 
 * Returns: (transfer full): FIXME
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
 * FIXME
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
 * FIXME
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
 * FIXME
 * 
 * Returns: (transfer full): FIXME
 */
NcmData *
ncm_data_dup (NcmData *data)
{
  g_assert (NCM_DATA_GET_CLASS (data)->dup != NULL);
  return NCM_DATA_GET_CLASS (data)->dup (data);
}

/**
 * ncm_data_get_length:
 * @data: a #NcmData.
 *
 * FIXME
 * 
 * Returns: FIXME
 */
guint
ncm_data_get_length (NcmData *data)
{
  g_assert (NCM_DATA_GET_CLASS (data)->dup != NULL);
  return NCM_DATA_GET_CLASS (data)->get_length (data);
}

/**
 * ncm_data_set_init:
 * @data: a #NcmData.
 *
 * FIXME
 * 
 */
void
ncm_data_set_init (NcmData *data)
{
  data->init  = TRUE;
  data->begin = FALSE;
}

/**
 * ncm_data_prepare:
 * @data: a #NcmData.
 * @mset: a #NcmMSet.
 *
 * FIXME
 *
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
 * FIXME
 *
 */
void
ncm_data_resample (NcmData *data, NcmMSet *mset)
{
  if (NCM_DATA_GET_CLASS (data)->resample == NULL)
    g_error ("ncm_data_resample: The data (%s) does not implement resample\n", 
             NCM_DATA_GET_CLASS (data)->name);

  ncm_data_prepare (data, mset);
  
  NCM_DATA_GET_CLASS (data)->resample (data, mset);
}

/**
 * ncm_data_leastsquares_f:
 * @data: a #NcmData.
 * @mset: a #NcmMSet.
 * @f: FIXME
 *
 * FIXME
 *
 */
void 
ncm_data_leastsquares_f (NcmData *data, NcmMSet *mset, NcmVector *f)
{
  ncm_data_prepare (data, mset);

  if (NCM_DATA_GET_CLASS (data)->leastsquares_f == NULL)
    g_error ("ncm_data_leastsquares_f: The data (%s) does not implement leastsquares_f\n", 
             NCM_DATA_GET_CLASS (data)->name);

  NCM_DATA_GET_CLASS (data)->leastsquares_f (data, mset, f);
}

/**
 * ncm_data_leastsquares_J:
 * @data: a #NcmData.
 * @mset: a #NcmMSet.
 * @J: FIXME
 *
 * FIXME
 *
 */
void 
ncm_data_leastsquares_J (NcmData *data, NcmMSet *mset, NcmMatrix *J)
{
  ncm_data_prepare (data, mset);

  if (NCM_DATA_GET_CLASS (data)->leastsquares_J == NULL)
    g_error ("ncm_data_leastsquares_J: The data (%s) does not implement leastsquares_J\n", 
             NCM_DATA_GET_CLASS (data)->name);

  NCM_DATA_GET_CLASS (data)->leastsquares_J (data, mset, J);
}

/**
 * ncm_data_leastsquares_f_J:
 * @data: a #NcmData.
 * @mset: a #NcmMSet.
 * @f: FIXME
 * @J: FIXME
 *
 * FIXME
 *
 */
void 
ncm_data_leastsquares_f_J (NcmData *data, NcmMSet *mset, NcmVector *f, NcmMatrix *J)
{
  ncm_data_prepare (data, mset);

  if (NCM_DATA_GET_CLASS (data)->leastsquares_f_J == NULL)
    g_error ("ncm_data_leastsquares_f_J: The data (%s) does not implement leastsquares_f_J\n", 
             NCM_DATA_GET_CLASS (data)->name);

  NCM_DATA_GET_CLASS (data)->leastsquares_f_J (data, mset, f, J);
}

/**
 * ncm_data_m2lnL_val:
 * @data: a #NcmData.
 * @mset: a #NcmMSet.
 * @m2lnL: (out): FIXME
 *
 * FIXME
 *
 */
void 
ncm_data_m2lnL_val (NcmData *data, NcmMSet *mset, gdouble *m2lnL)
{
  ncm_data_prepare (data, mset);

  if (NCM_DATA_GET_CLASS (data)->m2lnL_val == NULL)
    g_error ("ncm_data_m2lnL_val: The data (%s) does not implement m2lnL_val\n", 
             NCM_DATA_GET_CLASS (data)->name);

  NCM_DATA_GET_CLASS (data)->m2lnL_val (data, mset, m2lnL);
}

/**
 * ncm_data_m2lnL_grad:
 * @data: a #NcmData.
 * @mset: a #NcmMSet.
 * @grad: FIXME
 *
 * FIXME
 *
 */
void 
ncm_data_m2lnL_grad (NcmData *data, NcmMSet *mset, NcmVector *grad)
{
  ncm_data_prepare (data, mset);

  if (NCM_DATA_GET_CLASS (data)->m2lnL_grad == NULL)
    g_error ("ncm_data_m2lnL_grad: The data (%s) does not implement m2lnL_grad\n", 
             NCM_DATA_GET_CLASS (data)->name);

  NCM_DATA_GET_CLASS (data)->m2lnL_grad (data, mset, grad);
}

/**
 * ncm_data_m2lnL_val_grad:
 * @data: a #NcmData.
 * @mset: a #NcmMSet.
 * @m2lnL: (out): FIXME
 * @grad: FIXME
 *
 * FIXME
 *
 */
void ncm_data_m2lnL_val_grad (NcmData *data, NcmMSet *mset, gdouble *m2lnL, NcmVector *grad)
{
  ncm_data_prepare (data, mset);

  if (NCM_DATA_GET_CLASS (data)->m2lnL_val_grad == NULL)
    g_error ("ncm_data_m2lnL_val_grad: The data (%s) does not implement m2lnL_val_grad\n", 
             NCM_DATA_GET_CLASS (data)->name);

  NCM_DATA_GET_CLASS (data)->m2lnL_val_grad (data, mset, m2lnL, grad);
}

/**
 * ncm_data_get_desc:
 * @data: a #NcmData.
 *
 * FIXME
 * 
 * Returns: (transfer none): FIXME
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


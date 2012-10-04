/***************************************************************************
 *            data.c
 *
 *  Sat Mar 29 15:51:46 2008
 *  Copyright  2008  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <sandro@lapsandro>
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
 * SECTION:data
 * @title: Generic Data Object
 * @short_description: Object representing generic data
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include <numcosmo/numcosmo.h>

#include <glib.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_sys.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_combination.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_eigen.h>
#include <math.h>
#include <string.h>
#include <gmp.h>

G_DEFINE_BOXED_TYPE (NcData, nc_data, nc_data_copy, nc_data_free);
G_DEFINE_BOXED_TYPE (NcDataStruct, nc_data_struct, nc_data_struct_copy, nc_data_struct_free);

/**
 * nc_data_struct_new:
 *
 * FIXME
 *
 * Returns: FIXME
 */
NcDataStruct *
nc_data_struct_new (void)
{
  NcDataStruct *dts = g_slice_new0 (NcDataStruct);
  return dts;
}

/**
 * nc_data_struct_copy:
 * @dts: A #NcDataStruct
 *
 * FIXME
 *
 * Returns: FIXME
 */
NcDataStruct *
nc_data_struct_copy (NcDataStruct *dts)
{
  NcDataStruct *clone = nc_data_struct_new ();
  memcpy (clone, dts, sizeof (NcDataStruct));
  if (dts->dup)
    clone->data = dts->dup (dts->data);
  else
    g_error ("NcDataStruct do not implement duplication.");

  if (NC_DATA_STRUCT_HAS_BEGIN (dts))
    NC_DATA_STRUCT_BEGIN (dts);

  return clone;
}

/**
 * nc_data_struct_free:
 * @dts: A #NcDataStruct
 *
 * FIXME
 *
 * Returns: FIXME
 */
void
nc_data_struct_free (NcDataStruct *dts)
{
  if (NC_DATA_STRUCT_HAS_FREE(dts))
    dts->free (dts->data);
  g_slice_free (NcDataStruct, dts);
}

/**
 * nc_data_new:
 *
 * FIXME
 *
 * Returns: FIXME
 */
NcData *
nc_data_new (void)
{
  NcData *data = g_slice_new0 (NcData);

  data->name                  = NULL;
  data->type                  = 0;
  data->init                  = FALSE;
  data->model                 = NULL;
  data->model_free            = NULL;
  data->model_ref             = NULL;
  data->model_init            = NULL;
  data->dts                   = NULL;
  data->resample              = NULL;
  data->prepare               = NULL;
  data->calc_leastsquares_f   = NULL;
  data->calc_leastsquares_J   = NULL;
  data->calc_leastsquares_f_J = NULL;
  data->calc_m2lnL_val        = NULL;
  data->calc_m2lnL_grad       = NULL;
  data->calc_m2lnL_val_grad   = NULL;

  return data;
}

/**
 * nc_data_free:
 * @data: pointer to type defined by #NcData
 *
 * FIXME
 */
void
nc_data_free (NcData *data)
{
  if (NC_DATA_STRUCT_HAS_END (data->dts))
    NC_DATA_STRUCT_END (data->dts);

  nc_data_struct_free (data->dts);

  if ((data->model_free != NULL) && (data->model != NULL))
	data->model_free (data->model);

  if (data->orig != NULL)
    nc_data_struct_free (data->orig);

  g_slice_free (NcData, data);
  return;
}

/**
 * nc_data_model_init:
 * @data: pointer to type defined by #NcData
 *
 * FIXME
 *
 * Returns: FIXME
 */
gboolean
nc_data_model_init (NcData *data)
{
  if (NC_DATA_HAS_MODEL_INIT (data))
    NC_DATA_MODEL_INIT (data);
  return TRUE;
}


/**
 * nc_data_begin:
 * @data: pointer to type defined by #NcData
 *
 * FIXME
 *
 * Returns: FIXME
 */
gboolean
nc_data_begin (NcData *data)
{
  if (NC_DATA_STRUCT_HAS_BEGIN (data->dts))
    NC_DATA_STRUCT_BEGIN (data->dts);
  return TRUE;
}

/**
 * nc_data_init:
 * @data: pointer to type defined by #NcData
 *
 * This function must be called when a NcData is
 * loaded in order to call all household methods
 * and its sets data->init = TRUE
 *
 * Returns: TRUE is the initialization went ok
 */
gboolean
nc_data_init (NcData *data)
{
  nc_data_begin (data);
  nc_data_model_init (data);
  data->init = TRUE;
  return TRUE;
}

/**
 * nc_data_copyto:
 * @dest: pointer to type defined by #NcData
 * @data: pointer to type defined by #NcData
 *
 * FIXME
 *
 * Returns: FIXME
 */
gboolean
nc_data_copyto (NcData *dest, NcData *src)
{
  if (NC_DATA_STRUCT_HAS_COPY (src->dts))
    NC_DATA_STRUCT_COPY (dest->dts, src->dts);
  else
    g_error ("Cannot copy NcDataStruct (%s), function not implemented.", src->name);

  if (NC_DATA_STRUCT_HAS_BEGIN (dest->dts))
    NC_DATA_STRUCT_BEGIN (dest->dts);

  return TRUE;
}

/**
 * nc_data_copy:
 * @data: pointer to type defined by #NcData
 *
 * FIXME
 *
 * Returns: FIXME
 */
NcData *
nc_data_copy (NcData *data)
{
  NcData *clone = nc_data_new ();

  clone->name                  = data->name;
  clone->type                  = data->type;
  clone->init                  = data->init;
  clone->model                 = data->model_ref (data->model);
  clone->model_ref             = data->model_ref;
  clone->model_free            = data->model_free;
  clone->model_init            = data->model_init;
  clone->resample              = data->resample;
  clone->prepare               = data->prepare;
  clone->calc_leastsquares_f   = data->calc_leastsquares_f;
  clone->calc_leastsquares_J   = data->calc_leastsquares_J;
  clone->calc_leastsquares_f_J = data->calc_leastsquares_f_J;
  clone->calc_m2lnL_val        = data->calc_m2lnL_val;
  clone->calc_m2lnL_grad       = data->calc_m2lnL_grad;
  clone->calc_m2lnL_val_grad   = data->calc_m2lnL_val_grad;

  if (NC_DATA_STRUCT_HAS_DUP (data->dts))
    clone->dts = nc_data_struct_copy (data->dts);
  else
    g_error ("Cannot duplicate NcData (%s), function not implemented.", data->name);
  clone->clone = TRUE;
  return clone;
}

/**
 * nc_data_resample:
 * @data: a #NcData
 * @mset: a #NcmMSet
 * @save: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gboolean
nc_data_resample (NcData *data, NcmMSet *mset, gboolean save)
{
  if (!data->resample)
    g_error ("The pdf (%s) do not implement resample\n", data->name);

  if (save && (data->orig == NULL))
    data->orig = nc_data_struct_copy (data->dts);

  if (data->prepare)
    NC_DATA_PREPARE (data, mset);

  data->resample (mset, NC_DATA_MODEL (data), NC_DATA_DATA (data));

  if (NC_DATA_STRUCT_HAS_BEGIN (data->dts))
    NC_DATA_STRUCT_BEGIN (data->dts);

  return TRUE;
}

/**
 * nc_data_set_orig:
 * @data: pointer to type defined by #NcData
 *
 * FIXME
 *
 * Returns: FIXME
 */
gboolean
nc_data_set_orig (NcData *data)
{
  if (data->orig == NULL)
    return FALSE;

  if (NC_DATA_STRUCT_HAS_COPY (data->dts))
  {
    NC_DATA_STRUCT_COPY (data->dts, data->orig);
    if (NC_DATA_STRUCT_HAS_FREE (data->orig))
      nc_data_struct_free (data->orig);
    data->orig = NULL;
  }
  else
    g_error ("Cannot set original data, copy not implemented.");

  if (NC_DATA_STRUCT_HAS_BEGIN (data->dts))
    NC_DATA_STRUCT_BEGIN (data->dts);

  return TRUE;
}

/**
 * nc_data_bootstrap:
 * @data: pointer to type defined by #NcData
 *
 * FIXME
 *
 * Returns: FIXME
 */
gboolean
nc_data_bootstrap (NcData *data)
{
  g_assert_not_reached ();
  return TRUE;
}

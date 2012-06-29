/***************************************************************************
 *            data_onevardist.c
 *
 *  Thu Apr 15 11:16:11 2010
 *  Copyright  2010  Mariana Penna Lima
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Mariana Penna Lima 2012 <pennalima@gmail.com>
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
 * SECTION:data_onevardist
 * @title: One Variable Distribution Data
 * @short_description: #NcData object representing a one variable distribution data
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include <numcosmo/numcosmo.h>

#include <glib.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_sf_gamma.h>
#include <gmp.h>

/********************************************************************************************************/

static NcDataOneVarDist *
_onevardist_data_alloc (guint np)
{
  NcDataOneVarDist *onevardist = g_slice_new0 (NcDataOneVarDist);

  onevardist->np         = np;
  onevardist->x          = ncm_vector_new (np);
  onevardist->extra_data = NULL;

  return onevardist;
}

/**********************************************************************************************/

static void
_onevardist_data_copy (gpointer dest_ptr, gpointer src_ptr)
{
  NcDataOneVarDist *dest = (NcDataOneVarDist *) dest_ptr;
  NcDataOneVarDist *src = (NcDataOneVarDist *) src_ptr;

  ncm_vector_memcpy (dest->x, src->x);

  if (src->extra_data && dest->extra_data)
  {
	if (NC_DATA_STRUCT_HAS_COPY(src->extra_data))
	  NC_DATA_STRUCT_COPY (dest->extra_data, src->extra_data);
	else
	  g_error ("NcDataStruct extra do not implement copy.");
  }

}

static void
_onevardist_data_free (gpointer onevardist_ptr)
{
  NcDataOneVarDist *onevardist = (NcDataOneVarDist *) onevardist_ptr;

  ncm_vector_free (onevardist->x);

  if (onevardist->extra_data)
	nc_data_struct_free (onevardist->extra_data);

  g_slice_free (NcDataOneVarDist, onevardist);
}

static guint _onevardist_get_length (gpointer onevardist_ptr) { return ((NcDataOneVarDist *) onevardist_ptr)->np; }

/********************************************************************************************************/

static void
_onevardist_calc_m2lnL (NcmMSet *mset, gpointer model, gpointer data, gdouble *m2lnL)
{
  NcDataModelOneVarDist *dmovd = (NcDataModelOneVarDist *) model;
  NcDataOneVarDist *onevardist = (NcDataOneVarDist *) data;

  gint i;
  *m2lnL = 0.0;

  for (i = 0; i < onevardist->np; i++)
  {
	const gdouble x_i = ncm_vector_get (onevardist->x, i);
	*m2lnL += ncm_mset_func_eval1 (dmovd->dist, mset, x_i);
  }

  return;
}

/********************************************************************************************************/

static gpointer
_onevardist_data_dup (gpointer onevardist_ptr)
{
  NcDataOneVarDist *onevardist = (NcDataOneVarDist *) onevardist_ptr;
  NcDataOneVarDist *clone;
  clone = _onevardist_data_alloc (onevardist->np);
  _onevardist_data_copy (clone, onevardist);

  return clone;
}

/********************************************************************************************************/

static void
_onevardist_data_end (gpointer onevardist_ptr)
{
  NcDataOneVarDist *onevardist = (NcDataOneVarDist *) onevardist_ptr;

  if (onevardist->extra_data)
  {
	if (NC_DATA_STRUCT_HAS_END (onevardist->extra_data))
	  NC_DATA_STRUCT_END (onevardist->extra_data);
  }
}

/********************************************************************************************************/

static void
_onevardist_data_begin (gpointer onevardist_ptr)
{
  NcDataOneVarDist *onevardist = (NcDataOneVarDist *) onevardist_ptr;

  if (onevardist->extra_data)
  {
	if (NC_DATA_STRUCT_HAS_BEGIN (onevardist->extra_data))
	  NC_DATA_STRUCT_BEGIN (onevardist->extra_data);
  }
}

/********************************************************************************************************/

static void
_onevardist_resample (NcmMSet *mset, gpointer model, gpointer data)
{
  NcDataModelOneVarDist *dmovd = (NcDataModelOneVarDist *) model;
  NcDataOneVarDist *onevardist = (NcDataOneVarDist *) data;
  gsl_rng *rng = ncm_get_rng ();
  gint i;

  if (dmovd->inv_pdf == NULL)
	g_error ("This data do not implement the inverse of the pdf.");

  for (i = 0; i < onevardist->np; i++)
  {
	const gdouble u_i = gsl_rng_uniform (rng);
	const gdouble x_i = ncm_mset_func_eval1 (dmovd->inv_pdf, mset, u_i);
	ncm_vector_set (onevardist->x, i, x_i);
  }
}

/********************************************************************************************************/

static gpointer
_onevardist_model_ref (gpointer model)
{
  NcDataModelOneVarDist *dmovd = (NcDataModelOneVarDist *) model;
  NcDataModelOneVarDist *dmovd_ref = g_slice_new (NcDataModelOneVarDist);

  dmovd_ref->dist    = ncm_mset_func_ref (dmovd->dist);
  dmovd_ref->inv_pdf = ncm_mset_func_ref (dmovd->inv_pdf);

  return dmovd_ref;
}

/********************************************************************************************************/

static void
_onevardist_model_free (gpointer model)
{
  NcDataModelOneVarDist *dmovd = (NcDataModelOneVarDist *) model;
  ncm_mset_func_free (dmovd->dist);
  ncm_mset_func_free (dmovd->inv_pdf);

  g_slice_free (NcDataModelOneVarDist, dmovd);
}

/********************************************************************************************************/

static NcDataStruct *
_nc_data_struct_onevardist_new (void)
{
  NcDataStruct *dts = nc_data_struct_new ();
  dts->data       = NULL;
  dts->dup        = &_onevardist_data_dup;
  dts->copy       = &_onevardist_data_copy;
  dts->free       = &_onevardist_data_free;
  dts->begin      = &_onevardist_data_begin;
  dts->end        = &_onevardist_data_end;
  dts->get_length = &_onevardist_get_length;

  return dts;
}

/********************************************************************************************************/

/**
 * nc_data_onevardist_new:
 * @dist: FIXME
 * @inv_pdf: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
NcData *
nc_data_onevardist_new (NcmMSetFunc *dist, NcmMSetFunc *inv_pdf)
{
  NcData *data = nc_data_new ();
  NcDataModelOneVarDist *dmovd = g_slice_new (NcDataModelOneVarDist);
  static gchar *name = "One dimension distribution";

  g_assert (dist != NULL);

  dmovd->dist = dist;
  dmovd->inv_pdf = inv_pdf;

  data->name                  = name;
  data->type                  = 0;
  data->init                  = FALSE;
  data->model                 = dmovd;
  data->model_init            = NULL;
  data->model_ref             = &_onevardist_model_ref;
  data->model_free            = &_onevardist_model_free;
  data->dts                   = _nc_data_struct_onevardist_new ();
  data->resample              = &_onevardist_resample;
  data->prepare               = NULL;
  data->calc_leastsquares_f   = NULL;
  data->calc_leastsquares_J   = NULL;
  data->calc_leastsquares_f_J = NULL;
  data->calc_m2lnL_val        = &_onevardist_calc_m2lnL;
  data->calc_m2lnL_grad       = NULL;
  data->calc_m2lnL_val_grad   = NULL;

  return data;
}

/********************************************************************************************************/

/**
 * nc_data_onevardist_init_from_vector: (skip)
 * @data: a #NcData
 * @x: a #NcmVector
 * @extra_data: a #NcDataStruct
 *
 * FIXME
 *
 * Returns: FIXME
 */
void
nc_data_onevardist_init_from_vector (NcData *data, NcmVector *x, NcDataStruct *extra_data)
{
  NcDataOneVarDist *onevardist = g_slice_new0 (NcDataOneVarDist);
  onevardist->np = ncm_vector_len (x);
  onevardist->x = x;
  g_object_ref (x);

  NC_DATA_STRUCT_DATA (data->dts) = onevardist;
  onevardist->extra_data = extra_data;

  nc_data_init (data);
}

/**
 * nc_data_onevardist_set_prepare:
 * @data: FIXME
 * @prepare: (scope call): FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
void
nc_data_onevardist_set_prepare (NcData *data, NcDataPrepare prepare)
{
  data->prepare = prepare;
}

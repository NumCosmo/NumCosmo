/***************************************************************************
 *            data_gaussian.c
 *
 *  Fri Mar 19 14:57:35 2010
 *  Copyright  2010  Sandro Dias Pinto Vitenti
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
 * SECTION:data_gaussian
 * @title: Gaussian Data
 * @short_description: Generic Gaussian object
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "data/data_gaussian.h"
#include "math/util.h"

#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_blas.h>
#ifdef NUMCOSMO_HAVE_SQLITE3
#include <sqlite3.h>
#endif

/********************************************************************************************************/

static NcDataGaussian *
_gauss_data_alloc (guint np, gboolean has_x)
{
  NcDataGaussian *gauss = g_slice_new0 (NcDataGaussian);

  gauss->wt = 0.0;
  gauss->np = np;
  if (has_x)
	gauss->x = ncm_vector_new (np);
  gauss->y = ncm_vector_new (np);
  gauss->sigma = ncm_vector_new (np);
  gauss->weight = ncm_vector_new (np);
  gauss->extra_data = NULL;

  return gauss;
}

static NcDataGaussian *
_gauss_cov_data_alloc (guint np, gboolean has_x)
{
  NcDataGaussian *gauss = g_slice_new0 (NcDataGaussian);

  gauss->wt = 0.0;
  gauss->np = np;
  if (has_x)
	gauss->x = ncm_vector_new (np);
  gauss->y = ncm_vector_new (np);
  gauss->sigma = ncm_vector_new (np);
  gauss->weight = ncm_vector_new (np);
  gauss->inv_cov = ncm_matrix_new (np, np);
  gauss->orto = ncm_matrix_new (np, np);
  gauss->vp = ncm_vector_new (np);
  gauss->vw = gsl_eigen_symmv_alloc (np);
  gauss->extra_data = NULL;

  return gauss;
}

/**********************************************************************************************/

static void
_gauss_cov_data_copy (gpointer dest_ptr, gpointer src_ptr)
{
  NcDataGaussian *dest = (NcDataGaussian *) dest_ptr;
  NcDataGaussian *src = (NcDataGaussian *) src_ptr;

  dest->wt = src->wt;
  if (src->x)
	ncm_vector_memcpy (dest->x, src->x);
  ncm_vector_memcpy (dest->y, src->y);
  ncm_vector_memcpy (dest->weight, src->weight);
  ncm_vector_memcpy (dest->sigma, src->sigma);
  ncm_matrix_memcpy (dest->inv_cov, src->inv_cov);
  ncm_matrix_memcpy (dest->orto, src->orto);

  if (src->extra_data && dest->extra_data)
  {
	if (NC_DATA_STRUCT_HAS_COPY(src->extra_data))
	  NC_DATA_STRUCT_COPY (dest->extra_data, src->extra_data);
	else
	  g_error ("NcDataStruct extra do not implement copy.");
  }
}

static void
_gauss_data_copy (gpointer dest_ptr, gpointer src_ptr)
{
  NcDataGaussian *dest = (NcDataGaussian *) dest_ptr;
  NcDataGaussian *src = (NcDataGaussian *) src_ptr;

  dest->wt = src->wt;
  if (src->x)
	ncm_vector_memcpy (dest->x, src->x);
  ncm_vector_memcpy (dest->y, src->y);
  ncm_vector_memcpy (dest->sigma, src->sigma);
  ncm_vector_memcpy (dest->weight, src->weight);

  if (src->extra_data && dest->extra_data)
  {
	if (NC_DATA_STRUCT_HAS_COPY(src->extra_data))
	  NC_DATA_STRUCT_COPY (dest->extra_data, src->extra_data);
	else
	  g_error ("NcDataStruct extra do not implement copy.");
  }
}

static void
_gauss_cov_data_free (gpointer gauss_ptr)
{
  NcDataGaussian *gauss = (NcDataGaussian *) gauss_ptr;

  if (gauss->x)
	ncm_vector_free (gauss->x);
  ncm_vector_free (gauss->y);
  ncm_vector_free (gauss->sigma);
  ncm_vector_free (gauss->weight);
  ncm_vector_free (gauss->vp);
  ncm_matrix_free (gauss->inv_cov);
  ncm_matrix_free (gauss->orto);
  gsl_eigen_symmv_free (gauss->vw);

  if (gauss->extra_data)
	nc_data_struct_free (gauss->extra_data);

  g_slice_free (NcDataGaussian, gauss);
}

static void _gauss_data_free (gpointer gauss_ptr)
{
  NcDataGaussian *gauss = (NcDataGaussian *)gauss_ptr;

  if (gauss->x)
	ncm_vector_free (gauss->x);
  ncm_vector_free (gauss->y);
  ncm_vector_free (gauss->sigma);
  ncm_vector_free (gauss->weight);

  if (gauss->extra_data)
	nc_data_struct_free (gauss->extra_data);

  g_slice_free (NcDataGaussian, gauss);
}

static void
_gauss_cov_data_begin (gpointer gauss_ptr)
{
  NcDataGaussian *gauss = (NcDataGaussian *) gauss_ptr;
  NcmMatrix *inv_cov = ncm_matrix_new (NCM_MATRIX_NROWS (gauss->inv_cov), NCM_MATRIX_NCOLS (gauss->inv_cov));
  gint i;

  ncm_matrix_memcpy (inv_cov, gauss->inv_cov);
  gsl_eigen_symmv (NCM_MATRIX_GSL (inv_cov), ncm_vector_gsl(gauss->weight), NCM_MATRIX_GSL (gauss->orto), gauss->vw);
  ncm_matrix_transpose (gauss->orto);

  gauss->wt = 0.0;
  for (i = 0; i < gauss->np; i++)
  {
	const gdouble w = ncm_vector_get (gauss->weight, i);
	const gdouble sigma_i = 1.0 / sqrt(w);
	ncm_vector_set (gauss->sigma, i, sigma_i);
	gauss->wt += w;
  }

  if (gauss->extra_data)
  {
	if (NC_DATA_STRUCT_HAS_BEGIN (gauss->extra_data))
	  NC_DATA_STRUCT_BEGIN (gauss->extra_data);
  }

  ncm_matrix_free (inv_cov);
}

static void
_gauss_data_begin (gpointer gauss_ptr)
{
  NcDataGaussian *gauss = (NcDataGaussian *) gauss_ptr;
  gint i;

  gauss->wt = 0.0;
  for (i = 0; i < gauss->np; i++)
  {
	const gdouble sigma_i = ncm_vector_get (gauss->sigma, i);
	const gdouble w_i = 1.0 / (sigma_i * sigma_i);
	ncm_vector_set (gauss->weight, i, w_i);
	gauss->wt += w_i;
  }

  if (gauss->extra_data)
  {
	if (NC_DATA_STRUCT_HAS_BEGIN (gauss->extra_data))
	  NC_DATA_STRUCT_BEGIN (gauss->extra_data);
  }
}

/**********************************************************************************************/

static void
_gauss_data_end (gpointer gauss_ptr)
{
  NcDataGaussian *gauss = (NcDataGaussian *) gauss_ptr;

  if (gauss->extra_data)
  {
	if (NC_DATA_STRUCT_HAS_END (gauss->extra_data))
	  NC_DATA_STRUCT_END (gauss->extra_data);
  }
}

/**********************************************************************************************/

static gpointer
_gauss_data_dup (gpointer gauss_ptr)
{
  NcDataGaussian *gauss = (NcDataGaussian *) gauss_ptr;
  NcDataGaussian *clone;
  gdouble has_x = (gauss->x != NULL);
  clone = _gauss_data_alloc (gauss->np, has_x);
  _gauss_data_copy (clone, gauss);

  if (gauss->extra_data)
  {
	if (NC_DATA_STRUCT_HAS_DUP (gauss->extra_data))
	  clone->extra_data = nc_data_struct_copy (gauss->extra_data);
	else
	  g_error ("NcDataStruct extra do not implement duplication.");
  }

  return clone;
}

static gpointer
_gauss_cov_data_dup (gpointer gauss_ptr)
{
  NcDataGaussian *gauss = (NcDataGaussian *) gauss_ptr;
  NcDataGaussian *clone;
  gdouble has_x = (gauss->x != NULL);
  clone = _gauss_cov_data_alloc (gauss->np, has_x);
  _gauss_cov_data_copy (clone, gauss);

  if (gauss->extra_data)
  {
	if (NC_DATA_STRUCT_HAS_DUP (gauss->extra_data))
	  clone->extra_data = nc_data_struct_copy (gauss->extra_data);
	else
	  g_error ("NcDataStruct extra do not implement duplication.");
  }

  return clone;
}

/**********************************************************************************************/

static guint _gauss_get_length (gpointer gauss_ptr) { return ((NcDataGaussian *) gauss_ptr)->np; }

/********************************************************************************************************/

void
_gaussian_x_resample (NcmMSet *mset, gpointer model, gpointer data)
{
  GPtrArray *fa = (GPtrArray *) model;
  NcDataGaussian *gauss = (NcDataGaussian *) data;
  gint i;
  gsl_rng *rng = ncm_get_rng ();

  for (i = 0; i < gauss->np; i++)
  {
	NcmMSetFunc *mean_func = g_ptr_array_index (fa, i % fa->len);
	const gdouble x_i = ncm_vector_get (gauss->x, i);
	const gdouble yt_i = ncm_mset_func_eval1 (mean_func, mset, x_i);
	const gdouble sigma_i = ncm_vector_get (gauss->sigma, i);
	const gdouble y_i = yt_i + gsl_ran_gaussian (rng, sigma_i);
	ncm_vector_set (gauss->y, i, y_i);
  }
}

/********************************************************************************************************/

static void
_gauss_wmean_x_calc_leastsquares_f (NcmMSet *mset, gpointer model, gpointer data, NcmVector *v)
{
  GPtrArray *fa = (GPtrArray *) model;
  NcDataGaussian *gauss = (NcDataGaussian *) data;
  gdouble wmean;
  gint i;

  for (i = 0; i < gauss->np; i++)
  {
	NcmMSetFunc *mean_func = g_ptr_array_index (fa, i % fa->len);
	const gdouble x_i = ncm_vector_get (gauss->x, i);
	const gdouble y_i = ncm_vector_get (gauss->y, i);
	const gdouble yt_i = ncm_mset_func_eval1 (mean_func, mset, x_i);
	ncm_vector_set (v, i, (yt_i - y_i));
  }

  wmean = gsl_stats_wmean (ncm_vector_gsl (gauss->weight)->data,
                           ncm_vector_gsl (gauss->weight)->stride,
                           ncm_vector_gsl (v)->data,
                           ncm_vector_gsl (v)->stride,
                           ncm_vector_gsl (v)->size);

  for (i = 0; i < gauss->np; i++)
  {
	const gdouble sigma_i = ncm_vector_get (gauss->sigma, i);
	const gdouble v_i = ncm_vector_get (v, i);
	ncm_vector_set (v, i, (v_i - wmean) / sigma_i);
  }
}

static void
_gauss_wmean_x_calc_m2lnL_val (NcmMSet *mset, gpointer model, gpointer data, gdouble *m2lnL)
{
  GPtrArray *fa = (GPtrArray *) model;
  NcDataGaussian *gauss = (NcDataGaussian *) data;
  //  gdouble wmean = 0.0;
  gdouble wt = 0.0;
  gdouble tmp2 = 0.0;
  gint i;

  *m2lnL = 0.0;

  for (i = 0; i < gauss->np; i++)
  {
	NcmMSetFunc *mean_func = g_ptr_array_index (fa, i % fa->len);
	const gdouble x_i = ncm_vector_get (gauss->x, i);
	const gdouble y_i = ncm_vector_get (gauss->y, i);
	const gdouble yt_i = ncm_mset_func_eval1 (mean_func, mset, x_i);
	const gdouble r_i = (yt_i - y_i);
	const gdouble w_i = ncm_vector_get (gauss->weight, i);
	const gdouble r_i2 = r_i *  r_i;
	*m2lnL += r_i2 * w_i;
	tmp2 += r_i * w_i;
	wt += w_i;
  }

  *m2lnL -= tmp2 * tmp2 / wt;
}

/********************************************************************************************************/

static void
_gauss_x_calc_leastsquares_f (NcmMSet *mset, gpointer model, gpointer data, NcmVector *v)
{
  GPtrArray *fa = (GPtrArray *) model;
  NcDataGaussian *gauss = (NcDataGaussian *) data;
  gint i;

  for (i = 0; i < gauss->np; i++)
  {
	NcmMSetFunc *mean_func = g_ptr_array_index (fa, i % fa->len);
	const gdouble x_i = ncm_vector_get (gauss->x, i);
	const gdouble y_i = ncm_vector_get (gauss->y, i);
	const gdouble yt_i = ncm_mset_func_eval1 (mean_func, mset, x_i);
	const gdouble sigma_i = ncm_vector_get (gauss->sigma, i);
	ncm_vector_set (v, i, (yt_i - y_i) / sigma_i);
  }
}

static void
_gauss_x_calc_m2lnL_val (NcmMSet *mset, gpointer model, gpointer data, gdouble *m2lnL)
{
  GPtrArray *fa = (GPtrArray *) model;
  NcDataGaussian *gauss = (NcDataGaussian *) data;
  gint i;

  *m2lnL = 0.0;

  for (i = 0; i < gauss->np; i++)
  {
	NcmMSetFunc *mean_func = g_ptr_array_index (fa, i % fa->len);
	const gdouble x_i = ncm_vector_get (gauss->x, i);
	const gdouble y_i = ncm_vector_get (gauss->y, i);
	const gdouble yt_i = ncm_mset_func_eval1 (mean_func, mset, x_i);
	const gdouble sigma_i = ncm_vector_get (gauss->sigma, i);
	const gdouble r_i = (yt_i - y_i) / sigma_i;
	*m2lnL += r_i * r_i;
  }
}

/********************************************************************************************************/

static void
_gauss_calc_leastsquares_f (NcmMSet *mset, gpointer model, gpointer data, NcmVector *v)
{
  GPtrArray *fa = (GPtrArray *) model;
  NcDataGaussian *gauss = (NcDataGaussian *) data;
  gint i;

  for (i = 0; i < gauss->np; i++)
  {
	NcmMSetFunc *mean_func = g_ptr_array_index (fa, i % fa->len);
	const gdouble y_i = ncm_vector_get (gauss->y, i);
	const gdouble yt_i = ncm_mset_func_eval0 (mean_func, mset);
	const gdouble sigma_i = ncm_vector_get (gauss->sigma, i);

	ncm_vector_set (v, i, (yt_i - y_i) / sigma_i);
  }
}

static void
_gauss_calc_m2lnL_val (NcmMSet *mset, gpointer model, gpointer data, gdouble *m2lnL)
{
  GPtrArray *fa = (GPtrArray *) model;
  NcDataGaussian *gauss = (NcDataGaussian *) data;
  gint i;

  *m2lnL = 0.0;

  for (i = 0; i < gauss->np; i++)
  {
	NcmMSetFunc *mean_func = g_ptr_array_index (fa, i % fa->len);
	const gdouble y_i = ncm_vector_get (gauss->y, i);
	const gdouble yt_i = ncm_mset_func_eval0 (mean_func, mset);
	const gdouble sigma_i = ncm_vector_get (gauss->sigma, i);
	const gdouble r_i = (yt_i - y_i) / sigma_i;
	*m2lnL += r_i * r_i;
  }
}

/********************************************************************************************************/

static void
_gauss_cov_calc_leastsquares_f (NcmMSet *mset, gpointer model, gpointer data, NcmVector *v)
{
  GPtrArray *fa = (GPtrArray *) model;
  NcDataGaussian *gauss = (NcDataGaussian *) data;
  gint i;

  for (i = 0; i < gauss->np; i++)
  {
	NcmMSetFunc *mean_func = g_ptr_array_index (fa, i % fa->len);
	const gdouble y_i = ncm_vector_get (gauss->y, i);
	const gdouble yt_i = ncm_mset_func_eval0 (mean_func, mset);
	ncm_vector_set (gauss->vp, i, (yt_i - y_i));
  }

  {
	gint ret = gsl_blas_dgemv (CblasNoTrans, 1.0, NCM_MATRIX_GSL (gauss->orto), ncm_vector_gsl (gauss->vp), 0.0, ncm_vector_gsl (v));
	ncm_vector_div (v, gauss->sigma);
	NC_TEST_GSL_RESULT("_gauss_cov_calc_m2lnL_v", ret);
  }
}

static void
_gauss_cov_calc_m2lnL_val (NcmMSet *mset, gpointer model, gpointer data, gdouble *m2lnL)
{
  GPtrArray *fa = (GPtrArray *) model;
  NcDataGaussian *gauss = (NcDataGaussian *) data;
  gint i, j;

  *m2lnL = 0.0;

  for (i = 0; i < gauss->np; i++)
  {
	NcmMSetFunc *mean_func = g_ptr_array_index (fa, i % fa->len);
	const gdouble y_i = ncm_vector_get (gauss->y, i);
	const gdouble yt_i = ncm_mset_func_eval0 (mean_func, mset);
	ncm_vector_set (gauss->vp, i, (yt_i - y_i));
  }

  for (i = 0; i < gauss->np; i++)
  {
	const gdouble f_i = ncm_vector_get (gauss->vp, i);
	gdouble u_i = 0.0;
	for (j = 0; j < gauss->np; j++)
	  u_i += ncm_matrix_get (gauss->inv_cov, i, j) * ncm_vector_get (gauss->vp, j);
	*m2lnL += u_i * f_i;
  }
}

/********************************************************************************************************/

static void
_gauss_x_cov_calc_leastsquares_f (NcmMSet *mset, gpointer model, gpointer data, NcmVector *v)
{
  GPtrArray *fa = (GPtrArray *) model;
  NcDataGaussian *gauss = (NcDataGaussian *) data;
  gint i;

  for (i = 0; i < gauss->np; i++)
  {
	NcmMSetFunc *mean_func = g_ptr_array_index (fa, i % fa->len);
	const gdouble x_i = ncm_vector_get (gauss->x, i);
	const gdouble y_i = ncm_vector_get (gauss->y, i);
	const gdouble yt_i = ncm_mset_func_eval1 (mean_func, mset, x_i);
	ncm_vector_set (gauss->vp, i, (yt_i - y_i));
  }

  {
	gint ret = gsl_blas_dgemv (CblasNoTrans, 1.0, NCM_MATRIX_GSL (gauss->orto), ncm_vector_gsl (gauss->vp), 0.0, ncm_vector_gsl (v));
	ncm_vector_div (v, gauss->sigma);
	NC_TEST_GSL_RESULT("_gauss_cov_calc_m2lnL_v", ret);
  }
}

static void
_gauss_x_cov_calc_m2lnL_val (NcmMSet *mset, gpointer model, gpointer data, gdouble *m2lnL)
{
  GPtrArray *fa = (GPtrArray *) model;
  NcDataGaussian *gauss = (NcDataGaussian *) data;
  gint i, j;

  *m2lnL = 0.0;

  for (i = 0; i < gauss->np; i++)
  {
	NcmMSetFunc *mean_func = g_ptr_array_index (fa, i % fa->len);
	const gdouble x_i = ncm_vector_get (gauss->x, i);
	const gdouble y_i = ncm_vector_get (gauss->y, i);
	const gdouble yt_i = ncm_mset_func_eval1 (mean_func, mset, x_i);
	ncm_vector_set (gauss->vp, i, (yt_i - y_i));
  }

  for (i = 0; i < gauss->np; i++)
  {
	const gdouble f_i = ncm_vector_get (gauss->vp, i);
	gdouble u_i = 0.0;
	for (j = 0; j < gauss->np; j++)
	  u_i += ncm_matrix_get (gauss->inv_cov, i, j) * ncm_vector_get (gauss->vp, j);
	*m2lnL += u_i * f_i;
  }
}

/********************************************************************************************************/

static NcDataStruct *
_nc_data_struct_gaussian_new (void)
{
  NcDataStruct *dts = nc_data_struct_new ();
  dts->data       = NULL;
  dts->dup        = &_gauss_data_dup;
  dts->copy       = &_gauss_data_copy;
  dts->free       = &_gauss_data_free;
  dts->begin      = &_gauss_data_begin;
  dts->end        = &_gauss_data_end;
  dts->get_length = &_gauss_get_length;

  return dts;
}

static NcDataStruct *
_nc_data_struct_gaussian_cov_new (void)
{
  NcDataStruct *dts = nc_data_struct_new ();
  dts->data       = NULL;
  dts->dup        = &_gauss_cov_data_dup;
  dts->copy       = &_gauss_cov_data_copy;
  dts->free       = &_gauss_cov_data_free;
  dts->begin      = &_gauss_cov_data_begin;
  dts->end        = &_gauss_data_end;
  dts->get_length = &_gauss_get_length;

  return dts;
}

/********************************************************************************************************/

static NcData *
_nc_data_gaussian_new ()
{
  NcData *data = nc_data_new ();
  static gchar *name = "Gaussian";

  data->name                  = name;
  data->type                  = NC_DATA_GAUSSIAN_SIGMA;
  data->init                  = FALSE;
  data->model                 = NULL;
  data->model_init            = NULL;
  data->model_ref             = (NcDataRef)g_ptr_array_ref;
  data->model_free            = (NcDataFree)g_ptr_array_unref;
  data->dts                   = _nc_data_struct_gaussian_new ();
  data->resample              = NULL;
  data->prepare               = NULL;
  data->calc_leastsquares_f   = &_gauss_calc_leastsquares_f;
  data->calc_leastsquares_J   = NULL;
  data->calc_leastsquares_f_J = NULL;
  data->calc_m2lnL_val        = &_gauss_calc_m2lnL_val;
  data->calc_m2lnL_grad       = NULL;
  data->calc_m2lnL_val_grad   = NULL;

  return data;
}

static NcData *
_nc_data_gaussian_cov_new ()
{
  NcData *data = nc_data_new ();
  static gchar *name = "Gaussian - covariance matrix";

  data->name                  = name;
  data->type                  = NC_DATA_GAUSSIAN_COV;
  data->init                  = FALSE;
  data->model                 = NULL;
  data->model_init            = NULL;
  data->model_ref             = (NcDataRef)g_ptr_array_ref;
  data->model_free            = (NcDataFree)g_ptr_array_unref;
  data->dts                   = _nc_data_struct_gaussian_cov_new ();
  data->resample              = NULL;
  data->prepare               = NULL;
  data->calc_leastsquares_f   = &_gauss_cov_calc_leastsquares_f;
  data->calc_leastsquares_J   = NULL;
  data->calc_leastsquares_f_J = NULL;
  data->calc_m2lnL_val        = &_gauss_cov_calc_m2lnL_val;
  data->calc_m2lnL_grad       = NULL;
  data->calc_m2lnL_val_grad   = NULL;

  return data;
}

static NcData *
_nc_data_gaussian_x_new ()
{
  NcData *data = nc_data_new ();
  static gchar *name = "Gaussian - x";

  data->name                  = name;
  data->type                  = NC_DATA_GAUSSIAN_X_SIGMA;
  data->init                  = FALSE;
  data->model                 = NULL;
  data->model_init            = NULL;
  data->model_ref             = (NcDataRef)g_ptr_array_ref;
  data->model_free            = (NcDataFree)g_ptr_array_unref;
  data->dts                   = _nc_data_struct_gaussian_new ();
  data->resample              = &_gaussian_x_resample;
  data->prepare               = NULL;
  data->calc_leastsquares_f   = &_gauss_x_calc_leastsquares_f;
  data->calc_leastsquares_J   = NULL;
  data->calc_leastsquares_f_J = NULL;
  data->calc_m2lnL_val        = &_gauss_x_calc_m2lnL_val;
  data->calc_m2lnL_grad       = NULL;
  data->calc_m2lnL_val_grad   = NULL;

  return data;
}

static NcData *
_nc_data_gaussian_x_cov_new ()
{
  NcData *data = nc_data_new ();
  static gchar *name = "Gaussian x - covariance matrix";

  data->name                  = name;
  data->type                  = NC_DATA_GAUSSIAN_X_COV;
  data->init                  = FALSE;
  data->model                 = NULL;
  data->model_init            = NULL;
  data->model_ref             = (NcDataRef)g_ptr_array_ref;
  data->model_free            = (NcDataFree)g_ptr_array_unref;
  data->dts                   = _nc_data_struct_gaussian_cov_new ();
  data->resample              = NULL;
  data->prepare               = NULL;
  data->calc_leastsquares_f   = &_gauss_x_cov_calc_leastsquares_f;
  data->calc_leastsquares_J   = NULL;
  data->calc_leastsquares_f_J = NULL;
  data->calc_m2lnL_val        = &_gauss_x_cov_calc_m2lnL_val;
  data->calc_m2lnL_grad       = NULL;
  data->calc_m2lnL_val_grad   = NULL;
  return data;
}

static NcData *
_nc_data_gaussian_x_wmean_new ()
{
  NcData *data = nc_data_new ();
  static gchar *name = "Gaussian - x wmean";

  data->name                  = name;
  data->type                  = NC_DATA_GAUSSIAN_X_WMEAN;
  data->init                  = FALSE;
  data->model                 = NULL;
  data->model_init            = NULL;
  data->model_ref             = (NcDataRef)g_ptr_array_ref;
  data->model_free            = (NcDataFree)g_ptr_array_unref;
  data->dts                   = _nc_data_struct_gaussian_new ();
  data->resample              = &_gaussian_x_resample;
  data->prepare               = NULL;
  data->calc_leastsquares_f   = &_gauss_wmean_x_calc_leastsquares_f;
  data->calc_leastsquares_J   = NULL;
  data->calc_leastsquares_f_J = NULL;
  data->calc_m2lnL_val        = &_gauss_wmean_x_calc_m2lnL_val;
  data->calc_m2lnL_grad       = NULL;
  data->calc_m2lnL_val_grad   = NULL;

  return data;
}

/********************************************************************************************************/

/**
 * nc_data_gaussian_new:
 * @gauss_type: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
NcData *
nc_data_gaussian_new (NcDataGaussianType gauss_type)
{
  NcData *data;

  switch (gauss_type)
  {
	case NC_DATA_GAUSSIAN_SIGMA:
	  data = _nc_data_gaussian_new ();
	  break;
	case NC_DATA_GAUSSIAN_COV:
	  data = _nc_data_gaussian_cov_new ();
	  break;
	case NC_DATA_GAUSSIAN_X_SIGMA:
	  data = _nc_data_gaussian_x_new ();
	  break;
	case NC_DATA_GAUSSIAN_X_COV:
	  data = _nc_data_gaussian_x_cov_new ();
	  break;
	case NC_DATA_GAUSSIAN_X_WMEAN:
	  data = _nc_data_gaussian_x_wmean_new ();
	  break;
	default:
	  g_assert_not_reached ();
	  break;
  }
  data->type = gauss_type;
  return data;
}

/********************************************************************************************************/

/**
 * nc_data_gaussian_set_func_array:
 * @data: FIXME
 * @func: (array length=n): Array of functions
   * @n: FIXME
 *
 * FIXME
 *
 */
void
nc_data_gaussian_set_func_array (NcData *data, NcmMSetFunc **func, guint n)
{
  GPtrArray *fa = ncm_mset_func_array_new ();
  guint i;
  for (i = 0; i < n; i++)
	g_ptr_array_add (fa, (gpointer) ncm_mset_func_ref (func[i]));

  data->model = fa;
}

/**
 * nc_data_gaussian_set_func:
 * @data: a #NcData
 * @func: a #NcmMSetFunc
 *
 * FIXME
 *
 */
void
nc_data_gaussian_set_func (NcData *data, NcmMSetFunc *func)
{
  GPtrArray *fa = ncm_mset_func_array_new ();
  g_ptr_array_add (fa, (gpointer) ncm_mset_func_ref (func));
  data->model = fa;
}

/********************************************************************************************************/

/**
 * nc_data_gaussian_init_from_matrix:
 * @data: FIXME
 * @cm: FIXME
 * @extra_data: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
void
nc_data_gaussian_init_from_matrix (NcData *data, NcmMatrix *cm, NcDataStruct *extra_data)
{
  NcDataGaussian *gauss;
  guint ncol;
  guint nrow = NCM_MATRIX_NROWS (cm);

  switch (data->type)
  {
	case NC_DATA_GAUSSIAN_SIGMA:
	  gauss = _gauss_data_alloc (nrow, FALSE);
	  ncol = 2;
	  break;
	case NC_DATA_GAUSSIAN_COV:
	  gauss = _gauss_cov_data_alloc (nrow, FALSE);
	  ncol = 1 + nrow;
	  break;
	case NC_DATA_GAUSSIAN_X_SIGMA:
	  gauss = _gauss_data_alloc (nrow, TRUE);
	  ncol = 3;
	  break;
	case NC_DATA_GAUSSIAN_X_COV:
	  gauss = _gauss_cov_data_alloc (nrow, TRUE);
	  ncol = 2 + nrow;
	  break;
	case NC_DATA_GAUSSIAN_X_WMEAN:
	  gauss = _gauss_data_alloc (nrow, TRUE);
	  ncol = 3;
	  break;
	default:
	  g_assert_not_reached ();
	  break;
  }

  g_assert (NCM_MATRIX_NCOLS (cm) == ncol);
  NC_DATA_DATA (data) = gauss;
  gauss->extra_data = extra_data;

  ncol = 0;

  if ((data->type >= NC_DATA_GAUSSIAN_X_SIGMA) && (data->type <= NC_DATA_GAUSSIAN_X_WMEAN))
  {
	NcmVector *x_view = ncm_matrix_get_col (cm, ncol++);
	ncm_vector_memcpy (gauss->x, x_view);
	ncm_vector_free (x_view);
  }

  {
	NcmVector *y_view = ncm_matrix_get_col (cm, ncol++);
	ncm_vector_memcpy (gauss->y, y_view);
	ncm_vector_free (y_view);
  }

  if (data->type == NC_DATA_GAUSSIAN_COV || data->type == NC_DATA_GAUSSIAN_X_COV)
  {
	NcmMatrix *inv_cov_view = ncm_matrix_get_submatrix (cm, 0, ncol, gauss->np, gauss->np);
	ncm_matrix_memcpy (gauss->inv_cov, inv_cov_view);
	ncm_matrix_free (inv_cov_view);
  }
  else
  {
	NcmVector *sigma_view = ncm_matrix_get_col (cm, ncol++);
	ncm_vector_memcpy (gauss->sigma, sigma_view);
	ncm_vector_free (sigma_view);
  }

  nc_data_init (data);
}

#ifdef NUMCOSMO_HAVE_SQLITE3
/**
 * nc_data_gaussian_init_from_query: (skip)
 * @data: FIXME
 * @db: FIXME
 * @query: FIXME
 * @extra_data: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
void
nc_data_gaussian_init_from_query (NcData *data, sqlite3 *db, gchar *query, NcDataStruct *extra_data)
{
  NcDataGaussian *gauss;
  gint i, nrow, qncol, ncol, ret;
  gchar **res;
  gchar *err_str;

  g_assert (db != NULL);

  ret = sqlite3_get_table (db, query, &res, &nrow, &qncol, &err_str);
  if (ret != SQLITE_OK)
  {
	sqlite3_free_table (res);
	g_error ("Query error: %s", err_str);
  }

  switch (data->type)
  {
	case NC_DATA_GAUSSIAN_SIGMA:
	  gauss = _gauss_data_alloc (nrow, FALSE);
	  ncol = 2;
	  break;
	case NC_DATA_GAUSSIAN_COV:
	  gauss = _gauss_cov_data_alloc (nrow, FALSE);
	  ncol = 1 + nrow;
	  break;
	case NC_DATA_GAUSSIAN_X_SIGMA:
	  gauss = _gauss_data_alloc (nrow, TRUE);
	  ncol = 3;
	  break;
	case NC_DATA_GAUSSIAN_X_COV:
	  gauss = _gauss_cov_data_alloc (nrow, TRUE);
	  ncol = 2 + nrow;
	  break;
	case NC_DATA_GAUSSIAN_X_WMEAN:
	  gauss = _gauss_data_alloc (nrow, TRUE);
	  ncol = 3;
	  break;
	default:
	  g_assert_not_reached ();
	  break;
  }

  g_assert (ncol == qncol);
  NC_DATA_STRUCT_DATA (data->dts) = gauss;
  gauss->extra_data = extra_data;

  for (i = 0; i < nrow; i++)
  {
	gint j = 0;
	if ((data->type >= NC_DATA_GAUSSIAN_X_SIGMA) && (data->type <= NC_DATA_GAUSSIAN_X_WMEAN))
	  ncm_vector_set (gauss->x, i, atof (res[(i + 1) * qncol + j++]));
	ncm_vector_set (gauss->y, i, atof (res[(i + 1) * qncol + j++]));

	if (data->type == NC_DATA_GAUSSIAN_COV || data->type == NC_DATA_GAUSSIAN_X_COV)
	{
	  for (j = 0; j < nrow; j++)
		ncm_matrix_set (gauss->inv_cov, i, j, atof (res[(i + 1) * qncol + j]));
	}
	else
	  ncm_vector_set (gauss->sigma, i, atof (res[(i + 1) * qncol + j++]));
  }

  sqlite3_free_table (res);

  nc_data_init (data);

  return;
}
#endif /* NUMCOSMO_HAVE_SQLITE3 */

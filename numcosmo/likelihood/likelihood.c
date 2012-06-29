/***************************************************************************
 *            likelihood.c
 *
 *  Wed May 30 15:36:41 2007
 *  Copyright  2007  Sandro Dias Pinto Vitenti
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
 * SECTION:likelihood
 * @title: Likelihood
 * @short_description: FIXME
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include <numcosmo/numcosmo.h>

#include <string.h>
#include <glib.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_cdf.h>

G_DEFINE_BOXED_TYPE (NcLikelihood, nc_likelihood, nc_likelihood_copy, nc_likelihood_free);

/**
 * nc_likelihood_new:
 * @ds: a #NcDataSet.
 *
 * FIXME
 *
 * Returns: FIXME
 */
NcLikelihood *
nc_likelihood_new (NcDataSet *ds)
{
  NcLikelihood *lh;
  lh = g_slice_new (NcLikelihood);
  lh->ds = ds;
  lh->priors = NULL;
  lh->clone = FALSE;
  return lh;
}

/**
 * nc_likelihood_copy:
 * @lh_orig: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
NcLikelihood *
nc_likelihood_copy (NcLikelihood *lh_orig)
{
  NcLikelihood *lh = g_slice_new (NcLikelihood);
  lh->ds = nc_dataset_copy (lh_orig->ds);

  if (lh_orig->priors != NULL)
	lh->priors = g_list_copy (lh_orig->priors);
  else
	lh->priors = NULL;
  lh->clone = TRUE;

  return lh;
}

/**
 * nc_likelihood_free:
 * @lh: FIXME
 *
 * FIXME
 */
void
nc_likelihood_free (NcLikelihood *lh)
{
  if (lh->clone)
  {
	nc_dataset_free0 (lh->ds, FALSE);
  }
  if (lh->priors != NULL)
	g_list_free_full (lh->priors, (GDestroyNotify)&ncm_mset_func_free);
  g_slice_free (NcLikelihood, lh);
  return;
}

/**
 * nc_likelihood_priors_add:
 * @lh: FIXME
 * @prior: FIXME
 *
 * FIXME
 */
void
nc_likelihood_priors_add (NcLikelihood *lh, NcmMSetFunc *prior)
{
  lh->priors = g_list_append (lh->priors, ncm_mset_func_ref (prior));
}

/**
 * nc_likelihood_priors_length:
 * @lh: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gint
nc_likelihood_priors_length (NcLikelihood *lh)
{
  return g_list_length (lh->priors);
}

/**
 * nc_likelihood_has_leastsquares_J:
 * @lh: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gboolean
nc_likelihood_has_leastsquares_J (NcLikelihood *lh)
{
  GList *data_list;
  if (lh->ds->data_list != NULL)
  {
	data_list = g_list_first (lh->ds->data_list);
	while (data_list)
	{
	  NcData *data = (NcData *)data_list->data;
	  if (!data->calc_leastsquares_f)
		return FALSE;
	  data_list = g_list_next(data_list);
	}
  }
  else
	return FALSE;

  return TRUE;
}

/**
 * nc_likelihood_has_m2lnL_grad:
 * @lh: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gboolean
nc_likelihood_has_m2lnL_grad (NcLikelihood *lh)
{
  GList *data_list;
  if (lh->ds->data_list != NULL)
  {
	data_list = g_list_first (lh->ds->data_list);
	while (data_list)
	{
	  NcData *data = (NcData *)data_list->data;
	  if (!data->calc_m2lnL_grad)
		return FALSE;
	  data_list = g_list_next(data_list);
	}
  }
  else
	return FALSE;

  return TRUE;
}


/**
 * nc_likelihood_leastsquares_f: (skip)
 * @lh: a #NcLikelihood.
 * @mset: a #NcmMSet.
 * @f: a #NcmVector.
 *
 * FIXME
 */
void
nc_likelihood_leastsquares_f (NcLikelihood *lh, NcmMSet *mset, NcmVector *f)
{
  guint data_size = nc_dataset_get_n (lh->ds);
  guint priors_size = nc_likelihood_priors_length (lh);
  if (data_size)
  {
	NcmVector *data_f = ncm_vector_get_subvector (f, 0, data_size);
	nc_likelihood_data_leastsquares_f (lh, mset, data_f);
	ncm_vector_free (data_f);
  }
  if (priors_size)
  {
	NcmVector *priors_f = ncm_vector_get_subvector (f, data_size, priors_size);
	nc_likelihood_priors_leastsquares_f (lh, mset, priors_f);
	ncm_vector_free (priors_f);
  }
  return;
}

/**
 * nc_likelihood_data_leastsquares_f: (skip)
 * @lh: a #NcLikelihood.
 * @mset: a #NcmMSet.
 * @data_f: a #NcmVector.
 *
 * FIXME
 */
void
nc_likelihood_data_leastsquares_f (NcLikelihood *lh, NcmMSet *mset, NcmVector *data_f)
{
  GList *data_list;
  guint pos = 0;

  if (lh->ds->data_list != NULL)
  {
	data_list = g_list_first (lh->ds->data_list);

	while (data_list)
	{
	  NcData *data = (NcData *)data_list->data;
	  guint n = NC_DATA_LENGTH (data);

	  if (!data->calc_leastsquares_f)
		g_error ("Cannot calculate leastsquares vector f, NcData (%s) is not compatible", data->name);

	  NcmVector *data_f_i = ncm_vector_get_subvector (data_f, pos, n);

	  if (NC_DATA_HAS_PREPARE(data))
		NC_DATA_PREPARE (data, mset);

	  data->calc_leastsquares_f (mset, NC_DATA_MODEL(data), NC_DATA_DATA(data), data_f_i);

	  pos += n;
	  data_list = g_list_next(data_list);
	  ncm_vector_free (data_f_i);
	}
  }

  return;
}

/**
 * nc_likelihood_priors_leastsquares_f: (skip)
 * @lh: a #NcLikelihood.
 * @mset: a #NcmMSet.
 * @priors_f: a #NcmVector.
 *
 * FIXME
 */
void
nc_likelihood_priors_leastsquares_f (NcLikelihood *lh, NcmMSet *mset, NcmVector *priors_f)
{
  GList *plist;
  guint i = 0;

  plist = g_list_first (lh->priors);
  while (plist != NULL)
  {
	NcmMSetFunc *func = NCM_MSET_FUNC (plist->data);
	ncm_vector_set (priors_f, i++, ncm_mset_func_eval0 (func, mset));
	plist = g_list_next (plist);
  }

  return;
}

/**
 * nc_likelihood_leastsquares_J: (skip)
 * @lh: a #NcLikelihood.
 * @mset: a #NcmMSet.
 * @J: a #NcmMatrix.
 *
 * FIXME
 */
void
nc_likelihood_leastsquares_J (NcLikelihood *lh, NcmMSet *mset, NcmMatrix *J)
{
  gint pos = 0;
  GList *data_list;
  GList *plist;

  if (lh->ds->data_list != NULL)
  {
	data_list = g_list_first (lh->ds->data_list);

	while (data_list)
	{
	  NcData *data = (NcData *)data_list->data;
	  guint n = NC_DATA_LENGTH (data);
	  NcmMatrix *jac = ncm_matrix_get_submatrix (J, pos, 0, n, NCM_MATRIX_NCOLS (J));

	  if (NC_DATA_HAS_PREPARE(data))
		NC_DATA_PREPARE (data, mset);

	  if (!data->calc_leastsquares_f)
		g_error ("Cannot calculate leastsquares vector f, NcData (%s) is not compatible", data->name);

	  if (data->calc_leastsquares_J)
		data->calc_leastsquares_J (mset, NC_DATA_MODEL (data), NC_DATA_DATA (data), jac);
	  else
		g_error ("Cannot calculate leastsquares jacobian J, NcData (%s) is not compatible", data->name);

	  pos += n;
	  data_list = g_list_next(data_list);
	  ncm_matrix_free (jac);
	}
  }

  plist = g_list_first (lh->priors);
  while (plist != NULL)
  {
	g_assert_not_reached ();
	/*
	NcmMSetFunc *func = NCM_MSET_FUNC (plist->data);
	NcmVector *row = ncm_matrix_get_row (J, pos + i++);
	ncm_vector_free (row);
	*/
	plist = g_list_next (plist);
  }
  return;
}

/**
 * nc_likelihood_leastsquares_f_J: (skip)
 * @lh: a #NcLikelihood.
 * @mset: a #NcmMSet.
 * @f: a #NcmVector.
 * @J: a #NcmMatrix.
 *
 * FIXME
 */
void
nc_likelihood_leastsquares_f_J (NcLikelihood *lh, NcmMSet *mset, NcmVector *f, NcmMatrix *J)
{
  gint pos = 0;
  GList *data_list;
  GList *plist;

  if (lh->ds->data_list != NULL)
  {
	data_list = g_list_first (lh->ds->data_list);

	while (data_list)
	{
	  NcData *data = (NcData *)data_list->data;
	  guint n = NC_DATA_LENGTH (data);
	  NcmMatrix *jac = ncm_matrix_get_submatrix (J, pos, 0, n, NCM_MATRIX_NCOLS (J));
	  NcmVector *v = ncm_vector_get_subvector (f, pos, n);

	  if (NC_DATA_HAS_PREPARE(data))
		NC_DATA_PREPARE (data, mset);

	  if (!data->calc_leastsquares_f)
		g_error ("Cannot calculate leastsquares vector f, NcData (%s) is not compatible", data->name);

	  if (data->calc_leastsquares_f_J)
		data->calc_leastsquares_f_J (mset, NC_DATA_MODEL (data), NC_DATA_DATA (data), v, J);
	  else
	  {
		data->calc_leastsquares_f (mset, NC_DATA_MODEL (data), NC_DATA_DATA (data), v);
		if (data->calc_leastsquares_J)
		  data->calc_leastsquares_J (mset, NC_DATA_MODEL (data), NC_DATA_DATA (data), jac);
		else
		  g_error ("Cannot calculate leastsquares jacobian J, NcData (%s) is not compatible", data->name);
	  }

	  pos += n;
	  data_list = g_list_next (data_list);
	  ncm_matrix_free (jac);
	  ncm_vector_free (v);
	}
  }

  plist = g_list_first (lh->priors);
  while (plist != NULL)
  {
	g_assert_not_reached ();
	/*
	NcmMSetFunc *func = NCM_MSET_FUNC (plist->data);
	NcmVector *row = ncm_matrix_get_row (J, pos + i);
	ncm_vector_set (f, i + pos, ncm_mset_func_eval0 (func, mset));
	NCM_FUNC_CONST_DF (func, model, pt, row);
	i++;
	ncm_vector_free (row);
	*/
	plist = g_list_next (plist);
  }
  return;
}

/**
 * nc_likelihood_m2lnL_val:
 * @lh: a #NcLikelihood.
 * @mset: a #NcmMSet.
 * @m2lnL: (out): FIXME
   *
 * FIXME
 */
void
nc_likelihood_m2lnL_val (NcLikelihood *lh, NcmMSet *mset, gdouble *m2lnL)
{
  gdouble data_m2lnL;
  gdouble priors_m2lnL;
  nc_likelihood_data_m2lnL_val (lh, mset, &data_m2lnL);
  nc_likelihood_priors_m2lnL_val (lh, mset, &priors_m2lnL);
  *m2lnL = data_m2lnL + priors_m2lnL;
  return;
}

/**
 * nc_likelihood_data_m2lnL_val:
 * @lh: a #NcLikelihood.
 * @mset: a #NcmMSet.
 * @data_m2lnL: (out): FIXME
   *
 * FIXME
 */
void
nc_likelihood_data_m2lnL_val (NcLikelihood *lh, NcmMSet *mset, gdouble *data_m2lnL)
{
  GList *data_list;
  *data_m2lnL = 0.0;

  if (lh->ds->data_list != NULL)
  {
	data_list = g_list_first (lh->ds->data_list);

	while (data_list)
	{
	  NcData *data = (NcData *)data_list->data;
	  gdouble m2lnL_i;

	  if (NC_DATA_HAS_PREPARE (data))
		NC_DATA_PREPARE (data, mset);

	  if (data->calc_m2lnL_val)
		data->calc_m2lnL_val (mset, NC_DATA_MODEL (data), NC_DATA_DATA (data), &m2lnL_i);
	  else
		g_error ("Cannot calculate m2lnL, NcData (%s) is not compatible", data->name);

	  *data_m2lnL += m2lnL_i;
	  data_list = g_list_next (data_list);
	}
  }

  return;
}

/**
 * nc_likelihood_priors_m2lnL_val:
 * @lh: a #NcLikelihood.
 * @mset: a #NcmMSet.
 * @priors_m2lnL: (out): FIXME
   *
 * FIXME
 */
void
nc_likelihood_priors_m2lnL_val (NcLikelihood *lh, NcmMSet *mset, gdouble *priors_m2lnL)
{
  GList *plist;
  *priors_m2lnL = 0.0;

  plist = g_list_first (lh->priors);
  while (plist != NULL)
  {
	NcmMSetFunc *func = NCM_MSET_FUNC (plist->data);
	gdouble sqrt_m2lnL_i = ncm_mset_func_eval0 (func, mset);
	*priors_m2lnL += sqrt_m2lnL_i * sqrt_m2lnL_i;
	plist = g_list_next (plist);
  }

  return;
}

/**
 * nc_likelihood_m2lnL_grad: (skip)
 * @lh: a #NcLikelihood.
 * @mset: a #NcmMSet.
 * @grad: a #NcmVector.
 *
 * FIXME
 */
void
nc_likelihood_m2lnL_grad (NcLikelihood *lh, NcmMSet *mset, NcmVector *grad)
{
  GList *data_list;
  GList *plist;
  guint i;
  guint free_params_len = ncm_mset_fparams_len (mset);
  NcmVector *grad_i = ncm_vector_new (free_params_len);

  ncm_vector_set_zero (grad);

  if (lh->ds->data_list != NULL)
  {
	data_list = g_list_first (lh->ds->data_list);

	while (data_list)
	{
	  NcData *data = (NcData *)data_list->data;

	  if (NC_DATA_HAS_PREPARE(data))
		NC_DATA_PREPARE(data, mset);

	  if (data->calc_m2lnL_grad)
		data->calc_m2lnL_grad (mset, NC_DATA_MODEL (data), NC_DATA_DATA (data), grad_i);
	  else
		g_error ("Cannot calculate Grad m2lnL, NcData (%s) is not compatible", data->name);

	  ncm_vector_add (grad, grad_i);
	  data_list = g_list_next(data_list);
	}
  }

  plist = g_list_first (lh->priors);
  while (plist != NULL)
  {
	NcmMSetFunc *func = NCM_MSET_FUNC (plist->data);
	gdouble sqrt_m2lnL_i = ncm_mset_func_eval0 (func, mset);
	/*NCM_FUNC_CONST_DF (func, model, pt, grad_i);*/
	g_assert_not_reached ();

	for (i = 0; i < ncm_vector_len (grad); i++)
	  ncm_vector_set (grad, i, ncm_vector_get (grad, i) + 2.0 * sqrt_m2lnL_i * ncm_vector_get (grad_i, i));

	plist = g_list_next (plist);
  }

  ncm_vector_free (grad_i);

  return;
}

/**
 * nc_likelihood_m2lnL_val_grad: (skip)
 * @lh: a #NcLikelihood.
 * @mset: a #NcmMSet.
 * @m2lnL: (out): FIXME
 * @grad: a #NcmVector.
 *
 * FIXME
 */
void
nc_likelihood_m2lnL_val_grad (NcLikelihood *lh, NcmMSet *mset, gdouble *m2lnL, NcmVector *grad)
{
  GList *data_list;
  GList *plist;
  guint i;
  guint free_params_len = ncm_mset_fparams_len (mset);
  NcmVector *grad_i = ncm_vector_new (free_params_len);
  *m2lnL = 0.0;

  ncm_vector_set_zero (grad);

  if (lh->ds->data_list != NULL)
  {
	data_list = g_list_first (lh->ds->data_list);

	while (data_list)
	{
	  NcData *data = (NcData *)data_list->data;
	  gdouble m2lnL_i;

	  if (NC_DATA_HAS_PREPARE(data))
		NC_DATA_PREPARE (data, mset);

	  if (data->calc_m2lnL_val_grad)
		data->calc_m2lnL_val_grad (mset, NC_DATA_MODEL (data), NC_DATA_DATA (data), &m2lnL_i, grad_i);
	  else if (data->calc_m2lnL_val && data->calc_m2lnL_grad)
	  {
		data->calc_m2lnL_val (mset, NC_DATA_MODEL (data), NC_DATA_DATA (data), &m2lnL_i);
		data->calc_m2lnL_grad (mset, NC_DATA_MODEL (data), NC_DATA_DATA (data), grad_i);
	  }
	  else
		g_error ("Cannot calculate Grad m2lnL and m2lnL, NcData (%s) is not compatible", data->name);

	  *m2lnL += m2lnL_i;
	  ncm_vector_add (grad, grad_i);
	  data_list = g_list_next(data_list);
	}
  }

  plist = g_list_first (lh->priors);
  while (plist != NULL)
  {
	NcmMSetFunc *func = NCM_MSET_FUNC (plist->data);
	gdouble sqrt_m2lnL_i = ncm_mset_func_eval0 (func, mset);
	/* NCM_FUNC_CONST_DF (func, model, pt, grad_i); */
	g_assert_not_reached ();

	for (i = 0; i < ncm_vector_len (grad); i++)
	  ncm_vector_set (grad, i, ncm_vector_get (grad, i) + 2.0 * sqrt_m2lnL_i * ncm_vector_get (grad_i, i));
	*m2lnL += sqrt_m2lnL_i * sqrt_m2lnL_i;
	plist = g_list_next (plist);
  }

  ncm_vector_free (grad_i);

  return;
}

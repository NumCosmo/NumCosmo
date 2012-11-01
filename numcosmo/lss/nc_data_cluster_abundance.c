/***************************************************************************
 *            nc_data_cluster_abundance.c
 *
 *  Tue Apr  6 01:11:23 2010
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
 * SECTION:nc_data_cluster_abundance
 * @title: Cluster Abundance Data
 * @short_description: FIXME
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include <numcosmo/numcosmo.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <glib.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_const_mksa.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_histogram2d.h>
#include <glib/gstdio.h>
#ifdef NUMCOSMO_HAVE_CFITSIO
#include <fitsio.h>
#endif /* NUMCOSMO_HAVE_CFITSIO */

static void
_ca_data_copy (gpointer dest_ptr, gpointer src_ptr)
{
  NcDataClusterAbundance *dest = (NcDataClusterAbundance *) dest_ptr;
  NcDataClusterAbundance *src = (NcDataClusterAbundance *) src_ptr;

  dest->area_survey = src->area_survey;
  dest->np          = src->np;

  if (dest->z != NULL)
  {
	nc_cluster_redshift_free (dest->z);
	dest->z = NULL;
  }
  if (dest->m != NULL)
  {
	nc_cluster_mass_free (dest->m);
	dest->m = NULL;
  }

  if (src->z != NULL)
	dest->z = nc_cluster_redshift_ref (src->z);
  if (src->m != NULL)
	dest->m = nc_cluster_mass_ref (src->m);

#define _VECTOR_COPY(name) \
  if (src->name != NULL) \
  { \
	if (dest->name == NULL) \
	  dest->name = ncm_vector_copy (src->name); \
	if (ncm_vector_len (dest->name) != ncm_vector_len (src->name)) \
	{ \
	  ncm_vector_free (dest->name); \
	  dest->name = ncm_vector_copy (src->name); \
	} \
	ncm_vector_memcpy (dest->name, src->name); \
  }

#define _MATRIX_COPY(name) \
  if (src->name != NULL) \
  { \
	if (dest->name == NULL) \
	  dest->name = ncm_matrix_copy (src->name); \
	if ((NCM_MATRIX_NROWS (dest->name) != NCM_MATRIX_NROWS (src->name)) || (NCM_MATRIX_NCOLS (dest->name) != NCM_MATRIX_NCOLS (src->name))) \
	{ \
	  ncm_matrix_free (dest->name); \
	  dest->name = ncm_matrix_copy (src->name); \
	} \
	ncm_matrix_memcpy (dest->name, src->name); \
  }

  _VECTOR_COPY (lnM_true)
  _VECTOR_COPY (z_true)

  _MATRIX_COPY (z_obs)
  _MATRIX_COPY (z_obs_params)
  _MATRIX_COPY (lnM_obs)
  _MATRIX_COPY (lnM_obs_params)
}

static gpointer
_ca_data_dup (gpointer ca_ptr)
{
  NcDataClusterAbundance *ca = (NcDataClusterAbundance *) ca_ptr;
  NcDataClusterAbundance *clone = g_slice_new0 (NcDataClusterAbundance);
  _ca_data_copy (clone, ca);

  return clone;
}

static void
_ca_data_free (gpointer ca_ptr)
{
  NcDataClusterAbundance *ca = (NcDataClusterAbundance *) ca_ptr;

  if (ca->z != NULL)
	nc_cluster_redshift_free (ca->z);
  if (ca->m != NULL)
	nc_cluster_mass_free (ca->m);

  if (ca->lnM_true != NULL)
	ncm_vector_free (ca->lnM_true);
  if (ca->z_true != NULL)
	ncm_vector_free (ca->z_true);

  if (ca->z_obs != NULL)
	ncm_matrix_free (ca->z_obs);
  if (ca->z_obs_params != NULL)
	ncm_matrix_free (ca->z_obs_params);
  if (ca->lnM_obs != NULL)
	ncm_matrix_free (ca->lnM_obs);
  if (ca->lnM_obs_params != NULL)
	ncm_matrix_free (ca->lnM_obs_params);

  g_slice_free (NcDataClusterAbundance, ca_ptr);
}

static void
_ca_data_begin (gpointer ca_ptr)
{
  NcDataClusterAbundance *ca = (NcDataClusterAbundance *) ca_ptr;
  ca->log_np_fac = lgamma (ca->np + 1);
}

static guint _ca_data_get_length (gpointer ca_ptr) { return ((NcDataClusterAbundance *) ca_ptr)->np; }

/**
 * FIXME
 */
static NcDataStruct *
_nc_data_struct_cluster_abundance_new ()
{
  NcDataStruct *dts = nc_data_struct_new ();
  NcDataClusterAbundance *dca = g_slice_new0 (NcDataClusterAbundance);

  dts->data       = dca;
  dts->dup        = &_ca_data_dup;
  dts->copy       = &_ca_data_copy;
  dts->free       = &_ca_data_free;
  dts->begin      = &_ca_data_begin;
  dts->end        = NULL;
  dts->get_length = &_ca_data_get_length;

  return dts;
}

static void _nc_data_cluster_abundance_unbinned_prepare (NcmMSet *mset, gpointer model, gpointer data);
static void _nc_data_cluster_abundance_binned_resample (NcmMSet *mset, gpointer model, gpointer data);

static void
_nc_data_cluster_abundance_binned_prepare (NcmMSet *mset, gpointer model, gpointer data)
{
  NcDataPoisson *poisson = (NcDataPoisson *) data;
  NcmMSetFunc *int_func = NCM_MSET_FUNC (model);

  _nc_data_cluster_abundance_unbinned_prepare (mset, int_func->obj, NC_DATA_STRUCT_DATA (poisson->extra_data));
}

/**
 * nc_data_cluster_abundance_binned_new:
 * @cad: a #NcClusterAbundance.
 *
 * FIXME
 *
 * Returns: FIXME
 */
NcData *
nc_data_cluster_abundance_binned_new (NcClusterAbundance *cad)
{
  NcData *data = nc_data_poisson_new (NC_DATA_POISSON_INT);
  NC_DATA_MODEL (data) = nc_data_cluster_abundance_binned_new_function (cad);

  nc_data_poisson_set_prepare (data, &_nc_data_cluster_abundance_binned_prepare);
  nc_data_poisson_set_resample (data, &_nc_data_cluster_abundance_binned_resample);

  return data;
}

/**
 * nc_data_cluster_abundance_init_from_fits_file:
 * @data: a #NcData
 * @filename: name of the file
 *
 * FIXME
 */
void
nc_data_cluster_abundance_init_from_fits_file (NcData *data, gchar *filename)
{
  g_assert_not_reached ();
}

static void
_nc_data_cluster_abundance_binned_f (NcmMSet *mset, gpointer obj, const gdouble *x, gdouble *f)
{
  //NcHICosmo *model = NC_HICOSMO (ncm_mset_peek (mset, NC_HICOSMO_ID));
  //NcClusterAbundance *cad = NC_CLUSTER_ABUNDANCE (obj);
  //f[0] = nc_cluster_abundance_N_val (cad, model, cad->lnMi, cad->lnMf, cad->zi, x[0]);
  g_assert_not_reached ();
  return;
}

/**
 * nc_data_cluster_abundance_binned_new_function: (skip)
 * @cad: a #NcClusterAbundance
 *
 * FIXME
 *
 * Returns: FIXME
 */
NcmMSetFunc *
nc_data_cluster_abundance_binned_new_function (NcClusterAbundance *cad)
{
  NcmMSetFunc *func = ncm_mset_func_new (_nc_data_cluster_abundance_binned_f, 1, 1, cad, (GDestroyNotify)nc_cluster_abundance_free);
  return func;
}

/**
 * nc_data_cluster_abundance_binned_lnM_z_new:
 * @cad: a #NcClusterAbundance
 *
 * FIXME
 *
 * Returns: FIXME
 */
NcData *
nc_data_cluster_abundance_binned_lnM_z_new (NcClusterAbundance *cad)
{
  g_assert_not_reached ();
  return NULL;
}

/************************************************************************************************************
 * Unbinned abundance data                                                                                  *
 ************************************************************************************************************/

static void
_nc_data_cluster_abundance_unbinned_model_init (gpointer model, gpointer data)
{
  NcClusterAbundance *cad = NC_CLUSTER_ABUNDANCE (model);
  NcDataClusterAbundance *dca = (NcDataClusterAbundance *) data;

  if (dca->z == NULL || dca->m == NULL)
	g_error ("Cannot init NcClusterAbundance missing NcClusterRedshift or NcClusterMass object.");

  nc_cluster_abundance_set_redshift (cad, dca->z);
  nc_cluster_abundance_set_mass (cad, dca->m);

  cad->completeness  = dca->completeness;
  cad->purity        = dca->purity;
  cad->sd_lnM        = dca->sd_lnM;

  cad->mfp->area_survey = dca->area_survey;
}

static void
_nc_data_cluster_abundance_unbinned_prepare (NcmMSet *mset, gpointer model, gpointer data)
{
  NcClusterAbundance *cad = NC_CLUSTER_ABUNDANCE (model);
  NcHICosmo *hic = NC_HICOSMO (ncm_mset_peek (mset, NC_HICOSMO_ID));
  /* NcDataNClusterAbundance *dca = (NcDataNClusterAbundance *) data; FIXME */

  if (ncm_model_ctrl_update (cad->ctrl, NCM_MODEL (hic)))
  {
	nc_cluster_abundance_prepare (cad, hic);
  }
}

/**
 * _nc_data_cluster_abundance_resample
 * @model: a #NcmModel.
 * @model: FIXME
 * @data: FIXME
 *
 * This function generates random numbers which are used to obtain redshift
 * and mass (logarithm base e) values...
   *
 */
static void
_nc_data_cluster_abundance_resample (NcmMSet *mset, gpointer model, gpointer data)
{
  NcClusterAbundance *cad = NC_CLUSTER_ABUNDANCE (model);
  NcDataClusterAbundance *dca = (NcDataClusterAbundance *) data;
  NcHICosmo *mc = NC_HICOSMO (ncm_mset_peek (mset, NC_HICOSMO_ID));
  gsl_rng *rng = ncm_get_rng ();
  guint z_obs_len = nc_cluster_redshift_obs_len (cad->z);
  guint z_obs_params_len = nc_cluster_redshift_obs_params_len (cad->z);
  guint lnM_obs_len = nc_cluster_mass_obs_len (cad->m);
  guint lnM_obs_params_len = nc_cluster_mass_obs_params_len (cad->m);
  GArray *lnM_true_array = NULL;
  GArray *z_true_array = NULL;
  GArray *z_obs_array = NULL;
  GArray *z_obs_params_array = NULL;
  GArray *lnM_obs_array = NULL;
  GArray *lnM_obs_params_array = NULL;
  guint total_np;
  gint i;

  gdouble *zi_obs = g_new (gdouble, z_obs_len);
  gdouble *zi_obs_params = z_obs_params_len > 0 ? g_new (gdouble, z_obs_params_len) : NULL;
  gdouble *lnMi_obs = g_new (gdouble, lnM_obs_len);
  gdouble *lnMi_obs_params = lnM_obs_params_len > 0 ? g_new (gdouble, lnM_obs_params_len) : NULL;

  total_np = gsl_ran_poisson (rng, cad->norma);

  lnM_true_array = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), total_np);
  z_true_array = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), total_np);

  z_obs_array = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), total_np * z_obs_len);
  if (z_obs_params_len > 0)
	z_obs_params_array = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), total_np * z_obs_params_len);

  lnM_obs_array = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), total_np * lnM_obs_len);
  if (lnM_obs_params_len > 0)
	lnM_obs_params_array = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), total_np * lnM_obs_params_len);

  //printf ("# Generating unbinned %u (z,lnM) %g\n", total_np, cad->norma);
  //printf ("# Resampling in range [% 20.15g, % 20.15g] [% 20.15e, % 20.15e]\n", cad->zi, cad->zf, exp (cad->lnMi), exp (cad->lnMf));
  nc_cluster_abundance_prepare_inv_dNdz (cad, NC_HICOSMO (ncm_mset_peek (mset, NC_HICOSMO_ID)));

  for (i = 0; i < total_np; i++)
  {
	const gdouble u1 = _nc_cad_inv_dNdz_convergence_f (gsl_rng_uniform_pos (rng), cad->z_epsilon);
	const gdouble u2 = _nc_cad_inv_dNdz_convergence_f (gsl_rng_uniform_pos (rng), cad->lnM_epsilon);
	const gdouble z_true = ncm_spline_eval (cad->inv_z, u1);
	const gdouble lnM_true = ncm_spline2d_eval (cad->inv_lnM_z, u2, z_true);

    if ( nc_cluster_redshift_resample (cad->z, lnM_true, z_true, zi_obs, zi_obs_params) &&
         nc_cluster_mass_resample (cad->m, mc, lnM_true, z_true, lnMi_obs, lnMi_obs_params) )
	{
	  g_array_append_val (lnM_true_array, lnM_true);
	  g_array_append_val (z_true_array, z_true);

	  g_array_append_vals (z_obs_array, zi_obs, z_obs_len);
	  g_array_append_vals (lnM_obs_array, lnMi_obs, lnM_obs_len);

	  if (z_obs_params_len > 0)
		g_array_append_vals (z_obs_params_array, zi_obs_params, z_obs_params_len);
	  if (lnM_obs_params_len > 0)
		g_array_append_vals (lnM_obs_params_array, lnMi_obs_params, lnM_obs_params_len);
	}
	//printf ("% 20.15g % 20.15g\n", zi_real, lnMi_real);
  }

  if (dca->lnM_true != NULL)
	ncm_vector_free (dca->lnM_true);
  dca->lnM_true = ncm_vector_new_array (lnM_true_array);
  g_array_unref (lnM_true_array);

  if (dca->z_true != NULL)
	ncm_vector_free (dca->z_true);
  dca->z_true = ncm_vector_new_array (z_true_array);
  g_array_unref (z_true_array);

  if (dca->z_obs != NULL)
	ncm_matrix_free (dca->z_obs);
  dca->z_obs = ncm_matrix_new_array (z_obs_array, z_obs_len);
  g_array_unref (z_obs_array);

  if (dca->lnM_obs != NULL)
	ncm_matrix_free (dca->lnM_obs);
  dca->lnM_obs = ncm_matrix_new_array (lnM_obs_array, lnM_obs_len);
  g_array_unref (lnM_obs_array);

  if (z_obs_params_len > 0)
  {
	if (dca->z_obs_params != NULL)
	  ncm_matrix_free (dca->z_obs_params);
	dca->z_obs_params = ncm_matrix_new_array (z_obs_params_array, z_obs_params_len);
	g_array_unref (z_obs_params_array);
  }

  if (lnM_obs_params_len > 0)
  {
	if (dca->lnM_obs_params != NULL)
	  ncm_matrix_free (dca->lnM_obs_params);
	dca->lnM_obs_params = ncm_matrix_new_array (lnM_obs_params_array, lnM_obs_params_len);
	g_array_unref (lnM_obs_params_array);
  }

  dca->np = NCM_MATRIX_NROWS (dca->z_obs);
  //printf ("# Generated %ld | expected % 20.15g\n", dca->np, nc_cluster_abundance_n (cad, NC_HICOSMO (ncm_mset_peek (mset, NC_HICOSMO_ID))));

  g_free (zi_obs);
  g_free (zi_obs_params);
  g_free (lnMi_obs);
  g_free (lnMi_obs_params);
}

static void
_nc_data_cluster_abundance_binned_resample (NcmMSet *mset, gpointer model, gpointer data)
{
  NcmMSetFunc *int_func = NCM_MSET_FUNC (model);
  NcDataPoisson *poisson = (NcDataPoisson *) data;
  NcDataClusterAbundance *dca = (NcDataClusterAbundance *) NC_DATA_STRUCT_DATA (poisson->extra_data);
  NcClusterAbundance *cad = NC_CLUSTER_ABUNDANCE (int_func->obj);
  gsl_rng *rng = ncm_get_rng ();
  static gdouble u2t = 0.0;
  gint i;

  nc_cluster_abundance_prepare_inv_dNdz (cad, NC_HICOSMO (ncm_mset_peek (mset, NC_HICOSMO_ID)));
  gsl_histogram_reset (poisson->h);

  dca->np = gsl_ran_poisson (rng, cad->norma);

  //printf ("# Generating binned %ld (z,lnM) %g\n", dca->np, cad->norma);

  for (i = 0; i < dca->np; i++)
  {
	const gdouble u1 = _nc_cad_inv_dNdz_convergence_f (gsl_rng_uniform_pos (rng), cad->z_epsilon); //gsl_rng_uniform (rng);
	const gdouble u2 = gsl_rng_uniform_pos (rng);
	const gdouble zi = ncm_spline_eval (cad->inv_z, u1);
	//printf ("% 20.15g % 20.15g % 20.15g % 20.15g\n", u1, 0.0, zi, 0.0);
	gsl_histogram_increment (poisson->h, zi);
	u2t += u2;
  }
}

typedef struct
{
  NcClusterAbundance *cad;
  NcDataClusterAbundance *dca;
  NcHICosmo *m;
  gdouble *m2lnL;
} _Evald2N;

static void
_eval_z_p_lnm_p_d2n (glong i, glong f, gpointer data)
{
  _Evald2N *evald2n = (_Evald2N *) data;
  glong n;
  gdouble m2lnL = 0.0;
  G_LOCK_DEFINE_STATIC (save_m2lnL);

  for (n = i; n < f; n++)
  {
	gdouble *lnMn_obs = ncm_matrix_ptr (evald2n->dca->lnM_obs, n, 0);
	gdouble *lnMn_obs_params = evald2n->dca->lnM_obs_params != NULL ? ncm_matrix_ptr (evald2n->dca->lnM_obs_params, n, 0) : NULL;
	gdouble *zn_obs = ncm_matrix_ptr (evald2n->dca->z_obs, n, 0);
	gdouble *zn_obs_params = evald2n->dca->z_obs_params != NULL ? ncm_matrix_ptr (evald2n->dca->z_obs_params, n, 0) : NULL;
	const gdouble mlnLn = -log (nc_cluster_abundance_z_p_lnm_p_d2n (evald2n->cad, evald2n->m, lnMn_obs, lnMn_obs_params, zn_obs, zn_obs_params));
	m2lnL += mlnLn;
  }

  G_LOCK (save_m2lnL);
  *evald2n->m2lnL += m2lnL;
  G_UNLOCK (save_m2lnL);
}

static void
_eval_z_p_d2n (glong i, glong f, gpointer data)
{
  _Evald2N *evald2n = (_Evald2N *) data;
  glong n;
  gdouble m2lnL = 0.0;
  G_LOCK_DEFINE_STATIC (save_m2lnL);

  for (n = i; n < f; n++)
  {
	const gdouble lnMn = ncm_vector_get (evald2n->dca->lnM_true, n);
	gdouble *zn_obs = ncm_matrix_ptr (evald2n->dca->z_obs, n, 0);
	gdouble *zn_obs_params = evald2n->dca->z_obs_params != NULL ? ncm_matrix_ptr (evald2n->dca->z_obs_params, n, 0) : NULL;
	const gdouble mlnLn = -log (nc_cluster_abundance_z_p_d2n (evald2n->cad, evald2n->m, lnMn, zn_obs, zn_obs_params));
	m2lnL += mlnLn;
  }

  G_LOCK (save_m2lnL);
  *evald2n->m2lnL += m2lnL;
  G_UNLOCK (save_m2lnL);
}

static void
_eval_lnm_p_d2n (glong i, glong f, gpointer data)
{
  _Evald2N *evald2n = (_Evald2N *) data;
  glong n;
  gdouble m2lnL = 0.0;
  G_LOCK_DEFINE_STATIC (save_m2lnL);

  for (n = i; n < f; n++)
  {
	const gdouble zn = ncm_vector_get (evald2n->dca->z_true, n);
	gdouble *lnMn_obs = ncm_matrix_ptr (evald2n->dca->lnM_obs, n, 0);
	gdouble *lnMn_obs_params = evald2n->dca->lnM_obs_params != NULL ? ncm_matrix_ptr (evald2n->dca->lnM_obs_params, n, 0) : NULL;
	const gdouble mlnLn = -log (nc_cluster_abundance_lnm_p_d2n (evald2n->cad, evald2n->m, lnMn_obs, lnMn_obs_params, zn));
	m2lnL += mlnLn;
  }

  G_LOCK (save_m2lnL);
  *evald2n->m2lnL += m2lnL;
  G_UNLOCK (save_m2lnL);
}

static void
_eval_d2n (glong i, glong f, gpointer data)
{
  _Evald2N *evald2n = (_Evald2N *) data;
  glong n;
  gdouble m2lnL = 0.0;
  G_LOCK_DEFINE_STATIC (save_m2lnL);

  for (n = i; n < f; n++)
  {
	const gdouble lnMn = ncm_vector_get (evald2n->dca->lnM_true, n);
	const gdouble zn = ncm_vector_get (evald2n->dca->z_true, n);
	const gdouble mlnLn = -log (nc_cluster_abundance_d2n (evald2n->cad, evald2n->m, lnMn, zn));
	m2lnL += mlnLn;
  }

  G_LOCK (save_m2lnL);
  *evald2n->m2lnL += m2lnL;
  G_UNLOCK (save_m2lnL);
}

static void
_eval_intp_d2n (glong i, glong f, gpointer data)
{
  _Evald2N *evald2n = (_Evald2N *) data;
  glong n;
  gdouble m2lnL = 0.0;
  G_LOCK_DEFINE_STATIC (save_m2lnL);

  for (n = i; n < f; n++)
  {
	const gdouble lnMn = ncm_vector_get (evald2n->dca->lnM_true, n);
	const gdouble zn = ncm_vector_get (evald2n->dca->z_true, n);
	const gdouble mlnLn = -log (nc_cluster_abundance_intp_d2n (evald2n->cad, evald2n->m, lnMn, zn));
	m2lnL += mlnLn;
  }

  G_LOCK (save_m2lnL);
  *evald2n->m2lnL += m2lnL;
  G_UNLOCK (save_m2lnL);
}

static void
_nc_data_cluster_abundance_calc_m2lnL (NcmMSet *mset, gpointer model, gpointer data, gdouble *m2lnL)
{
  NcClusterAbundance *cad = NC_CLUSTER_ABUNDANCE (model);
  NcDataClusterAbundance *dca = (NcDataClusterAbundance *) data;
  NcHICosmo *m = NC_HICOSMO (ncm_mset_peek (mset, NC_HICOSMO_ID));
  GTimer *bench = g_timer_new ();

  *m2lnL = 0.0;

  if (dca->use_true_data)
  {
	_Evald2N evald2n = {cad, dca, m, m2lnL};
	g_assert (dca->z_true);
	g_assert (dca->lnM_true);
	ncm_func_eval_threaded_loop (&_eval_intp_d2n, 0, dca->np, &evald2n);
  }
  else
  {
	NcClusterRedshiftImpl z_impl = nc_cluster_redshift_impl (dca->z);
	NcClusterMassImpl lnM_impl = nc_cluster_mass_impl (dca->m);
	gboolean z_p = z_impl & NC_CLUSTER_REDSHIFT_P;
	gboolean lnM_p = lnM_impl & NC_CLUSTER_MASS_P;
	if (z_p && lnM_p)
	{
	  _Evald2N evald2n = {cad, dca, m, m2lnL};
	  ncm_func_eval_threaded_loop (&_eval_z_p_lnm_p_d2n, 0, dca->np, &evald2n);
	}
	else if (z_p && !lnM_p)
	{
	  g_assert (dca->lnM_true);
	  _Evald2N evald2n = {cad, dca, m, m2lnL};
	  ncm_func_eval_threaded_loop (&_eval_z_p_d2n, 0, dca->np, &evald2n);
	}
	else if (!z_p && lnM_p)
	{
	  g_assert (dca->z_true);
	  _Evald2N evald2n = {cad, dca, m, m2lnL};
	  ncm_func_eval_threaded_loop (&_eval_lnm_p_d2n, 0, dca->np, &evald2n);
	}
	else
	{
	  g_assert (dca->z_true);
	  g_assert (dca->lnM_true);
       _Evald2N evald2n = {cad, dca, m, m2lnL};
	  ncm_func_eval_threaded_loop (&_eval_d2n, 0, dca->np, &evald2n);

	}
  }

  *m2lnL += (dca->log_np_fac + nc_cluster_abundance_n (cad, m));

  *m2lnL *= 2.0;

  g_timer_destroy (bench);
}

/**
 * nc_data_cluster_abundance_unbinned_new:
 * @cad: a #NcClusterAbundance
 *
 * FIXME
 *
 */
NcData *
nc_data_cluster_abundance_unbinned_new (NcClusterAbundance *cad)
{
  static gchar *name = "Cluster abundance unbinned";
  NcData *data = nc_data_new ();

  data->name                  = name;
  data->type                  = 0;
  data->init                  = FALSE;
  data->model                 = nc_cluster_abundance_ref (cad);
  data->model_init            = &_nc_data_cluster_abundance_unbinned_model_init;
  data->model_ref             = (NcDataRef) &nc_cluster_abundance_ref;
  data->model_free            = (NcDataFree) &nc_cluster_abundance_free;
  data->dts                   = _nc_data_struct_cluster_abundance_new ();
  data->prepare               = &_nc_data_cluster_abundance_unbinned_prepare;

  data->calc_leastsquares_f   = NULL;
  data->calc_leastsquares_J   = NULL;
  data->calc_leastsquares_f_J = NULL;

  data->resample              = &_nc_data_cluster_abundance_resample;
  data->calc_m2lnL_val        = &_nc_data_cluster_abundance_calc_m2lnL;

  data->calc_m2lnL_grad       = NULL;
  data->calc_m2lnL_val_grad   = NULL;

  return data;
}

/**
 * nc_data_cluster_abundance_true_data:
 * @data: a #NcData.
 * @use_true_data: FIXME
 *
 * FIXME
 *
 */
void
nc_data_cluster_abundance_true_data (NcData *data, gboolean use_true_data)
{
  NcDataClusterAbundance *dca = (NcDataClusterAbundance *) NC_DATA_DATA (data);
  if (use_true_data)
  {
	g_assert (dca->lnM_true != NULL);
	g_assert (dca->z_true != NULL);
  }
  dca->use_true_data = use_true_data;
}

/**
 * nc_data_cluster_abundance_unbinned_init_from_sampling:
 * @data: a #NcData.
 * @mset: a #NcmMSet.
 * @clusterz: a #NcClusterRedshift.
 * @clusterm: a #NcClusterMass.
 * @area_survey: area in units of square degrees.
 *
 * FIXME
 *
 */
void
nc_data_cluster_abundance_unbinned_init_from_sampling (NcData *data, NcmMSet *mset, NcClusterRedshift *clusterz, NcClusterMass *clusterm, gdouble area_survey)
{
  NcDataClusterAbundance *dca = (NcDataClusterAbundance *) NC_DATA_DATA (data);

  dca->area_survey = area_survey;

  dca->z = nc_cluster_redshift_ref (clusterz);
  dca->m = nc_cluster_mass_ref (clusterm);

  nc_data_model_init (data);
  nc_data_resample (data, mset, FALSE);

  nc_data_init (data);
}

/**
 * nc_data_cluster_abundance_unbinned_bin_data: (skip)
 * @ca_unbinned: a #NcData
 * @nodes: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
NcData *
nc_data_cluster_abundance_unbinned_bin_data (NcData *ca_unbinned, gsl_vector *nodes)
{
  NcDataClusterAbundance *dca = (NcDataClusterAbundance *) NC_DATA_DATA (ca_unbinned);
  NcClusterAbundance *cad = NC_CLUSTER_ABUNDANCE (NC_DATA_MODEL (ca_unbinned));
  NcData *ca_binned;
  NcDataStruct *dts_ca = _nc_data_struct_cluster_abundance_new ();
  gsl_histogram *hist;
  gint i;

  g_assert (ca_unbinned->init);
  g_assert (nodes->size > 1);
  g_assert (nodes->stride == 1);

  ca_binned = nc_data_cluster_abundance_binned_new (cad);

  hist = gsl_histogram_alloc (nodes->size - 1);
  gsl_histogram_set_ranges (hist, nodes->data, nodes->size);

  {
	for (i = 0; i < dca->np; i++)
	{
	  const gdouble z_i = 0.0;//gsl_matrix_get (dca->real.z_lnM, i, 0);
	  gsl_histogram_increment (hist, z_i);
	  g_assert_not_reached ();
	}
  }

  //nc_data_cluster_abundance_binned_init_from_hist (ca_binned, hist, dca->opt, dca->area_survey, dca->lnMi, dca->lnMf, dca->photoz_sigma0, dca->lnM_sigma0);
  nc_data_poisson_init_from_histogram (ca_binned, hist, FALSE, dts_ca);

  gsl_histogram_free (hist);

  return ca_binned;
}

/**
 * nc_data_cluster_abundance_hist_lnM_z: (skip)
 * @ca_unbinned: a #NcData.
 * @lnM_nodes: FIXME
 * @z_nodes: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gsl_histogram2d *
nc_data_cluster_abundance_hist_lnM_z (NcData *ca_unbinned, gsl_vector *lnM_nodes, gsl_vector *z_nodes)
{
  NcDataClusterAbundance *dca = (NcDataClusterAbundance *) NC_DATA_DATA (ca_unbinned);
  //NcClusterAbundance *cad = NC_CLUSTER_ABUNDANCE (NC_DATA_MODEL (ca_unbinned));
  gsl_histogram2d *hist;
  gint i;

  g_assert (ca_unbinned->init);
  g_assert (lnM_nodes->size > 1);
  g_assert (lnM_nodes->stride == 1);
  g_assert (z_nodes->size > 1);
  g_assert (z_nodes->stride == 1);

  //ca_binned = nc_data_cluster_abundance_binned_lnM_z_new (cad); /* I have to make this function. */
  //ca_binned = nc_data_cluster_abundance_unbinned_new (cad);

  hist = gsl_histogram2d_alloc (lnM_nodes->size - 1, z_nodes->size - 1);
  gsl_histogram2d_set_ranges (hist, lnM_nodes->data, lnM_nodes->size, z_nodes->data, z_nodes->size);

  {
	for (i = 0; i < dca->np; i++)
	{
	  const gdouble zi_real = 0.0;//gsl_matrix_get (dca->real.z_lnM, i, 0);
	  const gdouble lnMi_real = 0.0;//gsl_matrix_get (dca->real.z_lnM, i, 1);
	  g_assert_not_reached ();

	  gsl_histogram2d_increment (hist, lnMi_real, zi_real);
	}
  }

  return hist;
}

#define NC_FITS_ERROR(status) \
do { \
if (status) \
{ \
gchar errormsg[30]; \
fits_get_errstatus (status, errormsg); \
g_error ("FITS: %s", errormsg); \
} \
} while (FALSE);

#ifdef NUMCOSMO_HAVE_CFITSIO

/**
 * nc_cluster_abundance_catalog_save:
 * @data: a #NcData
 * @filename: name of the file
 * @overwrite: FIXME
 *
 * FIXME
 *
 */
void
nc_cluster_abundance_catalog_save (NcData *data, gchar *filename, gboolean overwrite)
{
  NcDataClusterAbundance *dca = (NcDataClusterAbundance *) NC_DATA_DATA (data);
  /*******************************************************************/
  /* Create a binary table extension                                 */
  /*******************************************************************/
  fitsfile *fptr;       /* pointer to the FITS file, defined in fitsio.h */
  gint status;
  gint tfields;

  gchar extname[] = "ClusterAbundance";   /* extension name */
  GPtrArray *ttype_array = g_ptr_array_sized_new (10);
  GPtrArray *tform_array = g_ptr_array_sized_new (10);
  GPtrArray *tunit_array = g_ptr_array_sized_new (10);

  guint z_obs_len = nc_cluster_redshift_obs_len (dca->z);
  guint z_obs_params_len = nc_cluster_redshift_obs_params_len (dca->z);
  guint lnM_obs_len = nc_cluster_mass_obs_len (dca->m);
  guint lnM_obs_params_len = nc_cluster_mass_obs_params_len (dca->m);

  g_ptr_array_set_free_func (tform_array, g_free);

  g_ptr_array_add (ttype_array, "Z_OBS");
  g_ptr_array_add (tform_array, g_strdup_printf ("%dD", z_obs_len));
  g_ptr_array_add (tunit_array, "REDSHIFT OBS");

  g_ptr_array_add (ttype_array, "LNM_OBS");
  g_ptr_array_add (tform_array, g_strdup_printf ("%dD", lnM_obs_len));
  g_ptr_array_add (tunit_array, "MASS OBS");

  if (dca->z_true != NULL)
  {
	g_ptr_array_add (ttype_array, "Z_TRUE");
	g_ptr_array_add (tform_array, g_strdup ("1D"));
	g_ptr_array_add (tunit_array, "TRUE REDSHIFT");
  }

  if (dca->lnM_true != NULL)
  {
	g_ptr_array_add (ttype_array, "LNM_TRUE");
	g_ptr_array_add (tform_array, g_strdup ("1D"));
	g_ptr_array_add (tunit_array, "TRUE LNM");
  }

  if (z_obs_params_len > 0)
  {
	g_ptr_array_add (ttype_array, "Z_OBS_PARAMS");
	g_ptr_array_add (tform_array, g_strdup_printf ("%dD", z_obs_params_len));
	g_ptr_array_add (tunit_array, "REDSHIFT OBS PARAMS");
  }

  if (lnM_obs_params_len > 0)
  {
	g_ptr_array_add (ttype_array, "LNM_OBS_PARAMS");
	g_ptr_array_add (tform_array, g_strdup_printf ("%dD", lnM_obs_params_len));
	g_ptr_array_add (tunit_array, "LNM OBS PARAMS");
  }

  tfields = ttype_array->len;

  /* initialize status before calling fitsio routines */
  status = 0;

  if (overwrite && g_file_test (filename, G_FILE_TEST_EXISTS))
	g_unlink (filename);

  /* create new FITS file */
  if (fits_create_file (&fptr, filename, &status))
	NC_FITS_ERROR (status);

  /* append a new empty binary table onto the FITS file */
  if (fits_create_tbl (fptr, BINARY_TBL, dca->np, tfields, (gchar **)ttype_array->pdata, (gchar **)tform_array->pdata,
                       (gchar **)tunit_array->pdata, extname, &status) )
	NC_FITS_ERROR (status);

  {
	gchar *z_ser = ncm_cfg_serialize_to_string (G_OBJECT (dca->z), FALSE);
	gchar *lnM_ser = ncm_cfg_serialize_to_string (G_OBJECT (dca->m), FALSE);

	if (fits_write_key_longstr (fptr, "Z_OBJ", z_ser, "Serialized redshift object", &status))
	  NC_FITS_ERROR(status);
	if (fits_write_key_longstr (fptr, "LNM_OBJ", lnM_ser, "Serialized mass object", &status))
	  NC_FITS_ERROR(status);

	g_free (z_ser);
	g_free (lnM_ser);
  }

  {
	gdouble sarea_d = dca->area_survey / gsl_pow_2 (M_PI / 180.0);
	if (fits_write_key(fptr, TDOUBLE, "AREA", &sarea_d, "Survey area in degree square", &status))
	  NC_FITS_ERROR(status);
  }

  {
	guint colnum = 1;
	if (fits_write_col (fptr, TDOUBLE, colnum, 1, 1, dca->np, ncm_matrix_ptr (dca->z_obs, 0, 0), &status))
	  NC_FITS_ERROR(status);

	colnum++;
	if (fits_write_col (fptr, TDOUBLE, colnum, 1, 1, dca->np, ncm_matrix_ptr (dca->lnM_obs, 0, 0), &status))
	  NC_FITS_ERROR(status);

	if (dca->z_true != NULL)
	{
	  colnum++;
	  if (fits_write_col (fptr, TDOUBLE, colnum, 1, 1, dca->np, ncm_vector_ptr (dca->z_true, 0), &status))
		NC_FITS_ERROR(status);
	}

	if (dca->lnM_true != NULL)
	{
	  colnum++;
	  if (fits_write_col (fptr, TDOUBLE, colnum, 1, 1, dca->np, ncm_vector_ptr (dca->lnM_true, 0), &status))
		NC_FITS_ERROR(status);
	}

	if (z_obs_params_len > 0)
	{
	  colnum++;
	  if (fits_write_col (fptr, TDOUBLE, colnum, 1, 1, dca->np, ncm_matrix_ptr (dca->z_obs_params, 0, 0), &status))
		NC_FITS_ERROR(status);
	}

	if (lnM_obs_params_len > 0)
	{
	  colnum++;
	  if (fits_write_col (fptr, TDOUBLE, colnum, 1, 1, dca->np, ncm_matrix_ptr (dca->lnM_obs_params, 0, 0), &status))
		NC_FITS_ERROR(status);
	}

  }

  if ( fits_close_file(fptr, &status) )
	NC_FITS_ERROR(status);

  g_ptr_array_unref (ttype_array);
  g_ptr_array_unref (tform_array);
  g_ptr_array_unref (tunit_array);

  return;
}

/**
 * nc_cluster_abundance_catalog_load:
 * @data: a #NcData.
 * @filename: name of the file
 *
 * FIXME
 *
 */
void
nc_cluster_abundance_catalog_load (NcData *data, gchar *filename)
{
  NcDataClusterAbundance *dca = (NcDataClusterAbundance *) NC_DATA_DATA (data);
  gint status, hdutype;
  gchar comment[FLEN_COMMENT];
  fitsfile *fptr;

  status = 0;

  if (filename == NULL)
	g_error ("nc_cluster_abundance_catalog_load: null filename");

  if (fits_open_file (&fptr, filename, READONLY, &status))
	NC_FITS_ERROR(status);

  if (fits_movabs_hdu (fptr, 2, &hdutype, &status))
	NC_FITS_ERROR (status);

  if (hdutype != BINARY_TBL)
	g_error ("%s (%d): Ncuster catalog is not binary!\n", __FILE__, __LINE__);

  if (dca->z != NULL)
	nc_cluster_redshift_free (dca->z);
  if (dca->m != NULL)
	nc_cluster_mass_free (dca->m);

  {
	gchar *z_ser = NULL;
	gchar *lnM_ser = NULL;

	if (fits_read_key_longstr (fptr, "Z_OBJ", &z_ser, comment, &status))
	  NC_FITS_ERROR (status);
	if (fits_read_key_longstr (fptr, "LNM_OBJ", &lnM_ser, comment, &status))
	  NC_FITS_ERROR (status);

    dca->z = nc_cluster_redshift_new_from_name (z_ser);
    dca->m = nc_cluster_mass_new_from_name (lnM_ser);

	g_free (z_ser);
	g_free (lnM_ser);
  }

  {
	glong nrows;
	if (fits_get_num_rows (fptr, &nrows, &status))
	  NC_FITS_ERROR(status);
	dca->np = nrows;
  }

  if (fits_read_key_dbl (fptr, "AREA", &dca->area_survey, comment, &status))
	g_error ("Fits file does not contain AREA in the header indicating the survey area (degree square). Use [col #AREA=...] to add this information.");
  dca->area_survey *= gsl_pow_2 (M_PI / 180.0);

  {
	gint z_obs_i, lnM_obs_i;
	gint z_obs_tc, lnM_obs_tc;
	glong z_obs_rp, lnM_obs_rp;
	glong z_obs_w, lnM_obs_w;

	if (fits_get_colnum (fptr, CASESEN, "Z_OBS", &z_obs_i, &status))
	  g_error ("Column Z_OBS not found, invalid fits file.");

	if (fits_get_colnum (fptr, CASESEN, "LNM_OBS", &lnM_obs_i, &status))
	  g_error ("Column LNM_OBS not found, invalid fits file.");

	if (fits_get_coltype (fptr, z_obs_i, &z_obs_tc, &z_obs_rp, &z_obs_w, &status))
	  g_error ("Column Z_OBS info not found, invalid fits file.");

	if (fits_get_coltype (fptr, lnM_obs_i, &lnM_obs_tc, &lnM_obs_rp, &lnM_obs_w, &status))
	  g_error ("Column LNM_OBS info not found, invalid fits file.");

	if (nc_cluster_redshift_obs_len (dca->z) != z_obs_rp)
	  g_error ("NcClusterRedshift object has observables length %d but fits has %ld.",
	           nc_cluster_redshift_obs_len (dca->z), z_obs_rp);

	if (nc_cluster_mass_obs_len (dca->m) != lnM_obs_rp)
	  g_error ("NcClusterMass object has observables length %d but fits has %ld.",
	           nc_cluster_mass_obs_len (dca->m), lnM_obs_rp);

	if (dca->z_obs)
	  ncm_matrix_free (dca->z_obs);
    dca->z_obs = ncm_matrix_new (dca->np, z_obs_rp);

	if (dca->lnM_obs)
	  ncm_matrix_free (dca->lnM_obs);
	dca->lnM_obs = ncm_matrix_new (dca->np, lnM_obs_rp);

	if (fits_read_col (fptr, TDOUBLE, z_obs_i, 1, 1, dca->np, NULL, ncm_matrix_ptr (dca->z_obs, 0, 0), NULL, &status))
	  NC_FITS_ERROR(status);

	if (fits_read_col (fptr, TDOUBLE, lnM_obs_i, 1, 1, dca->np, NULL, ncm_matrix_ptr (dca->lnM_obs, 0, 0), NULL, &status))
	  NC_FITS_ERROR(status);
  }

  {
	gint z_obs_params_i, lnM_obs_params_i;
	gint z_obs_params_tc, lnM_obs_params_tc;
	glong z_obs_params_rp, lnM_obs_params_rp;
	glong z_obs_params_w, lnM_obs_params_w;

	if (fits_get_colnum (fptr, CASESEN, "Z_OBS_PARAMS", &z_obs_params_i, &status))
	{
      if (nc_cluster_redshift_obs_params_len (dca->z) > 0)
		g_error ("NcClusterRedshift object has observable parameters length %d but fits has 0.",
		         nc_cluster_redshift_obs_params_len (dca->z));
	}
	else
	{
	  if (fits_get_coltype (fptr, z_obs_params_i, &z_obs_params_tc, &z_obs_params_rp, &z_obs_params_w, &status))
		g_error ("Column Z_OBS_PARAMS info not found, invalid fits file.");

	  if (nc_cluster_redshift_obs_params_len (dca->z) != z_obs_params_rp)
		g_error ("NcClusterRedshift object has observable parameters length %d but fits has %ld.",
		         nc_cluster_redshift_obs_params_len (dca->z), z_obs_params_rp);

	  if (dca->z_obs_params)
		ncm_matrix_free (dca->z_obs_params);
	  dca->z_obs_params = ncm_matrix_new (dca->np, z_obs_params_rp);

	  if (fits_read_col (fptr, TDOUBLE, z_obs_params_i, 1, 1, dca->np, NULL, ncm_matrix_ptr (dca->z_obs_params, 0, 0), NULL, &status))
		NC_FITS_ERROR(status);
	}

	if (fits_get_colnum (fptr, CASESEN, "LNM_OBS_PARAMS", &lnM_obs_params_i, &status))
	{
      if (nc_cluster_mass_obs_params_len (dca->m) > 0)
		g_error ("NcClusterMass object has observable parameters length %d but fits has 0.",
		         nc_cluster_mass_obs_params_len (dca->m));
	}
	else
	{
	  if (fits_get_coltype (fptr, lnM_obs_params_i, &lnM_obs_params_tc, &lnM_obs_params_rp, &lnM_obs_params_w, &status))
		g_error ("Column LNM_OBS_PARAMS info not found, invalid fits file.");

	  if (nc_cluster_mass_obs_params_len (dca->m) != lnM_obs_params_rp)
		g_error ("NcClusterMass object has observable parameters length %d but fits has %ld.",
		         nc_cluster_mass_obs_params_len (dca->m), lnM_obs_params_rp);

	  if (dca->lnM_obs_params)
		ncm_matrix_free (dca->lnM_obs_params);
	  dca->lnM_obs_params = ncm_matrix_new (dca->np, lnM_obs_params_rp);


	  if (fits_read_col (fptr, TDOUBLE, lnM_obs_params_i, 1, 1, dca->np, NULL, ncm_matrix_ptr (dca->lnM_obs_params, 0, 0), NULL, &status))
		NC_FITS_ERROR(status);
	}
  }

  {
	gint z_true_i, lnM_true_i;
	if (dca->z_true != NULL)
	{
	  ncm_vector_free (dca->z_true);
	  dca->z_true = NULL;
	}
	if (!fits_get_colnum (fptr, CASESEN, "Z_TRUE", &z_true_i, &status))
	{
	  dca->z_true = ncm_vector_new (dca->np);
	  if (fits_read_col (fptr, TDOUBLE, z_true_i, 1, 1, dca->np, NULL, ncm_vector_ptr (dca->z_true, 0), NULL, &status))
		NC_FITS_ERROR(status);
	}

	if (dca->lnM_true != NULL)
	{
	  ncm_vector_free (dca->lnM_true);
	  dca->lnM_true = NULL;
	}
	if (!fits_get_colnum (fptr, CASESEN, "LNM_TRUE", &lnM_true_i, &status))
	{
	  dca->lnM_true = ncm_vector_new (dca->np);
	  if (fits_read_col (fptr, TDOUBLE, lnM_true_i, 1, 1, dca->np, NULL, ncm_vector_ptr (dca->lnM_true, 0), NULL, &status))
		NC_FITS_ERROR(status);

	}
  }

  if ( fits_close_file(fptr, &status) )
	NC_FITS_ERROR(status);

  nc_data_init (data);

  return;
}
#endif /* NUMCOSMO_HAVE_CFITSIO */

/***************************************************************************
 *            data_cluster_abundance.c
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
 * SECTION:data_cluster_abundance
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
#include <fitsio.h>

static void
_ca_data_copy (gpointer dest_ptr, gpointer src_ptr)
{
  NcDataClusterAbundance *dest = (NcDataClusterAbundance *) dest_ptr;
  NcDataClusterAbundance *src = (NcDataClusterAbundance *) src_ptr;

  dest->area_survey        = src->area_survey;
  dest->real.lnMi          = src->real.lnMi;
  dest->real.lnMf          = src->real.lnMf;
  dest->real.lnM_sigma0    = src->real.lnM_sigma0;
  dest->real.lnM_bias      = src->real.lnM_bias;
  dest->obs.lnMi           = src->obs.lnMi;
  dest->obs.lnMf           = src->obs.lnMf;
  dest->obs.lnM_sigma0     = src->obs.lnM_sigma0;
  dest->obs.lnM_bias       = src->obs.lnM_bias;
  dest->np                 = src->np;
  dest->real.zi            = src->real.zi;
  dest->real.zf            = src->real.zf;
  dest->real.photoz_sigma0 = src->real.photoz_sigma0;
  dest->real.photoz_bias   = src->real.photoz_bias;
  dest->obs.zi             = src->obs.zi;
  dest->obs.zf             = src->obs.zf;
  dest->obs.photoz_sigma0  = src->obs.photoz_sigma0;
  dest->obs.photoz_bias    = src->obs.photoz_bias;

  dest->opt                = src->opt;

  if (src->real.z_lnM != NULL)
  {
	if (dest->real.z_lnM == NULL)
	  dest->real.z_lnM = gsl_matrix_alloc (src->real.z_lnM->size1, src->real.z_lnM->size2);
	if ((dest->real.z_lnM->size1 != src->real.z_lnM->size1) || (dest->real.z_lnM->size2 != src->real.z_lnM->size2))
	{
	  gsl_matrix_free (dest->real.z_lnM);
	  dest->real.z_lnM = gsl_matrix_alloc (src->real.z_lnM->size1, src->real.z_lnM->size2);
	}
	gsl_matrix_memcpy (dest->real.z_lnM, src->real.z_lnM);
  }

  if (src->obs.z_lnM != NULL)
  {
	if (dest->obs.z_lnM == NULL)
	  dest->obs.z_lnM = gsl_matrix_alloc (src->obs.z_lnM->size1, src->obs.z_lnM->size2);
	if ((dest->obs.z_lnM->size1 != src->obs.z_lnM->size1) || (dest->obs.z_lnM->size2 != src->obs.z_lnM->size2))
	{
	  gsl_matrix_free (dest->obs.z_lnM);
	  dest->obs.z_lnM = gsl_matrix_alloc (src->obs.z_lnM->size1, src->obs.z_lnM->size2);
	}
	gsl_matrix_memcpy (dest->obs.z_lnM, src->obs.z_lnM);
  }
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
  if (ca->real.z_lnM != NULL)
	gsl_matrix_free (ca->real.z_lnM);
  if (ca->obs.z_lnM != NULL)
	gsl_matrix_free (ca->obs.z_lnM);
  g_slice_free (NcDataClusterAbundance, ca_ptr);
}

static void
_ca_data_begin (gpointer ca_ptr)
{
  NcDataClusterAbundance *ca = (NcDataClusterAbundance *) ca_ptr;
  ca->log_np_fac = lgamma (ca->np + 1);
}

static guint _ca_data_get_length (gpointer ca_ptr) { return ((NcDataClusterAbundance *) ca_ptr)->np + 1; }

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

#define NC_DATA_CLUSTER_ABUNDANCE_GROUP "Cluster Abundance Data"
#define NC_DATA_CLUSTER_ABUNDANCE_NKEYS 5

static gchar *NC_DATA_CLUSTER_ABUNDANCE_KEYS[] = {"SurveyArea", "lnMi", "lnMf", "RedshiftNodes", "NumberCounts"};

/**
 * nc_data_cluster_abundance_binned_init_from_text_file_gkey:
 * @data: a #NcData
 * @obs: TRUE if the histogram provides observational data, FALSE if it provides real values of z and mass.
 * @filename: name of the file
 *
 * FIXME
 *
 */
void
nc_data_cluster_abundance_binned_init_from_text_file_gkey (NcData *data, gboolean obs, gchar *filename)
{
  GKeyFile *kf = g_key_file_new ();
  GError *err = NULL;
  gint i;

  g_key_file_load_from_file (kf, filename, G_KEY_FILE_NONE, &err);

  if (!g_key_file_has_group (kf, NC_DATA_CLUSTER_ABUNDANCE_GROUP))
	g_error ("File (%s) do not contain %s group.", filename, NC_DATA_CLUSTER_ABUNDANCE_GROUP);
  for (i = 0; i < NC_DATA_CLUSTER_ABUNDANCE_NKEYS; i++)
  {
	if (!g_key_file_has_key (kf, NC_DATA_CLUSTER_ABUNDANCE_GROUP, NC_DATA_CLUSTER_ABUNDANCE_KEYS[i], &err))
	  g_error ("File (%s) do not contain the key %s.", filename, NC_DATA_CLUSTER_ABUNDANCE_KEYS[i]);
  }

  {
	NcDataStruct *dts_ca = _nc_data_struct_cluster_abundance_new ();
	NcDataClusterAbundance *ca = NC_DATA_STRUCT_DATA (dts_ca);
	gsize nodes_length;
	gsize N_length;
	gdouble *nodes = g_key_file_get_double_list (kf, NC_DATA_CLUSTER_ABUNDANCE_GROUP, NC_DATA_CLUSTER_ABUNDANCE_KEYS[3], &nodes_length, &err);
	gdouble *bin = g_key_file_get_double_list (kf, NC_DATA_CLUSTER_ABUNDANCE_GROUP, NC_DATA_CLUSTER_ABUNDANCE_KEYS[4], &N_length, &err);
	if (N_length + 1 != nodes_length)
	  g_error ("Incompatible number of nodes and number counts [%zd %zd]\n", nodes_length, N_length);

	ca->area_survey = g_key_file_get_double (kf, NC_DATA_CLUSTER_ABUNDANCE_GROUP, NC_DATA_CLUSTER_ABUNDANCE_KEYS[0], &err);
	if (obs)
	{
	  ca->obs.lnMi = g_key_file_get_double (kf, NC_DATA_CLUSTER_ABUNDANCE_GROUP, NC_DATA_CLUSTER_ABUNDANCE_KEYS[1], &err);
	  ca->obs.lnMi = g_key_file_get_double (kf, NC_DATA_CLUSTER_ABUNDANCE_GROUP, NC_DATA_CLUSTER_ABUNDANCE_KEYS[2], &err);
	}
	else
	{
	  ca->real.lnMi = g_key_file_get_double (kf, NC_DATA_CLUSTER_ABUNDANCE_GROUP, NC_DATA_CLUSTER_ABUNDANCE_KEYS[1], &err);
	  ca->real.lnMi = g_key_file_get_double (kf, NC_DATA_CLUSTER_ABUNDANCE_GROUP, NC_DATA_CLUSTER_ABUNDANCE_KEYS[2], &err);
	}

	{
	  gsl_histogram *h = gsl_histogram_alloc (N_length);

	  for (i = 0; i < N_length; i++)
	  {
		h->range[i] = nodes[i];
		h->bin[i] = bin[i];
	  }
	  h->range[i] = nodes[i];

	  g_free (nodes);
	  g_free (bin);

	  if (obs)
	  {
		ca->obs.zi = h->range[0];
		ca->obs.zf = h->range[h->n];
	  }
	  else
	  {
		ca->real.zi = h->range[0];
		ca->real.zf = h->range[h->n];
	  }

	  nc_data_poisson_init_from_histogram (data, h, TRUE, dts_ca);
	}
  }

  g_key_file_free (kf);
}

/**
 * nc_data_cluster_abundance_binned_init_from_sampling:
 * @data: a #NcData.
 * @mset: a #NcmMSet.
 * @nodes: a #NcmVector.
 * @opt: a #NcClusterAbundanceOpt.
 * @obs: TRUE if the sample observational data, FALSE if it provides real values of z and mass.
 * @area_survey: area in units of square degrees.
 * @lnMi: logarithm base e of the minimum mass.
 * @lnMf: logarithm base e of the maximum mass.
 * @photoz_sigma0: FIXME
 * @photoz_bias: FIXME
 * @lnM_sigma0: FIXME
 * @lnM_bias: FIXME
 *
 * FIXME
 */
void
nc_data_cluster_abundance_binned_init_from_sampling (NcData *data, NcmMSet *mset, NcmVector *nodes, NcClusterAbundanceOpt opt, gboolean obs, gdouble area_survey, gdouble lnMi, gdouble lnMf, gdouble photoz_sigma0, gdouble photoz_bias, gdouble lnM_sigma0, gdouble lnM_bias)
{
  NcDataStruct *dts_ca = _nc_data_struct_cluster_abundance_new ();
  NcDataClusterAbundance *ca = NC_DATA_STRUCT_DATA (dts_ca);

  ca->area_survey = area_survey;
  ca->opt = opt;
  if (obs)
  {
	ca->obs.lnMi = lnMi;
	ca->obs.lnMf = lnMf;
	ca->obs.zi = ncm_vector_get (nodes, 0);
	ca->obs.zf = ncm_vector_get (nodes, ncm_vector_len (nodes) - 1);
	ca->obs.photoz_sigma0 = photoz_sigma0;
	ca->obs.photoz_bias = photoz_bias;
	ca->obs.lnM_sigma0 = lnM_sigma0;
	ca->obs.lnM_bias = lnM_bias;
  }
  else
  {
	ca->real.lnMi = lnMi;
	ca->real.lnMf = lnMf;
	ca->real.zi = ncm_vector_get (nodes, 0);
	ca->real.zf = ncm_vector_get (nodes, ncm_vector_len (nodes) - 1);
	ca->real.photoz_sigma0 = photoz_sigma0;
	ca->real.photoz_bias = photoz_bias;
	ca->real.lnM_sigma0 = lnM_sigma0;
	ca->real.lnM_bias = lnM_bias;
  }

  nc_data_poisson_init_zero (data, nodes, dts_ca);
  nc_data_resample (data, mset, FALSE);
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

/**
 * nc_data_cluster_abundance_binned_save:
 * @data: a #NcData.
 * @filename: name of the file.
 *
 * FIXME
 */
void
nc_data_cluster_abundance_binned_save (NcData *data, gchar *filename)
{
  GKeyFile *kf = g_key_file_new ();
  GError *err = NULL;
  NcDataPoisson *poisson = NC_DATA_STRUCT_DATA (data->dts);
  NcDataClusterAbundance *ca = NC_DATA_STRUCT_DATA (poisson->extra_data);
  gchar *file_data;
  gsize data_size;

  g_assert (data->init);

  g_key_file_set_double (kf, NC_DATA_CLUSTER_ABUNDANCE_GROUP, NC_DATA_CLUSTER_ABUNDANCE_KEYS[0], ca->area_survey);

  g_key_file_set_double (kf, NC_DATA_CLUSTER_ABUNDANCE_GROUP, NC_DATA_CLUSTER_ABUNDANCE_KEYS[1], ca->obs.lnMi);
  g_key_file_set_double (kf, NC_DATA_CLUSTER_ABUNDANCE_GROUP, NC_DATA_CLUSTER_ABUNDANCE_KEYS[2], ca->obs.lnMf);

  g_key_file_set_double (kf, NC_DATA_CLUSTER_ABUNDANCE_GROUP, NC_DATA_CLUSTER_ABUNDANCE_KEYS[1], ca->real.lnMi);
  g_key_file_set_double (kf, NC_DATA_CLUSTER_ABUNDANCE_GROUP, NC_DATA_CLUSTER_ABUNDANCE_KEYS[2], ca->real.lnMf);

  g_key_file_set_double_list (kf, NC_DATA_CLUSTER_ABUNDANCE_GROUP, NC_DATA_CLUSTER_ABUNDANCE_KEYS[3], poisson->h->range, poisson->h->n + 1);
  g_key_file_set_double_list (kf, NC_DATA_CLUSTER_ABUNDANCE_GROUP, NC_DATA_CLUSTER_ABUNDANCE_KEYS[4], poisson->h->bin, poisson->h->n);

  file_data = g_key_file_to_data (kf, &data_size, &err);
  ncm_cfg_init ();

  {
	gchar *full_filename = ncm_cfg_get_fullpath ("%s", filename);
	g_file_set_contents (full_filename, file_data, data_size, &err);
	g_free (full_filename);
  }
}

static void
_nc_data_cluster_abundance_binned_f (NcmMSet *mset, gpointer obj, const gdouble *x, gdouble *f)
{
  NcHICosmo *model = NC_HICOSMO (ncm_mset_peek (mset, NC_HICOSMO_ID));
  NcClusterAbundance *cad = NC_CLUSTER_ABUNDANCE (obj);
  f[0] = nc_cluster_abundance_N_val (cad, model, cad->lnMi, cad->lnMf, cad->zi, x[0]);
  //printf ("% 20.15g % 20.15g % 20.15g\n", cad->zi, z, mi);
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

/**
 * nc_data_cluster_abundance_binned_lnM_z_init_from_hist: (skip)
 * @data: a #NcData.
 * @obs: TRUE if the histogram provides observational data, FALSE if it provides real values of z and mass.
 * @hist: FIXME
 * @opt: a #NcClusterAbundanceOpt.
 * @area_survey: area in units of square degrees.
 * @photoz_sigma0: FIXME
 * @photoz_bias: FIXME
 * @lnM_sigma0: FIXME
 * @lnM_bias: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
void
nc_data_cluster_abundance_binned_lnM_z_init_from_hist (NcData *data, gboolean obs, gsl_histogram2d *hist, NcClusterAbundanceOpt opt, gdouble area_survey, gdouble photoz_sigma0, gdouble photoz_bias, gdouble lnM_sigma0, gdouble lnM_bias)
{
  NcDataStruct *dts_ca = _nc_data_struct_cluster_abundance_new ();
  NcDataClusterAbundance *ca = NC_DATA_STRUCT_DATA (dts_ca);
  ca->area_survey = area_survey;
  if (obs)
  {
	ca->obs.lnMi = hist->xrange[0];
	ca->obs.lnMf = hist->xrange[hist->nx];
	ca->obs.zi = hist->yrange[0];
	ca->obs.zf = hist->yrange[hist->ny];
	ca->obs.photoz_sigma0 = photoz_sigma0;
	ca->obs.photoz_bias = photoz_bias;
	ca->obs.lnM_sigma0 = lnM_sigma0;
	ca->obs.lnM_bias = lnM_bias;
  }
  else
  {
	ca->real.lnMi = hist->xrange[0];
	ca->real.lnMf = hist->xrange[hist->nx];
	ca->real.zi = hist->yrange[0];
	ca->real.zf = hist->yrange[hist->ny];
	ca->real.photoz_sigma0 = photoz_sigma0;
	ca->real.photoz_bias = photoz_bias;
	ca->real.lnM_sigma0 = lnM_sigma0;
	ca->real.lnM_bias = lnM_bias;
  }
  ca->opt = opt;

  //nc_data_poisson_init_from_histogram (data, hist, FALSE, dts_ca);
}

/************************************************************************************************************
 * Unbinned abundance data                                                                                  *
 ************************************************************************************************************/

static void
_nc_data_cluster_abundance_unbinned_model_init (gpointer model, gpointer data)
{
  NcClusterAbundance *cad = NC_CLUSTER_ABUNDANCE (model);
  NcDataClusterAbundance *dca = (NcDataClusterAbundance *) data;

  cad->zi = dca->real.zi;
  cad->zf = dca->real.zf;
  cad->lnMi = dca->real.lnMi;
  cad->lnMf = dca->real.lnMf;
//  cad->photoz_sigma0 = dca->real.photoz_sigma0;
  g_assert_not_reached ();
  cad->lnM_sigma0    = dca->real.lnM_sigma0;
  printf ("Real: % 15.5g % 15.5g\n", cad->zi, cad->zf);

  cad->zi = dca->obs.zi;
  cad->zf = dca->obs.zf;
  cad->lnMi = dca->obs.lnMi;
  cad->lnMf = dca->obs.lnMf;
  // cad->photoz_sigma0 = dca->obs.photoz_sigma0;
  g_assert_not_reached ();
  cad->lnM_sigma0    = dca->obs.lnM_sigma0;
  printf ("Obs: % 15.5g % 15.5g\n", cad->zi, cad->zf);

  cad->completeness  = dca->completeness;
  cad->purity        = dca->purity;
  cad->sd_lnM        = dca->sd_lnM;

  nc_cluster_abundance_set_options (cad, dca->opt);

  cad->mfp->area_survey = dca->area_survey;
}

static void
_nc_data_cluster_abundance_unbinned_prepare (NcmMSet *mset, gpointer model, gpointer data)
{
  NcClusterAbundance *cad = NC_CLUSTER_ABUNDANCE (model);
  /* NcDataNClusterAbundance *dca = (NcDataNClusterAbundance *) data; FIXME */

  if (ncm_model_ctrl_update (cad->ctrl, ncm_mset_peek (mset, NC_HICOSMO_ID)))
  {
	nc_cluster_abundance_prepare (cad, model);
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
  gsl_rng *rng = ncm_get_rng ();
  guint np;
  gint i;

  const gint test_mask = NC_CLUSTER_ABUNDANCE_PHOTOZ |
	NC_CLUSTER_ABUNDANCE_MOBS |
	NC_CLUSTER_ABUNDANCE_REAL_ZM |
	NC_CLUSTER_ABUNDANCE_OBS_ZM;
  const gint opt_mask = cad->opt & test_mask;

  np = gsl_ran_poisson (rng, cad->norma);

  printf ("Dentro do resample\n");
  switch (opt_mask)
  {
	case NC_CLUSTER_ABUNDANCE_REAL_ZM:
	{
	  if (dca->real.z_lnM == NULL)
		dca->real.z_lnM = gsl_matrix_alloc (np, 2);
	  else if (dca->np != np)
	  {
		gsl_matrix_free (dca->real.z_lnM);
		dca->real.z_lnM = gsl_matrix_alloc (np, 2);
	  }
	}
	  break;
	case NC_CLUSTER_ABUNDANCE_OBS_ZM:
	{
	  if (dca->obs.z_lnM == NULL)
		dca->obs.z_lnM = gsl_matrix_alloc (np, 2);
	  else if (dca->np != np)
	  {
		gsl_matrix_free (dca->obs.z_lnM);
		dca->obs.z_lnM = gsl_matrix_alloc (np, 2);
	  }
	}
	  break;
  }

  dca->np = np;

  nc_cluster_abundance_prepare_inv_dNdz_no_obs (cad, NC_HICOSMO (ncm_mset_peek (mset, NC_HICOSMO_ID)));
  //nc_cluster_abundance_prepare_inv_dNdz (cad, NC_HICOSMO(model));

  //printf ("# Generating unbinned %ld (z,lnM) %g\n", dca->np, cad->norma);

  switch (opt_mask)
  {
	case NC_CLUSTER_ABUNDANCE_PHOTOZ:
	{
	  printf ("Chegou para gerar!\n");
	  for (i = 0; i < dca->np; i++)
	  {
		gdouble zi_obs;
		const gdouble u1 = _nc_cad_inv_dNdz_convergence_f (gsl_rng_uniform_pos (rng), cad->z_epsilon);
		const gdouble u2 = _nc_cad_inv_dNdz_convergence_f (gsl_rng_uniform_pos (rng), cad->lnM_epsilon);
		const gdouble zi_real = ncm_spline_eval (cad->inv_z, u1);
		const gdouble lnMi_real = ncm_spline2d_eval (cad->inv_lnM_z, u2, zi_real);
		//const gdouble sigma_z = cad->photoz_sigma0 * (1.0 + zi_real);
		g_assert_not_reached ();

		do {
		  //zi_obs = zi_real + gsl_ran_gaussian (rng, sigma_z);
		  g_assert_not_reached ();
		} while (zi_obs < 0);

		gsl_matrix_set (dca->real.z_lnM, i, 0, zi_real);
		gsl_matrix_set (dca->real.z_lnM, i, 1, lnMi_real);
		gsl_matrix_set (dca->obs.z_lnM, i, 0, zi_obs);
		gsl_matrix_set (dca->obs.z_lnM, i, 1, lnMi_real);
		printf ("% 20.15g % 20.15g % 20.15g\n", zi_real, lnMi_real, zi_obs);
	  }
	}
	  break;
	case NC_CLUSTER_ABUNDANCE_MOBS:
	{
	  for (i = 0; i < dca->np; i++)
	  {
		gdouble lnMi_obs;
		const gdouble u1 = _nc_cad_inv_dNdz_convergence_f (gsl_rng_uniform_pos (rng), cad->z_epsilon);
		const gdouble u2 = _nc_cad_inv_dNdz_convergence_f (gsl_rng_uniform_pos (rng), cad->lnM_epsilon);
		const gdouble zi_real = ncm_spline_eval (cad->inv_z, u1);
		const gdouble lnMi_real = ncm_spline2d_eval (cad->inv_lnM_z, u2, zi_real);
		const gdouble sigma_M = cad->lnM_sigma0;
		//printf ("% 20.15g % 20.15g % 20.15g % 20.15g\n", u1, u2, zi, lnMi);

		lnMi_obs = lnMi_real + gsl_ran_gaussian (rng, sigma_M);

		gsl_matrix_set (dca->real.z_lnM, i, 0, zi_real);
		gsl_matrix_set (dca->real.z_lnM, i, 1, lnMi_real);
		gsl_matrix_set (dca->obs.z_lnM, i, 0, zi_real);
		gsl_matrix_set (dca->obs.z_lnM, i, 1, lnMi_obs);
	  }
	}
	  break;
	case NC_CLUSTER_ABUNDANCE_PHOTOZ | NC_CLUSTER_ABUNDANCE_MOBS:
	{
	  for (i = 0; i < dca->np; i++)
	  {
		gdouble zi_obs, lnMi_obs;
		const gdouble u1 = _nc_cad_inv_dNdz_convergence_f (gsl_rng_uniform_pos (rng), cad->z_epsilon);
		const gdouble u2 = _nc_cad_inv_dNdz_convergence_f (gsl_rng_uniform_pos (rng), cad->lnM_epsilon);
		const gdouble zi_real = ncm_spline_eval (cad->inv_z, u1);
		const gdouble lnMi_real = ncm_spline2d_eval (cad->inv_lnM_z, u2, zi_real);
		//const gdouble sigma_z = cad->photoz_sigma0 * (1.0 + zi_real);
		g_assert_not_reached ();
		const gdouble sigma_M = cad->lnM_sigma0;
		//printf ("% 20.15g % 20.15g % 20.15g % 20.15g\n", u1, u2, zi, lnMi);

		do {
		  //zi_obs = zi_real + gsl_ran_gaussian (rng, sigma_z);
		  g_assert_not_reached ();
		} while (zi_obs < 0);

		lnMi_obs = lnMi_real + gsl_ran_gaussian (rng, sigma_M);

		gsl_matrix_set (dca->real.z_lnM, i, 0, zi_real);
		gsl_matrix_set (dca->real.z_lnM, i, 1, lnMi_real);
		gsl_matrix_set (dca->obs.z_lnM, i, 0, zi_obs);
		gsl_matrix_set (dca->obs.z_lnM, i, 1, lnMi_obs);
	  }
	}
	  break;
	default:
	{
	  for (i = 0; i < dca->np; i++)
	  {
		const gdouble u1 = _nc_cad_inv_dNdz_convergence_f (gsl_rng_uniform_pos (rng), cad->z_epsilon);//gsl_rng_uniform (rng);
		const gdouble u2 = _nc_cad_inv_dNdz_convergence_f (gsl_rng_uniform_pos (rng), cad->lnM_epsilon);//gsl_rng_uniform (rng);
		const gdouble zi_real = ncm_spline_eval (cad->inv_z, u1);
		const gdouble lnMi_real = ncm_spline2d_eval (cad->inv_lnM_z, u2, zi_real);

		//printf ("% 20.15g % 20.15g % 20.15g % 20.15g\n", u1, u2, zi, lnMi);

		gsl_matrix_set (dca->real.z_lnM, i, 0, zi_real);
		gsl_matrix_set (dca->real.z_lnM, i, 1, lnMi_real);
	  }
	}
  }
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

static void
_nc_data_cluster_abundance_calc_m2lnL (NcmMSet *mset, gpointer model, gpointer data, gdouble *m2lnL)
{
  NcClusterAbundance *cad = NC_CLUSTER_ABUNDANCE (model);
  NcDataClusterAbundance *dca = (NcDataClusterAbundance *) data;
  NcHICosmo *m = NC_HICOSMO (ncm_mset_peek (mset, NC_HICOSMO_ID));
  GTimer *bench = g_timer_new ();
  gint i;

  const gint test_mask = NC_CLUSTER_ABUNDANCE_PHOTOZ |
	NC_CLUSTER_ABUNDANCE_MOBS;
  const gint opt_mask = cad->opt & test_mask;

  *m2lnL = 0.0;

  switch (opt_mask)
  {
	case NC_CLUSTER_ABUNDANCE_PHOTOZ:
	{
	  if (dca->opt & NC_CLUSTER_ABUNDANCE_BINMASS)
	  {
		for (i = 0; i < dca->np; i++)
		{
		  const gdouble zi = gsl_matrix_get (dca->obs.z_lnM, i, 0);
		  *m2lnL -= log (nc_cluster_abundance_dNdz_val (cad, m, cad->lnMi, cad->lnMf, zi));
		}
	  }
	  else
	  {
		for (i = 0; i < dca->np; i++)
		{
		  const gdouble zi = gsl_matrix_get (dca->obs.z_lnM, i, 0);
		  const gdouble lnMi = gsl_matrix_get (dca->real.z_lnM, i, 1);
		  const gdouble mlnL_i = -log (nc_cluster_abundance_d2NdzdlnM_val (cad, m, lnMi, zi));
		  *m2lnL += mlnL_i;
		}
	  }
	}
	case NC_CLUSTER_ABUNDANCE_MOBS:
	{
	  if (dca->opt & NC_CLUSTER_ABUNDANCE_BINMASS)
	  {
		for (i = 0; i < dca->np; i++)
		{
		  const gdouble zi = gsl_matrix_get (dca->real.z_lnM, i, 0);
		  *m2lnL -= log (nc_cluster_abundance_dNdz_val (cad, m, cad->lnMi, cad->lnMf, zi));
		}
	  }
	  else
	  {
		for (i = 0; i < dca->np; i++)
		{
		  const gdouble zi = gsl_matrix_get (dca->real.z_lnM, i, 0);
		  const gdouble lnMi = gsl_matrix_get (dca->obs.z_lnM, i, 1);
		  const gdouble mlnL_i = -log (nc_cluster_abundance_d2NdzdlnM_val (cad, m, lnMi, zi));
		  *m2lnL += mlnL_i;
		}
	  }
	}
	case NC_CLUSTER_ABUNDANCE_PHOTOZ | NC_CLUSTER_ABUNDANCE_MOBS:
	{
	  if (dca->opt & NC_CLUSTER_ABUNDANCE_BINMASS)
	  {
		for (i = 0; i < dca->np; i++)
		{
		  const gdouble zi = gsl_matrix_get (dca->obs.z_lnM, i, 0);
		  *m2lnL -= log (nc_cluster_abundance_dNdz_val (cad, m, cad->lnMi, cad->lnMf, zi));
		}
	  }
	  else
	  {
		for (i = 0; i < dca->np; i++)
		{
		  const gdouble zi = gsl_matrix_get (dca->obs.z_lnM, i, 0);
		  const gdouble lnMi = gsl_matrix_get (dca->obs.z_lnM, i, 1);
		  const gdouble mlnL_i = -log (nc_cluster_abundance_d2NdzdlnM_val (cad, m, lnMi, zi));
		  *m2lnL += mlnL_i;
		}
	  }
	}
	case 0:
	{
	  if (dca->opt & NC_CLUSTER_ABUNDANCE_BINMASS)
	  {
		for (i = 0; i < dca->np; i++)
		{
		  const gdouble zi = gsl_matrix_get (dca->real.z_lnM, i, 0);
		  *m2lnL -= log (nc_cluster_abundance_dNdz_val (cad, m, cad->lnMi, cad->lnMf, zi));
		}
	  }
	  else
	  {
		for (i = 0; i < dca->np; i++)
		{
		  const gdouble zi = gsl_matrix_get (dca->real.z_lnM, i, 0);
		  const gdouble lnMi = gsl_matrix_get (dca->real.z_lnM, i, 1);
		  const gdouble mlnL_i = -log (nc_cluster_abundance_d2NdzdlnM_val (cad, m, lnMi, zi));
		  *m2lnL += mlnL_i;
		}
	  }
	}
  }

  //printf ("% 20.15g\n", 2.0 * ((*m2lnL) + dca->np * log (cad->norma)));
  *m2lnL += (dca->log_np_fac + cad->norma * cad->completeness_factor);

  *m2lnL *= 2.0;
  //printf ("# nhoc % 20.15g | took %fs\n", *m2lnL, g_timer_elapsed (bench, NULL));
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
  data->model_init            = &_nc_data_cluster_abundance_unbinned_model_init;
  data->model                 = cad;
  data->model_free            = NULL;
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
 * nc_data_cluster_abundance_unbinned_init_from_sampling:
 * @data: a #NcData.
 * @mset: a #NcmMSet.
 * @opt: a #NcClusterAbundanceOpt.
 * @area_survey: area in units of square degrees.
 * @lnMi: logarithm base e of the minimum mass.
 * @lnMf: logarithm base e of the maximum mass.
 * @z_initial: minimum redshift.
 * @z_final: maximum redshift.
 * @photoz_sigma0: FIXME
 * @photoz_bias: FIXME
 * @lnM_sigma0: FIXME
 * @lnM_bias: FIXME
 *
 * FIXME
 *
 */
void
nc_data_cluster_abundance_unbinned_init_from_sampling (NcData *data, NcmMSet *mset, NcClusterAbundanceOpt opt, gdouble area_survey, gdouble lnMi, gdouble lnMf, gdouble z_initial, gdouble z_final, gdouble photoz_sigma0, gdouble photoz_bias, gdouble lnM_sigma0, gdouble lnM_bias)
{
  NcDataClusterAbundance *dca = (NcDataClusterAbundance *) NC_DATA_DATA (data);

  if (z_initial == 0.0)
  {
	g_message ("nc_data_cluster_abundance_unbinned_init_from_sampling: User requested zi == 0.0. Using zi = %e\n", _NC_CLUSTER_ABUNDANCE_MIN_Z);
	z_initial = _NC_CLUSTER_ABUNDANCE_MIN_Z;
  }

  dca->area_survey = area_survey;
  dca->opt = opt;

  if (opt & NC_CLUSTER_ABUNDANCE_OBS_ZM)
  {
	dca->obs.lnMi = lnMi;
	dca->obs.lnMf = lnMf;
	dca->obs.zi = z_initial;
	dca->obs.zf = z_final;

	dca->obs.photoz_sigma0 = photoz_sigma0;
	dca->obs.photoz_bias = photoz_bias;
	dca->obs.lnM_sigma0 = lnM_sigma0;
	dca->obs.lnM_bias = lnM_bias;
  }
  else
  {
	printf ("init_from_sampling\n");
	dca->real.lnMi = lnMi;
	dca->real.lnMf = lnMf;
	dca->real.zi = z_initial;
	dca->real.zf = z_final;

	dca->real.photoz_sigma0 = photoz_sigma0;
	dca->real.photoz_bias = photoz_bias;
	dca->real.lnM_sigma0 = lnM_sigma0;
	dca->real.lnM_bias = lnM_bias;
	printf ("z bias = % 15.5g lnM bias = % 15.5g\n", dca->real.photoz_bias, dca->real.lnM_bias);
  }

  nc_data_model_init (data);
  nc_data_resample (data, mset, FALSE);

  nc_data_init (data);
}

/**
 * nc_data_cluster_abundance_unbinned_init_from_text_file:
 * @data: a #NcData
 * @filename: name of the file
 * @opt: a #NcClusterAbundanceOpt
 * @area_survey: area in units of square degrees
 * @lnMi: logarithm base e of the minimum mass
 * @lnMf: logarithm base e of the maximum mass
 * @z_initial: minimum redshift
 * @z_final: maximum redshift
 * @photoz_sigma0: FIXME
 * @photoz_bias: FIXME
 * @lnM_sigma0: FIXME
 * @lnM_bias: FIXME
 *
 * FIXME
 */
void
nc_data_cluster_abundance_unbinned_init_from_text_file (NcData *data, gchar *filename, NcClusterAbundanceOpt opt, gdouble area_survey, gdouble lnMi, gdouble lnMf, gdouble z_initial, gdouble z_final, gdouble photoz_sigma0, gdouble photoz_bias, gdouble lnM_sigma0, gdouble lnM_bias)
{
  NcDataClusterAbundance *dca = (NcDataClusterAbundance *) NC_DATA_DATA (data);
  FILE *catalog;
  gint c, nlines = 0, i = 0, j = 0;
  gdouble lnM, z;
  glong np;

  if (z_initial == 0.0)
  {
	g_message ("nc_data_cluster_abundance_unbinned_init_from_text_file: User requested zi == 0.0. Using zi = %e\n", _NC_CLUSTER_ABUNDANCE_MIN_Z);
	z_initial = _NC_CLUSTER_ABUNDANCE_MIN_Z;
  }

  dca->area_survey = area_survey;
  dca->opt = opt;

  if (opt & NC_CLUSTER_ABUNDANCE_OBS_ZM)
  {
	dca->obs.lnMi = lnMi;
	dca->obs.lnMf = lnMf;
	dca->obs.zi = z_initial;
	dca->obs.zf = z_final;

	dca->obs.photoz_sigma0 = photoz_sigma0;
	dca->obs.photoz_bias = photoz_bias;
	dca->obs.lnM_sigma0 = lnM_sigma0;
	dca->obs.lnM_bias = lnM_bias;
  }
  else
  {
	dca->real.lnMi = lnMi;
	dca->real.lnMf = lnMf;
	dca->real.zi = z_initial;
	dca->real.zf = z_final;

	dca->real.photoz_sigma0 = photoz_sigma0;
	dca->real.photoz_bias = photoz_bias;
	dca->real.lnM_sigma0 = lnM_sigma0;
	dca->real.lnM_bias = lnM_bias;
  }

  if (filename == NULL)
	g_error ("It must pass the name of the file.");
  catalog = fopen (filename, "r");

  while ((c = fgetc(catalog)) != EOF)
	if (c == '\n') nlines++;
  rewind (catalog);

  np = 0;
  for (i = 0; i < nlines; i++)
  {
	gint nread = fscanf(catalog, " %lg %lg \n", &z, &lnM);
	if (!nread)
	  g_error ("nc_data_cluster_abundance_unbinned_init_from_text_file[fscanf]: cant find data");

	//printf("z = %.5g lnM = %.5g fscanf = %d\n", z, lnM, nread);

	if ((lnM >= lnMi) && (z <= z_final))
	{
	  np++;
	}
  }
  rewind (catalog);

  if (opt & NC_CLUSTER_ABUNDANCE_OBS_ZM)
  {
	if (dca->obs.z_lnM == NULL)
	  dca->obs.z_lnM = gsl_matrix_alloc (np, 2);
	else if (dca->obs.z_lnM->size1 != np)
	{
	  gsl_matrix_free (dca->obs.z_lnM);
	  dca->obs.z_lnM = gsl_matrix_alloc (np, 2);
	}
	//gsl_vector_ulong *v_mass = gsl_vector_alloc (np);
	dca->np = np;

	for (i = 0; i < nlines; i++)
	{
	  fscanf(catalog, "%lg %lg\n", &z, &lnM);
	  //printf ("lnM = %g z = %g [%d]\n", lnM, z, j);

	  if ((lnM >= lnMi) && (lnM <= lnMf) && (z <= z_final) && (z >= z_initial))
	  {
		//gsl_vector_set (v_mass, i, mass);
		gsl_matrix_set (dca->obs.z_lnM, j, 0, z);
		gsl_matrix_set (dca->obs.z_lnM, j, 1, lnM);
		j++;

		//printf ("%g %g\n", z, lnM);
	  }
	}
  }
  if (opt & NC_CLUSTER_ABUNDANCE_REAL_ZM)
  {
	if (dca->real.z_lnM == NULL)
	  dca->real.z_lnM = gsl_matrix_alloc (np, 2);
	else if (dca->obs.z_lnM->size1 != np)
	{
	  gsl_matrix_free (dca->real.z_lnM);
	  dca->real.z_lnM = gsl_matrix_alloc (np, 2);
	}
	//gsl_vector_ulong *v_mass = gsl_vector_alloc (np);
	dca->np = np;

	for (i = 0; i < nlines; i++)
	{
	  fscanf(catalog, "%lg %lg\n", &z, &lnM);
	  //printf ("lnM = %g z = %g [%d]\n", lnM, z, j);

	  if ((lnM >= lnMi) && (lnM <= lnMf) && (z <= z_final) && (z >= z_initial))
	  {
		//gsl_vector_set (v_mass, i, mass);
		gsl_matrix_set (dca->real.z_lnM, j, 0, z);
		gsl_matrix_set (dca->real.z_lnM, j, 1, lnM);
		j++;

		//printf ("%g %g\n", z, lnM);
	  }
	}
  }
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

  if (dca->opt & NC_CLUSTER_ABUNDANCE_PHOTOZ)
  {
	for (i = 0; i < dca->np; i++)
	{
	  const gdouble z_i = gsl_matrix_get (dca->obs.z_lnM, i, 0);
	  gsl_histogram_increment (hist, z_i);
	}
  }
  else
  {
	for (i = 0; i < dca->np; i++)
	{
	  const gdouble z_i = gsl_matrix_get (dca->real.z_lnM, i, 0);
	  gsl_histogram_increment (hist, z_i);
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

  if (dca->opt & NC_CLUSTER_ABUNDANCE_OBS_ZM)
  {
	for (i = 0; i < dca->np; i++)
	{
	  const gdouble zi_obs = gsl_matrix_get (dca->obs.z_lnM, i, 0);
	  const gdouble lnMi_obs = gsl_matrix_get (dca->obs.z_lnM, i, 1);

	  gsl_histogram2d_increment (hist, lnMi_obs, zi_obs);
	}
  }
  else
  {
	for (i = 0; i < dca->np; i++)
	{
	  const gdouble zi_real = gsl_matrix_get (dca->real.z_lnM, i, 0);
	  const gdouble lnMi_real = gsl_matrix_get (dca->real.z_lnM, i, 1);

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
  glong i;
  gint tfields = 4;

  gchar extname[] = "ClusterAbundance";   /* extension name */
  gchar *ttype[] = { "Z_REAL", "lnM_REAL", "Z_OBS", "lnM_OBS" };
  gchar *tform[] = { "1D", "1D", "1D", "1D" };
  gchar *tunit[] = { "REAL REDSHIFT", "REAL LN_MASS", "OBSERVATIONAL REDSHIFT", "OBSERVATIONAL LN_MASS" };

  /* initialize status before calling fitsio routines */
  status = 0;

  if (overwrite && g_file_test (filename, G_FILE_TEST_EXISTS))
	g_unlink (filename);

  /* create new FITS file */
  if (fits_create_file (&fptr, filename, &status))
	NC_FITS_ERROR (status);
  /*
   if (fits_create_img (fptr,  bitpix, naxis, naxes, &status))
   NC_FITS_ERROR (status);

   if (fits_write_date (fptr, &status))
   NC_FITS_ERROR (status);
   */
  /* move to 1st HDU  */
  /*
   if (fits_movabs_hdu (fptr, 2, &hdutype, &status))
   NC_FITS_ERROR (status);
   */
  /* append a new empty binary table onto the FITS file */
  if (fits_create_tbl (fptr, BINARY_TBL, dca->np, tfields, ttype, tform,
                       tunit, extname, &status) )
	NC_FITS_ERROR (status);

  if (fits_write_key(fptr, TDOUBLE, "ZI_T", &dca->real.zi, "Minimum true redshift", &status))
	NC_FITS_ERROR(status);
  if (fits_write_key(fptr, TDOUBLE, "ZF_T", &dca->real.zf, "Maximum true redshift", &status))
	NC_FITS_ERROR(status);
  if (fits_write_key(fptr, TDOUBLE, "ZI_O", &dca->obs.zi, "Minimum observational redshift", &status))
	NC_FITS_ERROR(status);
  if (fits_write_key(fptr, TDOUBLE, "ZF_O", &dca->obs.zf, "Maximum observational redshift", &status))
	NC_FITS_ERROR(status);

  {
	gdouble sarea_d = dca->area_survey / gsl_pow_2 (M_PI / 180.0);
	gdouble Mi_real = exp (dca->real.lnMi);
	gdouble Mf_real = exp (dca->real.lnMf);
	gdouble Mi_obs = exp (dca->obs.lnMi);
	gdouble Mf_obs = exp (dca->obs.lnMf);
	if (fits_write_key(fptr, TDOUBLE, "MI_T", &Mi_real, "Minimum true mass in h^{-1} * M_sun", &status))
	  NC_FITS_ERROR(status);
	if (fits_write_key(fptr, TDOUBLE, "MF_T", &Mf_real, "Maximum true mass in h^{-1} * M_sun", &status))
	  NC_FITS_ERROR(status);
	if (fits_write_key(fptr, TDOUBLE, "MI_O", &Mi_obs, "Minimum observational mass in h^{-1} * M_sun", &status))
	  NC_FITS_ERROR(status);
	if (fits_write_key(fptr, TDOUBLE, "MF_O", &Mf_obs, "Maximum observational mass in h^{-1} * M_sun", &status))
	  NC_FITS_ERROR(status);
	if (fits_write_key(fptr, TDOUBLE, "AREA", &sarea_d, "Survey area in degree square", &status))
	  NC_FITS_ERROR(status);
  }

  {
	gboolean Mobs, photoz;

	Mobs = dca->opt & NC_CLUSTER_ABUNDANCE_MOBS;
	photoz = dca->opt & NC_CLUSTER_ABUNDANCE_PHOTOZ;

	if (fits_write_key(fptr, TLOGICAL, "M_OBS", &Mobs, NULL, &status))
	  NC_FITS_ERROR(status);
	if (fits_write_key(fptr, TLOGICAL, "PHOTOZ", &photoz, NULL, &status))
	  NC_FITS_ERROR(status);

	if (Mobs)
	{
	  if (fits_write_key(fptr, TDOUBLE, "M_OBS_S0", &dca->obs.lnM_sigma0, NULL, &status))
		NC_FITS_ERROR(status);
	}
	if (photoz)
	{
	  if (fits_write_key(fptr, TDOUBLE, "PHOTOZS0", &dca->obs.photoz_sigma0, NULL, &status))
		NC_FITS_ERROR(status);
	}
  }

  if (dca->opt & NC_CLUSTER_ABUNDANCE_REAL_ZM)
  {
	for (i = 0; i < dca->np; i++)
	{
	  if (fits_write_col (fptr, TDOUBLE, 1, i + 1, 1, 1, gsl_matrix_ptr (dca->real.z_lnM, i, 0), &status))
		NC_FITS_ERROR(status);
	  if (fits_write_col (fptr, TDOUBLE, 2, i + 1, 1, 1, gsl_matrix_ptr (dca->real.z_lnM, i, 1), &status))
		NC_FITS_ERROR(status);
	  //printf ("Cat fit %g %g\n", gsl_matrix_get (dca->z_lnM, i, 0), gsl_matrix_get (dca->z_lnM, i, 1));
	}
  }
  else if (dca->opt & NC_CLUSTER_ABUNDANCE_OBS_ZM)
  {
	for (i = 0; i < dca->np; i++)
	{
	  if (fits_write_col (fptr, TDOUBLE, 3, i + 1, 1, 1, gsl_matrix_ptr (dca->obs.z_lnM, i, 0), &status))
		NC_FITS_ERROR(status);
	  if (fits_write_col (fptr, TDOUBLE, 4, i + 1, 1, 1, gsl_matrix_ptr (dca->obs.z_lnM, i, 1), &status))
		NC_FITS_ERROR(status);
	  //printf ("Cat fit %g %g\n", gsl_matrix_get (dca->z_lnM, i, 0), gsl_matrix_get (dca->z_lnM, i, 1));
	}
  }

  if ( fits_close_file(fptr, &status) )
	NC_FITS_ERROR(status);

  return;
}

/**
 * nc_cluster_abundance_catalog_load:
 * @data: a #NcData.
 * @filename: name of the file
 * @opt: a #NcClusterAbundanceOpt.
 *
 * FIXME
 *
 */
void
nc_cluster_abundance_catalog_load (NcData *data, gchar *filename, NcClusterAbundanceOpt opt)
{
  NcDataClusterAbundance *dca = (NcDataClusterAbundance *) NC_DATA_DATA (data);
  gint  status, hdutype, anynul, i;
  gchar comment[FLEN_COMMENT], ZT_name[100], MT_name[100], ZO_name[100], MO_name[100];
  gint zt_index, Mt_index, zo_index, Mo_index;
  gboolean is_lnM = TRUE;
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

  if (fits_read_key_dbl (fptr, "ZI_T", &dca->real.zi, comment, &status))
	g_error ("Fits file does not contain ZI_T in the header indicating the initial true redshift. Use [col #ZI_T=...] to add this information.");

  if (dca->real.zi == 0.0)
  {
	g_message ("nc_cluster_abundance_catalog_load: Catalog requested zi_true == 0.0. Using zi_true = %e\n", _NC_CLUSTER_ABUNDANCE_MIN_Z);
	dca->real.zi = _NC_CLUSTER_ABUNDANCE_MIN_Z;
  }

  if (fits_read_key_dbl (fptr, "ZF_T", &dca->real.zf, comment, &status))
	g_error ("Fits file does not contain ZF_T in the header indicating the final true redshift. Use [col #ZF_T=...] to add this information.");

  if (fits_read_key_dbl (fptr, "ZI_O", &dca->obs.zi, comment, &status))
	g_error ("Fits file does not contain ZI_O in the header indicating the initial observational redshift. Use [col #ZI_O=...] to add this information.");

  if (dca->obs.zi == 0.0)
  {
	g_message ("nc_cluster_abundance_catalog_load: Catalog requested zi_obs == 0.0. Using zi_obs = %e\n", _NC_CLUSTER_ABUNDANCE_MIN_Z);
	dca->obs.zi = _NC_CLUSTER_ABUNDANCE_MIN_Z;
  }

  if (fits_read_key_dbl (fptr, "ZF_O", &dca->obs.zf, comment, &status))
	g_error ("Fits file does not contain ZF_O in the header indicating the final observational redshift. Use [col #ZF_O=...] to add this information.");

  {
	gdouble Mi_real, Mf_real, Mi_obs, Mf_obs;
	if (fits_read_key_dbl (fptr, "MI_T", &Mi_real, comment, &status))
	  g_error ("Fits file does not contain MI_T in the header indicating the minimum true mass. Use [col #MI_T=...] to add this information.");
	if (fits_read_key_dbl (fptr, "MF_T", &Mf_real, comment, &status))
	  g_error ("Fits file does not contain MF_T in the header indicating the maximum true mass. Use [col #MF_T=...] to add this information.");

	dca->real.lnMi = log (Mi_real);
	dca->real.lnMf = log (Mf_real);

	if (fits_read_key_dbl (fptr, "MI_O", &Mi_obs, comment, &status))
	  g_error ("Fits file does not contain MI_O in the header indicating the minimum observational mass. Use [col #MI_O=...] to add this information.");
	if (fits_read_key_dbl (fptr, "MF_O", &Mf_obs, comment, &status))
	  g_error ("Fits file does not contain MF_O in the header indicating the maximum observational mass. Use [col #MF_O=...] to add this information.");

	dca->obs.lnMi = log (Mi_obs);
	dca->obs.lnMf = log (Mf_obs);
  }

  if (fits_read_key_dbl (fptr, "AREA", &dca->area_survey, comment, &status))
	g_error ("Fits file does not contain AREA in the header indicating the survey area (degree square). Use [col #AREA=...] to add this information.");
  dca->area_survey *= gsl_pow_2 (M_PI / 180.0);

  {
	gboolean Mobs, photoz;
	if (fits_read_key_log (fptr, "M_OBS", &Mobs, comment, &status))
	  g_error ("Fits file does not contain M_OBS in the header indicating the usage of mass observable relations. Use [col #M_OBS=T|F] to set it to TRUE or FALSE.");
	if (Mobs)
	{
	  if (fits_read_key_dbl (fptr, "M_OBS_S0", &dca->obs.lnM_sigma0, comment, &status))
		g_error ("Fits file does not contain M_OBS_S0 although it has M_OBS set true. Use [col #M_OBS_S0=...] to set its value.");
	  dca->opt = dca->opt | NC_CLUSTER_ABUNDANCE_MOBS;
	}
	else
	  dca->obs.lnM_sigma0 = 0.0;

	if (fits_read_key_log (fptr, "PHOTOZ", &photoz, comment, &status))
	  g_error ("Fits file does not contain PHOTOZ in the header indicating the usage of photometric redshift error. Use [col #PHOTOZ=T|F] to set it to TRUE or FALSE.");
	if (photoz)
	{
	  dca->opt = dca->opt | NC_CLUSTER_ABUNDANCE_PHOTOZ;
	  if (fits_read_key_dbl (fptr, "PHOTOZS0", &dca->obs.photoz_sigma0, comment, &status))
		g_error ("Fits file does not contain PHOTOZS0 although it has PHOTOZ set true. Use [col #PHOTOZS0=...] to set its value.");
	}
	else
	  dca->obs.photoz_sigma0 = 0.0;
  }
  dca->opt = dca->opt | opt;

  if (fits_get_colname (fptr, CASEINSEN, "ZT", ZT_name, &zt_index, &status))
	g_error ("Column ZT|zt not found, invalid fits file.");
  if (fits_get_colname (fptr, CASEINSEN, "ZO", ZO_name, &zo_index, &status))
	g_error ("Column ZO|zo not found, invalid fits file.");


  fits_get_colname (fptr, CASEINSEN, "lnM_T", MT_name, &Mt_index, &status);
  if ((status == COL_NOT_FOUND) || (status == COL_NOT_UNIQUE))
  {
	status = 0;
	fits_get_colname (fptr, CASEINSEN, "M200_T", MT_name, &Mt_index, &status);
	if ((status == COL_NOT_FOUND) || (status == COL_NOT_UNIQUE))
	{
	  status = 0;
	  fits_get_colname (fptr, CASEINSEN, "MT", MT_name, &Mt_index, &status);
	  if ((status == COL_NOT_FOUND) || (status == COL_NOT_UNIQUE))
		g_error ("Column (lnM_T|M200_T|MT) not found, invalid fits file.");
	}
	is_lnM = FALSE;
  }

  fits_get_colname (fptr, CASEINSEN, "lnM_O", MO_name, &Mo_index, &status);
  if ((status == COL_NOT_FOUND) || (status == COL_NOT_UNIQUE))
  {
	status = 0;
	fits_get_colname (fptr, CASEINSEN, "M200", MO_name, &Mo_index, &status);
	if ((status == COL_NOT_FOUND) || (status == COL_NOT_UNIQUE))
	{
	  status = 0;
	  fits_get_colname (fptr, CASEINSEN, "MO", MO_name, &Mo_index, &status);
	  if ((status == COL_NOT_FOUND) || (status == COL_NOT_UNIQUE))
		g_error ("Column (lnM_O|M200|MO) not found, invalid fits file.");
	}
	is_lnM = FALSE;
  }

  {
	gchar *select = g_strdup_printf (
	                                 "%s >= % 20.15g && %s <= % 20.15g && %s >= % 20.15g && %s <= % 20.15g %s >= % 20.15g && %s <= % 20.15g && %s >= % 20.15g && %s <= % 20.15g",
	                                 ZT_name, dca->real.zi, ZT_name, dca->real.zf,
	                                 MT_name, is_lnM ? dca->real.lnMi : exp(dca->real.lnMi),
	                                 MT_name, is_lnM ? dca->real.lnMf : exp(dca->real.lnMf),
	                                 ZO_name, dca->obs.zi, ZO_name, dca->obs.zf,
	                                 MO_name, is_lnM ? dca->obs.lnMi : exp(dca->obs.lnMi),
	                                 MO_name, is_lnM ? dca->obs.lnMf : exp(dca->obs.lnMf));
	fits_select_rows (fptr, fptr, select, &status);
	NC_FITS_ERROR (status);
	g_free (select);
  }

  if ( fits_read_key_lng (fptr, "NAXIS2", &dca->np, comment, &status) )
	NC_FITS_ERROR(status);

  if (dca->opt & NC_CLUSTER_ABUNDANCE_REAL_ZM)
  {
	if (dca->real.z_lnM != NULL)
	{
	  if (dca->real.z_lnM->size1 != dca->np)
	  {
		gsl_matrix_free (dca->real.z_lnM);
		dca->real.z_lnM = gsl_matrix_alloc (dca->np, 2);
	  }
	}
	else
	  dca->real.z_lnM = gsl_matrix_alloc (dca->np, 2);

	if (!is_lnM)
	{
	  for (i = 0; i < dca->np; i++)
	  {
		if (fits_read_col_dbl (fptr, zt_index, i + 1, 1, 1, 0.0, gsl_matrix_ptr (dca->real.z_lnM, i, 0), &anynul, &status))
		  NC_FITS_ERROR(status);
		if (fits_read_col_dbl (fptr, Mt_index, i + 1, 1, 1, 0.0, gsl_matrix_ptr (dca->real.z_lnM, i, 1), &anynul, &status))
		  NC_FITS_ERROR(status);
		*gsl_matrix_ptr (dca->real.z_lnM, i, 1) = log (*gsl_matrix_ptr (dca->real.z_lnM, i, 1));
	  }
	}
	else
	{
	  for (i = 0; i < dca->np; i++)
	  {
		if (fits_read_col_dbl (fptr, zt_index, i + 1, 1, 1, 0.0, gsl_matrix_ptr (dca->real.z_lnM, i, 0), &anynul, &status))
		  NC_FITS_ERROR(status);
		if (fits_read_col_dbl (fptr, Mt_index, i + 1, 1, 1, 0.0, gsl_matrix_ptr (dca->real.z_lnM, i, 1), &anynul, &status))
		  NC_FITS_ERROR(status);
	  }
	}
  }

  if (dca->opt & NC_CLUSTER_ABUNDANCE_OBS_ZM)
  {
	if (dca->obs.z_lnM != NULL)
	{
	  if (dca->obs.z_lnM->size1 != dca->np)
	  {
		gsl_matrix_free (dca->obs.z_lnM);
		dca->obs.z_lnM = gsl_matrix_alloc (dca->np, 2);
	  }
	}
	else
	  dca->obs.z_lnM = gsl_matrix_alloc (dca->np, 2);

	if (!is_lnM)
	{
	  for (i = 0; i < dca->np; i++)
	  {
		if (fits_read_col_dbl (fptr, zo_index, i + 1, 1, 1, 0.0, gsl_matrix_ptr (dca->obs.z_lnM, i, 0), &anynul, &status))
		  NC_FITS_ERROR(status);
		if (fits_read_col_dbl (fptr, Mo_index, i + 1, 1, 1, 0.0, gsl_matrix_ptr (dca->obs.z_lnM, i, 1), &anynul, &status))
		  NC_FITS_ERROR(status);
		*gsl_matrix_ptr (dca->obs.z_lnM, i, 1) = log (*gsl_matrix_ptr (dca->obs.z_lnM, i, 1));
	  }
	}
	else
	{
	  for (i = 0; i < dca->np; i++)
	  {
		if (fits_read_col_dbl (fptr, zo_index, i + 1, 1, 1, 0.0, gsl_matrix_ptr (dca->obs.z_lnM, i, 0), &anynul, &status))
		  NC_FITS_ERROR(status);
		if (fits_read_col_dbl (fptr, Mo_index, i + 1, 1, 1, 0.0, gsl_matrix_ptr (dca->obs.z_lnM, i, 1), &anynul, &status))
		  NC_FITS_ERROR(status);
	  }
	}
  }

  if ( fits_close_file(fptr, &status) )
	NC_FITS_ERROR(status);

  nc_data_init (data);

  return;
}

/**
 * nc_cluster_matching_catalog_save:
 * @data: a #NcData
 * @filename: name of the file
 * @overwrite: FIXME
 *
 * FIXME
 *
 */
void
nc_cluster_matching_catalog_save (NcData *data, gchar *filename, gboolean overwrite)
{
  NcDataClusterAbundance *dca = (NcDataClusterAbundance *) NC_DATA_DATA (data);
  /*******************************************************************/
  /* Create a binary table extension                                 */
  /*******************************************************************/
  fitsfile *fptr;       /* pointer to the FITS file, defined in fitsio.h */
  gint status, hdutype;
  glong i;
  gint n_zbins = gsl_histogram2d_nx (dca->completeness);
  /* gint n_logMbins = gsl_histogram2d_ny (dca->completeness); FIXME */

  gint bitpix   =  SHORT_IMG;
  glong naxis   =   0;
  glong naxes[] = { 0, 0 };

  gint tfields   = 3;

  gchar extname[] = "NcusterAbundance";   /* extension name */
  gchar *ttype[] = { "Z", "MASS_RANK" , "U_TWOWAY"};
  gchar *tform[] = { "1D", "1E", "1I"};
  gchar *tunit[] = { "REDSHIFT", "MASS", "U_TWOWAY"};

  gint tfields_2   = 6;

  gchar extname_2[] = "Statistics";   /* extension name */
  gchar *ttype_2[] = { "Z_MIN", "Z_MAX", "MBINS" , "C", "P", "MDIFF_SIGMA"};
  gchar *tform_2[] = { "1E", "1E", "E", "E", "E", "E"};
  gchar *tunit_2[] = { "Z_MIN", "Z_MAX", "MBINS" , "COMPLETENESS", "PURITY", "MDIFF_SIGMA"};

  /* initialize status before calling fitsio routines */
  status = 0;

  if (overwrite && g_file_test (filename, G_FILE_TEST_EXISTS))
	g_unlink (filename);

  /* create new FITS file */
  if (fits_create_file (&fptr, filename, &status))
	NC_FITS_ERROR (status);

  if (fits_create_img (fptr,  bitpix, naxis, naxes, &status))
	NC_FITS_ERROR (status);

  if (fits_write_date (fptr, &status))
	NC_FITS_ERROR (status);

  /* move to 2nd HDU  */
  if ( fits_movabs_hdu (fptr, 2, &hdutype, &status))
	NC_FITS_ERROR (status);

  /* append a new empty binary table onto the FITS file */
  if ( fits_create_tbl (fptr, BINARY_TBL, dca->np, tfields, ttype, tform,
                        tunit, extname, &status) )
	NC_FITS_ERROR (status);

  {
	gdouble sarea_d = dca->area_survey / gsl_pow_2 (M_PI / 180.0);
	if (fits_write_key(fptr, TDOUBLE, "AREA", &sarea_d, "Survey area in degree square", &status))
	  NC_FITS_ERROR(status);
  }

  {
	gdouble M;
	gint u_twoway = 1;
	for (i = 0; i < dca->np; i++)
	{
	  if (fits_write_col (fptr, TDOUBLE, 1, i + 1, 1, 1, gsl_matrix_ptr (dca->real.z_lnM, i, 0), &status))
		NC_FITS_ERROR(status);
	  M = exp (gsl_matrix_get (dca->real.z_lnM, i, 1));
	  if (fits_write_col (fptr, TDOUBLE, 2, i + 1, 1, 1, &M, &status))
		NC_FITS_ERROR(status);
	  if (fits_write_col (fptr, TINT, 3, i + 1, 1, 1, &u_twoway, &status))
		NC_FITS_ERROR(status);
	}
  }

  /* move to 3rd HDU  */
  if ( fits_movabs_hdu (fptr, 4, &hdutype, &status))
	NC_FITS_ERROR (status);

  if ( fits_create_tbl (fptr, BINARY_TBL, n_zbins, tfields_2, ttype_2, tform_2, tunit_2, extname_2, &status) )
	NC_FITS_ERROR (status);

  {
	gdouble z_min, z_max, lnM_min, lnM_max, M;
	gsl_vector *M_nodes = NULL;
	gint m_bins = gsl_histogram2d_ny (dca->completeness);
	M_nodes = gsl_vector_alloc (m_bins + 1);
	for(i = 0; i < m_bins; i++)
	{
	  gsl_histogram2d_get_yrange (dca->completeness, i, &lnM_min, &lnM_max);
	  M = exp (lnM_min); /* corrigir!!! quero transformar em log_10 (M), FIXME*/
	  gsl_vector_set (M_nodes, i, M);
	}
	M = exp (lnM_max);
	gsl_vector_set (M_nodes, m_bins, M);
	for (i = 0; i < n_zbins; i++)
	{
	  gsl_histogram2d_get_xrange (dca->completeness, i, &z_min, &z_max);
	  if (fits_write_col (fptr, TDOUBLE, 1, i + 1, 1, 1, &z_min, &status))
		NC_FITS_ERROR(status);
	  if (fits_write_col (fptr, TDOUBLE, 2, i + 1, 1, 1, &z_max, &status))
		NC_FITS_ERROR(status);
	  if (fits_write_col (fptr, TDOUBLE, 3, i + 1, 1, 1, gsl_vector_ptr (M_nodes, 0), &status))
		NC_FITS_ERROR(status);
	  //if (fits_write_col (fptr, TDOUBLE, 4, i + 1, 1, 1, gsl_matrix_ptr (complet_matrix, i, 0), &status))
	  //NC_FITS_ERROR(status);
	  //if (fits_write_col (fptr, TDOUBLE, 5, i + 1, 1, 1, gsl_matrix_ptr (pur_matrix, i, 0), &status))
	  //NC_FITS_ERROR(status);
	  //if (fits_write_col (fptr, TDOUBLE, 6, i + 1, 1, 1, gsl_matrix_ptr (sdM_matrix, i, 0), &status))
	  //NC_FITS_ERROR(status);
	}
  }

  if ( fits_close_file(fptr, &status) )
	NC_FITS_ERROR(status);

  return;
}

/**
 * nc_cluster_matching_catalog_load:
 * @data: a #NcData.
 * @filename: name of the file.
 * @opt: a #NcClusterAbundanceOpt.
 *
 * FIXME
 * P.S. This function was not adapted to have columns with true and observable values of z and lnM.
 *
 */
void
nc_cluster_matching_catalog_load (NcData *data, gchar *filename, NcClusterAbundanceOpt opt)
{
  NcDataClusterAbundance *dca = (NcDataClusterAbundance *) NC_DATA_DATA (data);
  gint  status, hdutype, anynul, i;
  gchar comment[FLEN_COMMENT]; //Z_name[100], M_name[100];
  gint z_index, Mobs_index, utwoway_index;
  //gboolean is_lnM = TRUE;
  fitsfile *fptr;

  gint j, typecode, zmin_index, zmax_index, mbins_index, completeness_index, purity_index, sigmaM_index;
  glong nbins_z, nbins_lnM, width, np_complet;
  gsl_vector *z_range = NULL;
  gsl_vector *mass_range = NULL;
  gsl_matrix *complet_matrix, *pur_matrix, *sdM_matrix;

  status = 0;

  if (filename == NULL)
	g_error ("nc_cluster_matching_catalog_load: null filename");

  if (fits_open_file (&fptr, filename, READONLY, &status))
	NC_FITS_ERROR(status);

  if (fits_movabs_hdu (fptr, 4, &hdutype, &status))
	NC_FITS_ERROR (status);
  if (fits_read_key_lng (fptr, "NAXIS2", &nbins_z, comment, &status))
	NC_FITS_ERROR(status);
  if (fits_get_colnum (fptr, CASEINSEN, "ZMIN", &zmin_index, &status))
	NC_FITS_ERROR(status);
  if (fits_get_colnum (fptr, CASEINSEN, "ZMAX", &zmax_index, &status))
	NC_FITS_ERROR(status);
  if (fits_get_colnum (fptr, CASEINSEN, "MBINS", &mbins_index, &status))
	NC_FITS_ERROR(status);
  if (fits_get_colnum (fptr, CASEINSEN, "C", &completeness_index, &status))
	NC_FITS_ERROR(status);
  if (fits_get_coltype (fptr, completeness_index, &typecode, &np_complet, &width, &status))
	NC_FITS_ERROR(status);
  //printf ("np = %ld", np_complet);
  if (fits_get_colnum (fptr, CASEINSEN, "P", &purity_index, &status))
	NC_FITS_ERROR(status);
  if (fits_get_colnum (fptr, CASEINSEN, "MDIFF_SIGMA", &sigmaM_index, &status))
	NC_FITS_ERROR(status);

  if (fits_get_coltype (fptr, mbins_index, &typecode, &nbins_lnM, &width, &status))
	NC_FITS_ERROR(status);
  nbins_lnM--; /* The function above return the number of nodes and we need the number of bins. */
  //printf ("# typecode = %d nbins_logM = %ld width = %ld\n", typecode, nbins_logM, width);

  z_range = gsl_vector_alloc (nbins_z + 1);
  mass_range = gsl_vector_alloc (nbins_lnM + 1);
  complet_matrix = gsl_matrix_alloc (nbins_z, nbins_lnM);
  pur_matrix = gsl_matrix_alloc (nbins_z, nbins_lnM);
  sdM_matrix = gsl_matrix_alloc (nbins_z, nbins_lnM);
  dca->completeness = gsl_histogram2d_alloc (nbins_z, nbins_lnM);
  dca->purity = gsl_histogram2d_alloc (nbins_z, nbins_lnM);
  dca->sd_lnM = gsl_histogram2d_alloc (nbins_z, nbins_lnM);

  for (j = 0; j < nbins_z; j++)
  {
	if (fits_read_col_dbl (fptr, zmin_index, j + 1, 1, 1, 0.0, gsl_vector_ptr (z_range, j), &anynul, &status))
	  NC_FITS_ERROR(status);
	if (fits_read_col_dbl (fptr, completeness_index, j + 1, 1, np_complet, 0.0, gsl_matrix_ptr (complet_matrix, j, 0), &anynul, &status))
	  NC_FITS_ERROR(status);
	if (fits_read_col_dbl (fptr, purity_index, j + 1, 1, np_complet, 0.0, gsl_matrix_ptr (pur_matrix, j, 0), &anynul, &status))
	  NC_FITS_ERROR(status);
	if (fits_read_col_dbl (fptr, sigmaM_index, j + 1, 1, np_complet, 0.0, gsl_matrix_ptr (sdM_matrix, j, 0), &anynul, &status))
	  NC_FITS_ERROR(status);
	//printf ("%.3g\n", gsl_vector_get (z_range, j));
  }

  if (fits_read_col_dbl (fptr, zmax_index, nbins_z, 1, 1, 0.0, gsl_vector_ptr (z_range, nbins_z), &anynul, &status))
	NC_FITS_ERROR(status);
  if (fits_read_col_dbl (fptr, mbins_index, 1, 1, nbins_lnM + 1, 0.0, gsl_vector_ptr (mass_range, 0), &anynul, &status))
	NC_FITS_ERROR(status);
  //printf ("%.3g\n", gsl_vector_get (z_range, nbins_z));

  gsl_vector_scale (mass_range, M_LN10); /* Matching catalog provides logM and we need lnM. */
  gsl_matrix_scale (sdM_matrix, M_LN10); /* Standard deviation is obtained for (logMobs - logMtrue) so we have to multiply by ln10. */

  gsl_histogram2d_set_ranges (dca->completeness, gsl_vector_ptr (z_range, 0), nbins_z + 1, gsl_vector_ptr (mass_range, 0), nbins_lnM + 1);
  gsl_histogram2d_set_ranges (dca->purity, gsl_vector_ptr (z_range, 0), nbins_z + 1, gsl_vector_ptr (mass_range, 0), nbins_lnM + 1);
  gsl_histogram2d_set_ranges (dca->sd_lnM, gsl_vector_ptr (z_range, 0), nbins_z + 1, gsl_vector_ptr (mass_range, 0), nbins_lnM + 1);
  //printf ("# zbin = %.3g mbin = %.3g | %ld %ld\n", gsl_vector_get (z_range, 0), gsl_vector_get (mass_range, 0), nbins_z, nbins_logM);

  for (j = 0; j < nbins_z; j++)
  {
	//gdouble xlower, xupper;
	gdouble x = gsl_vector_get (z_range, j);
	//gsl_histogram2d_get_xrange (completeness, j, &xlower, &xupper);
	//printf ("x = %.3g xlower = %.3g xupper = %.3g\n", x, xlower, xupper);
	for (i = 0; i < nbins_lnM; i++)
	{
	  gdouble y = gsl_vector_get (mass_range, i);
	  //printf ("y = %.3g\n", y);
	  gsl_histogram2d_accumulate (dca->completeness, x, y, gsl_matrix_get (complet_matrix, j, i));
	  gsl_histogram2d_accumulate (dca->purity, x, y, gsl_matrix_get (pur_matrix, j, i));
	  gsl_histogram2d_accumulate (dca->sd_lnM, x, y, gsl_matrix_get (sdM_matrix, j, i));
	  //printf ("p_matrix = %.3g p = %.3g\n", gsl_matrix_get (complet_matrix, j, i), gsl_histogram2d_get (completeness, j, i));
	  //printf ("p_matrix = %.3g p = %.3g\n", gsl_matrix_get (pur_matrix, j, i), gsl_histogram2d_get (purity, j, i));
	  //printf ("sd_matrix = %.3g sd_M = %.3g\n", gsl_matrix_get (sdM_matrix, j, i), gsl_histogram2d_get (sd_M, j, i));
	}
  }

  gsl_vector_free (z_range);
  gsl_vector_free (mass_range);
  gsl_matrix_free (complet_matrix);
  gsl_matrix_free (pur_matrix);
  gsl_matrix_free (sdM_matrix);

  dca->real.zi = gsl_histogram2d_xmin (dca->sd_lnM);
  dca->real.zf = gsl_histogram2d_xmax (dca->sd_lnM);
  dca->real.lnMi = gsl_histogram2d_ymin (dca->sd_lnM);
  dca->real.lnMf = gsl_histogram2d_ymax (dca->sd_lnM);

  printf ("zi = %.3g zf = %.3g lnMi = %.3g lnMf = %.3g\n", dca->real.zi, dca->real.zf, dca->real.lnMi, dca->real.lnMf);

  if (fits_movabs_hdu (fptr, 2, &hdutype, &status))
	NC_FITS_ERROR (status);

  if (hdutype != BINARY_TBL)
	g_error ("%s (%d): Matching catalog is not binary!\n", __FILE__, __LINE__);

  if (fits_read_key_dbl (fptr, "AREA", &dca->area_survey, comment, &status))
	g_error ("Fits file do not contain AREA in the header indicating the survey area (degree square). Use [col #AREA=...] to add this information.");
  dca->area_survey *= gsl_pow_2 (M_PI / 180.0);

  {
	gboolean Mobs_local, selection, photoz;
	if (fits_read_key_log (fptr, "M_OBS_LOCAL", &Mobs_local, comment, &status))
	  g_error ("Fits file does not contain M_OBS_LOCAL in the header indicating the usage of mass observable relations. Use [col #M_OBS_LOCAL=T|F] to set it to TRUE or FALSE.");
	if (Mobs_local)
	  dca->opt = dca->opt | NC_CLUSTER_ABUNDANCE_MOBS_LOCAL;

	if (fits_read_key_log (fptr, "SELECTION", &selection, comment, &status))
	  g_error ("Fits file does not contain SELECTION in the header indicating the usage of completeness and purity. Use [col #SELECTION=T|F] to set it to TRUE or FALSE.");
	if (selection)
	  dca->opt = dca->opt | (NC_CLUSTER_ABUNDANCE_COMPLETENESS | NC_CLUSTER_ABUNDANCE_PURITY);

	if (fits_read_key_log (fptr, "PHOTOZ", &photoz, comment, &status))
	  g_error ("Fits file does not contain PHOTOZ in the header indicating the usage of photometric redshift error. Use [col #PHOTOZ=T|F] to set it to TRUE or FALSE.");
	if (photoz)
	{
	  dca->opt = dca->opt | NC_CLUSTER_ABUNDANCE_PHOTOZ;
	  if (fits_read_key_dbl (fptr, "PHOTOZS0", &dca->real.photoz_sigma0, comment, &status))
		g_error ("Fits file does not contain PHOTOZS0 although it has PHOTOZ set true. Use [col #PHOTOZS0=...] to set its value.");
	}
	else
	  dca->real.photoz_sigma0 = 0.0;
  }
  dca->opt = dca->opt | opt;
  printf ("opt = %d\n", dca->opt);
  //if (fits_get_colnum (fptr, CASEINSEN, "TWOWAY", &twoway_index, &status))
  //NC_FITS_ERROR(status);
  if (fits_get_colnum (fptr, CASEINSEN, "U_TWOWAY", &utwoway_index, &status))
	NC_FITS_ERROR(status);
  if (fits_get_colnum (fptr, CASEINSEN, "Z", &z_index, &status))
	NC_FITS_ERROR(status);
  if (fits_get_colnum (fptr, CASEINSEN, "MASS_RANK", &Mobs_index, &status))
	NC_FITS_ERROR(status);
  //if (fits_get_colnum (fptr, CASEINSEN, "U_MASS_MATCH", &Mtrue_index, &status))
  //NC_FITS_ERROR(status);

  //if (fits_select_rows (fptr, fptr, "(TWOWAY == 1) || (U_TWOWAY == 1)", &status))
  //NC_FITS_ERROR(status);

  if (fits_select_rows (fptr, fptr, "(U_TWOWAY == 1)", &status))
	NC_FITS_ERROR(status);

  if (fits_read_key_lng (fptr, "NAXIS2", &dca->np, comment, &status))
	NC_FITS_ERROR(status);
  dca->real.z_lnM = gsl_matrix_alloc (dca->np, 2);
  //printf ("%d %d %d %ld\n", z_index, Mobs_index, Mtrue_index, np);

  for (i = 0; i < dca->np; i++)
  {
	gdouble M;
	if (fits_read_col_dbl (fptr, z_index, i + 1, 1, 1, 0.0, gsl_matrix_ptr (dca->real.z_lnM, i, 0), &anynul, &status))
	  NC_FITS_ERROR(status);
	if (fits_read_col_dbl (fptr, Mobs_index, i + 1, 1, 1, 0.0, &M, &anynul, &status))
	  NC_FITS_ERROR(status);
	gsl_matrix_set (dca->real.z_lnM, i, 1, log (M));
	//if (fits_read_col_dbl (fptr, Mtrue_index, i + 1, 1, 1, 0.0, gsl_matrix_ptr (dca->z_lnM, i, 2), &anynul, &status))
	//NC_FITS_ERROR(status);
	//printf ("z = %.3g Mobs = %.3g\n", gsl_matrix_get (dca->z_lnM, i, 0), gsl_matrix_get (dca->z_lnM, i, 1));
  }

  if ( fits_close_file(fptr, &status) )
	NC_FITS_ERROR(status);

  nc_data_init (data);

  return;
}

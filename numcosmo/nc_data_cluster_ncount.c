/***************************************************************************
 *            nc_data_cluster_ncount.c
 *
 *  Tue Apr  6 01:11:23 2010
 *  Copyright  2010  Mariana Penna Lima
 *  <pennalima@gmail.com>
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
 * SECTION:nc_data_cluster_ncount
 * @title: Cluster number count data
 * @short_description: FIXME
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc_data_cluster_ncount.h"
#include "nc_data_cluster_poisson.h"

#include "math/util.h"
#include "math/ncm_func_eval.h"
#include "math/ncm_cfg.h"

#include <glib/gstdio.h>
#include <gsl/gsl_randist.h>
#ifdef NUMCOSMO_HAVE_CFITSIO
#include <fitsio.h>
#endif /* NUMCOSMO_HAVE_CFITSIO */

enum
{
  PROP_0,
  PROP_CAD,
  PROP_SIZE,
};

G_DEFINE_TYPE (NcDataClusterNCount, nc_data_cluster_ncount, NCM_TYPE_DATA);

static void
nc_data_cluster_ncount_init (NcDataClusterNCount *ncount)
{
  ncount->cad            = NULL;
  ncount->z              = NULL;
  ncount->m              = NULL;
  ncount->lnM_true       = NULL;
  ncount->z_true         = NULL;
  ncount->z_obs          = NULL;
  ncount->z_obs_params   = NULL;
  ncount->lnM_obs        = NULL;
  ncount->lnM_obs_params = NULL;
  ncount->area_survey    = 0.0;
  ncount->np             = 0.0;
  ncount->log_np_fac     = 0.0;
  ncount->use_true_data  = FALSE;
  ncount->completeness   = NULL;
  ncount->purity         = NULL;
  ncount->sd_lnM         = NULL;
  ncount->fiducial       = FALSE;
  ncount->seed           = 0;
  ncount->rnd_name       = NULL;

}

static void
nc_data_cluster_ncount_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcDataClusterNCount *ncount = NC_DATA_CLUSTER_NCOUNT (object);
  g_return_if_fail (NC_IS_DATA_CLUSTER_NCOUNT (object));

  switch (prop_id)
  {
    case PROP_CAD:
      ncount->cad = g_value_dup_object (value);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_data_cluster_ncount_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcDataClusterNCount *ncount = NC_DATA_CLUSTER_NCOUNT (object);

  g_return_if_fail (NC_IS_DATA_CLUSTER_NCOUNT (object));

  switch (prop_id)
  {
    case PROP_CAD:
      g_value_set_object (value, ncount->cad);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_data_cluster_ncount_dispose (GObject *object)
{
  NcDataClusterNCount *ncount = NC_DATA_CLUSTER_NCOUNT (object);
  
  nc_cluster_redshift_clear (&ncount->z);
  nc_cluster_mass_clear (&ncount->m);

  ncm_vector_clear (&ncount->lnM_true);
  ncm_vector_clear (&ncount->z_true);

  ncm_matrix_clear (&ncount->z_obs);
  ncm_matrix_clear (&ncount->z_obs_params);
  ncm_matrix_clear (&ncount->lnM_obs);
  ncm_matrix_clear (&ncount->lnM_obs_params);
  
  /* Chain up : end */
  G_OBJECT_CLASS (nc_data_cluster_ncount_parent_class)->dispose (object);
}

static void
nc_data_cluster_ncount_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_data_cluster_ncount_parent_class)->finalize (object);
}

static guint _nc_data_cluster_ncount_get_length (NcmData *data);
static NcmData *_nc_data_cluster_ncount_dup (NcmData *data);
static void _nc_data_cluster_ncount_begin (NcmData *data);
static void _nc_data_cluster_ncount_prepare (NcmData *data, NcmMSet *mset);
static void _nc_data_cluster_ncount_resample (NcmData *data, NcmMSet *mset);
static void _nc_data_cluster_ncount_m2lnL_val (NcmData *data, NcmMSet *mset, gdouble *m2lnL);

static void
nc_data_cluster_ncount_class_init (NcDataClusterNCountClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);
  NcmDataClass *data_class   = NCM_DATA_CLASS (klass);

  object_class->set_property = &nc_data_cluster_ncount_set_property;
  object_class->get_property = &nc_data_cluster_ncount_get_property;
  object_class->dispose      = &nc_data_cluster_ncount_dispose;
  object_class->finalize     = &nc_data_cluster_ncount_finalize;

  g_object_class_install_property (object_class,
                                   PROP_CAD,
                                   g_param_spec_object ("cluster-abundance",
                                                        NULL,
                                                        "Cluster Abundance",
                                                        NC_TYPE_CLUSTER_ABUNDANCE,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  
  data_class->name = "Cluster abundance unbinned";

  data_class->get_length = &_nc_data_cluster_ncount_get_length;
  data_class->dup        = &_nc_data_cluster_ncount_dup;
  data_class->begin      = &_nc_data_cluster_ncount_begin;
  data_class->prepare    = &_nc_data_cluster_ncount_prepare;
  data_class->resample   = &_nc_data_cluster_ncount_resample;
  data_class->m2lnL_val  = &_nc_data_cluster_ncount_m2lnL_val;
}

static void
_ca_data_copy (NcDataClusterNCount *dest, NcDataClusterNCount *src)
{
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
do { if (src->name != NULL) \
{ \
if (dest->name == NULL) \
dest->name = ncm_vector_dup (src->name); \
if (ncm_vector_len (dest->name) != ncm_vector_len (src->name)) \
{ \
ncm_vector_free (dest->name); \
dest->name = ncm_vector_dup (src->name); \
} \
ncm_vector_memcpy (dest->name, src->name); \
} } while (FALSE)

#define _MATRIX_COPY(name) \
do { if (src->name != NULL) \
{ \
if (dest->name == NULL) \
dest->name = ncm_matrix_dup (src->name); \
if ((NCM_MATRIX_NROWS (dest->name) != NCM_MATRIX_NROWS (src->name)) || (NCM_MATRIX_NCOLS (dest->name) != NCM_MATRIX_NCOLS (src->name))) \
{ \
ncm_matrix_free (dest->name); \
dest->name = ncm_matrix_dup (src->name); \
} \
ncm_matrix_memcpy (dest->name, src->name); \
} } while (FALSE)

  _VECTOR_COPY (lnM_true);
  _VECTOR_COPY (z_true);

  _MATRIX_COPY (z_obs);
  _MATRIX_COPY (z_obs_params);
  _MATRIX_COPY (lnM_obs);
  _MATRIX_COPY (lnM_obs_params);
}

static NcmData *
_nc_data_cluster_ncount_dup (NcmData *data)
{
  NcDataClusterNCount *ncount = NC_DATA_CLUSTER_NCOUNT (data);
  NcDataClusterNCount *ncount_dup = NC_DATA_CLUSTER_NCOUNT (nc_data_cluster_ncount_new (ncount->cad));
  _ca_data_copy (ncount_dup, ncount);

  return NCM_DATA (ncount_dup);
}

static guint 
_nc_data_cluster_ncount_get_length (NcmData *data) 
{ 
  return NC_DATA_CLUSTER_NCOUNT (data)->np; 
}

static void
_nc_data_cluster_ncount_begin (NcmData *data)
{
  NcDataClusterNCount *ncount = NC_DATA_CLUSTER_NCOUNT (data);
  ncount->log_np_fac = lgamma (ncount->np + 1);
}

/**
 * nc_data_cluster_ncount_new:
 * @cad: FIXME
 *
 * FIXME
 * 
 * Returns: FIXME
 */
NcmData *
nc_data_cluster_ncount_new (NcClusterAbundance *cad)
{
  NcDataClusterNCount *ncount = g_object_new (NC_TYPE_DATA_CLUSTER_NCOUNT,
                                              "custer-abundance", cad,
                                              NULL);

  return NCM_DATA (ncount);
}

/**
 * nc_data_cluster_ncount_ref:
 * @ncount: FIXME
 *
 * FIXME
 * 
 * Returns: (transfer full): FIXME
 */
NcDataClusterNCount *
nc_data_cluster_ncount_ref (NcDataClusterNCount *ncount)
{
  return g_object_ref (ncount);
}

/**
 * nc_data_cluster_ncount_free:
 * @ncount: FIXME
 *
 * FIXME
 * 
 */
void
nc_data_cluster_ncount_free (NcDataClusterNCount *ncount)
{
  g_object_unref (ncount);
}

/**
 * nc_data_cluster_ncount_init_from_fits_file:
 * @data: a #NcmData
 * @filename: name of the file
 *
 * FIXME
 * 
 */
void
nc_data_cluster_ncount_init_from_fits_file (NcmData *data, gchar *filename)
{
  g_assert_not_reached ();
}

static void
_nc_data_cluster_ncount_binned_f (NcmMSet *mset, gpointer obj, const gdouble *x, gdouble *f)
{
  //NcHICosmo *cosmo = NC_HICOSMO (ncm_mset_peek (mset, NC_HICOSMO_ID));
  //NcClusterAbundance *cad = NC_CLUSTER_ABUNDANCE (obj);
  //f[0] = nc_cluster_abundance_N_val (cad, model, cad->lnMi, cad->lnMf, cad->zi, x[0]);
  g_assert_not_reached ();
  return;
}

/**
 * nc_data_cluster_ncount_binned_new_function: (skip)
 * @cad: a #NcClusterAbundance
 *
 * FIXME
 *
 * Returns: FIXME
 */
NcmMSetFunc *
nc_data_cluster_ncount_binned_new_function (NcClusterAbundance *cad)
{
  NcmMSetFunc *func = ncm_mset_func_new (_nc_data_cluster_ncount_binned_f, 1, 1, cad, (GDestroyNotify) nc_cluster_abundance_free);
  return func;
}

/**
 * nc_data_cluster_ncount_binned_lnM_z_new:
 * @cad: a #NcClusterAbundance
 *
 * FIXME
 *
 * Returns: FIXME
 */
NcmData *
nc_data_cluster_ncount_binned_lnM_z_new (NcClusterAbundance *cad)
{
  g_assert_not_reached ();
  return NULL;
}

/************************************************************************************************************
 * Unbinned number count data                                                                               *
 ************************************************************************************************************/

static void
_nc_data_cluster_ncount_model_init (NcDataClusterNCount *ncount)
{
  NcClusterAbundance *cad = ncount->cad;

  if (ncount->z == NULL || ncount->m == NULL)
    g_error ("_nc_data_cluster_ncount_model_init: Cannot init NcClusterAbundance missing NcClusterRedshift or NcClusterMass object.");

  nc_cluster_abundance_set_redshift (cad, ncount->z);
  nc_cluster_abundance_set_mass (cad, ncount->m);

  cad->completeness  = ncount->completeness;
  cad->purity        = ncount->purity;
  cad->sd_lnM        = ncount->sd_lnM;

  cad->mfp->area_survey = ncount->area_survey;
}

static void
_nc_data_cluster_ncount_prepare (NcmData *data, NcmMSet *mset)
{
  NcDataClusterNCount *ncount = NC_DATA_CLUSTER_NCOUNT (data);
  NcHICosmo *cosmo = NC_HICOSMO (ncm_mset_peek (mset, NC_HICOSMO_ID));

  if (ncm_model_ctrl_update (ncount->cad->ctrl, NCM_MODEL (cosmo)))
  {
    nc_cluster_abundance_prepare (ncount->cad, cosmo);
  }
}

/**
 * _nc_data_cluster_ncount_resample:
 * @mset: a #NcmModel.
 * @model: FIXME
 * @data: FIXME
 *
 * This function generates random numbers which are used to obtain redshift
 * and mass (logarithm base e) values...
 *
 */
static void
_nc_data_cluster_ncount_resample (NcmData *data, NcmMSet *mset)
{
  NcDataClusterNCount *ncount = NC_DATA_CLUSTER_NCOUNT (data);
  NcClusterAbundance *cad = ncount->cad;
  NcHICosmo *cosmo = NC_HICOSMO (ncm_mset_peek (mset, NC_HICOSMO_ID));
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
        nc_cluster_mass_resample (cad->m, cosmo, lnM_true, z_true, lnMi_obs, lnMi_obs_params) )
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

  if (ncount->lnM_true != NULL)
    ncm_vector_free (ncount->lnM_true);
  ncount->lnM_true = ncm_vector_new_array (lnM_true_array);
  g_array_unref (lnM_true_array);

  if (ncount->z_true != NULL)
    ncm_vector_free (ncount->z_true);
  ncount->z_true = ncm_vector_new_array (z_true_array);
  g_array_unref (z_true_array);

  if (ncount->z_obs != NULL)
    ncm_matrix_free (ncount->z_obs);
  ncount->z_obs = ncm_matrix_new_array (z_obs_array, z_obs_len);
  g_array_unref (z_obs_array);

  if (ncount->lnM_obs != NULL)
    ncm_matrix_free (ncount->lnM_obs);
  ncount->lnM_obs = ncm_matrix_new_array (lnM_obs_array, lnM_obs_len);
  g_array_unref (lnM_obs_array);

  if (z_obs_params_len > 0)
  {
    if (ncount->z_obs_params != NULL)
      ncm_matrix_free (ncount->z_obs_params);
    ncount->z_obs_params = ncm_matrix_new_array (z_obs_params_array, z_obs_params_len);
    g_array_unref (z_obs_params_array);
  }

  if (lnM_obs_params_len > 0)
  {
    if (ncount->lnM_obs_params != NULL)
      ncm_matrix_free (ncount->lnM_obs_params);
    ncount->lnM_obs_params = ncm_matrix_new_array (lnM_obs_params_array, lnM_obs_params_len);
    g_array_unref (lnM_obs_params_array);
  }

  ncount->np = NCM_MATRIX_NROWS (ncount->z_obs);
  //printf ("# Generated %ld | expected % 20.15g\n", ncount->np, nc_cluster_abundance_n (cad, NC_HICOSMO (ncm_mset_peek (mset, NC_HICOSMO_ID))));

  g_free (zi_obs);
  g_free (zi_obs_params);
  g_free (lnMi_obs);
  g_free (lnMi_obs_params);
}

typedef struct
{
  NcClusterAbundance *cad;
  NcDataClusterNCount *ncount;
  NcHICosmo *cosmo;
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
    gdouble *lnMn_obs = ncm_matrix_ptr (evald2n->ncount->lnM_obs, n, 0);
    gdouble *lnMn_obs_params = evald2n->ncount->lnM_obs_params != NULL ? ncm_matrix_ptr (evald2n->ncount->lnM_obs_params, n, 0) : NULL;
    gdouble *zn_obs = ncm_matrix_ptr (evald2n->ncount->z_obs, n, 0);
    gdouble *zn_obs_params = evald2n->ncount->z_obs_params != NULL ? ncm_matrix_ptr (evald2n->ncount->z_obs_params, n, 0) : NULL;
    const gdouble mlnLn = -log (nc_cluster_abundance_z_p_lnm_p_d2n (evald2n->cad, evald2n->cosmo, lnMn_obs, lnMn_obs_params, zn_obs, zn_obs_params));
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
    const gdouble lnMn = ncm_vector_get (evald2n->ncount->lnM_true, n);
    gdouble *zn_obs = ncm_matrix_ptr (evald2n->ncount->z_obs, n, 0);
    gdouble *zn_obs_params = evald2n->ncount->z_obs_params != NULL ? ncm_matrix_ptr (evald2n->ncount->z_obs_params, n, 0) : NULL;
    const gdouble mlnLn = -log (nc_cluster_abundance_z_p_d2n (evald2n->cad, evald2n->cosmo, lnMn, zn_obs, zn_obs_params));
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
    const gdouble zn = ncm_vector_get (evald2n->ncount->z_true, n);
    gdouble *lnMn_obs = ncm_matrix_ptr (evald2n->ncount->lnM_obs, n, 0);
    gdouble *lnMn_obs_params = evald2n->ncount->lnM_obs_params != NULL ? ncm_matrix_ptr (evald2n->ncount->lnM_obs_params, n, 0) : NULL;
    const gdouble mlnLn = -log (nc_cluster_abundance_lnm_p_d2n (evald2n->cad, evald2n->cosmo, lnMn_obs, lnMn_obs_params, zn));
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
    const gdouble lnMn = ncm_vector_get (evald2n->ncount->lnM_true, n);
    const gdouble zn = ncm_vector_get (evald2n->ncount->z_true, n);
    const gdouble mlnLn = -log (nc_cluster_abundance_d2n (evald2n->cad, evald2n->cosmo, lnMn, zn));
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
    const gdouble lnMn = ncm_vector_get (evald2n->ncount->lnM_true, n);
    const gdouble zn = ncm_vector_get (evald2n->ncount->z_true, n);
    const gdouble mlnLn = -log (nc_cluster_abundance_intp_d2n (evald2n->cad, evald2n->cosmo, lnMn, zn));
    m2lnL += mlnLn;
  }

  G_LOCK (save_m2lnL);
  *evald2n->m2lnL += m2lnL;
  G_UNLOCK (save_m2lnL);
}

static void
_nc_data_cluster_ncount_m2lnL_val (NcmData *data, NcmMSet *mset, gdouble *m2lnL)
{
  NcDataClusterNCount *ncount = NC_DATA_CLUSTER_NCOUNT (data);
  NcClusterAbundance *cad = ncount->cad;
  NcHICosmo *cosmo = NC_HICOSMO (ncm_mset_peek (mset, NC_HICOSMO_ID));
  GTimer *bench = g_timer_new ();

  *m2lnL = 0.0;

  if (ncount->use_true_data)
  {
    _Evald2N evald2n = {cad, ncount, cosmo, m2lnL};
    g_assert (ncount->z_true);
    g_assert (ncount->lnM_true);
    ncm_func_eval_threaded_loop (&_eval_intp_d2n, 0, ncount->np, &evald2n);
  }
  else
  {
    NcClusterRedshiftImpl z_impl = nc_cluster_redshift_impl (ncount->z);
    NcClusterMassImpl lnM_impl = nc_cluster_mass_impl (ncount->m);
    gboolean z_p = z_impl & NC_CLUSTER_REDSHIFT_P;
    gboolean lnM_p = lnM_impl & NC_CLUSTER_MASS_P;
    if (z_p && lnM_p)
    {
      _Evald2N evald2n = {cad, ncount, cosmo, m2lnL};
      ncm_func_eval_threaded_loop (&_eval_z_p_lnm_p_d2n, 0, ncount->np, &evald2n);
    }
    else if (z_p && !lnM_p)
    {
      g_assert (ncount->lnM_true);
      _Evald2N evald2n = {cad, ncount, cosmo, m2lnL};
      ncm_func_eval_threaded_loop (&_eval_z_p_d2n, 0, ncount->np, &evald2n);
    }
    else if (!z_p && lnM_p)
    {
      g_assert (ncount->z_true);
      _Evald2N evald2n = {cad, ncount, cosmo, m2lnL};
      ncm_func_eval_threaded_loop (&_eval_lnm_p_d2n, 0, ncount->np, &evald2n);
    }
    else
    {
      g_assert (ncount->z_true);
      g_assert (ncount->lnM_true);
      _Evald2N evald2n = {cad, ncount, cosmo, m2lnL};
      ncm_func_eval_threaded_loop (&_eval_d2n, 0, ncount->np, &evald2n);

    }
  }

  *m2lnL += (ncount->log_np_fac + nc_cluster_abundance_n (cad, cosmo));

  *m2lnL *= 2.0;

  g_timer_destroy (bench);
}

/**
 * nc_data_cluster_ncount_true_data:
 * @data: a #NcmData.
 * @use_true_data: FIXME
 *
 * FIXME
 *
 */
void
nc_data_cluster_ncount_true_data (NcmData *data, gboolean use_true_data)
{
  NcDataClusterNCount *ncount = NC_DATA_CLUSTER_NCOUNT (data);
  if (use_true_data)
  {
    g_assert (ncount->lnM_true != NULL);
    g_assert (ncount->z_true != NULL);
  }
  ncount->use_true_data = use_true_data;
}

static void _nc_data_cluster_ncount_model_init (NcDataClusterNCount *ncount);

/**
 * nc_data_cluster_ncount_init_from_sampling:
 * @data: a #NcmData.
 * @mset: a #NcmMSet.
 * @clusterz: a #NcClusterRedshift.
 * @clusterm: a #NcClusterMass.
 * @area_survey: area in units of square degrees.
 *
 * FIXME
 *
 */
void
nc_data_cluster_ncount_init_from_sampling (NcmData *data, NcmMSet *mset, NcClusterRedshift *clusterz, NcClusterMass *clusterm, gdouble area_survey)
{
  NcDataClusterNCount *ncount = NC_DATA_CLUSTER_NCOUNT (data);

  ncount->area_survey = area_survey;

  ncount->z = nc_cluster_redshift_ref (clusterz);
  ncount->m = nc_cluster_mass_ref (clusterm);

  _nc_data_cluster_ncount_model_init (ncount);
  ncm_data_resample (data, mset);

  ncm_data_set_init (data);
}

/**
 * nc_data_cluster_ncount_bin_data: (skip)
 * @data: a #NcmData
 * @nodes: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
NcmData *
nc_data_cluster_ncount_bin_data (NcmData *data, gsl_vector *nodes)
{
  NcDataClusterNCount *ncount = NC_DATA_CLUSTER_NCOUNT (data);
  NcmData *data_cpoisson;
  gsl_histogram *hist;
  gint i;

  g_assert (nodes->size > 1);
  g_assert (nodes->stride == 1);

  data_cpoisson = nc_data_cluster_poisson_new (ncount);

  hist = gsl_histogram_alloc (nodes->size - 1);
  gsl_histogram_set_ranges (hist, nodes->data, nodes->size);

  {
    for (i = 0; i < ncount->np; i++)
    {
      const gdouble z_i = 0.0;//gsl_matrix_get (ncount->real.z_lnM, i, 0);
      gsl_histogram_increment (hist, z_i);
      g_assert_not_reached ();
    }
  }

  nc_data_poisson_init_from_histogram (data_cpoisson, hist);

  gsl_histogram_free (hist);

  return data_cpoisson;
}

/**
 * nc_data_cluster_ncount_hist_lnM_z: (skip)
 * @data: a #NcmData.
 * @lnM_nodes: FIXME
 * @z_nodes: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gsl_histogram2d *
nc_data_cluster_ncount_hist_lnM_z (NcmData *data, gsl_vector *lnM_nodes, gsl_vector *z_nodes)
{
  NcDataClusterNCount *ncount = NC_DATA_CLUSTER_NCOUNT (data);
  gsl_histogram2d *hist;
  gint i;

  g_assert (data->init);
  g_assert (lnM_nodes->size > 1);
  g_assert (lnM_nodes->stride == 1);
  g_assert (z_nodes->size > 1);
  g_assert (z_nodes->stride == 1);

  //ca_binned = nc_data_cluster_ncount_binned_lnM_z_new (cad); /* I have to make this function. */
  //ca_binned = nc_data_cluster_ncount_new (cad);

  hist = gsl_histogram2d_alloc (lnM_nodes->size - 1, z_nodes->size - 1);
  gsl_histogram2d_set_ranges (hist, lnM_nodes->data, lnM_nodes->size, z_nodes->data, z_nodes->size);

  {
    for (i = 0; i < ncount->np; i++)
    {
      const gdouble zi_real = 0.0;//gsl_matrix_get (ncount->real.z_lnM, i, 0);
      const gdouble lnMi_real = 0.0;//gsl_matrix_get (ncount->real.z_lnM, i, 1);
      g_assert_not_reached ();

      gsl_histogram2d_increment (hist, lnMi_real, zi_real);
    }
  }

  return hist;
}

#ifdef NUMCOSMO_HAVE_CFITSIO

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
 * nc_data_cluster_ncount_catalog_save:
 * @data: a #NcmData
 * @filename: name of the file
 * @overwrite: FIXME
 *
 * FIXME
 *
 */
void
nc_data_cluster_ncount_catalog_save (NcmData *data, gchar *filename, gboolean overwrite)
{
  NcDataClusterNCount *ncount = NC_DATA_CLUSTER_NCOUNT (data);
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

  guint z_obs_len = nc_cluster_redshift_obs_len (ncount->z);
  guint z_obs_params_len = nc_cluster_redshift_obs_params_len (ncount->z);
  guint lnM_obs_len = nc_cluster_mass_obs_len (ncount->m);
  guint lnM_obs_params_len = nc_cluster_mass_obs_params_len (ncount->m);

  g_ptr_array_set_free_func (tform_array, g_free);

  g_ptr_array_add (ttype_array, "Z_OBS");
  g_ptr_array_add (tform_array, g_strdup_printf ("%dD", z_obs_len));
  g_ptr_array_add (tunit_array, "REDSHIFT OBS");

  g_ptr_array_add (ttype_array, "LNM_OBS");
  g_ptr_array_add (tform_array, g_strdup_printf ("%dD", lnM_obs_len));
  g_ptr_array_add (tunit_array, "MASS OBS");

  if (ncount->z_true != NULL)
  {
    g_ptr_array_add (ttype_array, "Z_TRUE");
    g_ptr_array_add (tform_array, g_strdup ("1D"));
    g_ptr_array_add (tunit_array, "TRUE REDSHIFT");
  }

  if (ncount->lnM_true != NULL)
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
  if (fits_create_tbl (fptr, BINARY_TBL, ncount->np, tfields, (gchar **)ttype_array->pdata, (gchar **)tform_array->pdata,
                       (gchar **)tunit_array->pdata, extname, &status) )
    NC_FITS_ERROR (status);

  {
    gchar *z_ser = ncm_cfg_serialize_to_string (G_OBJECT (ncount->z), FALSE);
    gchar *lnM_ser = ncm_cfg_serialize_to_string (G_OBJECT (ncount->m), FALSE);

    if (fits_write_key_longstr (fptr, "Z_OBJ", z_ser, "Serialized redshift object", &status))
      NC_FITS_ERROR(status);
    if (fits_write_key_longstr (fptr, "LNM_OBJ", lnM_ser, "Serialized mass object", &status))
      NC_FITS_ERROR(status);

    g_free (z_ser);
    g_free (lnM_ser);
  }

  {
    gdouble sarea_d = ncount->area_survey / gsl_pow_2 (M_PI / 180.0);
    if (fits_write_key(fptr, TDOUBLE, "AREA", &sarea_d, "Survey area in degree square", &status))
      NC_FITS_ERROR(status);
  }

  {
    guint colnum = 1;
    if (fits_write_col (fptr, TDOUBLE, colnum, 1, 1, ncount->np, ncm_matrix_ptr (ncount->z_obs, 0, 0), &status))
      NC_FITS_ERROR(status);

    colnum++;
    if (fits_write_col (fptr, TDOUBLE, colnum, 1, 1, ncount->np, ncm_matrix_ptr (ncount->lnM_obs, 0, 0), &status))
      NC_FITS_ERROR(status);

    if (ncount->z_true != NULL)
    {
      colnum++;
      if (fits_write_col (fptr, TDOUBLE, colnum, 1, 1, ncount->np, ncm_vector_ptr (ncount->z_true, 0), &status))
        NC_FITS_ERROR(status);
    }

    if (ncount->lnM_true != NULL)
    {
      colnum++;
      if (fits_write_col (fptr, TDOUBLE, colnum, 1, 1, ncount->np, ncm_vector_ptr (ncount->lnM_true, 0), &status))
        NC_FITS_ERROR(status);
    }

    if (z_obs_params_len > 0)
    {
      colnum++;
      if (fits_write_col (fptr, TDOUBLE, colnum, 1, 1, ncount->np, ncm_matrix_ptr (ncount->z_obs_params, 0, 0), &status))
        NC_FITS_ERROR(status);
    }

    if (lnM_obs_params_len > 0)
    {
      colnum++;
      if (fits_write_col (fptr, TDOUBLE, colnum, 1, 1, ncount->np, ncm_matrix_ptr (ncount->lnM_obs_params, 0, 0), &status))
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
 * nc_data_cluster_ncount_catalog_load:
 * @data: a #NcmData.
 * @filename: name of the file
 *
 * FIXME
 *
 */
void
nc_data_cluster_ncount_catalog_load (NcmData *data, gchar *filename)
{
  NcDataClusterNCount *ncount = NC_DATA_CLUSTER_NCOUNT (data);
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

  if (ncount->z != NULL)
    nc_cluster_redshift_free (ncount->z);
  if (ncount->m != NULL)
    nc_cluster_mass_free (ncount->m);

  {
    gchar *z_ser = NULL;
    gchar *lnM_ser = NULL;

    if (fits_read_key_longstr (fptr, "Z_OBJ", &z_ser, comment, &status))
      NC_FITS_ERROR (status);
    if (fits_read_key_longstr (fptr, "LNM_OBJ", &lnM_ser, comment, &status))
      NC_FITS_ERROR (status);

    {
      gint z_ser_len = strlen (z_ser);
      gint lnM_ser_len = strlen (lnM_ser);

      if (z_ser[z_ser_len - 1] == '&')
        z_ser[z_ser_len - 1] = ' ';

      if (lnM_ser[lnM_ser_len - 1] == '&')
        lnM_ser[lnM_ser_len - 1] = ' ';

    }

    ncount->z = nc_cluster_redshift_new_from_name (z_ser);
    ncount->m = nc_cluster_mass_new_from_name (lnM_ser);

    g_free (z_ser);
    g_free (lnM_ser);
  }

  {
    glong nrows;
    if (fits_get_num_rows (fptr, &nrows, &status))
      NC_FITS_ERROR(status);
    ncount->np = nrows;
  }

  if (fits_read_key_dbl (fptr, "AREA", &ncount->area_survey, comment, &status))
    g_error ("Fits file does not contain AREA in the header indicating the survey area (degree square). Use [col #AREA=...] to add this information.");
  ncount->area_survey *= gsl_pow_2 (M_PI / 180.0);

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

    if (nc_cluster_redshift_obs_len (ncount->z) != z_obs_rp)
      g_error ("NcClusterRedshift object has observables length %d but fits has %ld.",
               nc_cluster_redshift_obs_len (ncount->z), z_obs_rp);

    if (nc_cluster_mass_obs_len (ncount->m) != lnM_obs_rp)
      g_error ("NcClusterMass object has observables length %d but fits has %ld.",
               nc_cluster_mass_obs_len (ncount->m), lnM_obs_rp);

    if (ncount->z_obs)
      ncm_matrix_free (ncount->z_obs);
    ncount->z_obs = ncm_matrix_new (ncount->np, z_obs_rp);

    if (ncount->lnM_obs)
      ncm_matrix_free (ncount->lnM_obs);
    ncount->lnM_obs = ncm_matrix_new (ncount->np, lnM_obs_rp);

    if (fits_read_col (fptr, TDOUBLE, z_obs_i, 1, 1, ncount->np, NULL, ncm_matrix_ptr (ncount->z_obs, 0, 0), NULL, &status))
      NC_FITS_ERROR(status);

    if (fits_read_col (fptr, TDOUBLE, lnM_obs_i, 1, 1, ncount->np, NULL, ncm_matrix_ptr (ncount->lnM_obs, 0, 0), NULL, &status))
      NC_FITS_ERROR(status);
  }

  {
    gint z_obs_params_i, lnM_obs_params_i;
    gint z_obs_params_tc, lnM_obs_params_tc;
    glong z_obs_params_rp, lnM_obs_params_rp;
    glong z_obs_params_w, lnM_obs_params_w;

    if (fits_get_colnum (fptr, CASESEN, "Z_OBS_PARAMS", &z_obs_params_i, &status))
    {
      if (nc_cluster_redshift_obs_params_len (ncount->z) > 0)
        g_error ("NcClusterRedshift object has observable parameters length %d but fits has 0.",
                 nc_cluster_redshift_obs_params_len (ncount->z));
    }
    else
    {
      if (fits_get_coltype (fptr, z_obs_params_i, &z_obs_params_tc, &z_obs_params_rp, &z_obs_params_w, &status))
        g_error ("Column Z_OBS_PARAMS info not found, invalid fits file.");

      if (nc_cluster_redshift_obs_params_len (ncount->z) != z_obs_params_rp)
        g_error ("NcClusterRedshift object has observable parameters length %d but fits has %ld.",
                 nc_cluster_redshift_obs_params_len (ncount->z), z_obs_params_rp);

      if (ncount->z_obs_params)
        ncm_matrix_free (ncount->z_obs_params);
      ncount->z_obs_params = ncm_matrix_new (ncount->np, z_obs_params_rp);

      if (fits_read_col (fptr, TDOUBLE, z_obs_params_i, 1, 1, ncount->np, NULL, ncm_matrix_ptr (ncount->z_obs_params, 0, 0), NULL, &status))
        NC_FITS_ERROR(status);
    }

    if (fits_get_colnum (fptr, CASESEN, "LNM_OBS_PARAMS", &lnM_obs_params_i, &status))
    {
      if (nc_cluster_mass_obs_params_len (ncount->m) > 0)
        g_error ("NcClusterMass object has observable parameters length %d but fits has 0.",
                 nc_cluster_mass_obs_params_len (ncount->m));
      status = 0;
    }
    else
    {
      if (fits_get_coltype (fptr, lnM_obs_params_i, &lnM_obs_params_tc, &lnM_obs_params_rp, &lnM_obs_params_w, &status))
        g_error ("Column LNM_OBS_PARAMS info not found, invalid fits file.");

      if (nc_cluster_mass_obs_params_len (ncount->m) != lnM_obs_params_rp)
        g_error ("NcClusterMass object has observable parameters length %d but fits has %ld.",
                 nc_cluster_mass_obs_params_len (ncount->m), lnM_obs_params_rp);

      if (ncount->lnM_obs_params)
        ncm_matrix_free (ncount->lnM_obs_params);
      ncount->lnM_obs_params = ncm_matrix_new (ncount->np, lnM_obs_params_rp);


      if (fits_read_col (fptr, TDOUBLE, lnM_obs_params_i, 1, 1, ncount->np, NULL, ncm_matrix_ptr (ncount->lnM_obs_params, 0, 0), NULL, &status))
        NC_FITS_ERROR(status);
    }
  }

  {
    gint z_true_i, lnM_true_i;
    if (ncount->z_true != NULL)
    {
      ncm_vector_free (ncount->z_true);
      ncount->z_true = NULL;
    }
    if (!fits_get_colnum (fptr, CASESEN, "Z_TRUE", &z_true_i, &status))
    {
      ncount->z_true = ncm_vector_new (ncount->np);
      if (fits_read_col (fptr, TDOUBLE, z_true_i, 1, 1, ncount->np, NULL, ncm_vector_ptr (ncount->z_true, 0), NULL, &status))
        NC_FITS_ERROR(status);
    }
    else
      status = 0;

    if (ncount->lnM_true != NULL)
    {
      ncm_vector_free (ncount->lnM_true);
      ncount->lnM_true = NULL;
    }

    if (!fits_get_colnum (fptr, CASESEN, "LNM_TRUE", &lnM_true_i, &status))
    {
      ncount->lnM_true = ncm_vector_new (ncount->np);
      if (fits_read_col (fptr, TDOUBLE, lnM_true_i, 1, 1, ncount->np, NULL, ncm_vector_ptr (ncount->lnM_true, 0), NULL, &status))
        NC_FITS_ERROR(status);
    }
    else
      status = 0;

  }

  if (fits_close_file (fptr, &status) )
    NC_FITS_ERROR(status);

  _nc_data_cluster_ncount_model_init (ncount);

  ncm_data_set_init (data);

  return;
}

#endif /* NUMCOSMO_HAVE_CFITSIO */

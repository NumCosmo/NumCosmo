/***************************************************************************
 *            nc_data_cluster_ncount.c
 *
 *  Tue Apr  6 01:11:23 2010
 *  Copyright  2010  Mariana Penna Lima
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) 2012 Mariana Penna Lima <pennalima@gmail.com>
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
 * @title: NcDataClusterNCount
 * @short_description: Cluster number count data.
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "data/nc_data_cluster_ncount.h"
#include "data/nc_data_cluster_poisson.h"
#include "nc_hireion.h"

#include "math/ncm_func_eval.h"
#include "math/ncm_serialize.h"
#include "math/ncm_cfg.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <glib/gstdio.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_statistics_double.h>
#ifdef NUMCOSMO_HAVE_CFITSIO
#include <fitsio.h>
#endif /* NUMCOSMO_HAVE_CFITSIO */
#endif /* NUMCOSMO_GIR_SCAN */

enum
{
  PROP_0,
  PROP_CAD,
  PROP_N_Z_OBS,
  PROP_N_Z_OBS_PARAMS,
  PROP_N_M_OBS,
  PROP_N_M_OBS_PARAMS,
  PROP_LNM_TRUE,
  PROP_Z_TRUE,
  PROP_Z_OBS,
  PROP_Z_OBS_PARAMS,
  PROP_LNM_OBS,
  PROP_LNM_OBS_PARAMS,
  PROP_SURVEY_AREA,
  PROP_USE_TRUE,
  PROP_BINNED,
  PROP_Z_NODES,
  PROP_LNM_NODES,
  PROP_FIDUCIAL,
  PROP_SEED,
  PROP_RNG_NAME,
  PROP_SIZE,
};

G_DEFINE_TYPE (NcDataClusterNCount, nc_data_cluster_ncount, NCM_TYPE_DATA);

static void
nc_data_cluster_ncount_init (NcDataClusterNCount *ncount)
{
  ncount->cad            = NULL;
  ncount->n_z_obs        = 0;
  ncount->n_z_obs_params = 0;
  ncount->n_M_obs        = 0;
  ncount->n_M_obs_params = 0;
  ncount->lnM_true       = NULL;
  ncount->z_true         = NULL;
  ncount->z_obs          = NULL;
  ncount->z_obs_params   = NULL;
  ncount->lnM_obs        = NULL;
  ncount->lnM_obs_params = NULL;
  ncount->m2lnL_a        = g_array_new (FALSE, FALSE, sizeof (gdouble));
  ncount->area_survey    = 0.0;
  ncount->np             = 0.0;
  ncount->log_np_fac     = 0.0;
  ncount->use_true_data  = FALSE;
  ncount->binned         = FALSE;
  ncount->z_nodes        = NULL;
  ncount->lnM_nodes      = NULL;
  ncount->purity         = NULL;
  ncount->sd_lnM         = NULL;
  ncount->z_lnM          = NULL;
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
      nc_cluster_abundance_clear (&ncount->cad);
      ncount->cad = g_value_dup_object (value);
      break;
    case PROP_N_Z_OBS:
      nc_data_cluster_ncount_set_n_z_obs (ncount, g_value_get_uint (value));
      break;
    case PROP_N_Z_OBS_PARAMS:
      nc_data_cluster_ncount_set_n_z_obs_params (ncount, g_value_get_uint (value));
      break;
    case PROP_N_M_OBS:
      nc_data_cluster_ncount_set_n_M_obs (ncount, g_value_get_uint (value));
      break;
    case PROP_N_M_OBS_PARAMS:
      nc_data_cluster_ncount_set_n_M_obs_params (ncount, g_value_get_uint (value));
      break;   
    case PROP_LNM_TRUE:
      nc_data_cluster_ncount_set_lnM_true (ncount, g_value_get_object (value));
      break;
    case PROP_Z_TRUE:
      nc_data_cluster_ncount_set_z_true (ncount, g_value_get_object (value));
      break;
    case PROP_Z_OBS:
      nc_data_cluster_ncount_set_z_obs (ncount, g_value_get_object (value));
      break;
    case PROP_Z_OBS_PARAMS:
      nc_data_cluster_ncount_set_z_obs_params (ncount, g_value_get_object (value));
      break;
    case PROP_LNM_OBS:
      nc_data_cluster_ncount_set_lnM_obs (ncount, g_value_get_object (value));
      break;
    case PROP_LNM_OBS_PARAMS:
      nc_data_cluster_ncount_set_lnM_obs_params (ncount, g_value_get_object (value));
      break;
    case PROP_SURVEY_AREA:
      ncount->area_survey = g_value_get_double (value);
      break;
    case PROP_USE_TRUE:
      ncount->use_true_data = g_value_get_boolean (value);
      break;
    case PROP_BINNED:
      nc_data_cluster_ncount_set_binned (ncount, g_value_get_boolean (value));
      break;
    case PROP_Z_NODES:
    {
      ncm_vector_clear (&ncount->z_nodes);
      ncount->z_nodes = g_value_dup_object (value);
      if (ncount->z_nodes != NULL)
        ncount->binned = FALSE;
      break;
    }
    case PROP_LNM_NODES:
    {
      ncm_vector_clear (&ncount->lnM_nodes);
      ncount->lnM_nodes = g_value_dup_object (value);
      if (ncount->lnM_nodes != NULL)
        ncount->binned = FALSE;
      break;
    }
    case PROP_FIDUCIAL:
      ncount->fiducial = g_value_get_boolean (value);
      break;
    case PROP_SEED:
      ncount->seed = g_value_get_uint64 (value);
      break;
    case PROP_RNG_NAME:
      g_clear_pointer (&ncount->rnd_name, g_free);
      ncount->rnd_name = g_value_dup_string (value);
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
    case PROP_N_Z_OBS:
      g_value_set_uint (value, ncount->n_z_obs);
      break;
    case PROP_N_Z_OBS_PARAMS:
      g_value_set_uint (value, ncount->n_z_obs_params);
      break;
    case PROP_N_M_OBS:
      g_value_set_uint (value, ncount->n_M_obs);
      break;
    case PROP_N_M_OBS_PARAMS:
      g_value_set_uint (value, ncount->n_M_obs_params);
      break;
    case PROP_LNM_TRUE:
      g_value_set_object (value, ncount->lnM_true);
      break;
    case PROP_Z_TRUE:
      g_value_set_object (value, ncount->z_true);
      break;
    case PROP_Z_OBS:
      g_value_set_object (value, ncount->z_obs);
      break;
    case PROP_Z_OBS_PARAMS:
      g_value_set_object (value, ncount->z_obs_params);
      break;
    case PROP_LNM_OBS:
      g_value_set_object (value, ncount->lnM_obs);
      break;
    case PROP_LNM_OBS_PARAMS:
      g_value_set_object (value, ncount->lnM_obs_params);
      break;
    case PROP_SURVEY_AREA:
      g_value_set_double (value, ncount->area_survey);
      break;
    case PROP_USE_TRUE:
      g_value_set_boolean (value, ncount->use_true_data);
      break;
    case PROP_BINNED:
      g_value_set_boolean (value, ncount->binned);
      break;
    case PROP_Z_NODES:
      g_value_set_object (value, ncount->z_nodes);
      break;
    case PROP_LNM_NODES:
      g_value_set_object (value, ncount->lnM_nodes);
      break;
    case PROP_FIDUCIAL:
      g_value_set_boolean (value, ncount->fiducial);
      break;
    case PROP_SEED:
      g_value_set_uint64 (value, ncount->seed);
      break;
    case PROP_RNG_NAME:
      g_value_set_string (value, ncount->rnd_name);
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

  nc_cluster_abundance_clear (&ncount->cad);

  ncm_vector_clear (&ncount->lnM_true);
  ncm_vector_clear (&ncount->z_true);

  ncm_matrix_clear (&ncount->z_obs);
  ncm_matrix_clear (&ncount->z_obs_params);
  ncm_matrix_clear (&ncount->lnM_obs);
  ncm_matrix_clear (&ncount->lnM_obs_params);

  ncm_vector_clear (&ncount->lnM_nodes);
  ncm_vector_clear (&ncount->z_nodes);  

  g_clear_pointer (&ncount->m2lnL_a, g_array_unref);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_data_cluster_ncount_parent_class)->dispose (object);
}

static void
nc_data_cluster_ncount_finalize (GObject *object)
{
  NcDataClusterNCount *ncount = NC_DATA_CLUSTER_NCOUNT (object);

  g_clear_pointer (&ncount->rnd_name, g_free);
  g_clear_pointer (&ncount->purity, gsl_histogram2d_free);
  g_clear_pointer (&ncount->sd_lnM, gsl_histogram2d_free);
  g_clear_pointer (&ncount->z_lnM, gsl_histogram2d_free);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_data_cluster_ncount_parent_class)->finalize (object);
}

static guint _nc_data_cluster_ncount_get_length (NcmData *data);
static void _nc_data_cluster_ncount_begin (NcmData *data);
static void _nc_data_cluster_ncount_prepare (NcmData *data, NcmMSet *mset);
static void _nc_data_cluster_ncount_resample (NcmData *data, NcmMSet *mset, NcmRNG *rng);
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
                                                        "Cluster abundance",
                                                        NC_TYPE_CLUSTER_ABUNDANCE,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_N_Z_OBS,
                                   g_param_spec_uint ("n-z-obs",
                                                      NULL,
                                                      "Number of redshift observables",
                                                      1, G_MAXUINT32, 1,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_N_Z_OBS_PARAMS,
                                   g_param_spec_uint ("n-z-obs-params",
                                                      NULL,
                                                      "Number of redshift observables parameters",
                                                      0, G_MAXUINT32, 0,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_N_M_OBS,
                                   g_param_spec_uint ("n-M-obs",
                                                      NULL,
                                                      "Number of mass observables",
                                                      1, G_MAXUINT32, 1,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_N_M_OBS_PARAMS,
                                   g_param_spec_uint ("n-M-obs-params",
                                                      NULL,
                                                      "Number of mass observables parameters",
                                                      0, G_MAXUINT32, 0,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_LNM_TRUE,
                                   g_param_spec_object ("lnM-true",
                                                        NULL,
                                                        "Clusters true masses",
                                                        NCM_TYPE_VECTOR,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));  
  g_object_class_install_property (object_class,
                                   PROP_Z_TRUE,
                                   g_param_spec_object ("z-true",
                                                        NULL,
                                                        "Clusters true redshifts",
                                                        NCM_TYPE_VECTOR,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_LNM_OBS,
                                   g_param_spec_object ("lnM-obs",
                                                        NULL,
                                                        "Clusters mass observables",
                                                        NCM_TYPE_MATRIX,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_LNM_OBS_PARAMS,
                                   g_param_spec_object ("lnM-obs-params",
                                                        NULL,
                                                        "Clusters mass observables parameters",
                                                        NCM_TYPE_MATRIX,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_Z_OBS,
                                   g_param_spec_object ("z-obs",
                                                        NULL,
                                                        "Clusters redshift observables",
                                                        NCM_TYPE_MATRIX,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_Z_OBS_PARAMS,
                                   g_param_spec_object ("z-obs-params",
                                                        NULL,
                                                        "Clusters redshift observables parameters",
                                                        NCM_TYPE_MATRIX,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_SURVEY_AREA,
                                   g_param_spec_double ("area",
                                                        NULL,
                                                        "Cluster observation area",
                                                        0, G_MAXDOUBLE, 0,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_USE_TRUE,
                                   g_param_spec_boolean ("use-true",
                                                         NULL,
                                                         "If the true data must be used",
                                                         FALSE,
                                                         G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_BINNED,
                                   g_param_spec_boolean ("binned",
                                                         NULL,
                                                         "Whether use binned data",
                                                         FALSE,
                                                         G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_Z_NODES,
                                   g_param_spec_object ("z-nodes",
                                                        NULL,
                                                        "Clusters redshifts nodes for binning",
                                                        NCM_TYPE_VECTOR,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_LNM_NODES,
                                   g_param_spec_object ("lnM-nodes",
                                                        NULL,
                                                        "Clusters mass nodes for binning",
                                                        NCM_TYPE_VECTOR,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_FIDUCIAL,
                                   g_param_spec_boolean ("fiducial",
                                                         NULL,
                                                         "If it is fiducial data",
                                                         FALSE,
                                                         G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_RNG_NAME,
                                   g_param_spec_string ("rng-name",
                                                        NULL,
                                                        "Random number generator name",
                                                        NULL,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_SEED,
                                   g_param_spec_uint64 ("rng-seed",
                                                        NULL,
                                                        "Random number generator seed",
                                                        0, G_MAXUINT64, 0,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  data_class->name = "Cluster abundance unbinned";

  data_class->get_length = &_nc_data_cluster_ncount_get_length;
  data_class->begin      = &_nc_data_cluster_ncount_begin;
  data_class->prepare    = &_nc_data_cluster_ncount_prepare;
  data_class->resample   = &_nc_data_cluster_ncount_resample;
  data_class->m2lnL_val  = &_nc_data_cluster_ncount_m2lnL_val;
}

static guint 
_nc_data_cluster_ncount_get_length (NcmData *data) 
{ 
  return NC_DATA_CLUSTER_NCOUNT (data)->np; 
}

static void
_nc_data_cluster_ncount_begin (NcmData *data)
{
  gint signp = 0;
  NcDataClusterNCount *ncount = NC_DATA_CLUSTER_NCOUNT (data);   
  ncount->log_np_fac = lgamma_r (ncount->np + 1, &signp);
}

/**
 * nc_data_cluster_ncount_new:
 * @cad: a #NcClusterAbundance
 *
 * FIXME
 * 
 * Returns: FIXME
 */
NcDataClusterNCount *
nc_data_cluster_ncount_new (NcClusterAbundance *cad)
{
  NcDataClusterNCount *ncount = g_object_new (NC_TYPE_DATA_CLUSTER_NCOUNT,
                                              "cluster-abundance", cad,
                                              NULL);

  return ncount;
}

/**
 * nc_data_cluster_ncount_ref:
 * @ncount: a #NcDataClusterNCount
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
 * @ncount: a #NcDataClusterNCount
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
 * nc_data_cluster_ncount_clear:
 * @ncount: a #NcDataClusterNCount
 *
 * FIXME
 * 
 */
void
nc_data_cluster_ncount_clear (NcDataClusterNCount **ncount)
{
  g_clear_object (ncount);
}

/**
 * nc_data_cluster_ncount_set_n_z_obs:
 * @ncount: a #NcDataClusterNCount
 * @n_z_obs: FIXME
 *
 * FIXME
 * 
 */
void 
nc_data_cluster_ncount_set_n_z_obs (NcDataClusterNCount *ncount, guint n_z_obs)
{
  if (ncount->n_z_obs != n_z_obs)
  {
    if (ncount->z_obs)
      g_error ("nc_data_cluster_ncount_set_n_z_obs: cannot change the number of redshift observables in a non-empty data object.");
    ncount->n_z_obs = n_z_obs;
  }
}

/**
 * nc_data_cluster_ncount_set_n_z_obs_params:
 * @ncount: a #NcDataClusterNCount
 * @n_z_obs_params: FIXME
 *
 * FIXME
 * 
 */
void 
nc_data_cluster_ncount_set_n_z_obs_params (NcDataClusterNCount *ncount, guint n_z_obs_params)
{
  if (ncount->n_z_obs_params != n_z_obs_params)
  {
    if (ncount->z_obs_params)
      g_error ("nc_data_cluster_ncount_set_n_z_obs_params: cannot change the number of redshift observable parameters in a non-empty data object.");
    ncount->n_z_obs_params = n_z_obs_params;
  }
}

/**
 * nc_data_cluster_ncount_set_n_M_obs:
 * @ncount: a #NcDataClusterNCount
 * @n_M_obs: FIXME
 *
 * FIXME
 * 
 */
void 
nc_data_cluster_ncount_set_n_M_obs (NcDataClusterNCount *ncount, guint n_M_obs)
{
  if (ncount->n_M_obs != n_M_obs)
  {
    if (ncount->lnM_obs)
      g_error ("nc_data_cluster_ncount_set_n_M_obs: cannot change the number of mass observables in a non-empty data object.");
    ncount->n_M_obs = n_M_obs;
  }
}

/**
 * nc_data_cluster_ncount_set_n_M_obs_params:
 * @ncount: a #NcDataClusterNCount
 * @n_M_obs_params: FIXME
 *
 * FIXME
 * 
 */
void 
nc_data_cluster_ncount_set_n_M_obs_params (NcDataClusterNCount *ncount, guint n_M_obs_params)
{
  if (ncount->n_M_obs_params != n_M_obs_params)
  {
    if (ncount->lnM_obs_params)
      g_error ("nc_data_cluster_ncount_set_n_M_obs_params: cannot change the number of mass observable parameters in a non-empty data object.");
    ncount->n_M_obs_params = n_M_obs_params;
  }
}

/**
 * nc_data_cluster_ncount_set_lnM_true:
 * @ncount: a #NcDataClusterNCount
 * @v: a #NcmVector
 *
 * Sets the vector @v containing the values of the clusters true masses.
 * 
 */
void 
nc_data_cluster_ncount_set_lnM_true (NcDataClusterNCount *ncount, const NcmVector *v)
{
  if (ncount->np > 0)
  {
    if (ncm_vector_len (v) != ncount->np)
      g_error ("nc_data_cluster_ncount_set_lnM_true: incompatible vector, the data has %u clusters and the vector length is %u.",
               ncount->np, ncm_vector_len (v));
    if (ncount->lnM_true != NULL)
      ncm_vector_memcpy (ncount->lnM_true, v);
    else
      ncount->lnM_true = ncm_vector_dup (v);
  }
  else
  {
    ncount->np = ncm_vector_len (v);
    ncount->lnM_true = ncm_vector_dup (v);
  }
}

/**
 * nc_data_cluster_ncount_set_z_true:
 * @ncount: a #NcDataClusterNCount
 * @v: a #NcmVector 
 *
 * Sets the vector @v containing the true redshift values of the clusters.
 * 
 */
void 
nc_data_cluster_ncount_set_z_true (NcDataClusterNCount *ncount, const NcmVector *v)
{
  if (ncount->np > 0)
  {
    if (ncm_vector_len (v) != ncount->np)
      g_error ("nc_data_cluster_ncount_set_z_true: incompatible vector, the data has %u clusters and the vector length is %u.",
               ncount->np, ncm_vector_len (v));
    if (ncount->z_true != NULL)
      ncm_vector_memcpy (ncount->z_true, v);
    else
      ncount->z_true = ncm_vector_dup (v);
  }
  else
  {
    ncount->np = ncm_vector_len (v);
    ncount->z_true = ncm_vector_dup (v);
  }
}

/**
 * nc_data_cluster_ncount_set_lnM_obs:
 * @ncount: a #NcDataClusterNCount
 * @m: a #NcmMatrix
 *
 * Sets the matrix @m representing the cluster mass observables.
 * 
 */
void 
nc_data_cluster_ncount_set_lnM_obs (NcDataClusterNCount *ncount, const NcmMatrix *m)
{
  if (ncount->np > 0)
  {
    if (ncm_matrix_nrows (m) != ncount->np)
      g_error ("nc_data_cluster_ncount_set_lnM_obs: incompatible matrix, the data has %u clusters and the matrix has %u rows.",
               ncount->np, ncm_matrix_nrows (m));
    if (ncount->lnM_obs != NULL)
    {
      ncm_matrix_memcpy (ncount->lnM_obs, m);
      return;
    }
  }
  else
    ncount->np = ncm_matrix_nrows (m);

  if (ncount->n_M_obs != ncm_matrix_ncols (m))
    g_error ("nc_data_cluster_ncount_set_lnM_obs: incompatible matrix, NcmClusterMass object has %u points per observation and the matrix has %u cols.",
             ncount->n_M_obs, ncm_matrix_ncols (m));
  ncount->lnM_obs = ncm_matrix_dup (m);
}

/**
 * nc_data_cluster_ncount_set_lnM_obs_params:
 * @ncount: a #NcDataClusterNCount
 * @m: a #NcmMatrix
 *
 * Sets the matrix @m representing the cluster mass-observable parameters.
 * 
 */
void 
nc_data_cluster_ncount_set_lnM_obs_params (NcDataClusterNCount *ncount, const NcmMatrix *m)
{
  if (ncount->np > 0)
  {
    if (ncm_matrix_nrows (m) != ncount->np)
      g_error ("nc_data_cluster_ncount_set_lnM_obs_params: incompatible matrix, the data has %u clusters and the matrix has %u rows.",
               ncount->np, ncm_matrix_nrows (m));
    if (ncount->lnM_obs_params != NULL)
    {
      ncm_matrix_memcpy (ncount->lnM_obs_params, m);
      return;
    }
  }
  else
    ncount->np = ncm_matrix_nrows (m);

  if (ncount->n_M_obs != ncm_matrix_ncols (m))
    g_error ("nc_data_cluster_ncount_set_lnM_obs_params: incompatible matrix, NcmClusterMass object has %u parameters per observation and the matrix has %u cols.",
             ncount->n_M_obs_params, ncm_matrix_ncols (m));
  ncount->lnM_obs_params = ncm_matrix_dup (m);
}

/**
 * nc_data_cluster_ncount_set_z_obs:
 * @ncount: a #NcDataClusterNCount
 * @m: a #NcmMatrix
 *
 * Sets the matrix @m representing the cluster redshift observables.
 * 
 */
void 
nc_data_cluster_ncount_set_z_obs (NcDataClusterNCount *ncount, const NcmMatrix *m)
{
  if (ncount->np > 0)
  {
    if (ncm_matrix_nrows (m) != ncount->np)
      g_error ("nc_data_cluster_ncount_set_z_obs: incompatible matrix, the data has %u clusters and the matrix has %u rows.",
               ncount->np, ncm_matrix_nrows (m));
    if (ncount->z_obs != NULL)
    {
      ncm_matrix_memcpy (ncount->z_obs, m);
      return;
    }
  }
  else
    ncount->np = ncm_matrix_nrows (m);

  if (ncount->n_z_obs != ncm_matrix_ncols (m))
    g_error ("nc_data_cluster_ncount_set_z_obs: incompatible matrix, NcmClusterRedshift object has %u points per observation and the matrix has %u cols.",
             ncount->n_z_obs, ncm_matrix_ncols (m));
  ncount->z_obs = ncm_matrix_dup (m);
}

/**
 * nc_data_cluster_ncount_set_z_obs_params:
 * @ncount: a #NcDataClusterNCount
 * @m: a #NcmMatrix
 *
 * Sets the matrix @m representing the cluster redshift observable parameters.
 * 
 */
void 
nc_data_cluster_ncount_set_z_obs_params (NcDataClusterNCount *ncount, const NcmMatrix *m)
{
  if (ncount->np > 0)
  {
    if (ncm_matrix_nrows (m) != ncount->np)
      g_error ("nc_data_cluster_ncount_set_z_obs_params: incompatible matrix, the data has %u clusters and the matrix has %u rows.",
               ncount->np, ncm_matrix_nrows (m));
    if (ncount->z_obs_params != NULL)
    {
      ncm_matrix_memcpy (ncount->z_obs_params, m);
      return;
    }
  }
  else
    ncount->np = ncm_matrix_nrows (m);

  if (ncount->n_z_obs_params != ncm_matrix_ncols (m))
    g_error ("nc_data_cluster_ncount_set_lnM_obs: incompatible matrix, NcmClusterRedshift object has %u parameters per observation and the matrix has %u cols.",
             ncount->n_z_obs_params, ncm_matrix_ncols (m));
  ncount->z_obs_params = ncm_matrix_dup (m);
}

/**
 * nc_data_cluster_ncount_has_lnM_true:
 * @ncount: a #NcDataClusterNCount
 *
 * Returns: TRUE if it contains the lnM truth table.
 * 
 */
gboolean 
nc_data_cluster_ncount_has_lnM_true (NcDataClusterNCount *ncount)
{
  return ncount->lnM_true != NULL;
}

/**
 * nc_data_cluster_ncount_has_z_true:
 * @ncount: a #NcDataClusterNCount
 *
 * Returns: TRUE if it contains the redshift truth table.
 * 
 */
gboolean 
nc_data_cluster_ncount_has_z_true (NcDataClusterNCount *ncount)
{
  return ncount->z_true != NULL;
}

/**
 * nc_data_cluster_ncount_get_len:
 * @ncount: a #NcDataClusterNCount
 * 
 * Gets the total number of objects.
 * 
 * Returns: Number of objects in @ncount.
 */
guint 
nc_data_cluster_ncount_get_len (NcDataClusterNCount *ncount)
{
  return ncount->np;
}

/**
 * nc_data_cluster_ncount_lnM_obs_len:
 * @ncount: a #NcDataClusterNCount
 * 
 * Number of doubles to describe the observational data related to the mass
 * of each object.
 * 
 * Returns: Length of each row describing the mass proxy.
 */
guint 
nc_data_cluster_ncount_lnM_obs_len (NcDataClusterNCount *ncount)
{
  g_assert (ncount->lnM_obs != NULL);
  return ncm_matrix_row_len (ncount->lnM_obs);
}

/**
 * nc_data_cluster_ncount_lnM_obs_params_len:
 * @ncount: a #NcDataClusterNCount
 * 
 * Number of doubles to describe the observational data parameters related 
 * mass to each object.
 * 
 * Returns: Length of each row describing the mass proxy, it can be zero.
 */
guint 
nc_data_cluster_ncount_lnM_obs_params_len (NcDataClusterNCount *ncount)
{
  if (ncount->lnM_obs_params == NULL)
    return 0;
  else
    return ncm_matrix_row_len (ncount->lnM_obs_params);
}

/**
 * nc_data_cluster_ncount_z_obs_len:
 * @ncount: a #NcDataClusterNCount
 * 
 * Number of doubles to describe the observational data related to
 * the redshift of each object.
 * 
 * Returns: Length of each row describing the mass proxy.
 */
guint 
nc_data_cluster_ncount_z_obs_len (NcDataClusterNCount *ncount)
{
  g_assert (ncount->z_obs != NULL);
  return ncm_matrix_row_len (ncount->z_obs);
}

/**
 * nc_data_cluster_ncount_z_params_len:
 * @ncount: a #NcDataClusterNCount
 * 
 * Number of doubles to describe the observational data parameters related 
 * to the redshift of each object.
 * 
 * Returns: Length of each row describing the mass proxy, it can be zero.
 */
guint 
nc_data_cluster_ncount_z_obs_params_len (NcDataClusterNCount *ncount)
{
  if (ncount->z_obs_params == NULL)
    return 0;
  else
    return ncm_matrix_row_len (ncount->z_obs_params);
}

/**
 * nc_data_cluster_ncount_get_lnM_true:
 * @ncount: a #NcDataClusterNCount
 * 
 * Gets the vector containing the true values of the masses.
 * 
 * Returns: (transfer full): True masses #NcmVector.
 */
NcmVector *
nc_data_cluster_ncount_get_lnM_true (NcDataClusterNCount *ncount)
{
  if (ncount->lnM_true != NULL)
    return ncm_vector_ref (ncount->lnM_true);
  else
    return NULL;
}

/**
 * nc_data_cluster_ncount_get_z_true:
 * @ncount: a #NcDataClusterNCount
 * 
 * Gets the vector containing the true values of the redshifts.
 * 
 * Returns: (transfer full): True redshift #NcmVector.
 */
NcmVector *
nc_data_cluster_ncount_get_z_true (NcDataClusterNCount *ncount)
{
  if (ncount->z_true != NULL)
    return ncm_vector_ref (ncount->z_true);
  else
    return NULL;
}

/**
 * nc_data_cluster_ncount_get_lnM_obs:
 * @ncount: a #NcDataClusterNCount
 * 
 * Gets the matrix containing the mass observables.
 * 
 * Returns: (transfer full): Mass observable #NcmMatrix.
 */
NcmMatrix *
nc_data_cluster_ncount_get_lnM_obs (NcDataClusterNCount *ncount)
{
  if (ncount->lnM_obs != NULL)
    return ncm_matrix_ref (ncount->lnM_obs);
  else
    return NULL;
}

/**
 * nc_data_cluster_ncount_get_lnM_obs_params:
 * @ncount: a #NcDataClusterNCount
 * 
 * Gets the matrix containing the mass observables parameters.
 * 
 * Returns: (transfer full): Mass observable parameters #NcmMatrix.
 */
NcmMatrix *
nc_data_cluster_ncount_get_lnM_obs_params (NcDataClusterNCount *ncount)
{
  if (ncount->lnM_obs_params != NULL)
    return ncm_matrix_ref (ncount->lnM_obs_params);
  else
    return NULL;
}

/**
 * nc_data_cluster_ncount_get_z_obs:
 * @ncount: a #NcDataClusterNCount
 * 
 * Gets the matrix containing the redshift observables.
 * 
 * Returns: (transfer full): Redshift observable #NcmMatrix.
 */
NcmMatrix *
nc_data_cluster_ncount_get_z_obs (NcDataClusterNCount *ncount)
{
  if (ncount->z_obs != NULL)
    return ncm_matrix_ref (ncount->z_obs);
  else
    return NULL;
}

/**
 * nc_data_cluster_ncount_get_z_obs_params:
 * @ncount: a #NcDataClusterNCount
 * 
 * Gets the matrix containing the redshift observables parameters.
 * 
 * Returns: (transfer full): Redshift observable parameters #NcmMatrix.
 */
NcmMatrix *
nc_data_cluster_ncount_get_z_obs_params (NcDataClusterNCount *ncount)
{
  if (ncount->z_obs_params != NULL)
    return ncm_matrix_ref (ncount->z_obs_params);
  else
    return NULL;
}

/************************************************************************************************************
 * Unbinned number count data                                                                               *
 ************************************************************************************************************/

static void
_nc_data_cluster_ncount_model_init (NcDataClusterNCount *ncount)
{
  NcClusterAbundance *cad = ncount->cad;

  cad->purity        = ncount->purity;
  cad->sd_lnM        = ncount->sd_lnM;

  nc_halo_mass_function_set_area (cad->mfp, ncount->area_survey);
}

static void
_nc_data_cluster_ncount_prepare (NcmData *data, NcmMSet *mset)
{
  NcDataClusterNCount *ncount = NC_DATA_CLUSTER_NCOUNT (data);
  NcHICosmo *cosmo            = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));
  NcClusterRedshift *clusterz = NC_CLUSTER_REDSHIFT (ncm_mset_peek (mset, nc_cluster_redshift_id ()));
  NcClusterMass *clusterm     = NC_CLUSTER_MASS (ncm_mset_peek (mset, nc_cluster_mass_id ()));

  g_assert ((cosmo != NULL) && (clusterz != NULL) && (clusterm != NULL));

  g_assert_cmpuint (nc_cluster_mass_obs_len (clusterm), ==, ncount->n_z_obs);
  g_assert_cmpuint (nc_cluster_mass_obs_params_len (clusterm), ==, ncount->n_z_obs_params);
  
  g_assert_cmpuint (nc_cluster_redshift_obs_len (clusterz), ==, ncount->n_M_obs);
  g_assert_cmpuint (nc_cluster_redshift_obs_params_len (clusterz), ==, ncount->n_M_obs_params);
    
  nc_cluster_abundance_prepare_if_needed (ncount->cad, cosmo, clusterz, clusterm);
}

static gchar *
_nc_data_cluster_ncount_desc (NcDataClusterNCount *ncount, NcHICosmo *cosmo)
{
  NcClusterAbundance *cad = ncount->cad;
  return g_strdup_printf ("Cluster NCount resample %s. Generated %u from mean %10.5g (full). Resampled in range [%8.4f, %8.4f] [%1.8e, %1.8e] and area %8.4f degrees square",
                          ncount->binned ? "binned" : "unbinned",
                          ncount->np, cad->norma /*nc_cluster_abundance_n (cad, cosmo)*/, 
                          cad->zi, cad->zf, 
                          exp (cad->lnMi), exp (cad->lnMf), 
                          ncount->area_survey / gsl_pow_2 (M_PI / 180.0));
}

void _nc_data_cluster_ncount_bin_data (NcDataClusterNCount *ncount);

/**
 * _nc_data_cluster_ncount_resample:
 * @data: a #NcmData.
 * @mset: a #NcmMSet.
 *
 * This function generates random numbers which are used to obtain redshift
 * and mass (logarithm base e) values...
 *
 */
static void
_nc_data_cluster_ncount_resample (NcmData *data, NcmMSet *mset, NcmRNG *rng)
{
  NcDataClusterNCount *ncount = NC_DATA_CLUSTER_NCOUNT (data);
  NcClusterAbundance *cad     = ncount->cad;
  NcHICosmo *cosmo            = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));
  NcClusterRedshift *clusterz = NC_CLUSTER_REDSHIFT (ncm_mset_peek (mset, nc_cluster_redshift_id ()));
  NcClusterMass *clusterm     = NC_CLUSTER_MASS (ncm_mset_peek (mset, nc_cluster_mass_id ()));

  guint z_obs_len          = nc_cluster_redshift_obs_len (clusterz);
  guint z_obs_params_len   = nc_cluster_redshift_obs_params_len (clusterz);
  guint lnM_obs_len        = nc_cluster_mass_obs_len (clusterm);
  guint lnM_obs_params_len = nc_cluster_mass_obs_params_len (clusterm);

  GArray *lnM_true_array       = NULL;
  GArray *z_true_array         = NULL;
  GArray *z_obs_array          = NULL;
  GArray *z_obs_params_array   = NULL;
  GArray *lnM_obs_array        = NULL;
  GArray *lnM_obs_params_array = NULL;
  guint total_np;
  guint i;

  gdouble *zi_obs = g_new (gdouble, z_obs_len);
  gdouble *zi_obs_params = z_obs_params_len > 0 ? g_new (gdouble, z_obs_params_len) : NULL;
  gdouble *lnMi_obs = g_new (gdouble, lnM_obs_len);
  gdouble *lnMi_obs_params = lnM_obs_params_len > 0 ? g_new (gdouble, lnM_obs_params_len) : NULL;

  ncm_rng_lock (rng);
  total_np = gsl_ran_poisson (rng->r, cad->norma);
  ncm_rng_unlock (rng);

  nc_data_cluster_ncount_set_n_z_obs (ncount, z_obs_len);
  nc_data_cluster_ncount_set_n_z_obs_params (ncount, z_obs_params_len);

  nc_data_cluster_ncount_set_n_M_obs (ncount, lnM_obs_len);
  nc_data_cluster_ncount_set_n_M_obs_params (ncount, lnM_obs_params_len);
  
  if (total_np == 0)
  {
    ncount->np = 0;
    g_free (zi_obs);
    g_free (lnMi_obs);
    if (z_obs_params_len > 0)
      g_free (zi_obs_params);
    if (lnM_obs_params_len > 0)
      g_free (lnMi_obs_params);

    if (ncount->binned)
    {
      _nc_data_cluster_ncount_bin_data (ncount);
    }

    ncm_data_take_desc (data, _nc_data_cluster_ncount_desc (ncount, cosmo));
    return;
  }

  lnM_true_array = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), total_np);
  z_true_array   = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), total_np);

  z_obs_array = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), total_np * z_obs_len);
  if (z_obs_params_len > 0)
    z_obs_params_array = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), total_np * z_obs_params_len);

  lnM_obs_array = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), total_np * lnM_obs_len);
  if (lnM_obs_params_len > 0)
    lnM_obs_params_array = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), total_np * lnM_obs_params_len);

  nc_cluster_abundance_prepare_inv_dNdz (cad, NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ())),
                                         cad->lnMi);

  for (i = 0; i < total_np; i++)
  {
    ncm_rng_lock (rng);
    {
      const gdouble u1 = _nc_cad_inv_dNdz_convergence_f (gsl_rng_uniform_pos (rng->r), cad->z_epsilon);
      const gdouble u2 = _nc_cad_inv_dNdz_convergence_f (gsl_rng_uniform_pos (rng->r), cad->lnM_epsilon);
      const gdouble z_true = ncm_spline_eval (cad->inv_z, u1);
      const gdouble lnM_true = ncm_spline2d_eval (cad->inv_lnM_z, u2, z_true);
      ncm_rng_unlock (rng);

      if ( nc_cluster_redshift_resample (clusterz, lnM_true, z_true, zi_obs, zi_obs_params, rng) &&
          nc_cluster_mass_resample (clusterm, cosmo, lnM_true, z_true, lnMi_obs, lnMi_obs_params, rng) )
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
    }
  }

  if (z_obs_array->len == 0 || lnM_obs_array->len == 0)
  {
    ncount->np = 0;
    g_free (zi_obs);
    g_free (lnMi_obs);
    if (z_obs_params_len > 0)
      g_free (zi_obs_params);
    if (lnM_obs_params_len > 0)
      g_free (lnMi_obs_params);

    ncm_vector_clear (&ncount->lnM_true);
    g_array_unref (lnM_true_array);
    ncm_vector_clear (&ncount->z_true);
    g_array_unref (z_true_array);
    ncm_matrix_clear (&ncount->z_obs);
    g_array_unref (z_obs_array);
    ncm_matrix_clear (&ncount->lnM_obs);
    g_array_unref (lnM_obs_array);

    if (ncount->binned)
    {
      _nc_data_cluster_ncount_bin_data (ncount);
    }
    
    ncm_data_take_desc (data, _nc_data_cluster_ncount_desc (ncount, cosmo));
    return;
  }

  ncm_vector_clear (&ncount->lnM_true);
  ncount->lnM_true = ncm_vector_new_array (lnM_true_array);
  g_array_unref (lnM_true_array);

  ncm_vector_clear (&ncount->z_true);
  ncount->z_true = ncm_vector_new_array (z_true_array);
  g_array_unref (z_true_array);

  ncm_matrix_clear (&ncount->z_obs);
  ncount->z_obs = ncm_matrix_new_array (z_obs_array, z_obs_len);
  g_array_unref (z_obs_array);

  ncm_matrix_clear (&ncount->lnM_obs);
  ncount->lnM_obs = ncm_matrix_new_array (lnM_obs_array, lnM_obs_len);
  g_array_unref (lnM_obs_array);

  if (z_obs_params_len > 0)
  {
    ncm_matrix_clear (&ncount->z_obs_params);
    ncount->z_obs_params = ncm_matrix_new_array (z_obs_params_array, z_obs_params_len);
    g_array_unref (z_obs_params_array);
  }

  if (lnM_obs_params_len > 0)
  {
    ncm_matrix_clear (&ncount->lnM_obs_params);
    ncount->lnM_obs_params = ncm_matrix_new_array (lnM_obs_params_array, lnM_obs_params_len);
    g_array_unref (lnM_obs_params_array);
  }

  ncount->np = ncm_matrix_nrows (ncount->z_obs);

  /* printf ("Generated %u, Expected %10.5g\n", ncount->np, nc_cluster_abundance_n (cad, cosmo)); */
  ncm_data_take_desc (data, _nc_data_cluster_ncount_desc (ncount, cosmo));

  g_free (zi_obs);
  g_free (zi_obs_params);
  g_free (lnMi_obs);
  g_free (lnMi_obs_params);

  if (ncount->binned)
  {
    _nc_data_cluster_ncount_bin_data (ncount);
  }
}

typedef struct
{
  NcClusterAbundance *cad;
  NcDataClusterNCount *ncount;
  NcClusterRedshift *clusterz;
  NcClusterMass *clusterm;
  NcHICosmo *cosmo;
} _Evald2N;

static void
_eval_z_p_lnM_p_d2n (glong i, glong f, gpointer data)
{
  _Evald2N *evald2n = (_Evald2N *) data;
  glong n;
  
  for (n = i; n < f; n++)
  {
    gdouble *lnMn_obs = ncm_matrix_ptr (evald2n->ncount->lnM_obs, n, 0);
    gdouble *lnMn_obs_params = evald2n->ncount->lnM_obs_params != NULL ? ncm_matrix_ptr (evald2n->ncount->lnM_obs_params, n, 0) : NULL;
    gdouble *zn_obs = ncm_matrix_ptr (evald2n->ncount->z_obs, n, 0);
    gdouble *zn_obs_params = evald2n->ncount->z_obs_params != NULL ? ncm_matrix_ptr (evald2n->ncount->z_obs_params, n, 0) : NULL;
    const gdouble mlnLn = -log (nc_cluster_abundance_z_p_lnM_p_d2n (evald2n->cad, evald2n->cosmo, evald2n->clusterz, evald2n->clusterm, lnMn_obs, lnMn_obs_params, zn_obs, zn_obs_params));
    g_array_index (evald2n->ncount->m2lnL_a, gdouble, n) = mlnLn;
  }
}

static void
_eval_z_p_d2n (glong i, glong f, gpointer data)
{
  _Evald2N *evald2n = (_Evald2N *) data;
  glong n;

  for (n = i; n < f; n++)
  {
    const gdouble lnMn = ncm_vector_get (evald2n->ncount->lnM_true, n);
    gdouble *zn_obs = ncm_matrix_ptr (evald2n->ncount->z_obs, n, 0);
    gdouble *zn_obs_params = evald2n->ncount->z_obs_params != NULL ? ncm_matrix_ptr (evald2n->ncount->z_obs_params, n, 0) : NULL;
    const gdouble mlnLn = -log (nc_cluster_abundance_z_p_d2n (evald2n->cad, evald2n->cosmo, evald2n->clusterz, evald2n->clusterm, lnMn, zn_obs, zn_obs_params));
    g_array_index (evald2n->ncount->m2lnL_a, gdouble, n) = mlnLn;
  }
}

static void
_eval_lnM_p_d2n (glong i, glong f, gpointer data)
{
  _Evald2N *evald2n = (_Evald2N *) data;
  glong n;

  for (n = i; n < f; n++)
  {
    const gdouble zn = ncm_vector_get (evald2n->ncount->z_true, n);
    gdouble *lnMn_obs = ncm_matrix_ptr (evald2n->ncount->lnM_obs, n, 0);
    gdouble *lnMn_obs_params = evald2n->ncount->lnM_obs_params != NULL ? ncm_matrix_ptr (evald2n->ncount->lnM_obs_params, n, 0) : NULL;
    const gdouble mlnLn = -log (nc_cluster_abundance_lnM_p_d2n (evald2n->cad, evald2n->cosmo, evald2n->clusterz, evald2n->clusterm, lnMn_obs, lnMn_obs_params, zn));
    g_array_index (evald2n->ncount->m2lnL_a, gdouble, n) = mlnLn;
  }
}

static void
_eval_d2n (glong i, glong f, gpointer data)
{
  _Evald2N *evald2n = (_Evald2N *) data;
  glong n;

  for (n = i; n < f; n++)
  {
    const gdouble lnMn = ncm_vector_get (evald2n->ncount->lnM_true, n);
    const gdouble zn = ncm_vector_get (evald2n->ncount->z_true, n);
    const gdouble mlnLn = -log (nc_cluster_abundance_d2n (evald2n->cad, evald2n->cosmo, evald2n->clusterz, evald2n->clusterm, lnMn, zn));
    g_array_index (evald2n->ncount->m2lnL_a, gdouble, n) = mlnLn;
  }
}

static void
_eval_intp_d2n (glong i, glong f, gpointer data)
{
  _Evald2N *evald2n = (_Evald2N *) data;
  glong n;

  for (n = i; n < f; n++)
  {
    const gdouble lnMn = ncm_vector_get (evald2n->ncount->lnM_true, n);
    const gdouble zn = ncm_vector_get (evald2n->ncount->z_true, n);
    const gdouble mlnLn = -log (nc_cluster_abundance_intp_d2n (evald2n->cad, evald2n->cosmo, evald2n->clusterz, evald2n->clusterm, lnMn, zn));
    g_array_index (evald2n->ncount->m2lnL_a, gdouble, n) = mlnLn;
  }
}

static void
_nc_data_cluster_ncount_m2lnL_val (NcmData *data, NcmMSet *mset, gdouble *m2lnL)
{
  NcDataClusterNCount *ncount = NC_DATA_CLUSTER_NCOUNT (data);
  NcClusterAbundance *cad     = ncount->cad;
  NcHICosmo *cosmo            = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));
  NcClusterRedshift *clusterz = NC_CLUSTER_REDSHIFT (ncm_mset_peek (mset, nc_cluster_redshift_id ()));
  NcClusterMass *clusterm     = NC_CLUSTER_MASS (ncm_mset_peek (mset, nc_cluster_mass_id ()));
  
  *m2lnL = 0.0;

  if (ncount->binned)
    g_error ("_nc_data_cluster_ncount_m2lnL_val: don't support binned likelihood yet.");

  if (ncount->np == 0)
  {
    const gdouble n_th = nc_cluster_abundance_n (cad, cosmo, clusterz, clusterm);    
    *m2lnL = 2.0 * (ncount->log_np_fac + n_th);
    return;
  }

  g_array_set_size (ncount->m2lnL_a, ncount->np);

  if (ncount->use_true_data)
  {
    _Evald2N evald2n = {cad, ncount, clusterz, clusterm, cosmo};
    g_assert (ncount->z_true);
    g_assert (ncount->lnM_true);
    ncm_func_eval_threaded_loop_full (&_eval_intp_d2n, 0, ncount->np, &evald2n);
  }
  else
  {
    gboolean z_p   = ncm_model_check_impl_opt (NCM_MODEL (clusterz), NC_CLUSTER_REDSHIFT_P);
    gboolean lnM_p = ncm_model_check_impl_opt (NCM_MODEL (clusterm), NC_CLUSTER_MASS_P);
    
    if (z_p && lnM_p)
    {
      _Evald2N evald2n = {cad, ncount, clusterz, clusterm, cosmo};
      ncm_func_eval_threaded_loop_full (&_eval_z_p_lnM_p_d2n, 0, ncount->np, &evald2n);
    }
    else if (z_p && !lnM_p)
    {
      g_assert (ncount->lnM_true);
      _Evald2N evald2n = {cad, ncount, clusterz, clusterm, cosmo};
      ncm_func_eval_threaded_loop_full (&_eval_z_p_d2n, 0, ncount->np, &evald2n);
    }
    else if (!z_p && lnM_p)
    {
      g_assert (ncount->z_true);
      _Evald2N evald2n = {cad, ncount, clusterz, clusterm, cosmo};
      ncm_func_eval_threaded_loop_full (&_eval_lnM_p_d2n, 0, ncount->np, &evald2n);
    }
    else
    {
      g_assert (ncount->z_true);
      g_assert (ncount->lnM_true);
      _Evald2N evald2n = {cad, ncount, clusterz, clusterm, cosmo};
      ncm_func_eval_threaded_loop_full (&_eval_d2n, 0, ncount->np, &evald2n);

    }
  }

  {
    const gdouble n_th = nc_cluster_abundance_n (cad, cosmo, clusterz, clusterm);
    guint i;
    for (i = 0 ; i < ncount->np; i++) 
    { 
      *m2lnL += g_array_index (ncount->m2lnL_a, gdouble, i); 
    }
    *m2lnL += (ncount->log_np_fac + n_th);
  }

  *m2lnL *= 2.0;
}

/**
 * nc_data_cluster_ncount_true_data:
 * @ncount: a #NcDataClusterNCount
 * @use_true_data: FIXME
 *
 * FIXME
 *
 */
void
nc_data_cluster_ncount_true_data (NcDataClusterNCount *ncount, gboolean use_true_data)
{
  if (use_true_data)
  {
    g_assert (ncount->lnM_true != NULL);
    g_assert (ncount->z_true != NULL);
  }
  ncount->use_true_data = use_true_data;
}

/**
 * nc_data_cluster_ncount_using_true_data:
 * @ncount: a #NcDataClusterNCount
 *
 * Returns: if it is using true data.
 */
gboolean 
nc_data_cluster_ncount_using_true_data (NcDataClusterNCount *ncount)
{
  return ncount->use_true_data;
}

static void _nc_data_cluster_ncount_model_init (NcDataClusterNCount *ncount);

/**
 * nc_data_cluster_ncount_init_from_sampling:
 * @ncount: a #NcDataClusterNCount
 * @mset: a #NcmMSet
 * @area_survey: area in units of square degrees
 * @rng: a #NcmRNG
 *
 * FIXME
 *
 */
void
nc_data_cluster_ncount_init_from_sampling (NcDataClusterNCount *ncount, NcmMSet *mset, gdouble area_survey, NcmRNG *rng)
{
  ncount->area_survey = area_survey;
  _nc_data_cluster_ncount_model_init (ncount);

  ncm_data_resample (NCM_DATA (ncount), mset, rng);
  ncm_data_set_init (NCM_DATA (ncount), TRUE);
}

/**
 * nc_data_cluster_ncount_print:
 * @ncount: a #NcDataClusterNCount
 * @cosmo: a #NcHICosmo
 * @out: FIXME
 * @header: FIXME
 *
 * FIXME
 *
 */
void
nc_data_cluster_ncount_print (NcDataClusterNCount *ncount, NcHICosmo *cosmo, FILE *out, gchar *header)
{
  NcClusterAbundance *cad = ncount->cad;
  gint i, j;
  gint nbins_M = 100;
  gint nbins_z = 20;
  gsl_vector *lnM_nodes = gsl_vector_alloc (nbins_M);
  gsl_vector *z_nodes = gsl_vector_alloc (nbins_z);

  g_assert (ncount->np > 0);
  
  for (i = 0; i < nbins_M; i++)
  {
    gdouble lnM = cad->lnMi + (cad->lnMf - cad->lnMi) / (nbins_M - 1.0) * i;
    gsl_vector_set (lnM_nodes, i, lnM);

    printf ("lnM = %5.5g Me = %5.5g M10 = %5.5g\n", gsl_vector_get (lnM_nodes, i), exp(gsl_vector_get (lnM_nodes, i)), pow(10, gsl_vector_get (lnM_nodes, i)));
  }

  for (j = 0; j < nbins_z; j++)
  {
    gdouble z = cad->zi + (cad->zf - cad->zi) / (nbins_z - 1.0) * j;
    gsl_vector_set (z_nodes, j, z);

    printf ("z = %5.5g\n", gsl_vector_get (z_nodes, j));
  }

  if (header != NULL)
    fprintf (out, "# %s\n# ", header);
  else
    fprintf (out, "# ");

  g_assert_not_reached ();
  gsl_histogram2d *h = NULL;
  /*gsl_histogram2d *h = nc_data_cluster_ncount_hist_lnM_z (ncount, lnM_nodes, z_nodes);*/

  fprintf (out, "# z M N/(logM * V) (catalog) dn/dlog10M (theory) Nmi(catalog)(abundance in bins of redshift and mass) Nmi(theory)\n");
  for (j = 0; j < nbins_z - 1; j++)
  {
    gdouble zm, dz, zl, zu, V;
    gsl_histogram2d_get_yrange (h, j, &zl, &zu);
    zm = (zu + zl) / 2.0;
    dz = (zu - zl);
    V = nc_halo_mass_function_dv_dzdomega (cad->mfp, cosmo, zm) * ncount->area_survey * dz;
    for (i = 0; i < nbins_M; i++)
    {
      gdouble ln_ml, ln_mu, Mm, lnMm, log_mu, log_ml;
      gdouble Nmi = gsl_histogram2d_get (h, i, j);
      gdouble dndlog10M, ca_M;
      gsl_histogram2d_get_xrange (h, i, &ln_ml, &ln_mu);
      lnMm = (ln_ml + ln_mu) / 2.0;
      Mm = exp (lnMm);
      log_mu = log10 (exp(ln_mu));
      log_ml = log10 (exp(ln_ml));
      dndlog10M = M_LN10 * nc_halo_mass_function_dn_dlnM (cad->mfp, cosmo, lnMm, zm);
      ca_M = (log_mu - log_ml) * V * dndlog10M;

      //printf ("log-mu = %5.5g log-ml = %5.5g\n", log_mu, log_ml);
      fprintf (out, "% 6.6g % 6.6e % 6.6g % 6.6g % 6.6g % 6.6g\n", zm, Mm, Nmi / ((log_mu - log_ml) * V), dndlog10M, Nmi, ca_M);
    }
    fprintf (out, "\n\n");
  }

  gsl_vector_free (lnM_nodes);
  gsl_vector_free (z_nodes);
  gsl_histogram2d_free (h);
}

void
_nc_data_cluster_ncount_bin_data (NcDataClusterNCount *ncount)
{
  gsl_histogram2d_reset (ncount->z_lnM);

  if (ncount->np == 0)
    return;

  if (ncount->use_true_data)
  {
    guint i;

    for (i = 0; i < ncount->np; i++)
    {
      const gdouble lnM = ncm_vector_get (ncount->lnM_true, i);
      const gdouble z   = ncm_vector_get (ncount->z_true, i);
      gsl_histogram2d_increment (ncount->z_lnM, z, lnM);
    }    
  }
  else
  {
    guint i;
    
    for (i = 0; i < ncount->np; i++)
    {
      const gdouble lnM = ncm_matrix_get (ncount->lnM_obs, i, 0);
      const gdouble z   = ncm_matrix_get (ncount->z_obs, i, 0);
      gsl_histogram2d_increment (ncount->z_lnM, z, lnM);
    }
  }
}

static void
_nc_data_cluster_ncount_bin_alloc (NcDataClusterNCount *ncount, guint z_bins, guint lnM_bins)
{
  g_assert_cmpint (z_bins, >=, 1);
  g_assert_cmpint (lnM_bins, >=, 1);

  if (ncount->z_lnM != NULL)
  {
    if (ncount->z_lnM->nx != z_bins || ncount->z_lnM->ny != lnM_bins)
    {
      g_clear_pointer (&ncount->z_lnM, gsl_histogram2d_free);
      ncount->z_lnM = gsl_histogram2d_alloc (z_bins, lnM_bins);
    }
  }
  else
    ncount->z_lnM = gsl_histogram2d_alloc (z_bins, lnM_bins);

}

/**
 * nc_data_cluster_ncount_set_bin_by_nodes:
 * @ncount: a #NcDataClusterNCount
 * @z_nodes: a #NcmVector
 * @lnM_nodes: a #NcmVector
 *
 * FIXME
 *
 */
void
nc_data_cluster_ncount_set_bin_by_nodes (NcDataClusterNCount *ncount, NcmVector *z_nodes, NcmVector *lnM_nodes)
{
  guint z_bins   = ncm_vector_len (z_nodes) - 1;
  guint lnM_bins = ncm_vector_len (lnM_nodes) - 1;
  _nc_data_cluster_ncount_bin_alloc (ncount, z_bins, lnM_bins);

  g_assert_cmpint (ncm_vector_stride (z_nodes), ==, 1);
  g_assert_cmpint (ncm_vector_stride (lnM_nodes), ==, 1);
  
  gsl_histogram2d_set_ranges (ncount->z_lnM, 
                              ncm_vector_ptr (z_nodes, 0), z_bins + 1, 
                              ncm_vector_ptr (lnM_nodes, 0), lnM_bins + 1);

  _nc_data_cluster_ncount_bin_data (ncount);
  ncount->binned = TRUE;
/*
  ncm_vector_log_vals (z_nodes,   "# z   ", "% 20.15g");
  ncm_vector_log_vals (lnM_nodes, "# lnM ", "% 20.15g");
*/  
  if (ncount->z_nodes != z_nodes)
  {
    ncm_vector_clear (&ncount->z_nodes);
    ncount->z_nodes = ncm_vector_ref (z_nodes);
  }

  if (ncount->lnM_nodes != lnM_nodes)
  {
    ncm_vector_clear (&ncount->lnM_nodes);
    ncount->lnM_nodes = ncm_vector_ref (lnM_nodes);
  }
}

/**
 * nc_data_cluster_ncount_set_bin_by_minmax:
 * @ncount: a #NcDataClusterNCount.
 * @z_nbins: number of bins in z.
 * @lnM_nbins: number of bins in $\ln(M)$. 
 *
 * Creates a uniform binning with @z_nbins and @lnM_nbins.
 *
 */
void
nc_data_cluster_ncount_set_bin_by_minmax (NcDataClusterNCount *ncount, guint z_nbins, guint lnM_nbins)
{
  gdouble z_min = 0.0, z_max = 0.0, lnM_min = 0.0, lnM_max = 0.0;
  _nc_data_cluster_ncount_bin_alloc (ncount, z_nbins, lnM_nbins);
  if (ncount->use_true_data)
  {
    gsl_vector_minmax (ncm_vector_gsl (ncount->z_true), &z_min, &z_max);
    gsl_vector_minmax (ncm_vector_gsl (ncount->lnM_true), &lnM_min, &lnM_max);
  }
  else
  {
    gsl_matrix_minmax (ncm_matrix_gsl (ncount->z_obs), &z_min, &z_max);
    gsl_matrix_minmax (ncm_matrix_gsl (ncount->lnM_obs), &lnM_min, &lnM_max);
  }

  z_min = 0.0;
  z_max = z_max * 1.5;

  lnM_min = lnM_min + (lnM_min > 0.0 ? -0.5 * lnM_min : 0.5 * lnM_min);
  lnM_max = lnM_max + (lnM_max > 0.0 ? 0.5 * lnM_min : - 0.5 * lnM_min);

  gsl_histogram2d_set_ranges_uniform (ncount->z_lnM, z_min, z_max, lnM_min, lnM_max);
  _nc_data_cluster_ncount_bin_data (ncount);
  ncount->binned = TRUE;

  ncm_vector_clear (&ncount->z_nodes);
  ncm_vector_clear (&ncount->lnM_nodes);
  ncount->z_nodes   = ncm_vector_new_data_dup (ncount->z_lnM->xrange, ncount->z_lnM->nx + 1, 1);
  ncount->lnM_nodes = ncm_vector_new_data_dup (ncount->z_lnM->yrange, ncount->z_lnM->ny + 1, 1);
/*
  ncm_vector_log_vals (ncount->z_nodes,   "minmax-z  ", "%8.4g");
  ncm_vector_log_vals (ncount->lnM_nodes, "minmax-lnM", "%8.4g");
*/
}

/**
 * nc_data_cluster_ncount_set_bin_by_quantile:
 * @ncount: a #NcDataClusterNCount.
 * @z_quantiles: FIXME
 * @lnM_quantiles: FIXME
 *
 * FIXME
 *
 */
void 
nc_data_cluster_ncount_set_bin_by_quantile (NcDataClusterNCount *ncount, NcmVector *z_quantiles, NcmVector *lnM_quantiles)
{
  guint z_bins   = ncm_vector_len (z_quantiles) + 1;
  guint lnM_bins = ncm_vector_len (lnM_quantiles) + 1;
  NcmVector *z_nodes = ncm_vector_new (z_bins + 1);
  NcmVector *lnM_nodes = ncm_vector_new (lnM_bins + 1);
  guint i;

  if ((ncm_vector_get (z_quantiles, 0) == 0.0) || 
      (ncm_vector_get (lnM_quantiles, 0) == 0.0) ||
      (ncm_vector_get (z_quantiles, ncm_vector_len (z_quantiles) - 1) == 1.0) ||
      (ncm_vector_get (lnM_quantiles, ncm_vector_len (lnM_quantiles) - 1) == 1.0))
  {
    g_error ("nc_data_cluster_ncount_set_bin_by_quantile: quantiles must be > 0.0 and < 1.0");
  }


  {
    gdouble *z_data, *lnM_data;
    guint z_stride, lnM_stride;
    guint z_len, lnM_len;
    NcmVector *z_dup = NULL;
    NcmVector *lnM_dup = NULL;
    gdouble z_min = 0.0, z_max = 0.0, lnM_min = 0.0, lnM_max = 0.0;
    
    if (ncount->use_true_data)
    {
      z_dup   = ncm_vector_dup (ncount->z_true);
      lnM_dup = ncm_vector_dup (ncount->lnM_true);
    }
    else
    {
      NcmVector *col = ncm_matrix_get_col (ncount->z_obs, 0);

      z_dup   = ncm_vector_dup (col);
      ncm_vector_free (col);

      col = ncm_matrix_get_col (ncount->lnM_obs, 0);
      lnM_dup = ncm_vector_dup (col);
      
      ncm_vector_free (col);
    }

    z_data     = ncm_vector_ptr (z_dup, 0);
    lnM_data   = ncm_vector_ptr (lnM_dup, 0);
    z_stride   = ncm_vector_stride (z_dup);
    lnM_stride = ncm_vector_stride (lnM_dup);
    z_len      = ncm_vector_len (z_dup);
    lnM_len    = ncm_vector_len (lnM_dup);

    gsl_sort (z_data, z_stride, z_len);
    gsl_sort (lnM_data, lnM_stride, lnM_len);
    
    z_min = 0.0;
    z_max = z_data[z_stride * (z_len - 1)] * 1.5;

    lnM_min = lnM_data[0];
    lnM_max = lnM_data[lnM_stride * (lnM_len - 1)];
    
    lnM_min = lnM_min + (lnM_min > 0.0 ? -0.5 * lnM_min : 0.5 * lnM_min);
    lnM_max = lnM_max + (lnM_max > 0.0 ? 0.5 * lnM_min : - 0.5 * lnM_min);    
    
    ncm_vector_set (z_nodes,   0, z_min);
    ncm_vector_set (lnM_nodes, 0, lnM_min);

    for (i = 1; i < z_bins; i++)
    {
      const gdouble z_node = gsl_stats_quantile_from_sorted_data (z_data, z_stride, z_len,
                                                                  ncm_vector_get (z_quantiles, i - 1)); 
      ncm_vector_set (z_nodes, i, z_node);
    }
    for (i = 1; i < lnM_bins; i++)
    {
      const gdouble lnM_node = gsl_stats_quantile_from_sorted_data (lnM_data, lnM_stride, lnM_len,
                                                                    ncm_vector_get (lnM_quantiles, i - 1)); 
      ncm_vector_set (lnM_nodes, i, lnM_node);
    }

    ncm_vector_set (z_nodes,   z_bins,   z_max);
    ncm_vector_set (lnM_nodes, lnM_bins, lnM_max);
    ncm_vector_free (z_dup);
    ncm_vector_free (lnM_dup);
  }

  nc_data_cluster_ncount_set_bin_by_nodes (ncount, z_nodes, lnM_nodes);

  ncm_vector_free (z_nodes);
  ncm_vector_free (lnM_nodes);
}

/**
 * nc_data_cluster_ncount_set_binned:
 * @ncount: a #NcDataClusterNCount.
 * @on: FIXME
 *
 * FIXME
 *
 */
void 
nc_data_cluster_ncount_set_binned (NcDataClusterNCount *ncount, gboolean on)
{
  if (on)
  {
    if (ncount->binned)
      return;
    else
    {
      if (ncount->z_nodes == NULL || ncount->lnM_nodes == NULL)
      {
        g_error ("nc_data_cluster_ncount_set_binned: cannot turn on, no nodes were defined.");
      }
      else
      {
        nc_data_cluster_ncount_set_bin_by_nodes (ncount, ncount->z_nodes, ncount->lnM_nodes);
      }
    }
  }
  else
    ncount->binned = FALSE;  
}

#ifdef NUMCOSMO_HAVE_CFITSIO

/**
 * nc_data_cluster_ncount_catalog_save:
 * @ncount: a #NcDataClusterNCount.
 * @filename: name of the file
 * @overwrite: FIXME
 *
 * FIXME
 *
 */
void
nc_data_cluster_ncount_catalog_save (NcDataClusterNCount *ncount, gchar *filename, gboolean overwrite)
{
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

  g_assert (ncount->np > 0);
  
  g_ptr_array_set_free_func (tform_array, g_free);

  g_ptr_array_add (ttype_array, "Z_OBS");
  g_ptr_array_add (tform_array, g_strdup_printf ("%dD", ncount->n_z_obs));
  g_ptr_array_add (tunit_array, "REDSHIFT OBS");

  g_ptr_array_add (ttype_array, "LNM_OBS");
  g_ptr_array_add (tform_array, g_strdup_printf ("%dD", ncount->n_M_obs));
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

  if (ncount->n_z_obs_params > 0)
  {
    g_ptr_array_add (ttype_array, "Z_OBS_PARAMS");
    g_ptr_array_add (tform_array, g_strdup_printf ("%dD", ncount->n_z_obs_params));
    g_ptr_array_add (tunit_array, "REDSHIFT OBS PARAMS");
  }

  if (ncount->n_M_obs_params > 0)
  {
    g_ptr_array_add (ttype_array, "LNM_OBS_PARAMS");
    g_ptr_array_add (tform_array, g_strdup_printf ("%dD", ncount->n_M_obs_params));
    g_ptr_array_add (tunit_array, "LNM OBS PARAMS");
  }

  tfields = ttype_array->len;

  /* initialize status before calling fitsio routines */
  status = 0;

  if (overwrite && g_file_test (filename, G_FILE_TEST_EXISTS))
    g_unlink (filename);

  /* create new FITS file */
  fits_create_file (&fptr, filename, &status);
  NCM_FITS_ERROR (status);

  /* append a new empty binary table onto the FITS file */
  fits_create_tbl (fptr, BINARY_TBL, ncount->np, tfields, (gchar **)ttype_array->pdata, (gchar **)tform_array->pdata,
                       (gchar **)tunit_array->pdata, extname, &status);
  NCM_FITS_ERROR (status);

  {
    gdouble sarea_d = ncount->area_survey / gsl_pow_2 (M_PI / 180.0);
    fits_write_key(fptr, TDOUBLE, "AREA", &sarea_d, "Survey area in degree square", &status);
    NCM_FITS_ERROR (status);
  }

  {
    guint colnum = 1;
    fits_write_col (fptr, TDOUBLE, colnum, 1, 1, ncount->np, ncm_matrix_ptr (ncount->z_obs, 0, 0), &status);
    NCM_FITS_ERROR (status);

    colnum++;
    fits_write_col (fptr, TDOUBLE, colnum, 1, 1, ncount->np, ncm_matrix_ptr (ncount->lnM_obs, 0, 0), &status);
    NCM_FITS_ERROR (status);

    if (ncount->z_true != NULL)
    {
      colnum++;
      fits_write_col (fptr, TDOUBLE, colnum, 1, 1, ncount->np, ncm_vector_ptr (ncount->z_true, 0), &status);
      NCM_FITS_ERROR (status);
    }

    if (ncount->lnM_true != NULL)
    {
      colnum++;
      fits_write_col (fptr, TDOUBLE, colnum, 1, 1, ncount->np, ncm_vector_ptr (ncount->lnM_true, 0), &status);
      NCM_FITS_ERROR (status);
    }

    if (ncount->n_z_obs_params > 0)
    {
      colnum++;
      fits_write_col (fptr, TDOUBLE, colnum, 1, 1, ncount->np, ncm_matrix_ptr (ncount->z_obs_params, 0, 0), &status);
      NCM_FITS_ERROR (status);
    }

    if (ncount->n_M_obs_params > 0)
    {
      colnum++;
      fits_write_col (fptr, TDOUBLE, colnum, 1, 1, ncount->np, ncm_matrix_ptr (ncount->lnM_obs_params, 0, 0), &status);
      NCM_FITS_ERROR (status);
    }

  }

  fits_close_file(fptr, &status);
  NCM_FITS_ERROR (status);

  g_ptr_array_unref (ttype_array);
  g_ptr_array_unref (tform_array);
  g_ptr_array_unref (tunit_array);

  return;
}

/**
 * nc_data_cluster_ncount_catalog_load:
 * @ncount: a #NcDataClusterNCount.
 * @filename: name of the file
 *
 * FIXME
 *
 */
void
nc_data_cluster_ncount_catalog_load (NcDataClusterNCount *ncount, gchar *filename)
{
  gint status, hdutype;
  gchar comment[FLEN_COMMENT];
  fitsfile *fptr;

  status = 0;

  if (filename == NULL)
    g_error ("nc_cluster_abundance_catalog_load: null filename");

  fits_open_file (&fptr, filename, READONLY, &status);
  NCM_FITS_ERROR (status);

  fits_movabs_hdu (fptr, 2, &hdutype, &status);
  NCM_FITS_ERROR (status);

  if (hdutype != BINARY_TBL)
    g_error ("%s (%d): NcDataClusterNCount catalog is not binary!", __FILE__, __LINE__);
  
  {
    glong nrows;
    fits_get_num_rows (fptr, &nrows, &status);
    NCM_FITS_ERROR (status);
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

    nc_data_cluster_ncount_set_n_z_obs (ncount, z_obs_rp);
    nc_data_cluster_ncount_set_n_M_obs (ncount, lnM_obs_rp);

    ncm_matrix_clear (&ncount->z_obs);
    ncount->z_obs = ncm_matrix_new (ncount->np, z_obs_rp);

    ncm_matrix_clear (&ncount->lnM_obs);
    ncount->lnM_obs = ncm_matrix_new (ncount->np, lnM_obs_rp);

    fits_read_col (fptr, TDOUBLE, z_obs_i, 1, 1, ncount->np * z_obs_rp, NULL, ncm_matrix_ptr (ncount->z_obs, 0, 0), NULL, &status);
    NCM_FITS_ERROR (status);

    fits_read_col (fptr, TDOUBLE, lnM_obs_i, 1, 1, ncount->np * lnM_obs_rp, NULL, ncm_matrix_ptr (ncount->lnM_obs, 0, 0), NULL, &status);
    NCM_FITS_ERROR (status);
  }

  {
    gint z_obs_params_i, lnM_obs_params_i;
    gint z_obs_params_tc, lnM_obs_params_tc;
    glong z_obs_params_rp, lnM_obs_params_rp;
    glong z_obs_params_w, lnM_obs_params_w;

    if (fits_get_colnum (fptr, CASESEN, "Z_OBS_PARAMS", &z_obs_params_i, &status))
    {
      nc_data_cluster_ncount_set_n_z_obs_params (ncount, 0);
    }
    else
    {
      if (fits_get_coltype (fptr, z_obs_params_i, &z_obs_params_tc, &z_obs_params_rp, &z_obs_params_w, &status))
        g_error ("Column Z_OBS_PARAMS info not found, invalid fits file.");

      nc_data_cluster_ncount_set_n_z_obs_params (ncount, z_obs_params_rp);

      ncm_matrix_clear (&ncount->z_obs_params);
      ncount->z_obs_params = ncm_matrix_new (ncount->np, z_obs_params_rp);

      fits_read_col (fptr, TDOUBLE, z_obs_params_i, 1, 1, ncount->np * z_obs_params_rp, NULL, ncm_matrix_ptr (ncount->z_obs_params, 0, 0), NULL, &status);
      NCM_FITS_ERROR (status);
    }

    if (fits_get_colnum (fptr, CASESEN, "LNM_OBS_PARAMS", &lnM_obs_params_i, &status))
    {
      nc_data_cluster_ncount_set_n_M_obs_params (ncount, 0);
    }
    else
    {
      if (fits_get_coltype (fptr, lnM_obs_params_i, &lnM_obs_params_tc, &lnM_obs_params_rp, &lnM_obs_params_w, &status))
        g_error ("Column LNM_OBS_PARAMS info not found, invalid fits file.");

      nc_data_cluster_ncount_set_n_M_obs_params (ncount, lnM_obs_params_rp);

      ncm_matrix_clear (&ncount->lnM_obs_params);
      ncount->lnM_obs_params = ncm_matrix_new (ncount->np, lnM_obs_params_rp);

      fits_read_col (fptr, TDOUBLE, lnM_obs_params_i, 1, 1, ncount->np * lnM_obs_params_rp, NULL, ncm_matrix_ptr (ncount->lnM_obs_params, 0, 0), NULL, &status);
      NCM_FITS_ERROR (status);
    }
  }

  {
    gint z_true_i, lnM_true_i;
    
    ncm_vector_clear (&ncount->z_true);
    if (!fits_get_colnum (fptr, CASESEN, "Z_TRUE", &z_true_i, &status))
    {
      ncount->z_true = ncm_vector_new (ncount->np);
      fits_read_col (fptr, TDOUBLE, z_true_i, 1, 1, ncount->np, NULL, ncm_vector_ptr (ncount->z_true, 0), NULL, &status);
      NCM_FITS_ERROR (status);
    }
    else
      status = 0;

    ncm_vector_clear (&ncount->lnM_true);
    if (!fits_get_colnum (fptr, CASESEN, "LNM_TRUE", &lnM_true_i, &status))
    {
      ncount->lnM_true = ncm_vector_new (ncount->np);
      fits_read_col (fptr, TDOUBLE, lnM_true_i, 1, 1, ncount->np, NULL, ncm_vector_ptr (ncount->lnM_true, 0), NULL, &status);
      NCM_FITS_ERROR (status);
    }
    else
      status = 0;
  }

  fits_close_file (fptr, &status);
  NCM_FITS_ERROR (status);

  _nc_data_cluster_ncount_model_init (ncount);

  ncm_data_set_init (NCM_DATA (ncount), TRUE);
  
  return;
}

#endif /* NUMCOSMO_HAVE_CFITSIO */

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

#include "math/ncm_func_eval.h"
#include "math/ncm_serialize.h"
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
  PROP_CLUSTERZ,
  PROP_CLUSTERM,
  PROP_LNM_TRUE,
  PROP_Z_TRUE,
  PROP_Z_OBS,
  PROP_Z_OBS_PARAMS,
  PROP_LNM_OBS,
  PROP_LNM_OBS_PARAMS,
  PROP_SURVEY_AREA,
  PROP_USE_TRUE,
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
      nc_cluster_abundance_clear (&ncount->cad);
      ncount->cad = g_value_dup_object (value);
      break;
    case PROP_CLUSTERZ:
      nc_data_cluster_ncount_set_redshift (ncount, g_value_get_object (value));
      break;
    case PROP_CLUSTERM:
      nc_data_cluster_ncount_set_mass (ncount, g_value_get_object (value));
      break;
    case PROP_LNM_TRUE:
    {
      GVariant *var = g_value_get_variant (value);
      const NcmVector *v = ncm_vector_const_new_variant (var);
      nc_data_cluster_ncount_set_lnM_true (ncount, v);
      ncm_vector_const_free (v);
      break;
    }
    case PROP_Z_TRUE:
    {
      GVariant *var = g_value_get_variant (value);
      const NcmVector *v = ncm_vector_const_new_variant (var);
      nc_data_cluster_ncount_set_z_true (ncount, v);
      ncm_vector_const_free (v);
      break;
    }
    case PROP_Z_OBS:
    {
      GVariant *var = g_value_get_variant (value);
      const NcmMatrix *m = ncm_matrix_const_new_variant (var);
      nc_data_cluster_ncount_set_z_obs (ncount, m);
      ncm_matrix_const_free (m);
      break;
    }
    case PROP_Z_OBS_PARAMS:
    {
      GVariant *var = g_value_get_variant (value);
      const NcmMatrix *m = ncm_matrix_const_new_variant (var);
      nc_data_cluster_ncount_set_z_obs_params (ncount, m);
      ncm_matrix_const_free (m);
      break;
    }
    case PROP_LNM_OBS:
    {
      GVariant *var = g_value_get_variant (value);
      const NcmMatrix *m = ncm_matrix_const_new_variant (var);
      nc_data_cluster_ncount_set_lnM_obs (ncount, m);
      ncm_matrix_const_free (m);
      break;
    }
    case PROP_LNM_OBS_PARAMS:
    {
      GVariant *var = g_value_get_variant (value);
      const NcmMatrix *m = ncm_matrix_const_new_variant (var);
      nc_data_cluster_ncount_set_lnM_obs_params (ncount, m);
      ncm_matrix_const_free (m);
      break;
    }
    case PROP_SURVEY_AREA:
      ncount->area_survey = g_value_get_double (value);
      break;
    case PROP_USE_TRUE:
      ncount->use_true_data = g_value_get_boolean (value);
      break;
    case PROP_FIDUCIAL:
      ncount->fiducial = g_value_get_boolean (value);
      break;
    case PROP_SEED:
      ncount->seed = g_value_get_ulong (value);
      break;
    case PROP_RNG_NAME:
      g_clear_pointer (&ncount->rnd_name, &g_free);
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
    case PROP_CLUSTERZ:
      g_value_set_object (value, ncount->z);
      break;
    case PROP_CLUSTERM:
      g_value_set_object (value, ncount->m);
      break;
    case PROP_LNM_TRUE:
    {
      GVariant *var = ncm_vector_peek_variant (ncount->lnM_true); 
      g_value_take_variant (value, var);
      break;
    }
    case PROP_Z_TRUE:
    {
      GVariant *var = ncm_vector_peek_variant (ncount->z_true); 
      g_value_take_variant (value, var);
      break;
    }
    case PROP_Z_OBS:
    {
      GVariant *var = ncm_matrix_peek_variant (ncount->z_obs); 
      g_value_take_variant (value, var);
      break;
    }
    case PROP_Z_OBS_PARAMS:
    {
      GVariant *var = ncm_matrix_peek_variant (ncount->z_obs_params); 
      g_value_take_variant (value, var);
      break;
    }
    case PROP_LNM_OBS:
    {
      GVariant *var = ncm_matrix_peek_variant (ncount->lnM_obs); 
      g_value_take_variant (value, var);
      break;
    }
    case PROP_LNM_OBS_PARAMS:
    {
      GVariant *var = ncm_matrix_peek_variant (ncount->lnM_obs_params); 
      g_value_take_variant (value, var);
      break;
    }
    case PROP_SURVEY_AREA:
      g_value_set_double (value, ncount->area_survey);
      break;
    case PROP_USE_TRUE:
      g_value_set_boolean (value, ncount->use_true_data);
      break;
    case PROP_FIDUCIAL:
      g_value_set_boolean (value, ncount->fiducial);
      break;
    case PROP_SEED:
      g_value_set_ulong (value, ncount->seed);
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
                                   PROP_CLUSTERM,
                                   g_param_spec_object ("mass",
                                                        NULL,
                                                        "Cluster mass observable",
                                                        NC_TYPE_CLUSTER_MASS,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_CLUSTERZ,
                                   g_param_spec_object ("redshift",
                                                        NULL,
                                                        "Cluster redshift observable",
                                                        NC_TYPE_CLUSTER_REDSHIFT,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_LNM_TRUE,
                                   g_param_spec_variant ("lnM-true",
                                                         NULL,
                                                         "Clusters true masses",
                                                         G_VARIANT_TYPE ("ad"), NULL,
                                                         G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));  
  g_object_class_install_property (object_class,
                                   PROP_Z_TRUE,
                                   g_param_spec_variant ("z-true",
                                                         NULL,
                                                         "Clusters true redshifts",
                                                         G_VARIANT_TYPE ("ad"), NULL,
                                                         G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_LNM_OBS,
                                   g_param_spec_variant ("lnM-obs",
                                                         NULL,
                                                         "Clusters mass observables",
                                                         G_VARIANT_TYPE ("aad"), NULL,
                                                         G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_LNM_OBS_PARAMS,
                                   g_param_spec_variant ("lnM-obs-params",
                                                         NULL,
                                                         "Clusters mass observables parameters",
                                                         G_VARIANT_TYPE ("aad"), NULL,
                                                         G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_Z_OBS,
                                   g_param_spec_variant ("z-obs",
                                                         NULL,
                                                         "Clusters redshift observables",
                                                         G_VARIANT_TYPE ("aad"), NULL,
                                                         G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_Z_OBS_PARAMS,
                                   g_param_spec_variant ("z-obs-params",
                                                         NULL,
                                                         "Clusters redshift observables parameters",
                                                         G_VARIANT_TYPE ("aad"), NULL,
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
                                   g_param_spec_ulong ("rng-seed",
                                                        NULL,
                                                        "Random number generator seed",
                                                        0, G_MAXULONG, 0,
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
                                              "cluster-abundance", cad,
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
 * nc_data_cluster_ncount_clear:
 * @ncount: FIXME
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
 * nc_data_cluster_ncount_set_mass:
 * @ncount: a #NcDataClusterNCount.
 * @m: the #NcClusterMass object.
 *
 * Sets the mass observable object.
 * 
 */
void 
nc_data_cluster_ncount_set_mass (NcDataClusterNCount *ncount, NcClusterMass *m)
{
  if (ncount->lnM_obs != NULL)
    g_assert_cmpuint (ncm_matrix_ncols (ncount->lnM_obs), ==, nc_cluster_mass_obs_len (m));
  if (ncount->lnM_obs_params != NULL)
    g_assert_cmpuint (ncm_matrix_ncols (ncount->lnM_obs_params), ==, nc_cluster_mass_obs_params_len (m));
  
  if (ncount->m == NULL)
    ncount->m = nc_cluster_mass_ref (m);
  else
  {
    nc_cluster_mass_clear (&ncount->m);
    ncount->m = nc_cluster_mass_ref (m);
  }
}

/**
 * nc_data_cluster_ncount_set_redshift:
 * @ncount: a #NcDataClusterNCount.
 * @z: the #NcClusterRedshift object.
 *
 * Sets the redshift observable object.
 * 
 */
void 
nc_data_cluster_ncount_set_redshift (NcDataClusterNCount *ncount, NcClusterRedshift *z)
{
  if (ncount->lnM_obs != NULL)
    g_assert_cmpuint (ncm_matrix_ncols (ncount->z_obs), ==, nc_cluster_redshift_obs_len (z));
  if (ncount->lnM_obs_params != NULL)
    g_assert_cmpuint (ncm_matrix_ncols (ncount->z_obs_params), ==, nc_cluster_redshift_obs_params_len (z));

  if (ncount->z == NULL)
    ncount->z = nc_cluster_redshift_ref (z);
  else
  {
    nc_cluster_redshift_clear (&ncount->z);
    ncount->z = nc_cluster_redshift_ref (z);
  }
}

/**
 * nc_data_cluster_ncount_set_lnM_true:
 * @ncount: a #NcDataClusterNCount.
 * @v: the value of the masses.
 *
 * Sets the vector representing the clusters true masses.
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
 * @ncount: a #NcDataClusterNCount.
 * @v: the value of the redshifts.
 *
 * Sets the vector representing the clusters true redshifts.
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
 * @ncount: a #NcDataClusterNCount.
 * @m: the value of the mass observables.
 *
 * Sets the matrix representing the clusters mass observables.
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

  if (ncount->m != NULL && nc_cluster_mass_obs_len (ncount->m) != ncm_matrix_ncols (m))
    g_error ("nc_data_cluster_ncount_set_lnM_obs: incompatible matrix, NcmClusterMass object has %u points per observation and the matrix has %u cols.",
             nc_cluster_mass_obs_len (ncount->m), ncm_matrix_ncols (m));
  ncount->lnM_obs = ncm_matrix_dup (m);
}

/**
 * nc_data_cluster_ncount_set_lnM_obs_params:
 * @ncount: a #NcDataClusterNCount.
 * @m: the mass observables parameters.
 *
 * Sets the matrix representing the clusters mass observables parameters.
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

  if (ncount->m != NULL && nc_cluster_mass_obs_params_len (ncount->m) != ncm_matrix_ncols (m))
    g_error ("nc_data_cluster_ncount_set_lnM_obs_params: incompatible matrix, NcmClusterMass object has %u parameters per observation and the matrix has %u cols.",
             nc_cluster_mass_obs_params_len (ncount->m), ncm_matrix_ncols (m));
  ncount->lnM_obs_params = ncm_matrix_dup (m);
}

/**
 * nc_data_cluster_ncount_set_z_obs:
 * @ncount: a #NcDataClusterNCount.
 * @m: the value of the redshift observables.
 *
 * Sets the matrix representing the clusters redshift observables.
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

  if (ncount->m != NULL && nc_cluster_redshift_obs_len (ncount->z) != ncm_matrix_ncols (m))
    g_error ("nc_data_cluster_ncount_set_z_obs: incompatible matrix, NcmClusterRedshift object has %u points per observation and the matrix has %u cols.",
             nc_cluster_redshift_obs_len (ncount->z), ncm_matrix_ncols (m));
  ncount->z_obs = ncm_matrix_dup (m);
}

/**
 * nc_data_cluster_ncount_set_z_obs_params:
 * @ncount: a #NcDataClusterNCount.
 * @m: the redshift observables parameters.
 *
 * Sets the matrix representing the clusters redshift observables parameters.
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

  if (ncount->m != NULL && nc_cluster_redshift_obs_params_len (ncount->z) != ncm_matrix_ncols (m))
    g_error ("nc_data_cluster_ncount_set_lnM_obs: incompatible matrix, NcmClusterRedshift object has %u parameters per observation and the matrix has %u cols.",
             nc_cluster_redshift_obs_params_len (ncount->z), ncm_matrix_ncols (m));
  ncount->z_obs_params = ncm_matrix_dup (m);
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
  NCM_UNUSED (data);
  NCM_UNUSED (filename);
  g_assert_not_reached ();
}

static void
_nc_data_cluster_ncount_binned_f (NcmMSet *mset, gpointer obj, const gdouble *x, gdouble *f)
{
  //NcHICosmo *cosmo = NC_HICOSMO (ncm_mset_peek (mset, NC_HICOSMO_ID));
  //NcClusterAbundance *cad = NC_CLUSTER_ABUNDANCE (obj);
  //f[0] = nc_cluster_abundance_N_val (cad, model, cad->lnMi, cad->lnMf, cad->zi, x[0]);
  NCM_UNUSED (mset);
  NCM_UNUSED (obj);
  NCM_UNUSED (x);
  NCM_UNUSED (f);
  g_assert_not_reached ();
  return;
}

/**
 * nc_data_cluster_ncount_binned_create_func: (skip)
 * @cad: a #NcClusterAbundance
 *
 * FIXME
 *
 * Returns: FIXME
 */
NcmMSetFunc *
nc_data_cluster_ncount_binned_create_func (NcClusterAbundance *cad)
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
  NCM_UNUSED (cad);
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
  NcHICosmo *cosmo = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));

  if (ncm_model_ctrl_update (ncount->cad->ctrl, NCM_MODEL (cosmo)))
  {
    nc_cluster_abundance_prepare (ncount->cad, cosmo);
  }
}

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
  NcClusterAbundance *cad = ncount->cad;
  NcHICosmo *cosmo = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));
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
  guint i;

  gdouble *zi_obs = g_new (gdouble, z_obs_len);
  gdouble *zi_obs_params = z_obs_params_len > 0 ? g_new (gdouble, z_obs_params_len) : NULL;
  gdouble *lnMi_obs = g_new (gdouble, lnM_obs_len);
  gdouble *lnMi_obs_params = lnM_obs_params_len > 0 ? g_new (gdouble, lnM_obs_params_len) : NULL;

  ncm_rng_lock (rng);
  total_np = gsl_ran_poisson (rng->r, cad->norma);
  ncm_rng_unlock (rng);

  lnM_true_array = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), total_np);
  z_true_array = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), total_np);

  z_obs_array = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), total_np * z_obs_len);
  if (z_obs_params_len > 0)
    z_obs_params_array = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), total_np * z_obs_params_len);

  lnM_obs_array = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), total_np * lnM_obs_len);
  if (lnM_obs_params_len > 0)
    lnM_obs_params_array = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), total_np * lnM_obs_params_len);
  
  nc_cluster_abundance_prepare_inv_dNdz (cad, NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ())));

  for (i = 0; i < total_np; i++)
  {
    ncm_rng_lock (rng);
    {
      const gdouble u1 = _nc_cad_inv_dNdz_convergence_f (gsl_rng_uniform_pos (rng->r), cad->z_epsilon);
      const gdouble u2 = _nc_cad_inv_dNdz_convergence_f (gsl_rng_uniform_pos (rng->r), cad->lnM_epsilon);
      const gdouble z_true = ncm_spline_eval (cad->inv_z, u1);
      const gdouble lnM_true = ncm_spline2d_eval (cad->inv_lnM_z, u2, z_true);
      ncm_rng_unlock (rng);

      if ( nc_cluster_redshift_resample (cad->z, lnM_true, z_true, zi_obs, zi_obs_params, rng) &&
          nc_cluster_mass_resample (cad->m, cosmo, lnM_true, z_true, lnMi_obs, lnMi_obs_params, rng) )
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
  
  ncm_vector_clear (&ncount->lnM_true);
  ncount->lnM_true = ncm_vector_new_array (lnM_true_array);
  g_array_unref (lnM_true_array);

  ncm_vector_clear (&ncount->z_true);
  ncount->z_true = ncm_vector_new_array (z_true_array);
  g_array_unref (z_true_array);

  if (z_obs_array->len == 0 || lnM_obs_array->len == 0)
    g_error ("_nc_data_cluster_ncount_resample: error generating sample zero objects generated (%u, %u)", z_obs_array->len, lnM_obs_array->len);

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

  ncm_data_take_desc (data, 
                      g_strdup_printf ("Cluster NCount resample unbinned. Generated %u from mean %10.5g. Resampled in range [%8.4f, %8.4f] [%1.8e, %1.8e] and area %8.4f degrees square", 
                                       ncount->np, nc_cluster_abundance_n (cad, cosmo), 
                                       cad->zi, cad->zf, 
                                       exp (cad->lnMi), exp (cad->lnMf), 
                                       ncount->area_survey / gsl_pow_2 (M_PI / 180.0)));
  ncm_rng_unlock (rng);
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
  NcHICosmo *cosmo = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));
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

  {
    const gdouble n_th = nc_cluster_abundance_n (cad, cosmo);    
    *m2lnL += (ncount->log_np_fac + n_th);
  }

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
 * @rng: a #NcmRNG.
 *
 * FIXME
 *
 */
void
nc_data_cluster_ncount_init_from_sampling (NcmData *data, NcmMSet *mset, NcClusterRedshift *clusterz, NcClusterMass *clusterm, gdouble area_survey, NcmRNG *rng)
{
  NcDataClusterNCount *ncount = NC_DATA_CLUSTER_NCOUNT (data);

  nc_data_cluster_ncount_set_redshift (ncount, clusterz);
  nc_data_cluster_ncount_set_mass (ncount, clusterm);
  ncount->area_survey = area_survey;

  _nc_data_cluster_ncount_model_init (ncount);

  ncm_data_set_init (data, TRUE);

  ncm_data_resample (data, mset, rng);
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
  guint i;

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

  ncm_data_poisson_init_from_histogram (data_cpoisson, hist);

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
  guint i;

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

/**
 * nc_data_cluster_ncount_print:
 * @data: FIXME
 * @cosmo: a NcHICosmo
 * @out: FIXME
 * @header: FIXME
 *
 * FIXME
 *
 */
void
nc_data_cluster_ncount_print (NcmData *data, NcHICosmo *cosmo, FILE *out, gchar *header)
{
  NcDataClusterNCount *ncount = NC_DATA_CLUSTER_NCOUNT (data);
  NcClusterAbundance *cad = ncount->cad;
  gint i, j;
  gint nbins_M = 100;
  gint nbins_z = 20;
  gsl_vector *lnM_nodes = gsl_vector_alloc (nbins_M);
  gsl_vector *z_nodes = gsl_vector_alloc (nbins_z);

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

  gsl_histogram2d *h = nc_data_cluster_ncount_hist_lnM_z (data, lnM_nodes, z_nodes);


  fprintf (out, "# z M N/(logM * V) (catalog) dn/dlog10M (theory) Nmi(catalog)(abundance in bins of redshift and mass) Nmi(theory)\n");
  for (j = 0; j < nbins_z - 1; j++)
  {
    gdouble zm, dz, zl, zu, V;
    gsl_histogram2d_get_yrange (h, j, &zl, &zu);
    zm = (zu + zl) / 2.0;
    dz = (zu - zl);
    V = nc_mass_function_dv_dzdomega (cad->mfp, cosmo, zm) * cad->mfp->area_survey * dz;
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
      dndlog10M = M_LN10 * nc_mass_function_dn_dlnm (cad->mfp, cosmo, lnMm, zm);
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

#ifdef NUMCOSMO_HAVE_CFITSIO

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
  fits_create_file (&fptr, filename, &status);
  NCM_FITS_ERROR (status);

  /* append a new empty binary table onto the FITS file */
  fits_create_tbl (fptr, BINARY_TBL, ncount->np, tfields, (gchar **)ttype_array->pdata, (gchar **)tform_array->pdata,
                       (gchar **)tunit_array->pdata, extname, &status);
  NCM_FITS_ERROR (status);

  {
    gchar *z_ser = ncm_serialize_global_to_string (G_OBJECT (ncount->z), FALSE);
    gchar *lnM_ser = ncm_serialize_global_to_string (G_OBJECT (ncount->m), FALSE);

    fits_write_key_longstr (fptr, "Z_OBJ", z_ser, "Serialized redshift object", &status);
    NCM_FITS_ERROR (status);
    fits_write_key_longstr (fptr, "LNM_OBJ", lnM_ser, "Serialized mass object", &status);
    NCM_FITS_ERROR (status);

    g_free (z_ser);
    g_free (lnM_ser);
  }

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

    if (z_obs_params_len > 0)
    {
      colnum++;
      fits_write_col (fptr, TDOUBLE, colnum, 1, 1, ncount->np, ncm_matrix_ptr (ncount->z_obs_params, 0, 0), &status);
      NCM_FITS_ERROR (status);
    }

    if (lnM_obs_params_len > 0)
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

  fits_open_file (&fptr, filename, READONLY, &status);
  NCM_FITS_ERROR (status);

  fits_movabs_hdu (fptr, 2, &hdutype, &status);
  NCM_FITS_ERROR (status);

  if (hdutype != BINARY_TBL)
    g_error ("%s (%d): NcDataClusterNCount catalog is not binary!\n", __FILE__, __LINE__);

  if (ncount->z != NULL)
    nc_cluster_redshift_free (ncount->z);
  if (ncount->m != NULL)
    nc_cluster_mass_free (ncount->m);

  {
    gchar *z_ser = NULL;
    gchar *lnM_ser = NULL;

    fits_read_key_longstr (fptr, "Z_OBJ", &z_ser, comment, &status);
    NCM_FITS_ERROR (status);
    fits_read_key_longstr (fptr, "LNM_OBJ", &lnM_ser, comment, &status);
    NCM_FITS_ERROR (status);

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

    fits_free_memory (z_ser, &status);
    NCM_FITS_ERROR (status);
    fits_free_memory (lnM_ser, &status);
    NCM_FITS_ERROR (status);
  }

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

    if (nc_cluster_redshift_obs_len (ncount->z) != z_obs_rp)
      g_error ("NcClusterRedshift object has observables length %d but fits has %ld.",
               nc_cluster_redshift_obs_len (ncount->z), z_obs_rp);

    if (nc_cluster_mass_obs_len (ncount->m) != lnM_obs_rp)
      g_error ("NcClusterMass object has observables length %d but fits has %ld.",
               nc_cluster_mass_obs_len (ncount->m), lnM_obs_rp);

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

      ncm_matrix_clear (&ncount->z_obs_params);
      ncount->z_obs_params = ncm_matrix_new (ncount->np, z_obs_params_rp);

      fits_read_col (fptr, TDOUBLE, z_obs_params_i, 1, 1, ncount->np * z_obs_params_rp, NULL, ncm_matrix_ptr (ncount->z_obs_params, 0, 0), NULL, &status);
      NCM_FITS_ERROR (status);
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

  ncm_data_set_init (data, TRUE);
  
  return;
}

#endif /* NUMCOSMO_HAVE_CFITSIO */

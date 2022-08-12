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
#include "nc_hireion.h"

#include "math/ncm_func_eval.h"
#include "math/ncm_serialize.h"
#include "math/ncm_obj_array.h"
#include "math/ncm_cfg.h"

#include "misc/cubature.h"

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
  PROP_MASS_TYPE,
  PROP_REDSHIFT_TYPE,
  PROP_LNM_TRUE,
  PROP_Z_TRUE,
  PROP_Z_OBS,
  PROP_Z_OBS_PARAMS,
  PROP_LNM_OBS,
  PROP_LNM_OBS_PARAMS,
  PROP_SURVEY_AREA,
  PROP_USE_TRUE,
  PROP_BINNED,
  PROP_Z_OBS_BINS,
  PROP_LNM_OBS_BINS,
  PROP_BIN_COUNT,
  PROP_FIDUCIAL,
  PROP_SEED,
  PROP_RNG_NAME,
  PROP_SIZE,
};

struct _NcDataClusterNCountPrivate
{
  NcClusterAbundance *cad;
  GType mass_type;
  GType redshift_type;
  NcmVector *lnM_true;
  NcmVector *z_true;
  NcmMatrix *z_obs;
  NcmMatrix *z_obs_params;
  NcmMatrix *lnM_obs;
  NcmMatrix *lnM_obs_params;
  GArray *m2lnL_a;
  GArray *m2lnL_err_a;
  GArray *z_order;
  GArray *lnM_order;
  GArray *p_z;
  GArray *p_lnM;
  GArray *d2n;
  gdouble area_survey;
  guint np;
  guint z_obs_len;
  guint z_obs_params_len;
  guint lnM_obs_len;
  guint lnM_obs_params_len;
  gdouble log_np_fac;
  gboolean use_true_data;
  gboolean binned;
  gsl_histogram2d *purity;
  gsl_histogram2d *sd_lnM;
  NcmObjArray *z_obs_bins;
  NcmObjArray *lnM_obs_bins;
  NcmVector *bin_count;
  gboolean fiducial;
  guint64 seed;
  gchar *rnd_name;
};

G_DEFINE_TYPE_WITH_PRIVATE (NcDataClusterNCount, nc_data_cluster_ncount, NCM_TYPE_DATA);

static void
nc_data_cluster_ncount_init (NcDataClusterNCount *ncount)
{
  NcDataClusterNCountPrivate * const self = ncount->priv = nc_data_cluster_ncount_get_instance_private (ncount);

  self->cad              = NULL;
  self->mass_type        = G_TYPE_INVALID;
  self->redshift_type    = G_TYPE_INVALID;
  self->z_obs_len        = 0;
  self->z_obs_params_len = 0;
  self->lnM_obs_len        = 0;
  self->lnM_obs_params_len = 0;
  self->lnM_true         = NULL;
  self->z_true           = NULL;
  self->z_obs            = NULL;
  self->z_obs_params     = NULL;
  self->lnM_obs          = NULL;
  self->lnM_obs_params   = NULL;
  self->m2lnL_a          = g_array_new (FALSE, FALSE, sizeof (gdouble));
  self->m2lnL_err_a      = g_array_new (FALSE, FALSE, sizeof (gdouble));
  self->z_order          = g_array_new (FALSE, FALSE, sizeof (size_t));
  self->lnM_order        = g_array_new (FALSE, FALSE, sizeof (size_t));
  self->p_z              = g_array_new (FALSE, FALSE, sizeof (gdouble));
  self->p_lnM            = g_array_new (FALSE, FALSE, sizeof (gdouble));
  self->d2n              = g_array_new (FALSE, FALSE, sizeof (gdouble));
  self->area_survey      = 0.0;
  self->np               = 0.0;
  self->log_np_fac       = 0.0;
  self->use_true_data    = FALSE;
  self->binned           = FALSE;
  self->purity           = NULL;
  self->sd_lnM           = NULL;
  self->z_obs_bins       = ncm_obj_array_new ();
  self->lnM_obs_bins     = ncm_obj_array_new ();
  self->bin_count        = NULL;
  self->fiducial         = FALSE;
  self->seed             = 0;
  self->rnd_name         = NULL;
}

static void
nc_data_cluster_ncount_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcDataClusterNCount *ncount = NC_DATA_CLUSTER_NCOUNT (object);
  NcDataClusterNCountPrivate * const self = ncount->priv;
  g_return_if_fail (NC_IS_DATA_CLUSTER_NCOUNT (object));

  switch (prop_id)
  {
    case PROP_CAD:
      nc_cluster_abundance_clear (&self->cad);
      self->cad = g_value_dup_object (value);
      break;
    case PROP_MASS_TYPE:
    {
      NcClusterMassClass *clusterm_class;

      self->mass_type = g_type_from_name (g_value_get_string (value));
      if (self->mass_type == G_TYPE_INVALID)
        g_error ("nc_data_cluster_ncount_set_property: GType `%s' unregistered or invalid.", g_value_get_string (value));

      clusterm_class         = g_type_class_ref (self->mass_type);
      self->lnM_obs_len        = nc_cluster_mass_class_obs_len (clusterm_class);
      self->lnM_obs_params_len = nc_cluster_mass_class_obs_params_len (clusterm_class);

      g_assert (self->lnM_obs_len > 0);

      g_type_class_unref (clusterm_class);
      break;
    }
    case PROP_REDSHIFT_TYPE:
    {
      NcClusterRedshiftClass *clusterz_class;

      self->redshift_type = g_type_from_name (g_value_get_string (value));
      if (self->redshift_type == G_TYPE_INVALID)
        g_error ("nc_data_cluster_ncount_set_property: GType `%s' unregistered or invalid.", g_value_get_string (value));

      clusterz_class         = g_type_class_ref (self->redshift_type);
      self->z_obs_len        = nc_cluster_redshift_class_obs_len (clusterz_class);
      self->z_obs_params_len = nc_cluster_redshift_class_obs_params_len (clusterz_class);

      g_assert (self->z_obs_len > 0);

      g_type_class_unref (clusterz_class);
      break;
    }
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
      self->area_survey = g_value_get_double (value);
      break;
    case PROP_USE_TRUE:
      self->use_true_data = g_value_get_boolean (value);
      break;
    case PROP_BINNED:
      nc_data_cluster_ncount_set_binned (ncount, g_value_get_boolean (value));
      break;
    case PROP_Z_OBS_BINS:
      nc_data_cluster_ncount_set_z_obs_bins (ncount, g_value_get_boxed (value));
      break;
    case PROP_LNM_OBS_BINS:
      nc_data_cluster_ncount_set_lnM_obs_bins (ncount, g_value_get_boxed (value));
      break;
    case PROP_BIN_COUNT:
      nc_data_cluster_ncount_set_bin_count (ncount, g_value_get_object (value));
      break;
    case PROP_FIDUCIAL:
      self->fiducial = g_value_get_boolean (value);
      break;
    case PROP_SEED:
      self->seed = g_value_get_uint64 (value);
      break;
    case PROP_RNG_NAME:
      g_clear_pointer (&self->rnd_name, g_free);
      self->rnd_name = g_value_dup_string (value);
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
  NcDataClusterNCountPrivate * const self = ncount->priv;

  g_return_if_fail (NC_IS_DATA_CLUSTER_NCOUNT (object));

  switch (prop_id)
  {
    case PROP_CAD:
      g_value_set_object (value, self->cad);
      break;
    case PROP_REDSHIFT_TYPE:
      g_value_set_string (value, g_type_name (self->redshift_type));
      break;
    case PROP_MASS_TYPE:
      g_value_set_string (value, g_type_name (self->mass_type));
      break;
    case PROP_LNM_TRUE:
      g_value_set_object (value, self->lnM_true);
      break;
    case PROP_Z_TRUE:
      g_value_set_object (value, self->z_true);
      break;
    case PROP_Z_OBS:
      g_value_set_object (value, self->z_obs);
      break;
    case PROP_Z_OBS_PARAMS:
      g_value_set_object (value, self->z_obs_params);
      break;
    case PROP_LNM_OBS:
      g_value_set_object (value, self->lnM_obs);
      break;
    case PROP_LNM_OBS_PARAMS:
      g_value_set_object (value, self->lnM_obs_params);
      break;
    case PROP_SURVEY_AREA:
      g_value_set_double (value, self->area_survey);
      break;
    case PROP_USE_TRUE:
      g_value_set_boolean (value, self->use_true_data);
      break;
    case PROP_BINNED:
      g_value_set_boolean (value, self->binned);
      break;
    case PROP_Z_OBS_BINS:
      g_value_set_boxed (value, self->z_obs_bins);
      break;
    case PROP_LNM_OBS_BINS:
      g_value_set_boxed (value, self->lnM_obs_bins);
      break;
    case PROP_BIN_COUNT:
      g_value_set_object (value, self->bin_count);
      break;
    case PROP_FIDUCIAL:
      g_value_set_boolean (value, self->fiducial);
      break;
    case PROP_SEED:
      g_value_set_uint64 (value, self->seed);
      break;
    case PROP_RNG_NAME:
      g_value_set_string (value, self->rnd_name);
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
  NcDataClusterNCountPrivate * const self = ncount->priv;

  nc_cluster_abundance_clear (&self->cad);

  ncm_vector_clear (&self->lnM_true);
  ncm_vector_clear (&self->z_true);

  ncm_matrix_clear (&self->z_obs);
  ncm_matrix_clear (&self->z_obs_params);
  ncm_matrix_clear (&self->lnM_obs);
  ncm_matrix_clear (&self->lnM_obs_params);

  ncm_obj_array_clear (&self->lnM_obs_bins);
  ncm_obj_array_clear (&self->z_obs_bins);
  ncm_vector_clear (&self->bin_count);

  g_clear_pointer (&self->m2lnL_a, g_array_unref);
  g_clear_pointer (&self->m2lnL_err_a, g_array_unref);
  g_clear_pointer (&self->z_order, g_array_unref);
  g_clear_pointer (&self->lnM_order, g_array_unref);

  g_clear_pointer (&self->p_z, g_array_unref);
  g_clear_pointer (&self->p_lnM, g_array_unref);
  g_clear_pointer (&self->d2n, g_array_unref);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_data_cluster_ncount_parent_class)->dispose (object);
}

static void
nc_data_cluster_ncount_finalize (GObject *object)
{
  NcDataClusterNCount *ncount = NC_DATA_CLUSTER_NCOUNT (object);
  NcDataClusterNCountPrivate * const self = ncount->priv;

  g_clear_pointer (&self->rnd_name, g_free);
  g_clear_pointer (&self->purity, gsl_histogram2d_free);
  g_clear_pointer (&self->sd_lnM, gsl_histogram2d_free);

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
                                   PROP_REDSHIFT_TYPE,
                                   g_param_spec_string ("redshift-type",
                                                        NULL,
                                                        "Cluster redshift proxy type",
                                                        g_type_name (NC_TYPE_CLUSTER_REDSHIFT),
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_MASS_TYPE,
                                   g_param_spec_string ("mass-type",
                                                        NULL,
                                                        "Cluster mass proxy type",
                                                        g_type_name (NC_TYPE_CLUSTER_MASS),
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
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
                                   PROP_Z_OBS_BINS,
                                   g_param_spec_boxed ("z-obs-bins",
                                                        NULL,
                                                        "Clusters redshifts bins",
                                                        NCM_TYPE_OBJ_ARRAY,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_LNM_OBS_BINS,
                                   g_param_spec_boxed ("lnM-obs-bins",
                                                        NULL,
                                                        "Clusters mass bins",
                                                        NCM_TYPE_OBJ_ARRAY,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_BIN_COUNT,
                                   g_param_spec_object ("bin-count",
                                                        NULL,
                                                        "Bin count",
                                                        NCM_TYPE_VECTOR,
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
  NcDataClusterNCount *ncount = NC_DATA_CLUSTER_NCOUNT (data);
  NcDataClusterNCountPrivate * const self = ncount->priv;

  return self->np;
}

static void
_nc_data_cluster_ncount_begin (NcmData *data)
{
  NcDataClusterNCount *ncount = NC_DATA_CLUSTER_NCOUNT (data);
  NcDataClusterNCountPrivate * const self = ncount->priv;
  gint signp = 0;

  self->log_np_fac = lgamma_r (self->np + 1, &signp);
}

static void
_nc_data_cluster_ncount_prepare (NcmData *data, NcmMSet *mset)
{
  NcDataClusterNCount *ncount = NC_DATA_CLUSTER_NCOUNT (data);
  NcDataClusterNCountPrivate * const self = ncount->priv;
  NcHICosmo *cosmo            = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));
  NcClusterRedshift *clusterz = NC_CLUSTER_REDSHIFT (ncm_mset_peek (mset, nc_cluster_redshift_id ()));
  NcClusterMass *clusterm     = NC_CLUSTER_MASS (ncm_mset_peek (mset, nc_cluster_mass_id ()));

  g_assert ((cosmo != NULL) && (clusterz != NULL) && (clusterm != NULL));

  g_assert (G_OBJECT_TYPE (clusterz) == self->redshift_type);
  g_assert (G_OBJECT_TYPE (clusterm) == self->mass_type);

  nc_cluster_abundance_prepare_if_needed (self->cad, cosmo, clusterz, clusterm);
}

static gchar *_nc_data_cluster_ncount_desc (NcDataClusterNCount *ncount, NcHICosmo *cosmo);

static void
_nc_data_cluster_ncount_resample (NcmData *data, NcmMSet *mset, NcmRNG *rng)
{
  NcDataClusterNCount *ncount = NC_DATA_CLUSTER_NCOUNT (data);
  NcDataClusterNCountPrivate * const self = ncount->priv;
  NcClusterAbundance *cad      = self->cad;
  NcHICosmo *cosmo             = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));
  NcClusterRedshift *clusterz  = NC_CLUSTER_REDSHIFT (ncm_mset_peek (mset, nc_cluster_redshift_id ()));
  NcClusterMass *clusterm      = NC_CLUSTER_MASS (ncm_mset_peek (mset, nc_cluster_mass_id ()));

  GArray *lnM_true_array       = NULL;
  GArray *z_true_array         = NULL;
  GArray *z_obs_array          = NULL;
  GArray *z_obs_params_array   = NULL;
  GArray *lnM_obs_array        = NULL;
  GArray *lnM_obs_params_array = NULL;
  guint total_np;
  guint i;

  gdouble *zi_obs = g_new (gdouble, self->z_obs_len);
  gdouble *zi_obs_params = self->z_obs_params_len > 0 ? g_new (gdouble, self->z_obs_params_len) : NULL;
  gdouble *lnMi_obs = g_new (gdouble, self->lnM_obs_len);
  gdouble *lnMi_obs_params = self->lnM_obs_params_len > 0 ? g_new (gdouble, self->lnM_obs_params_len) : NULL;

  g_assert (G_OBJECT_TYPE (clusterz) == self->redshift_type);
  g_assert (G_OBJECT_TYPE (clusterm) == self->mass_type);

  ncm_rng_lock (rng);
  total_np = gsl_ran_poisson (rng->r, cad->norma);
  ncm_rng_unlock (rng);

  if (total_np == 0)
  {
    self->np = 0;
    g_free (zi_obs);
    g_free (lnMi_obs);
    if (self->z_obs_params_len > 0)
      g_free (zi_obs_params);
    if (self->lnM_obs_params_len > 0)
      g_free (lnMi_obs_params);

    ncm_data_take_desc (data, _nc_data_cluster_ncount_desc (ncount, cosmo));
    return;
  }

  lnM_true_array = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), total_np);
  z_true_array   = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), total_np);

  z_obs_array = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), total_np * self->z_obs_len);
  if (self->z_obs_params_len > 0)
    z_obs_params_array = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), total_np * self->z_obs_params_len);

  lnM_obs_array = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), total_np * self->lnM_obs_len);
  if (self->lnM_obs_params_len > 0)
    lnM_obs_params_array = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), total_np * self->lnM_obs_params_len);

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

      if ( nc_cluster_redshift_resample (clusterz, cosmo, lnM_true, z_true, zi_obs, zi_obs_params, rng) &&
          nc_cluster_mass_resample (clusterm, cosmo, lnM_true, z_true, lnMi_obs, lnMi_obs_params, rng) )
      {
        g_array_append_val (lnM_true_array, lnM_true);
        g_array_append_val (z_true_array, z_true);

        g_array_append_vals (z_obs_array, zi_obs, self->z_obs_len);
        g_array_append_vals (lnM_obs_array, lnMi_obs, self->lnM_obs_len);

        if (self->z_obs_params_len > 0)
          g_array_append_vals (z_obs_params_array, zi_obs_params, self->z_obs_params_len);
        if (self->lnM_obs_params_len > 0)
          g_array_append_vals (lnM_obs_params_array, lnMi_obs_params, self->lnM_obs_params_len);
      }
    }
  }

  if (z_obs_array->len == 0 || lnM_obs_array->len == 0)
  {
    self->np = 0;
    g_free (zi_obs);
    g_free (lnMi_obs);
    if (self->z_obs_params_len > 0)
      g_free (zi_obs_params);
    if (self->lnM_obs_params_len > 0)
      g_free (lnMi_obs_params);

    ncm_vector_clear (&self->lnM_true);
    g_array_unref (lnM_true_array);
    ncm_vector_clear (&self->z_true);
    g_array_unref (z_true_array);
    ncm_matrix_clear (&self->z_obs);
    g_array_unref (z_obs_array);
    ncm_matrix_clear (&self->lnM_obs);
    g_array_unref (lnM_obs_array);

    ncm_data_take_desc (data, _nc_data_cluster_ncount_desc (ncount, cosmo));
    return;
  }

  ncm_vector_clear (&self->lnM_true);
  self->lnM_true = ncm_vector_new_array (lnM_true_array);
  g_array_unref (lnM_true_array);

  ncm_vector_clear (&self->z_true);
  self->z_true = ncm_vector_new_array (z_true_array);
  g_array_unref (z_true_array);

  ncm_matrix_clear (&self->z_obs);
  self->z_obs = ncm_matrix_new_array (z_obs_array, self->z_obs_len);
  g_array_unref (z_obs_array);

  ncm_matrix_clear (&self->lnM_obs);
  self->lnM_obs = ncm_matrix_new_array (lnM_obs_array, self->lnM_obs_len);
  g_array_unref (lnM_obs_array);

  if (self->z_obs_params_len > 0)
  {
    ncm_matrix_clear (&self->z_obs_params);
    self->z_obs_params = ncm_matrix_new_array (z_obs_params_array, self->z_obs_params_len);
    g_array_unref (z_obs_params_array);
  }

  if (self->lnM_obs_params_len > 0)
  {
    ncm_matrix_clear (&self->lnM_obs_params);
    self->lnM_obs_params = ncm_matrix_new_array (lnM_obs_params_array, self->lnM_obs_params_len);
    g_array_unref (lnM_obs_params_array);
  }

  self->np = ncm_matrix_nrows (self->z_obs);

  /* printf ("Generated %u, Expected %10.5g\n", self->np, nc_cluster_abundance_n (cad, cosmo)); */
  ncm_data_take_desc (data, _nc_data_cluster_ncount_desc (ncount, cosmo));

  if (self->binned)
    nc_data_cluster_ncount_bin_data (ncount);

  g_free (zi_obs);
  g_free (zi_obs_params);
  g_free (lnMi_obs);
  g_free (lnMi_obs_params);
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
nc_data_cluster_ncount_new (NcClusterAbundance *cad, const gchar *redshift_type, const gchar *mass_type)
{
  NcDataClusterNCount *ncount = g_object_new (NC_TYPE_DATA_CLUSTER_NCOUNT,
                                              "cluster-abundance", cad,
                                              "redshift-type", redshift_type,
                                              "mass-type", mass_type,
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
  NcDataClusterNCountPrivate * const self = ncount->priv;

  if (self->np > 0)
  {
    if (ncm_vector_len (v) != self->np)
      g_error ("nc_data_cluster_ncount_set_lnM_true: incompatible vector, the data has %u clusters and the vector length is %u.",
               self->np, ncm_vector_len (v));
    if (self->lnM_true != NULL)
      ncm_vector_memcpy (self->lnM_true, v);
    else
      self->lnM_true = ncm_vector_dup (v);
  }
  else
  {
    self->np = ncm_vector_len (v);
    self->lnM_true = ncm_vector_dup (v);
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
  NcDataClusterNCountPrivate * const self = ncount->priv;

  if (self->np > 0)
  {
    if (ncm_vector_len (v) != self->np)
      g_error ("nc_data_cluster_ncount_set_z_true: incompatible vector, the data has %u clusters and the vector length is %u.",
               self->np, ncm_vector_len (v));
    if (self->z_true != NULL)
      ncm_vector_memcpy (self->z_true, v);
    else
      self->z_true = ncm_vector_dup (v);
  }
  else
  {
    self->np = ncm_vector_len (v);
    self->z_true = ncm_vector_dup (v);
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
  NcDataClusterNCountPrivate * const self = ncount->priv;

  if (self->np > 0)
  {
    if (ncm_matrix_nrows (m) != self->np)
      g_error ("nc_data_cluster_ncount_set_lnM_obs: incompatible matrix, the data has %u clusters and the matrix has %u rows.",
               self->np, ncm_matrix_nrows (m));
    if (self->lnM_obs != NULL)
    {
      ncm_matrix_memcpy (self->lnM_obs, m);
      return;
    }
  }
  else
    self->np = ncm_matrix_nrows (m);

  if (self->lnM_obs_len != ncm_matrix_ncols (m))
    g_error ("nc_data_cluster_ncount_set_lnM_obs: incompatible matrix, NcmClusterMass object has %u points per observation and the matrix has %u cols.",
             self->lnM_obs_len, ncm_matrix_ncols (m));
  self->lnM_obs = ncm_matrix_dup (m);
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
  NcDataClusterNCountPrivate * const self = ncount->priv;

  if (self->np > 0)
  {
    if (ncm_matrix_nrows (m) != self->np)
      g_error ("nc_data_cluster_ncount_set_lnM_obs_params: incompatible matrix, the data has %u clusters and the matrix has %u rows.",
               self->np, ncm_matrix_nrows (m));
    if (self->lnM_obs_params != NULL)
    {
      ncm_matrix_memcpy (self->lnM_obs_params, m);
      return;
    }
  }
  else
    self->np = ncm_matrix_nrows (m);

  if (self->lnM_obs_params_len != ncm_matrix_ncols (m))
    g_error ("nc_data_cluster_ncount_set_lnM_obs_params: incompatible matrix, NcmClusterMass object has %u parameters per observation and the matrix has %u cols.",
             self->lnM_obs_params_len, ncm_matrix_ncols (m));
  self->lnM_obs_params = ncm_matrix_dup (m);
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
  NcDataClusterNCountPrivate * const self = ncount->priv;

  if (self->np > 0)
  {
    if (ncm_matrix_nrows (m) != self->np)
      g_error ("nc_data_cluster_ncount_set_z_obs: incompatible matrix, the data has %u clusters and the matrix has %u rows.",
               self->np, ncm_matrix_nrows (m));
    if (self->z_obs != NULL)
    {
      ncm_matrix_memcpy (self->z_obs, m);
      return;
    }
  }
  else
    self->np = ncm_matrix_nrows (m);

  if (self->z_obs_len != ncm_matrix_ncols (m))
    g_error ("nc_data_cluster_ncount_set_z_obs: incompatible matrix, NcmClusterRedshift object has %u points per observation and the matrix has %u cols.",
             self->z_obs_len, ncm_matrix_ncols (m));
  self->z_obs = ncm_matrix_dup (m);
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
  NcDataClusterNCountPrivate * const self = ncount->priv;

  if (self->np > 0)
  {
    if (ncm_matrix_nrows (m) != self->np)
      g_error ("nc_data_cluster_ncount_set_z_obs_params: incompatible matrix, the data has %u clusters and the matrix has %u rows.",
               self->np, ncm_matrix_nrows (m));
    if (self->z_obs_params != NULL)
    {
      ncm_matrix_memcpy (self->z_obs_params, m);
      return;
    }
  }
  else
    self->np = ncm_matrix_nrows (m);

  if (self->z_obs_params_len != ncm_matrix_ncols (m))
    g_error ("nc_data_cluster_ncount_set_lnM_obs: incompatible matrix, NcmClusterRedshift object has %u parameters per observation and the matrix has %u cols.",
             self->z_obs_params_len, ncm_matrix_ncols (m));
  self->z_obs_params = ncm_matrix_dup (m);
}

/**
 * nc_data_cluster_ncount_set_lnM_obs_bins:
 * @ncount: a #NcDataClusterNCount
 * @lnM_obs_bins: a #NcmObsArray
 *
 * Sets array of #NcmVector's representing the lower and upper bounds
 * of each bin.
 *
 */
void
nc_data_cluster_ncount_set_lnM_obs_bins (NcDataClusterNCount *ncount, NcmObjArray *lnM_obs_bins)
{
  NcDataClusterNCountPrivate * const self = ncount->priv;

  ncm_obj_array_clear (&self->lnM_obs_bins);
  self->lnM_obs_bins = ncm_obj_array_ref (lnM_obs_bins);
}

/**
 * nc_data_cluster_ncount_set_z_obs_bins:
 * @ncount: a #NcDataClusterNCount
 * @z_obs_bins: a #NcmObsArray
 *
 * Sets array of #NcmVector's representing the lower and upper bounds
 * of each bin.
 *
 */
void
nc_data_cluster_ncount_set_z_obs_bins (NcDataClusterNCount *ncount, NcmObjArray *z_obs_bins)
{
  NcDataClusterNCountPrivate * const self = ncount->priv;

  ncm_obj_array_clear (&self->z_obs_bins);
  self->z_obs_bins = ncm_obj_array_ref (z_obs_bins);
}

/**
 * nc_data_cluster_ncount_set_bin_count:
 * @ncount: a #NcDataClusterNCount
 * @bin_count: a #NcmVector
 *
 * Sets vector containing the bin counts.
 *
 */
void
nc_data_cluster_ncount_set_bin_count (NcDataClusterNCount *ncount, NcmVector *bin_count)
{
  NcDataClusterNCountPrivate * const self = ncount->priv;

  ncm_vector_clear (&self->bin_count);
  self->bin_count = ncm_vector_ref (bin_count);
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
  NcDataClusterNCountPrivate * const self = ncount->priv;

  return self->lnM_true != NULL;
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
  NcDataClusterNCountPrivate * const self = ncount->priv;

  return self->z_true != NULL;
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
  NcDataClusterNCountPrivate * const self = ncount->priv;

  return self->np;
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
  NcDataClusterNCountPrivate * const self = ncount->priv;

  return self->lnM_obs_len;
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
  NcDataClusterNCountPrivate * const self = ncount->priv;

  return self->lnM_obs_params_len;
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
  NcDataClusterNCountPrivate * const self = ncount->priv;

  return self->z_obs_len;
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
  NcDataClusterNCountPrivate * const self = ncount->priv;

  return self->z_obs_params_len;
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
  NcDataClusterNCountPrivate * const self = ncount->priv;

  if (self->lnM_true != NULL)
    return ncm_vector_ref (self->lnM_true);
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
  NcDataClusterNCountPrivate * const self = ncount->priv;

  if (self->z_true != NULL)
    return ncm_vector_ref (self->z_true);
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
  NcDataClusterNCountPrivate * const self = ncount->priv;

  if (self->lnM_obs != NULL)
    return ncm_matrix_ref (self->lnM_obs);
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
  NcDataClusterNCountPrivate * const self = ncount->priv;

  if (self->lnM_obs_params != NULL)
    return ncm_matrix_ref (self->lnM_obs_params);
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
  NcDataClusterNCountPrivate * const self = ncount->priv;

  if (self->z_obs != NULL)
    return ncm_matrix_ref (self->z_obs);
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
  NcDataClusterNCountPrivate * const self = ncount->priv;

  if (self->z_obs_params != NULL)
    return ncm_matrix_ref (self->z_obs_params);
  else
    return NULL;
}

/************************************************************************************************************
 * Unbinned number count data                                                                               *
 ************************************************************************************************************/

static void
_nc_data_cluster_ncount_model_init (NcDataClusterNCount *ncount)
{
  NcDataClusterNCountPrivate * const self = ncount->priv;
  NcClusterAbundance *cad = self->cad;

  cad->purity = self->purity;
  cad->sd_lnM = self->sd_lnM;

  nc_cluster_abundance_set_area (cad, self->area_survey);
}

static gchar *
_nc_data_cluster_ncount_desc (NcDataClusterNCount *ncount, NcHICosmo *cosmo)
{
  NcDataClusterNCountPrivate * const self = ncount->priv;
  NcClusterAbundance *cad = self->cad;

  return g_strdup_printf ("Cluster NCount resample %s. Generated %u from mean %10.5g (full). "
                          "Mass proxy type `%s', redshift proxy type `%s'. "
                          "Resampled in range [%8.4f, %8.4f] [%1.8e, %1.8e] and area %8.4f degrees square.",
                          self->binned ? "binned" : "unbinned",
                          self->np, cad->norma,
                          g_type_name (self->mass_type),
                          g_type_name (self->redshift_type),
                          cad->zi, cad->zf, 
                          exp (cad->lnMi), exp (cad->lnMf), 
                          self->area_survey / gsl_pow_2 (M_PI / 180.0));
}

typedef struct
{
  NcClusterAbundance *cad;
  NcDataClusterNCountPrivate *self;
  NcClusterRedshift *clusterz;
  NcClusterMass *clusterm;
  NcHICosmo *cosmo;
  gdouble v_pp;
} _Evald2N;

static void
_eval_z_p_lnM_p_d2n (glong i, glong f, gpointer data)
{
  _Evald2N *evald2n = (_Evald2N *) data;
  glong n;
  
  for (n = i; n < f; n++)
  {
    gdouble *lnMn_obs = ncm_matrix_ptr (evald2n->self->lnM_obs, n, 0);
    gdouble *lnMn_obs_params = evald2n->self->lnM_obs_params != NULL ? ncm_matrix_ptr (evald2n->self->lnM_obs_params, n, 0) : NULL;
    gdouble *zn_obs = ncm_matrix_ptr (evald2n->self->z_obs, n, 0);
    gdouble *zn_obs_params = evald2n->self->z_obs_params != NULL ? ncm_matrix_ptr (evald2n->self->z_obs_params, n, 0) : NULL;
    const gdouble mlnLn = -log (nc_cluster_abundance_z_p_lnM_p_d2n (evald2n->cad, evald2n->cosmo, evald2n->clusterz, evald2n->clusterm, lnMn_obs, lnMn_obs_params, zn_obs, zn_obs_params) * evald2n->v_pp);
    g_array_index (evald2n->self->m2lnL_a, gdouble, n) = mlnLn;
  }
}

static gint
func_eval_lnM_p_d2n (unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval)
{
  _Evald2N *evald2n = (_Evald2N *) fdata;
  NcClusterAbundance *cad = evald2n->cad;
  gint n;

  nc_cluster_mass_p_vec_z_lnMobs (evald2n->clusterm, evald2n->cosmo, x[0], evald2n->self->z_true, evald2n->self->lnM_obs, NULL, evald2n->self->p_z);
  ncm_spline2d_eval_vec_y (cad->mfp->d2NdzdlnM, x[0], evald2n->self->z_true, evald2n->self->z_order, evald2n->self->d2n);

  for (n = 0; n < fdim; n++)
  {
    const gdouble d2NdzdlnM = g_array_index (evald2n->self->d2n, gdouble, n);
    const gdouble p_M_Mobs  = g_array_index (evald2n->self->p_z, gdouble, n);
    fval[n] = d2NdzdlnM * p_M_Mobs;
  }

  return 0;
}

static void
_eval_z_p_d2n (glong i, glong f, gpointer data)
{
  _Evald2N *evald2n = (_Evald2N *) data;
  glong n;

  for (n = i; n < f; n++)
  {
    const gdouble lnMn = ncm_vector_get (evald2n->self->lnM_true, n);
    gdouble *zn_obs = ncm_matrix_ptr (evald2n->self->z_obs, n, 0);
    gdouble *zn_obs_params = evald2n->self->z_obs_params != NULL ? ncm_matrix_ptr (evald2n->self->z_obs_params, n, 0) : NULL;
    const gdouble mlnLn = -log (nc_cluster_abundance_z_p_d2n (evald2n->cad, evald2n->cosmo, evald2n->clusterz, evald2n->clusterm, lnMn, zn_obs, zn_obs_params) * evald2n->v_pp);
    g_array_index (evald2n->self->m2lnL_a, gdouble, n) = mlnLn;
  }
}

static void
_eval_lnM_p_d2n (glong i, glong f, gpointer data)
{
  _Evald2N *evald2n = (_Evald2N *) data;
  const gboolean vec_int = TRUE;
  glong n;

  if (!vec_int)
  {
    for (n = i; n < f; n++)
    {
      const gdouble zn = ncm_vector_get (evald2n->self->z_true, n);
      gdouble *lnMn_obs = ncm_matrix_ptr (evald2n->self->lnM_obs, n, 0);
      gdouble *lnMn_obs_params = evald2n->self->lnM_obs_params != NULL ? ncm_matrix_ptr (evald2n->self->lnM_obs_params, n, 0) : NULL;
      const gdouble mlnLn = -log (nc_cluster_abundance_lnM_p_d2n (evald2n->cad, evald2n->cosmo, evald2n->clusterz, evald2n->clusterm, lnMn_obs, lnMn_obs_params, zn) * evald2n->v_pp);
      g_array_index (evald2n->self->m2lnL_a, gdouble, n) = mlnLn;
    }
  }
  else
  {
    gulong len = f - i;
    gdouble *val = &g_array_index (evald2n->self->m2lnL_a, gdouble, 0);
    gdouble *err = &g_array_index (evald2n->self->m2lnL_err_a, gdouble, 0);
    gdouble lnMl, lnMu;
    gint ret;

    g_assert (i == 0);

    nc_cluster_mass_n_limits (evald2n->clusterm, evald2n->cosmo, &lnMl, &lnMu);

    gsl_sort_index ((size_t *)evald2n->self->z_order->data, ncm_vector_data (evald2n->self->z_true), 1, evald2n->self->np);
    ret = pcubature (len, func_eval_lnM_p_d2n, evald2n, 1, &lnMl, &lnMu, 0, 0.0, 1.0e-7, ERROR_INDIVIDUAL, val, err);
    g_assert (ret == 0);

    for (i = 0; i < len; i++)
    {
      g_array_index (evald2n->self->m2lnL_a, gdouble, i) = -log (g_array_index (evald2n->self->m2lnL_a, gdouble, i) * evald2n->v_pp);
    }
  }

}

static void
_eval_d2n (glong i, glong f, gpointer data)
{
  _Evald2N *evald2n = (_Evald2N *) data;
  glong n;

  for (n = i; n < f; n++)
  {
    const gdouble lnMn = ncm_vector_get (evald2n->self->lnM_true, n);
    const gdouble zn = ncm_vector_get (evald2n->self->z_true, n);
    const gdouble mlnLn = -log (nc_cluster_abundance_d2n (evald2n->cad, evald2n->cosmo, evald2n->clusterz, evald2n->clusterm, lnMn, zn) * evald2n->v_pp);
    g_array_index (evald2n->self->m2lnL_a, gdouble, n) = mlnLn;
  }
}

static void
_eval_intp_d2n (glong i, glong f, gpointer data)
{
  _Evald2N *evald2n = (_Evald2N *) data;
  glong n;

  for (n = i; n < f; n++)
  {
    const gdouble lnMn = ncm_vector_get (evald2n->self->lnM_true, n);
    const gdouble zn = ncm_vector_get (evald2n->self->z_true, n);
    const gdouble mlnLn = -log (nc_cluster_abundance_intp_d2n (evald2n->cad, evald2n->cosmo, evald2n->clusterz, evald2n->clusterm, lnMn, zn) * evald2n->v_pp);
    g_array_index (evald2n->self->m2lnL_a, gdouble, n) = mlnLn;
  }
}

static void
_nc_data_cluster_ncount_m2lnL_val (NcmData *data, NcmMSet *mset, gdouble *m2lnL)
{
  NcDataClusterNCount *ncount = NC_DATA_CLUSTER_NCOUNT (data);
  NcDataClusterNCountPrivate * const self = ncount->priv;
  NcClusterAbundance *cad     = self->cad;
  NcHICosmo *cosmo            = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));
  NcClusterRedshift *clusterz = NC_CLUSTER_REDSHIFT (ncm_mset_peek (mset, nc_cluster_redshift_id ()));
  NcClusterMass *clusterm     = NC_CLUSTER_MASS (ncm_mset_peek (mset, nc_cluster_mass_id ()));
  const gdouble volume        = nc_cluster_mass_volume (clusterm) * nc_cluster_redshift_volume (clusterz);
  _Evald2N evald2n = {cad, self, clusterz, clusterm, cosmo, volume / self->np};
  
  *m2lnL = 0.0;

  if (self->binned)
  {
    const guint len   = self->z_obs_bins->len;
    const guint nbins = len / 2;
    gdouble lambda = 0.0;
    gint i;

    if ((len != self->lnM_obs_bins->len) || (len == 0) || (len % 2 == 1))
      g_error ("_nc_data_cluster_ncount_m2lnL_val: cannot use bin data, inconsistent bins (%u %u).",
          len, self->lnM_obs_bins->len);

    g_assert_cmpuint (ncm_vector_len (self->bin_count), ==, nbins);

    if ((self->lnM_obs_params > 0) || (self->z_obs_params > 0))
      g_error ("_nc_data_cluster_ncount_m2lnL_val: binned distribution with parameters not supported yet.");

    for (i = 0; i < nbins; i++)
    {
      const gint j = 2 * i;

      NcmVector *lnM_obs_lb = NCM_VECTOR (ncm_obj_array_peek (self->lnM_obs_bins, j + 0));
      NcmVector *lnM_obs_ub = NCM_VECTOR (ncm_obj_array_peek (self->lnM_obs_bins, j + 1));
      NcmVector *z_obs_lb   = NCM_VECTOR (ncm_obj_array_peek (self->z_obs_bins, j + 0));
      NcmVector *z_obs_ub   = NCM_VECTOR (ncm_obj_array_peek (self->z_obs_bins, j + 1));
      const gdouble n_i     = ncm_vector_get (self->bin_count, i);

      const gdouble lambda_i = nc_cluster_abundance_intp_bin_d2n (self->cad, cosmo, clusterz, clusterm,
          ncm_vector_data (lnM_obs_lb),
          ncm_vector_data (lnM_obs_ub),
          NULL,
          ncm_vector_data (z_obs_lb),
          ncm_vector_data (z_obs_ub),
          NULL);

      if (n_i > 0.0)
        *m2lnL += n_i * log (lambda_i / n_i);

      lambda += lambda_i;
      *m2lnL += -(lambda_i - n_i);
    }

/*
    printf ("% 22.15g % 22.15g %e\n",
        nc_cluster_abundance_n (cad, cosmo, clusterz, clusterm),
        lambda,
        nc_cluster_abundance_n (cad, cosmo, clusterz, clusterm) / lambda - 1.0);
    fflush (stdout);
*/

    *m2lnL = -2.0 * (*m2lnL);

    return;
  }

  if (self->np == 0)
  {
    const gdouble lambda = nc_cluster_abundance_n (cad, cosmo, clusterz, clusterm);
    *m2lnL = -2.0 * lambda;
    return;
  }

  g_array_set_size (self->m2lnL_a, self->np);
  g_array_set_size (self->m2lnL_err_a, self->np);
  g_array_set_size (self->z_order, self->np);
  g_array_set_size (self->lnM_order, self->np);

  g_array_set_size (self->p_z, self->np);
  g_array_set_size (self->p_lnM, self->np);
  g_array_set_size (self->d2n, self->np);

  if (self->use_true_data)
  {
    g_assert (self->z_true);
    g_assert (self->lnM_true);
    /*ncm_func_eval_threaded_loop_full (&_eval_intp_d2n, 0, self->np, &evald2n);*/
    _eval_intp_d2n (0, self->np, &evald2n);
  }
  else
  {
    gboolean z_p   = ncm_model_check_impl_opt (NCM_MODEL (clusterz), NC_CLUSTER_REDSHIFT_P);
    gboolean lnM_p = ncm_model_check_impl_opt (NCM_MODEL (clusterm), NC_CLUSTER_MASS_P);

    if (z_p && lnM_p)
    {
      /*ncm_func_eval_threaded_loop_full (&_eval_z_p_lnM_p_d2n, 0, self->np, &evald2n);*/
      _eval_z_p_lnM_p_d2n (0, self->np, &evald2n);
    }
    else if (z_p && !lnM_p)
    {
      g_assert (self->lnM_true);
      /*ncm_func_eval_threaded_loop_full (&_eval_z_p_d2n, 0, self->np, &evald2n);*/
      _eval_z_p_d2n (0, self->np, &evald2n);
    }
    else if (!z_p && lnM_p)
    {
      g_assert (self->z_true);
      /*ncm_func_eval_threaded_loop_full (&_eval_lnM_p_d2n, 0, self->np, &evald2n);*/
      _eval_lnM_p_d2n (0, self->np, &evald2n);
    }
    else
    {
      g_assert (self->z_true);
      g_assert (self->lnM_true);
      /*ncm_func_eval_threaded_loop_full (&_eval_d2n, 0, self->np, &evald2n);*/
      _eval_d2n (0, self->np, &evald2n);
    }
  }

  {
    const gdouble lambda = nc_cluster_abundance_n (cad, cosmo, clusterz, clusterm);
    guint i;
    for (i = 0; i < self->np; i++)
    { 
      *m2lnL += g_array_index (self->m2lnL_a, gdouble, i);
    }
    *m2lnL += (lambda - self->np);
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
  NcDataClusterNCountPrivate * const self = ncount->priv;

  if (use_true_data)
  {
    g_assert (self->lnM_true != NULL);
    g_assert (self->z_true != NULL);
  }
  self->use_true_data = use_true_data;
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
  NcDataClusterNCountPrivate * const self = ncount->priv;

  return self->use_true_data;
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
  NcDataClusterNCountPrivate * const self = ncount->priv;

  self->area_survey = area_survey;
  _nc_data_cluster_ncount_model_init (ncount);

  ncm_data_resample (NCM_DATA (ncount), mset, rng);
  ncm_data_set_init (NCM_DATA (ncount), TRUE);
}

/**
 * nc_data_cluster_ncount_add_bin:
 * @ncount: a #NcDataClusterNCount
 * @lnM_obs_lb: a #NcmVector
 * @lnM_obs_ub: a #NcmVector
 * @z_obs_lb: a #NcmVector
 * @z_obs_ub: a #NcmVector
 *
 * Adds a new bin using the lower and upper bounds defined in
 * (@lnM_obs_lb, @lnM_obs_ub), (@z_obs_lb @z_obs_ub).
 *
 */
void
nc_data_cluster_ncount_add_bin (NcDataClusterNCount *ncount, NcmVector *lnM_obs_lb, NcmVector *lnM_obs_ub, NcmVector *z_obs_lb, NcmVector *z_obs_ub)
{
  NcDataClusterNCountPrivate * const self = ncount->priv;

  g_assert_cmpuint (ncm_vector_len (lnM_obs_lb), ==, self->lnM_obs_len);
  g_assert_cmpuint (ncm_vector_len (lnM_obs_ub), ==, self->lnM_obs_len);

  g_assert_cmpuint (ncm_vector_len (z_obs_ub), ==, self->z_obs_len);
  g_assert_cmpuint (ncm_vector_len (z_obs_ub), ==, self->z_obs_len);

  ncm_obj_array_add (self->lnM_obs_bins, G_OBJECT (ncm_vector_dup (lnM_obs_lb)));
  ncm_obj_array_add (self->lnM_obs_bins, G_OBJECT (ncm_vector_dup (lnM_obs_ub)));

  ncm_obj_array_add (self->z_obs_bins, G_OBJECT (ncm_vector_dup (z_obs_lb)));
  ncm_obj_array_add (self->z_obs_bins, G_OBJECT (ncm_vector_dup (z_obs_ub)));
}

/**
 * nc_data_cluster_ncount_del_bins:
 * @ncount: a #NcDataClusterNCount
 *
 * Remove all bins in @ncount.
 *
 */
void
nc_data_cluster_ncount_del_bins (NcDataClusterNCount *ncount)
{
  NcDataClusterNCountPrivate * const self = ncount->priv;

  g_ptr_array_set_size (self->lnM_obs_bins, 0);
  g_ptr_array_set_size (self->z_obs_bins, 0);

  ncm_vector_clear (&self->bin_count);
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
  NcDataClusterNCountPrivate * const self = ncount->priv;

  self->binned = on;
}

/**
 * nc_data_cluster_ncount_bin_data:
 * @ncount: a #NcDataClusterNCount
 *
 * Bin the current data present in @ncount. The bins must
 * be already set.
 *
 */
void
nc_data_cluster_ncount_bin_data (NcDataClusterNCount *ncount)
{
  NcDataClusterNCountPrivate * const self = ncount->priv;
  const guint len = self->z_obs_bins->len;
  guint i;

  if ((len != self->lnM_obs_bins->len) || (len == 0) || (len % 2 == 1))
    g_error ("_nc_data_cluster_ncount_bin_data: cannot bin data, inconsistent bins (%u %u).",
        len, self->lnM_obs_bins->len);

  if (self->np == 0)
    return;

  ncm_vector_clear (&self->bin_count);
  self->bin_count = ncm_vector_new (len / 2);

  ncm_vector_set_zero (self->bin_count);

  for (i = 0; i < self->np; i++)
  {
    NcmVector *lnM_obs_row = ncm_matrix_get_row (self->lnM_obs, i);
    NcmVector *z_obs_row   = ncm_matrix_get_row (self->z_obs, i);
    gint j;

    for (j = 0; j < len; j += 2)
    {
      NcmVector *lnM_obs_lb = NCM_VECTOR (ncm_obj_array_peek (self->lnM_obs_bins, j + 0));
      NcmVector *lnM_obs_ub = NCM_VECTOR (ncm_obj_array_peek (self->lnM_obs_bins, j + 1));
      NcmVector *z_obs_lb   = NCM_VECTOR (ncm_obj_array_peek (self->z_obs_bins, j + 0));
      NcmVector *z_obs_ub   = NCM_VECTOR (ncm_obj_array_peek (self->z_obs_bins, j + 1));
      const gint bin_index  = j / 2;

      if (ncm_vector_between (lnM_obs_row, lnM_obs_lb, lnM_obs_ub, 0) &&
          ncm_vector_between (z_obs_row, z_obs_lb, z_obs_ub, 0))
      {
        ncm_vector_addto (self->bin_count, bin_index, 1.0);
      }
    }

    ncm_vector_free (lnM_obs_row);
    ncm_vector_free (z_obs_row);
  }
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
  NcDataClusterNCountPrivate * const self = ncount->priv;
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

  g_assert (self->np > 0);
  
  g_ptr_array_set_free_func (tform_array, g_free);

  g_ptr_array_add (ttype_array, "Z_OBS");
  g_ptr_array_add (tform_array, g_strdup_printf ("%dD", self->z_obs_len));
  g_ptr_array_add (tunit_array, "REDSHIFT OBS");

  g_ptr_array_add (ttype_array, "LNM_OBS");
  g_ptr_array_add (tform_array, g_strdup_printf ("%dD", self->lnM_obs_len));
  g_ptr_array_add (tunit_array, "MASS OBS");

  if (self->z_true != NULL)
  {
    g_ptr_array_add (ttype_array, "Z_TRUE");
    g_ptr_array_add (tform_array, g_strdup ("1D"));
    g_ptr_array_add (tunit_array, "TRUE REDSHIFT");
  }

  if (self->lnM_true != NULL)
  {
    g_ptr_array_add (ttype_array, "LNM_TRUE");
    g_ptr_array_add (tform_array, g_strdup ("1D"));
    g_ptr_array_add (tunit_array, "TRUE LNM");
  }

  if (self->z_obs_params_len > 0)
  {
    g_ptr_array_add (ttype_array, "Z_OBS_PARAMS");
    g_ptr_array_add (tform_array, g_strdup_printf ("%dD", self->z_obs_params_len));
    g_ptr_array_add (tunit_array, "REDSHIFT OBS PARAMS");
  }

  if (self->lnM_obs_params_len > 0)
  {
    g_ptr_array_add (ttype_array, "LNM_OBS_PARAMS");
    g_ptr_array_add (tform_array, g_strdup_printf ("%dD", self->lnM_obs_params_len));
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
  fits_create_tbl (fptr, BINARY_TBL, self->np, tfields, (gchar **)ttype_array->pdata, (gchar **)tform_array->pdata,
                       (gchar **)tunit_array->pdata, extname, &status);
  NCM_FITS_ERROR (status);

  fits_write_key_str(fptr, "MTYPE", g_type_name (self->mass_type), "Mass proxy type", &status);
  NCM_FITS_ERROR (status);

  fits_write_key_str(fptr, "ZTYPE", g_type_name (self->redshift_type), "Redshift proxy type", &status);
  NCM_FITS_ERROR (status);

  fits_write_key(fptr, TLOGICAL, "USETDATA", &self->use_true_data, "Whether to used true data", &status);
  NCM_FITS_ERROR (status);

  {
    gdouble sarea_d = self->area_survey / gsl_pow_2 (M_PI / 180.0);
    fits_write_key(fptr, TDOUBLE, "AREA", &sarea_d, "Survey area in degree square", &status);
    NCM_FITS_ERROR (status);
  }

  {
    guint colnum = 1;
    fits_write_col (fptr, TDOUBLE, colnum, 1, 1, self->np, ncm_matrix_ptr (self->z_obs, 0, 0), &status);
    NCM_FITS_ERROR (status);

    colnum++;
    fits_write_col (fptr, TDOUBLE, colnum, 1, 1, self->np, ncm_matrix_ptr (self->lnM_obs, 0, 0), &status);
    NCM_FITS_ERROR (status);

    if (self->z_true != NULL)
    {
      colnum++;
      fits_write_col (fptr, TDOUBLE, colnum, 1, 1, self->np, ncm_vector_ptr (self->z_true, 0), &status);
      NCM_FITS_ERROR (status);
    }

    if (self->lnM_true != NULL)
    {
      colnum++;
      fits_write_col (fptr, TDOUBLE, colnum, 1, 1, self->np, ncm_vector_ptr (self->lnM_true, 0), &status);
      NCM_FITS_ERROR (status);
    }

    if (self->z_obs_params_len > 0)
    {
      colnum++;
      fits_write_col (fptr, TDOUBLE, colnum, 1, 1, self->np, ncm_matrix_ptr (self->z_obs_params, 0, 0), &status);
      NCM_FITS_ERROR (status);
    }

    if (self->lnM_obs_params_len > 0)
    {
      colnum++;
      fits_write_col (fptr, TDOUBLE, colnum, 1, 1, self->np, ncm_matrix_ptr (self->lnM_obs_params, 0, 0), &status);
      NCM_FITS_ERROR (status);
    }

  }

  {
    const guint len   = self->z_obs_bins->len;
    const guint nbins = len / 2;

    if (len > 0)
    {
      guint colnum;
      gint i;

      if ((len != self->lnM_obs_bins->len) || (len % 2 == 1))
        g_error ("nc_data_cluster_ncount_catalog_save: cannot save binning information, inconsistent bins (%u %u).",
            len, self->lnM_obs_bins->len);
      g_assert_cmpuint (ncm_vector_len (self->bin_count), ==, nbins);

      g_ptr_array_set_size (ttype_array, 0);
      g_ptr_array_set_size (tform_array, 0);
      g_ptr_array_set_size (tunit_array, 0);

      g_ptr_array_add (ttype_array, "Z_OBS_LB");
      g_ptr_array_add (tform_array, g_strdup_printf ("%dD", self->z_obs_len));
      g_ptr_array_add (tunit_array, "REDSHIFT OBS");

      g_ptr_array_add (ttype_array, "Z_OBS_UB");
      g_ptr_array_add (tform_array, g_strdup_printf ("%dD", self->z_obs_len));
      g_ptr_array_add (tunit_array, "REDSHIFT OBS");

      g_ptr_array_add (ttype_array, "LNM_OBS_LB");
      g_ptr_array_add (tform_array, g_strdup_printf ("%dD", self->lnM_obs_len));
      g_ptr_array_add (tunit_array, "MASS OBS");

      g_ptr_array_add (ttype_array, "LNM_OBS_UB");
      g_ptr_array_add (tform_array, g_strdup_printf ("%dD", self->lnM_obs_len));
      g_ptr_array_add (tunit_array, "MASS OBS");

      g_ptr_array_add (ttype_array, "BIN_COUNT");
      g_ptr_array_add (tform_array, g_strdup ("1D"));
      g_ptr_array_add (tunit_array, "COUNT");

      tfields = ttype_array->len;

      fits_create_tbl (fptr, BINARY_TBL, nbins, tfields, (gchar **)ttype_array->pdata, (gchar **)tform_array->pdata,
                           (gchar **)tunit_array->pdata, extname, &status);
      NCM_FITS_ERROR (status);

      fits_write_key(fptr, TLOGICAL, "BINNED", &self->binned, "Whether to used binned data", &status);
      NCM_FITS_ERROR (status);

      for (i = 0; i < nbins; i++)
      {
        const gint j = 2 * i;
        NcmVector *lnM_obs_lb = NCM_VECTOR (ncm_obj_array_peek (self->lnM_obs_bins, j + 0));
        NcmVector *lnM_obs_ub = NCM_VECTOR (ncm_obj_array_peek (self->lnM_obs_bins, j + 1));
        NcmVector *z_obs_lb   = NCM_VECTOR (ncm_obj_array_peek (self->z_obs_bins, j + 0));
        NcmVector *z_obs_ub   = NCM_VECTOR (ncm_obj_array_peek (self->z_obs_bins, j + 1));

        colnum = 1;

        fits_write_col (fptr, TDOUBLE, colnum, i + 1, 1, 1, ncm_vector_data (z_obs_lb), &status);
        NCM_FITS_ERROR (status);
        colnum++;

        fits_write_col (fptr, TDOUBLE, colnum, i + 1, 1, 1, ncm_vector_data (z_obs_ub), &status);
        NCM_FITS_ERROR (status);
        colnum++;

        fits_write_col (fptr, TDOUBLE, colnum, i + 1, 1, 1, ncm_vector_data (lnM_obs_lb), &status);
        NCM_FITS_ERROR (status);
        colnum++;

        fits_write_col (fptr, TDOUBLE, colnum, i + 1, 1, 1, ncm_vector_data (lnM_obs_ub), &status);
        NCM_FITS_ERROR (status);
        colnum++;
      }

      fits_write_col (fptr, TDOUBLE, colnum, 1, 1, nbins, ncm_vector_data (self->bin_count), &status);
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
  NcDataClusterNCountPrivate * const self = ncount->priv;
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
    self->np = nrows;
  }

  {
    gchar mass_type[FLEN_VALUE];
    if (fits_read_key_str (fptr, "MTYPE", mass_type, NULL, &status))
      g_error ("Fits file does not contain MTYPE in the header indicating the mass proxy object type. Use [col #ZTYPE=...] to add this information.");

    g_assert (g_type_from_name (mass_type) == self->mass_type);
  }

  {
    gchar redshift_type[FLEN_VALUE];
    if (fits_read_key_str (fptr, "ZTYPE", redshift_type, NULL, &status))
      g_error ("Fits file does not contain ZTYPE in the header indicating the redshift proxy object type. Use [col #ZTYPE=...] to add this information.");

    g_assert (g_type_from_name (redshift_type) == self->redshift_type);
  }

  if (fits_read_key (fptr, TLOGICAL, "USETDATA", &self->use_true_data, comment, &status))
    g_error ("Fits file does not contain USETDATA in the header indicating whether to used true data. Use [col #USETDATA=...] to add this information.");

  if (fits_read_key_dbl (fptr, "AREA", &self->area_survey, comment, &status))
    g_error ("Fits file does not contain AREA in the header indicating the survey area (degree square). Use [col #AREA=...] to add this information.");
  self->area_survey *= gsl_pow_2 (M_PI / 180.0);

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

    g_assert (self->z_obs_len == z_obs_rp);
    g_assert (self->lnM_obs_len == lnM_obs_rp);

    ncm_matrix_clear (&self->z_obs);
    self->z_obs = ncm_matrix_new (self->np, z_obs_rp);

    ncm_matrix_clear (&self->lnM_obs);
    self->lnM_obs = ncm_matrix_new (self->np, lnM_obs_rp);

    fits_read_col (fptr, TDOUBLE, z_obs_i, 1, 1, self->np * z_obs_rp, NULL, ncm_matrix_ptr (self->z_obs, 0, 0), NULL, &status);
    NCM_FITS_ERROR (status);

    fits_read_col (fptr, TDOUBLE, lnM_obs_i, 1, 1, self->np * lnM_obs_rp, NULL, ncm_matrix_ptr (self->lnM_obs, 0, 0), NULL, &status);
    NCM_FITS_ERROR (status);
  }

  {
    gint z_obs_params_i, lnM_obs_params_i;
    gint z_obs_params_tc, lnM_obs_params_tc;
    glong z_obs_params_rp, lnM_obs_params_rp;
    glong z_obs_params_w, lnM_obs_params_w;

    if (fits_get_colnum (fptr, CASESEN, "Z_OBS_PARAMS", &z_obs_params_i, &status))
    {
      g_assert (self->z_obs_params_len == 0);
      status = 0;
    }
    else
    {
      if (fits_get_coltype (fptr, z_obs_params_i, &z_obs_params_tc, &z_obs_params_rp, &z_obs_params_w, &status))
        g_error ("Column Z_OBS_PARAMS info not found, invalid fits file.");

      g_assert (self->z_obs_params_len == z_obs_params_rp);

      ncm_matrix_clear (&self->z_obs_params);
      self->z_obs_params = ncm_matrix_new (self->np, z_obs_params_rp);

      fits_read_col (fptr, TDOUBLE, z_obs_params_i, 1, 1, self->np * z_obs_params_rp, NULL, ncm_matrix_ptr (self->z_obs_params, 0, 0), NULL, &status);
      NCM_FITS_ERROR (status);
    }

    if (fits_get_colnum (fptr, CASESEN, "LNM_OBS_PARAMS", &lnM_obs_params_i, &status))
    {
      g_assert (self->lnM_obs_params_len == 0);
      status = 0;
    }
    else
    {
      if (fits_get_coltype (fptr, lnM_obs_params_i, &lnM_obs_params_tc, &lnM_obs_params_rp, &lnM_obs_params_w, &status))
        g_error ("Column LNM_OBS_PARAMS info not found, invalid fits file.");

      g_assert (self->lnM_obs_params_len == lnM_obs_params_rp);

      ncm_matrix_clear (&self->lnM_obs_params);
      self->lnM_obs_params = ncm_matrix_new (self->np, lnM_obs_params_rp);

      fits_read_col (fptr, TDOUBLE, lnM_obs_params_i, 1, 1, self->np * lnM_obs_params_rp, NULL, ncm_matrix_ptr (self->lnM_obs_params, 0, 0), NULL, &status);
      NCM_FITS_ERROR (status);
    }
  }

  {
    gint z_true_i, lnM_true_i;
    
    ncm_vector_clear (&self->z_true);
    if (!fits_get_colnum (fptr, CASESEN, "Z_TRUE*", &z_true_i, &status))
    {
      self->z_true = ncm_vector_new (self->np);
      fits_read_col (fptr, TDOUBLE, z_true_i, 1, 1, self->np, NULL, ncm_vector_ptr (self->z_true, 0), NULL, &status);
      NCM_FITS_ERROR (status);
    }
    else
      status = 0;

    ncm_vector_clear (&self->lnM_true);
    if (!fits_get_colnum (fptr, CASESEN, "LNM_TRUE", &lnM_true_i, &status))
    {
      self->lnM_true = ncm_vector_new (self->np);
      fits_read_col (fptr, TDOUBLE, lnM_true_i, 1, 1, self->np, NULL, ncm_vector_ptr (self->lnM_true, 0), NULL, &status);
      NCM_FITS_ERROR (status);
    }
    else
      status = 0;
  }

  {
    gint hdunum;

    fits_get_num_hdus (fptr, &hdunum, &status);
    NCM_FITS_ERROR (status);

    if (hdunum >= 3)
    {
      gint z_obs_lb_i, z_obs_ub_i, lnM_obs_lb_i, lnM_obs_ub_i, bin_count_i;
      gint z_obs_lb_tc, z_obs_ub_tc, lnM_obs_lb_tc, lnM_obs_ub_tc, bin_count_tc;
      glong z_obs_lb_rp, z_obs_ub_rp, lnM_obs_lb_rp, lnM_obs_ub_rp, bin_count_rp;
      glong z_obs_lb_w, z_obs_ub_w, lnM_obs_lb_w, lnM_obs_ub_w, bin_count_w;
      glong nbins;
      gint i;

      fits_movabs_hdu (fptr, 3, &hdutype, &status);
      NCM_FITS_ERROR (status);

      if (hdutype != BINARY_TBL)
        g_error ("%s (%d): NcDataClusterNCount catalog binning information is not binary!", __FILE__, __LINE__);

      fits_get_num_rows (fptr, &nbins, &status);
      NCM_FITS_ERROR (status);

      if (fits_read_key (fptr, TLOGICAL, "BINNED", &self->binned, comment, &status))
        g_error ("Fits file does not contain BINNED in the header indicating whether to used binned data. Use [col #BINNED=...] to add this information.");

      ncm_vector_clear (&self->bin_count);
      self->bin_count = ncm_vector_new (nbins);

      g_ptr_array_set_size (self->lnM_obs_bins, 0);
      g_ptr_array_set_size (self->z_obs_bins, 0);

      if (fits_get_colnum (fptr, CASESEN, "Z_OBS_LB", &z_obs_lb_i, &status))
        g_error ("Column Z_OBS_LB not found, invalid fits file.");

      if (fits_get_colnum (fptr, CASESEN, "Z_OBS_UB", &z_obs_ub_i, &status))
        g_error ("Column Z_OBS_UB not found, invalid fits file.");

      if (fits_get_colnum (fptr, CASESEN, "LNM_OBS_LB", &lnM_obs_lb_i, &status))
        g_error ("Column LNM_OBS_LB not found, invalid fits file.");

      if (fits_get_colnum (fptr, CASESEN, "LNM_OBS_UB", &lnM_obs_ub_i, &status))
        g_error ("Column LNM_OBS_LB not found, invalid fits file.");

      if (fits_get_colnum (fptr, CASESEN, "BIN_COUNT", &bin_count_i, &status))
        g_error ("Column LNM_OBS_LB not found, invalid fits file.");

      if (fits_get_coltype (fptr, z_obs_lb_i, &z_obs_lb_tc, &z_obs_lb_rp, &z_obs_lb_w, &status))
        g_error ("Column Z_OBS_LB info not found, invalid fits file.");

      if (fits_get_coltype (fptr, z_obs_ub_i, &z_obs_ub_tc, &z_obs_ub_rp, &z_obs_ub_w, &status))
        g_error ("Column Z_OBS_UB info not found, invalid fits file.");

      if (fits_get_coltype (fptr, lnM_obs_lb_i, &lnM_obs_lb_tc, &lnM_obs_lb_rp, &lnM_obs_lb_w, &status))
        g_error ("Column LNM_OBS_LB info not found, invalid fits file.");

      if (fits_get_coltype (fptr, lnM_obs_ub_i, &lnM_obs_ub_tc, &lnM_obs_ub_rp, &lnM_obs_ub_w, &status))
        g_error ("Column LNM_OBS_UB info not found, invalid fits file.");

      if (fits_get_coltype (fptr, bin_count_i, &bin_count_tc, &bin_count_rp, &bin_count_w, &status))
        g_error ("Column LNM_OBS_UB info not found, invalid fits file.");

      g_assert (self->z_obs_len == z_obs_lb_rp);
      g_assert (self->z_obs_len == z_obs_ub_rp);
      g_assert (self->lnM_obs_len == lnM_obs_lb_rp);
      g_assert (self->lnM_obs_len == lnM_obs_ub_rp);

      for (i = 0; i < nbins; i++)
      {
        NcmVector *z_obs_lb   = ncm_vector_new (self->z_obs_len);
        NcmVector *z_obs_ub   = ncm_vector_new (self->z_obs_len);
        NcmVector *lnM_obs_lb = ncm_vector_new (self->lnM_obs_len);
        NcmVector *lnM_obs_ub = ncm_vector_new (self->lnM_obs_len);

        fits_read_col (fptr, TDOUBLE, z_obs_lb_i, i + 1, 1, self->z_obs_len, NULL, ncm_vector_data (z_obs_lb), NULL, &status);
        NCM_FITS_ERROR (status);

        fits_read_col (fptr, TDOUBLE, z_obs_ub_i, i + 1, 1, self->z_obs_len, NULL, ncm_vector_data (z_obs_ub), NULL, &status);
        NCM_FITS_ERROR (status);

        fits_read_col (fptr, TDOUBLE, lnM_obs_lb_i, i + 1, 1, self->lnM_obs_len, NULL, ncm_vector_data (lnM_obs_lb), NULL, &status);
        NCM_FITS_ERROR (status);

        fits_read_col (fptr, TDOUBLE, lnM_obs_ub_i, i + 1, 1, self->lnM_obs_len, NULL, ncm_vector_data (lnM_obs_ub), NULL, &status);
        NCM_FITS_ERROR (status);

        ncm_obj_array_add (self->lnM_obs_bins, G_OBJECT (lnM_obs_lb));
        ncm_obj_array_add (self->lnM_obs_bins, G_OBJECT (lnM_obs_ub));

        ncm_obj_array_add (self->z_obs_bins, G_OBJECT (z_obs_lb));
        ncm_obj_array_add (self->z_obs_bins, G_OBJECT (z_obs_ub));
      }

      fits_read_col (fptr, TDOUBLE, bin_count_i, 1, 1, nbins, NULL, ncm_vector_data (self->bin_count), NULL, &status);
      NCM_FITS_ERROR (status);
    }
  }

  fits_close_file (fptr, &status);
  NCM_FITS_ERROR (status);

  _nc_data_cluster_ncount_model_init (ncount);

  ncm_data_set_init (NCM_DATA (ncount), TRUE);

  return;
}

#endif /* NUMCOSMO_HAVE_CFITSIO */

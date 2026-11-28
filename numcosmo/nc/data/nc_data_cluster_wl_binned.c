/***************************************************************************
 *            nc_data_cluster_wl_binned.c
 *
 *  Thu November 26 2026
 *  Copyright  2026  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * nc_data_cluster_wl_binned.c
 * Copyright (C) 2026 Caio Lima de Oliveira <caiolimadeoliveira@pm.me>
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
 * NcDataClusterWLBinned:
 *
 * Binned weak-lensing likelihood of a cluster.
 *
 * The likelihood of a binned weak-lensing profile (#NcGalaxyWLObsBinned) of a cluster. The
 * #NcGalaxyWLBinnedEstimator that built the profile provides its model, the classic one by
 * default or, when #NcDataClusterWLBinned:exact is %TRUE, the expected value of the estimator
 * computed galaxy by galaxy (nc_galaxy_wl_binned_estimator_eval_model_exact()). The
 * #NcGalaxyWLBinnedScatter of the #NcmMSet gives the distribution of the profile around the
 * model, for instance #NcGalaxyWLBinnedScatterGauss.
 *
 * The profile is built once, by nc_data_cluster_wl_binned_new_from_obs() or
 * nc_galaxy_wl_binned_estimator_build(), and stays fixed while the parameters change, as in
 * the classic analyses: its radii, weights and critical densities are those of the models
 * it was built with. A resample draws new profile values around the model of the current
 * parameters and keeps everything else.
 *
 * The #NcmMSet must hold a #NcHICosmo, a #NcHaloPosition, a #NcHaloDensityProfile, a
 * #NcWLSurfaceMassDensity and a #NcGalaxyWLBinnedScatter.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc/data/nc_data_cluster_wl_binned.h"
#include "nc/lss/galaxy/nc_galaxy_wl_binned_scatter.h"
#include "ncm/core/ncm_serialize.h"

struct _NcDataClusterWLBinned
{
  /*< private >*/
  NcmData parent_instance;
  NcGalaxyWLBinnedEstimator *est;
  NcGalaxyWLObsBinned *binned;
  gboolean exact;
  NcmVector *model;
};

enum
{
  PROP_0,
  PROP_ESTIMATOR,
  PROP_BINNED,
  PROP_EXACT,
  PROP_LEN,
};

G_DEFINE_TYPE (NcDataClusterWLBinned, nc_data_cluster_wl_binned, NCM_TYPE_DATA)

static void
nc_data_cluster_wl_binned_init (NcDataClusterWLBinned *dcwlb)
{
  dcwlb->est    = NULL;
  dcwlb->binned = NULL;
  dcwlb->exact  = FALSE;
  dcwlb->model  = NULL;
}

static void
_nc_data_cluster_wl_binned_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcDataClusterWLBinned *dcwlb = NC_DATA_CLUSTER_WL_BINNED (object);

  switch (prop_id)
  {
    case PROP_ESTIMATOR:
      nc_galaxy_wl_binned_estimator_clear (&dcwlb->est);
      dcwlb->est = g_value_dup_object (value);
      break;
    case PROP_BINNED:
    {
      NcGalaxyWLObsBinned *binned = g_value_get_object (value);

      if (binned != NULL)
        nc_data_cluster_wl_binned_set_binned (dcwlb, binned);

      break;
    }
    case PROP_EXACT:
      nc_data_cluster_wl_binned_set_exact (dcwlb, g_value_get_boolean (value));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_data_cluster_wl_binned_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcDataClusterWLBinned *dcwlb = NC_DATA_CLUSTER_WL_BINNED (object);

  switch (prop_id)
  {
    case PROP_ESTIMATOR:
      g_value_set_object (value, dcwlb->est);
      break;
    case PROP_BINNED:
      g_value_set_object (value, dcwlb->binned);
      break;
    case PROP_EXACT:
      g_value_set_boolean (value, dcwlb->exact);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_data_cluster_wl_binned_dispose (GObject *object)
{
  NcDataClusterWLBinned *dcwlb = NC_DATA_CLUSTER_WL_BINNED (object);

  nc_galaxy_wl_binned_estimator_clear (&dcwlb->est);
  nc_galaxy_wl_obs_binned_clear (&dcwlb->binned);
  ncm_vector_clear (&dcwlb->model);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_data_cluster_wl_binned_parent_class)->dispose (object);
}

static guint
_nc_data_cluster_wl_binned_get_len (NcmData *data)
{
  NcDataClusterWLBinned *dcwlb = NC_DATA_CLUSTER_WL_BINNED (data);

  return (dcwlb->binned != NULL) ? nc_galaxy_wl_obs_binned_len (dcwlb->binned) : 0;
}

static NcGalaxyWLBinnedScatter *
_nc_data_cluster_wl_binned_peek_scatter (NcmMSet *mset)
{
  NcmModel *scatter = ncm_mset_peek (mset, nc_galaxy_wl_binned_scatter_id ());

  if (scatter == NULL)
    g_error ("nc_data_cluster_wl_binned: the mset must hold a NcGalaxyWLBinnedScatter.");

  return NC_GALAXY_WL_BINNED_SCATTER (scatter);
}

static void
_nc_data_cluster_wl_binned_m2lnL_val (NcmData *data, NcmMSet *mset, gdouble *m2lnL)
{
  NcDataClusterWLBinned *dcwlb     = NC_DATA_CLUSTER_WL_BINNED (data);
  NcGalaxyWLBinnedScatter *scatter = _nc_data_cluster_wl_binned_peek_scatter (mset);

  nc_data_cluster_wl_binned_eval_model (dcwlb, mset, dcwlb->model);
  *m2lnL = nc_galaxy_wl_binned_scatter_m2lnL (scatter, dcwlb->binned, dcwlb->model);
}

static void
_nc_data_cluster_wl_binned_resample (NcmData *data, NcmMSet *mset, NcmRNG *rng)
{
  NcDataClusterWLBinned *dcwlb     = NC_DATA_CLUSTER_WL_BINNED (data);
  NcGalaxyWLBinnedScatter *scatter = _nc_data_cluster_wl_binned_peek_scatter (mset);
  NcmVector *value                 = ncm_vector_new (nc_galaxy_wl_obs_binned_len (dcwlb->binned));

  nc_data_cluster_wl_binned_eval_model (dcwlb, mset, dcwlb->model);
  nc_galaxy_wl_binned_scatter_gen (scatter, dcwlb->binned, dcwlb->model, rng, value);
  {
    guint b;

    for (b = 0; b < ncm_vector_len (value); b++)
      nc_galaxy_wl_obs_binned_set (dcwlb->binned, NC_GALAXY_WL_OBS_BINNED_COL_VALUE, b, ncm_vector_get (value, b));
  }

  ncm_vector_free (value);
}

static void
_nc_data_cluster_wl_binned_mean_vector (NcmData *data, NcmMSet *mset, NcmVector *mu)
{
  nc_data_cluster_wl_binned_eval_model (NC_DATA_CLUSTER_WL_BINNED (data), mset, mu);
}

/*
 * Registers the source catalog of the profile as shared, so that the copies of the dataset
 * made through @ser (for instance one per thread of a sampler) reference it instead of
 * copying it. The source catalog is only read, by the exact model; a resample writes the
 * profile, which every copy has on its own.
 */
static void
_nc_data_cluster_wl_binned_register_shared (NcmData *data, NcmSerialize *ser)
{
  NcDataClusterWLBinned *dcwlb = NC_DATA_CLUSTER_WL_BINNED (data);

  if ((dcwlb->binned != NULL) && nc_galaxy_wl_obs_binned_has_members (dcwlb->binned))
  {
    gchar *name = g_strdup_printf ("dcwlbSource%p", (void *) dcwlb);

    ncm_serialize_set (ser, nc_galaxy_wl_obs_binned_peek_source (dcwlb->binned), name, FALSE);
    g_free (name);
  }
}

static void
nc_data_cluster_wl_binned_class_init (NcDataClusterWLBinnedClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);
  NcmDataClass *data_class   = NCM_DATA_CLASS (klass);

  object_class->set_property = &_nc_data_cluster_wl_binned_set_property;
  object_class->get_property = &_nc_data_cluster_wl_binned_get_property;
  object_class->dispose      = &_nc_data_cluster_wl_binned_dispose;

  /**
   * NcDataClusterWLBinned:estimator:
   *
   * The estimator of the profile.
   */
  g_object_class_install_property (object_class,
                                   PROP_ESTIMATOR,
                                   g_param_spec_object ("estimator",
                                                        NULL,
                                                        "Estimator of the binned profile",
                                                        NC_TYPE_GALAXY_WL_BINNED_ESTIMATOR,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcDataClusterWLBinned:binned:
   *
   * The binned profile.
   */
  g_object_class_install_property (object_class,
                                   PROP_BINNED,
                                   g_param_spec_object ("binned",
                                                        NULL,
                                                        "Binned profile",
                                                        NC_TYPE_GALAXY_WL_OBS_BINNED,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcDataClusterWLBinned:exact:
   *
   * Whether the model is the galaxy-by-galaxy expected value of the estimator
   * (nc_galaxy_wl_binned_estimator_eval_model_exact()) instead of its classic model. It needs
   * a profile built from a per-galaxy catalog.
   */
  g_object_class_install_property (object_class,
                                   PROP_EXACT,
                                   g_param_spec_boolean ("exact",
                                                         NULL,
                                                         "Whether to use the exact expected value",
                                                         FALSE,
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  data_class->name            = "Cluster binned weak lensing likelihood";
  data_class->bootstrap       = FALSE;
  data_class->get_length      = &_nc_data_cluster_wl_binned_get_len;
  data_class->m2lnL_val       = &_nc_data_cluster_wl_binned_m2lnL_val;
  data_class->resample        = &_nc_data_cluster_wl_binned_resample;
  data_class->mean_vector     = &_nc_data_cluster_wl_binned_mean_vector;
  data_class->register_shared = &_nc_data_cluster_wl_binned_register_shared;
}

/**
 * nc_data_cluster_wl_binned_new:
 * @est: the estimator of the profile
 * @binned: the binned profile
 *
 * Creates the likelihood of @binned with the model of @est.
 *
 * Returns: (transfer full): a new #NcDataClusterWLBinned.
 */
NcDataClusterWLBinned *
nc_data_cluster_wl_binned_new (NcGalaxyWLBinnedEstimator *est, NcGalaxyWLObsBinned *binned)
{
  return g_object_new (NC_TYPE_DATA_CLUSTER_WL_BINNED,
                       "estimator", est,
                       "binned", binned,
                       NULL);
}

/**
 * nc_data_cluster_wl_binned_new_from_obs:
 * @est: the estimator of the profile
 * @mset: the #NcmMSet of the geometry
 * @obs: the per-galaxy catalog
 * @edges: (element-type gdouble): the bin edges (proper Mpc)
 *
 * Builds the profile of @obs with @est (nc_galaxy_wl_binned_estimator_build()) and creates
 * its likelihood.
 *
 * Returns: (transfer full): a new #NcDataClusterWLBinned.
 */
NcDataClusterWLBinned *
nc_data_cluster_wl_binned_new_from_obs (NcGalaxyWLBinnedEstimator *est, NcmMSet *mset, NcGalaxyWLObs *obs, GArray *edges)
{
  NcGalaxyWLObsBinned *binned  = nc_galaxy_wl_binned_estimator_build (est, mset, obs, edges);
  NcDataClusterWLBinned *dcwlb = nc_data_cluster_wl_binned_new (est, binned);

  nc_galaxy_wl_obs_binned_free (binned);

  return dcwlb;
}

/**
 * nc_data_cluster_wl_binned_ref:
 * @dcwlb: a #NcDataClusterWLBinned
 *
 * Increases the reference count of @dcwlb by one.
 *
 * Returns: (transfer full): @dcwlb.
 */
NcDataClusterWLBinned *
nc_data_cluster_wl_binned_ref (NcDataClusterWLBinned *dcwlb)
{
  return g_object_ref (dcwlb);
}

/**
 * nc_data_cluster_wl_binned_free:
 * @dcwlb: a #NcDataClusterWLBinned
 *
 * Decreases the reference count of @dcwlb by one.
 */
void
nc_data_cluster_wl_binned_free (NcDataClusterWLBinned *dcwlb)
{
  g_object_unref (dcwlb);
}

/**
 * nc_data_cluster_wl_binned_clear:
 * @dcwlb: a #NcDataClusterWLBinned
 *
 * Decreases the reference count of *@dcwlb by one and sets *@dcwlb to %NULL.
 */
void
nc_data_cluster_wl_binned_clear (NcDataClusterWLBinned **dcwlb)
{
  g_clear_object (dcwlb);
}

/**
 * nc_data_cluster_wl_binned_set_binned:
 * @dcwlb: a #NcDataClusterWLBinned
 * @binned: the binned profile
 *
 * Sets the binned profile, referencing @binned.
 */
void
nc_data_cluster_wl_binned_set_binned (NcDataClusterWLBinned *dcwlb, NcGalaxyWLObsBinned *binned)
{
  if (dcwlb->exact && !nc_galaxy_wl_obs_binned_has_members (binned))
    g_error ("nc_data_cluster_wl_binned_set_binned: the exact model needs a profile built from a per-galaxy catalog.");

  nc_galaxy_wl_obs_binned_clear (&dcwlb->binned);
  ncm_vector_clear (&dcwlb->model);

  dcwlb->binned = nc_galaxy_wl_obs_binned_ref (binned);
  dcwlb->model  = ncm_vector_new (nc_galaxy_wl_obs_binned_len (binned));

  ncm_data_set_init (NCM_DATA (dcwlb), TRUE);
}

/**
 * nc_data_cluster_wl_binned_peek_binned:
 * @dcwlb: a #NcDataClusterWLBinned
 *
 * Returns: (transfer none): the binned profile.
 */
NcGalaxyWLObsBinned *
nc_data_cluster_wl_binned_peek_binned (NcDataClusterWLBinned *dcwlb)
{
  return dcwlb->binned;
}

/**
 * nc_data_cluster_wl_binned_peek_estimator:
 * @dcwlb: a #NcDataClusterWLBinned
 *
 * Returns: (transfer none): the estimator of the profile.
 */
NcGalaxyWLBinnedEstimator *
nc_data_cluster_wl_binned_peek_estimator (NcDataClusterWLBinned *dcwlb)
{
  return dcwlb->est;
}

/**
 * nc_data_cluster_wl_binned_set_exact:
 * @dcwlb: a #NcDataClusterWLBinned
 * @exact: whether to use the exact expected value
 *
 * Sets #NcDataClusterWLBinned:exact.
 */
void
nc_data_cluster_wl_binned_set_exact (NcDataClusterWLBinned *dcwlb, gboolean exact)
{
  if (exact && (dcwlb->binned != NULL) && !nc_galaxy_wl_obs_binned_has_members (dcwlb->binned))
    g_error ("nc_data_cluster_wl_binned_set_exact: the exact model needs a profile built from a per-galaxy catalog.");

  dcwlb->exact = exact;
}

/**
 * nc_data_cluster_wl_binned_get_exact:
 * @dcwlb: a #NcDataClusterWLBinned
 *
 * Returns: #NcDataClusterWLBinned:exact.
 */
gboolean
nc_data_cluster_wl_binned_get_exact (NcDataClusterWLBinned *dcwlb)
{
  return dcwlb->exact;
}

/**
 * nc_data_cluster_wl_binned_eval_model:
 * @dcwlb: a #NcDataClusterWLBinned
 * @mset: a #NcmMSet
 * @model: the model, one value per bin, set on return
 *
 * Evaluates the model of the profile for the models of @mset: the classic model of the
 * estimator, or its exact expected value when #NcDataClusterWLBinned:exact is %TRUE.
 */
void
nc_data_cluster_wl_binned_eval_model (NcDataClusterWLBinned *dcwlb, NcmMSet *mset, NcmVector *model)
{
  g_assert_nonnull (dcwlb->est);
  g_assert_nonnull (dcwlb->binned);

  if (dcwlb->exact)
    nc_galaxy_wl_binned_estimator_eval_model_exact (dcwlb->est, mset, dcwlb->binned, model);
  else
    nc_galaxy_wl_binned_estimator_eval_model (dcwlb->est, mset, dcwlb->binned, model);
}


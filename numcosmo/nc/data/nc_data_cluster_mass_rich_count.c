/***************************************************************************
 *            nc_data_cluster_mass_rich_count.c
 *
 *  Fri Jul 10 00:00:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_data_cluster_mass_rich_count.c
 * Copyright (C) 2026 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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
 * NcDataClusterMassRichCount:
 *
 * Cluster mass-richness data object using true (integer) richness counts.
 *
 * #NcDataClusterMassRichCount is a data object analogous to
 * #NcDataClusterMassRich, but for the case where the richness is the actual
 * (integer) number of galaxies detected in a cluster, rather than a
 * continuous proxy already marginalized over projection effects. The
 * observed count $N$ is modeled as a Poisson realization of a rate $\lambda$
 * that is itself log-normally distributed, i.e.,
 * $$ N \sim \text{Poisson}(\lambda), \quad \ln\lambda \sim
 * \mathcal{N}(\mu(\ln M, z), \sigma(\ln M, z)), $$
 * where $\mu$ and $\sigma$ are provided by a #NcClusterMassRichness model.
 * The resulting Poisson-Lognormal probability $P(N \mid \mu,\sigma)$ is
 * evaluated using #NcmPLN1D.
 *
 * If the mass-richness model's `cut` parameter (see
 * nc_cluster_mass_richness_get_cut()) is non-zero, each cluster's
 * probability is conditioned on $N \geq N_{\rm cut}$, with
 * $N_{\rm cut} = \mathrm{round}(e^{\rm cut})$, mirroring the truncation
 * applied by #NcDataClusterMassRich for the continuous case.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc/data/nc_data_cluster_mass_rich_count.h"
#include "nc/lss/cluster/nc_cluster_mass_richness.h"

#include "ncm/core/ncm_cfg.h"
#include "ncm/core/ncm_util.h"
#include "ncm/core/ncm_rng.h"
#include "ncm/core/ncm_memory_pool.h"
#include "ncm/core/ncm_pln1d.h"

#include <math.h>

typedef struct _NcDataClusterMassRichCountPrivate
{
  NcmVector *z_cluster;
  NcmVector *lnM_cluster;
  NcmVector *N_cluster;
  NcmVector *lnM_original;
  NcmVector *z_original;
  NcmVector *N_original;
  NcmPLN1D *pln1d;
} NcDataClusterMassRichCountPrivate;

enum
{
  PROP_0,
  PROP_Z_CLUSTER,
  PROP_LNM_CLUSTER,
  PROP_N_CLUSTER,
  PROP_LNM_ORIGINAL,
  PROP_Z_ORIGINAL,
  PROP_N_ORIGINAL,
  PROP_SIZE,
};

struct _NcDataClusterMassRichCount
{
  NcmData parent_instance;
};

G_DEFINE_TYPE_WITH_PRIVATE (NcDataClusterMassRichCount, nc_data_cluster_mass_rich_count, NCM_TYPE_DATA);

static void
nc_data_cluster_mass_rich_count_init (NcDataClusterMassRichCount *dmrc)
{
  NcDataClusterMassRichCountPrivate * const self = nc_data_cluster_mass_rich_count_get_instance_private (dmrc);

  self->z_cluster    = NULL;
  self->lnM_cluster  = NULL;
  self->N_cluster    = NULL;
  self->lnM_original = NULL;
  self->z_original   = NULL;
  self->N_original   = NULL;
  self->pln1d        = ncm_pln1d_new (120);
}

static void
_nc_data_cluster_mass_rich_count_set_vector (NcmVector **target, const GValue *value)
{
  NcmVector *new_vector = g_value_get_object (value);

  ncm_vector_clear (target);

  if (new_vector == NULL)
  {
    return;
  }
  else
  {
    NcmVector *dup_vector = ncm_vector_dup (new_vector);

    *target = dup_vector;
  }
}

static void
nc_data_cluster_mass_rich_count_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcDataClusterMassRichCount *dmrc               = NC_DATA_CLUSTER_MASS_RICH_COUNT (object);
  NcDataClusterMassRichCountPrivate * const self = nc_data_cluster_mass_rich_count_get_instance_private (dmrc);

  g_return_if_fail (NC_IS_DATA_CLUSTER_MASS_RICH_COUNT (object));

  switch (prop_id)
  {
    case PROP_Z_CLUSTER:
      _nc_data_cluster_mass_rich_count_set_vector (&self->z_cluster, value);
      break;
    case PROP_LNM_CLUSTER:
      _nc_data_cluster_mass_rich_count_set_vector (&self->lnM_cluster, value);
      break;
    case PROP_N_CLUSTER:
      _nc_data_cluster_mass_rich_count_set_vector (&self->N_cluster, value);
      break;
    case PROP_LNM_ORIGINAL:
      _nc_data_cluster_mass_rich_count_set_vector (&self->lnM_original, value);
      break;
    case PROP_Z_ORIGINAL:
      _nc_data_cluster_mass_rich_count_set_vector (&self->z_original, value);
      break;
    case PROP_N_ORIGINAL:
      _nc_data_cluster_mass_rich_count_set_vector (&self->N_original, value);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
nc_data_cluster_mass_rich_count_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcDataClusterMassRichCount *dmrc               = NC_DATA_CLUSTER_MASS_RICH_COUNT (object);
  NcDataClusterMassRichCountPrivate * const self = nc_data_cluster_mass_rich_count_get_instance_private (dmrc);

  g_return_if_fail (NC_IS_DATA_CLUSTER_MASS_RICH_COUNT (object));

  switch (prop_id)
  {
    case PROP_Z_CLUSTER:
      g_value_set_object (value, self->z_cluster);
      break;
    case PROP_LNM_CLUSTER:
      g_value_set_object (value, self->lnM_cluster);
      break;
    case PROP_N_CLUSTER:
      g_value_set_object (value, self->N_cluster);
      break;
    case PROP_LNM_ORIGINAL:
      g_value_set_object (value, self->lnM_original);
      break;
    case PROP_Z_ORIGINAL:
      g_value_set_object (value, self->z_original);
      break;
    case PROP_N_ORIGINAL:
      g_value_set_object (value, self->N_original);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
nc_data_cluster_mass_rich_count_dispose (GObject *object)
{
  NcDataClusterMassRichCount *dmrc               = NC_DATA_CLUSTER_MASS_RICH_COUNT (object);
  NcDataClusterMassRichCountPrivate * const self = nc_data_cluster_mass_rich_count_get_instance_private (dmrc);

  ncm_vector_clear (&self->z_cluster);
  ncm_vector_clear (&self->lnM_cluster);
  ncm_vector_clear (&self->N_cluster);
  ncm_vector_clear (&self->lnM_original);
  ncm_vector_clear (&self->z_original);
  ncm_vector_clear (&self->N_original);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_data_cluster_mass_rich_count_parent_class)->dispose (object);
}

static void
nc_data_cluster_mass_rich_count_finalize (GObject *object)
{
  NcDataClusterMassRichCount *dmrc               = NC_DATA_CLUSTER_MASS_RICH_COUNT (object);
  NcDataClusterMassRichCountPrivate * const self = nc_data_cluster_mass_rich_count_get_instance_private (dmrc);

  ncm_pln1d_clear (&self->pln1d);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_data_cluster_mass_rich_count_parent_class)->finalize (object);
}

static guint _nc_data_cluster_mass_rich_count_get_length (NcmData *data);
static guint _nc_data_cluster_mass_rich_count_get_dof (NcmData *data);
static void _nc_data_cluster_mass_rich_count_m2lnL_val (NcmData *data, NcmMSet *mset, gdouble *m2lnL);
static void _nc_data_cluster_mass_rich_count_prepare (NcmData *data, NcmMSet *mset);
static void _nc_data_cluster_mass_rich_count_resample (NcmData *data, NcmMSet *mset, NcmRNG *rng);

static void
nc_data_cluster_mass_rich_count_class_init (NcDataClusterMassRichCountClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);
  NcmDataClass *data_class   = NCM_DATA_CLASS (klass);

  object_class->set_property = nc_data_cluster_mass_rich_count_set_property;
  object_class->get_property = nc_data_cluster_mass_rich_count_get_property;
  object_class->dispose      = nc_data_cluster_mass_rich_count_dispose;
  object_class->finalize     = nc_data_cluster_mass_rich_count_finalize;

  g_object_class_install_property (object_class,
                                   PROP_Z_CLUSTER,
                                   g_param_spec_object ("z-cluster",
                                                        NULL,
                                                        "Clusters (halo) redshift array",
                                                        NCM_TYPE_VECTOR,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_LNM_CLUSTER,
                                   g_param_spec_object ("lnM-cluster",
                                                        NULL,
                                                        "Clusters (halo) ln-mass array",
                                                        NCM_TYPE_VECTOR,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_N_CLUSTER,
                                   g_param_spec_object ("N-cluster",
                                                        NULL,
                                                        "Clusters (halo) richness (galaxy count) array",
                                                        NCM_TYPE_VECTOR,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_LNM_ORIGINAL,
                                   g_param_spec_object ("lnM-original",
                                                        NULL,
                                                        "Clusters (halo) ln-mass array with original data",
                                                        NCM_TYPE_VECTOR,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_Z_ORIGINAL,
                                   g_param_spec_object ("z-original",
                                                        NULL,
                                                        "Clusters (halo) redshift array with original data",
                                                        NCM_TYPE_VECTOR,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_N_ORIGINAL,
                                   g_param_spec_object ("N-original",
                                                        NULL,
                                                        "Clusters (halo) richness (galaxy count) array with original data",
                                                        NCM_TYPE_VECTOR,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  data_class->bootstrap  = TRUE;
  data_class->get_length = &_nc_data_cluster_mass_rich_count_get_length;
  data_class->get_dof    = &_nc_data_cluster_mass_rich_count_get_dof;
  data_class->m2lnL_val  = &_nc_data_cluster_mass_rich_count_m2lnL_val;
  data_class->prepare    = &_nc_data_cluster_mass_rich_count_prepare;
  data_class->resample   = &_nc_data_cluster_mass_rich_count_resample;
}

static guint
_nc_data_cluster_mass_rich_count_get_length (NcmData *data)
{
  NcDataClusterMassRichCount *dmrc               = NC_DATA_CLUSTER_MASS_RICH_COUNT (data);
  NcDataClusterMassRichCountPrivate * const self = nc_data_cluster_mass_rich_count_get_instance_private (dmrc);

  if (self->z_original == NULL)
    g_error ("nc_data_cluster_mass_rich_count_get_length: sample data not set, "
             "the object must be initialized calling nc_data_cluster_mass_rich_count_set_data().");

  return ncm_vector_len (self->z_cluster);
}

static guint
_nc_data_cluster_mass_rich_count_get_dof (NcmData *data)
{
  return _nc_data_cluster_mass_rich_count_get_length (data);
}

static inline gdouble
_nc_data_cluster_mass_rich_count_log1mexp (const gdouble a)
{
  /* Computes log(1 - exp(a)) for a <= 0, e.g. Machler (2012). */
  if (a > -M_LN2)
    return log (-expm1 (a));
  else
    return log1p (-exp (a));
}

static inline guint
_nc_data_cluster_mass_rich_count_get_N_cut (NcClusterMassRichness *mr)
{
  const gdouble cut   = nc_cluster_mass_richness_get_cut (mr);
  const gdouble N_cut = round (exp (cut));

  return (N_cut > 0.0) ? (guint) N_cut : 0;
}

static inline gdouble
_nc_data_cluster_mass_rich_count_m2lnL_single (NcmPLN1D *pln1d, const gdouble N_i, const gdouble mu_i, const gdouble sigma_i, const guint N_cut)
{
  const gdouble lnp = ncm_pln1d_eval_lnp (pln1d, N_i, mu_i, sigma_i);

  if (N_cut == 0)
    return -2.0 * lnp;

  {
    const gdouble ln_left_tail = ncm_pln1d_eval_range_sum_lnp (pln1d, 0, N_cut - 1, mu_i, sigma_i);
    const gdouble ln_norm      = _nc_data_cluster_mass_rich_count_log1mexp (ln_left_tail);

    return -2.0 * (lnp - ln_norm);
  }
}

static gdouble
_nc_data_cluster_mass_rich_count_m2lnL (NcDataClusterMassRichCountPrivate *self, NcClusterMassRichness *mr)
{
  const guint ncluster = ncm_vector_len (self->z_cluster);
  const guint N_cut    = _nc_data_cluster_mass_rich_count_get_N_cut (mr);
  gdouble local_m2lnL  = 0.0;
  guint i;

  for (i = 0; i < ncluster; i++)
  {
    const gdouble z_i   = ncm_vector_get (self->z_cluster, i);
    const gdouble lnM_i = ncm_vector_get (self->lnM_cluster, i);
    const gdouble N_i   = ncm_vector_get (self->N_cluster, i);
    gdouble mu_i, sigma_i;

    nc_cluster_mass_richness_mu_sigma (mr, lnM_i, z_i, &mu_i, &sigma_i);

    local_m2lnL += _nc_data_cluster_mass_rich_count_m2lnL_single (self->pln1d, N_i, mu_i, sigma_i, N_cut);
  }

  return local_m2lnL;
}

static gdouble
_nc_data_cluster_mass_rich_count_m2lnL_bootstrap (NcDataClusterMassRichCountPrivate *self, NcClusterMassRichness *mr, NcmBootstrap *bstrap)
{
  const guint bsize   = ncm_bootstrap_get_bsize (bstrap);
  const guint N_cut   = _nc_data_cluster_mass_rich_count_get_N_cut (mr);
  gdouble local_m2lnL = 0.0;
  guint k;

  for (k = 0; k < bsize; k++)
  {
    const guint i       = ncm_bootstrap_get (bstrap, k);
    const gdouble z_i   = ncm_vector_get (self->z_cluster, i);
    const gdouble lnM_i = ncm_vector_get (self->lnM_cluster, i);
    const gdouble N_i   = ncm_vector_get (self->N_cluster, i);
    gdouble mu_i, sigma_i;

    nc_cluster_mass_richness_mu_sigma (mr, lnM_i, z_i, &mu_i, &sigma_i);

    local_m2lnL += _nc_data_cluster_mass_rich_count_m2lnL_single (self->pln1d, N_i, mu_i, sigma_i, N_cut);
  }

  return local_m2lnL;
}

static void
_nc_data_cluster_mass_rich_count_m2lnL_val (NcmData *data, NcmMSet *mset, gdouble *m2lnL)
{
  NcDataClusterMassRichCount *dmrc               = NC_DATA_CLUSTER_MASS_RICH_COUNT (data);
  NcDataClusterMassRichCountPrivate * const self = nc_data_cluster_mass_rich_count_get_instance_private (dmrc);
  NcClusterMass *cluster_mass                    = NC_CLUSTER_MASS (ncm_mset_peek (mset, nc_cluster_mass_id ()));
  NcClusterMassRichness *mr                      = NC_CLUSTER_MASS_RICHNESS (cluster_mass);
  gdouble local_m2lnL                            = 0.0;

  if (!ncm_data_bootstrap_enabled (data))
  {
    local_m2lnL = _nc_data_cluster_mass_rich_count_m2lnL (self, mr);
  }
  else
  {
    NcmBootstrap *bstrap = ncm_data_peek_bootstrap (data);
    const guint bsize    = ncm_bootstrap_get_bsize (bstrap);

    if (bsize == 0)
    {
      *m2lnL = 0.0;

      return;
    }

    g_assert (ncm_bootstrap_is_init (bstrap));

    local_m2lnL = _nc_data_cluster_mass_rich_count_m2lnL_bootstrap (self, mr, bstrap);
  }

  *m2lnL = local_m2lnL;
}

static void
_nc_data_cluster_mass_rich_count_prepare (NcmData *data, NcmMSet *mset)
{
  NcClusterMass *cmass = NC_CLUSTER_MASS (ncm_mset_peek (mset, nc_cluster_mass_id ()));

  /* Currently only compatible with #NcClusterMassRichness subclasses */
  g_assert (NC_IS_CLUSTER_MASS_RICHNESS (cmass));
}

/**
 * nc_data_cluster_mass_rich_count_new:
 *
 * Creates a new #NcDataClusterMassRichCount.
 *
 * Returns: the newly created #NcDataClusterMassRichCount.
 */
NcDataClusterMassRichCount *
nc_data_cluster_mass_rich_count_new (void)
{
  NcDataClusterMassRichCount *dmrc = g_object_new (NC_TYPE_DATA_CLUSTER_MASS_RICH_COUNT,
                                                   NULL);

  return dmrc;
}

/**
 * nc_data_cluster_mass_rich_count_ref:
 * @dmrc: a #NcDataClusterMassRichCount
 *
 * Increases the reference count of @dmrc by one.
 *
 * Returns: (transfer full): @dmrc
 */
NcDataClusterMassRichCount *
nc_data_cluster_mass_rich_count_ref (NcDataClusterMassRichCount *dmrc)
{
  return g_object_ref (dmrc);
}

/**
 * nc_data_cluster_mass_rich_count_free:
 * @dmrc: a #NcDataClusterMassRichCount
 *
 * Atomically decrements the reference count of @dmrc by one. If the reference count drops to 0,
 * all memory allocated by @dmrc is released.
 *
 */
void
nc_data_cluster_mass_rich_count_free (NcDataClusterMassRichCount *dmrc)
{
  g_object_unref (dmrc);
}

/**
 * nc_data_cluster_mass_rich_count_clear:
 * @dmrc: a #NcDataClusterMassRichCount
 *
 * If @dmrc is not %NULL, unrefs it and sets the pointer to %NULL.
 *
 */
void
nc_data_cluster_mass_rich_count_clear (NcDataClusterMassRichCount **dmrc)
{
  g_clear_object (dmrc);
}

/**
 * nc_data_cluster_mass_rich_count_set_data:
 * @dmrc: a #NcDataClusterMassRichCount
 * @lnM: a #NcmVector
 * @z: a #NcmVector
 * @N: a #NcmVector containing the observed richness (galaxy count) of each cluster
 *
 * Sets the data of @dmrc.
 *
 */
void
nc_data_cluster_mass_rich_count_set_data (NcDataClusterMassRichCount *dmrc, NcmVector *lnM, NcmVector *z, NcmVector *N)
{
  NcDataClusterMassRichCountPrivate * const self = nc_data_cluster_mass_rich_count_get_instance_private (dmrc);
  const guint len                                = ncm_vector_len (lnM);

  if ((len != ncm_vector_len (z)) || (len != ncm_vector_len (N)))
  {
    g_error ("nc_data_cluster_mass_rich_count_set_data: lnM, z and N must have the same length");

    return;
  }

  ncm_vector_clear (&self->z_original);
  ncm_vector_clear (&self->lnM_original);
  ncm_vector_clear (&self->N_original);
  ncm_vector_clear (&self->z_cluster);
  ncm_vector_clear (&self->lnM_cluster);
  ncm_vector_clear (&self->N_cluster);

  self->z_original   = ncm_vector_dup (z);
  self->lnM_original = ncm_vector_dup (lnM);
  self->N_original   = ncm_vector_dup (N);

  self->z_cluster   = ncm_vector_dup (z);
  self->lnM_cluster = ncm_vector_dup (lnM);
  self->N_cluster   = ncm_vector_dup (N);

  ncm_data_set_init (NCM_DATA (dmrc), TRUE);
}

static void
_nc_data_cluster_mass_rich_count_resample (NcmData *data, NcmMSet *mset, NcmRNG *rng)
{
  NcDataClusterMassRichCount *dmrc               = NC_DATA_CLUSTER_MASS_RICH_COUNT (data);
  NcDataClusterMassRichCountPrivate * const self = nc_data_cluster_mass_rich_count_get_instance_private (dmrc);
  NcClusterMass *clusterm                        = NC_CLUSTER_MASS (ncm_mset_peek (mset, nc_cluster_mass_id ()));
  NcClusterMassRichness *mr                      = NC_CLUSTER_MASS_RICHNESS (clusterm);
  const guint np                                 = ncm_vector_len (self->z_original);
  const guint N_cut                              = _nc_data_cluster_mass_rich_count_get_N_cut (mr);
  const gboolean sample_full_dist                = nc_cluster_mass_richness_get_sample_full_dist (mr);
  GArray *lnM_array                              = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), np);
  GArray *z_array                                = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), np);
  GArray *N_array                                = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), np);
  guint i;

  /*
   * A cluster is a Poisson realization N ~ Poisson(lambda) of a rate
   * lambda = exp(mu + sigma * xi), xi ~ N(0, 1).
   *
   * When sample_full_dist is TRUE, we draw once from the full
   * Poisson-Lognormal distribution. Samples with N < N_cut are discarded,
   * resulting in fewer clusters in the resampled catalog (selection
   * effect), analogous to #NcDataClusterMassRich.
   *
   * When sample_full_dist is FALSE, we redraw until N >= N_cut, guaranteeing
   * that all input clusters produce exactly one output cluster.
   */
  for (i = 0; i < np; i++)
  {
    const gdouble lnM = ncm_vector_get (self->lnM_original, i);
    const gdouble z   = ncm_vector_get (self->z_original, i);
    gdouble mu, sigma;
    gdouble N;
    gboolean keep;

    nc_cluster_mass_richness_mu_sigma (mr, lnM, z, &mu, &sigma);

    do {
      const gdouble lambda = exp (ncm_rng_gaussian_gen (rng, mu, sigma));

      N    = ncm_rng_poisson_gen (rng, lambda);
      keep = (N >= N_cut);
    } while (!keep && !sample_full_dist);

    if (keep)
    {
      g_array_append_val (lnM_array, lnM);
      g_array_append_val (z_array, z);
      g_array_append_val (N_array, N);
    }
  }

  ncm_vector_clear (&self->lnM_cluster);
  self->lnM_cluster = ncm_vector_new_array (lnM_array);

  ncm_vector_clear (&self->z_cluster);
  self->z_cluster = ncm_vector_new_array (z_array);

  ncm_vector_clear (&self->N_cluster);
  self->N_cluster = ncm_vector_new_array (N_array);

  g_array_unref (lnM_array);
  g_array_unref (z_array);
  g_array_unref (N_array);
  ncm_data_set_init (NCM_DATA (dmrc), TRUE);
}

/**
 * nc_data_cluster_mass_rich_count_peek_lnM:
 * @dmrc: a #NcDataClusterMassRichCount
 *
 * Gets the vector containing the true mass.
 *
 * Returns: (transfer full): Mass  #NcmVector.
 */
NcmVector *
nc_data_cluster_mass_rich_count_peek_lnM (NcDataClusterMassRichCount *dmrc)
{
  NcDataClusterMassRichCountPrivate * const self = nc_data_cluster_mass_rich_count_get_instance_private (dmrc);

  return ncm_vector_ref (self->lnM_cluster);
}

/**
 * nc_data_cluster_mass_rich_count_peek_z:
 * @dmrc: a #NcDataClusterMassRichCount
 *
 * Gets the vector containing the true z.
 *
 * Returns: (transfer full): z  #NcmVector.
 */
NcmVector *
nc_data_cluster_mass_rich_count_peek_z (NcDataClusterMassRichCount *dmrc)
{
  NcDataClusterMassRichCountPrivate * const self = nc_data_cluster_mass_rich_count_get_instance_private (dmrc);

  return ncm_vector_ref (self->z_cluster);
}

/**
 * nc_data_cluster_mass_rich_count_peek_N:
 * @dmrc: a #NcDataClusterMassRichCount
 *
 * Gets the vector containing the observed richness (galaxy count).
 *
 * Returns: (transfer full): Richness count  #NcmVector.
 */
NcmVector *
nc_data_cluster_mass_rich_count_peek_N (NcDataClusterMassRichCount *dmrc)
{
  NcDataClusterMassRichCountPrivate * const self = nc_data_cluster_mass_rich_count_get_instance_private (dmrc);

  return ncm_vector_ref (self->N_cluster);
}

/**
 * nc_data_cluster_mass_rich_count_apply_cut:
 * @dmrc: a #NcDataClusterMassRichCount
 * @N_min: minimum richness (galaxy count) to keep
 *
 * Apply a cut to the data, it will discard all clusters with richness
 * (galaxy count) lower than @N_min.
 *
 */
void
nc_data_cluster_mass_rich_count_apply_cut (NcDataClusterMassRichCount *dmrc, guint N_min)
{
  NcDataClusterMassRichCountPrivate * const self = nc_data_cluster_mass_rich_count_get_instance_private (dmrc);
  guint np                                       = ncm_vector_len (self->z_cluster);
  GArray *lnM_array                              = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), np);
  GArray *z_array                                = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), np);
  GArray *N_array                                = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), np);
  guint i;

  for (i = 0; i < np; i++)
  {
    if (ncm_vector_get (self->N_cluster, i) >= N_min)
    {
      const gdouble lnM = ncm_vector_get (self->lnM_cluster, i);
      const gdouble z   = ncm_vector_get (self->z_cluster, i);
      const gdouble N   = ncm_vector_get (self->N_cluster, i);

      g_array_append_val (lnM_array, lnM);
      g_array_append_val (z_array, z);
      g_array_append_val (N_array, N);
    }
  }

  ncm_vector_clear (&self->lnM_cluster);
  self->lnM_cluster = ncm_vector_new_array (lnM_array);

  ncm_vector_clear (&self->z_cluster);
  self->z_cluster = ncm_vector_new_array (z_array);

  ncm_vector_clear (&self->N_cluster);
  self->N_cluster = ncm_vector_new_array (N_array);

  g_array_unref (lnM_array);
  g_array_unref (z_array);
  g_array_unref (N_array);
}


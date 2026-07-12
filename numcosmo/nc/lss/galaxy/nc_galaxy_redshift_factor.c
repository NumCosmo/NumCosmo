/***************************************************************************
 *            nc_galaxy_redshift_factor.c
 *
 *  Tue Jul 1 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_galaxy_redshift_factor.c
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
 * with this program. If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * NcGalaxyRedshiftFactor:
 *
 * Abstract calculator for the galaxy joint redshift distribution.
 *
 * A calculator (a plain #GObject, NOT an #NcmModel and NOT held in an #NcmMSet)
 * that produces the per-galaxy JOINT density $p(z_\mathrm{phot}, z \mid I)$ as a
 * function of the true redshift $z$. The calculator never integrates $z$ itself
 * for the likelihood: it hands out an integrand via
 * nc_galaxy_redshift_factor_integ(), and the orchestrator (#NcDataClusterWL)
 * integrates that against every other $z$-dependent factor.
 *
 * Concrete schemes live as subclasses following the #NcPowspecML pattern
 * (abstract base owns the shared machinery once; each scheme owns its observed
 * fragment, engine and sampler). The redshift schemes are the Composed scheme
 * (a #NcGalaxyRedshiftPop slot convolved with a
 * #NcGalaxyRedshiftObs slot) and the Joint scheme (a pre-tabulated
 * $p(z)$ spline).
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc/lss/galaxy/nc_galaxy_redshift_factor.h"
#include "ncm/algebra/ncm_vector.h"
#include "ncm/core/ncm_rng.h"

typedef struct _NcGalaxyRedshiftFactorPrivate
{
  gint placeholder;
} NcGalaxyRedshiftFactorPrivate;

G_DEFINE_ABSTRACT_TYPE_WITH_PRIVATE (NcGalaxyRedshiftFactor, nc_galaxy_redshift_factor, G_TYPE_OBJECT);
G_DEFINE_BOXED_TYPE (NcGalaxyRedshiftFactorData, nc_galaxy_redshift_factor_data, nc_galaxy_redshift_factor_data_ref, nc_galaxy_redshift_factor_data_unref);
NCM_UTIL_DEFINE_CALLBACK (NcGalaxyRedshiftFactorIntegrand,
                          NC_GALAXY_REDSHIFT_FACTOR_INTEGRAND,
                          nc_galaxy_redshift_factor_integrand,
                          gdouble,
                          NCM_UTIL_CALLBACK_ARGS (const gdouble z, NcGalaxyRedshiftFactorData * data),
                          NCM_UTIL_CALLBACK_ARGS (z, data))

static void
nc_galaxy_redshift_factor_init (NcGalaxyRedshiftFactor *gsdr)
{
}

static void
_nc_galaxy_redshift_factor_finalize (GObject *object)
{
  G_OBJECT_CLASS (nc_galaxy_redshift_factor_parent_class)->finalize (object);
}

/* LCOV_EXCL_START */
static void
_nc_galaxy_redshift_factor_data_init (NcGalaxyRedshiftFactor *gsdr, NcmMSet *mset, NcGalaxyRedshiftFactorData *data)
{
  g_error ("_nc_galaxy_redshift_factor_data_init: method not implemented");
}

static void
_nc_galaxy_redshift_factor_gen (NcGalaxyRedshiftFactor *gsdr, NcmMSet *mset, NcGalaxyRedshiftFactorData *data, NcmRNG *rng)
{
  g_error ("_nc_galaxy_redshift_factor_gen: method not implemented");
}

static gboolean
_nc_galaxy_redshift_factor_gen1 (NcGalaxyRedshiftFactor *gsdr, NcmMSet *mset, NcGalaxyRedshiftFactorData *data, NcmRNG *rng)
{
  g_error ("_nc_galaxy_redshift_factor_gen1: method not implemented");

  return FALSE;
}

static NcGalaxyRedshiftFactorIntegrand *
_nc_galaxy_redshift_factor_integ (NcGalaxyRedshiftFactor *gsdr, NcmMSet *mset, gboolean use_lnp)
{
  g_error ("_nc_galaxy_redshift_factor_integ: method not implemented");

  return NULL;
}

static void
_nc_galaxy_redshift_factor_get_integ_lim (NcGalaxyRedshiftFactor *gsdr, NcmMSet *mset, NcGalaxyRedshiftFactorData *data, gdouble *z_min, gdouble *z_max)
{
  g_error ("_nc_galaxy_redshift_factor_get_integ_lim: method not implemented");
}

static gdouble
_nc_galaxy_redshift_factor_norm (NcGalaxyRedshiftFactor *gsdr, NcmMSet *mset, NcGalaxyRedshiftFactorData *data)
{
  g_error ("_nc_galaxy_redshift_factor_norm: method not implemented");

  return 0.0;
}

static NcmIntegralFixed *
_nc_galaxy_redshift_factor_make_fixed_nodes (NcGalaxyRedshiftFactor *gsdr, NcmMSet *mset,
                                             NcGalaxyRedshiftFactorData *data,
                                             gdouble z_lo, gdouble z_hi,
                                             guint n_nodes, guint rule_n)
{
  g_error ("_nc_galaxy_redshift_factor_make_fixed_nodes: method not implemented");

  return NULL;
}

/* LCOV_EXCL_STOP */

static void
_nc_galaxy_redshift_factor_prepare (NcGalaxyRedshiftFactor *gsdr, NcmMSet *mset)
{
  /* Default: no models to resolve, nothing to cache. */
}

static guint64
_nc_galaxy_redshift_factor_get_hash (NcGalaxyRedshiftFactor *gsdr)
{
  /* Default: assumed to never change. */
  return 0;
}

static void
_nc_galaxy_redshift_factor_update_data (NcGalaxyRedshiftFactor *gsdr, NcGalaxyRedshiftFactorData *data)
{
  /* Default: nothing cached per-galaxy to refresh. */
}

static void
nc_galaxy_redshift_factor_class_init (NcGalaxyRedshiftFactorClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->finalize = &_nc_galaxy_redshift_factor_finalize;

  klass->data_init        = &_nc_galaxy_redshift_factor_data_init;
  klass->gen              = &_nc_galaxy_redshift_factor_gen;
  klass->gen1             = &_nc_galaxy_redshift_factor_gen1;
  klass->prepare          = &_nc_galaxy_redshift_factor_prepare;
  klass->integ            = &_nc_galaxy_redshift_factor_integ;
  klass->get_integ_lim    = &_nc_galaxy_redshift_factor_get_integ_lim;
  klass->norm             = &_nc_galaxy_redshift_factor_norm;
  klass->make_fixed_nodes = &_nc_galaxy_redshift_factor_make_fixed_nodes;
  klass->get_hash         = &_nc_galaxy_redshift_factor_get_hash;
  klass->update_data      = &_nc_galaxy_redshift_factor_update_data;
}

/**
 * nc_galaxy_redshift_factor_data_ref:
 * @data: a #NcGalaxyRedshiftFactorData
 *
 * Increments the reference count of @data by one.
 *
 * Returns: (transfer full): @data
 */
NcGalaxyRedshiftFactorData *
nc_galaxy_redshift_factor_data_ref (NcGalaxyRedshiftFactorData *data)
{
  g_atomic_ref_count_inc (&data->ref_count);

  return data;
}

/**
 * nc_galaxy_redshift_factor_data_unref:
 * @data: a #NcGalaxyRedshiftFactorData
 *
 * Decreases the reference count of @data by one. If the reference count reaches
 * 0, the data is freed.
 */
void
nc_galaxy_redshift_factor_data_unref (NcGalaxyRedshiftFactorData *data)
{
  if (g_atomic_ref_count_dec (&data->ref_count))
  {
    g_assert_nonnull (data->ldata_destroy);
    data->ldata_destroy (data->ldata);
    g_free (data);
  }
}

/**
 * nc_galaxy_redshift_factor_data_read_row:
 * @data: a #NcGalaxyRedshiftFactorData
 * @obs: a #NcGalaxyWLObs
 * @i: the row index
 *
 * Reads the galaxy redshift data from the observation.
 *
 */
void
nc_galaxy_redshift_factor_data_read_row (NcGalaxyRedshiftFactorData *data, NcGalaxyWLObs *obs, const guint i)
{
  data->z = nc_galaxy_wl_obs_get (obs, NC_GALAXY_REDSHIFT_FACTOR_COL_Z, i, NULL);
  data->ldata_read_row (data, obs, i);
}

/**
 * nc_galaxy_redshift_factor_data_write_row:
 * @data: a #NcGalaxyRedshiftFactorData
 * @obs: a #NcGalaxyWLObs
 * @i: the row index
 *
 * Writes the galaxy redshift data to the observation.
 *
 */
void
nc_galaxy_redshift_factor_data_write_row (NcGalaxyRedshiftFactorData *data, NcGalaxyWLObs *obs, const guint i)
{
  nc_galaxy_wl_obs_set (obs, NC_GALAXY_REDSHIFT_FACTOR_COL_Z, i, data->z, NULL);
  data->ldata_write_row (data, obs, i);
}

/**
 * nc_galaxy_redshift_factor_data_required_columns:
 * @data: a #NcGalaxyRedshiftFactorData
 *
 * Returns: (element-type utf8) (transfer full): the required columns for the galaxy redshift data.
 */
GList *
nc_galaxy_redshift_factor_data_required_columns (NcGalaxyRedshiftFactorData *data)
{
  GList *columns = NULL;

  columns = g_list_append (columns, g_strdup (NC_GALAXY_REDSHIFT_FACTOR_COL_Z));
  data->ldata_required_columns (data, &columns);

  return columns;
}

/**
 * nc_galaxy_redshift_factor_integrand_new: (skip)
 * @func: (scope async) (closure callback_data): a #NcGalaxyRedshiftFactorIntegrandFunc
 * @callback_data_free: (scope async) (closure callback_data): a #NcGalaxyRedshiftFactorIntegrandFreeData
 * @callback_data_copy: (scope async) (closure callback_data): a #NcGalaxyRedshiftFactorIntegrandCopyData
 * @callback_data_prepare: (scope async) (closure callback_data): a #NcGalaxyRedshiftFactorIntegrandPrepareData
 * @callback_data: a gpointer
 *
 * Creates a new integrand for the galaxy redshift data. The integrand takes the
 * true redshift @z and the per-galaxy @data and returns the joint density (or
 * its natural logarithm) at @z.
 *
 * Returns: (transfer full): a new #NcGalaxyRedshiftFactorIntegrand object.
 */

/**
 * nc_galaxy_redshift_factor_ref:
 * @gsdr: a #NcGalaxyRedshiftFactor
 *
 * Increases the reference count of @gsdr by one.
 *
 * Returns: (transfer full): @gsdr
 */
NcGalaxyRedshiftFactor *
nc_galaxy_redshift_factor_ref (NcGalaxyRedshiftFactor *gsdr)
{
  return g_object_ref (gsdr);
}

/**
 * nc_galaxy_redshift_factor_free:
 * @gsdr: a #NcGalaxyRedshiftFactor
 *
 * Decreases the reference count of @gsdr by one.
 *
 */
void
nc_galaxy_redshift_factor_free (NcGalaxyRedshiftFactor *gsdr)
{
  g_object_unref (gsdr);
}

/**
 * nc_galaxy_redshift_factor_clear:
 * @gsdr: a #NcGalaxyRedshiftFactor
 *
 * Decreases the reference count of @gsdr by one, and sets the pointer *@gsdr to
 * NULL.
 *
 */
void
nc_galaxy_redshift_factor_clear (NcGalaxyRedshiftFactor **gsdr)
{
  g_clear_object (gsdr);
}

/**
 * nc_galaxy_redshift_factor_data_new:
 * @gsdr: a #NcGalaxyRedshiftFactor
 * @mset: a #NcmMSet supplying the scheme's models
 *
 * Creates a new per-galaxy #NcGalaxyRedshiftFactorData for the scheme @gsdr,
 * delegating fragment allocation to the scheme via its data_init vfunc (which
 * reads the required models from @mset).
 *
 * Returns: (transfer full): a new #NcGalaxyRedshiftFactorData.
 */
NcGalaxyRedshiftFactorData *
nc_galaxy_redshift_factor_data_new (NcGalaxyRedshiftFactor *gsdr, NcmMSet *mset)
{
  NcGalaxyRedshiftFactorData *data = g_new0 (NcGalaxyRedshiftFactorData, 1);

  data->z                      = 0.0;
  data->ldata                  = NULL;
  data->ldata_destroy          = NULL;
  data->ldata_read_row         = NULL;
  data->ldata_write_row        = NULL;
  data->ldata_required_columns = NULL;

  g_atomic_ref_count_init (&data->ref_count);
  NC_GALAXY_REDSHIFT_FACTOR_GET_CLASS (gsdr)->data_init (gsdr, mset, data);

  g_assert_nonnull (data->ldata_destroy);
  g_assert_nonnull (data->ldata_read_row);
  g_assert_nonnull (data->ldata_write_row);
  g_assert_nonnull (data->ldata_required_columns);

  return data;
}

/**
 * nc_galaxy_redshift_factor_gen:
 * @gsdr: a #NcGalaxyRedshiftFactor
 * @mset: a #NcmMSet supplying the scheme's models
 * @data: a #NcGalaxyRedshiftFactorData
 * @rng: a #NcmRNG
 *
 * Generates a new galaxy redshift sample into @data (draws the true redshift and
 * the redshift observation) using the models in @mset. The @data object must be
 * initialized beforehand.
 *
 */
void
nc_galaxy_redshift_factor_gen (NcGalaxyRedshiftFactor *gsdr, NcmMSet *mset, NcGalaxyRedshiftFactorData *data, NcmRNG *rng)
{
  NC_GALAXY_REDSHIFT_FACTOR_GET_CLASS (gsdr)->gen (gsdr, mset, data, rng);
}

/**
 * nc_galaxy_redshift_factor_gen1:
 * @gsdr: a #NcGalaxyRedshiftFactor
 * @mset: a #NcmMSet supplying the scheme's models
 * @data: a pre-initialized #NcGalaxyRedshiftFactorData
 * @rng: a #NcmRNG
 *
 * Attempts to generate a single redshift sample subject to the scheme's
 * selection/limits, storing the result in @data. Used for fixed-total-count
 * subsampling where proposals violating the constraints are discarded.
 *
 * Returns: %TRUE if a valid redshift was generated; %FALSE otherwise.
 */
gboolean
nc_galaxy_redshift_factor_gen1 (NcGalaxyRedshiftFactor *gsdr, NcmMSet *mset, NcGalaxyRedshiftFactorData *data, NcmRNG *rng)
{
  return NC_GALAXY_REDSHIFT_FACTOR_GET_CLASS (gsdr)->gen1 (gsdr, mset, data, rng);
}

/**
 * nc_galaxy_redshift_factor_prepare:
 * @gsdr: a #NcGalaxyRedshiftFactor
 * @mset: a #NcmMSet supplying the scheme's models
 *
 * Factory-level prepare: validates the scheme's models in @mset and
 * refreshes whatever it caches for efficient subsequent
 * nc_galaxy_redshift_factor_update_data() calls.
 *
 */
void
nc_galaxy_redshift_factor_prepare (NcGalaxyRedshiftFactor *gsdr, NcmMSet *mset)
{
  NC_GALAXY_REDSHIFT_FACTOR_GET_CLASS (gsdr)->prepare (gsdr, mset);
}

/**
 * nc_galaxy_redshift_factor_get_hash:
 * @gsdr: a #NcGalaxyRedshiftFactor
 *
 * Returns: an opaque value that changes whenever the last
 * nc_galaxy_redshift_factor_prepare() call refreshed something relevant
 * (default: a constant).
 */
guint64
nc_galaxy_redshift_factor_get_hash (NcGalaxyRedshiftFactor *gsdr)
{
  return NC_GALAXY_REDSHIFT_FACTOR_GET_CLASS (gsdr)->get_hash (gsdr);
}

/**
 * nc_galaxy_redshift_factor_update_data:
 * @gsdr: a #NcGalaxyRedshiftFactor
 * @data: a #NcGalaxyRedshiftFactorData
 *
 * Unconditionally refreshes @data's cached state from what the last
 * nc_galaxy_redshift_factor_prepare() call resolved (default: no-op).
 */
void
nc_galaxy_redshift_factor_update_data (NcGalaxyRedshiftFactor *gsdr, NcGalaxyRedshiftFactorData *data)
{
  NC_GALAXY_REDSHIFT_FACTOR_GET_CLASS (gsdr)->update_data (gsdr, data);
}

/**
 * nc_galaxy_redshift_factor_integ:
 * @gsdr: a #NcGalaxyRedshiftFactor
 * @mset: a #NcmMSet supplying the scheme's models
 * @use_lnp: if %TRUE the integrand returns the natural logarithm of the density
 *
 * Builds the per-galaxy joint integrand $p(z_\mathrm{phot}, z \mid I)$ as a
 * function of the true redshift $z$, resolving the scheme's models from @mset.
 * The calculator does not integrate $z$: the returned integrand is consumed by
 * the orchestrator's single $z$-integral.
 *
 * Returns: (transfer full): a new #NcGalaxyRedshiftFactorIntegrand.
 */
NcGalaxyRedshiftFactorIntegrand *
nc_galaxy_redshift_factor_integ (NcGalaxyRedshiftFactor *gsdr, NcmMSet *mset, gboolean use_lnp)
{
  return NC_GALAXY_REDSHIFT_FACTOR_GET_CLASS (gsdr)->integ (gsdr, mset, use_lnp);
}

/**
 * nc_galaxy_redshift_factor_get_integ_lim:
 * @gsdr: a #NcGalaxyRedshiftFactor
 * @mset: a #NcmMSet
 * @data: a #NcGalaxyRedshiftFactorData
 * @z_min: (out): the minimum redshift for integration
 * @z_max: (out): the maximum redshift for integration
 *
 * Gets the effective redshift integration support for @data: the range over
 * which the per-galaxy $p(z)$ factor is non-negligible. This is the domain used
 * by both the fixed-node and adaptive quadratures.
 *
 */
void
nc_galaxy_redshift_factor_get_integ_lim (NcGalaxyRedshiftFactor *gsdr, NcmMSet *mset, NcGalaxyRedshiftFactorData *data, gdouble *z_min, gdouble *z_max)
{
  NC_GALAXY_REDSHIFT_FACTOR_GET_CLASS (gsdr)->get_integ_lim (gsdr, mset, data, z_min, z_max);
}

/**
 * nc_galaxy_redshift_factor_norm:
 * @gsdr: a #NcGalaxyRedshiftFactor
 * @mset: a #NcmMSet
 * @data: a #NcGalaxyRedshiftFactorData
 *
 * Computes the normalization $\int p(z_\mathrm{phot}, z \mid I)\, \mathrm{d}z$
 * over the integration support (standalone/population facet, not the per-galaxy
 * likelihood path).
 *
 * Returns: the normalization integral.
 */
gdouble
nc_galaxy_redshift_factor_norm (NcGalaxyRedshiftFactor *gsdr, NcmMSet *mset, NcGalaxyRedshiftFactorData *data)
{
  return NC_GALAXY_REDSHIFT_FACTOR_GET_CLASS (gsdr)->norm (gsdr, mset, data);
}

/**
 * nc_galaxy_redshift_factor_make_fixed_nodes: (skip)
 * @gsdr: a #NcGalaxyRedshiftFactor
 * @mset: a #NcmMSet
 * @data: a #NcGalaxyRedshiftFactorData
 * @z_lo: the lower integration bound
 * @z_hi: the upper integration bound
 * @n_nodes: the number of sub-intervals
 * @rule_n: the fixed rule order per sub-interval
 *
 * Builds the fixed integration nodes over [@z_lo, @z_hi] used to marginalize the
 * true redshift.
 *
 * Returns: (transfer full): a new #NcmIntegralFixed.
 */
NcmIntegralFixed *
nc_galaxy_redshift_factor_make_fixed_nodes (NcGalaxyRedshiftFactor *gsdr, NcmMSet *mset, NcGalaxyRedshiftFactorData *data, gdouble z_lo, gdouble z_hi, guint n_nodes, guint rule_n)
{
  return NC_GALAXY_REDSHIFT_FACTOR_GET_CLASS (gsdr)->make_fixed_nodes (gsdr, mset, data, z_lo, z_hi, n_nodes, rule_n);
}


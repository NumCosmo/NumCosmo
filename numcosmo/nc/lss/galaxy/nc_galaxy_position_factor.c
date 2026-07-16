/***************************************************************************
 *            nc_galaxy_position_factor.c
 *
 *  Wed Jul 2 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 *  Copyright  2026  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * nc_galaxy_position_factor.c
 * Copyright (C) 2026 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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
 * with this program. If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * NcGalaxyPositionFactor:
 *
 * Abstract calculator for the galaxy sky-position distribution.
 *
 * A calculator (a plain #GObject, NOT an #NcmModel and NOT held in an #NcmMSet)
 * that produces the per-galaxy position density $p(\mathrm{ra}, \mathrm{dec}
 * \mid I)$. The galaxy position is observed directly — there is no
 * true-vs-observed scatter split and (for the flat scheme) no fitted
 * parameters — so, unlike #NcGalaxyRedshiftFactor, the position collapses to a
 * single distribution with no companion #NcmModel: the scheme carries its own
 * survey geometry as configuration.
 *
 * Concrete schemes live as subclasses following the #NcPowspecML pattern (the
 * abstract base owns the shared machinery once; each scheme owns its observed
 * fragment, engine and sampler). The only scheme today is the flat scheme
 * (#NcGalaxyPositionFactorFlat, uniform over a sky footprint); a future
 * clustered scheme would gain fitted parameters and read a cosmology/halo model
 * from @mset, which is why @mset is threaded through the vtable even though the
 * flat scheme ignores it.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc/lss/galaxy/nc_galaxy_position_factor.h"
#include "ncm/core/ncm_rng.h"

typedef struct _NcGalaxyPositionFactorPrivate
{
  gint placeholder;
} NcGalaxyPositionFactorPrivate;

G_DEFINE_ABSTRACT_TYPE_WITH_PRIVATE (NcGalaxyPositionFactor, nc_galaxy_position_factor, G_TYPE_OBJECT);
G_DEFINE_BOXED_TYPE (NcGalaxyPositionFactorData, nc_galaxy_position_factor_data, nc_galaxy_position_factor_data_ref, nc_galaxy_position_factor_data_unref);
NCM_UTIL_DEFINE_CALLBACK (NcGalaxyPositionFactorIntegrand,
                          NC_GALAXY_POSITION_FACTOR_INTEGRAND,
                          nc_galaxy_position_factor_integrand,
                          gdouble,
                          NCM_UTIL_CALLBACK_ARGS (NcGalaxyPositionFactorData * data),
                          NCM_UTIL_CALLBACK_ARGS (data))

static void
nc_galaxy_position_factor_init (NcGalaxyPositionFactor *gspf)
{
}

static void
_nc_galaxy_position_factor_finalize (GObject *object)
{
  G_OBJECT_CLASS (nc_galaxy_position_factor_parent_class)->finalize (object);
}

/* LCOV_EXCL_START */
static void
_nc_galaxy_position_factor_data_init (NcGalaxyPositionFactor *gspf, NcmMSet *mset, NcGalaxyPositionFactorData *data)
{
  g_error ("_nc_galaxy_position_factor_data_init: method not implemented");
}

static void
_nc_galaxy_position_factor_gen (NcGalaxyPositionFactor *gspf, NcmMSet *mset, NcGalaxyPositionFactorData *data, NcmRNG *rng)
{
  g_error ("_nc_galaxy_position_factor_gen: method not implemented");
}

static NcGalaxyPositionFactorIntegrand *
_nc_galaxy_position_factor_integ (NcGalaxyPositionFactor *gspf, NcmMSet *mset, gboolean use_lnp)
{
  g_error ("_nc_galaxy_position_factor_integ: method not implemented");

  return NULL;
}

/* LCOV_EXCL_STOP */

static void
_nc_galaxy_position_factor_prepare (NcGalaxyPositionFactor *gspf, NcmMSet *mset)
{
  /* Default: no models to resolve, nothing to cache. */
}

static guint64
_nc_galaxy_position_factor_get_hash (NcGalaxyPositionFactor *gspf)
{
  /* Default: assumed to never change. */
  return 0;
}

static void
_nc_galaxy_position_factor_update_data (NcGalaxyPositionFactor *gspf, NcGalaxyPositionFactorData *data)
{
  /* Default: nothing cached per-galaxy to refresh. */
}

static void
nc_galaxy_position_factor_class_init (NcGalaxyPositionFactorClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->finalize = &_nc_galaxy_position_factor_finalize;

  klass->data_init   = &_nc_galaxy_position_factor_data_init;
  klass->gen         = &_nc_galaxy_position_factor_gen;
  klass->prepare     = &_nc_galaxy_position_factor_prepare;
  klass->integ       = &_nc_galaxy_position_factor_integ;
  klass->get_hash    = &_nc_galaxy_position_factor_get_hash;
  klass->update_data = &_nc_galaxy_position_factor_update_data;
}

/**
 * nc_galaxy_position_factor_data_ref:
 * @data: a #NcGalaxyPositionFactorData
 *
 * Increments the reference count of @data by one.
 *
 * Returns: (transfer full): @data
 */
NcGalaxyPositionFactorData *
nc_galaxy_position_factor_data_ref (NcGalaxyPositionFactorData *data)
{
  g_atomic_ref_count_inc (&data->ref_count);

  return data;
}

/**
 * nc_galaxy_position_factor_data_unref:
 * @data: a #NcGalaxyPositionFactorData
 *
 * Decreases the reference count of @data by one. If the reference count reaches
 * 0, the data is freed.
 */
void
nc_galaxy_position_factor_data_unref (NcGalaxyPositionFactorData *data)
{
  if (g_atomic_ref_count_dec (&data->ref_count))
  {
    g_assert_nonnull (data->ldata_destroy);
    data->ldata_destroy (data->ldata);
    g_free (data);
  }
}

/**
 * nc_galaxy_position_factor_data_read_row:
 * @data: a #NcGalaxyPositionFactorData
 * @obs: a #NcGalaxyWLObs
 * @i: the row index
 *
 * Reads the galaxy position data from the observation.
 *
 */
void
nc_galaxy_position_factor_data_read_row (NcGalaxyPositionFactorData *data, NcGalaxyWLObs *obs, const guint i)
{
  data->ra  = nc_galaxy_wl_obs_get (obs, NC_GALAXY_POSITION_FACTOR_COL_RA, i, NULL);
  data->dec = nc_galaxy_wl_obs_get (obs, NC_GALAXY_POSITION_FACTOR_COL_DEC, i, NULL);
  data->ldata_read_row (data, obs, i);
}

/**
 * nc_galaxy_position_factor_data_write_row:
 * @data: a #NcGalaxyPositionFactorData
 * @obs: a #NcGalaxyWLObs
 * @i: the row index
 *
 * Writes the galaxy position data to the observation.
 *
 */
void
nc_galaxy_position_factor_data_write_row (NcGalaxyPositionFactorData *data, NcGalaxyWLObs *obs, const guint i)
{
  nc_galaxy_wl_obs_set (obs, NC_GALAXY_POSITION_FACTOR_COL_RA, i, data->ra, NULL);
  nc_galaxy_wl_obs_set (obs, NC_GALAXY_POSITION_FACTOR_COL_DEC, i, data->dec, NULL);
  data->ldata_write_row (data, obs, i);
}

/**
 * nc_galaxy_position_factor_data_required_columns:
 * @data: a #NcGalaxyPositionFactorData
 *
 * Returns: (element-type utf8) (transfer full): the required columns for the galaxy position data.
 */
GList *
nc_galaxy_position_factor_data_required_columns (NcGalaxyPositionFactorData *data)
{
  GList *columns = NULL;

  columns = g_list_append (columns, g_strdup (NC_GALAXY_POSITION_FACTOR_COL_RA));
  columns = g_list_append (columns, g_strdup (NC_GALAXY_POSITION_FACTOR_COL_DEC));
  data->ldata_required_columns (data, &columns);

  return columns;
}

/**
 * nc_galaxy_position_factor_integrand_new: (skip)
 * @func: (scope async) (closure callback_data): a #NcGalaxyPositionFactorIntegrandFunc
 * @callback_data_free: (scope async) (closure callback_data): a #NcGalaxyPositionFactorIntegrandFreeData
 * @callback_data_copy: (scope async) (closure callback_data): a #NcGalaxyPositionFactorIntegrandCopyData
 * @callback_data_prepare: (scope async) (closure callback_data): a #NcGalaxyPositionFactorIntegrandPrepareData
 * @callback_data: a gpointer
 *
 * Creates a new integrand for the galaxy position data. The integrand takes the
 * per-galaxy @data and returns the position density (or its natural logarithm)
 * at the galaxy's stored (ra, dec).
 *
 * Returns: (transfer full): a new #NcGalaxyPositionFactorIntegrand object.
 */

/**
 * nc_galaxy_position_factor_ref:
 * @gspf: a #NcGalaxyPositionFactor
 *
 * Increases the reference count of @gspf by one.
 *
 * Returns: (transfer full): @gspf
 */
NcGalaxyPositionFactor *
nc_galaxy_position_factor_ref (NcGalaxyPositionFactor *gspf)
{
  return g_object_ref (gspf);
}

/**
 * nc_galaxy_position_factor_free:
 * @gspf: a #NcGalaxyPositionFactor
 *
 * Decreases the reference count of @gspf by one.
 *
 */
void
nc_galaxy_position_factor_free (NcGalaxyPositionFactor *gspf)
{
  g_object_unref (gspf);
}

/**
 * nc_galaxy_position_factor_clear:
 * @gspf: a #NcGalaxyPositionFactor
 *
 * Decreases the reference count of @gspf by one, and sets the pointer *@gspf to
 * NULL.
 *
 */
void
nc_galaxy_position_factor_clear (NcGalaxyPositionFactor **gspf)
{
  g_clear_object (gspf);
}

/**
 * nc_galaxy_position_factor_data_new:
 * @gspf: a #NcGalaxyPositionFactor
 * @mset: a #NcmMSet supplying the scheme's models
 *
 * Creates a new per-galaxy #NcGalaxyPositionFactorData for the scheme @gspf,
 * delegating fragment allocation to the scheme via its data_init vfunc.
 *
 * Returns: (transfer full): a new #NcGalaxyPositionFactorData.
 */
NcGalaxyPositionFactorData *
nc_galaxy_position_factor_data_new (NcGalaxyPositionFactor *gspf, NcmMSet *mset)
{
  NcGalaxyPositionFactorData *data = g_new0 (NcGalaxyPositionFactorData, 1);

  data->ra                     = 0.0;
  data->dec                    = 0.0;
  data->ldata                  = NULL;
  data->ldata_destroy          = NULL;
  data->ldata_read_row         = NULL;
  data->ldata_write_row        = NULL;
  data->ldata_required_columns = NULL;

  g_atomic_ref_count_init (&data->ref_count);
  NC_GALAXY_POSITION_FACTOR_GET_CLASS (gspf)->data_init (gspf, mset, data);

  g_assert_nonnull (data->ldata_destroy);
  g_assert_nonnull (data->ldata_read_row);
  g_assert_nonnull (data->ldata_write_row);
  g_assert_nonnull (data->ldata_required_columns);

  return data;
}

/**
 * nc_galaxy_position_factor_gen:
 * @gspf: a #NcGalaxyPositionFactor
 * @mset: a #NcmMSet supplying the scheme's models
 * @data: a #NcGalaxyPositionFactorData
 * @rng: a #NcmRNG
 *
 * Generates a new galaxy position sample into @data. The @data object must be
 * initialized beforehand.
 *
 */
void
nc_galaxy_position_factor_gen (NcGalaxyPositionFactor *gspf, NcmMSet *mset, NcGalaxyPositionFactorData *data, NcmRNG *rng)
{
  NC_GALAXY_POSITION_FACTOR_GET_CLASS (gspf)->gen (gspf, mset, data, rng);
}

/**
 * nc_galaxy_position_factor_prepare:
 * @gspf: a #NcGalaxyPositionFactor
 * @mset: a #NcmMSet supplying the scheme's models
 *
 * Factory-level prepare: validates the scheme's models in @mset and
 * refreshes whatever it caches for efficient subsequent
 * nc_galaxy_position_factor_update_data() calls. Cheap/no-op for the flat
 * scheme, which has no models to resolve.
 *
 */
void
nc_galaxy_position_factor_prepare (NcGalaxyPositionFactor *gspf, NcmMSet *mset)
{
  NC_GALAXY_POSITION_FACTOR_GET_CLASS (gspf)->prepare (gspf, mset);
}

/**
 * nc_galaxy_position_factor_get_hash:
 * @gspf: a #NcGalaxyPositionFactor
 *
 * Returns: an opaque value that changes whenever the last
 * nc_galaxy_position_factor_prepare() call refreshed something relevant
 * (default: a constant).
 */
guint64
nc_galaxy_position_factor_get_hash (NcGalaxyPositionFactor *gspf)
{
  return NC_GALAXY_POSITION_FACTOR_GET_CLASS (gspf)->get_hash (gspf);
}

/**
 * nc_galaxy_position_factor_update_data:
 * @gspf: a #NcGalaxyPositionFactor
 * @data: a #NcGalaxyPositionFactorData
 *
 * Unconditionally refreshes @data's cached state from what the last
 * nc_galaxy_position_factor_prepare() call resolved (default: no-op).
 */
void
nc_galaxy_position_factor_update_data (NcGalaxyPositionFactor *gspf, NcGalaxyPositionFactorData *data)
{
  NC_GALAXY_POSITION_FACTOR_GET_CLASS (gspf)->update_data (gspf, data);
}

/**
 * nc_galaxy_position_factor_integ:
 * @gspf: a #NcGalaxyPositionFactor
 * @mset: a #NcmMSet supplying the scheme's models
 * @use_lnp: if %TRUE the integrand returns the natural logarithm of the density
 *
 * Builds the per-galaxy position density $p(\mathrm{ra}, \mathrm{dec} \mid I)$
 * as a callback evaluated at the galaxy's stored position. The position is not
 * integrated (it is observed directly); the callback is a fixed-position
 * evaluation.
 *
 * Returns: (transfer full): a new #NcGalaxyPositionFactorIntegrand.
 */
NcGalaxyPositionFactorIntegrand *
nc_galaxy_position_factor_integ (NcGalaxyPositionFactor *gspf, NcmMSet *mset, gboolean use_lnp)
{
  return NC_GALAXY_POSITION_FACTOR_GET_CLASS (gspf)->integ (gspf, mset, use_lnp);
}


/***************************************************************************
 *            nc_galaxy_sd_obs_redshift.c
 *
 *  Thu Aug 1 00:45:32 2024
 *  Copyright  2024  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 *  Copyright  2024  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_galaxy_sd_obs_redshift.c
 * Copyright (C) 2024 Caio Lima de Oliveira <caiolimadeoliveira@pm.me>
 * Copyright (C) 2024 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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
 * NcGalaxySDObsRedshift:
 *
 * Class describing galaxy sample observed redshift distribution.
 *
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "galaxy/nc_galaxy_sd_obs_redshift.h"
#include "math/ncm_vector.h"
#include "math/ncm_rng.h"

typedef struct _NcGalaxySDObsRedshiftPrivate
{
  gint placeholder;
} NcGalaxySDObsRedshiftPrivate;

enum
{
  PROP_0,
  PROP_LEN,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcGalaxySDObsRedshift, nc_galaxy_sd_obs_redshift, NCM_TYPE_MODEL);
G_DEFINE_BOXED_TYPE (NcGalaxySDObsRedshiftData, nc_galaxy_sd_obs_redshift_data, nc_galaxy_sd_obs_redshift_data_ref, nc_galaxy_sd_obs_redshift_data_unref);
NCM_UTIL_DEFINE_CALLBACK (NcGalaxySDObsRedshiftIntegrand,
                          NC_GALAXY_SD_OBS_REDSHIFT_INTEGRAND,
                          nc_galaxy_sd_obs_redshift_integrand,
                          gdouble,
                          NCM_UTIL_CALLBACK_ARGS (const gdouble z, NcGalaxySDObsRedshiftData * data),
                          NCM_UTIL_CALLBACK_ARGS (z, data))

static void
nc_galaxy_sd_obs_redshift_init (NcGalaxySDObsRedshift *gsdor)
{
}

/* LCOV_EXCL_START */
static void
_nc_galaxy_sd_obs_redshift_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcGalaxySDObsRedshift *gsdor = NC_GALAXY_SD_OBS_REDSHIFT (object);

  g_return_if_fail (NC_IS_GALAXY_SD_OBS_REDSHIFT (gsdor));

  switch (prop_id)
  {
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_galaxy_sd_obs_redshift_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcGalaxySDObsRedshift *gsdor = NC_GALAXY_SD_OBS_REDSHIFT (object);

  g_return_if_fail (NC_IS_GALAXY_SD_OBS_REDSHIFT (gsdor));

  switch (prop_id)
  {
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

/* LCOV_EXCL_STOP */

static void
_nc_galaxy_sd_obs_redshift_finalize (GObject *object)
{
  G_OBJECT_CLASS (nc_galaxy_sd_obs_redshift_parent_class)->finalize (object);
}

NCM_MSET_MODEL_REGISTER_ID (nc_galaxy_sd_obs_redshift, NC_TYPE_GALAXY_SD_OBS_REDSHIFT)

/*  LCOV_EXCL_START */
static void
_nc_galaxy_sd_obs_redshift_gen (NcGalaxySDObsRedshift *gsdor, NcGalaxySDObsRedshiftData *data, NcmRNG *rng)
{
  g_error ("_nc_galaxy_sd_obs_redshift_gen: method not implemented");
}

static gboolean
_nc_galaxy_sd_obs_redshift_gen1 (NcGalaxySDObsRedshift *gsdor, NcGalaxySDObsRedshiftData *data, NcmRNG *rng)
{
  g_error ("_nc_galaxy_sd_obs_redshift_gen1: method not implemented");

  return FALSE;
}

static void
_nc_galaxy_sd_obs_redshift_prepare (NcGalaxySDObsRedshift *gsdor, NcGalaxySDObsRedshiftData *data)
{
  g_error ("_nc_galaxy_sd_obs_redshift_prepare: method not implemented");
}

static void
_nc_galaxy_sd_obs_redshift_get_integ_lim (NcGalaxySDObsRedshift *gsdor, NcmMSet *mset, NcGalaxySDObsRedshiftData *data, gdouble *z_min, gdouble *z_max)
{
  g_error ("_nc_galaxy_sd_obs_redshift_get_integ_lim: method not implemented");
}

static NcGalaxySDObsRedshiftIntegrand *
_nc_galaxy_sd_obs_redshift_integ (NcGalaxySDObsRedshift *gsdor, gboolean use_lnp)
{
  g_error ("_nc_galaxy_sd_obs_redshift_integ: method not implemented");

  return NULL;
}

static void
_nc_galaxy_sd_obs_redshift_data_init (NcGalaxySDObsRedshift *gsdor, NcGalaxySDObsRedshiftData *data)
{
  g_error ("_nc_galaxy_sd_obs_redshift_data_new: method not implemented");
}

static NcmSpline *
_nc_galaxy_sd_obs_redshift_compute_binned_dndz (NcGalaxySDObsRedshift *gsdor, NcmVector *z_array)
{
  g_error ("_nc_galaxy_sd_obs_redshift_compute_binned_dndz: method not implemented");

  return NULL;
}

static NcmIntegralFixed *
_nc_galaxy_sd_obs_redshift_make_fixed_nodes (NcGalaxySDObsRedshift *gsdor, NcmMSet *mset,
                                             NcGalaxySDObsRedshiftData *data,
                                             gdouble z_lo, gdouble z_hi,
                                             guint n_nodes, guint rule_n)
{
  g_error ("_nc_galaxy_sd_obs_redshift_make_fixed_nodes: method not implemented");

  return NULL;
}

/*  LCOV_EXCL_STOP */

static void
nc_galaxy_sd_obs_redshift_class_init (NcGalaxySDObsRedshiftClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);
  NcmModelClass *model_class = NCM_MODEL_CLASS (klass);

  model_class->set_property = &_nc_galaxy_sd_obs_redshift_set_property;
  model_class->get_property = &_nc_galaxy_sd_obs_redshift_get_property;
  object_class->finalize    = &_nc_galaxy_sd_obs_redshift_finalize;

  ncm_model_class_set_name_nick (model_class, "Galaxy sample observed redshift distribution", "GalaxySDObsRedshift");
  ncm_model_class_add_params (model_class, 0, 0, PROP_LEN);
  ncm_mset_model_register_id (model_class, "NcGalaxySDObsRedshift", "Galaxy sample observed redshift distribution", NULL, FALSE, NCM_MSET_MODEL_MAIN);
  ncm_model_class_check_params_info (model_class);

  klass->gen                   = &_nc_galaxy_sd_obs_redshift_gen;
  klass->gen1                  = &_nc_galaxy_sd_obs_redshift_gen1;
  klass->prepare               = &_nc_galaxy_sd_obs_redshift_prepare;
  klass->get_integ_lim         = &_nc_galaxy_sd_obs_redshift_get_integ_lim;
  klass->integ                 = &_nc_galaxy_sd_obs_redshift_integ;
  klass->data_init             = &_nc_galaxy_sd_obs_redshift_data_init;
  klass->compute_binned_dndz   = &_nc_galaxy_sd_obs_redshift_compute_binned_dndz;
  klass->make_fixed_nodes      = &_nc_galaxy_sd_obs_redshift_make_fixed_nodes;
}

/**
 * nc_galaxy_sd_obs_redshift_data_ref:
 * @data: a #NcGalaxySDObsRedshiftData
 *
 * Increments the reference count of @data by one.
 *
 * Returns: (transfer full): a copy of @data
 */
NcGalaxySDObsRedshiftData *
nc_galaxy_sd_obs_redshift_data_ref (NcGalaxySDObsRedshiftData *data)
{
  g_atomic_ref_count_inc (&data->ref_count);

  return data;
}

/**
 * nc_galaxy_sd_obs_redshift_data_unref:
 * @data: a #NcGalaxySDObsRedshiftData
 *
 * Decreases the reference count of @data by one. If the reference count reaches 0, the
 * data is freed.
 */
void
nc_galaxy_sd_obs_redshift_data_unref (NcGalaxySDObsRedshiftData *data)
{
  if (g_atomic_ref_count_dec (&data->ref_count))
  {
    g_assert_nonnull (data->ldata_destroy);
    data->ldata_destroy (data->ldata);
    g_free (data);
  }
}

/**
 * nc_galaxy_sd_obs_redshift_data_read_row:
 * @data: a #NcGalaxySDObsRedshiftData
 * @obs: a #NcGalaxyWLObs
 * @i: the row index
 *
 * Reads the galaxy redshift data from the observation.
 *
 */
void
nc_galaxy_sd_obs_redshift_data_read_row (NcGalaxySDObsRedshiftData *data, NcGalaxyWLObs *obs, const guint i)
{
  data->z = nc_galaxy_wl_obs_get (obs, NC_GALAXY_SD_OBS_REDSHIFT_COL_Z, i, NULL);
  data->ldata_read_row (data, obs, i);
}

/**
 * nc_galaxy_sd_obs_redshift_data_write_row:
 * @data: a #NcGalaxySDObsRedshiftData
 * @obs: a #NcGalaxyWLObs
 * @i: the row index
 *
 * Writes the galaxy redshift data to the observation.
 *
 */
void
nc_galaxy_sd_obs_redshift_data_write_row (NcGalaxySDObsRedshiftData *data, NcGalaxyWLObs *obs, const guint i)
{
  nc_galaxy_wl_obs_set (obs, NC_GALAXY_SD_OBS_REDSHIFT_COL_Z, i, data->z, NULL);
  data->ldata_write_row (data, obs, i);
}

/**
 * nc_galaxy_sd_obs_redshift_data_required_columns:
 * @data: a #NcGalaxySDObsRedshiftData
 *
 * Returns: (element-type utf8) (transfer full): the required columns for the galaxy redshift data.
 */
GList *
nc_galaxy_sd_obs_redshift_data_required_columns (NcGalaxySDObsRedshiftData *data)
{
  GList *columns = NULL;

  columns = g_list_append (columns, g_strdup (NC_GALAXY_SD_OBS_REDSHIFT_COL_Z));
  data->ldata_required_columns (data, columns);

  return columns;
}

/**
 * nc_galaxy_sd_obs_redshift_integrand_new:
 * @func: (scope async) (closure callback_data): a #NcGalaxySDObsRedshiftIntegrandFunc
 * @callback_data_free: (scope async) (closure callback_data): a #NcGalaxySDObsRedshiftIntegrandFreeData
 * @callback_data_copy: (scope async) (closure callback_data): a #NcGalaxySDObsRedshiftIntegrandCopyData
 * @callback_data_prepare: (scope async) (closure callback_data): a #NcGalaxySDObsRedshiftIntegrandPrepareData
 * @callback_data: a gpointer
 *
 * Creates a new integrand for the galaxy redshift data.
 * The integrand is a function that takes the redshift @z and the galaxy redshift data @data as arguments.
 * The function should return the integrand value at @z.
 *
 * Returns: (transfer full): a new #NcGalaxySDObsRedshiftIntegrand object.
 */
/**
 * nc_galaxy_sd_obs_redshift_integrand_copy:
 * @callback_obj: a #NcGalaxySDObsRedshiftIntegrand
 *
 * Copies the integrand for the galaxy redshift data.
 *
 * Returns: (transfer full): a copy of @callback_obj
 */
/**
 * nc_galaxy_sd_obs_redshift_integrand_free:
 * @callback_obj: a #NcGalaxySDObsRedshiftIntegrand
 *
 * Frees the integrand for the galaxy redshift data.
 *
 */
/**
 * nc_galaxy_sd_obs_redshift_integrand_prepare:
 * @callback_obj: a #NcGalaxySDObsRedshiftIntegrand
 * @mset: a #NcmMSet
 *
 * Prepares the integrand for the galaxy redshift data.
 *
 */

/**
 * nc_galaxy_sd_obs_redshift_ref:
 * @gsdor: a #NcGalaxySDObsRedshift
 *
 * Increases the reference count of @gsdor by one.
 *
 * Returns: (transfer full): @gsdor
 */
NcGalaxySDObsRedshift *
nc_galaxy_sd_obs_redshift_ref (NcGalaxySDObsRedshift *gsdor)
{
  return g_object_ref (gsdor);
}

/**
 * nc_galaxy_sd_obs_redshift_free:
 * @gsdor: a #NcGalaxySDObsRedshift
 *
 * Decreases the reference count of @gsdor by one.
 *
 */
void
nc_galaxy_sd_obs_redshift_free (NcGalaxySDObsRedshift *gsdor)
{
  g_object_unref (gsdor);
}

/**
 * nc_galaxy_sd_obs_redshift_clear:
 * @gsdor: a #NcGalaxySDObsRedshift
 *
 * Decreases the reference count of @gsdor by one, and sets the pointer *@gsdor to
 * NULL.
 *
 */
void
nc_galaxy_sd_obs_redshift_clear (NcGalaxySDObsRedshift **gsdor)
{
  g_clear_object (gsdor);
}

/**
 * nc_galaxy_sd_obs_redshift_gen:
 * @gsdor: a #NcGalaxySDObsRedshift
 * @data: a #NcGalaxySDObsRedshiftData
 * @rng: a #NcmRNG
 *
 * Generates a new galaxy redshift data. The #NcGalaxySDObsRedshiftData object @data must be
 * initialized before calling this method.
 *
 */
void
nc_galaxy_sd_obs_redshift_gen (NcGalaxySDObsRedshift *gsdor, NcGalaxySDObsRedshiftData *data, NcmRNG *rng)
{
  NC_GALAXY_SD_OBS_REDSHIFT_GET_CLASS (gsdor)->gen (gsdor, data, rng);
}

/**
 * nc_galaxy_sd_obs_redshift_gen1:
 * @gsdor: a #NcGalaxySDObsRedshift instance
 * @data: a pre-initialized #NcGalaxySDObsRedshiftData
 * @rng: a #NcmRNG random number generator
 *
 * Attempts to generate a single redshift sample consistent with the observational
 * constraints defined in @gsdor. The result is stored in @data.
 *
 * This method is typically used in scenarios where the total number of galaxies is
 * fixed, and we wish to construct a subsample that satisfies observational selection
 * criteria (e.g., redshift cuts or survey limitations). For each galaxy in the total
 * sample, a redshift is proposed, and if it violates the constraints, the galaxy is
 * discarded from the final subsample.
 *
 * This sampling approach avoids the need to compute the normalization or acceptance
 * fraction analytically or numerically, which may involve complex or model-dependent
 * integrals. It also offers flexibility for implementing more general or evolving
 * selection functions.
 *
 * Returns: %TRUE if a valid redshift was generated; %FALSE otherwise.
 */
gboolean
nc_galaxy_sd_obs_redshift_gen1 (NcGalaxySDObsRedshift *gsdor, NcGalaxySDObsRedshiftData *data, NcmRNG *rng)
{
  return NC_GALAXY_SD_OBS_REDSHIFT_GET_CLASS (gsdor)->gen1 (gsdor, data, rng);
}

/**
 * nc_galaxy_sd_obs_redshift_prepare:
 * @gsdor: a #NcGalaxySDObsRedshift
 * @data: a #NcGalaxySDObsRedshiftData
 *
 * Prepares the galaxy redshift data for generation.
 *
 */
void
nc_galaxy_sd_obs_redshift_prepare (NcGalaxySDObsRedshift *gsdor, NcGalaxySDObsRedshiftData *data)
{
  NC_GALAXY_SD_OBS_REDSHIFT_GET_CLASS (gsdor)->prepare (gsdor, data);
}

/**
 * nc_galaxy_sd_obs_redshift_get_integ_lim:
 * @gsdor: a #NcGalaxySDObsRedshift
 * @mset: a #NcmMSet
 * @data: a #NcGalaxySDObsRedshiftData
 * @z_min: (out): the minimum redshift for integration
 * @z_max: (out): the maximum redshift for integration
 *
 * Gets the effective redshift integration support for the galaxy redshift data:
 * the range over which the per-galaxy P(z) factor is non-negligible (e.g. the
 * photometric kernel restricted to its relevant range). This is the domain used
 * by every integration method - fixed Gauss-Legendre nodes as well as the
 * adaptive quadratures.
 *
 */
void
nc_galaxy_sd_obs_redshift_get_integ_lim (NcGalaxySDObsRedshift *gsdor, NcmMSet *mset, NcGalaxySDObsRedshiftData *data, gdouble *z_min, gdouble *z_max)
{
  NC_GALAXY_SD_OBS_REDSHIFT_GET_CLASS (gsdor)->get_integ_lim (gsdor, mset, data, z_min, z_max);
}

/**
 * nc_galaxy_sd_obs_redshift_integ:
 * @gsdor: a #NcGalaxySDObsRedshift
 * @use_lnp: if TRUE the integrand must return the natural logarithm of the probability density
 *
 * Prepares the integrand for the galaxy redshift data.
 *
 */
NcGalaxySDObsRedshiftIntegrand *
nc_galaxy_sd_obs_redshift_integ (NcGalaxySDObsRedshift *gsdor, gboolean use_lnp)
{
  return NC_GALAXY_SD_OBS_REDSHIFT_GET_CLASS (gsdor)->integ (gsdor, use_lnp);
}

/**
 * nc_galaxy_sd_obs_redshift_compute_binned_dndz:
 * @gsdor: a #NcGalaxySDObsRedshift
 * @z_array: (array) (element-type gdouble) (nullable): true redshift evaluation points
 *
 * Computes the binned true redshift distribution dndz(z) for the photometric redshift
 * bin defined by the observation model. The photo-z bin edges are specified when
 * creating the observation model (e.g., via nc_galaxy_sd_obs_redshift_gauss_new()).
 *
 * The returned distribution is normalized such that $\int \mathrm{d}n/\mathrm{d}z \,
 * \mathrm{d}z = 1$.
 *
 * This method works for both lens-type bins (fixed photo-z edges) and
 * source-type bins (equal-area photo-z edges). The distinction is only in
 * how the photo-z bin edges were initially determined.
 *
 * If @z_array is provided, the binned dndz will be evaluated at those redshift points.
 * If @z_array is NULL, the method will determine an appropriate set of redshift points
 * based on the integration limits and the desired resolution.
 *
 * Returns: (transfer full): a #NcmSpline containing the binned dndz(z)
 */
NcmSpline *
nc_galaxy_sd_obs_redshift_compute_binned_dndz (NcGalaxySDObsRedshift *gsdor, NcmVector *z_array)
{
  return NC_GALAXY_SD_OBS_REDSHIFT_GET_CLASS (gsdor)->compute_binned_dndz (gsdor, z_array);
}

/**
 * nc_galaxy_sd_obs_redshift_make_fixed_nodes: (skip)
 * @gsdor: a #NcGalaxySDObsRedshift
 * @mset: a #NcmMSet
 * @data: a #NcGalaxySDObsRedshiftData
 * @z_lo: lower integration limit
 * @z_hi: upper integration limit
 * @n_nodes: number of Gauss-Legendre intervals
 * @rule_n: number of GL points per interval
 *
 * Allocates and populates a #NcmIntegralFixed over [@z_lo, @z_hi] so that the
 * precomputed node weights absorb the per-galaxy P(z) factor. The returned
 * object is owned by the caller. Callers may build several panels over adjacent
 * sub-ranges so that a known non-smooth point (e.g. the lens redshift) falls on
 * a panel boundary rather than inside a Gauss-Legendre interval.
 *
 * Returns: (transfer full): a new #NcmIntegralFixed with P(z)-weighted nodes.
 */
NcmIntegralFixed *
nc_galaxy_sd_obs_redshift_make_fixed_nodes (NcGalaxySDObsRedshift *gsdor, NcmMSet *mset,
                                            NcGalaxySDObsRedshiftData *data,
                                            gdouble z_lo, gdouble z_hi,
                                            guint n_nodes, guint rule_n)
{
  return NC_GALAXY_SD_OBS_REDSHIFT_GET_CLASS (gsdor)->make_fixed_nodes (gsdor, mset, data, z_lo, z_hi, n_nodes, rule_n);
}

/**
 * nc_galaxy_sd_obs_redshift_data_new:
 * @gsdor: a #NcGalaxySDObsRedshift
 *
 * Creates a new galaxy redshift data.
 *
 * Returns: (transfer full): a new #NcGalaxySDObsRedshiftData object.
 */
NcGalaxySDObsRedshiftData *
nc_galaxy_sd_obs_redshift_data_new (NcGalaxySDObsRedshift *gsdor)
{
  NcGalaxySDObsRedshiftData *data = g_new0 (NcGalaxySDObsRedshiftData, 1);

  data->z                      = 0.0;
  data->ldata                  = NULL;
  data->ldata_destroy          = NULL;
  data->ldata_read_row         = NULL;
  data->ldata_write_row        = NULL;
  data->ldata_required_columns = NULL;

  g_atomic_ref_count_init (&data->ref_count);
  NC_GALAXY_SD_OBS_REDSHIFT_GET_CLASS (gsdor)->data_init (gsdor, data);

  g_assert_nonnull (data->ldata_destroy);
  g_assert_nonnull (data->ldata_read_row);
  g_assert_nonnull (data->ldata_write_row);
  g_assert_nonnull (data->ldata_required_columns);

  return data;
}


/***************************************************************************
 *            nc_galaxy_wl_obs_binned.c
 *
 *  Thu November 26 2026
 *  Copyright  2026  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * nc_galaxy_wl_obs_binned.c
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
 * NcGalaxyWLObsBinned:
 *
 * Binned weak-lensing profile of a cluster.
 *
 * A #NcmCatalog with one row per radial bin. Every binned profile has the columns
 * #NC_GALAXY_WL_OBS_BINNED_COL_R_LO and #NC_GALAXY_WL_OBS_BINNED_COL_R_HI (bin edges,
 * proper Mpc), #NC_GALAXY_WL_OBS_BINNED_COL_RADIUS (the radius at which the model is
 * evaluated), #NC_GALAXY_WL_OBS_BINNED_COL_VALUE (the measured profile),
 * #NC_GALAXY_WL_OBS_BINNED_COL_OFFSET (the data-only term of the expected value of the
 * estimator, see below) and #NC_GALAXY_WL_OBS_BINNED_COL_N_GAL (the number of galaxies of
 * the bin). An estimator adds its own columns, the bin-averaged quantities its classic model
 * needs (for instance the mean inverse critical density of #NcGalaxyWLBinnedEstimatorDeltaSigma).
 *
 * The object also holds the covariance of the measured values
 * (#NcGalaxyWLObsBinned:cov), the cluster redshift used to build it and the name of the
 * estimator that built it.
 *
 * A profile built from a #NcGalaxyWLObs by a #NcGalaxyWLBinnedEstimator also holds the
 * source catalog and a members catalog, with one row per galaxy that entered a bin: its row
 * in the source catalog, its bin, its projected radius and the factor $f_i$ such that the
 * expected value of the estimator in bin $b$ is
 * $$\langle\hat d_b\rangle = o_b + \sum_{i \in b} f_i\, \bar g_i,$$
 * where $\bar g_i$ is the reduced shear of galaxy $i$ averaged over its redshift density and
 * $o_b$ the offset column. A profile built from external arrays has neither.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc/lss/galaxy/nc_galaxy_wl_obs_binned.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_math.h>
#endif /* NUMCOSMO_GIR_SCAN */

struct _NcGalaxyWLObsBinned
{
  /*< private >*/
  NcmCatalog parent_instance;
  NcmMatrix *cov;
  gdouble z_cluster;
  gchar *estimator;
  NcGalaxyWLObs *source;
  NcmCatalog *members;
};

enum
{
  PROP_COV = 1000,
  PROP_Z_CLUSTER,
  PROP_ESTIMATOR,
  PROP_SOURCE,
  PROP_MEMBERS,
};

G_DEFINE_TYPE (NcGalaxyWLObsBinned, nc_galaxy_wl_obs_binned, NCM_TYPE_CATALOG)

static void
nc_galaxy_wl_obs_binned_init (NcGalaxyWLObsBinned *binned)
{
  binned->cov       = NULL;
  binned->z_cluster = 0.0;
  binned->estimator = NULL;
  binned->source    = NULL;
  binned->members   = NULL;
}

static void
_nc_galaxy_wl_obs_binned_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcGalaxyWLObsBinned *binned = NC_GALAXY_WL_OBS_BINNED (object);

  switch (prop_id)
  {
    case PROP_COV:
      ncm_matrix_clear (&binned->cov);
      binned->cov = g_value_dup_object (value);
      break;
    case PROP_Z_CLUSTER:
      binned->z_cluster = g_value_get_double (value);
      break;
    case PROP_ESTIMATOR:
      nc_galaxy_wl_obs_binned_set_estimator (binned, g_value_get_string (value));
      break;
    case PROP_SOURCE:
      nc_galaxy_wl_obs_clear (&binned->source);
      binned->source = g_value_dup_object (value);
      break;
    case PROP_MEMBERS:
      ncm_catalog_clear (&binned->members);
      binned->members = g_value_dup_object (value);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_galaxy_wl_obs_binned_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcGalaxyWLObsBinned *binned = NC_GALAXY_WL_OBS_BINNED (object);

  switch (prop_id)
  {
    case PROP_COV:
      g_value_set_object (value, binned->cov);
      break;
    case PROP_Z_CLUSTER:
      g_value_set_double (value, binned->z_cluster);
      break;
    case PROP_ESTIMATOR:
      g_value_set_string (value, binned->estimator);
      break;
    case PROP_SOURCE:
      g_value_set_object (value, binned->source);
      break;
    case PROP_MEMBERS:
      g_value_set_object (value, binned->members);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_galaxy_wl_obs_binned_dispose (GObject *object)
{
  NcGalaxyWLObsBinned *binned = NC_GALAXY_WL_OBS_BINNED (object);

  ncm_matrix_clear (&binned->cov);
  nc_galaxy_wl_obs_clear (&binned->source);
  ncm_catalog_clear (&binned->members);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_galaxy_wl_obs_binned_parent_class)->dispose (object);
}

static void
_nc_galaxy_wl_obs_binned_finalize (GObject *object)
{
  NcGalaxyWLObsBinned *binned = NC_GALAXY_WL_OBS_BINNED (object);

  g_clear_pointer (&binned->estimator, g_free);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_galaxy_wl_obs_binned_parent_class)->finalize (object);
}

static void
nc_galaxy_wl_obs_binned_class_init (NcGalaxyWLObsBinnedClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = &_nc_galaxy_wl_obs_binned_set_property;
  object_class->get_property = &_nc_galaxy_wl_obs_binned_get_property;
  object_class->dispose      = &_nc_galaxy_wl_obs_binned_dispose;
  object_class->finalize     = &_nc_galaxy_wl_obs_binned_finalize;

  /**
   * NcGalaxyWLObsBinned:cov:
   *
   * The covariance of the measured values (#NC_GALAXY_WL_OBS_BINNED_COL_VALUE).
   */
  g_object_class_install_property (object_class,
                                   PROP_COV,
                                   g_param_spec_object ("cov",
                                                        NULL,
                                                        "Covariance of the binned values",
                                                        NCM_TYPE_MATRIX,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcGalaxyWLObsBinned:z-cluster:
   *
   * The cluster redshift used to build the profile.
   */
  g_object_class_install_property (object_class,
                                   PROP_Z_CLUSTER,
                                   g_param_spec_double ("z-cluster",
                                                        NULL,
                                                        "Cluster redshift",
                                                        0.0, G_MAXDOUBLE, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcGalaxyWLObsBinned:estimator:
   *
   * The name of the estimator that built the profile, or %NULL for a profile built from
   * external arrays.
   */
  g_object_class_install_property (object_class,
                                   PROP_ESTIMATOR,
                                   g_param_spec_string ("estimator",
                                                        NULL,
                                                        "Name of the estimator",
                                                        NULL,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcGalaxyWLObsBinned:source:
   *
   * The per-galaxy catalog the profile was built from, or %NULL.
   */
  g_object_class_install_property (object_class,
                                   PROP_SOURCE,
                                   g_param_spec_object ("source",
                                                        NULL,
                                                        "Source per-galaxy catalog",
                                                        NC_TYPE_GALAXY_WL_OBS,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcGalaxyWLObsBinned:members:
   *
   * The members catalog: one row per galaxy of the source catalog that entered a bin, or
   * %NULL.
   */
  g_object_class_install_property (object_class,
                                   PROP_MEMBERS,
                                   g_param_spec_object ("members",
                                                        NULL,
                                                        "Galaxies of each bin",
                                                        NCM_TYPE_CATALOG,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

/**
 * nc_galaxy_wl_obs_binned_new:
 * @n_bins: the number of bins
 * @extra_cols: (array zero-terminated=1) (nullable): the estimator columns
 *
 * Creates a new #NcGalaxyWLObsBinned with @n_bins rows, the standard columns and
 * @extra_cols, all set to zero.
 *
 * Returns: (transfer full): a new #NcGalaxyWLObsBinned.
 */
NcGalaxyWLObsBinned *
nc_galaxy_wl_obs_binned_new (guint n_bins, GStrv extra_cols)
{
  GPtrArray *cols = g_ptr_array_new ();
  NcGalaxyWLObsBinned *binned;

  g_ptr_array_add (cols, (gpointer) NC_GALAXY_WL_OBS_BINNED_COL_R_LO);
  g_ptr_array_add (cols, (gpointer) NC_GALAXY_WL_OBS_BINNED_COL_R_HI);
  g_ptr_array_add (cols, (gpointer) NC_GALAXY_WL_OBS_BINNED_COL_RADIUS);
  g_ptr_array_add (cols, (gpointer) NC_GALAXY_WL_OBS_BINNED_COL_VALUE);
  g_ptr_array_add (cols, (gpointer) NC_GALAXY_WL_OBS_BINNED_COL_OFFSET);
  g_ptr_array_add (cols, (gpointer) NC_GALAXY_WL_OBS_BINNED_COL_N_GAL);

  if (extra_cols != NULL)
  {
    guint i;

    for (i = 0; extra_cols[i] != NULL; i++)
      g_ptr_array_add (cols, extra_cols[i]);
  }

  g_ptr_array_add (cols, NULL);

  binned = g_object_new (NC_TYPE_GALAXY_WL_OBS_BINNED,
                         "columns", (GStrv) cols->pdata,
                         "len", n_bins,
                         NULL);

  g_ptr_array_unref (cols);

  {
    NcmMatrix *cov = ncm_matrix_new (n_bins, n_bins);

    ncm_matrix_set_zero (cov);
    nc_galaxy_wl_obs_binned_set_cov (binned, cov);
    ncm_matrix_free (cov);
  }

  return binned;
}

/**
 * nc_galaxy_wl_obs_binned_new_from_arrays:
 * @radius: (element-type gdouble): the radii of the bins (proper Mpc)
 * @value: (element-type gdouble): the measured profile
 * @cov: the covariance of @value
 * @extra_cols: (array zero-terminated=1) (nullable): the estimator columns
 *
 * Creates a binned profile from external arrays, for instance a published profile. The bin
 * edges and the number of galaxies are set to NaN and the offsets to zero; the estimator
 * columns must be filled with nc_galaxy_wl_obs_binned_set_column(), and the estimator name
 * with nc_galaxy_wl_obs_binned_set_estimator(). Such a profile has no members, so it supports
 * only the classic model of the estimator.
 *
 * Returns: (transfer full): a new #NcGalaxyWLObsBinned.
 */
NcGalaxyWLObsBinned *
nc_galaxy_wl_obs_binned_new_from_arrays (GArray *radius, GArray *value, NcmMatrix *cov, GStrv extra_cols)
{
  const guint n_bins = radius->len;
  NcGalaxyWLObsBinned *binned;
  guint i;

  g_assert_cmpuint (value->len, ==, n_bins);
  g_assert_cmpuint (ncm_matrix_nrows (cov), ==, n_bins);
  g_assert_cmpuint (ncm_matrix_ncols (cov), ==, n_bins);

  binned = nc_galaxy_wl_obs_binned_new (n_bins, extra_cols);

  for (i = 0; i < n_bins; i++)
  {
    nc_galaxy_wl_obs_binned_set (binned, NC_GALAXY_WL_OBS_BINNED_COL_R_LO, i, GSL_NAN);
    nc_galaxy_wl_obs_binned_set (binned, NC_GALAXY_WL_OBS_BINNED_COL_R_HI, i, GSL_NAN);
    nc_galaxy_wl_obs_binned_set (binned, NC_GALAXY_WL_OBS_BINNED_COL_N_GAL, i, GSL_NAN);
  }

  nc_galaxy_wl_obs_binned_set_column (binned, NC_GALAXY_WL_OBS_BINNED_COL_RADIUS, radius);
  nc_galaxy_wl_obs_binned_set_column (binned, NC_GALAXY_WL_OBS_BINNED_COL_VALUE, value);

  {
    NcmMatrix *cov_dup = ncm_matrix_dup (cov);

    nc_galaxy_wl_obs_binned_set_cov (binned, cov_dup);
    ncm_matrix_free (cov_dup);
  }

  return binned;
}

/**
 * nc_galaxy_wl_obs_binned_ref:
 * @binned: a #NcGalaxyWLObsBinned
 *
 * Increases the reference count of @binned by one.
 *
 * Returns: (transfer full): @binned.
 */
NcGalaxyWLObsBinned *
nc_galaxy_wl_obs_binned_ref (NcGalaxyWLObsBinned *binned)
{
  return g_object_ref (binned);
}

/**
 * nc_galaxy_wl_obs_binned_free:
 * @binned: a #NcGalaxyWLObsBinned
 *
 * Decreases the reference count of @binned by one.
 */
void
nc_galaxy_wl_obs_binned_free (NcGalaxyWLObsBinned *binned)
{
  g_object_unref (binned);
}

/**
 * nc_galaxy_wl_obs_binned_clear:
 * @binned: a #NcGalaxyWLObsBinned
 *
 * Decreases the reference count of *@binned by one and sets *@binned to %NULL.
 */
void
nc_galaxy_wl_obs_binned_clear (NcGalaxyWLObsBinned **binned)
{
  g_clear_object (binned);
}

/**
 * nc_galaxy_wl_obs_binned_len:
 * @binned: a #NcGalaxyWLObsBinned
 *
 * Returns: the number of bins.
 */
guint
nc_galaxy_wl_obs_binned_len (NcGalaxyWLObsBinned *binned)
{
  return ncm_catalog_len (NCM_CATALOG (binned));
}

/**
 * nc_galaxy_wl_obs_binned_set:
 * @binned: a #NcGalaxyWLObsBinned
 * @col: the column name
 * @i: the bin index
 * @val: the value
 *
 * Sets the value of column @col in bin @i; an unknown column aborts.
 */
void
nc_galaxy_wl_obs_binned_set (NcGalaxyWLObsBinned *binned, const gchar *col, const guint i, const gdouble val)
{
  ncm_catalog_set (NCM_CATALOG (binned), col, i, val, NULL);
}

/**
 * nc_galaxy_wl_obs_binned_get:
 * @binned: a #NcGalaxyWLObsBinned
 * @col: the column name
 * @i: the bin index
 *
 * Returns: the value of column @col in bin @i; an unknown column aborts.
 */
gdouble
nc_galaxy_wl_obs_binned_get (NcGalaxyWLObsBinned *binned, const gchar *col, const guint i)
{
  return ncm_catalog_get (NCM_CATALOG (binned), col, i, NULL);
}

/**
 * nc_galaxy_wl_obs_binned_has_column:
 * @binned: a #NcGalaxyWLObsBinned
 * @col: the column name
 *
 * Returns: whether @binned has the column @col.
 */
gboolean
nc_galaxy_wl_obs_binned_has_column (NcGalaxyWLObsBinned *binned, const gchar *col)
{
  return ncm_catalog_has_column (NCM_CATALOG (binned), col);
}

/**
 * nc_galaxy_wl_obs_binned_set_column:
 * @binned: a #NcGalaxyWLObsBinned
 * @col: the column name
 * @v: (element-type gdouble): the values, one per bin
 *
 * Sets column @col of every bin.
 */
void
nc_galaxy_wl_obs_binned_set_column (NcGalaxyWLObsBinned *binned, const gchar *col, GArray *v)
{
  const guint n_bins = nc_galaxy_wl_obs_binned_len (binned);
  guint i;

  g_assert_cmpuint (v->len, ==, n_bins);

  for (i = 0; i < n_bins; i++)
    nc_galaxy_wl_obs_binned_set (binned, col, i, g_array_index (v, gdouble, i));
}

/**
 * nc_galaxy_wl_obs_binned_get_column:
 * @binned: a #NcGalaxyWLObsBinned
 * @col: the column name
 *
 * Returns: (transfer full): a new vector with column @col of every bin.
 */
NcmVector *
nc_galaxy_wl_obs_binned_get_column (NcGalaxyWLObsBinned *binned, const gchar *col)
{
  const guint n_bins = nc_galaxy_wl_obs_binned_len (binned);
  NcmVector *v       = ncm_vector_new (n_bins);
  guint i;

  for (i = 0; i < n_bins; i++)
    ncm_vector_set (v, i, nc_galaxy_wl_obs_binned_get (binned, col, i));

  return v;
}

/**
 * nc_galaxy_wl_obs_binned_set_cov:
 * @binned: a #NcGalaxyWLObsBinned
 * @cov: the covariance of the binned values
 *
 * Sets the covariance of the binned values, referencing @cov.
 */
void
nc_galaxy_wl_obs_binned_set_cov (NcGalaxyWLObsBinned *binned, NcmMatrix *cov)
{
  const guint n_bins = nc_galaxy_wl_obs_binned_len (binned);

  g_assert_cmpuint (ncm_matrix_nrows (cov), ==, n_bins);
  g_assert_cmpuint (ncm_matrix_ncols (cov), ==, n_bins);

  ncm_matrix_clear (&binned->cov);
  binned->cov = ncm_matrix_ref (cov);
}

/**
 * nc_galaxy_wl_obs_binned_peek_cov:
 * @binned: a #NcGalaxyWLObsBinned
 *
 * Returns: (transfer none): the covariance of the binned values.
 */
NcmMatrix *
nc_galaxy_wl_obs_binned_peek_cov (NcGalaxyWLObsBinned *binned)
{
  return binned->cov;
}

/**
 * nc_galaxy_wl_obs_binned_set_z_cluster:
 * @binned: a #NcGalaxyWLObsBinned
 * @z_cluster: the cluster redshift
 *
 * Sets the cluster redshift used to build the profile.
 */
void
nc_galaxy_wl_obs_binned_set_z_cluster (NcGalaxyWLObsBinned *binned, const gdouble z_cluster)
{
  binned->z_cluster = z_cluster;
}

/**
 * nc_galaxy_wl_obs_binned_get_z_cluster:
 * @binned: a #NcGalaxyWLObsBinned
 *
 * Returns: the cluster redshift used to build the profile.
 */
gdouble
nc_galaxy_wl_obs_binned_get_z_cluster (NcGalaxyWLObsBinned *binned)
{
  return binned->z_cluster;
}

/**
 * nc_galaxy_wl_obs_binned_set_estimator:
 * @binned: a #NcGalaxyWLObsBinned
 * @estimator: (nullable): the name of the estimator
 *
 * Sets the name of the estimator of the profile.
 */
void
nc_galaxy_wl_obs_binned_set_estimator (NcGalaxyWLObsBinned *binned, const gchar *estimator)
{
  g_clear_pointer (&binned->estimator, g_free);
  binned->estimator = g_strdup (estimator);
}

/**
 * nc_galaxy_wl_obs_binned_peek_estimator:
 * @binned: a #NcGalaxyWLObsBinned
 *
 * Returns: (transfer none) (nullable): the name of the estimator of the profile.
 */
const gchar *
nc_galaxy_wl_obs_binned_peek_estimator (NcGalaxyWLObsBinned *binned)
{
  return binned->estimator;
}

/**
 * nc_galaxy_wl_obs_binned_set_members:
 * @binned: a #NcGalaxyWLObsBinned
 * @source: (nullable): the per-galaxy catalog of the profile
 * @members: (nullable): the members catalog
 *
 * Sets the source catalog and the members catalog (see the class description), or removes
 * them when both are %NULL.
 */
void
nc_galaxy_wl_obs_binned_set_members (NcGalaxyWLObsBinned *binned, NcGalaxyWLObs *source, NcmCatalog *members)
{
  g_assert ((source == NULL) == (members == NULL));

  nc_galaxy_wl_obs_clear (&binned->source);
  ncm_catalog_clear (&binned->members);

  if (source != NULL)
  {
    binned->source  = nc_galaxy_wl_obs_ref (source);
    binned->members = ncm_catalog_ref (members);
  }
}

/**
 * nc_galaxy_wl_obs_binned_peek_source:
 * @binned: a #NcGalaxyWLObsBinned
 *
 * Returns: (transfer none) (nullable): the per-galaxy catalog of the profile.
 */
NcGalaxyWLObs *
nc_galaxy_wl_obs_binned_peek_source (NcGalaxyWLObsBinned *binned)
{
  return binned->source;
}

/**
 * nc_galaxy_wl_obs_binned_peek_members:
 * @binned: a #NcGalaxyWLObsBinned
 *
 * Returns: (transfer none) (nullable): the members catalog of the profile.
 */
NcmCatalog *
nc_galaxy_wl_obs_binned_peek_members (NcGalaxyWLObsBinned *binned)
{
  return binned->members;
}

/**
 * nc_galaxy_wl_obs_binned_has_members:
 * @binned: a #NcGalaxyWLObsBinned
 *
 * Returns: whether the profile holds its source and members catalogs.
 */
gboolean
nc_galaxy_wl_obs_binned_has_members (NcGalaxyWLObsBinned *binned)
{
  return (binned->source != NULL) && (binned->members != NULL);
}


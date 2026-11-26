/***************************************************************************
 *            nc_galaxy_wl_obs_cheb_grid.c
 *
 *  Thu November 26 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_galaxy_wl_obs_cheb_grid.c
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
 * NcGalaxyWLObsChebGrid:
 *
 * The CHEB_WEIGHTS redshift grid of the galaxies of a #NcGalaxyWLObs, with the record of
 * what built it.
 *
 * The grid of a galaxy is what #NcDataClusterWLFactor computes under
 * %NC_DATA_CLUSTER_WL_INTEG_METHOD_CHEB_WEIGHTS before any evaluation: the weights of the
 * nested Chebyshev levels, the redshift values (the lower end of the p(z) support followed
 * by the Chebyshev nodes) and the normalizations of p(z) (total, background and
 * foreground). It depends on the p(z) of the galaxy, on its integration limits and on the
 * cluster redshift only, so it is computed once, stored with the catalog
 * (#NcGalaxyWLObs:cheb-grid) and used by every dataset on that catalog.
 *
 * The record holds the number of galaxies, the cluster redshift, the type name of the
 * redshift factor that built the grid and the data hash of the rows
 * (ncm_catalog_get_data_hash()). Each galaxy has a row of #NcGalaxyWLObsChebGrid:n-weights
 * weights and one of #NcGalaxyWLObsChebGrid:n-nodes redshift values; a galaxy with
 * $n_\mathrm{total}$ nodes after the lower end uses the first $n_\mathrm{total} + 1$ values,
 * and its weights only when $n_\mathrm{total} > 0$.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc/lss/galaxy/nc_galaxy_wl_obs_cheb_grid.h"
#include "ncm/algebra/ncm_matrix.h"
#include "ncm/core/ncm_cfg.h"

struct _NcGalaxyWLObsChebGrid
{
  /*< private >*/
  GObject parent_instance;
  guint len;
  guint n_weights;
  guint n_nodes;
  gdouble z_cl;
  gchar *z_type;
  gchar *data_hash;
  NcmMatrix *weights; /* len x n_weights */
  NcmMatrix *z_nodes; /* len x n_nodes */
  NcmMatrix *norms;   /* len x 3: total, background, foreground */
  GArray *n_total;    /* guint, len */
  gboolean constructed;
};

enum
{
  PROP_0,
  PROP_LEN,
  PROP_N_WEIGHTS,
  PROP_N_NODES,
  PROP_Z_CL,
  PROP_Z_TYPE,
  PROP_DATA_HASH,
  PROP_WEIGHTS,
  PROP_Z_NODES,
  PROP_NORMS,
  PROP_N_TOTAL,
};

G_DEFINE_TYPE (NcGalaxyWLObsChebGrid, nc_galaxy_wl_obs_cheb_grid, G_TYPE_OBJECT)

static void
nc_galaxy_wl_obs_cheb_grid_init (NcGalaxyWLObsChebGrid *grid)
{
  grid->len         = 0;
  grid->n_weights   = 0;
  grid->n_nodes     = 0;
  grid->z_cl        = 0.0;
  grid->z_type      = NULL;
  grid->data_hash   = NULL;
  grid->weights     = NULL;
  grid->z_nodes     = NULL;
  grid->norms       = NULL;
  grid->n_total     = NULL;
  grid->constructed = FALSE;
}

static void _nc_galaxy_wl_obs_cheb_grid_set_matrix (NcGalaxyWLObsChebGrid *grid, NcmMatrix **slot, NcmMatrix *m, const guint ncols, const gchar *name);

static void
_nc_galaxy_wl_obs_cheb_grid_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcGalaxyWLObsChebGrid *grid = NC_GALAXY_WL_OBS_CHEB_GRID (object);

  switch (prop_id)
  {
    case PROP_LEN:
      grid->len = g_value_get_uint (value);
      break;
    case PROP_N_WEIGHTS:
      grid->n_weights = g_value_get_uint (value);
      break;
    case PROP_N_NODES:
      grid->n_nodes = g_value_get_uint (value);
      break;
    case PROP_Z_CL:
      grid->z_cl = g_value_get_double (value);
      break;
    case PROP_Z_TYPE:
      g_clear_pointer (&grid->z_type, g_free);
      grid->z_type = g_value_dup_string (value);
      break;
    case PROP_DATA_HASH:
      g_clear_pointer (&grid->data_hash, g_free);
      grid->data_hash = g_value_dup_string (value);
      break;
    case PROP_WEIGHTS:
      _nc_galaxy_wl_obs_cheb_grid_set_matrix (grid, &grid->weights, g_value_get_object (value), grid->n_weights, "weights");
      break;
    case PROP_Z_NODES:
      _nc_galaxy_wl_obs_cheb_grid_set_matrix (grid, &grid->z_nodes, g_value_get_object (value), grid->n_nodes, "z-nodes");
      break;
    case PROP_NORMS:
      _nc_galaxy_wl_obs_cheb_grid_set_matrix (grid, &grid->norms, g_value_get_object (value), 3, "norms");
      break;
    case PROP_N_TOTAL:
    {
      GVariant *var = g_value_get_variant (value);

      if (var != NULL)
      {
        GArray *n_total = g_array_new (FALSE, TRUE, sizeof (guint));

        ncm_cfg_array_set_variant (n_total, var);

        if (grid->constructed && (n_total->len != grid->len))
          g_error ("nc_galaxy_wl_obs_cheb_grid: %u node counts for %u galaxies.", n_total->len, grid->len);

        g_clear_pointer (&grid->n_total, g_array_unref);
        grid->n_total = n_total;
      }

      break;
    }
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_galaxy_wl_obs_cheb_grid_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcGalaxyWLObsChebGrid *grid = NC_GALAXY_WL_OBS_CHEB_GRID (object);

  switch (prop_id)
  {
    case PROP_LEN:
      g_value_set_uint (value, grid->len);
      break;
    case PROP_N_WEIGHTS:
      g_value_set_uint (value, grid->n_weights);
      break;
    case PROP_N_NODES:
      g_value_set_uint (value, grid->n_nodes);
      break;
    case PROP_Z_CL:
      g_value_set_double (value, grid->z_cl);
      break;
    case PROP_Z_TYPE:
      g_value_set_string (value, grid->z_type);
      break;
    case PROP_DATA_HASH:
      g_value_set_string (value, grid->data_hash);
      break;
    case PROP_WEIGHTS:
      g_value_set_object (value, grid->weights);
      break;
    case PROP_Z_NODES:
      g_value_set_object (value, grid->z_nodes);
      break;
    case PROP_NORMS:
      g_value_set_object (value, grid->norms);
      break;
    case PROP_N_TOTAL:
      g_value_take_variant (value, ncm_cfg_array_to_variant (grid->n_total, G_VARIANT_TYPE_UINT32));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

/*
 * Allocates a zero matrix of @ncols columns into an empty @slot, or checks the size of
 * the one set before construction.
 */
static void
_nc_galaxy_wl_obs_cheb_grid_alloc_matrix (NcGalaxyWLObsChebGrid *grid, NcmMatrix **slot, const guint ncols, const gchar *name)
{
  if (*slot == NULL)
  {
    *slot = ncm_matrix_new (grid->len, ncols);
    ncm_matrix_set_zero (*slot);
  }
  else if ((ncm_matrix_nrows (*slot) != grid->len) || (ncm_matrix_ncols (*slot) != ncols))
  {
    g_error ("nc_galaxy_wl_obs_cheb_grid: the %s matrix is %u x %u, not %u x %u.",
             name, ncm_matrix_nrows (*slot), ncm_matrix_ncols (*slot), grid->len, ncols);
  }
}

static void
_nc_galaxy_wl_obs_cheb_grid_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (nc_galaxy_wl_obs_cheb_grid_parent_class)->constructed (object);
  {
    NcGalaxyWLObsChebGrid *grid = NC_GALAXY_WL_OBS_CHEB_GRID (object);

    if (grid->n_nodes == 0)
      g_error ("nc_galaxy_wl_obs_cheb_grid: a grid needs at least one redshift value per galaxy.");

    _nc_galaxy_wl_obs_cheb_grid_alloc_matrix (grid, &grid->weights, grid->n_weights, "weights");
    _nc_galaxy_wl_obs_cheb_grid_alloc_matrix (grid, &grid->z_nodes, grid->n_nodes, "z-nodes");
    _nc_galaxy_wl_obs_cheb_grid_alloc_matrix (grid, &grid->norms, 3, "norms");

    if (grid->n_total == NULL)
    {
      grid->n_total = g_array_new (FALSE, TRUE, sizeof (guint));
      g_array_set_size (grid->n_total, grid->len);
    }
    else if (grid->n_total->len != grid->len)
    {
      g_error ("nc_galaxy_wl_obs_cheb_grid: %u node counts for %u galaxies.", grid->n_total->len, grid->len);
    }

    grid->constructed = TRUE;
  }
}

static void
_nc_galaxy_wl_obs_cheb_grid_dispose (GObject *object)
{
  NcGalaxyWLObsChebGrid *grid = NC_GALAXY_WL_OBS_CHEB_GRID (object);

  ncm_matrix_clear (&grid->weights);
  ncm_matrix_clear (&grid->z_nodes);
  ncm_matrix_clear (&grid->norms);
  g_clear_pointer (&grid->n_total, g_array_unref);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_galaxy_wl_obs_cheb_grid_parent_class)->dispose (object);
}

static void
_nc_galaxy_wl_obs_cheb_grid_finalize (GObject *object)
{
  NcGalaxyWLObsChebGrid *grid = NC_GALAXY_WL_OBS_CHEB_GRID (object);

  g_clear_pointer (&grid->z_type, g_free);
  g_clear_pointer (&grid->data_hash, g_free);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_galaxy_wl_obs_cheb_grid_parent_class)->finalize (object);
}

static void
nc_galaxy_wl_obs_cheb_grid_class_init (NcGalaxyWLObsChebGridClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = &_nc_galaxy_wl_obs_cheb_grid_set_property;
  object_class->get_property = &_nc_galaxy_wl_obs_cheb_grid_get_property;
  object_class->constructed  = &_nc_galaxy_wl_obs_cheb_grid_constructed;
  object_class->dispose      = &_nc_galaxy_wl_obs_cheb_grid_dispose;
  object_class->finalize     = &_nc_galaxy_wl_obs_cheb_grid_finalize;

  /**
   * NcGalaxyWLObsChebGrid:len:
   *
   * The number of galaxies.
   */
  g_object_class_install_property (object_class,
                                   PROP_LEN,
                                   g_param_spec_uint ("len",
                                                      "Length",
                                                      "Number of galaxies",
                                                      0, G_MAXUINT, 0,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_STRINGS));

  /**
   * NcGalaxyWLObsChebGrid:n-weights:
   *
   * The number of weights of a galaxy, the weights of all the nested levels.
   */
  g_object_class_install_property (object_class,
                                   PROP_N_WEIGHTS,
                                   g_param_spec_uint ("n-weights",
                                                      "Number of weights",
                                                      "Number of weights per galaxy",
                                                      0, G_MAXUINT, 0,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_STRINGS));

  /**
   * NcGalaxyWLObsChebGrid:n-nodes:
   *
   * The largest number of redshift values of a galaxy: the lower end of the p(z) support
   * and the Chebyshev nodes.
   */
  g_object_class_install_property (object_class,
                                   PROP_N_NODES,
                                   g_param_spec_uint ("n-nodes",
                                                      "Number of nodes",
                                                      "Number of redshift values per galaxy",
                                                      0, G_MAXUINT, 0,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_STRINGS));

  /**
   * NcGalaxyWLObsChebGrid:z-cl:
   *
   * The cluster redshift the grid was built for.
   */
  g_object_class_install_property (object_class,
                                   PROP_Z_CL,
                                   g_param_spec_double ("z-cl",
                                                        "Cluster redshift",
                                                        "Cluster redshift of the grid",
                                                        0.0, G_MAXDOUBLE, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_STRINGS));

  /**
   * NcGalaxyWLObsChebGrid:z-type:
   *
   * The type name of the #NcGalaxyRedshiftFactor that built the grid.
   */
  g_object_class_install_property (object_class,
                                   PROP_Z_TYPE,
                                   g_param_spec_string ("z-type",
                                                        "Redshift factor type",
                                                        "Type name of the redshift factor of the grid",
                                                        NULL,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_STRINGS));

  /**
   * NcGalaxyWLObsChebGrid:data-hash:
   *
   * The data hash (ncm_catalog_get_data_hash()) of the rows the grid was built from.
   */
  g_object_class_install_property (object_class,
                                   PROP_DATA_HASH,
                                   g_param_spec_string ("data-hash",
                                                        "Data hash",
                                                        "Data hash of the rows of the grid",
                                                        NULL,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_STRINGS));

  /**
   * NcGalaxyWLObsChebGrid:weights:
   *
   * The weights, one row per galaxy.
   */
  g_object_class_install_property (object_class,
                                   PROP_WEIGHTS,
                                   g_param_spec_object ("weights",
                                                        "Weights",
                                                        "Weights of the nested levels",
                                                        NCM_TYPE_MATRIX,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_STRINGS));

  /**
   * NcGalaxyWLObsChebGrid:z-nodes:
   *
   * The redshift values, one row per galaxy.
   */
  g_object_class_install_property (object_class,
                                   PROP_Z_NODES,
                                   g_param_spec_object ("z-nodes",
                                                        "Redshift values",
                                                        "Lower end and Chebyshev nodes",
                                                        NCM_TYPE_MATRIX,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_STRINGS));

  /**
   * NcGalaxyWLObsChebGrid:norms:
   *
   * The normalizations of p(z), one row of total, background and foreground per galaxy.
   */
  g_object_class_install_property (object_class,
                                   PROP_NORMS,
                                   g_param_spec_object ("norms",
                                                        "Normalizations",
                                                        "Total, background and foreground normalizations",
                                                        NCM_TYPE_MATRIX,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_STRINGS));

  /**
   * NcGalaxyWLObsChebGrid:n-total:
   *
   * The number of Chebyshev nodes of each galaxy, 0 for a galaxy without background
   * support.
   */
  g_object_class_install_property (object_class,
                                   PROP_N_TOTAL,
                                   g_param_spec_variant ("n-total",
                                                         "Node counts",
                                                         "Number of Chebyshev nodes per galaxy",
                                                         G_VARIANT_TYPE ("au"), NULL,
                                                         G_PARAM_READWRITE | G_PARAM_STATIC_STRINGS));
}

static void
_nc_galaxy_wl_obs_cheb_grid_set_matrix (NcGalaxyWLObsChebGrid *grid, NcmMatrix **slot, NcmMatrix *m, const guint ncols, const gchar *name)
{
  if (m == NULL)
    return;

  if (grid->constructed && ((ncm_matrix_nrows (m) != grid->len) || (ncm_matrix_ncols (m) != ncols)))
    g_error ("nc_galaxy_wl_obs_cheb_grid: the %s matrix is %u x %u, not %u x %u.",
             name, ncm_matrix_nrows (m), ncm_matrix_ncols (m), grid->len, ncols);

  ncm_matrix_ref (m);
  ncm_matrix_clear (slot);
  *slot = m;
}

/**
 * nc_galaxy_wl_obs_cheb_grid_new:
 * @len: the number of galaxies
 * @n_weights: the number of weights per galaxy
 * @n_nodes: the largest number of redshift values per galaxy
 * @z_cl: the cluster redshift
 * @z_type: the type name of the redshift factor
 * @data_hash: the data hash of the rows
 *
 * Creates a grid of @len galaxies with every value zero.
 *
 * Returns: (transfer full): a new #NcGalaxyWLObsChebGrid.
 */
NcGalaxyWLObsChebGrid *
nc_galaxy_wl_obs_cheb_grid_new (const guint len, const guint n_weights, const guint n_nodes, const gdouble z_cl, const gchar *z_type, const gchar *data_hash)
{
  return g_object_new (NC_TYPE_GALAXY_WL_OBS_CHEB_GRID,
                       "len", len,
                       "n-weights", n_weights,
                       "n-nodes", n_nodes,
                       "z-cl", z_cl,
                       "z-type", z_type,
                       "data-hash", data_hash,
                       NULL);
}

/**
 * nc_galaxy_wl_obs_cheb_grid_ref:
 * @grid: a #NcGalaxyWLObsChebGrid
 *
 * Increases the reference count of @grid by one.
 *
 * Returns: (transfer full): @grid.
 */
NcGalaxyWLObsChebGrid *
nc_galaxy_wl_obs_cheb_grid_ref (NcGalaxyWLObsChebGrid *grid)
{
  return g_object_ref (grid);
}

/**
 * nc_galaxy_wl_obs_cheb_grid_free:
 * @grid: a #NcGalaxyWLObsChebGrid
 *
 * Decreases the reference count of @grid by one.
 */
void
nc_galaxy_wl_obs_cheb_grid_free (NcGalaxyWLObsChebGrid *grid)
{
  g_object_unref (grid);
}

/**
 * nc_galaxy_wl_obs_cheb_grid_clear:
 * @grid: a #NcGalaxyWLObsChebGrid
 *
 * Decreases the reference count of *@grid by one and sets *@grid to %NULL, when it is
 * not %NULL.
 */
void
nc_galaxy_wl_obs_cheb_grid_clear (NcGalaxyWLObsChebGrid **grid)
{
  g_clear_object (grid);
}

/**
 * nc_galaxy_wl_obs_cheb_grid_len:
 * @grid: a #NcGalaxyWLObsChebGrid
 *
 * Returns: the number of galaxies.
 */
guint
nc_galaxy_wl_obs_cheb_grid_len (NcGalaxyWLObsChebGrid *grid)
{
  return grid->len;
}

/**
 * nc_galaxy_wl_obs_cheb_grid_get_n_weights:
 * @grid: a #NcGalaxyWLObsChebGrid
 *
 * Returns: the number of weights per galaxy.
 */
guint
nc_galaxy_wl_obs_cheb_grid_get_n_weights (NcGalaxyWLObsChebGrid *grid)
{
  return grid->n_weights;
}

/**
 * nc_galaxy_wl_obs_cheb_grid_get_n_nodes:
 * @grid: a #NcGalaxyWLObsChebGrid
 *
 * Returns: the largest number of redshift values per galaxy.
 */
guint
nc_galaxy_wl_obs_cheb_grid_get_n_nodes (NcGalaxyWLObsChebGrid *grid)
{
  return grid->n_nodes;
}

/**
 * nc_galaxy_wl_obs_cheb_grid_get_z_cl:
 * @grid: a #NcGalaxyWLObsChebGrid
 *
 * Returns: the cluster redshift the grid was built for.
 */
gdouble
nc_galaxy_wl_obs_cheb_grid_get_z_cl (NcGalaxyWLObsChebGrid *grid)
{
  return grid->z_cl;
}

/**
 * nc_galaxy_wl_obs_cheb_grid_peek_z_type:
 * @grid: a #NcGalaxyWLObsChebGrid
 *
 * Returns: (transfer none): the type name of the redshift factor that built the grid.
 */
const gchar *
nc_galaxy_wl_obs_cheb_grid_peek_z_type (NcGalaxyWLObsChebGrid *grid)
{
  return grid->z_type;
}

/**
 * nc_galaxy_wl_obs_cheb_grid_peek_data_hash:
 * @grid: a #NcGalaxyWLObsChebGrid
 *
 * Returns: (transfer none): the data hash of the rows the grid was built from.
 */
const gchar *
nc_galaxy_wl_obs_cheb_grid_peek_data_hash (NcGalaxyWLObsChebGrid *grid)
{
  return grid->data_hash;
}

/**
 * nc_galaxy_wl_obs_cheb_grid_set_galaxy:
 * @grid: a #NcGalaxyWLObsChebGrid
 * @i: the galaxy index
 * @n_total: the number of Chebyshev nodes of the galaxy
 * @weights: (nullable) (element-type gdouble): the weights, %NULL when @n_total is 0
 * @z_nodes: the lower end of the p(z) support followed by the @n_total nodes
 * @norm: the normalization of p(z)
 * @bg_norm: the background normalization
 * @fg_norm: the foreground normalization
 *
 * Sets the grid of galaxy @i. A galaxy index out of range, a length of @z_nodes other
 * than @n_total + 1 or above #NcGalaxyWLObsChebGrid:n-nodes, or a length of @weights
 * other than #NcGalaxyWLObsChebGrid:n-weights when @n_total is positive aborts.
 */
void
nc_galaxy_wl_obs_cheb_grid_set_galaxy (NcGalaxyWLObsChebGrid *grid, const guint i, const guint n_total, GArray *weights, NcmVector *z_nodes, const gdouble norm, const gdouble bg_norm, const gdouble fg_norm)
{
  guint k;

  if (i >= grid->len)
    g_error ("nc_galaxy_wl_obs_cheb_grid_set_galaxy: galaxy %u of %u.", i, grid->len);

  if ((ncm_vector_len (z_nodes) != n_total + 1) || (ncm_vector_len (z_nodes) > grid->n_nodes))
    g_error ("nc_galaxy_wl_obs_cheb_grid_set_galaxy: %u redshift values for %u nodes and at most %u values.",
             ncm_vector_len (z_nodes), n_total, grid->n_nodes);

  if ((n_total > 0) && ((weights == NULL) || (weights->len != grid->n_weights)))
    g_error ("nc_galaxy_wl_obs_cheb_grid_set_galaxy: %u weights, not %u.",
             (weights == NULL) ? 0 : weights->len, grid->n_weights);

  for (k = 0; k < grid->n_weights; k++)
    ncm_matrix_set (grid->weights, i, k, (n_total > 0) ? g_array_index (weights, gdouble, k) : 0.0);

  for (k = 0; k < grid->n_nodes; k++)
    ncm_matrix_set (grid->z_nodes, i, k, (k <= n_total) ? ncm_vector_get (z_nodes, k) : 0.0);

  ncm_matrix_set (grid->norms, i, 0, norm);
  ncm_matrix_set (grid->norms, i, 1, bg_norm);
  ncm_matrix_set (grid->norms, i, 2, fg_norm);

  g_array_index (grid->n_total, guint, i) = n_total;
}

/**
 * nc_galaxy_wl_obs_cheb_grid_get_galaxy:
 * @grid: a #NcGalaxyWLObsChebGrid
 * @i: the galaxy index
 * @weights: (out) (transfer full) (nullable) (element-type gdouble): the weights, %NULL when the galaxy has no node
 * @z_nodes: (out) (transfer full): the lower end of the p(z) support followed by the nodes
 * @norm: (out): the normalization of p(z)
 * @bg_norm: (out): the background normalization
 * @fg_norm: (out): the foreground normalization
 *
 * Gets the grid of galaxy @i in new arrays. A galaxy index out of range aborts.
 *
 * Returns: the number of Chebyshev nodes of the galaxy.
 */
guint
nc_galaxy_wl_obs_cheb_grid_get_galaxy (NcGalaxyWLObsChebGrid *grid, const guint i, GArray **weights, NcmVector **z_nodes, gdouble *norm, gdouble *bg_norm, gdouble *fg_norm)
{
  guint n_total, k;

  if (i >= grid->len)
    g_error ("nc_galaxy_wl_obs_cheb_grid_get_galaxy: galaxy %u of %u.", i, grid->len);

  n_total = g_array_index (grid->n_total, guint, i);

  if (n_total > 0)
  {
    *weights = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), grid->n_weights);
    g_array_append_vals (*weights, ncm_matrix_ptr (grid->weights, i, 0), grid->n_weights);
  }
  else
  {
    *weights = NULL;
  }

  *z_nodes = ncm_vector_new (n_total + 1);

  for (k = 0; k <= n_total; k++)
    ncm_vector_set (*z_nodes, k, ncm_matrix_get (grid->z_nodes, i, k));

  *norm    = ncm_matrix_get (grid->norms, i, 0);
  *bg_norm = ncm_matrix_get (grid->norms, i, 1);
  *fg_norm = ncm_matrix_get (grid->norms, i, 2);

  return n_total;
}


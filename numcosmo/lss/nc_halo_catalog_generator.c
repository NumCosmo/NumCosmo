/***************************************************************************
 *            nc_halo_catalog_generator.c
 *
 *  Sat Jun 13 17:30 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_halo_catalog_generator.c
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
 * NcHaloCatalogGenerator:
 *
 * Samples a #NcHaloCatalog from a cluster abundance model.
 *
 * A reusable object that draws a Poisson realization of the cluster number
 * counts from a #NcClusterAbundance and, for each drawn object, samples the
 * observed redshift and mass through the #NcClusterRedshift and #NcClusterMass
 * models found in the #NcmMSet. The result is returned as a #NcHaloCatalog of
 * kind %NC_HALO_CATALOG_KIND_CLUSTER.
 *
 * It assembles the existing sampling machinery rather than reimplementing it,
 * and is meant to be reused by the likelihood data objects (for example
 * #NcDataClusterNCount) and by higher-level mock generators. Richer mocks
 * (sky positions, member galaxies) are expected to be built by composition,
 * holding extra collaborators that augment the generated catalog.
 *
 * The generated catalog columns are, in order:
 *
 * - `z_true`, `lnM_true`: the true redshift and mass;
 * - `z_obs_0` .. `z_obs_{n-1}`: the observed redshift block (`n` = redshift
 *   observable length);
 * - `lnM_obs_0` .. `lnM_obs_{m-1}`: the observed mass block;
 * - `z_obs_param_0` .. and `lnM_obs_param_0` .. : the observable parameter
 *   blocks, when the models declare any;
 * - `ra`, `dec`: sky position (degrees), only when a #NcmSkyFootprint is set
 *   (see nc_halo_catalog_generator_set_footprint());
 * - `r_Delta`: the spherical-overdensity radius (Mpc) of the true halo, only when
 *   radius output is enabled (see nc_halo_catalog_generator_set_with_radius()).
 *   It is computed from $M_\Delta = \frac{4\pi}{3}\,\Delta\,\rho_\mathrm{bg}\,
 *   r_\Delta^3$ using the mass definition ($\Delta$ and background density)
 *   carried by the multiplicity function used to draw the counts.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"
#include "lss/nc_halo_catalog_generator.h"
#include "lss/nc_cluster_abundance.h"
#include "lss/nc_cluster_redshift.h"
#include "lss/nc_cluster_mass.h"
#include "lss/nc_halo_mass_function.h"
#include "lss/nc_multiplicity_func.h"
#include "nc_hicosmo.h"
#include "math/ncm_c.h"
#include "math/ncm_matrix.h"
#include "math/ncm_spline.h"
#include "math/ncm_spline2d.h"

typedef struct _NcHaloCatalogGeneratorPrivate
{
  NcClusterAbundance *cad;
  NcmSkyFootprint *footprint;
  gboolean with_radius;
} NcHaloCatalogGeneratorPrivate;

enum
{
  PROP_0,
  PROP_ABUNDANCE,
  PROP_FOOTPRINT,
  PROP_WITH_RADIUS,
};

struct _NcHaloCatalogGenerator
{
  GObject parent_instance;
};

G_DEFINE_TYPE_WITH_PRIVATE (NcHaloCatalogGenerator, nc_halo_catalog_generator, G_TYPE_OBJECT)

static void
nc_halo_catalog_generator_init (NcHaloCatalogGenerator *gen)
{
  NcHaloCatalogGeneratorPrivate * const self = nc_halo_catalog_generator_get_instance_private (gen);

  self->cad         = NULL;
  self->footprint   = NULL;
  self->with_radius = FALSE;
}

static void
_nc_halo_catalog_generator_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcHaloCatalogGenerator *gen                = NC_HALO_CATALOG_GENERATOR (object);
  NcHaloCatalogGeneratorPrivate * const self = nc_halo_catalog_generator_get_instance_private (gen);

  switch (prop_id)
  {
    case PROP_ABUNDANCE:
      nc_cluster_abundance_clear (&self->cad);
      self->cad = g_value_dup_object (value);
      break;
    case PROP_FOOTPRINT:
      nc_halo_catalog_generator_set_footprint (gen, g_value_get_object (value));
      break;
    case PROP_WITH_RADIUS:
      nc_halo_catalog_generator_set_with_radius (gen, g_value_get_boolean (value));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_halo_catalog_generator_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcHaloCatalogGenerator *gen                = NC_HALO_CATALOG_GENERATOR (object);
  NcHaloCatalogGeneratorPrivate * const self = nc_halo_catalog_generator_get_instance_private (gen);

  switch (prop_id)
  {
    case PROP_ABUNDANCE:
      g_value_set_object (value, self->cad);
      break;
    case PROP_FOOTPRINT:
      g_value_set_object (value, self->footprint);
      break;
    case PROP_WITH_RADIUS:
      g_value_set_boolean (value, self->with_radius);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_halo_catalog_generator_dispose (GObject *object)
{
  NcHaloCatalogGenerator *gen                = NC_HALO_CATALOG_GENERATOR (object);
  NcHaloCatalogGeneratorPrivate * const self = nc_halo_catalog_generator_get_instance_private (gen);

  nc_cluster_abundance_clear (&self->cad);
  ncm_sky_footprint_clear (&self->footprint);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_halo_catalog_generator_parent_class)->dispose (object);
}

static void
nc_halo_catalog_generator_class_init (NcHaloCatalogGeneratorClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = &_nc_halo_catalog_generator_set_property;
  object_class->get_property = &_nc_halo_catalog_generator_get_property;
  object_class->dispose      = &_nc_halo_catalog_generator_dispose;

  /**
   * NcHaloCatalogGenerator:abundance:
   *
   * The cluster abundance model used to draw the number counts.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_ABUNDANCE,
                                   g_param_spec_object ("abundance",
                                                        "Abundance",
                                                        "Cluster abundance model",
                                                        NC_TYPE_CLUSTER_ABUNDANCE,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_STRINGS));

  /**
   * NcHaloCatalogGenerator:footprint:
   *
   * Optional sky footprint. When set, the generated catalog gains `ra` and `dec`
   * columns sampled over the footprint; when unset, no positions are generated.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_FOOTPRINT,
                                   g_param_spec_object ("footprint",
                                                        "Footprint",
                                                        "Optional sky footprint for position sampling",
                                                        NCM_TYPE_SKY_FOOTPRINT,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_STRINGS));

  /**
   * NcHaloCatalogGenerator:with-radius:
   *
   * Whether to append an `r_Delta` column with the spherical-overdensity radius
   * (Mpc) of each true halo, computed from the mass definition carried by the
   * multiplicity function. Disabled by default.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_WITH_RADIUS,
                                   g_param_spec_boolean ("with-radius",
                                                         "With radius",
                                                         "Whether to output the halo spherical-overdensity radius",
                                                         FALSE,
                                                         G_PARAM_READWRITE | G_PARAM_STATIC_STRINGS));
}

/**
 * nc_halo_catalog_generator_new:
 * @cad: a #NcClusterAbundance
 *
 * Creates a new #NcHaloCatalogGenerator using @cad to draw the number counts.
 *
 * Returns: (transfer full): a new #NcHaloCatalogGenerator.
 */
NcHaloCatalogGenerator *
nc_halo_catalog_generator_new (NcClusterAbundance *cad)
{
  return g_object_new (NC_TYPE_HALO_CATALOG_GENERATOR,
                       "abundance", cad,
                       NULL);
}

/**
 * nc_halo_catalog_generator_ref:
 * @gen: a #NcHaloCatalogGenerator
 *
 * Increases the reference count of @gen by one.
 *
 * Returns: (transfer full): @gen.
 */
NcHaloCatalogGenerator *
nc_halo_catalog_generator_ref (NcHaloCatalogGenerator *gen)
{
  return g_object_ref (gen);
}

/**
 * nc_halo_catalog_generator_free:
 * @gen: a #NcHaloCatalogGenerator
 *
 * Decreases the reference count of @gen by one.
 *
 */
void
nc_halo_catalog_generator_free (NcHaloCatalogGenerator *gen)
{
  g_object_unref (gen);
}

/**
 * nc_halo_catalog_generator_clear:
 * @gen: a #NcHaloCatalogGenerator
 *
 * If *@gen is different from %NULL, decreases its reference count and sets
 * *@gen to %NULL.
 *
 */
void
nc_halo_catalog_generator_clear (NcHaloCatalogGenerator **gen)
{
  g_clear_object (gen);
}

/**
 * nc_halo_catalog_generator_peek_abundance:
 * @gen: a #NcHaloCatalogGenerator
 *
 * Gets the cluster abundance model used by @gen.
 *
 * Returns: (transfer none): the #NcClusterAbundance.
 */
NcClusterAbundance *
nc_halo_catalog_generator_peek_abundance (NcHaloCatalogGenerator *gen)
{
  NcHaloCatalogGeneratorPrivate * const self = nc_halo_catalog_generator_get_instance_private (gen);

  return self->cad;
}

/**
 * nc_halo_catalog_generator_set_footprint:
 * @gen: a #NcHaloCatalogGenerator
 * @footprint: (nullable): a #NcmSkyFootprint, or %NULL to disable position sampling
 *
 * Sets the sky footprint used to sample positions. When a footprint is set, the
 * generated catalog gains `ra` and `dec` columns; when %NULL, no positions are
 * generated.
 *
 */
void
nc_halo_catalog_generator_set_footprint (NcHaloCatalogGenerator *gen, NcmSkyFootprint *footprint)
{
  NcHaloCatalogGeneratorPrivate * const self = nc_halo_catalog_generator_get_instance_private (gen);

  ncm_sky_footprint_clear (&self->footprint);

  if (footprint != NULL)
    self->footprint = ncm_sky_footprint_ref (footprint);
}

/**
 * nc_halo_catalog_generator_peek_footprint:
 * @gen: a #NcHaloCatalogGenerator
 *
 * Gets the sky footprint used to sample positions, or %NULL when none is set.
 *
 * Returns: (transfer none) (nullable): the #NcmSkyFootprint.
 */
NcmSkyFootprint *
nc_halo_catalog_generator_peek_footprint (NcHaloCatalogGenerator *gen)
{
  NcHaloCatalogGeneratorPrivate * const self = nc_halo_catalog_generator_get_instance_private (gen);

  return self->footprint;
}

/**
 * nc_halo_catalog_generator_set_with_radius:
 * @gen: a #NcHaloCatalogGenerator
 * @with_radius: whether to output the halo radius
 *
 * Sets whether the generated catalog gains an `r_Delta` column with the
 * spherical-overdensity radius (Mpc) of each true halo.
 *
 */
void
nc_halo_catalog_generator_set_with_radius (NcHaloCatalogGenerator *gen, gboolean with_radius)
{
  NcHaloCatalogGeneratorPrivate * const self = nc_halo_catalog_generator_get_instance_private (gen);

  self->with_radius = with_radius;
}

/**
 * nc_halo_catalog_generator_get_with_radius:
 * @gen: a #NcHaloCatalogGenerator
 *
 * Gets whether the generated catalog includes the `r_Delta` radius column.
 *
 * Returns: %TRUE when the radius column is generated.
 */
gboolean
nc_halo_catalog_generator_get_with_radius (NcHaloCatalogGenerator *gen)
{
  NcHaloCatalogGeneratorPrivate * const self = nc_halo_catalog_generator_get_instance_private (gen);

  return self->with_radius;
}

static gdouble
_nc_halo_catalog_generator_radius (NcMultiplicityFunc *mulf, NcHICosmo *cosmo, gdouble lnM, gdouble z)
{
  const gdouble rho_crit0  = ncm_c_crit_mass_density_h2_solar_mass_Mpc3 () * nc_hicosmo_h2 (cosmo);
  const gdouble Delta      = nc_multiplicity_func_get_Delta (mulf);
  const gdouble M          = exp (lnM);
  gdouble Delta_rho_bg;

  switch (nc_multiplicity_func_get_mdef (mulf))
  {
    case NC_MULTIPLICITY_FUNC_MASS_DEF_MEAN:
      Delta_rho_bg = Delta * rho_crit0 * nc_hicosmo_E2Omega_m (cosmo, z);
      break;
    case NC_MULTIPLICITY_FUNC_MASS_DEF_CRITICAL:
      Delta_rho_bg = Delta * rho_crit0 * nc_hicosmo_E2 (cosmo, z);
      break;
    default:
      g_error ("nc_halo_catalog_generator: radius output unsupported for this mass definition");
      return 0.0;
  }

  return cbrt (3.0 * M / (4.0 * M_PI * Delta_rho_bg));
}

static GStrv
_nc_halo_catalog_generator_build_columns (guint z_obs_len, guint z_obs_params_len, guint lnM_obs_len, guint lnM_obs_params_len, gboolean with_position, gboolean with_radius)
{
  GPtrArray *cols = g_ptr_array_new ();
  guint k;

  g_ptr_array_add (cols, g_strdup ("z_true"));
  g_ptr_array_add (cols, g_strdup ("lnM_true"));

  for (k = 0; k < z_obs_len; k++)
    g_ptr_array_add (cols, g_strdup_printf ("z_obs_%u", k));

  for (k = 0; k < lnM_obs_len; k++)
    g_ptr_array_add (cols, g_strdup_printf ("lnM_obs_%u", k));

  for (k = 0; k < z_obs_params_len; k++)
    g_ptr_array_add (cols, g_strdup_printf ("z_obs_param_%u", k));

  for (k = 0; k < lnM_obs_params_len; k++)
    g_ptr_array_add (cols, g_strdup_printf ("lnM_obs_param_%u", k));

  if (with_position)
  {
    g_ptr_array_add (cols, g_strdup ("ra"));
    g_ptr_array_add (cols, g_strdup ("dec"));
  }

  if (with_radius)
    g_ptr_array_add (cols, g_strdup ("r_Delta"));

  g_ptr_array_add (cols, NULL);

  return (GStrv) g_ptr_array_free (cols, FALSE);
}

/**
 * nc_halo_catalog_generator_generate:
 * @gen: a #NcHaloCatalogGenerator
 * @mset: a #NcmMSet holding the cosmology and the cluster redshift/mass models
 * @rng: a #NcmRNG
 *
 * Draws a Poisson realization of the cluster number counts and samples the
 * observed redshift and mass of each object, returning them as a
 * #NcHaloCatalog of kind %NC_HALO_CATALOG_KIND_CLUSTER. See the class
 * description for the column layout.
 *
 * Returns: (transfer full): a new #NcHaloCatalog with the sampled objects.
 */
NcHaloCatalog *
nc_halo_catalog_generator_generate (NcHaloCatalogGenerator *gen, NcmMSet *mset, NcmRNG *rng)
{
  NcHaloCatalogGeneratorPrivate * const self = nc_halo_catalog_generator_get_instance_private (gen);
  NcClusterAbundance *cad                    = self->cad;
  NcHICosmo *cosmo                           = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));
  NcClusterRedshift *clusterz                = NC_CLUSTER_REDSHIFT (ncm_mset_peek (mset, nc_cluster_redshift_id ()));
  NcClusterMass *clusterm                    = NC_CLUSTER_MASS (ncm_mset_peek (mset, nc_cluster_mass_id ()));

  const guint z_obs_len          = nc_cluster_redshift_class_obs_len (NC_CLUSTER_REDSHIFT_GET_CLASS (clusterz));
  const guint z_obs_params_len   = nc_cluster_redshift_class_obs_params_len (NC_CLUSTER_REDSHIFT_GET_CLASS (clusterz));
  const guint lnM_obs_len        = nc_cluster_mass_class_obs_len (NC_CLUSTER_MASS_GET_CLASS (clusterm));
  const guint lnM_obs_params_len = nc_cluster_mass_class_obs_params_len (NC_CLUSTER_MASS_GET_CLASS (clusterm));
  const gboolean with_position   = self->footprint != NULL;
  const gboolean with_radius     = self->with_radius;
  NcMultiplicityFunc *mulf       = with_radius ? nc_halo_mass_function_peek_multiplicity_function (cad->mfp) : NULL;

  gdouble *zi_obs          = g_new (gdouble, z_obs_len);
  gdouble *zi_obs_params   = z_obs_params_len > 0 ? g_new (gdouble, z_obs_params_len) : NULL;
  gdouble *lnMi_obs        = g_new (gdouble, lnM_obs_len);
  gdouble *lnMi_obs_params = lnM_obs_params_len > 0 ? g_new (gdouble, lnM_obs_params_len) : NULL;

  GArray *rows = g_array_new (FALSE, FALSE, sizeof (gdouble));
  GStrv col_names;
  NcHaloCatalog *hcat;
  NcmMatrix *data;
  guint total_np;
  guint np = 0;
  guint i;

  ncm_rng_lock (rng);
  total_np = ncm_rng_poisson_gen (rng, cad->norma);
  ncm_rng_unlock (rng);

  nc_cluster_abundance_prepare_inv_dNdz (cad, cosmo, cad->lnMi);

  for (i = 0; i < total_np; i++)
  {
    ncm_rng_lock (rng);
    {
      const gdouble u1       = _nc_cad_inv_dNdz_convergence_f (ncm_rng_uniform01_pos_gen (rng), cad->z_epsilon);
      const gdouble u2       = _nc_cad_inv_dNdz_convergence_f (ncm_rng_uniform01_pos_gen (rng), cad->lnM_epsilon);
      const gdouble z_true   = ncm_spline_eval (cad->inv_z, u1);
      const gdouble lnM_true = ncm_spline2d_eval (cad->inv_lnM_z, u2, z_true);

      ncm_rng_unlock (rng);

      if (nc_cluster_redshift_resample (clusterz, cosmo, lnM_true, z_true, zi_obs, zi_obs_params, rng) &&
          nc_cluster_mass_resample (clusterm, cosmo, lnM_true, z_true, lnMi_obs, lnMi_obs_params, rng))
      {
        guint k;

        g_array_append_val (rows, z_true);
        g_array_append_val (rows, lnM_true);

        for (k = 0; k < z_obs_len; k++)
          g_array_append_val (rows, zi_obs[k]);

        for (k = 0; k < lnM_obs_len; k++)
          g_array_append_val (rows, lnMi_obs[k]);

        for (k = 0; k < z_obs_params_len; k++)
          g_array_append_val (rows, zi_obs_params[k]);

        for (k = 0; k < lnM_obs_params_len; k++)
          g_array_append_val (rows, lnMi_obs_params[k]);

        if (with_position)
        {
          gdouble ra, dec;

          ncm_sky_footprint_gen_ra_dec (self->footprint, rng, &ra, &dec);
          g_array_append_val (rows, ra);
          g_array_append_val (rows, dec);
        }

        if (with_radius)
        {
          const gdouble r_Delta = _nc_halo_catalog_generator_radius (mulf, cosmo, lnM_true, z_true);

          g_array_append_val (rows, r_Delta);
        }

        np++;
      }
    }
  }

  col_names = _nc_halo_catalog_generator_build_columns (z_obs_len, z_obs_params_len, lnM_obs_len, lnM_obs_params_len, with_position, with_radius);
  hcat      = nc_halo_catalog_new (NC_HALO_CATALOG_KIND_CLUSTER, NULL, NULL, np, col_names, NULL, 0);

  if (np > 0)
  {
    data = ncm_catalog_peek_data (NCM_CATALOG (hcat));
    ncm_matrix_set_from_array (data, rows);
  }

  g_strfreev (col_names);
  g_array_unref (rows);
  g_free (zi_obs);
  g_free (lnMi_obs);
  g_free (zi_obs_params);
  g_free (lnMi_obs_params);

  return hcat;
}

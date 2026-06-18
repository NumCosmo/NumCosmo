/***************************************************************************
 *            nc_halo_catalog_member_generator.c
 *
 *  Sun Jun 14 14:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_halo_catalog_member_generator.c
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
 * NcHaloCatalogMemberGenerator:
 *
 * Populates a host #NcHaloCatalog with member galaxies through a #NcGalaxyHOD.
 *
 * A standalone, reusable object that consumes an existing host catalog (halos or
 * clusters) and, for each host, draws its central and satellite galaxies from a
 * halo occupation distribution. The members are placed in 3D within the host and
 * projected back to sky coordinates and redshift, and returned as a linked
 * #NcHaloCatalog of kind %NC_HALO_CATALOG_KIND_MEMBER.
 *
 * The host catalog must provide the columns produced by
 * #NcHaloCatalogGenerator with positions and radius enabled:
 *
 * - `ra`, `dec`: host sky position (degrees);
 * - `z_true`: host redshift;
 * - `lnM_true`: host natural-log mass (fed to the HOD);
 * - `r_Delta`: host spherical-overdensity radius (Mpc).
 *
 * For each host the central galaxy sits at the host position; satellites are
 * placed uniformly in volume within `r_Delta` ($R = r_\Delta\,U^{1/3}$),
 * isotropically in direction ($\cos\theta$ uniform, $\phi$ uniform). The
 * displacement is decomposed about the line of sight: the radial part
 * $R\cos\theta$ becomes a Hubble-flow redshift offset $\delta z = H(z)\,
 * R\cos\theta / c$, and the transverse part $R\sin\theta$ becomes an angular
 * separation $R\sin\theta / d_A(z)$ applied to the host position with the
 * spherical law of cosines. (Two inconsistencies of the reference Python
 * prototype are not reproduced: it used $c$ in m/s while $H$ is in km/s/Mpc,
 * and it used the full $R$ instead of $R\sin\theta$ for the transverse
 * separation; here $c$ is in km/s and the projection uses $R\sin\theta$.)
 *
 * The generated member catalog has columns `galaxy_id` (id), `parent_id`
 * (linking back to the host id, or the host row index when the host has no id
 * column), `is_central`, `ra`, `dec` and `z`.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"
#include "nc/lss/halo/nc_halo_catalog_member_generator.h"
#include "nc/lss/galaxy/nc_galaxy_hod.h"
#include "nc/background/nc_hicosmo.h"
#include "nc/background/nc_distance.h"
#include "math/ncm_c.h"
#include "math/ncm_matrix.h"

#include <math.h>

/* *INDENT-OFF* */
G_DEFINE_QUARK (nc-halo-catalog-member-generator-error, nc_halo_catalog_member_generator_error)
/* *INDENT-ON* */

typedef struct _NcHaloCatalogMemberGeneratorPrivate
{
  NcGalaxyHOD *hod;
  NcDistance *dist;
} NcHaloCatalogMemberGeneratorPrivate;

enum
{
  PROP_0,
  PROP_HOD,
  PROP_DISTANCE,
};

struct _NcHaloCatalogMemberGenerator
{
  GObject parent_instance;
};

G_DEFINE_TYPE_WITH_PRIVATE (NcHaloCatalogMemberGenerator, nc_halo_catalog_member_generator, G_TYPE_OBJECT)

static void
nc_halo_catalog_member_generator_init (NcHaloCatalogMemberGenerator *memgen)
{
  NcHaloCatalogMemberGeneratorPrivate * const self = nc_halo_catalog_member_generator_get_instance_private (memgen);

  self->hod  = NULL;
  self->dist = NULL;
}

static void
_nc_halo_catalog_member_generator_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcHaloCatalogMemberGenerator *memgen             = NC_HALO_CATALOG_MEMBER_GENERATOR (object);
  NcHaloCatalogMemberGeneratorPrivate * const self = nc_halo_catalog_member_generator_get_instance_private (memgen);

  switch (prop_id)
  {
    case PROP_HOD:
      nc_galaxy_hod_clear (&self->hod);
      self->hod = g_value_dup_object (value);
      break;
    case PROP_DISTANCE:
      nc_halo_catalog_member_generator_set_distance (memgen, g_value_get_object (value));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_halo_catalog_member_generator_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcHaloCatalogMemberGenerator *memgen             = NC_HALO_CATALOG_MEMBER_GENERATOR (object);
  NcHaloCatalogMemberGeneratorPrivate * const self = nc_halo_catalog_member_generator_get_instance_private (memgen);

  switch (prop_id)
  {
    case PROP_HOD:
      g_value_set_object (value, self->hod);
      break;
    case PROP_DISTANCE:
      g_value_set_object (value, self->dist);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_halo_catalog_member_generator_dispose (GObject *object)
{
  NcHaloCatalogMemberGenerator *memgen             = NC_HALO_CATALOG_MEMBER_GENERATOR (object);
  NcHaloCatalogMemberGeneratorPrivate * const self = nc_halo_catalog_member_generator_get_instance_private (memgen);

  nc_galaxy_hod_clear (&self->hod);
  nc_distance_clear (&self->dist);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_halo_catalog_member_generator_parent_class)->dispose (object);
}

static void
nc_halo_catalog_member_generator_class_init (NcHaloCatalogMemberGeneratorClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = &_nc_halo_catalog_member_generator_set_property;
  object_class->get_property = &_nc_halo_catalog_member_generator_get_property;
  object_class->dispose      = &_nc_halo_catalog_member_generator_dispose;

  /**
   * NcHaloCatalogMemberGenerator:hod:
   *
   * The halo occupation distribution used to draw the member galaxies.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_HOD,
                                   g_param_spec_object ("hod",
                                                        "HOD",
                                                        "Halo occupation distribution",
                                                        NC_TYPE_GALAXY_HOD,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_STRINGS));

  /**
   * NcHaloCatalogMemberGenerator:distance:
   *
   * Optional #NcDistance used to convert the transverse member displacement into
   * an angular separation. When unset, a temporary one covering the host
   * redshift range is built for each generation.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_DISTANCE,
                                   g_param_spec_object ("distance",
                                                        "Distance",
                                                        "Distance object for angular-diameter distances",
                                                        NC_TYPE_DISTANCE,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_STRINGS));
}

/**
 * nc_halo_catalog_member_generator_new:
 * @hod: a #NcGalaxyHOD
 *
 * Creates a new #NcHaloCatalogMemberGenerator drawing members through @hod.
 *
 * Returns: (transfer full): a new #NcHaloCatalogMemberGenerator.
 */
NcHaloCatalogMemberGenerator *
nc_halo_catalog_member_generator_new (NcGalaxyHOD *hod)
{
  return g_object_new (NC_TYPE_HALO_CATALOG_MEMBER_GENERATOR,
                       "hod", hod,
                       NULL);
}

/**
 * nc_halo_catalog_member_generator_ref:
 * @memgen: a #NcHaloCatalogMemberGenerator
 *
 * Increases the reference count of @memgen by one.
 *
 * Returns: (transfer full): @memgen.
 */
NcHaloCatalogMemberGenerator *
nc_halo_catalog_member_generator_ref (NcHaloCatalogMemberGenerator *memgen)
{
  return g_object_ref (memgen);
}

/**
 * nc_halo_catalog_member_generator_free:
 * @memgen: a #NcHaloCatalogMemberGenerator
 *
 * Decreases the reference count of @memgen by one.
 *
 */
void
nc_halo_catalog_member_generator_free (NcHaloCatalogMemberGenerator *memgen)
{
  g_object_unref (memgen);
}

/**
 * nc_halo_catalog_member_generator_clear:
 * @memgen: a #NcHaloCatalogMemberGenerator
 *
 * If *@memgen is different from %NULL, decreases its reference count and sets
 * *@memgen to %NULL.
 *
 */
void
nc_halo_catalog_member_generator_clear (NcHaloCatalogMemberGenerator **memgen)
{
  g_clear_object (memgen);
}

/**
 * nc_halo_catalog_member_generator_peek_hod:
 * @memgen: a #NcHaloCatalogMemberGenerator
 *
 * Gets the halo occupation distribution used by @memgen.
 *
 * Returns: (transfer none): the #NcGalaxyHOD.
 */
NcGalaxyHOD *
nc_halo_catalog_member_generator_peek_hod (NcHaloCatalogMemberGenerator *memgen)
{
  NcHaloCatalogMemberGeneratorPrivate * const self = nc_halo_catalog_member_generator_get_instance_private (memgen);

  return self->hod;
}

/**
 * nc_halo_catalog_member_generator_set_distance:
 * @memgen: a #NcHaloCatalogMemberGenerator
 * @dist: (nullable): a #NcDistance, or %NULL to build one on demand
 *
 * Sets the distance object used to compute angular-diameter distances.
 *
 */
void
nc_halo_catalog_member_generator_set_distance (NcHaloCatalogMemberGenerator *memgen, NcDistance *dist)
{
  NcHaloCatalogMemberGeneratorPrivate * const self = nc_halo_catalog_member_generator_get_instance_private (memgen);

  nc_distance_clear (&self->dist);

  if (dist != NULL)
    self->dist = nc_distance_ref (dist);
}

/**
 * nc_halo_catalog_member_generator_peek_distance:
 * @memgen: a #NcHaloCatalogMemberGenerator
 *
 * Gets the distance object used by @memgen, or %NULL when none is set.
 *
 * Returns: (transfer none) (nullable): the #NcDistance.
 */
NcDistance *
nc_halo_catalog_member_generator_peek_distance (NcHaloCatalogMemberGenerator *memgen)
{
  NcHaloCatalogMemberGeneratorPrivate * const self = nc_halo_catalog_member_generator_get_instance_private (memgen);

  return self->dist;
}

static void
_nc_halo_catalog_member_generator_sky_offset (gdouble ra_c, gdouble dec_c, gdouble sep, gdouble phi, gdouble *ra_out, gdouble *dec_out)
{
  const gdouble ra_c_rad  = ra_c * (M_PI / 180.0);
  const gdouble dec_c_rad = dec_c * (M_PI / 180.0);
  const gdouble sin_dec   = sin (dec_c_rad) * cos (sep) + cos (dec_c_rad) * sin (sep) * cos (phi);
  const gdouble dec_rad   = asin (CLAMP (sin_dec, -1.0, 1.0));
  const gdouble y         = sin (phi) * sin (sep) * cos (dec_c_rad);
  const gdouble x         = cos (sep) - sin (dec_c_rad) * sin (dec_rad);
  const gdouble ra_rad    = ra_c_rad + atan2 (y, x);
  gdouble ra_deg          = fmod (ra_rad * (180.0 / M_PI), 360.0);

  if (ra_deg < 0.0)
    ra_deg += 360.0;

  *ra_out  = ra_deg;
  *dec_out = dec_rad * (180.0 / M_PI);
}

static gboolean
_nc_halo_catalog_member_generator_check_columns (NcmCatalog *host, GError **error)
{
  const gchar *const required[] = {"ra", "dec", "z_true", "lnM_true", "r_Delta"};
  guint i;

  for (i = 0; i < G_N_ELEMENTS (required); i++)
  {
    if (!ncm_catalog_has_column (host, required[i]))
    {
      g_set_error (error, NC_HALO_CATALOG_MEMBER_GENERATOR_ERROR, NC_HALO_CATALOG_MEMBER_GENERATOR_ERROR_MISSING_COLUMN,
                   "Host catalog is missing the required column '%s'.", required[i]);

      return FALSE;
    }
  }

  return TRUE;
}

/**
 * nc_halo_catalog_member_generator_generate:
 * @memgen: a #NcHaloCatalogMemberGenerator
 * @host: the host #NcHaloCatalog to populate
 * @mset: a #NcmMSet holding the cosmology
 * @rng: a #NcmRNG
 * @error: a #GError for error reporting
 *
 * Draws member galaxies for every object in @host through the HOD and places
 * them in 3D, returning them as a #NcHaloCatalog of kind
 * %NC_HALO_CATALOG_KIND_MEMBER. See the class description for the host schema and
 * the column layout of the result.
 *
 * Returns: (transfer full) (nullable): a new member #NcHaloCatalog, or %NULL on
 *     error.
 */
NcHaloCatalog *
nc_halo_catalog_member_generator_generate (NcHaloCatalogMemberGenerator *memgen, NcHaloCatalog *host, NcmMSet *mset, NcmRNG *rng, GError **error)
{
  NcHaloCatalogMemberGeneratorPrivate * const self = nc_halo_catalog_member_generator_get_instance_private (memgen);
  NcmCatalog *hostc                                = NCM_CATALOG (host);
  NcHICosmo *cosmo                                 = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));
  const gchar *id_col                              = nc_halo_catalog_peek_id_col (host);
  const gdouble c_kms                              = ncm_c_c () / 1.0e3;
  const guint nhost                                = ncm_catalog_len (hostc);

  GArray *rows = g_array_new (FALSE, FALSE, sizeof (gdouble));
  NcDistance *dist;
  gboolean local_dist = FALSE;
  gdouble zmax        = 0.0;
  guint galaxy_id     = 0;
  NcHaloCatalog *members;
  GStrv col_names;
  NcmCatalogColType col_types[6];
  guint i;

  g_return_val_if_fail (error == NULL || *error == NULL, NULL);

  if (!_nc_halo_catalog_member_generator_check_columns (hostc, error))
  {
    g_array_unref (rows);

    return NULL;
  }

  for (i = 0; i < nhost; i++)
  {
    const gdouble z_c = ncm_catalog_get (hostc, "z_true", i, NULL);

    zmax = MAX (zmax, z_c);
  }

  dist = self->dist;

  if (dist == NULL)
  {
    dist       = nc_distance_new (zmax + 0.5);
    local_dist = TRUE;
  }

  nc_distance_prepare_if_needed (dist, cosmo);

  for (i = 0; i < nhost; i++)
  {
    const gdouble ra_c  = ncm_catalog_get (hostc, "ra", i, NULL);
    const gdouble dec_c = ncm_catalog_get (hostc, "dec", i, NULL);
    const gdouble z_c   = ncm_catalog_get (hostc, "z_true", i, NULL);
    const gdouble lnM   = ncm_catalog_get (hostc, "lnM_true", i, NULL);
    const gdouble r_D   = ncm_catalog_get (hostc, "r_Delta", i, NULL);
    const gint64 parent = id_col != NULL ? nc_halo_catalog_get_id (host, i, NULL) : (gint64) i;
    const gdouble Hz    = nc_hicosmo_H (cosmo, z_c);
    const gdouble dA    = nc_distance_angular_diameter (dist, cosmo, z_c) * nc_hicosmo_RH_Mpc (cosmo);

    gint n_cen, n_sat;
    guint total, j;

    nc_galaxy_hod_gen (self->hod, lnM, rng, &n_cen, &n_sat);
    total = (guint) (n_cen + n_sat);

    for (j = 0; j < total; j++)
    {
      const gboolean is_central = j < (guint) n_cen;
      gdouble R, cos_theta, phi;
      gdouble z_gal, sep, ra_gal, dec_gal;
      gdouble gid, pid, cen;

      gdouble sin_theta;

      ncm_rng_lock (rng);
      R         = is_central ? 0.0 : r_D * cbrt (ncm_rng_uniform01_gen (rng));
      cos_theta = ncm_rng_uniform_gen (rng, -1.0, 1.0);
      phi       = ncm_rng_uniform_gen (rng, 0.0, 2.0 * M_PI);
      ncm_rng_unlock (rng);

      sin_theta = sqrt (1.0 - cos_theta * cos_theta);

      /* Decompose the 3D displacement into a line-of-sight part R cos(theta),
       * turned into a Hubble-flow redshift offset, and a transverse part
       * R sin(theta), turned into an angular separation. */
      z_gal = z_c + (Hz / c_kms) * (R * cos_theta);
      sep   = R * sin_theta / dA;

      _nc_halo_catalog_member_generator_sky_offset (ra_c, dec_c, sep, phi, &ra_gal, &dec_gal);

      gid = (gdouble) galaxy_id;
      pid = (gdouble) parent;
      cen = is_central ? 1.0 : 0.0;

      g_array_append_val (rows, gid);
      g_array_append_val (rows, pid);
      g_array_append_val (rows, cen);
      g_array_append_val (rows, ra_gal);
      g_array_append_val (rows, dec_gal);
      g_array_append_val (rows, z_gal);

      galaxy_id++;
    }
  }

  col_names    = g_new0 (gchar *, 7);
  col_names[0] = g_strdup ("galaxy_id");
  col_names[1] = g_strdup ("parent_id");
  col_names[2] = g_strdup ("is_central");
  col_names[3] = g_strdup ("ra");
  col_names[4] = g_strdup ("dec");
  col_names[5] = g_strdup ("z");

  col_types[0] = NCM_CATALOG_COL_TYPE_INT;
  col_types[1] = NCM_CATALOG_COL_TYPE_INT;
  col_types[2] = NCM_CATALOG_COL_TYPE_BOOL;
  col_types[3] = NCM_CATALOG_COL_TYPE_DOUBLE;
  col_types[4] = NCM_CATALOG_COL_TYPE_DOUBLE;
  col_types[5] = NCM_CATALOG_COL_TYPE_DOUBLE;

  members = nc_halo_catalog_new (NC_HALO_CATALOG_KIND_MEMBER, "galaxy_id", "parent_id", galaxy_id, col_names, col_types, 6);

  if (galaxy_id > 0)
    ncm_matrix_set_from_array (ncm_catalog_peek_data (NCM_CATALOG (members)), rows);

  g_strfreev (col_names);
  g_array_unref (rows);

  if (local_dist)
    nc_distance_free (dist);

  return members;
}

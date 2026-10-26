/***************************************************************************
 *            nc_dist_mu_obs.c
 *
 *  Sun Oct 18 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_dist_mu_obs.c
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
 * NcDistMuObs:
 *
 * Catalog of distance moduli with sky positions.
 *
 * A #NcmCatalog with one row per object and the columns #NC_DIST_MU_OBS_COL_RA and
 * #NC_DIST_MU_OBS_COL_DEC (equatorial J2000, degrees), #NC_DIST_MU_OBS_COL_Z (redshift),
 * #NC_DIST_MU_OBS_COL_MU and #NC_DIST_MU_OBS_COL_SIGMA_MU (distance modulus and its
 * uncertainty, magnitudes). Any other column is kept as given, so alternative moduli
 * or identifiers of the source table travel with the object; the integer column
 * #NC_DIST_MU_OBS_COL_SAMPLE, when present, labels the sample of each object, see
 * #NcDistMuNuisance. Creating a catalog without one of the five columns is an error. The provenance of the table goes in
 * the catalog metadata, see ncm_catalog_set_meta().
 *
 * nc_dist_mu_obs_get_direction() gives the unit vector of a row in the equatorial
 * frame, $(\cos\delta\cos\alpha, \cos\delta\sin\alpha, \sin\delta)$, and
 * nc_dist_mu_obs_set_direction() writes the position of a row from a unit vector.
 *
 * nc_dist_mu_obs_new_from_catalog_id() loads a catalog of the NumCosmo data release,
 * see #NcDistMuObsCatalogId.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc/data/nc_dist_mu_obs.h"
#include "ncm/core/ncm_cfg.h"
#include "ncm/core/ncm_serialize.h"
#include "nc/data/nc_data_download.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <math.h>
#endif /* NUMCOSMO_GIR_SCAN */

struct _NcDistMuObs
{
  /*< private >*/
  NcmCatalog parent_instance;
};

G_DEFINE_TYPE (NcDistMuObs, nc_dist_mu_obs, NCM_TYPE_CATALOG)

static const gchar *_nc_dist_mu_obs_required[] = {
  NC_DIST_MU_OBS_COL_RA,
  NC_DIST_MU_OBS_COL_DEC,
  NC_DIST_MU_OBS_COL_Z,
  NC_DIST_MU_OBS_COL_MU,
  NC_DIST_MU_OBS_COL_SIGMA_MU,
  NULL,
};

static void
nc_dist_mu_obs_init (NcDistMuObs *obs)
{
}

static void
_nc_dist_mu_obs_constructed (GObject *object)
{
  NcmCatalog *catalog = NCM_CATALOG (object);
  guint i;

  /* Chain up : start */
  G_OBJECT_CLASS (nc_dist_mu_obs_parent_class)->constructed (object);

  for (i = 0; _nc_dist_mu_obs_required[i] != NULL; i++)
  {
    if (!ncm_catalog_has_column (catalog, _nc_dist_mu_obs_required[i]))
      g_error ("nc_dist_mu_obs: the column `%s' is required.", _nc_dist_mu_obs_required[i]);
  }
}

static void
_nc_dist_mu_obs_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_dist_mu_obs_parent_class)->finalize (object);
}

static void
nc_dist_mu_obs_class_init (NcDistMuObsClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->constructed = &_nc_dist_mu_obs_constructed;
  object_class->finalize    = &_nc_dist_mu_obs_finalize;
}

/**
 * nc_dist_mu_obs_new:
 * @nrows: number of objects
 * @col_names: (nullable) (array zero-terminated=1): column names
 *
 * Creates a new #NcDistMuObs with @nrows rows and the columns @col_names, which
 * must include the five required ones; %NULL gives the five required columns only.
 *
 * Returns: (transfer full): a new #NcDistMuObs
 */
NcDistMuObs *
nc_dist_mu_obs_new (guint nrows, GStrv col_names)
{
  NcDistMuObs *obs = g_object_new (NC_TYPE_DIST_MU_OBS,
                                   "columns", (col_names != NULL) ? col_names : (GStrv) _nc_dist_mu_obs_required,
                                   "len", nrows,
                                   NULL);

  return obs;
}

/**
 * nc_dist_mu_obs_new_full:
 * @nrows: number of objects
 * @col_names: (array zero-terminated=1): column names, including the five required ones
 * @col_types: (array length=n_types): the #NcmCatalogColType of each column
 * @n_types: number of columns
 *
 * Creates a new #NcDistMuObs with @nrows rows, the columns @col_names and their types,
 * see ncm_catalog_new_full().
 *
 * Returns: (transfer full): a new #NcDistMuObs
 */
NcDistMuObs *
nc_dist_mu_obs_new_full (guint nrows, GStrv col_names, const NcmCatalogColType *col_types, guint n_types)
{
  GArray *types = g_array_sized_new (FALSE, TRUE, sizeof (guint32), n_types);
  NcDistMuObs *obs;
  GVariant *var;
  guint i;

  for (i = 0; i < n_types; i++)
  {
    const guint32 t = col_types[i];

    g_array_append_val (types, t);
  }

  var = ncm_cfg_array_to_variant (types, G_VARIANT_TYPE ("u"));
  obs = g_object_new (NC_TYPE_DIST_MU_OBS,
                      "columns", col_names,
                      "col-types", var,
                      "len", nrows,
                      NULL);

  g_variant_unref (var);
  g_array_unref (types);

  return obs;
}

/**
 * nc_dist_mu_obs_catalog_id_get_filename:
 * @id: a #NcDistMuObsCatalogId
 *
 * Gives the local file of catalog @id, downloaded from the NumCosmo data release into
 * the user data directory when it is not there yet, see nc_data_download_release_file().
 *
 * Returns: (transfer full): the full path of the file
 */
gchar *
nc_dist_mu_obs_catalog_id_get_filename (NcDistMuObsCatalogId id)
{
  switch (id)
  {
    case NC_DIST_MU_OBS_CATALOG_CF4_SNIA:

      return nc_data_download_release_file ("nc_dist_mu_obs_cf4_snia.gvar");

      break;
    default:
      g_error ("nc_dist_mu_obs_catalog_id_get_filename: id %d not recognized.", id);

      return NULL;

      break;
  }
}

/**
 * nc_dist_mu_obs_new_from_catalog_id:
 * @id: a #NcDistMuObsCatalogId
 *
 * Loads the catalog @id from the NumCosmo data release, see
 * nc_dist_mu_obs_catalog_id_get_filename(); the provenance of its table is in its
 * metadata, see ncm_catalog_peek_meta().
 *
 * Returns: (transfer full): the catalog
 */
NcDistMuObs *
nc_dist_mu_obs_new_from_catalog_id (NcDistMuObsCatalogId id)
{
  gchar *filename   = nc_dist_mu_obs_catalog_id_get_filename (id);
  NcmSerialize *ser = ncm_serialize_new (NCM_SERIALIZE_OPT_KEEP_SHARED);
  GObject *obj      = ncm_serialize_from_binfile (ser, filename);

  g_free (filename);
  ncm_serialize_free (ser);

  if (!NC_IS_DIST_MU_OBS (obj))
    g_error ("nc_dist_mu_obs_new_from_catalog_id: the file of catalog id %d does not contain a NcDistMuObs.", id);

  return NC_DIST_MU_OBS (obj);
}

/**
 * nc_dist_mu_obs_ref:
 * @obs: a #NcDistMuObs
 *
 * Increases the reference count of @obs by one.
 *
 * Returns: (transfer full): @obs
 */
NcDistMuObs *
nc_dist_mu_obs_ref (NcDistMuObs *obs)
{
  return g_object_ref (obs);
}

/**
 * nc_dist_mu_obs_free:
 * @obs: a #NcDistMuObs
 *
 * Decreases the reference count of @obs by one.
 *
 */
void
nc_dist_mu_obs_free (NcDistMuObs *obs)
{
  g_object_unref (obs);
}

/**
 * nc_dist_mu_obs_clear:
 * @obs: a #NcDistMuObs
 *
 * If *@obs is not %NULL, decreases its reference count and sets *@obs to %NULL.
 *
 */
void
nc_dist_mu_obs_clear (NcDistMuObs **obs)
{
  g_clear_object (obs);
}

/**
 * nc_dist_mu_obs_set_direction:
 * @obs: a #NcDistMuObs
 * @i: row
 * @n: a unit #NcmTriVec in the equatorial frame
 *
 * Writes the right ascension and declination of row @i from @n.
 *
 */
void
nc_dist_mu_obs_set_direction (NcDistMuObs *obs, const guint i, const NcmTriVec *n)
{
  NcmCatalog *catalog = NCM_CATALOG (obs);
  NcmTriVec v;
  gdouble r, ra, dec;

  ncm_trivec_memcpy (&v, n);
  ncm_trivec_get_astro_ra_dec (&v, &r, &ra, &dec);

  ncm_catalog_set (catalog, NC_DIST_MU_OBS_COL_RA, i, ra, NULL);
  ncm_catalog_set (catalog, NC_DIST_MU_OBS_COL_DEC, i, dec, NULL);
}

/**
 * nc_dist_mu_obs_get_direction:
 * @obs: a #NcDistMuObs
 * @i: row
 * @n: (out caller-allocates): a #NcmTriVec
 *
 * Sets @n to the unit vector of row @i in the equatorial frame.
 *
 */
void
nc_dist_mu_obs_get_direction (NcDistMuObs *obs, const guint i, NcmTriVec *n)
{
  NcmCatalog *catalog = NCM_CATALOG (obs);
  const gdouble ra    = ncm_catalog_get (catalog, NC_DIST_MU_OBS_COL_RA, i, NULL);
  const gdouble dec   = ncm_catalog_get (catalog, NC_DIST_MU_OBS_COL_DEC, i, NULL);

  ncm_trivec_set_astro_ra_dec (n, 1.0, ra, dec);
}


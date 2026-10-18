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
 * or identifiers of the source table travel with the object. Creating a catalog
 * without one of the five columns is an error. The provenance of the table goes in
 * the catalog metadata, see ncm_catalog_set_meta().
 *
 * nc_dist_mu_obs_get_direction() gives the unit vector of a row in the equatorial
 * frame, $(\cos\delta\cos\alpha, \cos\delta\sin\alpha, \sin\delta)$, and
 * nc_dist_mu_obs_set_direction() writes the position of a row from a unit vector.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc/data/nc_dist_mu_obs.h"

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


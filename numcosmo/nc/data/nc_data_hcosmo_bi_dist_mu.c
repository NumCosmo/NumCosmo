/***************************************************************************
 *            nc_data_hcosmo_bi_dist_mu.c
 *
 *  Sat Oct 17 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_data_hcosmo_bi_dist_mu.c
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
 * NcDataHCosmoBIDistMu:
 *
 * Distance moduli with directions, for a #NcHCosmoBI.
 *
 * Each object has a redshift $z$ (#NcDataHCosmoBIDistMu:z), a unit vector $\hat p$ in the
 * equatorial frame (row of #NcDataHCosmoBIDistMu:p), a distance modulus and its
 * uncertainty (the mean and standard deviation of #NcmDataGaussDiag). The model value is
 * nc_hcosmo_bi_dist_dmodulus_sky() plus $5\log_{10}$ nc_hcosmo_bi_RH_Mpc(), computed with
 * #NcDataHCosmoBIDistMu:dist, whose #NcHCosmoBIDist:zmax must cover the redshifts.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc/data/nc_data_hcosmo_bi_dist_mu.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <math.h>
#endif /* NUMCOSMO_GIR_SCAN */

struct _NcDataHCosmoBIDistMu
{
  /*< private >*/
  NcmDataGaussDiag parent_instance;
  NcHCosmoBIDist *dist;
  NcmVector *z;
  NcmMatrix *p;
};

enum
{
  PROP_0,
  PROP_DIST,
  PROP_Z,
  PROP_P,
  PROP_SIZE,
};

G_DEFINE_TYPE (NcDataHCosmoBIDistMu, nc_data_hcosmo_bi_dist_mu, NCM_TYPE_DATA_GAUSS_DIAG)

static void
nc_data_hcosmo_bi_dist_mu_init (NcDataHCosmoBIDistMu *dmu)
{
  dmu->dist = NULL;
  dmu->z    = NULL;
  dmu->p    = NULL;
}

static void
_nc_data_hcosmo_bi_dist_mu_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcDataHCosmoBIDistMu *dmu = NC_DATA_HCOSMO_BI_DIST_MU (object);

  g_return_if_fail (NC_IS_DATA_HCOSMO_BI_DIST_MU (object));

  switch (prop_id)
  {
    case PROP_DIST:
      nc_data_hcosmo_bi_dist_mu_set_dist (dmu, g_value_get_object (value));
      break;
    case PROP_Z:
      ncm_vector_substitute (&dmu->z, g_value_get_object (value), TRUE);
      break;
    case PROP_P:
      ncm_matrix_substitute (&dmu->p, g_value_get_object (value), TRUE);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_data_hcosmo_bi_dist_mu_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcDataHCosmoBIDistMu *dmu = NC_DATA_HCOSMO_BI_DIST_MU (object);

  g_return_if_fail (NC_IS_DATA_HCOSMO_BI_DIST_MU (object));

  switch (prop_id)
  {
    case PROP_DIST:
      g_value_set_object (value, dmu->dist);
      break;
    case PROP_Z:
      g_value_set_object (value, dmu->z);
      break;
    case PROP_P:
      g_value_set_object (value, dmu->p);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_data_hcosmo_bi_dist_mu_dispose (GObject *object)
{
  NcDataHCosmoBIDistMu *dmu = NC_DATA_HCOSMO_BI_DIST_MU (object);

  nc_hcosmo_bi_dist_clear (&dmu->dist);
  ncm_vector_clear (&dmu->z);
  ncm_matrix_clear (&dmu->p);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_data_hcosmo_bi_dist_mu_parent_class)->dispose (object);
}

static void
_nc_data_hcosmo_bi_dist_mu_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_data_hcosmo_bi_dist_mu_parent_class)->finalize (object);
}

static void _nc_data_hcosmo_bi_dist_mu_prepare (NcmData *data, NcmMSet *mset);
static void _nc_data_hcosmo_bi_dist_mu_mean_func (NcmDataGaussDiag *diag, NcmMSet *mset, NcmVector *vp);
static void _nc_data_hcosmo_bi_dist_mu_set_size (NcmDataGaussDiag *diag, guint np);

static void
nc_data_hcosmo_bi_dist_mu_class_init (NcDataHCosmoBIDistMuClass *klass)
{
  GObjectClass *object_class        = G_OBJECT_CLASS (klass);
  NcmDataClass *data_class          = NCM_DATA_CLASS (klass);
  NcmDataGaussDiagClass *diag_class = NCM_DATA_GAUSS_DIAG_CLASS (klass);

  object_class->set_property = &_nc_data_hcosmo_bi_dist_mu_set_property;
  object_class->get_property = &_nc_data_hcosmo_bi_dist_mu_get_property;
  object_class->dispose      = &_nc_data_hcosmo_bi_dist_mu_dispose;
  object_class->finalize     = &_nc_data_hcosmo_bi_dist_mu_finalize;

  /**
   * NcDataHCosmoBIDistMu:dist:
   *
   * The #NcHCosmoBIDist used for the model values.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_DIST,
                                   g_param_spec_object ("dist",
                                                        NULL,
                                                        "Distance object",
                                                        NC_TYPE_HCOSMO_BI_DIST,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcDataHCosmoBIDistMu:z:
   *
   * Redshifts of the objects.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_Z,
                                   g_param_spec_object ("z",
                                                        NULL,
                                                        "Data redshifts",
                                                        NCM_TYPE_VECTOR,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcDataHCosmoBIDistMu:p:
   *
   * Unit vectors of the objects in the equatorial frame, one per row.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_P,
                                   g_param_spec_object ("p",
                                                        NULL,
                                                        "Data directions",
                                                        NCM_TYPE_MATRIX,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  data_class->prepare   = &_nc_data_hcosmo_bi_dist_mu_prepare;
  diag_class->mean_func = &_nc_data_hcosmo_bi_dist_mu_mean_func;
  diag_class->set_size  = &_nc_data_hcosmo_bi_dist_mu_set_size;
}

static void
_nc_data_hcosmo_bi_dist_mu_prepare (NcmData *data, NcmMSet *mset)
{
  NcDataHCosmoBIDistMu *dmu = NC_DATA_HCOSMO_BI_DIST_MU (data);
  NcHCosmoBI *cosmo         = NC_HCOSMO_BI (ncm_mset_peek (mset, nc_hcosmo_bi_id ()));

  nc_hcosmo_bi_dist_prepare_if_needed (dmu->dist, cosmo);
}

static void
_nc_data_hcosmo_bi_dist_mu_mean_func (NcmDataGaussDiag *diag, NcmMSet *mset, NcmVector *vp)
{
  NcDataHCosmoBIDistMu *dmu = NC_DATA_HCOSMO_BI_DIST_MU (diag);
  NcHCosmoBI *cosmo         = NC_HCOSMO_BI (ncm_mset_peek (mset, nc_hcosmo_bi_id ()));
  const gdouble f_lnRH      = 5.0 * log10 (nc_hcosmo_bi_RH_Mpc (cosmo));
  const guint np            = ncm_data_gauss_diag_get_size (diag);
  guint i;

  for (i = 0; i < np; i++)
  {
    const gdouble z     = ncm_vector_get (dmu->z, i);
    const gdouble dmu_i = nc_hcosmo_bi_dist_dmodulus_sky (dmu->dist, z,
                                                          ncm_matrix_get (dmu->p, i, 0),
                                                          ncm_matrix_get (dmu->p, i, 1),
                                                          ncm_matrix_get (dmu->p, i, 2));

    ncm_vector_set (vp, i, dmu_i + f_lnRH);
  }
}

static void
_nc_data_hcosmo_bi_dist_mu_set_size (NcmDataGaussDiag *diag, guint np)
{
  NcDataHCosmoBIDistMu *dmu = NC_DATA_HCOSMO_BI_DIST_MU (diag);
  const guint cnp           = ncm_data_gauss_diag_get_size (diag);

  if ((np == 0) || (np != cnp))
  {
    ncm_vector_clear (&dmu->z);
    ncm_matrix_clear (&dmu->p);
  }

  if ((np != 0) && (np != cnp))
  {
    dmu->z = ncm_vector_new (np);
    dmu->p = ncm_matrix_new (np, 3);
  }

  /* Chain up : end */
  NCM_DATA_GAUSS_DIAG_CLASS (nc_data_hcosmo_bi_dist_mu_parent_class)->set_size (diag, np);
}

/**
 * nc_data_hcosmo_bi_dist_mu_new:
 * @dist: (nullable): a #NcHCosmoBIDist
 *
 * Creates an empty #NcDataHCosmoBIDistMu; nc_data_hcosmo_bi_dist_mu_set_data() or
 * nc_data_hcosmo_bi_dist_mu_set_obs() fills it, and nc_data_hcosmo_bi_dist_mu_set_dist()
 * sets the calculator when @dist is %NULL.
 *
 * Returns: (transfer full): a new #NcDataHCosmoBIDistMu
 */
NcDataHCosmoBIDistMu *
nc_data_hcosmo_bi_dist_mu_new (NcHCosmoBIDist *dist)
{
  NcDataHCosmoBIDistMu *dmu = g_object_new (NC_TYPE_DATA_HCOSMO_BI_DIST_MU,
                                            "dist", dist,
                                            NULL);

  return dmu;
}

/**
 * nc_data_hcosmo_bi_dist_mu_ref:
 * @dmu: a #NcDataHCosmoBIDistMu
 *
 * Increases the reference count of @dmu by one.
 *
 * Returns: (transfer full): @dmu
 */
NcDataHCosmoBIDistMu *
nc_data_hcosmo_bi_dist_mu_ref (NcDataHCosmoBIDistMu *dmu)
{
  return g_object_ref (dmu);
}

/**
 * nc_data_hcosmo_bi_dist_mu_free:
 * @dmu: a #NcDataHCosmoBIDistMu
 *
 * Decreases the reference count of @dmu by one.
 *
 */
void
nc_data_hcosmo_bi_dist_mu_free (NcDataHCosmoBIDistMu *dmu)
{
  g_object_unref (dmu);
}

/**
 * nc_data_hcosmo_bi_dist_mu_clear:
 * @dmu: a #NcDataHCosmoBIDistMu
 *
 * If *@dmu is not %NULL, decreases its reference count and sets *@dmu to %NULL.
 *
 */
void
nc_data_hcosmo_bi_dist_mu_clear (NcDataHCosmoBIDistMu **dmu)
{
  g_clear_object (dmu);
}

/**
 * nc_data_hcosmo_bi_dist_mu_set_dist:
 * @dmu: a #NcDataHCosmoBIDistMu
 * @dist: (nullable): a #NcHCosmoBIDist
 *
 * Sets #NcDataHCosmoBIDistMu:dist.
 *
 */
void
nc_data_hcosmo_bi_dist_mu_set_dist (NcDataHCosmoBIDistMu *dmu, NcHCosmoBIDist *dist)
{
  nc_hcosmo_bi_dist_clear (&dmu->dist);

  if (dist != NULL)
    dmu->dist = nc_hcosmo_bi_dist_ref (dist);
}

/**
 * nc_data_hcosmo_bi_dist_mu_set_data:
 * @dmu: a #NcDataHCosmoBIDistMu
 * @z: redshifts
 * @p: unit vectors in the equatorial frame, one per row
 * @mu: distance moduli
 * @sigma: uncertainties of the distance moduli
 *
 * Copies the data into @dmu and marks it initialized. The four arguments must have the
 * same number of objects, and @p three columns.
 *
 */
void
nc_data_hcosmo_bi_dist_mu_set_data (NcDataHCosmoBIDistMu *dmu, NcmVector *z, NcmMatrix *p, NcmVector *mu, NcmVector *sigma)
{
  NcmDataGaussDiag *diag = NCM_DATA_GAUSS_DIAG (dmu);
  const guint np         = ncm_vector_len (z);

  g_assert_cmpuint (ncm_matrix_nrows (p), ==, np);
  g_assert_cmpuint (ncm_matrix_ncols (p), ==, 3);
  g_assert_cmpuint (ncm_vector_len (mu), ==, np);
  g_assert_cmpuint (ncm_vector_len (sigma), ==, np);

  ncm_data_gauss_diag_set_size (diag, np);

  ncm_vector_memcpy (dmu->z, z);
  ncm_matrix_memcpy (dmu->p, p);
  ncm_vector_memcpy (ncm_data_gauss_diag_peek_mean (diag), mu);
  ncm_vector_memcpy (ncm_data_gauss_diag_peek_std (diag), sigma);

  ncm_data_set_init (NCM_DATA (dmu), TRUE);
}

/**
 * nc_data_hcosmo_bi_dist_mu_set_obs:
 * @dmu: a #NcDataHCosmoBIDistMu
 * @obs: a #NcDistMuObs
 * @mu_col: (nullable): column of the distance moduli
 * @sigma_col: (nullable): column of their uncertainties
 * @z_min: lower redshift limit
 * @z_max: upper redshift limit
 * @error: a #GError
 *
 * Fills @dmu with the objects of @obs whose redshift lies in [@z_min, @z_max], taking
 * the moduli from @mu_col and their uncertainties from @sigma_col; %NULL selects
 * #NC_DIST_MU_OBS_COL_MU and #NC_DIST_MU_OBS_COL_SIGMA_MU. A column missing from @obs
 * sets @error. The directions come from nc_dist_mu_obs_get_direction().
 *
 */
void
nc_data_hcosmo_bi_dist_mu_set_obs (NcDataHCosmoBIDistMu *dmu, NcDistMuObs *obs, const gchar *mu_col, const gchar *sigma_col, const gdouble z_min, const gdouble z_max, GError **error)
{
  NcmCatalog *catalog = NCM_CATALOG (obs);
  const guint nrows   = ncm_catalog_len (catalog);
  const gchar *mcol   = (mu_col != NULL) ? mu_col : NC_DIST_MU_OBS_COL_MU;
  const gchar *scol   = (sigma_col != NULL) ? sigma_col : NC_DIST_MU_OBS_COL_SIGMA_MU;
  guint np            = 0;
  NcmVector *z, *mu, *sigma;
  NcmMatrix *p;
  guint i, k;

  if (!ncm_catalog_has_column (catalog, mcol))
  {
    ncm_util_set_or_call_error (error, NCM_CATALOG_ERROR, NCM_CATALOG_ERROR_COLUMN_NOT_FOUND,
                                "nc_data_hcosmo_bi_dist_mu_set_obs: column `%s' not found.", mcol);

    return;
  }

  if (!ncm_catalog_has_column (catalog, scol))
  {
    ncm_util_set_or_call_error (error, NCM_CATALOG_ERROR, NCM_CATALOG_ERROR_COLUMN_NOT_FOUND,
                                "nc_data_hcosmo_bi_dist_mu_set_obs: column `%s' not found.", scol);

    return;
  }

  for (i = 0; i < nrows; i++)
  {
    const gdouble zi = ncm_catalog_get (catalog, NC_DIST_MU_OBS_COL_Z, i, NULL);

    if ((zi >= z_min) && (zi <= z_max))
      np++;
  }

  z     = ncm_vector_new (np);
  mu    = ncm_vector_new (np);
  sigma = ncm_vector_new (np);
  p     = ncm_matrix_new (np, 3);

  for (i = 0, k = 0; i < nrows; i++)
  {
    const gdouble zi = ncm_catalog_get (catalog, NC_DIST_MU_OBS_COL_Z, i, NULL);
    NcmTriVec n;

    if (!((zi >= z_min) && (zi <= z_max)))
      continue;

    nc_dist_mu_obs_get_direction (obs, i, &n);
    ncm_vector_set (z, k, zi);
    ncm_vector_set (mu, k, ncm_catalog_get (catalog, mcol, i, NULL));
    ncm_vector_set (sigma, k, ncm_catalog_get (catalog, scol, i, NULL));
    ncm_matrix_set (p, k, 0, n.c[0]);
    ncm_matrix_set (p, k, 1, n.c[1]);
    ncm_matrix_set (p, k, 2, n.c[2]);
    k++;
  }

  nc_data_hcosmo_bi_dist_mu_set_data (dmu, z, p, mu, sigma);

  ncm_vector_free (z);
  ncm_vector_free (mu);
  ncm_vector_free (sigma);
  ncm_matrix_free (p);
}

/**
 * nc_data_hcosmo_bi_dist_mu_get_obs:
 * @dmu: a #NcDataHCosmoBIDistMu
 *
 * Builds a #NcDistMuObs with the five required columns from the current data of
 * @dmu: the moduli are the current mean, so a resampled @dmu gives a mock catalog.
 *
 * Returns: (transfer full): a new #NcDistMuObs
 */
NcDistMuObs *
nc_data_hcosmo_bi_dist_mu_get_obs (NcDataHCosmoBIDistMu *dmu)
{
  NcmDataGaussDiag *diag = NCM_DATA_GAUSS_DIAG (dmu);
  const guint np         = ncm_data_gauss_diag_get_size (diag);
  NcDistMuObs *obs       = nc_dist_mu_obs_new (np, NULL);
  NcmCatalog *catalog    = NCM_CATALOG (obs);
  NcmVector *mean        = ncm_data_gauss_diag_peek_mean (diag);
  NcmVector *std         = ncm_data_gauss_diag_peek_std (diag);
  guint i;

  for (i = 0; i < np; i++)
  {
    NcmTriVec n = {
      {
        ncm_matrix_get (dmu->p, i, 0), ncm_matrix_get (dmu->p, i, 1), ncm_matrix_get (dmu->p, i, 2)
      }
    };

    nc_dist_mu_obs_set_direction (obs, i, &n);
    ncm_catalog_set (catalog, NC_DIST_MU_OBS_COL_Z, i, ncm_vector_get (dmu->z, i), NULL);
    ncm_catalog_set (catalog, NC_DIST_MU_OBS_COL_MU, i, ncm_vector_get (mean, i), NULL);
    ncm_catalog_set (catalog, NC_DIST_MU_OBS_COL_SIGMA_MU, i, ncm_vector_get (std, i), NULL);
  }

  return obs;
}

/**
 * nc_data_hcosmo_bi_dist_mu_peek_z:
 * @dmu: a #NcDataHCosmoBIDistMu
 *
 * Returns: (transfer none): the redshifts
 */
NcmVector *
nc_data_hcosmo_bi_dist_mu_peek_z (NcDataHCosmoBIDistMu *dmu)
{
  return dmu->z;
}

/**
 * nc_data_hcosmo_bi_dist_mu_peek_p:
 * @dmu: a #NcDataHCosmoBIDistMu
 *
 * Returns: (transfer none): the unit vectors in the equatorial frame
 */
NcmMatrix *
nc_data_hcosmo_bi_dist_mu_peek_p (NcDataHCosmoBIDistMu *dmu)
{
  return dmu->p;
}


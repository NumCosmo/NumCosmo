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
 * equatorial frame (row of #NcDataHCosmoBIDistMu:p), a distance modulus (the mean of
 * #NcmDataGaussDiag), its measured uncertainty $\sigma_\mu$ (#NcDataHCosmoBIDistMu:sigma-mu)
 * and a sample index (#NcDataHCosmoBIDistMu:sample, $-1$ for none). The model value is
 * nc_hcosmo_bi_dist_dmodulus_sky() plus $5\log_{10}$ nc_hcosmo_bi_RH_Mpc(), computed with
 * #NcDataHCosmoBIDistMu:dist, whose #NcHCosmoBIDist:zmax must cover the redshifts.
 *
 * When the model set holds a #NcDistMuNuisance, the model value of an object of sample
 * $s$ gains the offset $\delta\mu_s$ (objects whose sample index is $-1$ or at least the
 * number of offsets get none), and the standard deviation becomes
 * $\sqrt{\sigma_\mu^2 + \sigma_{\mu,v}(z)^2}$, see nc_dist_mu_nuisance_sigma_mu_v().
 * Without it the standard deviation is $\sigma_\mu$.
 *
 * When the model set holds a #NcObserverVelocity, the model value of an object is that of
 * the rest frame at nc_observer_velocity_z_rest() minus $5\log_{10}$ of
 * nc_observer_velocity_doppler(), the observer's motion relative to the rest frame of the
 * model; the velocity term of the standard deviation keeps the observed redshift. The
 * rest-frame redshift of every object must stay positive and within
 * #NcHCosmoBIDist:zmax, which bounds the velocity by the smallest redshift of the data.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc/data/nc_data_hcosmo_bi_dist_mu.h"
#include "nc/data/nc_observer_velocity.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <math.h>
#include <gsl/gsl_math.h>
#endif /* NUMCOSMO_GIR_SCAN */

struct _NcDataHCosmoBIDistMu
{
  /*< private >*/
  NcmDataGaussDiag parent_instance;
  NcHCosmoBIDist *dist;
  NcmVector *z;
  NcmMatrix *p;
  NcmVector *sigma_mu;
  NcmVector *sample;
  gdouble sigma_v_last;
  gboolean sigma_stale;
};

enum
{
  PROP_0,
  PROP_DIST,
  PROP_Z,
  PROP_P,
  PROP_SIGMA_MU,
  PROP_SAMPLE,
  PROP_SIZE,
};

G_DEFINE_TYPE (NcDataHCosmoBIDistMu, nc_data_hcosmo_bi_dist_mu, NCM_TYPE_DATA_GAUSS_DIAG)

static void
nc_data_hcosmo_bi_dist_mu_init (NcDataHCosmoBIDistMu *dmu)
{
  dmu->dist         = NULL;
  dmu->z            = NULL;
  dmu->p            = NULL;
  dmu->sigma_mu     = NULL;
  dmu->sample       = NULL;
  dmu->sigma_v_last = 0.0;
  dmu->sigma_stale  = TRUE;
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
    case PROP_SIGMA_MU:
      ncm_vector_substitute (&dmu->sigma_mu, g_value_get_object (value), TRUE);
      dmu->sigma_stale = TRUE;
      break;
    case PROP_SAMPLE:
      ncm_vector_substitute (&dmu->sample, g_value_get_object (value), TRUE);
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
    case PROP_SIGMA_MU:
      g_value_set_object (value, dmu->sigma_mu);
      break;
    case PROP_SAMPLE:
      g_value_set_object (value, dmu->sample);
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
  ncm_vector_clear (&dmu->sigma_mu);
  ncm_vector_clear (&dmu->sample);

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
static gboolean _nc_data_hcosmo_bi_dist_mu_sigma_func (NcmDataGaussDiag *diag, NcmMSet *mset, NcmVector *sigma);
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

  /**
   * NcDataHCosmoBIDistMu:sigma-mu:
   *
   * Measured uncertainties of the distance moduli.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_SIGMA_MU,
                                   g_param_spec_object ("sigma-mu",
                                                        NULL,
                                                        "Measured uncertainties of the moduli",
                                                        NCM_TYPE_VECTOR,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcDataHCosmoBIDistMu:sample:
   *
   * Sample index of each object, $-1$ for none.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_SAMPLE,
                                   g_param_spec_object ("sample",
                                                        NULL,
                                                        "Sample index of each object",
                                                        NCM_TYPE_VECTOR,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  data_class->prepare    = &_nc_data_hcosmo_bi_dist_mu_prepare;
  diag_class->mean_func  = &_nc_data_hcosmo_bi_dist_mu_mean_func;
  diag_class->sigma_func = &_nc_data_hcosmo_bi_dist_mu_sigma_func;
  diag_class->set_size   = &_nc_data_hcosmo_bi_dist_mu_set_size;
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
  NcDataHCosmoBIDistMu *dmu  = NC_DATA_HCOSMO_BI_DIST_MU (diag);
  NcHCosmoBI *cosmo          = NC_HCOSMO_BI (ncm_mset_peek (mset, nc_hcosmo_bi_id ()));
  NcDistMuNuisance *nuisance = NC_DIST_MU_NUISANCE (ncm_mset_peek (mset, nc_dist_mu_nuisance_id ()));
  NcObserverVelocity *ov     = NC_OBSERVER_VELOCITY (ncm_mset_peek (mset, nc_observer_velocity_id ()));
  const gdouble f_lnRH       = 5.0 * log10 (nc_hcosmo_bi_RH_Mpc (cosmo));
  const guint np             = ncm_data_gauss_diag_get_size (diag);
  const gint nsamples        = (nuisance != NULL) ? (gint) nc_dist_mu_nuisance_get_nsamples (nuisance) : 0;
  guint i;

  for (i = 0; i < np; i++)
  {
    const gdouble px  = ncm_matrix_get (dmu->p, i, 0);
    const gdouble py  = ncm_matrix_get (dmu->p, i, 1);
    const gdouble pz  = ncm_matrix_get (dmu->p, i, 2);
    const gint s      = (gint) ncm_vector_get (dmu->sample, i);
    const gdouble off = ((s >= 0) && (s < nsamples)) ? nc_dist_mu_nuisance_get_dmu (nuisance, s) : 0.0;
    gdouble z         = ncm_vector_get (dmu->z, i);
    gdouble f_lnD     = 0.0;

    if (ov != NULL)
    {
      const gdouble D = nc_observer_velocity_doppler (ov, px, py, pz);

      z     = (1.0 + z) * D - 1.0;
      f_lnD = 5.0 * log10 (D);

      if (z <= 0.0)
        g_error ("_nc_data_hcosmo_bi_dist_mu_mean_func: object %u has a rest-frame redshift %g <= 0 for the "
                 "observer velocity %g km/s; raise the lower redshift limit of the data or bound the velocity.",
                 i, z, nc_observer_velocity_speed (ov));
    }

    ncm_vector_set (vp, i, nc_hcosmo_bi_dist_dmodulus_sky (dmu->dist, z, px, py, pz) - f_lnD + f_lnRH + off);
  }
}

/*
 * Fills @sigma with sqrt(sigma_mu^2 + sigma_mu_v(z)^2) when the model set holds a
 * NcDistMuNuisance and with sigma_mu otherwise; returns TRUE when it wrote @sigma, which
 * happens when sigma_v changed since the last call or the data changed since.
 */
static gboolean
_nc_data_hcosmo_bi_dist_mu_sigma_func (NcmDataGaussDiag *diag, NcmMSet *mset, NcmVector *sigma)
{
  NcDataHCosmoBIDistMu *dmu  = NC_DATA_HCOSMO_BI_DIST_MU (diag);
  NcDistMuNuisance *nuisance = NC_DIST_MU_NUISANCE (ncm_mset_peek (mset, nc_dist_mu_nuisance_id ()));
  const gdouble sigma_v      = (nuisance != NULL) ? nc_dist_mu_nuisance_get_sigma_v (nuisance) : 0.0;
  const guint np             = ncm_data_gauss_diag_get_size (diag);
  guint i;

  if (!dmu->sigma_stale && (sigma_v == dmu->sigma_v_last))
    return FALSE;

  for (i = 0; i < np; i++)
  {
    const gdouble sm = ncm_vector_get (dmu->sigma_mu, i);

    if (sigma_v == 0.0)
    {
      ncm_vector_set (sigma, i, sm);
    }
    else
    {
      const gdouble sv = nc_dist_mu_nuisance_sigma_mu_v (nuisance, ncm_vector_get (dmu->z, i));

      ncm_vector_set (sigma, i, hypot (sm, sv));
    }
  }

  dmu->sigma_v_last = sigma_v;
  dmu->sigma_stale  = FALSE;

  return TRUE;
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
    ncm_vector_clear (&dmu->sigma_mu);
    ncm_vector_clear (&dmu->sample);
  }

  if ((np != 0) && (np != cnp))
  {
    dmu->z        = ncm_vector_new (np);
    dmu->p        = ncm_matrix_new (np, 3);
    dmu->sigma_mu = ncm_vector_new (np);
    dmu->sample   = ncm_vector_new (np);
    ncm_vector_set_all (dmu->sample, -1.0);
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
 * Copies the data into @dmu, sets every sample index to $-1$ and marks @dmu
 * initialized. The four arguments must have the same number of objects, and @p three
 * columns.
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
  ncm_vector_memcpy (dmu->sigma_mu, sigma);
  ncm_vector_set_all (dmu->sample, -1.0);
  ncm_vector_memcpy (ncm_data_gauss_diag_peek_mean (diag), mu);
  ncm_vector_memcpy (ncm_data_gauss_diag_peek_std (diag), sigma);
  dmu->sigma_stale = TRUE;

  ncm_data_set_init (NCM_DATA (dmu), TRUE);
}

/**
 * nc_data_hcosmo_bi_dist_mu_set_samples:
 * @dmu: a #NcDataHCosmoBIDistMu
 * @sample: sample index of each object, $-1$ for none
 *
 * Copies @sample into @dmu, which must already hold the data; @sample must have one
 * entry per object.
 *
 */
void
nc_data_hcosmo_bi_dist_mu_set_samples (NcDataHCosmoBIDistMu *dmu, NcmVector *sample)
{
  g_assert_cmpuint (ncm_vector_len (sample), ==, ncm_data_gauss_diag_get_size (NCM_DATA_GAUSS_DIAG (dmu)));

  ncm_vector_memcpy (dmu->sample, sample);
}

/**
 * nc_data_hcosmo_bi_dist_mu_set_obs:
 * @dmu: a #NcDataHCosmoBIDistMu
 * @obs: a #NcDistMuObs
 * @mu_col: (nullable): column of the distance moduli
 * @sigma_col: (nullable): column of their uncertainties
 * @sample_col: (nullable): integer column of the sample indices
 * @z_min: lower redshift limit
 * @z_max: upper redshift limit
 * @error: a #GError
 *
 * Fills @dmu with the objects of @obs whose redshift lies in [@z_min, @z_max] and whose
 * modulus and uncertainty in @mu_col and @sigma_col are finite; %NULL selects
 * #NC_DIST_MU_OBS_COL_MU and #NC_DIST_MU_OBS_COL_SIGMA_MU, so a per-survey column with
 * empty entries selects the objects of that survey. The sample indices come from
 * @sample_col, or are $-1$ when it is %NULL. A column missing from @obs sets @error.
 * The directions come from nc_dist_mu_obs_get_direction().
 *
 */
void
nc_data_hcosmo_bi_dist_mu_set_obs (NcDataHCosmoBIDistMu *dmu, NcDistMuObs *obs, const gchar *mu_col, const gchar *sigma_col, const gchar *sample_col, const gdouble z_min, const gdouble z_max, GError **error)
{
  NcmCatalog *catalog = NCM_CATALOG (obs);
  const guint nrows   = ncm_catalog_len (catalog);
  const gchar *mcol   = (mu_col != NULL) ? mu_col : NC_DIST_MU_OBS_COL_MU;
  const gchar *scol   = (sigma_col != NULL) ? sigma_col : NC_DIST_MU_OBS_COL_SIGMA_MU;
  guint np            = 0;
  NcmVector *z, *mu, *sigma, *sample;
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

  if ((sample_col != NULL) && !ncm_catalog_has_column (catalog, sample_col))
  {
    ncm_util_set_or_call_error (error, NCM_CATALOG_ERROR, NCM_CATALOG_ERROR_COLUMN_NOT_FOUND,
                                "nc_data_hcosmo_bi_dist_mu_set_obs: column `%s' not found.", sample_col);

    return;
  }

  for (i = 0; i < nrows; i++)
  {
    const gdouble zi  = ncm_catalog_get (catalog, NC_DIST_MU_OBS_COL_Z, i, NULL);
    const gdouble mui = ncm_catalog_get (catalog, mcol, i, NULL);
    const gdouble si  = ncm_catalog_get (catalog, scol, i, NULL);

    if ((zi >= z_min) && (zi <= z_max) && gsl_finite (mui) && gsl_finite (si))
      np++;
  }

  z      = ncm_vector_new (np);
  mu     = ncm_vector_new (np);
  sigma  = ncm_vector_new (np);
  sample = ncm_vector_new (np);
  p      = ncm_matrix_new (np, 3);

  for (i = 0, k = 0; i < nrows; i++)
  {
    const gdouble zi  = ncm_catalog_get (catalog, NC_DIST_MU_OBS_COL_Z, i, NULL);
    const gdouble mui = ncm_catalog_get (catalog, mcol, i, NULL);
    const gdouble si  = ncm_catalog_get (catalog, scol, i, NULL);
    NcmTriVec n;

    if (!((zi >= z_min) && (zi <= z_max) && gsl_finite (mui) && gsl_finite (si)))
      continue;

    nc_dist_mu_obs_get_direction (obs, i, &n);
    ncm_vector_set (z, k, zi);
    ncm_vector_set (mu, k, mui);
    ncm_vector_set (sigma, k, si);
    ncm_vector_set (sample, k, (sample_col != NULL) ? (gdouble) ncm_catalog_get_int (catalog, sample_col, i, NULL) : -1.0);
    ncm_matrix_set (p, k, 0, n.c[0]);
    ncm_matrix_set (p, k, 1, n.c[1]);
    ncm_matrix_set (p, k, 2, n.c[2]);
    k++;
  }

  nc_data_hcosmo_bi_dist_mu_set_data (dmu, z, p, mu, sigma);
  nc_data_hcosmo_bi_dist_mu_set_samples (dmu, sample);

  ncm_vector_free (z);
  ncm_vector_free (mu);
  ncm_vector_free (sigma);
  ncm_vector_free (sample);
  ncm_matrix_free (p);
}

/**
 * nc_data_hcosmo_bi_dist_mu_get_obs:
 * @dmu: a #NcDataHCosmoBIDistMu
 *
 * Builds a #NcDistMuObs with the five required columns and #NC_DIST_MU_OBS_COL_SAMPLE
 * from the current data of @dmu: the moduli are the current mean, so a resampled @dmu
 * gives a mock catalog; the uncertainties are the measured ones, #NcDataHCosmoBIDistMu:sigma-mu.
 *
 * Returns: (transfer full): a new #NcDistMuObs
 */
NcDistMuObs *
nc_data_hcosmo_bi_dist_mu_get_obs (NcDataHCosmoBIDistMu *dmu)
{
  NcmDataGaussDiag *diag          = NCM_DATA_GAUSS_DIAG (dmu);
  const guint np                  = ncm_data_gauss_diag_get_size (diag);
  const gchar *cols[]             = {NC_DIST_MU_OBS_COL_RA, NC_DIST_MU_OBS_COL_DEC, NC_DIST_MU_OBS_COL_Z, NC_DIST_MU_OBS_COL_MU, NC_DIST_MU_OBS_COL_SIGMA_MU, NC_DIST_MU_OBS_COL_SAMPLE, NULL};
  const NcmCatalogColType types[] = {NCM_CATALOG_COL_TYPE_DOUBLE, NCM_CATALOG_COL_TYPE_DOUBLE, NCM_CATALOG_COL_TYPE_DOUBLE, NCM_CATALOG_COL_TYPE_DOUBLE, NCM_CATALOG_COL_TYPE_DOUBLE, NCM_CATALOG_COL_TYPE_INT};
  NcDistMuObs *obs                = nc_dist_mu_obs_new_full (np, (GStrv) cols, types, 6);
  NcmCatalog *catalog             = NCM_CATALOG (obs);
  NcmVector *mean                 = ncm_data_gauss_diag_peek_mean (diag);
  NcmVector *std                  = dmu->sigma_mu;
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
    ncm_catalog_set_int (catalog, NC_DIST_MU_OBS_COL_SAMPLE, i, (gint64) ncm_vector_get (dmu->sample, i), NULL);
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

/**
 * nc_data_hcosmo_bi_dist_mu_peek_sigma_mu:
 * @dmu: a #NcDataHCosmoBIDistMu
 *
 * Returns: (transfer none): the measured uncertainties of the moduli
 */
NcmVector *
nc_data_hcosmo_bi_dist_mu_peek_sigma_mu (NcDataHCosmoBIDistMu *dmu)
{
  return dmu->sigma_mu;
}

/**
 * nc_data_hcosmo_bi_dist_mu_peek_samples:
 * @dmu: a #NcDataHCosmoBIDistMu
 *
 * Returns: (transfer none): the sample index of each object, $-1$ for none
 */
NcmVector *
nc_data_hcosmo_bi_dist_mu_peek_samples (NcDataHCosmoBIDistMu *dmu)
{
  return dmu->sample;
}


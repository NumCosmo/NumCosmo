/***************************************************************************
 *            nc_galaxy_sd_obs_redshift_gauss.c
 *
 *  Thu Aug 1 20:03:55 2024
 *  Copyright  2024  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 *  Copyright  2024  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_galaxy_sd_obs_redshift_gauss.c
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
 * NcGalaxySDObsRedshiftGauss:
 *
 * Class describing photometric redshift observations with gaussian errors.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "gsl/gsl_randist.h"
#include "gsl/gsl_sf_erf.h"
#include "galaxy/nc_galaxy_sd_obs_redshift.h"
#include "galaxy/nc_galaxy_sd_obs_redshift_gauss.h"
#include "galaxy/nc_galaxy_sd_true_redshift.h"
#include "galaxy/nc_galaxy_sd_true_redshift_lsst_srd.h"
#include "math/ncm_vector.h"
#include "math/ncm_rng.h"
#include "math/ncm_dtuple.h"
#include "math/ncm_spline_cubic_notaknot.h"
#include "math/ncm_integral1d_ptr.h"
#include "math/ncm_stats_dist1d_spline.h"
#include "math/ncm_spline_func.h"
#include "math/ncm_util.h"

typedef struct _NcGalaxySDObsRedshiftGaussPrivate
{
  NcGalaxySDTrueRedshift *sdz;
  gdouble zp_min;
  gdouble zp_max;
  gboolean use_true_z;
  gdouble bin_sigma0;
  gdouble reltol;
  gdouble zp_support_max;
  /* Cache for p(zp) distribution */
  NcmSpline *pzp_spline;
  NcmStatsDist1d *stats;
  /* Cache for p(z|zp_min,zp_max) distribution */
  NcmSpline *pz_given_zp_spline;
} NcGalaxySDObsRedshiftGaussPrivate;

struct _NcGalaxySDObsRedshiftGauss
{
  NcGalaxySDObsRedshift parent_instance;
};

typedef struct _NcGalaxySDObsRedshiftGaussData
{
  gdouble zp;
  gdouble sigma;
  gdouble sigma0;
} NcGalaxySDObsRedshiftGaussData;

enum
{
  PROP_0,
  PROP_ZP_LIM,
  PROP_USE_TRUE_Z,
  PROP_BIN_SIGMA0,
  PROP_RELTOL,
  PROP_ZP_SUPPORT_MAX,
  PROP_LEN,
};

enum
{
  NC_GALAXY_SD_OBS_REDSHIFT_GAUSS_LSTATE_PZP = 0,
  NC_GALAXY_SD_OBS_REDSHIFT_GAUSS_LSTATE_PZ_GIVEN_ZP,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcGalaxySDObsRedshiftGauss, nc_galaxy_sd_obs_redshift_gauss, NC_TYPE_GALAXY_SD_OBS_REDSHIFT);

static void
nc_galaxy_sd_obs_redshift_gauss_init (NcGalaxySDObsRedshiftGauss *gsdorgauss)
{
  NcGalaxySDObsRedshiftGaussPrivate * const self = nc_galaxy_sd_obs_redshift_gauss_get_instance_private (gsdorgauss);

  self->sdz                = NULL;
  self->bin_sigma0         = 0.0;
  self->reltol             = 0.0;
  self->zp_min             = 0.0;
  self->zp_max             = 0.0;
  self->zp_support_max     = 0.0;
  self->pzp_spline         = NULL;
  self->stats              = NULL;
  self->pz_given_zp_spline = NULL;
}

static void
_nc_galaxy_sd_obs_redshift_gauss_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcGalaxySDObsRedshiftGauss *gsdorgauss = NC_GALAXY_SD_OBS_REDSHIFT_GAUSS (object);

  g_return_if_fail (NC_IS_GALAXY_SD_OBS_REDSHIFT_GAUSS (gsdorgauss));

  switch (prop_id)
  {
    case PROP_ZP_LIM:
    {
      NcmDTuple2 *lim = g_value_get_boxed (value);

      if (lim == NULL)
        g_error ("_nc_galaxy_sd_obs_redshift_gauss_set_property: lim is NULL");

      nc_galaxy_sd_obs_redshift_gauss_set_zp_lim (gsdorgauss, lim->elements[0], lim->elements[1]);
      break;
    }
    case PROP_USE_TRUE_Z:
      nc_galaxy_sd_obs_redshift_gauss_set_use_true_z (gsdorgauss, g_value_get_boolean (value));
      break;
    case PROP_BIN_SIGMA0:
      nc_galaxy_sd_obs_redshift_gauss_set_bin_sigma0 (gsdorgauss, g_value_get_double (value));
      break;
    case PROP_RELTOL:
      nc_galaxy_sd_obs_redshift_gauss_set_reltol (gsdorgauss, g_value_get_double (value));
      break;
    case PROP_ZP_SUPPORT_MAX:
      nc_galaxy_sd_obs_redshift_gauss_set_zp_support_max (gsdorgauss, g_value_get_double (value));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_galaxy_sd_obs_redshift_gauss_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcGalaxySDObsRedshiftGauss *gsdorgauss         = NC_GALAXY_SD_OBS_REDSHIFT_GAUSS (object);
  NcGalaxySDObsRedshiftGaussPrivate * const self = nc_galaxy_sd_obs_redshift_gauss_get_instance_private (gsdorgauss);

  g_return_if_fail (NC_IS_GALAXY_SD_OBS_REDSHIFT_GAUSS (gsdorgauss));

  switch (prop_id)
  {
    case PROP_ZP_LIM:
    {
      gdouble zp_min = self->zp_min;
      gdouble zp_max = self->zp_max;

      g_value_take_boxed (value, ncm_dtuple2_new (zp_min, zp_max));
      break;
    }
    case PROP_USE_TRUE_Z:
      g_value_set_boolean (value, nc_galaxy_sd_obs_redshift_gauss_get_use_true_z (gsdorgauss));
      break;
    case PROP_BIN_SIGMA0:
      g_value_set_double (value, nc_galaxy_sd_obs_redshift_gauss_get_bin_sigma0 (gsdorgauss));
      break;
    case PROP_RELTOL:
      g_value_set_double (value, nc_galaxy_sd_obs_redshift_gauss_get_reltol (gsdorgauss));
      break;
    case PROP_ZP_SUPPORT_MAX:
      g_value_set_double (value, nc_galaxy_sd_obs_redshift_gauss_get_zp_support_max (gsdorgauss));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_galaxy_sd_obs_redshift_gauss_dispose (GObject *object)
{
  NcGalaxySDObsRedshiftGauss *gsdorgauss         = NC_GALAXY_SD_OBS_REDSHIFT_GAUSS (object);
  NcGalaxySDObsRedshiftGaussPrivate * const self = nc_galaxy_sd_obs_redshift_gauss_get_instance_private (gsdorgauss);

  ncm_spline_clear (&self->pzp_spline);
  ncm_spline_clear (&self->pz_given_zp_spline);
  ncm_stats_dist1d_clear (&self->stats);

  /* Chain up: end */
  G_OBJECT_CLASS (nc_galaxy_sd_obs_redshift_gauss_parent_class)->dispose (object);
}

static void
_nc_galaxy_sd_obs_redshift_gauss_finalize (GObject *object)
{
  /* Chain up: end */
  G_OBJECT_CLASS (nc_galaxy_sd_obs_redshift_gauss_parent_class)->finalize (object);
}

static void _nc_galaxy_sd_obs_redshift_gauss_gen (NcGalaxySDObsRedshift *gsdor, NcGalaxySDObsRedshiftData *data, NcmRNG *rng);
static gboolean _nc_galaxy_sd_obs_redshift_gauss_gen1 (NcGalaxySDObsRedshift *gsdor, NcGalaxySDObsRedshiftData *data, NcmRNG *rng);
static void _nc_galaxy_sd_obs_redshift_gauss_prepare (NcGalaxySDObsRedshift *gsdor, NcGalaxySDObsRedshiftData *data);
static void _nc_galaxy_sd_obs_redshift_gauss_get_integ_lim (NcGalaxySDObsRedshift *gsdor, NcGalaxySDObsRedshiftData *data, gdouble *z_min, gdouble *z_max);
static NcGalaxySDObsRedshiftIntegrand *_nc_galaxy_sd_obs_redshift_gauss_integ (NcGalaxySDObsRedshift *gsdor, gboolean use_lnp);
static void _nc_galaxy_sd_obs_redshift_gauss_data_init (NcGalaxySDObsRedshift *gsdor, NcGalaxySDObsRedshiftData *data);
static NcmSpline *_nc_galaxy_sd_obs_redshift_gauss_compute_binned_dndz (NcGalaxySDObsRedshift *gsdor, NcmVector *z_array);
static void _nc_galaxy_sd_obs_redshift_gauss_add_submodel (NcmModel *model, NcmModel *submodel);
static NcmIntegralFixed *_nc_galaxy_sd_obs_redshift_gauss_prepare_fixed_nodes (NcGalaxySDObsRedshift *gsdor, NcmMSet *mset, NcGalaxySDObsRedshiftData *data, guint n_nodes, guint rule_n);

static void
nc_galaxy_sd_obs_redshift_gauss_class_init (NcGalaxySDObsRedshiftGaussClass *klass)
{
  NcGalaxySDObsRedshiftClass *gsdor_class = NC_GALAXY_SD_OBS_REDSHIFT_CLASS (klass);
  GObjectClass *object_class              = G_OBJECT_CLASS (klass);
  NcmModelClass *model_class              = NCM_MODEL_CLASS (klass);

  model_class->set_property = &_nc_galaxy_sd_obs_redshift_gauss_set_property;
  model_class->get_property = &_nc_galaxy_sd_obs_redshift_gauss_get_property;
  object_class->dispose     = &_nc_galaxy_sd_obs_redshift_gauss_dispose;
  object_class->finalize    = &_nc_galaxy_sd_obs_redshift_gauss_finalize;

  ncm_model_class_set_name_nick (model_class, "Gaussian Observed Redshift", "GalaxySDObsRedshiftGauss");
  ncm_model_class_add_params (model_class, 0, 0, PROP_LEN);

  /**
   * NcGalaxySDObsRedshiftGauss:zp_lim:
   *
   * Galaxy sample photometric redshift limits.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_ZP_LIM,
                                   g_param_spec_boxed ("zp-lim",
                                                       NULL,
                                                       "Galaxy sample photometric redshift limits",
                                                       NCM_TYPE_DTUPLE2,
                                                       G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcGalaxySDObsRedshiftGauss:use_true_z:
   *
   * Determines how the redshift distribution is modeled.
   *
   * If set to true, both the true redshift distribution and the observed redshift distribution are used.
   * In this case, the standard deviation is set to sigma0 * (1 + z), where z is the true redshift.
   *
   * If set to false, a simplified Gaussian redshift distribution is assumed, centered at zp with
   * a fixed standard deviation sigmaz.
   */
  g_object_class_install_property (object_class,
                                   PROP_USE_TRUE_Z,
                                   g_param_spec_boolean ("use-true-z",
                                                         NULL,
                                                         "Use the true redshift distribution",
                                                         TRUE,
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcGalaxySDObsRedshiftGauss:bin_sigma0:
   *
   * Base photometric redshift scatter parameter for binned analyses.
   *
   * This property stores the characteristic $\sigma_0$ value used in computing binned
   * redshift distributions. It allows methods like compute_binned_dndz to use the
   * stored value instead of requiring it as a parameter.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_BIN_SIGMA0,
                                   g_param_spec_double ("bin-sigma0",
                                                        NULL,
                                                        "Base photometric redshift scatter for binned analyses",
                                                        1.0e-6, G_MAXDOUBLE, 0.03,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcGalaxySDObsRedshiftGauss:reltol:
   *
   * Relative tolerance for numerical integration.
   *
   * This property sets the default relative error tolerance used in numerical
   * integrations when computing photometric redshift distributions.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_RELTOL,
                                   g_param_spec_double ("reltol",
                                                        NULL,
                                                        "Relative tolerance for numerical integration",
                                                        DBL_EPSILON, 1.0e-2, 1.0e-7,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcGalaxySDObsRedshiftGauss:zp-support-max:
   *
   * Maximum photometric redshift for support.
   *
   * This property defines the maximum photometric redshift value to consider
   * when computing distributions and bin edges.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_ZP_SUPPORT_MAX,
                                   g_param_spec_double ("zp-support-max",
                                                        NULL,
                                                        "Maximum photometric redshift for support",
                                                        0.1, 100.0, 20.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  ncm_model_class_check_params_info (model_class);

  gsdor_class->gen                 = &_nc_galaxy_sd_obs_redshift_gauss_gen;
  gsdor_class->gen1                = &_nc_galaxy_sd_obs_redshift_gauss_gen1;
  gsdor_class->prepare             = &_nc_galaxy_sd_obs_redshift_gauss_prepare;
  gsdor_class->get_integ_lim       = &_nc_galaxy_sd_obs_redshift_gauss_get_integ_lim;
  gsdor_class->integ               = &_nc_galaxy_sd_obs_redshift_gauss_integ;
  gsdor_class->data_init           = &_nc_galaxy_sd_obs_redshift_gauss_data_init;
  gsdor_class->compute_binned_dndz = &_nc_galaxy_sd_obs_redshift_gauss_compute_binned_dndz;
  gsdor_class->prepare_fixed_nodes = &_nc_galaxy_sd_obs_redshift_gauss_prepare_fixed_nodes;
  model_class->add_submodel        = &_nc_galaxy_sd_obs_redshift_gauss_add_submodel;
}

#define VECTOR (NCM_MODEL (gsdor))
#define SIGMA  (ncm_model_orig_param_get (VECTOR, NC_GALAXY_SD_OBS_REDSHIFT_GAUSS_SIGMA))

static void
_nc_galaxy_sd_obs_redshift_gauss_gen (NcGalaxySDObsRedshift *gsdor, NcGalaxySDObsRedshiftData *data, NcmRNG *rng)
{
  NcGalaxySDObsRedshiftGauss *gsdorgauss         = NC_GALAXY_SD_OBS_REDSHIFT_GAUSS (gsdor);
  NcGalaxySDObsRedshiftGaussPrivate * const self = nc_galaxy_sd_obs_redshift_gauss_get_instance_private (gsdorgauss);
  NcGalaxySDObsRedshiftGaussData * const ldata   = (NcGalaxySDObsRedshiftGaussData *) data->ldata;
  guint max_iter                                 = 1000;
  gdouble sigmaz;
  gdouble zp;
  gdouble z;

  do {
    z      = nc_galaxy_sd_true_redshift_gen (self->sdz, rng);
    sigmaz = ldata->sigma0 * (1.0 + z);
    zp     = ncm_rng_gaussian_gen (rng, z, sigmaz);

    if (max_iter-- == 0)
      g_error ("nc_galaxy_sd_obs_redshift_gauss_gen: maximum number of iterations reached.");
  } while ((zp > self->zp_max) || (zp < self->zp_min));

  data->z      = z;
  ldata->zp    = zp;
  ldata->sigma = sigmaz;
}

static gboolean
_nc_galaxy_sd_obs_redshift_gauss_gen1 (NcGalaxySDObsRedshift *gsdor, NcGalaxySDObsRedshiftData *data, NcmRNG *rng)
{
  NcGalaxySDObsRedshiftGauss *gsdorgauss         = NC_GALAXY_SD_OBS_REDSHIFT_GAUSS (gsdor);
  NcGalaxySDObsRedshiftGaussPrivate * const self = nc_galaxy_sd_obs_redshift_gauss_get_instance_private (gsdorgauss);
  NcGalaxySDObsRedshiftGaussData * const ldata   = (NcGalaxySDObsRedshiftGaussData *) data->ldata;
  gdouble sigmaz;
  gdouble zp;
  gdouble z;

  z      = nc_galaxy_sd_true_redshift_gen (self->sdz, rng);
  sigmaz = ldata->sigma0 * (1.0 + z);
  zp     = ncm_rng_gaussian_gen (rng, z, sigmaz);

  data->z      = z;
  ldata->zp    = zp;
  ldata->sigma = sigmaz;

  return (zp >= self->zp_min) && (zp <= self->zp_max);
}

static void
_nc_galaxy_sd_obs_redshift_gauss_prepare (NcGalaxySDObsRedshift *gsdor, NcGalaxySDObsRedshiftData *data)
{
  /* Nothing to do */
}

static void
_nc_galaxy_sd_obs_redshift_gauss_get_integ_lim (NcGalaxySDObsRedshift *gsdor, NcGalaxySDObsRedshiftData *data, gdouble *z_min, gdouble *z_max)
{
  NcGalaxySDObsRedshiftGauss *gsdorgauss         = NC_GALAXY_SD_OBS_REDSHIFT_GAUSS (gsdor);
  NcGalaxySDObsRedshiftGaussPrivate * const self = nc_galaxy_sd_obs_redshift_gauss_get_instance_private (gsdorgauss);

  nc_galaxy_sd_true_redshift_get_lim (self->sdz, z_min, z_max);
}

struct _IntegData
{
  NcGalaxySDObsRedshiftGauss *gsdorgauss;
};

/* LCOV_EXCL_START */
static gpointer
_integ_data_copy (gpointer idata)
{
  struct _IntegData *new_idata = g_new0 (struct _IntegData, 1);

  *new_idata = *(struct _IntegData *) idata;

  return new_idata;
}

/* LCOV_EXCL_STOP */

static void
_integ_data_free (gpointer idata)
{
  g_free (idata);
}

static gdouble
_nc_galaxy_sd_obs_redshift_gauss_ln_integ_f (gpointer callback_data, const gdouble z, NcGalaxySDObsRedshiftData *data)
{
  const struct _IntegData *int_data              = (struct _IntegData *) callback_data;
  NcGalaxySDObsRedshiftGaussPrivate * const self = nc_galaxy_sd_obs_redshift_gauss_get_instance_private (int_data->gsdorgauss);
  NcGalaxySDObsRedshiftGaussData * const ldata   = (NcGalaxySDObsRedshiftGaussData *) data->ldata;
  const gdouble zp                               = ldata->zp;
  gdouble sign;

  if (self->use_true_z)
  {
    const gdouble sigmaz    = ldata->sigma0 * (1.0 + z);
    const gdouble norm      = sqrt (2.0 * M_PI) * sigmaz;
    const gdouble lognorm   = ncm_util_log_gaussian_integral (self->zp_min, self->zp_max, z, sigmaz, &sign);
    const gdouble ln_int_z  = nc_galaxy_sd_true_redshift_ln_integ (self->sdz, z);
    const gdouble ln_int_zp = -0.5 * gsl_pow_2 ((zp - z) / sigmaz) - lognorm - log (norm);

    return ln_int_z + ln_int_zp;
  }
  else
  {
    const gdouble norm      = sqrt (2.0 * M_PI) * ldata->sigma;
    const gdouble ln_int_zp = -0.5 * gsl_pow_2 ((zp - z) / ldata->sigma) - log (norm);

    return ln_int_zp;
  }
}

static gdouble
_nc_galaxy_sd_obs_redshift_gauss_integ_f (gpointer callback_data, const gdouble z, NcGalaxySDObsRedshiftData *data)
{
  const struct _IntegData *int_data              = (struct _IntegData *) callback_data;
  NcGalaxySDObsRedshiftGaussPrivate * const self = nc_galaxy_sd_obs_redshift_gauss_get_instance_private (int_data->gsdorgauss);
  NcGalaxySDObsRedshiftGaussData * const ldata   = (NcGalaxySDObsRedshiftGaussData *) data->ldata;
  const gdouble zp                               = ldata->zp;
  gdouble sign;

  if (self->use_true_z)
  {
    const gdouble sigmaz  = ldata->sigma0 * (1.0 + z);
    const gdouble norm    = sqrt (2.0 * M_PI) * sigmaz;
    const gdouble lognorm = ncm_util_log_gaussian_integral (self->zp_min, self->zp_max, z, sigmaz, &sign);
    const gdouble int_z   = nc_galaxy_sd_true_redshift_integ (self->sdz, z);
    const gdouble int_zp  = exp (-0.5 * gsl_pow_2 ((zp - z) / sigmaz) - lognorm) / norm;

    return int_z * int_zp;
  }
  else
  {
    const gdouble norm   = sqrt (2.0 * M_PI) * ldata->sigma;
    const gdouble int_zp = exp (-0.5 * gsl_pow_2 ((zp - z) / ldata->sigma)) / norm;

    return int_zp;
  }
}

static NcGalaxySDObsRedshiftIntegrand *
_nc_galaxy_sd_obs_redshift_gauss_integ (NcGalaxySDObsRedshift *gsdor, gboolean use_lnp)
{
  NcGalaxySDObsRedshiftGauss *gsdorgauss = NC_GALAXY_SD_OBS_REDSHIFT_GAUSS (gsdor);
  struct _IntegData *int_data            = g_new0 (struct _IntegData, 1);
  NcGalaxySDObsRedshiftIntegrand *integ  = nc_galaxy_sd_obs_redshift_integrand_new (use_lnp ? _nc_galaxy_sd_obs_redshift_gauss_ln_integ_f : _nc_galaxy_sd_obs_redshift_gauss_integ_f,
                                                                                    _integ_data_free,
                                                                                    _integ_data_copy,
                                                                                    NULL,
                                                                                    int_data);

  int_data->gsdorgauss = gsdorgauss;

  return integ;
}

static void
_nc_galaxy_sd_obs_redshift_gauss_ldata_free (gpointer ldata)
{
  g_free (ldata);
}

static void
_nc_galaxy_sd_obs_redshift_gauss_ldata_read_row (NcGalaxySDObsRedshiftData *data, NcGalaxyWLObs *obs, const guint i)
{
  NcGalaxySDObsRedshiftGaussData *ldata = (NcGalaxySDObsRedshiftGaussData *) data->ldata;

  ldata->zp     = nc_galaxy_wl_obs_get (obs, NC_GALAXY_SD_OBS_REDSHIFT_GAUSS_COL_ZP, i, NULL);
  ldata->sigma0 = nc_galaxy_wl_obs_get (obs, NC_GALAXY_SD_OBS_REDSHIFT_GAUSS_COL_SIGMA0, i, NULL);
  ldata->sigma  = nc_galaxy_wl_obs_get (obs, NC_GALAXY_SD_OBS_REDSHIFT_GAUSS_COL_SIGMA, i, NULL);
}

static void
_nc_galaxy_sd_obs_redshift_gauss_ldata_write_row (NcGalaxySDObsRedshiftData *data, NcGalaxyWLObs *obs, const guint i)
{
  NcGalaxySDObsRedshiftGaussData *ldata = (NcGalaxySDObsRedshiftGaussData *) data->ldata;

  nc_galaxy_wl_obs_set (obs, NC_GALAXY_SD_OBS_REDSHIFT_GAUSS_COL_ZP, i, ldata->zp, NULL);
  nc_galaxy_wl_obs_set (obs, NC_GALAXY_SD_OBS_REDSHIFT_GAUSS_COL_SIGMA0, i, ldata->sigma0, NULL);
  nc_galaxy_wl_obs_set (obs, NC_GALAXY_SD_OBS_REDSHIFT_GAUSS_COL_SIGMA, i, ldata->sigma, NULL);
}

static void
_nc_galaxy_sd_obs_redshift_gauss_ldata_required_columns (NcGalaxySDObsRedshiftData *data, GList *columns)
{
  columns = g_list_append (columns, g_strdup (NC_GALAXY_SD_OBS_REDSHIFT_GAUSS_COL_ZP));
  columns = g_list_append (columns, g_strdup (NC_GALAXY_SD_OBS_REDSHIFT_GAUSS_COL_SIGMA0));
  columns = g_list_append (columns, g_strdup (NC_GALAXY_SD_OBS_REDSHIFT_GAUSS_COL_SIGMA));
}

static void
_nc_galaxy_sd_obs_redshift_gauss_data_init (NcGalaxySDObsRedshift *gsdor, NcGalaxySDObsRedshiftData *data)
{
  NcGalaxySDObsRedshiftGaussData *ldata = g_new0 (NcGalaxySDObsRedshiftGaussData, 1);

  data->ldata                  = ldata;
  data->ldata_destroy          = &_nc_galaxy_sd_obs_redshift_gauss_ldata_free;
  data->ldata_read_row         = &_nc_galaxy_sd_obs_redshift_gauss_ldata_read_row;
  data->ldata_write_row        = &_nc_galaxy_sd_obs_redshift_gauss_ldata_write_row;
  data->ldata_required_columns = &_nc_galaxy_sd_obs_redshift_gauss_ldata_required_columns;
}

typedef struct _GaussFixedNodesGSLArg
{
  struct _IntegData int_data;
  NcGalaxySDObsRedshiftData *data;
} GaussFixedNodesGSLArg;

static gdouble
_gauss_fixed_nodes_gsl_f (gdouble z, gpointer user_data)
{
  GaussFixedNodesGSLArg *arg = (GaussFixedNodesGSLArg *) user_data;

  return _nc_galaxy_sd_obs_redshift_gauss_integ_f (&arg->int_data, z, arg->data);
}

static NcmIntegralFixed *
_nc_galaxy_sd_obs_redshift_gauss_prepare_fixed_nodes (
  NcGalaxySDObsRedshift     *gsdor,
  NcmMSet                   *mset,
  NcGalaxySDObsRedshiftData *data,
  guint                     n_nodes,
  guint                     rule_n
)
{
  NcGalaxySDObsRedshiftGauss *gsdorgauss         = NC_GALAXY_SD_OBS_REDSHIFT_GAUSS (gsdor);
  NcGalaxySDObsRedshiftGaussPrivate * const self = nc_galaxy_sd_obs_redshift_gauss_get_instance_private (gsdorgauss);
  NcGalaxySDObsRedshiftGaussData * const ldata   = (NcGalaxySDObsRedshiftGaussData *) data->ldata;
  GaussFixedNodesGSLArg arg;
  NcmIntegralFixed *intf;
  gsl_function F;
  gdouble z_min, z_max;

  g_assert_nonnull (self->sdz);
  nc_galaxy_sd_true_redshift_get_lim (self->sdz, &z_min, &z_max);

  /* Restricting the nodes to the Gaussian's effective support */
  {
    const gdouble sigma_max = self->use_true_z
      ? ldata->sigma0 * (
      (1.0 + 7.0 * ldata->sigma0) * (1.0 + ldata->zp)
      )
      : ldata->sigma;
    const gdouble half_width = 7.0 * sigma_max;

    z_min = MAX (z_min, ldata->zp - half_width);
    z_max = MIN (z_max, ldata->zp + half_width);
  }

  intf = ncm_integral_fixed_new (n_nodes, rule_n, z_min, z_max);

  arg.int_data.gsdorgauss = gsdorgauss;
  arg.data                = data;

  F.function = &_gauss_fixed_nodes_gsl_f;
  F.params   = &arg;

  ncm_integral_fixed_calc_nodes (intf, &F);

  return intf;
}

typedef struct _BinnedDndzIntegData
{
  NcGalaxySDTrueRedshift *gsdtr;
  gdouble zpl;
  gdouble zpu;
  gdouble sigma0;
} BinnedDndzIntegData;

static gdouble
_binned_dndz_integrand (gpointer user_data, const gdouble z, const gdouble w)
{
  BinnedDndzIntegData *data  = (BinnedDndzIntegData *) user_data;
  const gdouble Pz           = nc_galaxy_sd_true_redshift_integ (data->gsdtr, z);
  const gdouble sigmaz       = data->sigma0 * (1.0 + z);
  const gdouble sqrt2_sigmaz = M_SQRT2 * sigmaz;
  const gdouble W            = ncm_util_gaussian_integral (data->zpl, data->zpu, z, sigmaz);
  const gdouble N            = 0.5 * erfc (-z / sqrt2_sigmaz);

  return Pz * W / N;
}

static gdouble
_binned_dndz_integrand_for_gsl (const gdouble z, gpointer user_data)
{
  return _binned_dndz_integrand (user_data, z, 1.0);
}

/**
 * _find_effective_support:
 * @F: GSL function to evaluate
 * @z_min: minimum z value
 * @z_max: maximum z value
 * @threshold: minimum relative value to consider (e.g., 1e-20)
 * @z_min_eff: (out): effective minimum z
 * @z_max_eff: (out): effective maximum z
 *
 * Finds the effective support of a function by locating where it drops
 * below a threshold relative to its maximum value.
 */
static void
_find_effective_support (gsl_function *F, gdouble z_min, gdouble z_max, gdouble threshold, gdouble *z_min_eff, gdouble *z_max_eff)
{
  const guint n_search = 500;
  gdouble z_array[n_search];
  gdouble f_array[n_search];
  gdouble f_max = -G_MAXDOUBLE;
  guint i_max   = 0;
  guint i;

  /* Sample the function to find maximum */
  for (i = 0; i < n_search; i++)
  {
    z_array[i] = z_min + (z_max - z_min) * i / (n_search - 1.0);
    f_array[i] = GSL_FN_EVAL (F, z_array[i]);

    if (f_array[i] > f_max)
    {
      f_max = f_array[i];
      i_max = i;
    }
  }

  /* Find left tail: search backwards from maximum */
  *z_min_eff = z_min;

  for (i = i_max; i > 0; i--)
  {
    if (f_array[i] < threshold * f_max)
    {
      *z_min_eff = z_array[i];
      break;
    }
  }

  /* Find right tail: search forwards from maximum */
  *z_max_eff = z_max;

  for (i = i_max; i < n_search; i++)
  {
    if (f_array[i] < threshold * f_max)
    {
      *z_max_eff = z_array[i];
      break;
    }
  }

  /* Ensure we have a valid range */
  if (*z_min_eff >= *z_max_eff)
  {
    *z_min_eff = z_min;
    *z_max_eff = z_max;
  }
}

typedef struct _PhotozDistIntegData
{
  NcGalaxySDTrueRedshift *gsdtr;
  NcmIntegral1dPtr *integrator;
  gdouble zp;
  gdouble sigma0;
  gdouble z_min;
  gdouble z_max;
} PhotozDistIntegData;

static gdouble
_photoz_distribution_integrand (gpointer user_data, const gdouble z, const gdouble w)
{
  PhotozDistIntegData *data    = (PhotozDistIntegData *) user_data;
  const gdouble Pz             = nc_galaxy_sd_true_redshift_integ (data->gsdtr, z);
  const gdouble sigmaz         = data->sigma0 * (1.0 + z);
  const gdouble sqrt2pi_sigmaz = M_SQRT2 * M_SQRTPI * sigmaz;
  const gdouble gauss          = exp (-0.5 * gsl_pow_2 ((data->zp - z) / sigmaz));
  const gdouble norm           = 0.5 * erfc (-z / (M_SQRT2 * sigmaz));

  return Pz * gauss / (sqrt2pi_sigmaz * norm);
}

static gdouble
_photoz_ln_distribution_gsl (gdouble zp, gpointer user_data)
{
  PhotozDistIntegData *data = (PhotozDistIntegData *) user_data;
  gdouble error             = 0.0;
  gdouble result;

  data->zp = zp;
  ncm_integral1d_ptr_set_userdata (data->integrator, data);
  result = ncm_integral1d_eval (
    NCM_INTEGRAL1D (data->integrator),
    data->z_min,
    data->z_max,
    &error);

  g_assert_cmpfloat (result, >, 0.0);

  return -2.0 * log (fabs (result));
}

static void
_nc_galaxy_sd_obs_redshift_gauss_add_submodel (NcmModel *model, NcmModel *submodel)
{
  /* Chain up: start */
  NCM_MODEL_CLASS (nc_galaxy_sd_obs_redshift_gauss_parent_class)->add_submodel (model, submodel);
  {
    NcGalaxySDObsRedshiftGauss *gsdorgauss         = NC_GALAXY_SD_OBS_REDSHIFT_GAUSS (model);
    NcGalaxySDObsRedshiftGaussPrivate * const self = nc_galaxy_sd_obs_redshift_gauss_get_instance_private (gsdorgauss);

    g_assert (ncm_model_is_submodel (submodel));
    g_assert (NC_IS_GALAXY_SD_TRUE_REDSHIFT (submodel));

    self->sdz = NC_GALAXY_SD_TRUE_REDSHIFT (submodel);
  }
}

/**
 * _nc_galaxy_sd_obs_redshift_gauss_prepare_pzp:
 * @gsdorgauss: a #NcGalaxySDObsRedshiftGauss
 *
 * Prepares and caches the marginal photometric redshift distribution.
 *
 * This distribution is computed as:
 *
 * $$P(z_p) = \int_0^\infty P(z) \, \mathrm{Gauss}(z_p|z,\sigma_z) \, \mathrm{d}z$$
 *
 * where $\sigma_z = \sigma_0 (1 + z)$ and the normalization ensures $z_p > 0$.
 *
 * The spline is stored as $-2\ln P(z_p)$ to enable reuse by equal-area binning.
 * The result is cached and will be reused if called again with the same
 * parameters, unless the model state has been invalidated.
 */
static void
_nc_galaxy_sd_obs_redshift_gauss_prepare_pzp (NcGalaxySDObsRedshiftGauss *gsdorgauss)
{
  NcGalaxySDObsRedshiftGaussPrivate * const self = nc_galaxy_sd_obs_redshift_gauss_get_instance_private (gsdorgauss);

  if ((self->pzp_spline == NULL) ||
      !ncm_model_lstate_is_update (NCM_MODEL (gsdorgauss), NC_GALAXY_SD_OBS_REDSHIFT_GAUSS_LSTATE_PZP))
  {
    gdouble z_min, z_max;
    PhotozDistIntegData integ_data;
    gsl_function F;

    nc_galaxy_sd_true_redshift_get_lim (self->sdz, &z_min, &z_max);

    integ_data.gsdtr      = self->sdz;
    integ_data.integrator = ncm_integral1d_ptr_new (&_photoz_distribution_integrand, NULL);
    integ_data.zp         = 0.0;
    integ_data.sigma0     = self->bin_sigma0;
    integ_data.z_min      = z_min;
    integ_data.z_max      = z_max;

    F.function = &_photoz_ln_distribution_gsl;
    F.params   = &integ_data;

    /* Clear old spline if it exists */
    ncm_spline_clear (&self->pzp_spline);

    self->pzp_spline = NCM_SPLINE (ncm_spline_cubic_notaknot_new ());
    ncm_spline_set_func (
      self->pzp_spline,
      NCM_SPLINE_FUNCTION_SPLINE,
      &F, 0.0,
      self->zp_support_max,
      0,
      self->reltol);

    ncm_integral1d_ptr_free (integ_data.integrator);

    ncm_stats_dist1d_clear (&self->stats);
    self->stats = NCM_STATS_DIST1D (ncm_stats_dist1d_spline_new (self->pzp_spline));
    g_object_set (self->stats, "abstol", 1.0e-50, NULL);
    ncm_stats_dist1d_prepare (NCM_STATS_DIST1D (self->stats));

    ncm_model_lstate_set_update (NCM_MODEL (gsdorgauss), NC_GALAXY_SD_OBS_REDSHIFT_GAUSS_LSTATE_PZP);
  }
}

/**
 * _nc_galaxy_sd_obs_redshift_gauss_prepare_pz_given_zp:
 * @gsdorgauss: a #NcGalaxySDObsRedshiftGauss
 *
 * Prepares and caches the conditional true redshift distribution.
 *
 * This distribution gives the probability of true redshift z for galaxies observed
 * in the photometric redshift bin [zp_min, zp_max]:
 *
 * $$P(z|z_{p,\mathrm{min}}, z_{p,\mathrm{max}}) = \frac{P(z) \, W(z_{p,\mathrm{min}}, z_{p,\mathrm{max}}|z)}{N}$$
 *
 * where $W$ is the Gaussian integral over the bin and $N$ is a normalization factor.
 *
 * The result is cached as a spline and will be reused if called again with the same
 * parameters, unless the model state has been invalidated.
 */
static void
_nc_galaxy_sd_obs_redshift_gauss_prepare_pz_given_zp (NcGalaxySDObsRedshiftGauss *gsdorgauss)
{
  NcGalaxySDObsRedshiftGaussPrivate * const self = nc_galaxy_sd_obs_redshift_gauss_get_instance_private (gsdorgauss);

  if ((self->pz_given_zp_spline == NULL) ||
      !ncm_model_lstate_is_update (NCM_MODEL (gsdorgauss), NC_GALAXY_SD_OBS_REDSHIFT_GAUSS_LSTATE_PZ_GIVEN_ZP))
  {
    gdouble z_min, z_max;
    gdouble z_min_eff, z_max_eff;
    BinnedDndzIntegData integ_data;
    gsl_function F;

    nc_galaxy_sd_true_redshift_get_lim (self->sdz, &z_min, &z_max);

    integ_data.gsdtr  = self->sdz;
    integ_data.zpl    = self->zp_min;
    integ_data.zpu    = self->zp_max;
    integ_data.sigma0 = self->bin_sigma0;

    F.function = &_binned_dndz_integrand_for_gsl;
    F.params   = &integ_data;

    /* Find effective support where function is above 1e-20 of its maximum */
    _find_effective_support (&F, z_min, z_max, 1.0e-20, &z_min_eff, &z_max_eff);

    ncm_spline_clear (&self->pz_given_zp_spline);
    self->pz_given_zp_spline = NCM_SPLINE (ncm_spline_cubic_notaknot_new ());
    ncm_spline_set_func (self->pz_given_zp_spline, NCM_SPLINE_FUNCTION_SPLINE, &F, z_min_eff, z_max_eff, 0, self->reltol);

    {
      NcmIntegral1dPtr *integrator = ncm_integral1d_ptr_new (&_binned_dndz_integrand, NULL);
      gdouble norm_error           = 0.0;
      gdouble norm;

      ncm_integral1d_set_reltol (NCM_INTEGRAL1D (integrator), self->reltol);
      ncm_integral1d_ptr_set_userdata (integrator, &integ_data);

      norm = ncm_integral1d_eval (NCM_INTEGRAL1D (integrator), z_min, z_max, &norm_error);

      g_assert_cmpfloat (norm, >, 0.0);

      ncm_vector_scale (ncm_spline_peek_yv (self->pz_given_zp_spline), 1.0 / norm);
      ncm_spline_prepare (self->pz_given_zp_spline);

      ncm_integral1d_ptr_free (integrator);
    }

    /* Mark model local state as updated */
    ncm_model_lstate_set_update (NCM_MODEL (gsdorgauss), NC_GALAXY_SD_OBS_REDSHIFT_GAUSS_LSTATE_PZ_GIVEN_ZP);
  }
}

/**
 * nc_galaxy_sd_obs_redshift_gauss_new:
 * @sdz: a #NcGalaxySDTrueRedshift
 * @zp_min: the minimum photometric redshift
 * @zp_max: the maximum photometric redshift
 *
 * Creates a new #NcGalaxySDObsRedshiftGauss object.
 *
 * Returns: (transfer full): a new #NcGalaxySDObsRedshiftGauss object.
 */
NcGalaxySDObsRedshiftGauss *
nc_galaxy_sd_obs_redshift_gauss_new (NcGalaxySDTrueRedshift *sdz, const gdouble zp_min, const gdouble zp_max)
{
  NcmDTuple2 lim                         = NCM_DTUPLE2_STATIC_INIT (zp_min, zp_max);
  NcGalaxySDObsRedshiftGauss *gsdorgauss = g_object_new (NC_TYPE_GALAXY_SD_OBS_REDSHIFT_GAUSS,
                                                         "zp-lim", &lim,
                                                         NULL);

  ncm_model_add_submodel (NCM_MODEL (gsdorgauss), NCM_MODEL (sdz));

  return gsdorgauss;
}

/**
 * nc_galaxy_sd_obs_redshift_gauss_new_lsst_srd_bins:
 * @type: a #NcGalaxySDTrueRedshiftLSSTSRDType
 * @gsdtr_new: (out) (transfer full) (optional) (nullable): the true redshift distribution submodel
 *
 * Creates LSST SRD photometric redshift bins.
 *
 * Creates a #GPtrArray of #NcGalaxySDObsRedshiftGauss objects, one for each LSST SRD
 * bin of the specified type. All objects share the same true redshift distribution
 * submodel. The bin edges are computed based on the type:
 *
 * - For lens types (Y1 and Y10): Uses linearly spaced bins in photo-z
 *   - Y1: 5 bins from 0.2 to 1.2 with sigma_z = 0.03
 *   - Y10: 10 bins from 0.2 to 1.2 with sigma_z = 0.03
 * - For source types (Y1 and Y10): Uses equal-area bins up to photo-z = 3.5
 *   - Both Y1 and Y10: 5 bins with sigma_z = 0.05
 *
 * Returns: (transfer container) (element-type NcGalaxySDObsRedshiftGauss): a #GPtrArray
 * containing all bin objects.
 */
GPtrArray *
nc_galaxy_sd_obs_redshift_gauss_new_lsst_srd_bins (NcGalaxySDTrueRedshiftLSSTSRDType type, NcGalaxySDTrueRedshiftLSSTSRD **gsdtr_new)
{
  NcGalaxySDTrueRedshiftLSSTSRD *gsdtr = nc_galaxy_sd_true_redshift_lsst_srd_new_from_type (type);
  const gdouble zp_support_max         = 20.0;
  NcmVector *bin_edges                 = NULL;
  GPtrArray *bins;
  gdouble sigma_z;
  guint n_bins;

  /* Determine number of bins and sigma_z based on type */
  switch (type)
  {
    case NC_GALAXY_SD_TRUE_REDSHIFT_LSST_SRD_Y1_SOURCE:
      n_bins  = 5;
      sigma_z = 0.05;
      break;
    case NC_GALAXY_SD_TRUE_REDSHIFT_LSST_SRD_Y1_LENS:
      n_bins  = 5;
      sigma_z = 0.03;
      break;
    case NC_GALAXY_SD_TRUE_REDSHIFT_LSST_SRD_Y10_SOURCE:
      n_bins  = 5;
      sigma_z = 0.05;
      break;
    case NC_GALAXY_SD_TRUE_REDSHIFT_LSST_SRD_Y10_LENS:
      n_bins  = 10;
      sigma_z = 0.03;
      break;
    default:
      g_error ("nc_galaxy_sd_obs_redshift_gauss_new_lsst_srd_bins: invalid type %d", type);

      return NULL;
  }

  /* Compute bin edges based on whether it's lens or source */
  if ((type == NC_GALAXY_SD_TRUE_REDSHIFT_LSST_SRD_Y1_LENS) || (type == NC_GALAXY_SD_TRUE_REDSHIFT_LSST_SRD_Y10_LENS))
  {
    /* Linear bins for lens */
    const gdouble z_min = 0.2;
    const gdouble z_max = 1.2;

    bin_edges = ncm_vector_new (n_bins + 1);

    for (guint i = 0; i <= n_bins; i++)
    {
      const gdouble zp = z_min + i * (z_max - z_min) / n_bins;

      ncm_vector_set (bin_edges, i, zp);
    }
  }
  else
  {
    const gdouble zp_max_total                  = 3.5;
    NcGalaxySDObsRedshiftGauss *gsdorgauss_temp = g_object_new (
      NC_TYPE_GALAXY_SD_OBS_REDSHIFT_GAUSS,
      "zp-lim", &(NcmDTuple2) {
      {0.0, zp_max_total}
    },
      "bin-sigma0", sigma_z,
      "zp-support-max", 20.0,
      NULL);

    ncm_model_add_submodel (NCM_MODEL (gsdorgauss_temp), NCM_MODEL (gsdtr));

    bin_edges = nc_galaxy_sd_obs_redshift_gauss_compute_equal_area_photoz_bins (
      gsdorgauss_temp,
      n_bins,
      zp_max_total);
  }

  /* Create array to hold all bins */
  bins = g_ptr_array_new_full (n_bins, (GDestroyNotify) nc_galaxy_sd_obs_redshift_gauss_free);

  /* Create all bin objects sharing the same submodel */
  for (guint i = 0; i < n_bins; i++)
  {
    const gdouble zp_min = ncm_vector_get (bin_edges, i);
    const gdouble zp_max = ncm_vector_get (bin_edges, i + 1);
    NcGalaxySDObsRedshiftGauss *gsdorgauss;
    NcmDTuple2 lim = NCM_DTUPLE2_STATIC_INIT (zp_min, zp_max);

    gsdorgauss = g_object_new (NC_TYPE_GALAXY_SD_OBS_REDSHIFT_GAUSS,
                               "zp-lim", &lim,
                               "bin-sigma0", sigma_z,
                               "zp-support-max", zp_support_max,
                               NULL);
    ncm_model_add_submodel (NCM_MODEL (gsdorgauss), NCM_MODEL (gsdtr));

    g_ptr_array_add (bins, gsdorgauss);
  }

  /* Clean up */
  ncm_vector_free (bin_edges);

  if (gsdtr_new != NULL)
    *gsdtr_new = gsdtr;
  else
    nc_galaxy_sd_true_redshift_lsst_srd_free (gsdtr);

  return bins;
}

/**
 * nc_galaxy_sd_obs_redshift_gauss_ref:
 * @gsdorgauss: a #NcGalaxySDObsRedshiftGauss
 *
 * Increases the reference count of @gsdorgauss by one.
 *
 * Returns: (transfer full): @gsdorgauss.
 */
NcGalaxySDObsRedshiftGauss *
nc_galaxy_sd_obs_redshift_gauss_ref (NcGalaxySDObsRedshiftGauss *gsdorgauss)
{
  return g_object_ref (gsdorgauss);
}

/**
 * nc_galaxy_sd_obs_redshift_gauss_free:
 * @gsdorgauss: a #NcGalaxySDObsRedshiftGauss
 *
 * Decreases the reference count of @gsdorgauss by one.
 *
 */
void
nc_galaxy_sd_obs_redshift_gauss_free (NcGalaxySDObsRedshiftGauss *gsdorgauss)
{
  g_object_unref (gsdorgauss);
}

/**
 * nc_galaxy_sd_obs_redshift_gauss_clear:
 * @gsdorgauss: a #NcGalaxySDObsRedshiftGauss
 *
 * Decreases the reference count of @gsdorgauss by one, and sets the pointer *@gsdorgauss to
 * NULL.
 *
 */
void
nc_galaxy_sd_obs_redshift_gauss_clear (NcGalaxySDObsRedshiftGauss **gsdorgauss)
{
  g_clear_object (gsdorgauss);
}

/**
 * nc_galaxy_sd_obs_redshift_gauss_set_zp_lim:
 * @gsdorgauss: a #NcGalaxySDObsRedshiftGauss
 * @zp_min: the minimum photometric redshift
 * @zp_max: the maximum photometric redshift
 *
 * Sets the minimum and maximum photometric redshifts.
 *
 */
void
nc_galaxy_sd_obs_redshift_gauss_set_zp_lim (NcGalaxySDObsRedshiftGauss *gsdorgauss, const gdouble zp_min, const gdouble zp_max)
{
  NcGalaxySDObsRedshiftGaussPrivate * const self = nc_galaxy_sd_obs_redshift_gauss_get_instance_private (gsdorgauss);

  self->zp_min = zp_min;
  self->zp_max = zp_max;

  ncm_model_state_mark_outdated (NCM_MODEL (gsdorgauss));
}

/**
 * nc_galaxy_sd_obs_redshift_gauss_get_zp_lim:
 * @gsdorgauss: a #NcGalaxySDObsRedshiftGauss
 * @zp_min: (out): the minimum photometric redshift
 * @zp_max: (out): the maximum photometric redshift
 *
 * Gets the minimum and maximum photometric redshifts.
 *
 */
void
nc_galaxy_sd_obs_redshift_gauss_get_zp_lim (NcGalaxySDObsRedshiftGauss *gsdorgauss, gdouble *zp_min, gdouble *zp_max)
{
  NcGalaxySDObsRedshiftGaussPrivate * const self = nc_galaxy_sd_obs_redshift_gauss_get_instance_private (gsdorgauss);

  *zp_min = self->zp_min;
  *zp_max = self->zp_max;
}

/**
 * nc_galaxy_sd_obs_redshift_gauss_set_use_true_z:
 * @gsdorgauss: a #NcGalaxySDObsRedshiftGauss
 * @use_true_z: whether to use the true redshift
 *
 * Sets whether to use the true redshift
 *
 */
void
nc_galaxy_sd_obs_redshift_gauss_set_use_true_z (NcGalaxySDObsRedshiftGauss *gsdorgauss, const gboolean use_true_z)
{
  NcGalaxySDObsRedshiftGaussPrivate * const self = nc_galaxy_sd_obs_redshift_gauss_get_instance_private (gsdorgauss);

  self->use_true_z = use_true_z;

  ncm_model_state_mark_outdated (NCM_MODEL (gsdorgauss));
}

/**
 * nc_galaxy_sd_obs_redshift_gauss_get_use_true_z:
 * @gsdorgauss: a #NcGalaxySDObsRedshiftGauss
 *
 * Gets whether to use the true redshift
 *
 * Returns: whether to use the true redshift
 */
gboolean
nc_galaxy_sd_obs_redshift_gauss_get_use_true_z (NcGalaxySDObsRedshiftGauss *gsdorgauss)
{
  NcGalaxySDObsRedshiftGaussPrivate * const self = nc_galaxy_sd_obs_redshift_gauss_get_instance_private (gsdorgauss);

  return self->use_true_z;
}

/**
 * nc_galaxy_sd_obs_redshift_gauss_set_bin_sigma0:
 * @gsdorgauss: a #NcGalaxySDObsRedshiftGauss
 * @bin_sigma0: base photometric redshift scatter parameter for binned analyses
 *
 * Sets the base photometric redshift scatter parameter used in binned analyses.
 * This value will be used by compute_binned_dndz and related methods. Set to -1.0
 * to unset and require explicit parameter passing.
 *
 */
void
nc_galaxy_sd_obs_redshift_gauss_set_bin_sigma0 (NcGalaxySDObsRedshiftGauss *gsdorgauss, const gdouble bin_sigma0)
{
  NcGalaxySDObsRedshiftGaussPrivate * const self = nc_galaxy_sd_obs_redshift_gauss_get_instance_private (gsdorgauss);

  self->bin_sigma0 = bin_sigma0;

  ncm_model_state_mark_outdated (NCM_MODEL (gsdorgauss));
}

/**
 * nc_galaxy_sd_obs_redshift_gauss_get_bin_sigma0:
 * @gsdorgauss: a #NcGalaxySDObsRedshiftGauss
 *
 * Gets the base photometric redshift scatter parameter for binned analyses.
 * Returns -1.0 if not set.
 *
 * Returns: the bin_sigma0 value
 */
gdouble
nc_galaxy_sd_obs_redshift_gauss_get_bin_sigma0 (NcGalaxySDObsRedshiftGauss *gsdorgauss)
{
  NcGalaxySDObsRedshiftGaussPrivate * const self = nc_galaxy_sd_obs_redshift_gauss_get_instance_private (gsdorgauss);

  return self->bin_sigma0;
}

/**
 * nc_galaxy_sd_obs_redshift_gauss_set_reltol:
 * @gsdorgauss: a #NcGalaxySDObsRedshiftGauss
 * @reltol: relative tolerance for numerical integration
 *
 * Sets the relative tolerance used in numerical integrations when computing
 * photometric redshift distributions.
 *
 */
void
nc_galaxy_sd_obs_redshift_gauss_set_reltol (NcGalaxySDObsRedshiftGauss *gsdorgauss, const gdouble reltol)
{
  NcGalaxySDObsRedshiftGaussPrivate * const self = nc_galaxy_sd_obs_redshift_gauss_get_instance_private (gsdorgauss);

  self->reltol = reltol;

  ncm_model_state_mark_outdated (NCM_MODEL (gsdorgauss));
}

/**
 * nc_galaxy_sd_obs_redshift_gauss_get_reltol:
 * @gsdorgauss: a #NcGalaxySDObsRedshiftGauss
 *
 * Gets the relative tolerance for numerical integration.
 *
 * Returns: the reltol value
 */
gdouble
nc_galaxy_sd_obs_redshift_gauss_get_reltol (NcGalaxySDObsRedshiftGauss *gsdorgauss)
{
  NcGalaxySDObsRedshiftGaussPrivate * const self = nc_galaxy_sd_obs_redshift_gauss_get_instance_private (gsdorgauss);

  return self->reltol;
}

/**
 * nc_galaxy_sd_obs_redshift_gauss_set_zp_support_max:
 * @gsdorgauss: a #NcGalaxySDObsRedshiftGauss
 * @zp_support_max: maximum photometric redshift for support
 *
 * Sets the maximum photometric redshift value to consider when computing
 * distributions and bin edges.
 *
 */
void
nc_galaxy_sd_obs_redshift_gauss_set_zp_support_max (NcGalaxySDObsRedshiftGauss *gsdorgauss, const gdouble zp_support_max)
{
  NcGalaxySDObsRedshiftGaussPrivate * const self = nc_galaxy_sd_obs_redshift_gauss_get_instance_private (gsdorgauss);

  self->zp_support_max = zp_support_max;

  ncm_model_state_mark_outdated (NCM_MODEL (gsdorgauss));
}

/**
 * nc_galaxy_sd_obs_redshift_gauss_get_zp_support_max:
 * @gsdorgauss: a #NcGalaxySDObsRedshiftGauss
 *
 * Gets the maximum photometric redshift for support.
 *
 * Returns: the zp_support_max value
 */
gdouble
nc_galaxy_sd_obs_redshift_gauss_get_zp_support_max (NcGalaxySDObsRedshiftGauss *gsdorgauss)
{
  NcGalaxySDObsRedshiftGaussPrivate * const self = nc_galaxy_sd_obs_redshift_gauss_get_instance_private (gsdorgauss);

  return self->zp_support_max;
}

/**
 * nc_galaxy_sd_obs_redshift_gauss_gen:
 * @gsdorgauss: a #NcGalaxySDObsRedshiftGauss
 * @mset: a #NcmMSet
 * @data: a #NcGalaxySDObsRedshiftData
 * @sigma0: the standard deviation of the redshift errors
 * @rng: a #NcmRNG
 *
 * Sets the required columns for the data and generates a redshift observation.
 *
 */
void
nc_galaxy_sd_obs_redshift_gauss_gen (NcGalaxySDObsRedshiftGauss *gsdorgauss, NcmMSet *mset, NcGalaxySDObsRedshiftData *data, const gdouble sigma0, NcmRNG *rng)
{
  NcGalaxySDObsRedshiftClass *klass            = NC_GALAXY_SD_OBS_REDSHIFT_GET_CLASS (gsdorgauss);
  NcGalaxySDObsRedshiftGaussData * const ldata = (NcGalaxySDObsRedshiftGaussData *) data->ldata;

  g_assert_cmpfloat (sigma0, >=, 0.0);

  ldata->sigma0 = sigma0;
  klass->gen (NC_GALAXY_SD_OBS_REDSHIFT (gsdorgauss), data, rng);
}

/**
 * nc_galaxy_sd_obs_redshift_gauss_gen1:
 * @gsdorgauss: a #NcGalaxySDObsRedshiftGauss
 * @mset: a #NcmMSet
 * @data: a #NcGalaxySDObsRedshiftData
 * @sigma0: the standard deviation of the redshift errors
 * @rng: a #NcmRNG
 *
 * Sets the required columns for the data and generates a redshift observation. See
 * nc_galaxy_sd_obs_redshift_gen1() for details.
 *
 * Returns: whether the redshift observation is within the limits.
 */
gboolean
nc_galaxy_sd_obs_redshift_gauss_gen1 (NcGalaxySDObsRedshiftGauss *gsdorgauss, NcmMSet *mset, NcGalaxySDObsRedshiftData *data, const gdouble sigma0, NcmRNG *rng)
{
  NcGalaxySDObsRedshiftClass *klass            = NC_GALAXY_SD_OBS_REDSHIFT_GET_CLASS (gsdorgauss);
  NcGalaxySDObsRedshiftGaussData * const ldata = (NcGalaxySDObsRedshiftGaussData *) data->ldata;

  g_assert_cmpfloat (sigma0, >=, 0.0);

  ldata->sigma0 = sigma0;

  return klass->gen1 (NC_GALAXY_SD_OBS_REDSHIFT (gsdorgauss), data, rng);
}

/**
 * nc_galaxy_sd_obs_redshift_gauss_data_set:
 * @gsdorgauss: a #NcGalaxySDObsRedshiftGauss
 * @data: a #NcGalaxySDObsRedshiftData
 * @zp: the observed redshift
 * @sigma0: the base standard deviation of the redshift errors
 * @sigma_z: the standard deviation of the redshift errors
 *
 * Sets the observed redshift and redshift error parameters.
 */
void
nc_galaxy_sd_obs_redshift_gauss_data_set (NcGalaxySDObsRedshiftGauss *gsdorgauss, NcGalaxySDObsRedshiftData *data, const gdouble zp, const gdouble sigma0, const gdouble sigma_z)
{
  NcGalaxySDObsRedshiftGaussData * const ldata = (NcGalaxySDObsRedshiftGaussData *) data->ldata;

  g_assert_cmpfloat (sigma0, >=, 0.0);
  g_assert_cmpfloat (sigma_z, >=, 0.0);

  ldata->zp     = zp;
  ldata->sigma0 = sigma0;
  ldata->sigma  = sigma_z;
}

/**
 * nc_galaxy_sd_obs_redshift_gauss_data_get:
 * @gsdorgauss: a #NcGalaxySDObsRedshiftGauss
 * @data: a #NcGalaxySDObsRedshiftData
 * @zp: (out): the observed redshift
 * @sigma0: (out): the base standard deviation of the redshift errors
 * @sigma_z: (out): the standard deviation of the redshift errors
 *
 * Gets the observed redshift and redshift error parameters.
 */
void
nc_galaxy_sd_obs_redshift_gauss_data_get (NcGalaxySDObsRedshiftGauss *gsdorgauss, NcGalaxySDObsRedshiftData *data, gdouble *zp, gdouble *sigma0, gdouble *sigma_z)
{
  NcGalaxySDObsRedshiftGaussData * const ldata = (NcGalaxySDObsRedshiftGaussData *) data->ldata;

  *zp      = ldata->zp;
  *sigma0  = ldata->sigma0;
  *sigma_z = ldata->sigma;
}

/**
 * nc_galaxy_sd_obs_redshift_gauss_eval_pzp:
 * @gsdorgauss: a #NcGalaxySDObsRedshiftGauss
 * @zp: photometric redshift value
 *
 * Evaluates the marginal photometric redshift distribution.
 *
 * This distribution represents the probability of observing a galaxy at photometric
 * redshift @zp, integrating over all possible true redshifts:
 *
 * $$P(z_p) = \int_0^\infty P(z) \, \mathrm{Gauss}(z_p|z,\sigma_z(z)) \, \mathrm{d}z$$
 *
 * where $\sigma_z(z) = \sigma_0 (1 + z)$ and the normalization enforces $z_p > 0$.
 *
 * The distribution is computed once and cached, so subsequent calls with the same
 * parameters will be very fast. The cache is invalidated when the model state changes.
 *
 * Returns: the probability density at photometric redshift @zp
 */
gdouble
nc_galaxy_sd_obs_redshift_gauss_eval_pzp (NcGalaxySDObsRedshiftGauss *gsdorgauss, gdouble zp)
{
  NcGalaxySDObsRedshiftGaussPrivate * const self = nc_galaxy_sd_obs_redshift_gauss_get_instance_private (gsdorgauss);

  _nc_galaxy_sd_obs_redshift_gauss_prepare_pzp (gsdorgauss);

  return ncm_stats_dist1d_eval_p (self->stats, zp);
}

/**
 * nc_galaxy_sd_obs_redshift_gauss_eval_pz_given_zp:
 * @gsdorgauss: a #NcGalaxySDObsRedshiftGauss
 * @z: true redshift value
 *
 * Evaluates the conditional true redshift distribution for a photometric bin.
 *
 * This distribution gives the probability that a galaxy observed in the photometric
 * redshift bin $[z_{p,\mathrm{min}}, z_{p,\mathrm{max}}]$ has true redshift @z:
 * $$
 * P(z|z_{p,\mathrm{min}}, z_{p,\mathrm{max}}) = \frac{P(z) \, W(z_{p,\mathrm{min}},
 * z_{p,\mathrm{max}}|z)}{N}
 * $$
 * where $W$ is the Gaussian integral over the photometric redshift bin:
 * $$
 * W = \int_{z_{p,\mathrm{min}}}^{z_{p,\mathrm{max}}} \mathrm{Gauss}(z_p|z,\sigma_z(z))
 * \, \mathrm{d}z_p,
 * $$
 * with $\sigma_z(z) = \sigma_0 (1 + z)$, and $N$ is the normalization constant.
 *
 * The distribution is computed once and cached, so subsequent calls with the same
 * parameters will be very fast. The cache is invalidated when the model state changes.
 *
 * This is particularly useful for computing the redshift distribution of galaxies in
 * weak lensing tomographic bins defined by photometric redshift cuts.
 *
 * Returns: the probability density at true redshift @z for the given photo-z bin
 */
gdouble
nc_galaxy_sd_obs_redshift_gauss_eval_pz_given_zp (NcGalaxySDObsRedshiftGauss *gsdorgauss, gdouble z)
{
  NcGalaxySDObsRedshiftGaussPrivate * const self = nc_galaxy_sd_obs_redshift_gauss_get_instance_private (gsdorgauss);

  _nc_galaxy_sd_obs_redshift_gauss_prepare_pz_given_zp (gsdorgauss);

  {
    NcmVector *xv       = ncm_spline_peek_xv (self->pz_given_zp_spline);
    const guint len     = ncm_vector_len (xv);
    const gdouble z_min = ncm_vector_get (xv, 0);
    const gdouble z_max = ncm_vector_get (xv, len - 1);

    /* Return 0 for values outside the effective support */
    if ((z < z_min) || (z > z_max))
      return 0.0;

    return ncm_spline_eval (self->pz_given_zp_spline, z);
  }
}

static NcmSpline *
_nc_galaxy_sd_obs_redshift_gauss_compute_binned_dndz (NcGalaxySDObsRedshift *gsdor, NcmVector *z_array)
{
  NcGalaxySDObsRedshiftGauss *gsdorgauss         = NC_GALAXY_SD_OBS_REDSHIFT_GAUSS (gsdor);
  NcGalaxySDObsRedshiftGaussPrivate * const self = nc_galaxy_sd_obs_redshift_gauss_get_instance_private (gsdorgauss);

  _nc_galaxy_sd_obs_redshift_gauss_prepare_pz_given_zp (gsdorgauss);

  return ncm_spline_copy (self->pz_given_zp_spline);
}

/**
 * nc_galaxy_sd_obs_redshift_gauss_compute_equal_area_photoz_bins:
 * @gsdorgauss: a #NcGalaxySDObsRedshiftGauss
 * @n_bins: number of bins to create
 * @zp_max: maximum photometric redshift to consider
 *
 * Computes equal-area photometric redshift bin edges.
 *
 * This function computes bin edges with equal areas under the convolved photo-z
 * distribution, suitable for weak lensing source samples where equal-area
 * binning in photo-z space is desired.
 *
 * The algorithm:
 *
 * 1. Computes the convolved photo-z distribution $P(z_p) = \int P(z) \times
 *    \mathrm{Gauss}(z_p|z,\sigma_z) \, \mathrm{d}z$
 * 2. Builds the cumulative distribution function (CDF)
 * 3. Inverts the CDF to find n_bins+1 edges with equal integrated probability
 *
 * Returns: (transfer full): a #NcmVector with n_bins+1 photo-z bin edges
 */
NcmVector *
nc_galaxy_sd_obs_redshift_gauss_compute_equal_area_photoz_bins (NcGalaxySDObsRedshiftGauss *gsdorgauss, guint n_bins, gdouble zp_max)
{
  NcGalaxySDObsRedshiftGaussPrivate * const self = nc_galaxy_sd_obs_redshift_gauss_get_instance_private (gsdorgauss);
  NcmVector *bin_edges;
  gdouble total_area;

  if (zp_max > self->zp_support_max)
    g_error ("Requested zp_max=%.2f exceeds the configured zp_support_max=%.2f.", zp_max, self->zp_support_max);

  /* Prepare the cached spline */
  _nc_galaxy_sd_obs_redshift_gauss_prepare_pzp (gsdorgauss);

  bin_edges  = ncm_vector_new (n_bins + 1);
  total_area = ncm_stats_dist1d_eval_pdf (self->stats, zp_max);

  for (guint i = 0; i <= n_bins; i++)
  {
    const gdouble target_area = (i * total_area) / n_bins;
    const gdouble zp_edge     = ncm_stats_dist1d_eval_inv_pdf (self->stats, target_area);

    ncm_vector_set (bin_edges, i, zp_edge);
  }

  return bin_edges;
}


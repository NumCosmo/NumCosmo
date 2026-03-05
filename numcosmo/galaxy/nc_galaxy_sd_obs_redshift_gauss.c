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
  /* Cache for p(zp) distribution */
  NcmSpline *pzp_spline;
  gdouble pzp_sigma0;
  gdouble pzp_zp_max;
  /* Cache for p(z|zp_min,zp_max) distribution */
  NcmSpline *pz_given_zp_spline;
  gdouble pz_given_zp_sigma0;
  gdouble pz_given_zp_zp_min;
  gdouble pz_given_zp_zp_max;
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
  self->pzp_spline         = NULL;
  self->pzp_sigma0         = -1.0;
  self->pzp_zp_max         = -1.0;
  self->pz_given_zp_spline = NULL;
  self->pz_given_zp_sigma0 = -1.0;
  self->pz_given_zp_zp_min = -1.0;
  self->pz_given_zp_zp_max = -1.0;
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
static NcmSpline *_nc_galaxy_sd_obs_redshift_gauss_compute_binned_dndz (NcGalaxySDObsRedshift *gsdor, gdouble sigma0, NcmVector *z_array, gdouble rel_error);
static void _nc_galaxy_sd_obs_redshift_gauss_add_submodel (NcmModel *model, NcmModel *submodel);

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


  ncm_model_class_check_params_info (model_class);

  gsdor_class->gen                 = &_nc_galaxy_sd_obs_redshift_gauss_gen;
  gsdor_class->gen1                = &_nc_galaxy_sd_obs_redshift_gauss_gen1;
  gsdor_class->prepare             = &_nc_galaxy_sd_obs_redshift_gauss_prepare;
  gsdor_class->get_integ_lim       = &_nc_galaxy_sd_obs_redshift_gauss_get_integ_lim;
  gsdor_class->integ               = &_nc_galaxy_sd_obs_redshift_gauss_integ;
  gsdor_class->data_init           = &_nc_galaxy_sd_obs_redshift_gauss_data_init;
  gsdor_class->compute_binned_dndz = &_nc_galaxy_sd_obs_redshift_gauss_compute_binned_dndz;
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

  ldata->zp     = nc_galaxy_wl_obs_get (obs, NC_GALAXY_SD_OBS_REDSHIFT_GAUSS_COL_ZP, i);
  ldata->sigma0 = nc_galaxy_wl_obs_get (obs, NC_GALAXY_SD_OBS_REDSHIFT_GAUSS_COL_SIGMA0, i);
  ldata->sigma  = nc_galaxy_wl_obs_get (obs, NC_GALAXY_SD_OBS_REDSHIFT_GAUSS_COL_SIGMA, i);
}

static void
_nc_galaxy_sd_obs_redshift_gauss_ldata_write_row (NcGalaxySDObsRedshiftData *data, NcGalaxyWLObs *obs, const guint i)
{
  NcGalaxySDObsRedshiftGaussData *ldata = (NcGalaxySDObsRedshiftGaussData *) data->ldata;

  nc_galaxy_wl_obs_set (obs, NC_GALAXY_SD_OBS_REDSHIFT_GAUSS_COL_ZP, i, ldata->zp);
  nc_galaxy_wl_obs_set (obs, NC_GALAXY_SD_OBS_REDSHIFT_GAUSS_COL_SIGMA0, i, ldata->sigma0);
  nc_galaxy_wl_obs_set (obs, NC_GALAXY_SD_OBS_REDSHIFT_GAUSS_COL_SIGMA, i, ldata->sigma);
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

  /* The normalization N is due to the zp>0 constraint. */

  return Pz * W / N;
}

static gdouble
_binned_dndz_integrand_for_gsl (const gdouble z, gpointer user_data)
{
  return _binned_dndz_integrand (user_data, z, 1.0);
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
_photoz_distribution_gsl (gdouble zp, gpointer user_data)
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

  return result;
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
 * @sigma0: base photometric redshift scatter parameter
 * @zp_max: maximum photometric redshift
 * @rel_error: relative error tolerance for integration
 *
 * Prepares and caches the marginal photometric redshift distribution.
 *
 * This distribution is computed as:
 *
 * $$P(z_p) = \int_0^\infty P(z) \, \mathrm{Gauss}(z_p|z,\sigma_z) \, \mathrm{d}z$$
 *
 * where $\sigma_z = \sigma_0 (1 + z)$ and the normalization ensures $z_p > 0$.
 *
 * The result is cached as a spline and will be reused if called again with the same
 * parameters, unless the model state has been invalidated.
 */
static void
_nc_galaxy_sd_obs_redshift_gauss_prepare_pzp (NcGalaxySDObsRedshiftGauss *gsdorgauss, gdouble sigma0, gdouble zp_max, gdouble rel_error)
{
  NcGalaxySDObsRedshiftGaussPrivate * const self = nc_galaxy_sd_obs_redshift_gauss_get_instance_private (gsdorgauss);
  gboolean need_update                           = FALSE;

  g_assert_cmpfloat (sigma0, >, 0.0);
  g_assert_cmpfloat (zp_max, >, 0.0);
  g_assert_cmpfloat (rel_error, >, 0.0);
  g_assert_nonnull (self->sdz);

  /* Check if we need to recompute */
  if ((self->pzp_spline == NULL) ||
      !ncm_model_lstate_is_update (NCM_MODEL (gsdorgauss), NC_GALAXY_SD_OBS_REDSHIFT_GAUSS_LSTATE_PZP) ||
      (self->pzp_sigma0 != sigma0) ||
      (self->pzp_zp_max != zp_max))
    need_update = TRUE;

  if (need_update)
  {
    gdouble z_min, z_max;
    PhotozDistIntegData integ_data;
    gsl_function F;

    nc_galaxy_sd_true_redshift_get_lim (self->sdz, &z_min, &z_max);

    integ_data.gsdtr      = self->sdz;
    integ_data.integrator = ncm_integral1d_ptr_new (&_photoz_distribution_integrand, NULL);
    integ_data.zp         = 0.0;
    integ_data.sigma0     = sigma0;
    integ_data.z_min      = z_min;
    integ_data.z_max      = z_max;

    F.function = &_photoz_distribution_gsl;
    F.params   = &integ_data;

    /* Clear old spline if it exists */
    ncm_spline_clear (&self->pzp_spline);

    /* Create and populate new spline */
    self->pzp_spline = NCM_SPLINE (ncm_spline_cubic_notaknot_new ());
    ncm_spline_set_func (self->pzp_spline, NCM_SPLINE_FUNCTION_SPLINE, &F, 0.0, zp_max, 0, rel_error);

    /* Cache parameters */
    self->pzp_sigma0 = sigma0;
    self->pzp_zp_max = zp_max;

    /* Clean up */
    ncm_integral1d_ptr_free (integ_data.integrator);

    /* Mark model local state as updated */
    ncm_model_lstate_set_update (NCM_MODEL (gsdorgauss), NC_GALAXY_SD_OBS_REDSHIFT_GAUSS_LSTATE_PZP);
  }
}

/**
 * _nc_galaxy_sd_obs_redshift_gauss_prepare_pz_given_zp:
 * @gsdorgauss: a #NcGalaxySDObsRedshiftGauss
 * @zp_min: minimum photometric redshift of the bin
 * @zp_max: maximum photometric redshift of the bin
 * @sigma0: base photometric redshift scatter parameter
 * @rel_error: relative error tolerance for integration
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
_nc_galaxy_sd_obs_redshift_gauss_prepare_pz_given_zp (NcGalaxySDObsRedshiftGauss *gsdorgauss, gdouble zp_min, gdouble zp_max, gdouble sigma0, gdouble rel_error)
{
  NcGalaxySDObsRedshiftGaussPrivate * const self = nc_galaxy_sd_obs_redshift_gauss_get_instance_private (gsdorgauss);
  gboolean need_update                           = FALSE;

  g_assert_cmpfloat (zp_min, >=, 0.0);
  g_assert_cmpfloat (zp_max, >, zp_min);
  g_assert_cmpfloat (sigma0, >, 0.0);
  g_assert_cmpfloat (rel_error, >, 0.0);
  g_assert_nonnull (self->sdz);

  /* Check if we need to recompute */
  if ((self->pz_given_zp_spline == NULL) ||
      !ncm_model_lstate_is_update (NCM_MODEL (gsdorgauss), NC_GALAXY_SD_OBS_REDSHIFT_GAUSS_LSTATE_PZ_GIVEN_ZP) ||
      (self->pz_given_zp_sigma0 != sigma0) ||
      (self->pz_given_zp_zp_min != zp_min) ||
      (self->pz_given_zp_zp_max != zp_max))
    need_update = TRUE;

  if (need_update)
  {
    gdouble z_min, z_max;
    BinnedDndzIntegData integ_data;
    gsl_function F;

    nc_galaxy_sd_true_redshift_get_lim (self->sdz, &z_min, &z_max);

    integ_data.gsdtr  = self->sdz;
    integ_data.zpl    = zp_min;
    integ_data.zpu    = zp_max;
    integ_data.sigma0 = sigma0;

    F.function = &_binned_dndz_integrand_for_gsl;
    F.params   = &integ_data;

    /* Clear old spline if it exists */
    ncm_spline_clear (&self->pz_given_zp_spline);

    /* Create and populate new spline */
    self->pz_given_zp_spline = NCM_SPLINE (ncm_spline_cubic_notaknot_new ());
    ncm_spline_set_func (self->pz_given_zp_spline, NCM_SPLINE_FUNCTION_SPLINE, &F, z_min, z_max, 0, rel_error);

    /* Normalize the spline */
    {
      NcmIntegral1dPtr *integrator = ncm_integral1d_ptr_new (&_binned_dndz_integrand, NULL);
      gdouble norm_error           = 0.0;
      gdouble norm;

      ncm_integral1d_set_reltol (NCM_INTEGRAL1D (integrator), rel_error);
      ncm_integral1d_ptr_set_userdata (integrator, &integ_data);

      norm = ncm_integral1d_eval (NCM_INTEGRAL1D (integrator), z_min, z_max, &norm_error);

      g_assert_cmpfloat (norm, >, 0.0);

      /* Scale the spline by 1/norm */
      {
        NcmVector *y_vec = ncm_spline_peek_yv (self->pz_given_zp_spline);
        const guint len  = ncm_vector_len (y_vec);

        for (guint i = 0; i < len; i++)
        {
          const gdouble y = ncm_vector_get (y_vec, i);

          ncm_vector_set (y_vec, i, y / norm);
        }

        ncm_spline_prepare (self->pz_given_zp_spline);
      }

      ncm_integral1d_ptr_free (integrator);
    }

    /* Cache parameters */
    self->pz_given_zp_sigma0 = sigma0;
    self->pz_given_zp_zp_min = zp_min;
    self->pz_given_zp_zp_max = zp_max;

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
 * Creates a #GPtrArray of #NcGalaxySDObsRedshiftGauss objects, one for each
 * LSST SRD bin of the specified type. All objects share the same true redshift
 * distribution submodel. The bin edges are computed based on the type:
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
  guint n_bins;
  gdouble sigma_z;
  NcmVector *bin_edges = NULL;
  GPtrArray *bins;

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
    /* Equal area bins for source */
    const gdouble zp_max_total = 3.5;
    const gdouble rel_error    = 1.0e-4;

    bin_edges = nc_galaxy_sd_obs_redshift_gauss_compute_equal_area_photoz_bins (NC_GALAXY_SD_TRUE_REDSHIFT (gsdtr),
                                                                                n_bins,
                                                                                sigma_z,
                                                                                zp_max_total,
                                                                                rel_error);
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
                               NULL);

    /* All bins share the same true redshift distribution submodel */
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
 * @sigma0: base photometric redshift scatter parameter
 * @zp_max: maximum photometric redshift for the distribution
 * @rel_error: relative error tolerance for integration (typical value: 1e-4)
 *
 * Evaluates the marginal photometric redshift distribution.
 *
 * This distribution represents the probability of observing a galaxy at
 * photometric redshift @zp, integrating over all possible true redshifts:
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
nc_galaxy_sd_obs_redshift_gauss_eval_pzp (NcGalaxySDObsRedshiftGauss *gsdorgauss, gdouble zp, gdouble sigma0, gdouble zp_max, gdouble rel_error)
{
  NcGalaxySDObsRedshiftGaussPrivate * const self = nc_galaxy_sd_obs_redshift_gauss_get_instance_private (gsdorgauss);

  g_assert_cmpfloat (zp, >=, 0.0);
  g_assert (zp <= zp_max * (1.0 + 1.0e-10));

  /* Prepare the distribution (uses cache if possible) */
  _nc_galaxy_sd_obs_redshift_gauss_prepare_pzp (gsdorgauss, sigma0, zp_max, rel_error);

  /* Evaluate the spline at zp */
  return ncm_spline_eval (self->pzp_spline, zp);
}

/**
 * nc_galaxy_sd_obs_redshift_gauss_eval_pz_given_zp:
 * @gsdorgauss: a #NcGalaxySDObsRedshiftGauss
 * @z: true redshift value
 * @zp_min: minimum photometric redshift of the bin
 * @zp_max: maximum photometric redshift of the bin
 * @sigma0: base photometric redshift scatter parameter
 * @rel_error: relative error tolerance for integration (typical value: 1e-4)
 *
 * Evaluates the conditional true redshift distribution for a photometric bin.
 *
 * This distribution gives the probability that a galaxy observed in the
 * photometric redshift bin $[z_{p,\mathrm{min}}, z_{p,\mathrm{max}}]$ has
 * true redshift @z:
 *
 * $$P(z|z_{p,\mathrm{min}}, z_{p,\mathrm{max}}) = \frac{P(z) \, W(z_{p,\mathrm{min}}, z_{p,\mathrm{max}}|z)}{N}$$
 *
 * where $W$ is the Gaussian integral over the photometric redshift bin:
 *
 * $$W = \int_{z_{p,\mathrm{min}}}^{z_{p,\mathrm{max}}} \mathrm{Gauss}(z_p|z,\sigma_z(z)) \, \mathrm{d}z_p$$
 *
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
nc_galaxy_sd_obs_redshift_gauss_eval_pz_given_zp (NcGalaxySDObsRedshiftGauss *gsdorgauss, gdouble z, gdouble zp_min, gdouble zp_max, gdouble sigma0, gdouble rel_error)
{
  NcGalaxySDObsRedshiftGaussPrivate * const self = nc_galaxy_sd_obs_redshift_gauss_get_instance_private (gsdorgauss);
  gdouble z_min, z_max;

  nc_galaxy_sd_true_redshift_get_lim (self->sdz, &z_min, &z_max);

  g_assert_cmpfloat (z, >=, z_min);
  g_assert_cmpfloat (z, <=, z_max);

  /* Prepare the distribution (uses cache if possible) */
  _nc_galaxy_sd_obs_redshift_gauss_prepare_pz_given_zp (gsdorgauss, zp_min, zp_max, sigma0, rel_error);

  /* Evaluate the spline at z */
  return ncm_spline_eval (self->pz_given_zp_spline, z);
}

static NcmSpline *
_nc_galaxy_sd_obs_redshift_gauss_compute_binned_dndz (NcGalaxySDObsRedshift *gsdor, gdouble sigma0, NcmVector *z_array, gdouble rel_error)
{
  NcGalaxySDObsRedshiftGauss *gsdorgauss         = NC_GALAXY_SD_OBS_REDSHIFT_GAUSS (gsdor);
  NcGalaxySDObsRedshiftGaussPrivate * const self = nc_galaxy_sd_obs_redshift_gauss_get_instance_private (gsdorgauss);

  g_assert_cmpfloat (sigma0, >, 0.0);
  g_assert_nonnull (z_array);
  g_assert_cmpuint (ncm_vector_len (z_array), >, 1);

  {
    const gdouble zpl              = self->zp_min;
    const gdouble zpu              = self->zp_max;
    NcmIntegral1dPtr *integrator   = ncm_integral1d_ptr_new (&_binned_dndz_integrand, NULL);
    NcmSpline *spline              = NCM_SPLINE (ncm_spline_cubic_notaknot_new ());
    BinnedDndzIntegData integ_data = {
      .gsdtr  = self->sdz,
      .zpl    = zpl,
      .zpu    = zpu,
      .sigma0 = sigma0
    };
    gdouble z_min, z_max;
    gdouble norm_error = 0.0, norm;

    nc_galaxy_sd_true_redshift_get_lim (self->sdz, &z_min, &z_max);
    ncm_integral1d_set_reltol (NCM_INTEGRAL1D (integrator), rel_error);
    ncm_integral1d_ptr_set_userdata (integrator, &integ_data);

    norm = ncm_integral1d_eval (NCM_INTEGRAL1D (integrator), z_min, z_max, &norm_error);

    g_assert_cmpfloat (norm, >, 0.0);

    if (z_array != NULL)
    {
      const guint n_points = ncm_vector_len (z_array);
      NcmVector *dndz_vec  = ncm_vector_new (n_points);

      for (guint i = 0; i < n_points; i++)
      {
        const gdouble z    = ncm_vector_get (z_array, i);
        const gdouble dndz = _binned_dndz_integrand (&integ_data, z, 1.0) / norm;

        ncm_vector_set (dndz_vec, i, dndz);
      }

      ncm_spline_set (spline, z_array, dndz_vec, TRUE);
      ncm_vector_free (dndz_vec);
    }
    else
    {
      gsl_function F = {
        .function = &_binned_dndz_integrand_for_gsl,
        .params   = &integ_data
      };

      ncm_spline_set_func (spline, NCM_SPLINE_FUNCTION_SPLINE, &F, z_min, z_max, 0, rel_error);
    }

    ncm_integral1d_ptr_free (integrator);

    return spline;
  }
}

/**
 * nc_galaxy_sd_obs_redshift_gauss_compute_equal_area_photoz_bins:
 * @gsdtr: a #NcGalaxySDTrueRedshift
 * @n_bins: number of bins to create
 * @sigma0: base photometric redshift scatter parameter
 * @zp_max: maximum photometric redshift to consider
 * @rel_error: relative error tolerance for integration
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
nc_galaxy_sd_obs_redshift_gauss_compute_equal_area_photoz_bins (NcGalaxySDTrueRedshift *gsdtr, guint n_bins, gdouble sigma0, gdouble zp_max, gdouble rel_error)
{
  gdouble z_min, z_max;

  g_assert_cmpuint (n_bins, >, 0);
  g_assert_cmpfloat (sigma0, >, 0.0);
  g_assert_cmpfloat (zp_max, >, 0.0);

  nc_galaxy_sd_true_redshift_get_lim (gsdtr, &z_min, &z_max);
  {
    NcmSplineCubicNotaknot *spline = ncm_spline_cubic_notaknot_new ();
    PhotozDistIntegData integ_data = {
      .gsdtr      = gsdtr,
      .integrator = ncm_integral1d_ptr_new (&_photoz_distribution_integrand, NULL),
      .zp         = 0.0,
      .sigma0     = sigma0,
      .z_min      = z_min,
      .z_max      = z_max
    };
    gsl_function F;

    F.function = &_photoz_distribution_gsl;
    F.params   = &integ_data;

    ncm_spline_set_func (NCM_SPLINE (spline), NCM_SPLINE_FUNCTION_SPLINE, &F, 0.0, zp_max, 0, rel_error);

    {
      /* Transform P(zp) spline to -2*ln(P(zp)) for NcmStatsDist1dSpline */
      NcmVector *xv        = ncm_spline_peek_xv (NCM_SPLINE (spline));
      NcmVector *yv        = ncm_spline_peek_yv (NCM_SPLINE (spline));
      NcmVector *m2lnP_vec = ncm_vector_new (ncm_vector_len (yv));
      NcmSpline *m2lnP_spline;
      NcmStatsDist1dSpline *stats;
      NcmVector *bin_edges;
      const guint len = ncm_vector_len (yv);
      gdouble total_area;

      for (guint i = 0; i < len; i++)
      {
        const gdouble P = ncm_vector_get (yv, i);

        g_assert_cmpfloat (P, >, 0.0);
        ncm_vector_set (m2lnP_vec, i, -2.0 * log (P));
      }

      m2lnP_spline = NCM_SPLINE (ncm_spline_cubic_notaknot_new ());
      ncm_spline_set (m2lnP_spline, xv, m2lnP_vec, TRUE);

      stats     = ncm_stats_dist1d_spline_new (m2lnP_spline);
      bin_edges = ncm_vector_new (n_bins + 1);

      g_object_set (stats, "abstol", 1.0e-50, NULL);
      ncm_stats_dist1d_prepare (NCM_STATS_DIST1D (stats));
      total_area = ncm_stats_dist1d_eval_pdf (NCM_STATS_DIST1D (stats), zp_max);

      for (guint i = 0; i <= n_bins; i++)
      {
        const gdouble target_area = (i * total_area) / n_bins;
        const gdouble zp_edge     = ncm_stats_dist1d_eval_inv_pdf (NCM_STATS_DIST1D (stats), target_area);

        ncm_vector_set (bin_edges, i, zp_edge);
      }

      ncm_stats_dist1d_free (NCM_STATS_DIST1D (stats));
      ncm_spline_free (m2lnP_spline);
      ncm_vector_free (m2lnP_vec);
      ncm_spline_free (NCM_SPLINE (spline));
      ncm_integral1d_ptr_free (integ_data.integrator);

      return bin_edges;
    }
  }
}


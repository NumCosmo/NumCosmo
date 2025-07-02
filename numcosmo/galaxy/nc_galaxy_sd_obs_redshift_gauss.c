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
#include "galaxy/nc_galaxy_sd_obs_redshift.h"
#include "galaxy/nc_galaxy_sd_obs_redshift_gauss.h"
#include "galaxy/nc_galaxy_sd_true_redshift.h"
#include "math/ncm_vector.h"
#include "math/ncm_rng.h"
#include "math/ncm_dtuple.h"

typedef struct _NcGalaxySDObsRedshiftGaussPrivate
{
  NcGalaxySDTrueRedshift *sdz;
  gdouble zp_min;
  gdouble zp_max;
  gboolean use_true_z;
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
  PROP_LIM,
  PROP_USE_TRUE_Z,
  PROP_LEN,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcGalaxySDObsRedshiftGauss, nc_galaxy_sd_obs_redshift_gauss, NC_TYPE_GALAXY_SD_OBS_REDSHIFT);

static void
nc_galaxy_sd_obs_redshift_gauss_init (NcGalaxySDObsRedshiftGauss *gsdorgauss)
{
  NcGalaxySDObsRedshiftGaussPrivate * const self = nc_galaxy_sd_obs_redshift_gauss_get_instance_private (gsdorgauss);

  self->sdz = NULL;
}

static void
_nc_galaxy_sd_obs_redshift_gauss_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcGalaxySDObsRedshiftGauss *gsdorgauss = NC_GALAXY_SD_OBS_REDSHIFT_GAUSS (object);

  g_return_if_fail (NC_IS_GALAXY_SD_OBS_REDSHIFT_GAUSS (gsdorgauss));

  switch (prop_id)
  {
    case PROP_LIM:
    {
      NcmDTuple2 *lim = g_value_get_boxed (value);

      if (lim == NULL)
        g_error ("_nc_galaxy_sd_obs_redshift_gauss_set_property: lim is NULL");

      nc_galaxy_sd_obs_redshift_gauss_set_lim (gsdorgauss, lim->elements[0], lim->elements[1]);
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
  NcGalaxySDObsRedshiftGauss *gsdorgauss = NC_GALAXY_SD_OBS_REDSHIFT_GAUSS (object);

  g_return_if_fail (NC_IS_GALAXY_SD_OBS_REDSHIFT_GAUSS (gsdorgauss));

  switch (prop_id)
  {
    case PROP_LIM:
    {
      gdouble zp_min, zp_max;

      nc_galaxy_sd_obs_redshift_gauss_get_lim (gsdorgauss, &zp_min, &zp_max);

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
static void _nc_galaxy_sd_obs_redshift_gauss_get_lim (NcGalaxySDObsRedshift *gsdor, NcGalaxySDObsRedshiftData *data, gdouble *z_min, gdouble *z_max);
static NcGalaxySDObsRedshiftIntegrand *_nc_galaxy_sd_obs_redshift_gauss_integ (NcGalaxySDObsRedshift *gsdor);
static void _nc_galaxy_sd_obs_redshift_gauss_data_init (NcGalaxySDObsRedshift *gsdor, NcGalaxySDObsRedshiftData *data);
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
   * NcGalaxySDObsRedshiftGauss:lim:
   *
   * Galaxy sample photometric redshift distribution limits.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_LIM,
                                   g_param_spec_boxed ("lim",
                                                       NULL,
                                                       "Galaxy sample photometric redshift distribution limits",
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

  gsdor_class->gen          = &_nc_galaxy_sd_obs_redshift_gauss_gen;
  gsdor_class->gen1         = &_nc_galaxy_sd_obs_redshift_gauss_gen1;
  gsdor_class->prepare      = &_nc_galaxy_sd_obs_redshift_gauss_prepare;
  gsdor_class->get_lim      = &_nc_galaxy_sd_obs_redshift_gauss_get_lim;
  gsdor_class->integ        = &_nc_galaxy_sd_obs_redshift_gauss_integ;
  gsdor_class->data_init    = &_nc_galaxy_sd_obs_redshift_gauss_data_init;
  model_class->add_submodel = &_nc_galaxy_sd_obs_redshift_gauss_add_submodel;
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
_nc_galaxy_sd_obs_redshift_gauss_get_lim (NcGalaxySDObsRedshift *gsdor, NcGalaxySDObsRedshiftData *data, gdouble *z_min, gdouble *z_max)
{
  NcGalaxySDObsRedshiftGauss *gsdorgauss         = NC_GALAXY_SD_OBS_REDSHIFT_GAUSS (gsdor);
  NcGalaxySDObsRedshiftGaussPrivate * const self = nc_galaxy_sd_obs_redshift_gauss_get_instance_private (gsdorgauss);

  *z_min = self->zp_min;
  *z_max = self->zp_max;
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
_nc_galaxy_sd_obs_redshift_gauss_integ_f (gpointer callback_data, const gdouble z, NcGalaxySDObsRedshiftData *data)
{
  const struct _IntegData *int_data              = (struct _IntegData *) callback_data;
  NcGalaxySDObsRedshiftGaussPrivate * const self = nc_galaxy_sd_obs_redshift_gauss_get_instance_private (int_data->gsdorgauss);
  NcGalaxySDObsRedshiftGaussData * const ldata   = (NcGalaxySDObsRedshiftGaussData *) data->ldata;
  const gdouble zp                               = ldata->zp;

  if (self->use_true_z)
  {
    const gdouble sigmaz = ldata->sigma0 * (1.0 + z);
    const gdouble norm   = sqrt (2.0 * M_PI) * sigmaz * 0.5 * (1.0 + erf (z / (M_SQRT2 * sigmaz)));
    const gdouble int_z  = nc_galaxy_sd_true_redshift_integ (self->sdz, z);
    const gdouble int_zp = exp (-0.5 * gsl_pow_2 ((zp - z) / sigmaz)) / norm;

    return int_z * int_zp;
  }
  else
  {
    const gdouble norm   = sqrt (2.0 * M_PI) * ldata->sigma;
    const gdouble chi2   = exp (-0.5 * gsl_pow_2 ((zp - z) / ldata->sigma));
    const gdouble int_zp = chi2 / norm;

    return int_zp;
  }
}

static NcGalaxySDObsRedshiftIntegrand *
_nc_galaxy_sd_obs_redshift_gauss_integ (NcGalaxySDObsRedshift *gsdor)
{
  NcGalaxySDObsRedshiftGauss *gsdorgauss = NC_GALAXY_SD_OBS_REDSHIFT_GAUSS (gsdor);
  struct _IntegData *int_data            = g_new0 (struct _IntegData, 1);
  NcGalaxySDObsRedshiftIntegrand *integ  = nc_galaxy_sd_obs_redshift_integrand_new (&_nc_galaxy_sd_obs_redshift_gauss_integ_f,
                                                                                    &_integ_data_free,
                                                                                    &_integ_data_copy,
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
 * nc_galaxy_sd_obs_redshift_gauss_new:
 * @sdz: a #NcGalaxySDTrueRedshift.
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
                                                         "lim", &lim,
                                                         NULL);

  ncm_model_add_submodel (NCM_MODEL (gsdorgauss), NCM_MODEL (sdz));

  return gsdorgauss;
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
 * nc_galaxy_sd_obs_redshift_gauss_set_lim:
 * @gsdorgauss: a #NcGalaxySDObsRedshiftGauss
 * @zp_min: the minimum redshift
 * @zp_max: the maximum redshift
 *
 * Sets the minimum and maximum redshifts.
 */
void
nc_galaxy_sd_obs_redshift_gauss_set_lim (NcGalaxySDObsRedshiftGauss *gsdorgauss, const gdouble zp_min, const gdouble zp_max)
{
  NcGalaxySDObsRedshiftGaussPrivate * const self = nc_galaxy_sd_obs_redshift_gauss_get_instance_private (gsdorgauss);

  self->zp_min = zp_min;
  self->zp_max = zp_max;

  ncm_model_state_mark_outdated (NCM_MODEL (gsdorgauss));
}

/**
 * nc_galaxy_sd_obs_redshift_gauss_get_lim:
 * @gsdorgauss: a #NcGalaxySDObsRedshiftGauss
 * @zp_min: (out): the minimum redshift
 * @zp_max: (out): the maximum redshift
 *
 * Gets the minimum and maximum redshifts.
 *
 */
void
nc_galaxy_sd_obs_redshift_gauss_get_lim (NcGalaxySDObsRedshiftGauss *gsdorgauss, gdouble *zp_min, gdouble *zp_max)
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
 * Sets the observed redshift and the standard deviation of the redshift errors.
 *
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
 * @zp: the observed redshift
 * @sigma0: the base standard deviation of the redshift errors
 * @sigma_z: the standard deviation of the redshift errors
 *
 * Gets the observed redshift and the standard deviation of the redshift errors.
 *
 */
void
nc_galaxy_sd_obs_redshift_gauss_data_get (NcGalaxySDObsRedshiftGauss *gsdorgauss, NcGalaxySDObsRedshiftData *data, gdouble *zp, gdouble *sigma0, gdouble *sigma_z)
{
  NcGalaxySDObsRedshiftGaussData * const ldata = (NcGalaxySDObsRedshiftGaussData *) data->ldata;

  *zp      = ldata->zp;
  *sigma0  = ldata->sigma0;
  *sigma_z = ldata->sigma;
}


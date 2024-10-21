/***************************************************************************
 *            nc_galaxy_sd_obs_redshift_gauss.c
 *
 *  Thu Aug 1 20:03:55 2024
 *  Copyright  2024  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * nc_galaxy_sd_obs_redshift_gauss.c
 * Copyright (C) 2024 Caio Lima de Oliveira <caiolimadeoliveira@pm.me>
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
 * SECTION:nc_galaxy_sd_obs_redshift_gauss
 * @title: NcGalaxySDObsRedshiftGauss
 * @short_description: Class describing photometric redshift observations with gaussian errors.
 * @stability: Unstable
 *
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

typedef struct _NcGalaxySDObsRedshiftGaussPrivate
{
  NcGalaxySDTrueRedshift *sdz;
} NcGalaxySDObsRedshiftGaussPrivate;

struct _NcGalaxySDObsRedshiftGauss
{
  NcGalaxySDObsRedshift parent_instance;
};

typedef struct _NcGalaxySDObsRedshiftGaussData
{
  gdouble zp;
  gdouble sigma;
} NcGalaxySDObsRedshiftGaussData;

enum
{
  PROP_0,
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
_nc_galaxy_sd_obs_redshift_gauss_dispose (GObject *object)
{
  /* Chain up: end */
  G_OBJECT_CLASS (nc_galaxy_sd_obs_redshift_gauss_parent_class)->dispose (object);
}

static void
nc_galaxy_sd_obs_redshift_gauss_finalize (GObject *object)
{
  /* Chain up: end */
  G_OBJECT_CLASS (nc_galaxy_sd_obs_redshift_gauss_parent_class)->finalize (object);
}

static void _nc_galaxy_sd_obs_redshift_gauss_gen (NcGalaxySDObsRedshift *gsdor, NcGalaxySDObsRedshiftData *data, NcmRNG *rng);
static NcGalaxySDObsRedshiftIntegrand *_nc_galaxy_sd_obs_redshift_gauss_integ (NcGalaxySDObsRedshift *gsdor);
static void _nc_galaxy_sd_obs_redshift_gauss_data_init (NcGalaxySDObsRedshift *gsdor, NcGalaxySDObsRedshiftData *data);
static void _nc_galaxy_sd_obs_redshift_gauss_add_submodel (NcmModel *model, NcmModel *submodel);

static void
nc_galaxy_sd_obs_redshift_gauss_class_init (NcGalaxySDObsRedshiftGaussClass *klass)
{
  NcGalaxySDObsRedshiftClass *gsdor_class = NC_GALAXY_SD_OBS_REDSHIFT_CLASS (klass);
  GObjectClass *object_class              = G_OBJECT_CLASS (klass);
  NcmModelClass *model_class              = NCM_MODEL_CLASS (klass);

  object_class->dispose  = &_nc_galaxy_sd_obs_redshift_gauss_dispose;
  object_class->finalize = &nc_galaxy_sd_obs_redshift_gauss_finalize;

  ncm_model_class_set_name_nick (model_class, "Gaussian Observed Redshift", "GalaxySDObsRedshiftGauss");
  ncm_model_class_add_params (model_class, 0, 0, PROP_LEN);

  ncm_model_class_check_params_info (model_class);

  gsdor_class->gen          = &_nc_galaxy_sd_obs_redshift_gauss_gen;
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
  gdouble sigma                                  = ldata->sigma;
  gdouble z_min                                  = 0.0;
  gdouble z_max                                  = 0.0;
  gdouble zp;
  gdouble z;

  if (!nc_galaxy_sd_true_redshift_get_lim (self->sdz, &z_min, &z_max))
    g_error ("Failed to get redshift limits.");

  z = nc_galaxy_sd_true_redshift_gen (self->sdz, rng);

  {
    const gdouble sigmaz = sigma * (1.0 + z);

    do {
      zp = ncm_rng_gaussian_gen (rng, z, sigmaz);
    } while ((zp > z_max) || (zp < z_min));

    data->z   = z;
    ldata->zp = zp;
  }
}

struct _IntegData
{
  NcGalaxySDObsRedshiftGauss *gsdorgauss;
};

static gpointer
_integ_data_copy (gpointer idata)
{
  struct _IntegData *new_idata = g_new0 (struct _IntegData, 1);

  *new_idata = *(struct _IntegData *) idata;

  return new_idata;
}

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
  const gdouble sigma                            = ldata->sigma;
  const gdouble sigmaz                           = sigma * (1.0 + z);
  const gdouble norm                             =  (sqrt (2.0 * M_PI) * sigmaz);
  const gdouble int_z                            = nc_galaxy_sd_true_redshift_integ (self->sdz, z);
  const gdouble int_zp                           = exp (-0.5 * gsl_pow_2 ((zp - z) / sigmaz)) / norm;

  return int_z * int_zp;
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

  ldata->zp    = nc_galaxy_wl_obs_get (obs, NC_GALAXY_SD_OBS_REDSHIFT_GAUSS_COL_ZP, i);
  ldata->sigma = nc_galaxy_wl_obs_get (obs, NC_GALAXY_SD_OBS_REDSHIFT_GAUSS_COL_SIGMA, i);
}

static void
_nc_galaxy_sd_obs_redshift_gauss_ldata_write_row (NcGalaxySDObsRedshiftData *data, NcGalaxyWLObs *obs, const guint i)
{
  NcGalaxySDObsRedshiftGaussData *ldata = (NcGalaxySDObsRedshiftGaussData *) data->ldata;

  nc_galaxy_wl_obs_set (obs, NC_GALAXY_SD_OBS_REDSHIFT_GAUSS_COL_ZP, i, ldata->zp);
  nc_galaxy_wl_obs_set (obs, NC_GALAXY_SD_OBS_REDSHIFT_GAUSS_COL_SIGMA, i, ldata->sigma);
}

static void
_nc_galaxy_sd_obs_redshift_gauss_ldata_required_columns (NcGalaxySDObsRedshiftData *data, GList *columns)
{
  columns = g_list_append (columns, g_strdup (NC_GALAXY_SD_OBS_REDSHIFT_GAUSS_COL_ZP));
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
nc_galaxy_sd_obs_redshift_gauss_new (NcGalaxySDTrueRedshift *sdz)
{
  NcGalaxySDObsRedshiftGauss *gsdorgauss = g_object_new (NC_TYPE_GALAXY_SD_OBS_REDSHIFT_GAUSS, NULL);

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
 * nc_galaxy_sd_obs_redshift_gauss_gen:
 * @gsdorgauss: a #NcGalaxySDObsRedshiftGauss
 * @data: a #NcGalaxySDObsRedshiftData
 * @sigma_z: the standard deviation of the redshift errors
 * @rng: a #NcmRNG
 *
 * Sets the required columns for the data and generates a redshift observation.
 *
 */
void
nc_galaxy_sd_obs_redshift_gauss_gen (NcGalaxySDObsRedshiftGauss *gsdorgauss, NcmMSet *mset, NcGalaxySDObsRedshiftData *data, const gdouble sigma_z, NcmRNG *rng)
{
  NcGalaxySDObsRedshiftClass *klass            = NC_GALAXY_SD_OBS_REDSHIFT_GET_CLASS (gsdorgauss);
  NcGalaxySDObsRedshiftGaussData * const ldata = (NcGalaxySDObsRedshiftGaussData *) data->ldata;

  ldata->sigma = sigma_z;
  klass->gen (NC_GALAXY_SD_OBS_REDSHIFT (gsdorgauss), data, rng);
}

/**
 * nc_galaxy_sd_obs_redshift_gauss_data_set:
 * @gsdorgauss: a #NcGalaxySDObsRedshiftGauss
 * @data: a #NcGalaxySDObsRedshiftData
 * @zp: the observed redshift
 * @sigma_z: the standard deviation of the redshift errors
 *
 * Sets the observed redshift and the standard deviation of the redshift errors.
 *
 */
void
nc_galaxy_sd_obs_redshift_gauss_data_set (NcGalaxySDObsRedshiftGauss *gsdorgauss, NcGalaxySDObsRedshiftData *data, const gdouble zp, const gdouble sigma_z)
{
  NcGalaxySDObsRedshiftGaussData * const ldata = (NcGalaxySDObsRedshiftGaussData *) data->ldata;

  ldata->zp    = zp;
  ldata->sigma = sigma_z;
}

/**
 * nc_galaxy_sd_obs_redshift_gauss_data_get:
 * @gsdorgauss: a #NcGalaxySDObsRedshiftGauss
 * @data: a #NcGalaxySDObsRedshiftData
 * @zp: the observed redshift
 * @sigma_z: the standard deviation of the redshift errors
 *
 * Gets the observed redshift and the standard deviation of the redshift errors.
 *
 */
void
nc_galaxy_sd_obs_redshift_gauss_data_get (NcGalaxySDObsRedshiftGauss *gsdorgauss, NcGalaxySDObsRedshiftData *data, gdouble *zp, gdouble *sigma_z)
{
  NcGalaxySDObsRedshiftGaussData * const ldata = (NcGalaxySDObsRedshiftGaussData *) data->ldata;

  *zp      = ldata->zp;
  *sigma_z = ldata->sigma;
}


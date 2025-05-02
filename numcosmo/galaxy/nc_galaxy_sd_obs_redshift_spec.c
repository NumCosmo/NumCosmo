/***************************************************************************
 *            nc_galaxy_sd_obs_redshift_spec.c
 *
 *  Thu Aug 1 15:10:22 2024
 *  Copyright  2024  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 *  Copyright  2024  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_galaxy_sd_obs_redshift_spec.c
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
 * NcGalaxySDObsRedshiftSpec:
 *
 * Class describing spectroscopic redshift observations.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "galaxy/nc_galaxy_sd_obs_redshift.h"
#include "galaxy/nc_galaxy_sd_obs_redshift_spec.h"
#include "galaxy/nc_galaxy_sd_true_redshift.h"
#include "math/ncm_rng.h"
#include "math/ncm_dtuple.h"

typedef struct _NcGalaxySDObsRedshiftSpecPrivate
{
  NcGalaxySDTrueRedshift *sdz;
  gdouble z_min;
  gdouble z_max;
} NcGalaxySDObsRedshiftSpecPrivate;

struct _NcGalaxySDObsRedshiftSpec
{
  NcGalaxySDObsRedshift parent_instance;
};

typedef struct _NcGalaxySDObsRedshiftSpecData
{
  gint placeholder;
} NcGalaxySDObsRedshiftSpecData;


enum
{
  PROP_0,
  PROP_LIM,
  PROP_LEN,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcGalaxySDObsRedshiftSpec, nc_galaxy_sd_obs_redshift_spec, NC_TYPE_GALAXY_SD_OBS_REDSHIFT);

static void
nc_galaxy_sd_obs_redshift_spec_init (NcGalaxySDObsRedshiftSpec *gsdorspec)
{
  NcGalaxySDObsRedshiftSpecPrivate * const self = nc_galaxy_sd_obs_redshift_spec_get_instance_private (gsdorspec);

  self->sdz   = NULL;
  self->z_min = 0.0;
  self->z_max = 0.0;
}

static void
_nc_galaxy_sd_obs_redshift_spec_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcGalaxySDObsRedshiftSpec *gsdorspec = NC_GALAXY_SD_OBS_REDSHIFT_SPEC (object);

  g_return_if_fail (NC_IS_GALAXY_SD_OBS_REDSHIFT_SPEC (gsdorspec));

  switch (prop_id)
  {
    case PROP_LIM:
    {
      NcmDTuple2 *lim = g_value_get_boxed (value);

      if (lim == NULL)
        g_error ("_nc_galaxy_sd_obs_redshift_spec_set_property: lim is NULL");

      nc_galaxy_sd_obs_redshift_spec_set_lim (gsdorspec, lim->elements[0], lim->elements[1]);
      break;
    }
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_galaxy_sd_obs_redshift_spec_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcGalaxySDObsRedshiftSpec *gsdorspec = NC_GALAXY_SD_OBS_REDSHIFT_SPEC (object);

  g_return_if_fail (NC_IS_GALAXY_SD_OBS_REDSHIFT_SPEC (gsdorspec));

  switch (prop_id)
  {
    case PROP_LIM:
    {
      gdouble zp_min, zp_max;

      nc_galaxy_sd_obs_redshift_spec_get_lim (gsdorspec, &zp_min, &zp_max);

      g_value_take_boxed (value, ncm_dtuple2_new (zp_min, zp_max));
      break;
    }
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_galaxy_sd_obs_redshift_spec_dispose (GObject *object)
{
  /* Chain up: end */
  G_OBJECT_CLASS (nc_galaxy_sd_obs_redshift_spec_parent_class)->dispose (object);
}

static void
nc_galaxy_sd_obs_redshift_spec_finalize (GObject *object)
{
  /* Chain up: end */
  G_OBJECT_CLASS (nc_galaxy_sd_obs_redshift_spec_parent_class)->finalize (object);
}

static void _nc_galaxy_sd_obs_redshift_spec_gen (NcGalaxySDObsRedshift *gsdor, NcGalaxySDObsRedshiftData *data, NcmRNG *rng);
static gboolean _nc_galaxy_sd_obs_redshift_spec_gen1 (NcGalaxySDObsRedshift *gsdor, NcGalaxySDObsRedshiftData *data, NcmRNG *rng);
static void _nc_galaxy_sd_obs_redshift_spec_prepare (NcGalaxySDObsRedshift *gsdor, NcGalaxySDObsRedshiftData *data);
static void _nc_galaxy_sd_obs_redshift_spec_get_lim (NcGalaxySDObsRedshift *gsdor, NcGalaxySDObsRedshiftData *data, gdouble *z_min, gdouble *z_max);
static NcGalaxySDObsRedshiftIntegrand *_nc_galaxy_sd_obs_redshift_spec_integ (NcGalaxySDObsRedshift *gsdor);
static void _nc_galaxy_sd_obs_redshift_spec_data_init (NcGalaxySDObsRedshift *gsdor, NcGalaxySDObsRedshiftData *data);
static void _nc_galaxy_sd_obs_redshift_spec_add_submodel (NcmModel *model, NcmModel *submodel);

static void
nc_galaxy_sd_obs_redshift_spec_class_init (NcGalaxySDObsRedshiftSpecClass *klass)
{
  NcGalaxySDObsRedshiftClass *gsdor_class = NC_GALAXY_SD_OBS_REDSHIFT_CLASS (klass);
  GObjectClass *object_class              = G_OBJECT_CLASS (klass);
  NcmModelClass *model_class              = NCM_MODEL_CLASS (klass);

  model_class->set_property = &_nc_galaxy_sd_obs_redshift_spec_set_property;
  model_class->get_property = &_nc_galaxy_sd_obs_redshift_spec_get_property;
  object_class->dispose     = &_nc_galaxy_sd_obs_redshift_spec_dispose;
  object_class->finalize    = &nc_galaxy_sd_obs_redshift_spec_finalize;

  ncm_model_class_set_name_nick (model_class, "Spectroscopic Observed Redshift", "GalaxySDObsRedshiftSpec");
  ncm_model_class_add_params (model_class, 0, 0, PROP_LEN);

  /**
   * NcGalaxySDObsRedshiftSpec:lim:
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

  ncm_model_class_check_params_info (model_class);

  gsdor_class->gen          = &_nc_galaxy_sd_obs_redshift_spec_gen;
  gsdor_class->gen1         = &_nc_galaxy_sd_obs_redshift_spec_gen1;
  gsdor_class->prepare      = &_nc_galaxy_sd_obs_redshift_spec_prepare;
  gsdor_class->get_lim      = &_nc_galaxy_sd_obs_redshift_spec_get_lim;
  gsdor_class->integ        = &_nc_galaxy_sd_obs_redshift_spec_integ;
  gsdor_class->data_init    = &_nc_galaxy_sd_obs_redshift_spec_data_init;
  model_class->add_submodel = &_nc_galaxy_sd_obs_redshift_spec_add_submodel;
}

static void
_nc_galaxy_sd_obs_redshift_spec_gen (NcGalaxySDObsRedshift *gsdor, NcGalaxySDObsRedshiftData *data, NcmRNG *rng)
{
  NcGalaxySDObsRedshiftSpec *gsdorspec          = NC_GALAXY_SD_OBS_REDSHIFT_SPEC (gsdor);
  NcGalaxySDObsRedshiftSpecPrivate * const self = nc_galaxy_sd_obs_redshift_spec_get_instance_private (gsdorspec);
  guint max_tries                               = 1000;

  do {
    data->z = nc_galaxy_sd_true_redshift_gen (self->sdz, rng);

    if (max_tries-- == 0)
      g_error ("nc_galaxy_sd_obs_redshift_spec_gen: maximum number of iterations reached.");
  } while (data->z < self->z_min || data->z > self->z_max);
}

static gboolean
_nc_galaxy_sd_obs_redshift_spec_gen1 (NcGalaxySDObsRedshift *gsdor, NcGalaxySDObsRedshiftData *data, NcmRNG *rng)
{
  NcGalaxySDObsRedshiftSpec *gsdorspec          = NC_GALAXY_SD_OBS_REDSHIFT_SPEC (gsdor);
  NcGalaxySDObsRedshiftSpecPrivate * const self = nc_galaxy_sd_obs_redshift_spec_get_instance_private (gsdorspec);

  data->z = nc_galaxy_sd_true_redshift_gen (self->sdz, rng);

  return (data->z >= self->z_min) && (data->z <= self->z_max);
}

static void
_nc_galaxy_sd_obs_redshift_spec_prepare (NcGalaxySDObsRedshift *gsdor, NcGalaxySDObsRedshiftData *data)
{
  /* Nothing to do */
}

static void
_nc_galaxy_sd_obs_redshift_spec_get_lim (NcGalaxySDObsRedshift *gsdor, NcGalaxySDObsRedshiftData *data, gdouble *z_min, gdouble *z_max)
{
  NcGalaxySDObsRedshiftSpec *gsdorspec          = NC_GALAXY_SD_OBS_REDSHIFT_SPEC (gsdor);
  NcGalaxySDObsRedshiftSpecPrivate * const self = nc_galaxy_sd_obs_redshift_spec_get_instance_private (gsdorspec);

  *z_min = self->z_min;
  *z_max = self->z_max;
}

struct _IntegData
{
  NcGalaxySDObsRedshiftSpec *gsdorspec;
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
_nc_galaxy_sd_obs_redshift_spec_integ_f (gpointer callback_data, const gdouble z, NcGalaxySDObsRedshiftData *data)
{
  const struct _IntegData *int_data             = (struct _IntegData *) callback_data;
  NcGalaxySDObsRedshiftSpecPrivate * const self = nc_galaxy_sd_obs_redshift_spec_get_instance_private (int_data->gsdorspec);

  return nc_galaxy_sd_true_redshift_integ (self->sdz, z);
}

static NcGalaxySDObsRedshiftIntegrand *
_nc_galaxy_sd_obs_redshift_spec_integ (NcGalaxySDObsRedshift *gsdor)
{
  NcGalaxySDObsRedshiftSpec *gsdorspec  = NC_GALAXY_SD_OBS_REDSHIFT_SPEC (gsdor);
  struct _IntegData *int_data           = g_new0 (struct _IntegData, 1);
  NcGalaxySDObsRedshiftIntegrand *integ = nc_galaxy_sd_obs_redshift_integrand_new (_nc_galaxy_sd_obs_redshift_spec_integ_f,
                                                                                   _integ_data_free,
                                                                                   _integ_data_copy,
                                                                                   NULL,
                                                                                   int_data);

  int_data->gsdorspec = gsdorspec;

  return integ;
}

static void
_nc_galaxy_sd_obs_redshift_spec_ldata_free (gpointer ldata)
{
  g_free (ldata);
}

static void
_nc_galaxy_sd_obs_redshift_spec_ldata_read_row (NcGalaxySDObsRedshiftData *data, NcGalaxyWLObs *obs, const guint i)
{
  /* Nothing to do */
}

static void
_nc_galaxy_sd_obs_redshift_spec_ldata_write_row (NcGalaxySDObsRedshiftData *data, NcGalaxyWLObs *obs, const guint i)
{
  /* Nothing to do */
}

static void
_nc_galaxy_sd_obs_redshift_spec_ldata_required_columns (NcGalaxySDObsRedshiftData *data, GList *columns)
{
  /* Nothing to do */
}

static void
_nc_galaxy_sd_obs_redshift_spec_data_init (NcGalaxySDObsRedshift *gsdor, NcGalaxySDObsRedshiftData *data)
{
  NcGalaxySDObsRedshiftSpecData *ldata = g_new0 (NcGalaxySDObsRedshiftSpecData, 1);

  data->ldata                  = ldata;
  data->ldata_destroy          = &_nc_galaxy_sd_obs_redshift_spec_ldata_free;
  data->ldata_read_row         = &_nc_galaxy_sd_obs_redshift_spec_ldata_read_row;
  data->ldata_write_row        = &_nc_galaxy_sd_obs_redshift_spec_ldata_write_row;
  data->ldata_required_columns = &_nc_galaxy_sd_obs_redshift_spec_ldata_required_columns;
}

static void
_nc_galaxy_sd_obs_redshift_spec_add_submodel (NcmModel *model, NcmModel *submodel)
{
  /* Chain up: start */
  NCM_MODEL_CLASS (nc_galaxy_sd_obs_redshift_spec_parent_class)->add_submodel (model, submodel);
  {
    NcGalaxySDObsRedshiftSpec *gsdorspec          = NC_GALAXY_SD_OBS_REDSHIFT_SPEC (model);
    NcGalaxySDObsRedshiftSpecPrivate * const self = nc_galaxy_sd_obs_redshift_spec_get_instance_private (gsdorspec);

    g_assert (ncm_model_is_submodel (submodel));
    g_assert (NC_IS_GALAXY_SD_TRUE_REDSHIFT (submodel));

    self->sdz = NC_GALAXY_SD_TRUE_REDSHIFT (submodel);
  }
}

/**
 * nc_galaxy_sd_obs_redshift_spec_new:
 * @sdz: a #NcGalaxySDTrueRedshift
 * @z_min: the minimum redshift
 * @z_max: the maximum redshift
 *
 * Creates a new #NcGalaxySDObsRedshiftSpec object.
 *
 * Returns: (transfer full): a new #NcGalaxySDObsRedshiftSpec object.
 */
NcGalaxySDObsRedshiftSpec *
nc_galaxy_sd_obs_redshift_spec_new (NcGalaxySDTrueRedshift *sdz, const gdouble z_min, const gdouble z_max)
{
  NcmDTuple2 lim                       = NCM_DTUPLE2_STATIC_INIT (z_min, z_max);
  NcGalaxySDObsRedshiftSpec *gsdorspec = g_object_new (NC_TYPE_GALAXY_SD_OBS_REDSHIFT_SPEC,
                                                       "lim", &lim,
                                                       NULL);

  ncm_model_add_submodel (NCM_MODEL (gsdorspec), NCM_MODEL (sdz));

  return gsdorspec;
}

/**
 * nc_galaxy_sd_obs_redshift_spec_ref:
 * @gsdorspec: a #NcGalaxySDObsRedshiftSpec
 *
 * Increases the reference count of @gsdorspec by one.
 *
 * Returns: (transfer full): @gsdorspec.
 */
NcGalaxySDObsRedshiftSpec *
nc_galaxy_sd_obs_redshift_spec_ref (NcGalaxySDObsRedshiftSpec *gsdorspec)
{
  return g_object_ref (gsdorspec);
}

/**
 * nc_galaxy_sd_obs_redshift_spec_free:
 * @gsdorspec: a #NcGalaxySDObsRedshiftSpec
 *
 * Decreases the reference count of @gsdorspec by one.
 *
 */
void
nc_galaxy_sd_obs_redshift_spec_free (NcGalaxySDObsRedshiftSpec *gsdorspec)
{
  g_object_unref (gsdorspec);
}

/**
 * nc_galaxy_sd_obs_redshift_spec_clear:
 * @gsdorspec: a #NcGalaxySDObsRedshiftSpec
 *
 * Decreases the reference count of @gsdorspec by one, and sets the pointer *@gsdorspec to
 * NULL.
 *
 */
void
nc_galaxy_sd_obs_redshift_spec_clear (NcGalaxySDObsRedshiftSpec **gsdorspec)
{
  g_clear_object (gsdorspec);
}

/**
 * nc_galaxy_sd_obs_redshift_spec_set_lim:
 * @gsdorspec: a #NcGalaxySDObsRedshiftSpec
 * @z_min: the minimum redshift
 * @z_max: the maximum redshift
 *
 * Sets the redshift limits of the galaxy sample redshift distribution.
 *
 */
void
nc_galaxy_sd_obs_redshift_spec_set_lim (NcGalaxySDObsRedshiftSpec *gsdorspec, const gdouble z_min, const gdouble z_max)
{
  NcGalaxySDObsRedshiftSpecPrivate * const self = nc_galaxy_sd_obs_redshift_spec_get_instance_private (gsdorspec);

  g_assert_cmpfloat (z_min, <, z_max);

  self->z_min = z_min;
  self->z_max = z_max;
}

/**
 * nc_galaxy_sd_obs_redshift_spec_get_lim:
 * @gsdorspec: a #NcGalaxySDObsRedshiftSpec
 * @z_min: (out): the minimum redshift
 * @z_max: (out): the maximum redshift
 *
 * Returns the redshift limits of the galaxy sample redshift distribution.
 */
void
nc_galaxy_sd_obs_redshift_spec_get_lim (NcGalaxySDObsRedshiftSpec *gsdorspec, gdouble *z_min, gdouble *z_max)
{
  NcGalaxySDObsRedshiftSpecPrivate * const self = nc_galaxy_sd_obs_redshift_spec_get_instance_private (gsdorspec);

  *z_min = self->z_min;
  *z_max = self->z_max;
}

/**
 * nc_galaxy_sd_obs_redshift_spec_gen:
 * @gsdorspec: a #NcGalaxySDObsRedshiftSpec
 * @mset: a #NcmMSet
 * @data: a #NcGalaxySDObsRedshiftData
 * @rng: a #NcmRNG
 *
 * Generates a galaxy observed redshift.
 *
 */
void
nc_galaxy_sd_obs_redshift_spec_gen (NcGalaxySDObsRedshiftSpec *gsdorspec, NcmMSet *mset, NcGalaxySDObsRedshiftData *data, NcmRNG *rng)
{
  NcGalaxySDObsRedshiftClass *klass = NC_GALAXY_SD_OBS_REDSHIFT_GET_CLASS (gsdorspec);

  klass->gen (NC_GALAXY_SD_OBS_REDSHIFT (gsdorspec), data, rng);
}

/**
 * nc_galaxy_sd_obs_redshift_spec_gen1:
 * @gsdorspec: a #NcGalaxySDObsRedshiftSpec
 * @mset: a #NcmMSet
 * @data: a #NcGalaxySDObsRedshiftData
 * @rng: a #NcmRNG
 *
 * Generates a galaxy observed redshift. See nc_galaxy_sd_obs_redshift_spec_gen() for
 * details.
 *
 * Returns: %TRUE if @data->z is within the redshift limits
 */
gboolean
nc_galaxy_sd_obs_redshift_spec_gen1 (NcGalaxySDObsRedshiftSpec *gsdorspec, NcmMSet *mset, NcGalaxySDObsRedshiftData *data, NcmRNG *rng)
{
  NcGalaxySDObsRedshiftClass *klass = NC_GALAXY_SD_OBS_REDSHIFT_GET_CLASS (gsdorspec);

  return klass->gen1 (NC_GALAXY_SD_OBS_REDSHIFT (gsdorspec), data, rng);
}


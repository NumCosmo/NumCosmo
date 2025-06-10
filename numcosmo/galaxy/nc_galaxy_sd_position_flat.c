/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 2; tab-width: 2 -*-  */

/***************************************************************************
 *            nc_galaxy_sd_position_flat.c
 *
 *  Wed March 1 12:51:22 2023
 *  Copyright  2023  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 *  Copyright  2024  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_galaxy_sd_position_flat.c
 * Copyright (C) 2023 Caio Lima de Oliveira <caiolimadeoliveira@pm.me>
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
 * NcGalaxySDPositionFlat:
 *
 * Class describing galaxy sample position distributions with flat distribution.
 *
 * Class defining a galaxy sample position distribution with flat probability
 * distribution $P(\text{ra})P(\text{dec})$.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "galaxy/nc_galaxy_sd_position_flat.h"
#include "galaxy/nc_galaxy_sd_position.h"
#include "math/ncm_dtuple.h"
#include "math/ncm_model.h"
#include "math/ncm_vector.h"
#include "math/ncm_rng.h"
#include "math/ncm_c.h"

typedef struct _NcGalaxySDPositionFlatPrivate
{
  gdouble ra_min;
  gdouble ra_max;
  gdouble ra_norm;
  gdouble dec_min;
  gdouble dec_max;
  gdouble sin_dec_min;
  gdouble sin_dec_max;
  gdouble dec_norm;
} NcGalaxySDPositionFlatPrivate;

struct _NcGalaxySDPositionFlat
{
  NcGalaxySDPosition parent_instance;
};

typedef struct _NcGalaxySDPositionFlatData
{
  gint placeholder;
} NcGalaxySDPositionFlatData;

enum
{
  PROP_0,
  PROP_LEN,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcGalaxySDPositionFlat, nc_galaxy_sd_position_flat, NC_TYPE_GALAXY_SD_POSITION);

static void
nc_galaxy_sd_position_flat_init (NcGalaxySDPositionFlat *gsdpflat)
{
  NcGalaxySDPositionFlatPrivate * const self = nc_galaxy_sd_position_flat_get_instance_private (gsdpflat);

  self->ra_min      = 0.0;
  self->ra_max      = 0.0;
  self->ra_norm     = 0.0;
  self->dec_min     = 0.0;
  self->dec_max     = 0.0;
  self->sin_dec_min = 0.0;
  self->sin_dec_max = 0.0;
  self->dec_norm    = 0.0;
}

static void
_nc_galaxy_sd_position_flat_dispose (GObject *object)
{
  /* Chain up: end */
  G_OBJECT_CLASS (nc_galaxy_sd_position_flat_parent_class)->dispose (object);
}

static void
_nc_galaxy_sd_position_flat_finalize (GObject *object)
{
  /* Chain up: end */
  G_OBJECT_CLASS (nc_galaxy_sd_position_flat_parent_class)->finalize (object);
}

static void _nc_galaxy_sd_position_flat_gen (NcGalaxySDPosition *gsdp, NcGalaxySDPositionData *data, NcmRNG *rng);
static NcGalaxySDPositionIntegrand *_nc_galaxy_sd_position_flat_integ (NcGalaxySDPosition *gsdp);
static gboolean _nc_galaxy_sd_position_flat_set_ra_lim (NcGalaxySDPosition *gsdp, gdouble ra_min, gdouble ra_max);
static gboolean _nc_galaxy_sd_position_flat_get_ra_lim (NcGalaxySDPosition *gsdp, gdouble *ra_min, gdouble *ra_max);
static gboolean _nc_galaxy_sd_position_flat_set_dec_lim (NcGalaxySDPosition *gsdp, gdouble dec_min, gdouble dec_max);
static gboolean _nc_galaxy_sd_position_flat_get_dec_lim (NcGalaxySDPosition *gsdp, gdouble *dec_min, gdouble *dec_max);
static void _nc_galaxy_sd_position_flat_data_init (NcGalaxySDPosition *gsdp, NcGalaxySDObsRedshiftData *sdz_data, NcGalaxySDPositionData *data);

static void
nc_galaxy_sd_position_flat_class_init (NcGalaxySDPositionFlatClass *klass)
{
  NcGalaxySDPositionClass *sd_position_class = NC_GALAXY_SD_POSITION_CLASS (klass);
  GObjectClass *object_class                 = G_OBJECT_CLASS (klass);
  NcmModelClass *model_class                 = NCM_MODEL_CLASS (klass);

  object_class->dispose  = &_nc_galaxy_sd_position_flat_dispose;
  object_class->finalize = &_nc_galaxy_sd_position_flat_finalize;

  ncm_model_class_set_name_nick (model_class, "Flat Galaxy Distribution", "FLAT_GALAXY_DISTRIBUTION");
  ncm_model_class_add_params (model_class, 0, 0, PROP_LEN);

  ncm_model_class_check_params_info (model_class);

  sd_position_class->gen         = &_nc_galaxy_sd_position_flat_gen;
  sd_position_class->integ       = &_nc_galaxy_sd_position_flat_integ;
  sd_position_class->set_ra_lim  = &_nc_galaxy_sd_position_flat_set_ra_lim;
  sd_position_class->get_ra_lim  = &_nc_galaxy_sd_position_flat_get_ra_lim;
  sd_position_class->set_dec_lim = &_nc_galaxy_sd_position_flat_set_dec_lim;
  sd_position_class->get_dec_lim = &_nc_galaxy_sd_position_flat_get_dec_lim;
  sd_position_class->data_init   = &_nc_galaxy_sd_position_flat_data_init;
}

static void
_nc_galaxy_sd_position_flat_gen (NcGalaxySDPosition *gsdp, NcGalaxySDPositionData *data, NcmRNG *rng)
{
  NcGalaxySDPositionFlat *gsdpflat           = NC_GALAXY_SD_POSITION_FLAT (gsdp);
  NcGalaxySDPositionFlatPrivate * const self = nc_galaxy_sd_position_flat_get_instance_private (gsdpflat);
  gdouble sin_dec                            = ncm_rng_uniform_gen (rng, self->sin_dec_min, self->sin_dec_max);

  data->ra  = ncm_rng_uniform_gen (rng, self->ra_min, self->ra_max);
  data->dec = ncm_c_radian_to_degree (asin (sin_dec));
}

struct _IntegData
{
  NcGalaxySDPositionFlat *gsdpflat;
  NcGalaxySDPositionData *data;
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
_nc_galaxy_sd_position_flat_integ_f (gpointer callback_data, NcGalaxySDPositionData *data)
{
  const struct _IntegData *int_data          = (struct _IntegData *) callback_data;
  NcGalaxySDPositionFlatPrivate * const self = nc_galaxy_sd_position_flat_get_instance_private (int_data->gsdpflat);
  const gdouble ra_norm                      = self->ra_norm;
  const gdouble dec_norm                     = self->dec_norm;
  const gdouble ra                           = data->ra;
  const gdouble dec                          = data->dec;

  if ((ra >= self->ra_min) && (ra <= self->ra_max) && (dec >= self->dec_min) && (dec <= self->dec_max))
    return ra_norm * dec_norm * cos (ncm_c_degree_to_radian (dec));

  return 0.0;
}

static NcGalaxySDPositionIntegrand *
_nc_galaxy_sd_position_flat_integ (NcGalaxySDPosition *gsdp)
{
  NcGalaxySDPositionFlat *gsdpflat   = NC_GALAXY_SD_POSITION_FLAT (gsdp);
  struct _IntegData *int_data        = g_new0 (struct _IntegData, 1);
  NcGalaxySDPositionIntegrand *integ = nc_galaxy_sd_position_integrand_new (_nc_galaxy_sd_position_flat_integ_f,
                                                                            _integ_data_free,
                                                                            _integ_data_copy,
                                                                            NULL,
                                                                            int_data);

  int_data->gsdpflat = gsdpflat;

  return integ;
}

static gboolean
_nc_galaxy_sd_position_flat_set_ra_lim (NcGalaxySDPosition *gsdp, gdouble ra_min, gdouble ra_max)
{
  NcGalaxySDPositionFlat *gsdpflat           = NC_GALAXY_SD_POSITION_FLAT (gsdp);
  NcGalaxySDPositionFlatPrivate * const self = nc_galaxy_sd_position_flat_get_instance_private (gsdpflat);

  g_assert_cmpfloat (ra_min, <, ra_max);

  self->ra_min  = ra_min;
  self->ra_max  = ra_max;
  self->ra_norm = 1.0 / (self->ra_max - self->ra_min);

  return TRUE;
}

static gboolean
_nc_galaxy_sd_position_flat_get_ra_lim (NcGalaxySDPosition *gsdp, gdouble *ra_min, gdouble *ra_max)
{
  NcGalaxySDPositionFlat *gsdpflat           = NC_GALAXY_SD_POSITION_FLAT (gsdp);
  NcGalaxySDPositionFlatPrivate * const self = nc_galaxy_sd_position_flat_get_instance_private (gsdpflat);

  g_assert_nonnull (ra_min);
  g_assert_nonnull (ra_max);

  *ra_min = self->ra_min;
  *ra_max = self->ra_max;

  return TRUE;
}

static gboolean
_nc_galaxy_sd_position_flat_set_dec_lim (NcGalaxySDPosition *gsdp, gdouble dec_min, gdouble dec_max)
{
  NcGalaxySDPositionFlat *gsdpflat           = NC_GALAXY_SD_POSITION_FLAT (gsdp);
  NcGalaxySDPositionFlatPrivate * const self = nc_galaxy_sd_position_flat_get_instance_private (gsdpflat);

  g_assert_cmpfloat (dec_min, <, dec_max);

  self->dec_min     = dec_min;
  self->dec_max     = dec_max;
  self->sin_dec_min = sin (ncm_c_degree_to_radian (dec_min));
  self->sin_dec_max = sin (ncm_c_degree_to_radian (dec_max));
  self->dec_norm    = ncm_c_pi () / (180.0 * (self->sin_dec_max - self->sin_dec_min));

  return TRUE;
}

static gboolean
_nc_galaxy_sd_position_flat_get_dec_lim (NcGalaxySDPosition *gsdp, gdouble *dec_min, gdouble *dec_max)
{
  NcGalaxySDPositionFlat *gsdpflat           = NC_GALAXY_SD_POSITION_FLAT (gsdp);
  NcGalaxySDPositionFlatPrivate * const self = nc_galaxy_sd_position_flat_get_instance_private (gsdpflat);

  g_assert_nonnull (dec_min);
  g_assert_nonnull (dec_max);

  *dec_min = self->dec_min;
  *dec_max = self->dec_max;

  return TRUE;
}

static void
_nc_galaxy_sd_position_flat_ldata_free (gpointer ldata)
{
  g_free (ldata);
}

static void
_nc_galaxy_sd_position_flat_ldata_read_row (NcGalaxySDPositionData *data, NcGalaxyWLObs *obs, const guint i)
{
  /* Nothing to do */
}

static void
_nc_galaxy_sd_position_flat_ldata_write_row (NcGalaxySDPositionData *data, NcGalaxyWLObs *obs, const guint i)
{
  /* Nothing to do */
}

static void
_nc_galaxy_sd_position_flat_ldata_required_columns (NcGalaxySDPositionData *data, GList *columns)
{
  /* Nothing to do */
}

static void
_nc_galaxy_sd_position_flat_data_init (NcGalaxySDPosition *gsdp, NcGalaxySDObsRedshiftData *sdz_data, NcGalaxySDPositionData *data)
{
  NcGalaxySDPositionFlatData *ldata = g_new0 (NcGalaxySDPositionFlatData, 1);

  data->sdz_data               = sdz_data;
  data->ra                     = 0.0;
  data->dec                    = 0.0;
  data->ldata                  = ldata;
  data->ldata_destroy          = &_nc_galaxy_sd_position_flat_ldata_free;
  data->ldata_read_row         = &_nc_galaxy_sd_position_flat_ldata_read_row;
  data->ldata_write_row        = &_nc_galaxy_sd_position_flat_ldata_write_row;
  data->ldata_required_columns = &_nc_galaxy_sd_position_flat_ldata_required_columns;
}

/**
 * nc_galaxy_sd_position_flat_new:
 * @ra_min: minimum right ascension
 * @ra_max: maximum right ascension
 * @dec_min: minimum declination
 * @dec_max: maximum declination
 *
 * Creates a new #NcGalaxySDPositionFlat
 *
 * Returns: (transfer full): a new NcGalaxySDPositionFlat.
 */
NcGalaxySDPositionFlat *
nc_galaxy_sd_position_flat_new (const gdouble ra_min, const gdouble ra_max, const gdouble dec_min, const gdouble dec_max)
{
  NcmDTuple2 ra_lim  = NCM_DTUPLE2_STATIC_INIT (ra_min, ra_max);
  NcmDTuple2 dec_lim = NCM_DTUPLE2_STATIC_INIT (dec_min, dec_max);
  NcGalaxySDPositionFlat *gsdpflat;

  g_assert_cmpfloat (ra_min, <, ra_max);
  g_assert_cmpfloat (dec_min, <, dec_max);

  gsdpflat = g_object_new (NC_TYPE_GALAXY_SD_POSITION_FLAT,
                           "ra-lim", &ra_lim,
                           "dec-lim", &dec_lim,
                           NULL);

  return gsdpflat;
}

/**
 * nc_galaxy_sd_position_flat_ref:
 * @gsdpflat: a #NcGalaxySDPositionFlat
 *
 * Increases the reference count of @gsdpflat by one.
 *
 * Returns: (transfer full): @gsdpflat.
 */
NcGalaxySDPositionFlat *
nc_galaxy_sd_position_flat_ref (NcGalaxySDPositionFlat *gsdpflat)
{
  return g_object_ref (gsdpflat);
}

/**
 * nc_galaxy_sd_position_flat_free:
 * @gsdpflat: a #NcGalaxySDPositionFlat
 *
 * Decrease sthe reference count of @gsdpflat by one.
 *
 */
void
nc_galaxy_sd_position_flat_free (NcGalaxySDPositionFlat *gsdpflat)
{
  g_object_unref (gsdpflat);
}

/**
 * nc_galaxy_sd_position_flat_clear:
 * @gsdpflat: a #NcGalaxySDPositionFlat
 *
 * Decrease the reference count of @gsdpflat by one, and sets the pointer *@gsdpflat to
 * NULL.
 *
 */
void
nc_galaxy_sd_position_flat_clear (NcGalaxySDPositionFlat **gsdpflat)
{
  g_clear_object (gsdpflat);
}

/**
 * nc_galaxy_sd_position_flat_gen:
 * @gsdpflat: a #NcGalaxySDPositionFlat
 * @mset: a #NcmMSet
 * @data: a #NcGalaxySDPositionData
 * @rng: a #NcmRNG
 *
 * Generates a galaxy sample position.
 *
 */
void
nc_galaxy_sd_position_flat_gen (NcGalaxySDPositionFlat *gsdpflat, NcmMSet *mset, NcGalaxySDPositionData *data, NcmRNG *rng)
{
  NcGalaxySDPositionClass *sd_position_class = NC_GALAXY_SD_POSITION_GET_CLASS (gsdpflat);

  sd_position_class->gen (NC_GALAXY_SD_POSITION (gsdpflat), data, rng);
}


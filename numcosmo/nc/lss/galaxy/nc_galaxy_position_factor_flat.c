/***************************************************************************
 *            nc_galaxy_position_factor_flat.c
 *
 *  Wed Jul 2 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_galaxy_position_factor_flat.c
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
 * with this program. If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * NcGalaxyPositionFactorFlat:
 *
 * Flat (uniform) galaxy sky-position scheme.
 *
 * A #NcGalaxyPositionFactor scheme whose density $P(\mathrm{ra},
 * \mathrm{dec})$ is uniform over a rectangular sky footprint. It carries no
 * fitted parameters: the footprint (survey geometry) is held as configuration
 * via nc_galaxy_position_factor_flat_set_ra_lim() /
 * nc_galaxy_position_factor_flat_set_dec_lim(), not as an #NcmModel in an
 * #NcmMSet.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc/lss/galaxy/nc_galaxy_position_factor_flat.h"
#include "nc/lss/galaxy/nc_galaxy_position_factor.h"
#include "ncm/core/ncm_dtuple.h"
#include "ncm/core/ncm_rng.h"
#include "ncm/sphere/ncm_sky_footprint.h"
#include "ncm/sphere/ncm_sky_footprint_rectangular.h"

typedef struct _NcGalaxyPositionFactorFlatPrivate
{
  NcmSkyFootprintRectangular *footprint;
} NcGalaxyPositionFactorFlatPrivate;

struct _NcGalaxyPositionFactorFlat
{
  NcGalaxyPositionFactor parent_instance;
};

typedef struct _NcGalaxyPositionFactorFlatData
{
  gint placeholder;
} NcGalaxyPositionFactorFlatData;

enum
{
  PROP_0,
  PROP_RA_LIM,
  PROP_DEC_LIM,
  PROP_LEN,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcGalaxyPositionFactorFlat, nc_galaxy_position_factor_flat, NC_TYPE_GALAXY_POSITION_FACTOR);

static void
nc_galaxy_position_factor_flat_init (NcGalaxyPositionFactorFlat *gspfflat)
{
  NcGalaxyPositionFactorFlatPrivate * const self = nc_galaxy_position_factor_flat_get_instance_private (gspfflat);

  self->footprint = g_object_new (NCM_TYPE_SKY_FOOTPRINT_RECTANGULAR, NULL);
}

static void
_nc_galaxy_position_factor_flat_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcGalaxyPositionFactorFlat *gspfflat = NC_GALAXY_POSITION_FACTOR_FLAT (object);

  g_return_if_fail (NC_IS_GALAXY_POSITION_FACTOR_FLAT (object));

  switch (prop_id)
  {
    case PROP_RA_LIM:
    {
      NcmDTuple2 *ra_lim = g_value_get_boxed (value);

      if (ra_lim == NULL)
        g_error ("_nc_galaxy_position_factor_flat_set_property: ra_lim is NULL.");

      nc_galaxy_position_factor_flat_set_ra_lim (gspfflat, ra_lim->elements[0], ra_lim->elements[1]);
      break;
    }
    case PROP_DEC_LIM:
    {
      NcmDTuple2 *dec_lim = g_value_get_boxed (value);

      if (dec_lim == NULL)
        g_error ("_nc_galaxy_position_factor_flat_set_property: dec_lim is NULL.");

      nc_galaxy_position_factor_flat_set_dec_lim (gspfflat, dec_lim->elements[0], dec_lim->elements[1]);
      break;
    }
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_galaxy_position_factor_flat_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcGalaxyPositionFactorFlat *gspfflat = NC_GALAXY_POSITION_FACTOR_FLAT (object);

  g_return_if_fail (NC_IS_GALAXY_POSITION_FACTOR_FLAT (object));

  switch (prop_id)
  {
    case PROP_RA_LIM:
    {
      gdouble ra_min, ra_max;

      nc_galaxy_position_factor_flat_get_ra_lim (gspfflat, &ra_min, &ra_max);

      g_value_take_boxed (value, ncm_dtuple2_new (ra_min, ra_max));
      break;
    }
    case PROP_DEC_LIM:
    {
      gdouble dec_min, dec_max;

      nc_galaxy_position_factor_flat_get_dec_lim (gspfflat, &dec_min, &dec_max);

      g_value_take_boxed (value, ncm_dtuple2_new (dec_min, dec_max));
      break;
    }
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_galaxy_position_factor_flat_dispose (GObject *object)
{
  NcGalaxyPositionFactorFlat *gspfflat           = NC_GALAXY_POSITION_FACTOR_FLAT (object);
  NcGalaxyPositionFactorFlatPrivate * const self = nc_galaxy_position_factor_flat_get_instance_private (gspfflat);

  ncm_sky_footprint_rectangular_clear (&self->footprint);

  /* Chain up: end */
  G_OBJECT_CLASS (nc_galaxy_position_factor_flat_parent_class)->dispose (object);
}

static void
_nc_galaxy_position_factor_flat_finalize (GObject *object)
{
  /* Chain up: end */
  G_OBJECT_CLASS (nc_galaxy_position_factor_flat_parent_class)->finalize (object);
}

static void _nc_galaxy_position_factor_flat_data_init (NcGalaxyPositionFactor *gspf, NcmMSet *mset, NcGalaxyPositionFactorData *data);
static void _nc_galaxy_position_factor_flat_gen (NcGalaxyPositionFactor *gspf, NcmMSet *mset, NcGalaxyPositionFactorData *data, NcmRNG *rng);
static void _nc_galaxy_position_factor_flat_prepare (NcGalaxyPositionFactor *gspf, NcmMSet *mset, NcGalaxyPositionFactorData *data);
static NcGalaxyPositionFactorIntegrand *_nc_galaxy_position_factor_flat_integ (NcGalaxyPositionFactor *gspf, NcmMSet *mset, gboolean use_lnp);

static void
nc_galaxy_position_factor_flat_class_init (NcGalaxyPositionFactorFlatClass *klass)
{
  GObjectClass *object_class                   = G_OBJECT_CLASS (klass);
  NcGalaxyPositionFactorClass *position_factor_class = NC_GALAXY_POSITION_FACTOR_CLASS (klass);

  object_class->set_property = &_nc_galaxy_position_factor_flat_set_property;
  object_class->get_property = &_nc_galaxy_position_factor_flat_get_property;
  object_class->dispose      = &_nc_galaxy_position_factor_flat_dispose;
  object_class->finalize     = &_nc_galaxy_position_factor_flat_finalize;

  /**
   * NcGalaxyPositionFactorFlat:ra-lim:
   *
   * Galaxy sample right ascension footprint limits.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_RA_LIM,
                                   g_param_spec_boxed ("ra-lim",
                                                       NULL,
                                                       "Galaxy sample right ascension footprint limits",
                                                       NCM_TYPE_DTUPLE2,
                                                       G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcGalaxyPositionFactorFlat:dec-lim:
   *
   * Galaxy sample declination footprint limits.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_DEC_LIM,
                                   g_param_spec_boxed ("dec-lim",
                                                       NULL,
                                                       "Galaxy sample declination footprint limits",
                                                       NCM_TYPE_DTUPLE2,
                                                       G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  position_factor_class->data_init = &_nc_galaxy_position_factor_flat_data_init;
  position_factor_class->gen       = &_nc_galaxy_position_factor_flat_gen;
  position_factor_class->prepare   = &_nc_galaxy_position_factor_flat_prepare;
  position_factor_class->integ     = &_nc_galaxy_position_factor_flat_integ;
}

static void
_nc_galaxy_position_factor_flat_ldata_free (gpointer ldata)
{
  g_free (ldata);
}

static void
_nc_galaxy_position_factor_flat_ldata_read_row (NcGalaxyPositionFactorData *data, NcGalaxyWLObs *obs, const guint i)
{
  /* Nothing to do */
}

static void
_nc_galaxy_position_factor_flat_ldata_write_row (NcGalaxyPositionFactorData *data, NcGalaxyWLObs *obs, const guint i)
{
  /* Nothing to do */
}

static void
_nc_galaxy_position_factor_flat_ldata_required_columns (NcGalaxyPositionFactorData *data, GList **columns)
{
  /* Nothing to do */
}

static void
_nc_galaxy_position_factor_flat_data_init (NcGalaxyPositionFactor *gspf, NcmMSet *mset, NcGalaxyPositionFactorData *data)
{
  NcGalaxyPositionFactorFlatData *ldata = g_new0 (NcGalaxyPositionFactorFlatData, 1);

  data->ldata                  = ldata;
  data->ldata_destroy          = &_nc_galaxy_position_factor_flat_ldata_free;
  data->ldata_read_row         = &_nc_galaxy_position_factor_flat_ldata_read_row;
  data->ldata_write_row        = &_nc_galaxy_position_factor_flat_ldata_write_row;
  data->ldata_required_columns = &_nc_galaxy_position_factor_flat_ldata_required_columns;
}

static void
_nc_galaxy_position_factor_flat_gen (NcGalaxyPositionFactor *gspf, NcmMSet *mset, NcGalaxyPositionFactorData *data, NcmRNG *rng)
{
  NcGalaxyPositionFactorFlat *gspfflat           = NC_GALAXY_POSITION_FACTOR_FLAT (gspf);
  NcGalaxyPositionFactorFlatPrivate * const self = nc_galaxy_position_factor_flat_get_instance_private (gspfflat);

  ncm_sky_footprint_gen_ra_dec (NCM_SKY_FOOTPRINT (self->footprint), rng, &data->ra, &data->dec);
}

static void
_nc_galaxy_position_factor_flat_prepare (NcGalaxyPositionFactor *gspf, NcmMSet *mset, NcGalaxyPositionFactorData *data)
{
  /* No models to resolve: the flat scheme's footprint is held configuration. */
}

struct _IntegData
{
  NcGalaxyPositionFactorFlat *gspfflat;
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
_nc_galaxy_position_factor_flat_ln_integ_f (gpointer callback_data, NcGalaxyPositionFactorData *data)
{
  const struct _IntegData *int_data              = (struct _IntegData *) callback_data;
  NcGalaxyPositionFactorFlatPrivate * const self = nc_galaxy_position_factor_flat_get_instance_private (int_data->gspfflat);

  return ncm_sky_footprint_ln_density (NCM_SKY_FOOTPRINT (self->footprint), data->ra, data->dec);
}

static gdouble
_nc_galaxy_position_factor_flat_integ_f (gpointer callback_data, NcGalaxyPositionFactorData *data)
{
  const struct _IntegData *int_data              = (struct _IntegData *) callback_data;
  NcGalaxyPositionFactorFlatPrivate * const self = nc_galaxy_position_factor_flat_get_instance_private (int_data->gspfflat);

  return ncm_sky_footprint_density (NCM_SKY_FOOTPRINT (self->footprint), data->ra, data->dec);
}

static NcGalaxyPositionFactorIntegrand *
_nc_galaxy_position_factor_flat_integ (NcGalaxyPositionFactor *gspf, NcmMSet *mset, gboolean use_lnp)
{
  NcGalaxyPositionFactorFlat *gspfflat        = NC_GALAXY_POSITION_FACTOR_FLAT (gspf);
  struct _IntegData *int_data                 = g_new0 (struct _IntegData, 1);
  NcGalaxyPositionFactorIntegrand *integ      = nc_galaxy_position_factor_integrand_new (use_lnp ? _nc_galaxy_position_factor_flat_ln_integ_f : _nc_galaxy_position_factor_flat_integ_f,
                                                                                         _integ_data_free,
                                                                                         _integ_data_copy,
                                                                                         NULL,
                                                                                         int_data);

  int_data->gspfflat = gspfflat;

  return integ;
}

/**
 * nc_galaxy_position_factor_flat_new:
 * @ra_min: minimum right ascension
 * @ra_max: maximum right ascension
 * @dec_min: minimum declination
 * @dec_max: maximum declination
 *
 * Creates a new #NcGalaxyPositionFactorFlat uniform over the rectangular
 * footprint [@ra_min, @ra_max] x [@dec_min, @dec_max].
 *
 * Returns: (transfer full): a new #NcGalaxyPositionFactorFlat.
 */
NcGalaxyPositionFactorFlat *
nc_galaxy_position_factor_flat_new (const gdouble ra_min, const gdouble ra_max, const gdouble dec_min, const gdouble dec_max)
{
  NcmDTuple2 ra_lim  = NCM_DTUPLE2_STATIC_INIT (ra_min, ra_max);
  NcmDTuple2 dec_lim = NCM_DTUPLE2_STATIC_INIT (dec_min, dec_max);
  NcGalaxyPositionFactorFlat *gspfflat;

  g_assert_cmpfloat (ra_min, <, ra_max);
  g_assert_cmpfloat (dec_min, <, dec_max);

  gspfflat = g_object_new (NC_TYPE_GALAXY_POSITION_FACTOR_FLAT,
                           "ra-lim", &ra_lim,
                           "dec-lim", &dec_lim,
                           NULL);

  return gspfflat;
}

/**
 * nc_galaxy_position_factor_flat_ref:
 * @gspfflat: a #NcGalaxyPositionFactorFlat
 *
 * Increases the reference count of @gspfflat by one.
 *
 * Returns: (transfer full): @gspfflat.
 */
NcGalaxyPositionFactorFlat *
nc_galaxy_position_factor_flat_ref (NcGalaxyPositionFactorFlat *gspfflat)
{
  return g_object_ref (gspfflat);
}

/**
 * nc_galaxy_position_factor_flat_free:
 * @gspfflat: a #NcGalaxyPositionFactorFlat
 *
 * Decreases the reference count of @gspfflat by one.
 *
 */
void
nc_galaxy_position_factor_flat_free (NcGalaxyPositionFactorFlat *gspfflat)
{
  g_object_unref (gspfflat);
}

/**
 * nc_galaxy_position_factor_flat_clear:
 * @gspfflat: a #NcGalaxyPositionFactorFlat
 *
 * Decreases the reference count of @gspfflat by one, and sets the pointer
 * *@gspfflat to NULL.
 *
 */
void
nc_galaxy_position_factor_flat_clear (NcGalaxyPositionFactorFlat **gspfflat)
{
  g_clear_object (gspfflat);
}

/**
 * nc_galaxy_position_factor_flat_set_ra_lim:
 * @gspfflat: a #NcGalaxyPositionFactorFlat
 * @ra_min: the minimum right ascension
 * @ra_max: the maximum right ascension
 *
 * Sets the right ascension footprint limits.
 *
 */
void
nc_galaxy_position_factor_flat_set_ra_lim (NcGalaxyPositionFactorFlat *gspfflat, const gdouble ra_min, const gdouble ra_max)
{
  NcGalaxyPositionFactorFlatPrivate * const self = nc_galaxy_position_factor_flat_get_instance_private (gspfflat);

  ncm_sky_footprint_rectangular_set_ra_lim (self->footprint, ra_min, ra_max);
}

/**
 * nc_galaxy_position_factor_flat_get_ra_lim:
 * @gspfflat: a #NcGalaxyPositionFactorFlat
 * @ra_min: (out): the minimum right ascension
 * @ra_max: (out): the maximum right ascension
 *
 * Gets the right ascension footprint limits.
 *
 */
void
nc_galaxy_position_factor_flat_get_ra_lim (NcGalaxyPositionFactorFlat *gspfflat, gdouble *ra_min, gdouble *ra_max)
{
  NcGalaxyPositionFactorFlatPrivate * const self = nc_galaxy_position_factor_flat_get_instance_private (gspfflat);

  g_assert_nonnull (ra_min);
  g_assert_nonnull (ra_max);

  ncm_sky_footprint_rectangular_get_ra_lim (self->footprint, ra_min, ra_max);
}

/**
 * nc_galaxy_position_factor_flat_set_dec_lim:
 * @gspfflat: a #NcGalaxyPositionFactorFlat
 * @dec_min: the minimum declination
 * @dec_max: the maximum declination
 *
 * Sets the declination footprint limits.
 *
 */
void
nc_galaxy_position_factor_flat_set_dec_lim (NcGalaxyPositionFactorFlat *gspfflat, const gdouble dec_min, const gdouble dec_max)
{
  NcGalaxyPositionFactorFlatPrivate * const self = nc_galaxy_position_factor_flat_get_instance_private (gspfflat);

  ncm_sky_footprint_rectangular_set_dec_lim (self->footprint, dec_min, dec_max);
}

/**
 * nc_galaxy_position_factor_flat_get_dec_lim:
 * @gspfflat: a #NcGalaxyPositionFactorFlat
 * @dec_min: (out): the minimum declination
 * @dec_max: (out): the maximum declination
 *
 * Gets the declination footprint limits.
 *
 */
void
nc_galaxy_position_factor_flat_get_dec_lim (NcGalaxyPositionFactorFlat *gspfflat, gdouble *dec_min, gdouble *dec_max)
{
  NcGalaxyPositionFactorFlatPrivate * const self = nc_galaxy_position_factor_flat_get_instance_private (gspfflat);

  g_assert_nonnull (dec_min);
  g_assert_nonnull (dec_max);

  ncm_sky_footprint_rectangular_get_dec_lim (self->footprint, dec_min, dec_max);
}

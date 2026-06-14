/***************************************************************************
 *            ncm_sky_footprint_rectangular.c
 *
 *  Sat Jun 13 19:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_sky_footprint_rectangular.c
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
 * NcmSkyFootprintRectangular:
 *
 * Rectangular sky footprint in right ascension and declination.
 *
 * A #NcmSkyFootprint covering the rectangle
 * $[\mathrm{ra}_\mathrm{min}, \mathrm{ra}_\mathrm{max}] \times
 * [\mathrm{dec}_\mathrm{min}, \mathrm{dec}_\mathrm{max}]$ (degrees). Positions
 * are sampled uniformly on the sphere patch: right ascension uniform in its
 * range and declination uniform in $\sin(\mathrm{dec})$. The position density is
 * the corresponding normalized distribution, proportional to
 * $\cos(\mathrm{dec})$ inside the rectangle.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"
#include "math/ncm_sky_footprint_rectangular.h"
#include "math/ncm_dtuple.h"
#include "math/ncm_rng.h"
#include "math/ncm_c.h"

#include <math.h>
#include <gsl/gsl_math.h>

typedef struct _NcmSkyFootprintRectangularPrivate
{
  gdouble ra_min;
  gdouble ra_max;
  gdouble ra_norm;
  gdouble dec_min;
  gdouble dec_max;
  gdouble sin_dec_min;
  gdouble sin_dec_max;
  gdouble dec_norm;
} NcmSkyFootprintRectangularPrivate;

struct _NcmSkyFootprintRectangular
{
  NcmSkyFootprint parent_instance;
};

enum
{
  PROP_0,
  PROP_RA_LIM,
  PROP_DEC_LIM,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcmSkyFootprintRectangular, ncm_sky_footprint_rectangular, NCM_TYPE_SKY_FOOTPRINT)

static void
ncm_sky_footprint_rectangular_init (NcmSkyFootprintRectangular *rect)
{
  NcmSkyFootprintRectangularPrivate * const self = ncm_sky_footprint_rectangular_get_instance_private (rect);

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
_ncm_sky_footprint_rectangular_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmSkyFootprintRectangular *rect = NCM_SKY_FOOTPRINT_RECTANGULAR (object);

  switch (prop_id)
  {
    case PROP_RA_LIM:
    {
      NcmDTuple2 *ra_lim = g_value_get_boxed (value);

      if (ra_lim == NULL)
        g_error ("_ncm_sky_footprint_rectangular_set_property: ra_lim is NULL.");

      ncm_sky_footprint_rectangular_set_ra_lim (rect, ra_lim->elements[0], ra_lim->elements[1]);
      break;
    }
    case PROP_DEC_LIM:
    {
      NcmDTuple2 *dec_lim = g_value_get_boxed (value);

      if (dec_lim == NULL)
        g_error ("_ncm_sky_footprint_rectangular_set_property: dec_lim is NULL.");

      ncm_sky_footprint_rectangular_set_dec_lim (rect, dec_lim->elements[0], dec_lim->elements[1]);
      break;
    }
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_sky_footprint_rectangular_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmSkyFootprintRectangular *rect = NCM_SKY_FOOTPRINT_RECTANGULAR (object);

  switch (prop_id)
  {
    case PROP_RA_LIM:
    {
      gdouble ra_min, ra_max;

      ncm_sky_footprint_rectangular_get_ra_lim (rect, &ra_min, &ra_max);
      g_value_take_boxed (value, ncm_dtuple2_new (ra_min, ra_max));
      break;
    }
    case PROP_DEC_LIM:
    {
      gdouble dec_min, dec_max;

      ncm_sky_footprint_rectangular_get_dec_lim (rect, &dec_min, &dec_max);
      g_value_take_boxed (value, ncm_dtuple2_new (dec_min, dec_max));
      break;
    }
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void _ncm_sky_footprint_rectangular_gen_ra_dec (NcmSkyFootprint *footprint, NcmRNG *rng, gdouble *ra, gdouble *dec);
static gboolean _ncm_sky_footprint_rectangular_contains (NcmSkyFootprint *footprint, const gdouble ra, const gdouble dec);
static gdouble _ncm_sky_footprint_rectangular_get_area (NcmSkyFootprint *footprint);
static gdouble _ncm_sky_footprint_rectangular_density (NcmSkyFootprint *footprint, const gdouble ra, const gdouble dec);
static gdouble _ncm_sky_footprint_rectangular_ln_density (NcmSkyFootprint *footprint, const gdouble ra, const gdouble dec);

static void
ncm_sky_footprint_rectangular_class_init (NcmSkyFootprintRectangularClass *klass)
{
  GObjectClass *object_class           = G_OBJECT_CLASS (klass);
  NcmSkyFootprintClass *footprint_class = NCM_SKY_FOOTPRINT_CLASS (klass);

  object_class->set_property = &_ncm_sky_footprint_rectangular_set_property;
  object_class->get_property = &_ncm_sky_footprint_rectangular_get_property;

  /**
   * NcmSkyFootprintRectangular:ra-lim:
   *
   * The right ascension limits (min, max) in degrees.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_RA_LIM,
                                   g_param_spec_boxed ("ra-lim",
                                                       "RA limits",
                                                       "Right ascension limits (min, max) in degrees",
                                                       NCM_TYPE_DTUPLE2,
                                                       G_PARAM_READWRITE | G_PARAM_STATIC_STRINGS));

  /**
   * NcmSkyFootprintRectangular:dec-lim:
   *
   * The declination limits (min, max) in degrees.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_DEC_LIM,
                                   g_param_spec_boxed ("dec-lim",
                                                       "DEC limits",
                                                       "Declination limits (min, max) in degrees",
                                                       NCM_TYPE_DTUPLE2,
                                                       G_PARAM_READWRITE | G_PARAM_STATIC_STRINGS));

  footprint_class->gen_ra_dec = &_ncm_sky_footprint_rectangular_gen_ra_dec;
  footprint_class->contains   = &_ncm_sky_footprint_rectangular_contains;
  footprint_class->get_area   = &_ncm_sky_footprint_rectangular_get_area;
  footprint_class->density    = &_ncm_sky_footprint_rectangular_density;
  footprint_class->ln_density = &_ncm_sky_footprint_rectangular_ln_density;
}

static void
_ncm_sky_footprint_rectangular_gen_ra_dec (NcmSkyFootprint *footprint, NcmRNG *rng, gdouble *ra, gdouble *dec)
{
  NcmSkyFootprintRectangular *rect               = NCM_SKY_FOOTPRINT_RECTANGULAR (footprint);
  NcmSkyFootprintRectangularPrivate * const self = ncm_sky_footprint_rectangular_get_instance_private (rect);
  const gdouble sin_dec                          = ncm_rng_uniform_gen (rng, self->sin_dec_min, self->sin_dec_max);

  *ra  = ncm_rng_uniform_gen (rng, self->ra_min, self->ra_max);
  *dec = ncm_c_radian_to_degree (asin (sin_dec));
}

static gboolean
_ncm_sky_footprint_rectangular_contains (NcmSkyFootprint *footprint, const gdouble ra, const gdouble dec)
{
  NcmSkyFootprintRectangular *rect               = NCM_SKY_FOOTPRINT_RECTANGULAR (footprint);
  NcmSkyFootprintRectangularPrivate * const self = ncm_sky_footprint_rectangular_get_instance_private (rect);

  return (ra >= self->ra_min) && (ra <= self->ra_max) && (dec >= self->dec_min) && (dec <= self->dec_max);
}

static gdouble
_ncm_sky_footprint_rectangular_get_area (NcmSkyFootprint *footprint)
{
  NcmSkyFootprintRectangular *rect               = NCM_SKY_FOOTPRINT_RECTANGULAR (footprint);
  NcmSkyFootprintRectangularPrivate * const self = ncm_sky_footprint_rectangular_get_instance_private (rect);

  return ncm_c_degree_to_radian (self->ra_max - self->ra_min) * (self->sin_dec_max - self->sin_dec_min);
}

static gdouble
_ncm_sky_footprint_rectangular_density (NcmSkyFootprint *footprint, const gdouble ra, const gdouble dec)
{
  NcmSkyFootprintRectangular *rect               = NCM_SKY_FOOTPRINT_RECTANGULAR (footprint);
  NcmSkyFootprintRectangularPrivate * const self = ncm_sky_footprint_rectangular_get_instance_private (rect);

  if ((ra >= self->ra_min) && (ra <= self->ra_max) && (dec >= self->dec_min) && (dec <= self->dec_max))
    return self->ra_norm * self->dec_norm * cos (ncm_c_degree_to_radian (dec));

  return 0.0;
}

static gdouble
_ncm_sky_footprint_rectangular_ln_density (NcmSkyFootprint *footprint, const gdouble ra, const gdouble dec)
{
  NcmSkyFootprintRectangular *rect               = NCM_SKY_FOOTPRINT_RECTANGULAR (footprint);
  NcmSkyFootprintRectangularPrivate * const self = ncm_sky_footprint_rectangular_get_instance_private (rect);

  if ((ra >= self->ra_min) && (ra <= self->ra_max) && (dec >= self->dec_min) && (dec <= self->dec_max))
    return log (self->ra_norm * self->dec_norm * cos (ncm_c_degree_to_radian (dec)));

  return GSL_NEGINF;
}

/**
 * ncm_sky_footprint_rectangular_new:
 * @ra_min: minimum right ascension (degrees)
 * @ra_max: maximum right ascension (degrees)
 * @dec_min: minimum declination (degrees)
 * @dec_max: maximum declination (degrees)
 *
 * Creates a new #NcmSkyFootprintRectangular covering the given rectangle.
 *
 * Returns: (transfer full): a new #NcmSkyFootprintRectangular.
 */
NcmSkyFootprintRectangular *
ncm_sky_footprint_rectangular_new (const gdouble ra_min, const gdouble ra_max, const gdouble dec_min, const gdouble dec_max)
{
  NcmDTuple2 ra_lim  = NCM_DTUPLE2_STATIC_INIT (ra_min, ra_max);
  NcmDTuple2 dec_lim = NCM_DTUPLE2_STATIC_INIT (dec_min, dec_max);

  return g_object_new (NCM_TYPE_SKY_FOOTPRINT_RECTANGULAR,
                       "ra-lim", &ra_lim,
                       "dec-lim", &dec_lim,
                       NULL);
}

/**
 * ncm_sky_footprint_rectangular_ref:
 * @rect: a #NcmSkyFootprintRectangular
 *
 * Increases the reference count of @rect by one.
 *
 * Returns: (transfer full): @rect.
 */
NcmSkyFootprintRectangular *
ncm_sky_footprint_rectangular_ref (NcmSkyFootprintRectangular *rect)
{
  return g_object_ref (rect);
}

/**
 * ncm_sky_footprint_rectangular_free:
 * @rect: a #NcmSkyFootprintRectangular
 *
 * Decreases the reference count of @rect by one.
 *
 */
void
ncm_sky_footprint_rectangular_free (NcmSkyFootprintRectangular *rect)
{
  g_object_unref (rect);
}

/**
 * ncm_sky_footprint_rectangular_clear:
 * @rect: a #NcmSkyFootprintRectangular
 *
 * If *@rect is different from %NULL, decreases its reference count and sets
 * *@rect to %NULL.
 *
 */
void
ncm_sky_footprint_rectangular_clear (NcmSkyFootprintRectangular **rect)
{
  g_clear_object (rect);
}

/**
 * ncm_sky_footprint_rectangular_set_ra_lim:
 * @rect: a #NcmSkyFootprintRectangular
 * @ra_min: minimum right ascension (degrees)
 * @ra_max: maximum right ascension (degrees)
 *
 * Sets the right ascension limits of @rect.
 *
 */
void
ncm_sky_footprint_rectangular_set_ra_lim (NcmSkyFootprintRectangular *rect, const gdouble ra_min, const gdouble ra_max)
{
  NcmSkyFootprintRectangularPrivate * const self = ncm_sky_footprint_rectangular_get_instance_private (rect);

  g_assert_cmpfloat (ra_min, <, ra_max);

  self->ra_min  = ra_min;
  self->ra_max  = ra_max;
  self->ra_norm = 1.0 / (self->ra_max - self->ra_min);
}

/**
 * ncm_sky_footprint_rectangular_get_ra_lim:
 * @rect: a #NcmSkyFootprintRectangular
 * @ra_min: (out): the minimum right ascension (degrees)
 * @ra_max: (out): the maximum right ascension (degrees)
 *
 * Gets the right ascension limits of @rect.
 *
 */
void
ncm_sky_footprint_rectangular_get_ra_lim (NcmSkyFootprintRectangular *rect, gdouble *ra_min, gdouble *ra_max)
{
  NcmSkyFootprintRectangularPrivate * const self = ncm_sky_footprint_rectangular_get_instance_private (rect);

  g_assert_nonnull (ra_min);
  g_assert_nonnull (ra_max);

  *ra_min = self->ra_min;
  *ra_max = self->ra_max;
}

/**
 * ncm_sky_footprint_rectangular_set_dec_lim:
 * @rect: a #NcmSkyFootprintRectangular
 * @dec_min: minimum declination (degrees)
 * @dec_max: maximum declination (degrees)
 *
 * Sets the declination limits of @rect.
 *
 */
void
ncm_sky_footprint_rectangular_set_dec_lim (NcmSkyFootprintRectangular *rect, const gdouble dec_min, const gdouble dec_max)
{
  NcmSkyFootprintRectangularPrivate * const self = ncm_sky_footprint_rectangular_get_instance_private (rect);

  g_assert_cmpfloat (dec_min, <, dec_max);

  self->dec_min     = dec_min;
  self->dec_max     = dec_max;
  self->sin_dec_min = sin (ncm_c_degree_to_radian (dec_min));
  self->sin_dec_max = sin (ncm_c_degree_to_radian (dec_max));
  self->dec_norm    = ncm_c_pi () / (180.0 * (self->sin_dec_max - self->sin_dec_min));
}

/**
 * ncm_sky_footprint_rectangular_get_dec_lim:
 * @rect: a #NcmSkyFootprintRectangular
 * @dec_min: (out): the minimum declination (degrees)
 * @dec_max: (out): the maximum declination (degrees)
 *
 * Gets the declination limits of @rect.
 *
 */
void
ncm_sky_footprint_rectangular_get_dec_lim (NcmSkyFootprintRectangular *rect, gdouble *dec_min, gdouble *dec_max)
{
  NcmSkyFootprintRectangularPrivate * const self = ncm_sky_footprint_rectangular_get_instance_private (rect);

  g_assert_nonnull (dec_min);
  g_assert_nonnull (dec_max);

  *dec_min = self->dec_min;
  *dec_max = self->dec_max;
}

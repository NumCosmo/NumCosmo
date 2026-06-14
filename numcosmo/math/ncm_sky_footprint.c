/***************************************************************************
 *            ncm_sky_footprint.c
 *
 *  Sat Jun 13 19:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_sky_footprint.c
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
 * NcmSkyFootprint:
 *
 * Abstract sky region with sampling and density.
 *
 * A generic, cosmology-agnostic description of a region of the celestial sphere
 * in right ascension and declination. It provides uniform sampling of positions
 * within the region (ncm_sky_footprint_gen_ra_dec()), a membership test
 * (ncm_sky_footprint_contains()), the solid angle (ncm_sky_footprint_get_area())
 * and the normalized position density (ncm_sky_footprint_density() and its log).
 *
 * It is meant to be shared by composition: both the galaxy sample position
 * distribution used in weak-lensing likelihoods and the halo/cluster mock
 * generators hold a footprint and delegate their sky sampling and geometry to
 * it. Concrete shapes (for example #NcmSkyFootprintRectangular) implement the
 * virtual methods.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"
#include "math/ncm_sky_footprint.h"

G_DEFINE_ABSTRACT_TYPE (NcmSkyFootprint, ncm_sky_footprint, G_TYPE_OBJECT)

static void
ncm_sky_footprint_init (NcmSkyFootprint *footprint)
{
}

/* LCOV_EXCL_START */
static void
_ncm_sky_footprint_gen_ra_dec (NcmSkyFootprint *footprint, NcmRNG *rng, gdouble *ra, gdouble *dec)
{
  g_error ("_ncm_sky_footprint_gen_ra_dec: method not implemented by `%s'.", G_OBJECT_TYPE_NAME (footprint));
}

static gboolean
_ncm_sky_footprint_contains (NcmSkyFootprint *footprint, const gdouble ra, const gdouble dec)
{
  g_error ("_ncm_sky_footprint_contains: method not implemented by `%s'.", G_OBJECT_TYPE_NAME (footprint));

  return FALSE;
}

static gdouble
_ncm_sky_footprint_get_area (NcmSkyFootprint *footprint)
{
  g_error ("_ncm_sky_footprint_get_area: method not implemented by `%s'.", G_OBJECT_TYPE_NAME (footprint));

  return 0.0;
}

static gdouble
_ncm_sky_footprint_density (NcmSkyFootprint *footprint, const gdouble ra, const gdouble dec)
{
  g_error ("_ncm_sky_footprint_density: method not implemented by `%s'.", G_OBJECT_TYPE_NAME (footprint));

  return 0.0;
}

static gdouble
_ncm_sky_footprint_ln_density (NcmSkyFootprint *footprint, const gdouble ra, const gdouble dec)
{
  g_error ("_ncm_sky_footprint_ln_density: method not implemented by `%s'.", G_OBJECT_TYPE_NAME (footprint));

  return 0.0;
}

/* LCOV_EXCL_STOP */

static void
ncm_sky_footprint_class_init (NcmSkyFootprintClass *klass)
{
  klass->gen_ra_dec = &_ncm_sky_footprint_gen_ra_dec;
  klass->contains   = &_ncm_sky_footprint_contains;
  klass->get_area   = &_ncm_sky_footprint_get_area;
  klass->density    = &_ncm_sky_footprint_density;
  klass->ln_density = &_ncm_sky_footprint_ln_density;
}

/**
 * ncm_sky_footprint_ref:
 * @footprint: a #NcmSkyFootprint
 *
 * Increases the reference count of @footprint by one.
 *
 * Returns: (transfer full): @footprint.
 */
NcmSkyFootprint *
ncm_sky_footprint_ref (NcmSkyFootprint *footprint)
{
  return g_object_ref (footprint);
}

/**
 * ncm_sky_footprint_free:
 * @footprint: a #NcmSkyFootprint
 *
 * Decreases the reference count of @footprint by one.
 *
 */
void
ncm_sky_footprint_free (NcmSkyFootprint *footprint)
{
  g_object_unref (footprint);
}

/**
 * ncm_sky_footprint_clear:
 * @footprint: a #NcmSkyFootprint
 *
 * If *@footprint is different from %NULL, decreases its reference count and sets
 * *@footprint to %NULL.
 *
 */
void
ncm_sky_footprint_clear (NcmSkyFootprint **footprint)
{
  g_clear_object (footprint);
}

/**
 * ncm_sky_footprint_gen_ra_dec:
 * @footprint: a #NcmSkyFootprint
 * @rng: a #NcmRNG
 * @ra: (out): the sampled right ascension (degrees)
 * @dec: (out): the sampled declination (degrees)
 *
 * Samples a position uniformly distributed within @footprint.
 *
 */
void
ncm_sky_footprint_gen_ra_dec (NcmSkyFootprint *footprint, NcmRNG *rng, gdouble *ra, gdouble *dec)
{
  NCM_SKY_FOOTPRINT_GET_CLASS (footprint)->gen_ra_dec (footprint, rng, ra, dec);
}

/**
 * ncm_sky_footprint_contains:
 * @footprint: a #NcmSkyFootprint
 * @ra: a right ascension (degrees)
 * @dec: a declination (degrees)
 *
 * Checks whether the position (@ra, @dec) lies within @footprint.
 *
 * Returns: %TRUE if the position is inside @footprint, %FALSE otherwise.
 */
gboolean
ncm_sky_footprint_contains (NcmSkyFootprint *footprint, const gdouble ra, const gdouble dec)
{
  return NCM_SKY_FOOTPRINT_GET_CLASS (footprint)->contains (footprint, ra, dec);
}

/**
 * ncm_sky_footprint_get_area:
 * @footprint: a #NcmSkyFootprint
 *
 * Gets the solid angle covered by @footprint.
 *
 * Returns: the area of @footprint in steradians.
 */
gdouble
ncm_sky_footprint_get_area (NcmSkyFootprint *footprint)
{
  return NCM_SKY_FOOTPRINT_GET_CLASS (footprint)->get_area (footprint);
}

/**
 * ncm_sky_footprint_density:
 * @footprint: a #NcmSkyFootprint
 * @ra: a right ascension (degrees)
 * @dec: a declination (degrees)
 *
 * Computes the normalized position density at (@ra, @dec), in units of
 * inverse square degree, returning zero outside the footprint.
 *
 * Returns: the position density at (@ra, @dec).
 */
gdouble
ncm_sky_footprint_density (NcmSkyFootprint *footprint, const gdouble ra, const gdouble dec)
{
  return NCM_SKY_FOOTPRINT_GET_CLASS (footprint)->density (footprint, ra, dec);
}

/**
 * ncm_sky_footprint_ln_density:
 * @footprint: a #NcmSkyFootprint
 * @ra: a right ascension (degrees)
 * @dec: a declination (degrees)
 *
 * Computes the natural logarithm of the normalized position density at
 * (@ra, @dec), returning negative infinity outside the footprint.
 *
 * Returns: the log position density at (@ra, @dec).
 */
gdouble
ncm_sky_footprint_ln_density (NcmSkyFootprint *footprint, const gdouble ra, const gdouble dec)
{
  return NCM_SKY_FOOTPRINT_GET_CLASS (footprint)->ln_density (footprint, ra, dec);
}


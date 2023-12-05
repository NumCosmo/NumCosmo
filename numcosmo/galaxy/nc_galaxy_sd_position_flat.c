/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 2; tab-width: 2 -*-  */

/***************************************************************************
 *            nc_galaxy_sd_position_flat.c
 *
 *  Wed March 1 12:51:22 2023
 *  Copyright  2023  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * nc_galaxy_sd_position_flat.c
 * Copyright (C) 2023 Caio Lima de Oliveira <caiolimadeoliveira@pm.me>
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
 * SECTION:nc_galaxy_sd_position_flat
 * @title: NcGalaxySDPositionFlat
 * @short_description: Class describing galaxy sample position distributions with flat distribution
 * @stability: Unstable
 *
 *
 * Class defining a galaxy sample position distribution with flat
 * probability distribution $P(z)$.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

/* #include "nc_enum_types.h" */
#include "galaxy/nc_galaxy_sd_position_flat.h"
#include "galaxy/nc_galaxy_sd_position.h"
#include <math.h>
#include <gsl/gsl_math.h>

struct _NcGalaxySDPositionFlatPrivate
{
  NcmVector *z_lim;
  NcmVector *r_lim;
};

enum
{
  PROP_0,
  PROP_Z_LIM,
  PROP_R_LIM,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcGalaxySDPositionFlat, nc_galaxy_sd_position_flat, NC_TYPE_GALAXY_SD_POSITION);

static void
nc_galaxy_sd_position_flat_init (NcGalaxySDPositionFlat *gsdpflat)
{
  NcGalaxySDPositionFlatPrivate * const self = gsdpflat->priv = nc_galaxy_sd_position_flat_get_instance_private (gsdpflat);

  self->z_lim = NULL;
  self->r_lim = NULL;
}

static void
_nc_galaxy_sd_position_flat_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcGalaxySDPositionFlat *gsdpflat = NC_GALAXY_SD_POSITION_FLAT (object);

  g_return_if_fail (NC_IS_GALAXY_SD_POSITION_FLAT (object));

  switch (prop_id)
  {
    case PROP_Z_LIM:
      nc_galaxy_sd_position_flat_set_z_lim (gsdpflat, g_value_get_object (value));
      break;
    case PROP_R_LIM:
      nc_galaxy_sd_position_flat_set_r_lim (gsdpflat, g_value_get_object (value));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_galaxy_sd_position_flat_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcGalaxySDPositionFlat *gsdpflat = NC_GALAXY_SD_POSITION_FLAT (object);

  g_return_if_fail (NC_IS_GALAXY_SD_POSITION_FLAT (object));

  switch (prop_id)
  {
    case PROP_Z_LIM:
      g_value_set_object (value, nc_galaxy_sd_position_flat_peek_z_lim (gsdpflat));
      break;
    case PROP_R_LIM:
      g_value_set_object (value, nc_galaxy_sd_position_flat_peek_r_lim (gsdpflat));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_galaxy_sd_position_flat_dispose (GObject *object)
{
  NcGalaxySDPositionFlat *gsdpflat           = NC_GALAXY_SD_POSITION_FLAT (object);
  NcGalaxySDPositionFlatPrivate * const self = gsdpflat->priv;

  ncm_vector_clear (&self->z_lim);
  ncm_vector_clear (&self->r_lim);

  G_OBJECT_CLASS (nc_galaxy_sd_position_flat_parent_class)->dispose (object);
}

static void
_nc_galaxy_sd_position_flat_finalize (GObject *object)
{
  G_OBJECT_CLASS (nc_galaxy_sd_position_flat_parent_class)->finalize (object);
}

static gdouble _nc_galaxy_sd_position_flat_gen_r (NcGalaxySDPosition *gsdp, NcmRNG *rng);
static gdouble _nc_galaxy_sd_position_flat_gen_z (NcGalaxySDPosition *gsdp, NcmRNG *rng);
static gdouble _nc_galaxy_sd_position_flat_integ (NcGalaxySDPosition *gsdp, const gdouble r, const gdouble z);

static void
nc_galaxy_sd_position_flat_class_init (NcGalaxySDPositionFlatClass *klass)
{
  NcGalaxySDPositionClass *sd_position_class = NC_GALAXY_SD_POSITION_CLASS (klass);
  GObjectClass *object_class                 = G_OBJECT_CLASS (klass);

  object_class->set_property = &_nc_galaxy_sd_position_flat_set_property;
  object_class->get_property = &_nc_galaxy_sd_position_flat_get_property;
  object_class->dispose      = &_nc_galaxy_sd_position_flat_dispose;
  object_class->finalize     = &_nc_galaxy_sd_position_flat_finalize;

  /**
   * NcGalaxySDPositionFlat:z-lim:
   *
   * Galaxy sample redshift distribution limits.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_Z_LIM,
                                   g_param_spec_object ("z-lim",
                                                        NULL,
                                                        "Galaxy sample redshift distribution limits",
                                                        NCM_TYPE_VECTOR,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcGalaxySDPositionFlat:R-lim:
   *
   * Galaxy sample radius distribution limits.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_R_LIM,
                                   g_param_spec_object ("r-lim",
                                                        NULL,
                                                        "Galaxy sample radius distribution limits",
                                                        NCM_TYPE_VECTOR,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  sd_position_class->gen_r = &_nc_galaxy_sd_position_flat_gen_r;
  sd_position_class->gen_z = &_nc_galaxy_sd_position_flat_gen_z;
  sd_position_class->integ = &_nc_galaxy_sd_position_flat_integ;
}

static gdouble
_nc_galaxy_sd_position_flat_gen_r (NcGalaxySDPosition *gsdp, NcmRNG *rng)
{
  NcGalaxySDPositionFlat *gsdpflat           = NC_GALAXY_SD_POSITION_FLAT (gsdp);
  NcGalaxySDPositionFlatPrivate * const self = gsdpflat->priv;
  const gdouble r_lb                         = ncm_vector_get (self->r_lim, 0);
  const gdouble r_ub                         = ncm_vector_get (self->r_lim, 1);
  const gdouble r_lb2                        = r_lb * r_lb;
  const gdouble r_ub2                        = r_ub * r_ub;
  gdouble cumul_gen                          = ncm_rng_uniform_gen (rng, 0.0, 1.0);

  return sqrt (cumul_gen * (r_ub2 - r_lb2) + r_lb2);
}

static gdouble
_nc_galaxy_sd_position_flat_gen_z (NcGalaxySDPosition *gsdp, NcmRNG *rng)
{
  NcGalaxySDPositionFlat *gsdpflat           = NC_GALAXY_SD_POSITION_FLAT (gsdp);
  NcGalaxySDPositionFlatPrivate * const self = gsdpflat->priv;

  return ncm_rng_uniform_gen (rng, ncm_vector_get (self->z_lim, 0), ncm_vector_get (self->z_lim, 1));
}

static gdouble
_nc_galaxy_sd_position_flat_integ (NcGalaxySDPosition *gsdp, const gdouble r, const gdouble z)
{
  NcGalaxySDPositionFlat *gsdpflat           = NC_GALAXY_SD_POSITION_FLAT (gsdp);
  NcGalaxySDPositionFlatPrivate * const self = gsdpflat->priv;

  if ((z >= ncm_vector_get (self->z_lim, 0)) && (z <= ncm_vector_get (self->z_lim, 1)) &&
      (r >= ncm_vector_get (self->r_lim, 0)) && (r <= ncm_vector_get (self->r_lim, 1)))
    return r;

  g_assert_not_reached ();

  return 0.0;
}

/**
 * nc_galaxy_sd_position_flat_new:
 *
 * Creates a new #NcGalaxySDPositionFlat
 *
 * Returns: (transfer full): a new NcGalaxySDPositionFlat.
 */
NcGalaxySDPositionFlat *
nc_galaxy_sd_position_flat_new ()
{
  NcGalaxySDPositionFlat *gsdpflat = g_object_new (NC_TYPE_GALAXY_SD_POSITION_FLAT,
                                                   NULL);

  return gsdpflat;
}

/**
 * nc_galaxy_sd_position_flat_ref:
 * @gsdpflat: a #NcGalaxySDPositionFlat
 *
 * Increase the reference of @gsdpflat by one.
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
 * Decrease the reference count of @gsdpflat by one.
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
 * nc_galaxy_sd_position_flat_set_z_lim:
 * @gsdpflat: a #NcGalaxySDPositionFlat
 * @lim: a #NcmVector
 *
 * Sets the redshift limits @lim.
 */
void
nc_galaxy_sd_position_flat_set_z_lim (NcGalaxySDPositionFlat *gsdpflat, NcmVector *lim)
{
  NcGalaxySDPositionFlatPrivate * const self = gsdpflat->priv;

  g_assert_cmpuint (ncm_vector_len (lim), ==, 2);

  ncm_vector_clear (&self->z_lim);

  self->z_lim = ncm_vector_ref (lim);
}

/**
 * nc_galaxy_sd_position_flat_peek_z_lim:
 * @gsdpflat: a #NcGalaxySDPositionFlat
 *
 * Gets the redshift limits.
 *
 * Returns: (transfer none): the redshift limits.
 */
NcmVector *
nc_galaxy_sd_position_flat_peek_z_lim (NcGalaxySDPositionFlat *gsdpflat)
{
  NcGalaxySDPositionFlatPrivate * const self = gsdpflat->priv;

  return self->z_lim;
}

/**
 * nc_galaxy_sd_position_flat_set_r_lim:
 * @gsdpflat: a #NcGalaxySDPositionFlat
 * @lim: a #NcmVector
 *
 * Sets the radius limits @lim.
 */
void
nc_galaxy_sd_position_flat_set_r_lim (NcGalaxySDPositionFlat *gsdpflat, NcmVector *lim)
{
  NcGalaxySDPositionFlatPrivate * const self = gsdpflat->priv;

  g_assert_cmpuint (ncm_vector_len (lim), ==, 2);

  ncm_vector_clear (&self->r_lim);

  self->r_lim = ncm_vector_ref (lim);
}

/**
 * nc_galaxy_sd_position_flat_peek_r_lim:
 * @gsdpflat: a #NcGalaxySDPositionFlat
 *
 * Gets the radius limits.
 *
 * Returns: (transfer none): the radius limits.
 */
NcmVector *
nc_galaxy_sd_position_flat_peek_r_lim (NcGalaxySDPositionFlat *gsdpflat)
{
  NcGalaxySDPositionFlatPrivate * const self = gsdpflat->priv;

  return self->r_lim;
}


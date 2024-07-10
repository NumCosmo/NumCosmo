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
#include "math/ncm_dtuple.h"
#include <math.h>
#include <gsl/gsl_math.h>

typedef struct _NcGalaxySDPositionFlatPrivate
{
  gdouble z_lb;
  gdouble z_ub;
  gdouble theta_norm;
  gdouble theta_lb;
  gdouble theta_ub;
  gdouble theta_lb2;
  gdouble theta_ub2;
} NcGalaxySDPositionFlatPrivate;

struct _NcGalaxySDPositionFlat
{
  NcGalaxySDPosition parent_instance;
};

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

  self->z_lb       = 0.0;
  self->z_ub       = 0.0;
  self->theta_norm = 0.0;
  self->theta_lb   = 0.0;
  self->theta_ub   = 0.0;
  self->theta_lb2  = 0.0;
  self->theta_ub2  = 0.0;
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

static gdouble _nc_galaxy_sd_position_flat_gen_theta (NcGalaxySDPosition *gsdp, NcmRNG *rng);
static gdouble _nc_galaxy_sd_position_flat_gen_z (NcGalaxySDPosition *gsdp, NcmRNG *rng);
static gdouble _nc_galaxy_sd_position_flat_integ (NcGalaxySDPosition *gsdp, const gdouble r, const gdouble z);
static void _nc_galaxy_sd_position_flat_set_z_lim (NcGalaxySDPosition *gsdp, gdouble z_min, gdouble z_max);
static void _nc_galaxy_sd_position_flat_get_z_lim (NcGalaxySDPosition *gsdp, gdouble *z_min, gdouble *z_max);
static void _nc_galaxy_sd_position_flat_set_theta_lim (NcGalaxySDPosition *gsdp, gdouble theta_min, gdouble theta_max);
static void _nc_galaxy_sd_position_flat_get_theta_lim (NcGalaxySDPosition *gsdp, gdouble *theta_min, gdouble *theta_max);

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

  sd_position_class->gen_theta     = &_nc_galaxy_sd_position_flat_gen_theta;
  sd_position_class->gen_z         = &_nc_galaxy_sd_position_flat_gen_z;
  sd_position_class->integ         = &_nc_galaxy_sd_position_flat_integ;
  sd_position_class->set_z_lim     = &_nc_galaxy_sd_position_flat_set_z_lim;
  sd_position_class->get_z_lim     = &_nc_galaxy_sd_position_flat_get_z_lim;
  sd_position_class->set_theta_lim = &_nc_galaxy_sd_position_flat_set_theta_lim;
  sd_position_class->get_theta_lim = &_nc_galaxy_sd_position_flat_get_theta_lim;
}

static gdouble
_nc_galaxy_sd_position_flat_gen_theta (NcGalaxySDPosition *gsdp, NcmRNG *rng)
{
  NcGalaxySDPositionFlat *gsdpflat           = NC_GALAXY_SD_POSITION_FLAT (gsdp);
  NcGalaxySDPositionFlatPrivate * const self = nc_galaxy_sd_position_flat_get_instance_private (gsdpflat);
  const gdouble cumul_gen                    = ncm_rng_uniform_gen (rng, 0.0, 1.0);

  return sqrt (cumul_gen * 2.0 / self->theta_norm + self->theta_lb2);
}

static gdouble
_nc_galaxy_sd_position_flat_gen_z (NcGalaxySDPosition *gsdp, NcmRNG *rng)
{
  NcGalaxySDPositionFlat *gsdpflat           = NC_GALAXY_SD_POSITION_FLAT (gsdp);
  NcGalaxySDPositionFlatPrivate * const self = nc_galaxy_sd_position_flat_get_instance_private (gsdpflat);

  return ncm_rng_uniform_gen (rng, self->z_lb, self->z_ub);
}

static gdouble
_nc_galaxy_sd_position_flat_integ (NcGalaxySDPosition *gsdp, const gdouble r, const gdouble z)
{
  NcGalaxySDPositionFlat *gsdpflat           = NC_GALAXY_SD_POSITION_FLAT (gsdp);
  NcGalaxySDPositionFlatPrivate * const self = nc_galaxy_sd_position_flat_get_instance_private (gsdpflat);

  if ((z >= self->z_lb) && (z <= self->z_ub) && (r >= self->theta_lb) && (r <= self->theta_ub))
    return r;

  return 0.0;
}

static void
_nc_galaxy_sd_position_flat_set_z_lim (NcGalaxySDPosition *gsdp, gdouble z_min, gdouble z_max)
{
  NcGalaxySDPositionFlat *gsdpflat           = NC_GALAXY_SD_POSITION_FLAT (gsdp);
  NcGalaxySDPositionFlatPrivate * const self = nc_galaxy_sd_position_flat_get_instance_private (gsdpflat);

  g_assert_cmpfloat (z_min, <, z_max);

  self->z_lb = z_min;
  self->z_ub = z_max;
}

static void
_nc_galaxy_sd_position_flat_get_z_lim (NcGalaxySDPosition *gsdp, gdouble *z_min, gdouble *z_max)
{
  NcGalaxySDPositionFlat *gsdpflat           = NC_GALAXY_SD_POSITION_FLAT (gsdp);
  NcGalaxySDPositionFlatPrivate * const self = nc_galaxy_sd_position_flat_get_instance_private (gsdpflat);

  g_assert_nonnull (z_min);
  g_assert_nonnull (z_max);

  *z_min = self->z_lb;
  *z_max = self->z_ub;
}

static void
_nc_galaxy_sd_position_flat_set_theta_lim (NcGalaxySDPosition *gsdp, gdouble theta_min, gdouble theta_max)
{
  NcGalaxySDPositionFlat *gsdpflat           = NC_GALAXY_SD_POSITION_FLAT (gsdp);
  NcGalaxySDPositionFlatPrivate * const self = nc_galaxy_sd_position_flat_get_instance_private (gsdpflat);

  g_assert_cmpfloat (theta_min, <, theta_max);

  self->theta_lb = theta_min;
  self->theta_ub = theta_max;

  self->theta_lb2  = self->theta_lb * self->theta_lb;
  self->theta_ub2  = self->theta_ub * self->theta_ub;
  self->theta_norm = 2.0 / (self->theta_ub2 - self->theta_lb2);
}

static void
_nc_galaxy_sd_position_flat_get_theta_lim (NcGalaxySDPosition *gsdp, gdouble *theta_min, gdouble *theta_max)
{
  NcGalaxySDPositionFlat *gsdpflat           = NC_GALAXY_SD_POSITION_FLAT (gsdp);
  NcGalaxySDPositionFlatPrivate * const self = nc_galaxy_sd_position_flat_get_instance_private (gsdpflat);

  g_assert_nonnull (theta_min);
  g_assert_nonnull (theta_max);

  *theta_min = self->theta_lb;
  *theta_max = self->theta_ub;
}

/**
 * nc_galaxy_sd_position_flat_new:
 * @z_min: minimum redshift
 * @z_max: maximum redshift
 * @theta_min: minimum radius
 * @theta_max: maximum radius
 *
 * Creates a new #NcGalaxySDPositionFlat
 *
 * Returns: (transfer full): a new NcGalaxySDPositionFlat.
 */
NcGalaxySDPositionFlat *
nc_galaxy_sd_position_flat_new (const gdouble z_min, const gdouble z_max, const gdouble theta_min, const gdouble theta_max)
{
  NcmDTuple2 z_lim                 = NCM_DTUPLE2_STATIC_INIT (z_min, z_max);
  NcmDTuple2 theta_lim             = NCM_DTUPLE2_STATIC_INIT (theta_min, theta_max);
  NcGalaxySDPositionFlat *gsdpflat = g_object_new (NC_TYPE_GALAXY_SD_POSITION_FLAT,
                                                   "z-lim", &z_lim,
                                                   "r-lim", &theta_lim,
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


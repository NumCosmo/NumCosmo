/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 2; tab-width: 2 -*-  */

/***************************************************************************
 *            nc_galaxy_sd_position.c
 *
 *  Sat May 20 17:52:48 2023
 *  Copyright  2023  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * nc_galaxy_sd_position.c
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
 * SECTION: nc_galaxy_sd_position
 * @title: NcGalaxySDPosition
 * @short_description: Class describing galaxy sample position distributions.
 * @stability: Unstable
 *
 *
 * This class describes a galaxy sample position distributions.
 * It is composed by two distributions: an image position distribution $P(r)$ and a redshift distribution $P(z)$.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "galaxy/nc_galaxy_sd_position.h"
#include <math.h>
#include <gsl/gsl_math.h>


struct _NcGalaxySDPositionPrivate
{
  gint placeholder;
};

G_DEFINE_TYPE_WITH_PRIVATE (NcGalaxySDPosition, nc_galaxy_sd_position, G_TYPE_OBJECT);

static void
nc_galaxy_sd_position_init (NcGalaxySDPosition *gsdp)
{
  gsdp->priv = nc_galaxy_sd_position_get_instance_private (gsdp);
}

static void
_nc_galaxy_sd_position_finalize (GObject *object)
{
  G_OBJECT_CLASS (nc_galaxy_sd_position_parent_class)->finalize (object);
}

static void
_nc_galaxy_sd_position_gen_r (NcGalaxySDPosition *gsdp, NcmRNG *rng, gdouble *gen_r)
{
  g_error ("_nc_galaxy_sd_position_gen_r: method not implemented.");
}

static void
_nc_galaxy_sd_position_gen_z (NcGalaxySDPosition *gsdp, NcmRNG *rng, gdouble *gen_z)
{
  g_error ("_nc_galaxy_sd_position_gen_z: method not implemented.");
}

static gdouble
_nc_galaxy_sd_position_integ (NcGalaxySDPosition *gsdp, NcmVector *pos)
{
  g_error ("_nc_galaxy_sd_position_integ: method not implemented.");

  return 0.0;
}

static void
nc_galaxy_sd_position_class_init (NcGalaxySDPositionClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->finalize = &_nc_galaxy_sd_position_finalize;

  klass->gen_r = &_nc_galaxy_sd_position_gen_r;
  klass->gen_z = &_nc_galaxy_sd_position_gen_z;
  klass->integ = &_nc_galaxy_sd_position_integ;
}

/**
 * nc_galaxy_sd_position_ref:
 * @gsdp: a #NcGalaxySDPosition
 *
 * Increase the reference of @gsdp by one.
 *
 * Returns: (transfer full): @gsdp.
 */
NcGalaxySDPosition *
nc_galaxy_sd_position_ref (NcGalaxySDPosition *gsdp)
{
  return g_object_ref (gsdp);
}

/**
 * nc_galaxy_sd_position_free:
 * @gsdp: a #NcGalaxySDPosition
 *
 * Decrease the reference count of @gsdp by one.
 *
 */
void
nc_galaxy_sd_position_free (NcGalaxySDPosition *gsdp)
{
  g_object_unref (gsdp);
}

/**
 * nc_galaxy_sd_position_clear:
 * @gsdp: a #NcGalaxySDPosition
 *
 * Decrease the reference count of @gsdp by one, and sets the pointer *@gsdp to
 * NULL.
 *
 */
void
nc_galaxy_sd_position_clear (NcGalaxySDPosition **gsdp)
{
  g_clear_object (gsdp);
}

/**
 * nc_galaxy_sd_position_gen_r: (virtual gen_r)
 * @gsdp: a #NcGalaxySDPosition
 * @rng: a #NcmRNG
 * @gen_r: a #gdouble
 *
 * Generates a $r$ value from the distribution using @rng
 * and saves it in @r.
 *
 */
/**
 * nc_galaxy_sd_position_gen_z: (virtual gen_z)
 * @gsdp: a #NcGalaxySDPosition
 * @rng: a #NcmRNG
 * @gen_z: a #gdouble
 *
 * Generates a $z$ value from the distribution using @rng
 * and saves it in @z.
 *
 */
/**
 * nc_galaxy_sd_position_integ: (virtual integ)
 * @gsdp: a #NcGalaxySDPosition
 * @pos: a #NcmVector stating position and redshift
 *
 * Computes the probability density of the observables $r$ and $z$ given the redshift.
 * The probability density is given by $P(z)P(r)$.
 *
 * Returns: the probability density at $(r, z)$, $P(z)P(r)$.
 */


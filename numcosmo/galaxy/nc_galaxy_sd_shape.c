/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 2; tab-width: 2 -*-  */

/***************************************************************************
 *            nc_galaxy_sd_shape.c
 *
 *  Sat May 21 20:43:32 2023
 *  Copyright  2023  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * nc_galaxy_sd_shape.c
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
 * SECTION: nc_galaxy_sd_shape
 * @title: NcGalaxySDShape
 * @short_description: Class describing galaxy sample shape distribution.
 * @stability: Unstable
 *
 *
 * This class describes a galaxy sample shape distribution.
 * It is composed by a distribution $P(s)$.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "galaxy/nc_galaxy_sd_shape.h"
#include <math.h>
#include <gsl/gsl_math.h>


struct _NcGalaxySDShapePrivate
{
  gint placeholder;
};

G_DEFINE_TYPE_WITH_PRIVATE (NcGalaxySDShape, nc_galaxy_sd_shape, G_TYPE_OBJECT);

static void
nc_galaxy_sd_shape_init (NcGalaxySDShape *gsds)
{
  gsds->priv = nc_galaxy_sd_shape_get_instance_private (gsds);
}

static void
_nc_galaxy_sd_shape_finalize (GObject *object)
{
  G_OBJECT_CLASS (nc_galaxy_sd_shape_parent_class)->finalize (object);
}

static void
_nc_galaxy_sd_shape_gen (NcGalaxySDShape *gsds, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, const gdouble z_cluster, NcmRNG *rng, const gdouble r, const gdouble z, gdouble *et, gdouble *ex)
{
  g_error ("_nc_galaxy_sd_shape_gen: method not implemented.");
}

static gdouble
_nc_galaxy_sd_shape_integ (NcGalaxySDShape *gsds, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, const gdouble z_cluster, const gdouble r, const gdouble z, const gdouble et, const gdouble ex)
{
  g_error ("_nc_galaxy_sd_shape_integ: method not implemented.");

  return 0.0;
}

void
_nc_galaxy_sd_shape_integ_optzs_prep (NcGalaxySDShape *gsds, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, const gdouble z_cluster, const gdouble R)
{
  g_error ("_nc_galaxy_sd_shape_integ_optzs_prep: method not implemented.");
}

static gdouble
_nc_galaxy_sd_shape_integ_optzs (NcGalaxySDShape *gsds, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, const gdouble z_cluster, const gdouble z_source, const gdouble et, const gdouble ex)
{
  g_error ("_nc_galaxy_sd_shape_integ_optzs: method not implemented.");

  return 0.0;
}

static void
nc_galaxy_sd_shape_class_init (NcGalaxySDShapeClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->finalize = &_nc_galaxy_sd_shape_finalize;

  klass->gen   = &_nc_galaxy_sd_shape_gen;
  klass->integ = &_nc_galaxy_sd_shape_integ;
}

/**
 * nc_galaxy_sd_shape_ref:
 * @gsds: a #NcGalaxySDShape
 *
 * Increase the reference of @gsds by one.
 *
 * Returns: (transfer full): @gsds.
 */
NcGalaxySDShape *
nc_galaxy_sd_shape_ref (NcGalaxySDShape *gsds)
{
  return g_object_ref (gsds);
}

/**
 * nc_galaxy_sd_shape_free:
 * @gsds: a #NcGalaxySDShape
 *
 * Decrease the reference count of @gsds by one.
 *
 */
void
nc_galaxy_sd_shape_free (NcGalaxySDShape *gsds)
{
  g_object_unref (gsds);
}

/**
 * nc_galaxy_sd_shape_clear:
 * @gsds: a #NcGalaxySDShape
 *
 * Decrease the reference count of @gsds by one, and sets the pointer *@gsds to
 * NULL.
 *
 */
void
nc_galaxy_sd_shape_clear (NcGalaxySDShape **gsds)
{
  g_clear_object (gsds);
}

/**
 * nc_galaxy_sd_shape_gen: (virtual gen)
 * @gsds: a #NcGalaxySDShape
 * @cosmo: a #NcHICosmo
 * @dp: a #NcHaloDensityProfile
 * @smd: a #NcWLSurfaceMassDensity
 * @z_cluster: cluster redshift $z_\mathrm{cl}
 * @rng: a #NcmRNG
 * @r: a #gdouble
 * @z: a #gdouble
 * @et: (out): the generated tangential ellipticity component $e_\mathrm{t}$
 * @ex: (out): the generated cross ellipticity component $e_\mathrm{x}$
 *
 * Generates a shape value from the position using @rng.
 *
 */
void
nc_galaxy_sd_shape_gen (NcGalaxySDShape *gsds, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, const gdouble z_cluster, NcmRNG *rng, gdouble r, gdouble z, gdouble *et, gdouble *ex)
{
  NC_GALAXY_SD_SHAPE_GET_CLASS (gsds)->gen (gsds, cosmo, dp, smd, z_cluster, rng, r, z, et, ex);
}

/**
 * nc_galaxy_sd_shape_integ: (virtual integ)
 * @gsds: a #NcGalaxySDShape
 * @cosmo: a #NcHICosmo
 * @dp: a #NcHaloDensityProfile
 * @smd: a #NcWLSurfaceMassDensity
 * @z_cluster: cluster redshift $z_\mathrm{cl}$
 * @r: the projected radius $r$
 * @z: the redshift $z$
 * @et: the tangential ellipticity component $e_\mathrm{t}$
 * @ex: the cross ellipticity component $e_\mathrm{x}$
 *
 *
 * Computes the probability density of the observable shape given the position.
 * The probability density is given by $P(s)P$.
 *
 * Returns: the probability density of observable shape, $P(s)$.
 */
gdouble
nc_galaxy_sd_shape_integ (NcGalaxySDShape *gsds, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, const gdouble z_cluster, const gdouble r, const gdouble z, const gdouble et, const gdouble ex)
{
  return NC_GALAXY_SD_SHAPE_GET_CLASS (gsds)->integ (gsds, cosmo, dp, smd, z_cluster, r, z, et, ex);
}

/**
 * nc_galaxy_sd_shape_integ_optzs_prep: (virtual integ_optzs_prep)
 * @gsds: a #NcGalaxySDShape
 * @cosmo: a #NcHICosmo
 * @dp: a #NcHaloDensityProfile
 * @smd: a #NcWLSurfaceMassDensity
 * @z_cluster: cluster redshift $z_\mathrm{cl}$
 * @R: the projected radius $R$
 *
 *
 * Prepares $gsds$ to compute the probability density of the observable shape given the position.
 *
 */
void
nc_galaxy_sd_shape_integ_optzs_prep (NcGalaxySDShape *gsds, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, const gdouble z_cluster, const gdouble R)
{
  NC_GALAXY_SD_SHAPE_GET_CLASS (gsds)->integ_optzs_prep (gsds, cosmo, dp, smd, z_cluster, R);
}

/**
 * nc_galaxy_sd_shape_integ_optzs: (virtual integ_optzs)
 * @gsds: a #NcGalaxySDShape
 * @cosmo: a #NcHICosmo
 * @dp: a #NcHaloDensityProfile
 * @smd: a #NcWLSurfaceMassDensity
 * @z_cluster: cluster redshift $z_\mathrm{cl}$
 * @z_source: the redshift $z$
 * @et: the tangential ellipticity component $e_\mathrm{t}$
 * @ex: the cross ellipticity component $e_\mathrm{x}$
 *
 *
 * Computes the probability density of the observable shape given the position.
 * The probability density is given by $P(s)P$.
 *
 * Returns: the probability density of observable shape, $P(s)$.
 */
gdouble
nc_galaxy_sd_shape_integ_optzs (NcGalaxySDShape *gsds, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, const gdouble z_cluster, const gdouble z_source, const gdouble et, const gdouble ex)
{
  return NC_GALAXY_SD_SHAPE_GET_CLASS (gsds)->integ_optzs (gsds, cosmo, dp, smd, z_cluster, z_source, et, ex);
}


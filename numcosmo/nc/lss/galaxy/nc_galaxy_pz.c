/***************************************************************************
 *            nc_galaxy_pz.c
 *
 *  Tue Nov 17 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_galaxy_pz.c
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
 * NcGalaxyPz:
 *
 * Abstract class for the redshift probability density $p(z)$ of one galaxy.
 *
 * An object is data, not a model: it has no parameters and does not change after
 * construction, so one object may be shared by several catalogs and read concurrently.
 * The density is not normalized: nc_galaxy_pz_get_norm() gives its integral over the
 * range of nc_galaxy_pz_get_range(), and nc_galaxy_pz_gen() draws from the normalized
 * density. The evaluation functions check nothing: the caller keeps $z$ inside the
 * range.
 *
 * A subclass implements all five virtual functions.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc/lss/galaxy/nc_galaxy_pz.h"

G_DEFINE_ABSTRACT_TYPE (NcGalaxyPz, nc_galaxy_pz, G_TYPE_OBJECT)

static void
nc_galaxy_pz_init (NcGalaxyPz *gpz)
{
}

static void
nc_galaxy_pz_class_init (NcGalaxyPzClass *klass)
{
}

/**
 * nc_galaxy_pz_ref:
 * @gpz: a #NcGalaxyPz
 *
 * Increases the reference count of @gpz by one.
 *
 * Returns: (transfer full): @gpz.
 */
NcGalaxyPz *
nc_galaxy_pz_ref (NcGalaxyPz *gpz)
{
  return g_object_ref (gpz);
}

/**
 * nc_galaxy_pz_free:
 * @gpz: a #NcGalaxyPz
 *
 * Decreases the reference count of @gpz by one.
 */
void
nc_galaxy_pz_free (NcGalaxyPz *gpz)
{
  g_object_unref (gpz);
}

/**
 * nc_galaxy_pz_clear:
 * @gpz: a #NcGalaxyPz
 *
 * Decreases the reference count of *@gpz by one and sets *@gpz to %NULL; does nothing
 * when *@gpz is %NULL.
 */
void
nc_galaxy_pz_clear (NcGalaxyPz **gpz)
{
  g_clear_object (gpz);
}

/**
 * nc_galaxy_pz_get_range:
 * @gpz: a #NcGalaxyPz
 * @z_min: (out): the lower end of the range
 * @z_max: (out): the upper end of the range
 *
 * Gets the redshift range where the density is defined.
 */
void
nc_galaxy_pz_get_range (NcGalaxyPz *gpz, gdouble *z_min, gdouble *z_max)
{
  NC_GALAXY_PZ_GET_CLASS (gpz)->get_range (gpz, z_min, z_max);
}

/**
 * nc_galaxy_pz_eval_p:
 * @gpz: a #NcGalaxyPz
 * @z: the redshift
 *
 * Evaluates the density at @z, not normalized. The redshift is not checked against the
 * range.
 *
 * Returns: $p(z)$.
 */
gdouble
nc_galaxy_pz_eval_p (NcGalaxyPz *gpz, const gdouble z)
{
  return NC_GALAXY_PZ_GET_CLASS (gpz)->eval_p (gpz, z);
}

/**
 * nc_galaxy_pz_eval_ln_p:
 * @gpz: a #NcGalaxyPz
 * @z: the redshift
 *
 * Evaluates the logarithm of the density at @z, not normalized. The redshift is not
 * checked against the range.
 *
 * Returns: $\ln p(z)$.
 */
gdouble
nc_galaxy_pz_eval_ln_p (NcGalaxyPz *gpz, const gdouble z)
{
  return NC_GALAXY_PZ_GET_CLASS (gpz)->eval_ln_p (gpz, z);
}

/**
 * nc_galaxy_pz_get_norm:
 * @gpz: a #NcGalaxyPz
 *
 * Gets the integral of the density over its range.
 *
 * Returns: the normalization.
 */
gdouble
nc_galaxy_pz_get_norm (NcGalaxyPz *gpz)
{
  return NC_GALAXY_PZ_GET_CLASS (gpz)->get_norm (gpz);
}

/**
 * nc_galaxy_pz_gen:
 * @gpz: a #NcGalaxyPz
 * @rng: a #NcmRNG
 *
 * Draws a redshift from the normalized density with @rng.
 *
 * Returns: the redshift drawn.
 */
gdouble
nc_galaxy_pz_gen (NcGalaxyPz *gpz, NcmRNG *rng)
{
  return NC_GALAXY_PZ_GET_CLASS (gpz)->gen (gpz, rng);
}


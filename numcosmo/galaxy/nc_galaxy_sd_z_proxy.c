/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 2; tab-width: 2 -*-  */

/***************************************************************************
 *            nc_galaxy_sd_z_proxy.c
 *
 *  Sat May 20 23:08:38 2023
 *  Copyright  2023  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * nc_galaxy_sd_z_proxy.c
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
 * SECTION: nc_galaxy_sd_z_proxy
 * @title: NcGalaxySDZProxy
 * @short_description: Class describing galaxy sample position distributions.
 * @stability: Unstable
 *
 *
 * This class describes a galaxy sample proxy redshift distribution.
 * It is composed by a $P(z_p)$ distribution.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "galaxy/nc_galaxy_sd_z_proxy.h"
#include <math.h>
#include <gsl/gsl_math.h>


struct _NcGalaxySDZProxyPrivate
{
  gint placeholder;
};

G_DEFINE_TYPE_WITH_PRIVATE (NcGalaxySDZProxy, nc_galaxy_sd_z_proxy, G_TYPE_OBJECT);

static void
nc_galaxy_sd_z_proxy_init (NcGalaxySDZProxy *gsdzp)
{
  gsdzp->priv = nc_galaxy_sd_z_proxy_get_instance_private (gsdzp);
}

static void
_nc_galaxy_sd_z_proxy_finalize (GObject *object)
{
  G_OBJECT_CLASS (nc_galaxy_sd_z_proxy_parent_class)->finalize (object);
}

static gboolean
_nc_galaxy_sd_z_proxy_gen (NcGalaxySDZProxy *gsdzp, NcmRNG *rng, const gdouble z, gdouble *gen_zp)
{
  g_error ("_nc_galaxy_sd_z_proxy_gen: method not implemented.");

  return FALSE;
}

static gdouble
_nc_galaxy_sd_z_proxy_integ (NcGalaxySDZProxy *gsdzp, const gdouble z)
{
  g_error ("_nc_galaxy_sd_z_proxy_integ: method not implemented.");

  return 0.0;
}

static void
nc_galaxy_sd_z_proxy_class_init (NcGalaxySDZProxyClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->finalize = &_nc_galaxy_sd_z_proxy_finalize;

  klass->gen   = &_nc_galaxy_sd_z_proxy_gen;
  klass->integ = &_nc_galaxy_sd_z_proxy_integ;
}

/**
 * nc_galaxy_sd_z_proxy_ref:
 * @gsdzp: a #NcGalaxySDZProxy
 *
 * Increase the reference of @gsdzp by one.
 *
 * Returns: (transfer full): @gsdzp.
 */
NcGalaxySDZProxy *
nc_galaxy_sd_z_proxy_ref (NcGalaxySDZProxy *gsdzp)
{
  return g_object_ref (gsdzp);
}

/**
 * nc_galaxy_sd_z_proxy_free:
 * @gsdzp: a #NcGalaxySDZProxy
 *
 * Decrease the reference count of @gsdzp by one.
 *
 */
void
nc_galaxy_sd_z_proxy_free (NcGalaxySDZProxy *gsdzp)
{
  g_object_unref (gsdzp);
}

/**
 * nc_galaxy_sd_z_proxy_clear:
 * @gsdzp: a #NcGalaxySDZProxy
 *
 * Decrease the reference count of @gsdzp by one, and sets the pointer *@gsdzp to
 * NULL.
 *
 */
void
nc_galaxy_sd_z_proxy_clear (NcGalaxySDZProxy **gsdzp)
{
  g_clear_object (gsdzp);
}

/**
 * nc_galaxy_sd_z_proxy_gen: (virtual gen)
 * @gsdzp: a #NcGalaxySDZProxy
 * @rng: a #NcmRNG
 * @z: source redshift $z$
 *
 * Generates a $z_p$ value from the distribution using @rng.
 *
 * Returns: the generated value $z_p$.
 */
/**
 * nc_galaxy_sd_z_proxy_integ: (virtual integ)
 * @gsdzp: a #NcGalaxySDZProxy
 * @z: source redshift $z$
 *
 * Computes the probability density of the observable $z_p$ given the redshift.
 * The probability density is given by $P(z_p)P$.
 *
 * Returns: the probability density at $z_p$, $P(z_p)$.
 */


/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 2; tab-width: 2 -*-  */

/***************************************************************************
 *            nc_galaxy_sd_z_proxy_dirac.c
 *
 *  Wed June 21 19:40:56 2023
 *  Copyright  2023  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * nc_galaxy_sd_z_proxy_dirac.c
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
 * SECTION:nc_galaxy_sd_z_proxy_dirac
 * @title: NcGalaxySDZProxyDirac
 * @short_description: Class describing galaxy sample proxy redshift distributions with a dirac delta distribution
 * @stability: Unstable
 *
 *
 * Class defining a galaxy sample proxy redshift distribution with a dirac delta
 * probability distribution $P(z_p)$.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

/* #include "nc_enum_types.h" */
#include "galaxy/nc_galaxy_sd_z_proxy_dirac.h"
#include "galaxy/nc_galaxy_sd_z_proxy.h"
#include <math.h>
#include <gsl/gsl_math.h>

struct _NcGalaxySDZProxyDiracPrivate
{
  gint placeholder;
};

enum
{
  PROP_0,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcGalaxySDZProxyDirac, nc_galaxy_sd_z_proxy_dirac, NC_TYPE_GALAXY_SD_Z_PROXY);

static void
nc_galaxy_sd_z_proxy_dirac_init (NcGalaxySDZProxyDirac *gsdzpdirac)
{
}

static void
_nc_galaxy_sd_z_proxy_dirac_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  g_return_if_fail (NC_IS_GALAXY_SD_Z_PROXY_DIRAC (object));

  switch (prop_id)
  {
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_galaxy_sd_z_proxy_dirac_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  g_return_if_fail (NC_IS_GALAXY_SD_Z_PROXY_DIRAC (object));

  switch (prop_id)
  {
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_galaxy_sd_z_proxy_dirac_dispose (GObject *object)
{
  G_OBJECT_CLASS (nc_galaxy_sd_z_proxy_dirac_parent_class)->dispose (object);
}

static void
_nc_galaxy_sd_z_proxy_dirac_finalize (GObject *object)
{
  G_OBJECT_CLASS (nc_galaxy_sd_z_proxy_dirac_parent_class)->finalize (object);
}

static gboolean _nc_galaxy_sd_z_proxy_dirac_gen (NcGalaxySDZProxy *gsdp, NcmRNG *rng, const gdouble z, gdouble *gen_zp);
static gdouble _nc_galaxy_sd_z_proxy_dirac_integ (NcGalaxySDZProxy *gsdp, const gdouble z);

static void
nc_galaxy_sd_z_proxy_dirac_class_init (NcGalaxySDZProxyDiracClass *klass)
{
  NcGalaxySDZProxyClass *sd_position_class = NC_GALAXY_SD_Z_PROXY_CLASS (klass);
  GObjectClass *object_class               = G_OBJECT_CLASS (klass);

  object_class->set_property = &_nc_galaxy_sd_z_proxy_dirac_set_property;
  object_class->get_property = &_nc_galaxy_sd_z_proxy_dirac_get_property;
  object_class->dispose      = &_nc_galaxy_sd_z_proxy_dirac_dispose;
  object_class->finalize     = &_nc_galaxy_sd_z_proxy_dirac_finalize;

  sd_position_class->gen   = &_nc_galaxy_sd_z_proxy_dirac_gen;
  sd_position_class->integ = &_nc_galaxy_sd_z_proxy_dirac_integ;
}

static gboolean
_nc_galaxy_sd_z_proxy_dirac_gen (NcGalaxySDZProxy *gsdp, NcmRNG *rng, const gdouble z, gdouble *gen_zp)
{
  gen_zp[0] = z;

  return TRUE;
}

static gdouble
_nc_galaxy_sd_z_proxy_dirac_integ (NcGalaxySDZProxy *gsdp, gdouble z)
{
  return 0.0;
}

/**
 * nc_galaxy_sd_z_proxy_dirac_new:
 *
 * Creates a new #NcGalaxySDZProxyDirac
 *
 * Returns: (transfer full): a new NcGalaxySDZProxyDirac.
 */
NcGalaxySDZProxyDirac *
nc_galaxy_sd_z_proxy_dirac_new ()
{
  NcGalaxySDZProxyDirac *gsdzpdirac = g_object_new (NC_TYPE_GALAXY_SD_Z_PROXY_DIRAC,
                                                    NULL);

  return gsdzpdirac;
}

/**
 * nc_galaxy_sd_z_proxy_dirac_ref:
 * @gsdzpdirac: a #NcGalaxySDZProxyDirac
 *
 * Increase the reference of @gsdzpdirac by one.
 *
 * Returns: (transfer full): @gsdzpdirac.
 */
NcGalaxySDZProxyDirac *
nc_galaxy_sd_z_proxy_dirac_ref (NcGalaxySDZProxyDirac *gsdzpdirac)
{
  return g_object_ref (gsdzpdirac);
}

/**
 * nc_galaxy_sd_z_proxy_dirac_free:
 * @gsdzpdirac: a #NcGalaxySDZProxyDirac
 *
 * Decrease the reference count of @gsdzpdirac by one.
 *
 */
void
nc_galaxy_sd_z_proxy_dirac_free (NcGalaxySDZProxyDirac *gsdzpdirac)
{
  g_object_unref (gsdzpdirac);
}

/**
 * nc_galaxy_sd_z_proxy_dirac_clear:
 * @gsdzpdirac: a #NcGalaxySDZProxyDirac
 *
 * Decrease the reference count of @gsdzpdirac by one, and sets the pointer *@gsdzpdirac to
 * NULL.
 *
 */
void
nc_galaxy_sd_z_proxy_dirac_clear (NcGalaxySDZProxyDirac **gsdzpdirac)
{
  g_clear_object (gsdzpdirac);
}


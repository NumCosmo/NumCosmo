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


typedef struct _NcGalaxySDZProxyPrivate
{
  gint placeholder;
} NcGalaxySDZProxyPrivate;

G_DEFINE_TYPE_WITH_PRIVATE (NcGalaxySDZProxy, nc_galaxy_sd_z_proxy, NCM_TYPE_MODEL);

static void
nc_galaxy_sd_z_proxy_init (NcGalaxySDZProxy *gsdzp)
{
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
_nc_galaxy_sd_z_proxy_integ (NcGalaxySDZProxy *gsdzp, const gdouble z, const gdouble zp)
{
  g_error ("_nc_galaxy_sd_z_proxy_integ: method not implemented.");

  return 0.0;
}

static void
_nc_galaxy_sd_z_proxy_get_true_z_lim (NcGalaxySDZProxy *gsdzp, const gdouble zp, gdouble *z_min, gdouble *z_max)
{
  g_error ("_nc_galaxy_sd_z_proxy_get_true_z_lim: method not implemented.");
}

static void
nc_galaxy_sd_z_proxy_class_init (NcGalaxySDZProxyClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);
  NcmModelClass *model_class = NCM_MODEL_CLASS (klass);

  object_class->finalize = &_nc_galaxy_sd_z_proxy_finalize;

  ncm_model_class_set_name_nick (model_class, "Galaxy redshift proxy distribution", "NcGalaxySDZProxy");
  ncm_model_class_add_params (model_class, 0, 0, 1);

  ncm_mset_model_register_id (model_class,
                              "NcGalaxySDZProxy",
                              "Galaxy redshift proxy distribution.",
                              NULL,
                              TRUE,
                              NCM_MSET_MODEL_MAIN);

  ncm_model_class_check_params_info (NCM_MODEL_CLASS (klass));

  klass->gen            = &_nc_galaxy_sd_z_proxy_gen;
  klass->integ          = &_nc_galaxy_sd_z_proxy_integ;
  klass->get_true_z_lim = &_nc_galaxy_sd_z_proxy_get_true_z_lim;
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
 * @gen_zp: (out): generated proxy redshift $z_p$
 *
 * Generates a $z_p$ value from the distribution using @rng.
 *
 * Returns: the generated value $z_p$.
 */

gboolean
nc_galaxy_sd_z_proxy_gen (NcGalaxySDZProxy *gsdzp, NcmRNG *rng, const gdouble z, gdouble *gen_zp)
{
  return NC_GALAXY_SD_Z_PROXY_GET_CLASS (gsdzp)->gen (gsdzp, rng, z, gen_zp);
}

/**
 * nc_galaxy_sd_z_proxy_integ: (virtual integ)
 * @gsdzp: a #NcGalaxySDZProxy
 * @z: source redshift $z$
 * @zp: proxy redshift $z_p$
 *
 * Computes the probability density of the observable $z_p$ given the redshift.
 * The probability density is given by $P(z_p)P$.
 *
 * Returns: the probability density at $z_p$, $P(z_p)$.
 */
gdouble
nc_galaxy_sd_z_proxy_integ (NcGalaxySDZProxy *gsdzp, const gdouble z, const gdouble zp)
{
  return NC_GALAXY_SD_Z_PROXY_GET_CLASS (gsdzp)->integ (gsdzp, z, zp);
}

/**
 * nc_galaxy_sd_z_proxy_get_true_z_lim: (virtual get_true_z_lim)
 * @gsdzp: a #NcGalaxySDZProxy
 * @zp: proxy redshift $z_p$
 * @z_min: (out): minimum redshift
 * @z_max: (out): maximum redshift
 *
 * Computes the true redshift limits for a given proxy redshift $z_p$.
 *
 */
void
nc_galaxy_sd_z_proxy_get_true_z_lim (NcGalaxySDZProxy *gsdzp, const gdouble zp, gdouble *z_min, gdouble *z_max)
{
  NC_GALAXY_SD_Z_PROXY_GET_CLASS (gsdzp)->get_true_z_lim (gsdzp, zp, z_min, z_max);
}


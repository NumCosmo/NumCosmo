/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 2; tab-width: 2 -*-  */

/***************************************************************************
 *            nc_galaxy_sd_z_proxy_gauss.c
 *
 *  Tue June 1 17:18:06 2023
 *  Copyright  2023  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * nc_galaxy_sd_z_proxy_gauss.c
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
 * SECTION:nc_galaxy_sd_z_proxy_gauss
 * @title: NcGalaxySDZProxyGauss
 * @short_description: Class describing galaxy sample proxy redshift distributions with gaussian distribution
 * @stability: Unstable
 *
 *
 * Class defining a galaxy sample proxy redshift distribution with gaussian
 * probability distribution $P(z_p)$.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

// #include "nc_enum_types.h"
#include "galaxy/nc_galaxy_sd_z_proxy_gauss.h"
#include "galaxy/nc_galaxy_sd_z_proxy.h"
#include <math.h>
#include <gsl/gsl_math.h>

struct _NcGalaxySDZProxyGaussPrivate
{
  NcmVector *z_lim;
  gdouble sigma;
};

enum
{
  PROP_0,
  PROP_Z_LIM,
  PROP_SIGMA,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcGalaxySDZProxyGauss, nc_galaxy_sd_z_proxy_gauss, NC_TYPE_GALAXY_SD_Z_PROXY);

static void
nc_galaxy_sd_z_proxy_gauss_init (NcGalaxySDZProxyGauss *gsdzpgauss)
{
  NcGalaxySDZProxyGaussPrivate * const self = gsdzpgauss->priv = nc_galaxy_sd_z_proxy_gauss_get_instance_private (gsdzpgauss);

  self->z_lim = NULL;
  self->sigma = 0.0;
}

static void
_nc_galaxy_sd_z_proxy_gauss_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcGalaxySDZProxyGauss *gsdzpgauss = NC_GALAXY_SD_Z_PROXY_GAUSS (object);

  g_return_if_fail (NC_IS_GALAXY_SD_Z_PROXY_GAUSS (object));

  switch (prop_id)
  {
    case PROP_Z_LIM:
      nc_galaxy_sd_z_proxy_gauss_set_z_lim (gsdzpgauss, g_value_get_object (value));
      break;
    case PROP_SIGMA:
      nc_galaxy_sd_z_proxy_gauss_set_sigma (gsdzpgauss, g_value_get_double (value));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_galaxy_sd_z_proxy_gauss_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcGalaxySDZProxyGauss *gsdzpgauss = NC_GALAXY_SD_Z_PROXY_GAUSS (object);
  NcGalaxySDZProxyGaussPrivate * const self = gsdzpgauss->priv;

  g_return_if_fail (NC_IS_GALAXY_SD_Z_PROXY_GAUSS (object));

  switch (prop_id)
  {
    case PROP_Z_LIM:
      g_value_set_object (value, nc_galaxy_sd_z_proxy_gauss_peek_z_lim (gsdzpgauss));
      break;
    case PROP_SIGMA:
      g_value_set_double (value, self->sigma);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}


static void
_nc_galaxy_sd_z_proxy_gauss_dispose (GObject *object)
{
  NcGalaxySDZProxyGauss *gsdzpgauss = NC_GALAXY_SD_Z_PROXY_GAUSS (object);
  NcGalaxySDZProxyGaussPrivate * const self = gsdzpgauss->priv;

  ncm_vector_clear (&self->z_lim);

  G_OBJECT_CLASS (nc_galaxy_sd_z_proxy_gauss_parent_class)->dispose (object);
}

static void
_nc_galaxy_sd_z_proxy_gauss_finalize (GObject *object)
{
  G_OBJECT_CLASS (nc_galaxy_sd_z_proxy_gauss_parent_class)->finalize (object);
}

static gdouble _nc_galaxy_sd_z_proxy_gauss_gen (NcGalaxySDZProxy *gsdp, NcmRNG *rng, const gdouble z);
static gdouble _nc_galaxy_sd_z_proxy_gauss_integ (NcGalaxySDZProxy *gsdp, const gdouble z);

static void
nc_galaxy_sd_z_proxy_gauss_class_init (NcGalaxySDZProxyGaussClass *klass)
{
  NcGalaxySDZProxyClass *sd_position_class = NC_GALAXY_SD_Z_PROXY_CLASS (klass);
  GObjectClass *object_class         = G_OBJECT_CLASS (klass);

  object_class->set_property = &_nc_galaxy_sd_z_proxy_gauss_set_property;
  object_class->get_property = &_nc_galaxy_sd_z_proxy_gauss_get_property;
  object_class->dispose      = &_nc_galaxy_sd_z_proxy_gauss_dispose;
  object_class->finalize     = &_nc_galaxy_sd_z_proxy_gauss_finalize;

  /**
   * NcGalaxySDZProxyGauss:z-lim:
   *
   * Galaxy sample redshift distribution limits.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_Z_LIM,
                                   g_param_spec_object ("z-lim",
                                                        NULL,
                                                        "Galaxy sample photometric redshift distribution limits",
                                                        NCM_TYPE_VECTOR,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));



  sd_position_class->gen     = &_nc_galaxy_sd_z_proxy_gauss_gen;
  sd_position_class->integ   = &_nc_galaxy_sd_z_proxy_gauss_integ;
}


static gdouble
_nc_galaxy_sd_z_proxy_gauss_gen (NcGalaxySDZProxy *gsdp, NcmRNG *rng, const gdouble z)
{
  NcGalaxySDZProxyGauss *gsdzpgauss = NC_GALAXY_SD_Z_PROXY_GAUSS (gsdp);
  NcGalaxySDZProxyGaussPrivate * const self = gsdzpgauss->priv;
  gdouble z_p_gen = 0.0;
  gdouble z_ll = ncm_vector_get (self->z_lim, 0);
  gdouble z_ul = ncm_vector_get (self->z_lim, 1);

  while (z_p_gen < z_ll || z_p_gen > z_ul)
  {
    z_p_gen = ncm_rng_gaussian_gen (rng, z, self->sigma * (1+z));
  }

  return z_p_gen;
}

static gdouble 
_nc_galaxy_sd_z_proxy_gauss_integ (NcGalaxySDZProxy *gsdp, gdouble z)
{
  return 0.0;
}

/**
 * nc_galaxy_sd_z_proxy_gauss_new:
 *
 * Creates a new #NcGalaxySDZProxyGauss
 *
 * Returns: (transfer full): a new NcGalaxySDZProxyGauss.
 */
NcGalaxySDZProxyGauss *
nc_galaxy_sd_z_proxy_gauss_new ()
{
  NcGalaxySDZProxyGauss *gsdzpgauss = g_object_new (NC_TYPE_GALAXY_SD_Z_PROXY_GAUSS,
                                                  NULL);

  return gsdzpgauss;
}

/**
 * nc_galaxy_sd_z_proxy_gauss_ref:
 * @gsdzpgauss: a #NcGalaxySDZProxyGauss
 *
 * Increase the reference of @gsdzpgauss by one.
 *
 * Returns: (transfer full): @gsdzpgauss.
 */
NcGalaxySDZProxyGauss *
nc_galaxy_sd_z_proxy_gauss_ref (NcGalaxySDZProxyGauss *gsdzpgauss)
{
  return g_object_ref (gsdzpgauss);
}

/**
 * nc_galaxy_sd_z_proxy_gauss_free:
 * @gsdzpgauss: a #NcGalaxySDZProxyGauss
 *
 * Decrease the reference count of @gsdzpgauss by one.
 *
 */
void
nc_galaxy_sd_z_proxy_gauss_free (NcGalaxySDZProxyGauss *gsdzpgauss)
{
  g_object_unref (gsdzpgauss);
}

/**
 * nc_galaxy_sd_z_proxy_gauss_clear:
 * @gsdzpgauss: a #NcGalaxySDZProxyGauss
 *
 * Decrease the reference count of @gsdzpgauss by one, and sets the pointer *@gsdzpgauss to
 * NULL.
 *
 */
void
nc_galaxy_sd_z_proxy_gauss_clear (NcGalaxySDZProxyGauss **gsdzpgauss)
{
  g_clear_object (gsdzpgauss);
}

/**
 * nc_galaxy_sd_z_proxy_gauss_set_z_lim:
 * @gsdzpgauss: a #NcGalaxySDZProxyGauss
 * @lim: a #NcmVector
 *
 * Sets the redshift limits @lim.
 */
void
nc_galaxy_sd_z_proxy_gauss_set_z_lim (NcGalaxySDZProxyGauss *gsdzpgauss, NcmVector *lim)
{
  NcGalaxySDZProxyGaussPrivate * const self = gsdzpgauss->priv;

  g_assert_cmpuint (ncm_vector_len (lim), ==, 2);

  ncm_vector_clear (&self->z_lim);

  self->z_lim = ncm_vector_ref (lim);
}


/**
 * nc_galaxy_sd_z_proxy_gauss_peek_z_lim:
 * @gsdzpgauss: a #NcGalaxySDZProxyGauss
 *
 * Gets the redshift limits.
 *
 * Returns: (transfer none): the redshift limits.
 */
NcmVector *
nc_galaxy_sd_z_proxy_gauss_peek_z_lim (NcGalaxySDZProxyGauss *gsdzpgauss)
{
  NcGalaxySDZProxyGaussPrivate * const self = gsdzpgauss->priv;

  return self->z_lim;
}

/**
 * nc_galaxy_sd_z_proxy_gauss_set_sigma:
 * @gsdzpgauss: a #NcGalaxySDZProxyGauss
 * @sigma: the standard deviation
 *
 * Sets the standard deviation $\sigma_0$.
 */
void
nc_galaxy_sd_z_proxy_gauss_set_sigma (NcGalaxySDZProxyGauss *gsdzpgauss, gdouble sigma)
{
  NcGalaxySDZProxyGaussPrivate * const self = gsdzpgauss->priv;

  self->sigma = sigma;
}


/**
 * nc_galaxy_sd_z_proxy_gauss_peek_sigma:
 * @gsdzpgauss: a #NcGalaxySDZProxyGauss
 *
 * Gets the standard deviation $\sigma_0$.
 *
 * Returns: the standard deviation $\sigma_0$.
 */
gdouble
nc_galaxy_sd_z_proxy_gauss_get_sigma (NcGalaxySDZProxyGauss *gsdzpgauss)
{
  NcGalaxySDZProxyGaussPrivate * const self = gsdzpgauss->priv;

  return self->sigma;
}
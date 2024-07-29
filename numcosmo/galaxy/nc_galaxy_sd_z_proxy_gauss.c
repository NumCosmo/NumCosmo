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

/* #include "nc_enum_types.h" */
#include "galaxy/nc_galaxy_sd_z_proxy_gauss.h"
#include "galaxy/nc_galaxy_sd_z_proxy.h"
#include "math/ncm_dtuple.h"
#include <math.h>
#include <gsl/gsl_math.h>

typedef struct _NcGalaxySDZProxyGaussPrivate
{
  gdouble true_z_min;
  gdouble z_min;
  gdouble z_max;
  gdouble sigma;
} NcGalaxySDZProxyGaussPrivate;

struct _NcGalaxySDZProxyGauss
{
  NcGalaxySDZProxy parent_instance;
};

enum
{
  PROP_0,
  PROP_TRUE_Z_MIN,
  PROP_Z_LIM,
  PROP_SIGMA,
  PROP_LEN,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcGalaxySDZProxyGauss, nc_galaxy_sd_z_proxy_gauss, NC_TYPE_GALAXY_SD_Z_PROXY);

static void
nc_galaxy_sd_z_proxy_gauss_init (NcGalaxySDZProxyGauss *gsdzpgauss)
{
  NcGalaxySDZProxyGaussPrivate * const self = nc_galaxy_sd_z_proxy_gauss_get_instance_private (gsdzpgauss);

  self->true_z_min = 1e-11;
  self->z_min      = 0.0;
  self->z_max      = 0.0;
  self->sigma      = 0.0;
}

static void
_nc_galaxy_sd_z_proxy_gauss_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcGalaxySDZProxyGauss *gsdzpgauss = NC_GALAXY_SD_Z_PROXY_GAUSS (object);

  g_return_if_fail (NC_IS_GALAXY_SD_Z_PROXY_GAUSS (object));

  switch (prop_id)
  {
    case PROP_TRUE_Z_MIN:
    {
      if (value == 0)
        g_error ("_nc_galaxy_sd_position_set_property: value is zero.");

      nc_galaxy_sd_z_proxy_gauss_set_true_z_min (gsdzpgauss, g_value_get_double (value));
      break;
    }
    case PROP_Z_LIM:
    {
      NcmDTuple2 *z_lim = g_value_get_boxed (value);

      if (z_lim == NULL)
        g_error ("_nc_galaxy_sd_position_set_property: z_lim is NULL.");

      nc_galaxy_sd_z_proxy_gauss_set_z_lim (gsdzpgauss, z_lim->elements[0], z_lim->elements[1]);
      break;
    }
    case PROP_SIGMA:
      nc_galaxy_sd_z_proxy_gauss_set_sigma (gsdzpgauss, g_value_get_double (value));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_galaxy_sd_z_proxy_gauss_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcGalaxySDZProxyGauss *gsdzpgauss         = NC_GALAXY_SD_Z_PROXY_GAUSS (object);
  NcGalaxySDZProxyGaussPrivate * const self = nc_galaxy_sd_z_proxy_gauss_get_instance_private (gsdzpgauss);

  g_return_if_fail (NC_IS_GALAXY_SD_Z_PROXY_GAUSS (object));

  switch (prop_id)
  {
    case PROP_TRUE_Z_MIN:
      g_value_set_double (value, self->true_z_min);
      break;
    case PROP_Z_LIM:
    {
      gdouble z_min, z_max;

      nc_galaxy_sd_z_proxy_gauss_get_z_lim (gsdzpgauss, &z_min, &z_max);

      g_value_take_boxed (value, ncm_dtuple2_new (z_min, z_max));
      break;
    }
    case PROP_SIGMA:
      g_value_set_double (value, self->sigma);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_galaxy_sd_z_proxy_gauss_dispose (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_galaxy_sd_z_proxy_gauss_parent_class)->dispose (object);
}

static void
_nc_galaxy_sd_z_proxy_gauss_finalize (GObject *object)
{
  G_OBJECT_CLASS (nc_galaxy_sd_z_proxy_gauss_parent_class)->finalize (object);
}

static gboolean _nc_galaxy_sd_z_proxy_gauss_gen (NcGalaxySDZProxy *gsdzp, NcmRNG *rng, const gdouble z, gdouble *gen_zp);
static gdouble _nc_galaxy_sd_z_proxy_gauss_integ (NcGalaxySDZProxy *gsdzp, const gdouble z, const gdouble zp);
static void _nc_galaxy_sd_z_proxy_gauss_get_true_z_lim (NcGalaxySDZProxy *gsdzp, const gdouble zp, gdouble *z_min, gdouble *z_max);

static void
nc_galaxy_sd_z_proxy_gauss_class_init (NcGalaxySDZProxyGaussClass *klass)
{
  NcGalaxySDZProxyClass *sd_position_class = NC_GALAXY_SD_Z_PROXY_CLASS (klass);
  NcmModelClass *model_class               = NCM_MODEL_CLASS (klass);
  GObjectClass *object_class               = G_OBJECT_CLASS (klass);

  model_class->set_property = &_nc_galaxy_sd_z_proxy_gauss_set_property;
  model_class->get_property = &_nc_galaxy_sd_z_proxy_gauss_get_property;
  object_class->dispose     = &_nc_galaxy_sd_z_proxy_gauss_dispose;
  object_class->finalize    = &_nc_galaxy_sd_z_proxy_gauss_finalize;

  ncm_model_class_set_name_nick (model_class, "Gaussian photometric redshift distribution", "GAUSS_PHOTO_Z");
  ncm_model_class_add_params (model_class, 0, 0, PROP_LEN);

  /**
   * NcGalaxySDZProxyGauss:true-z-min:
   *
   * Minimum true redshift.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_TRUE_Z_MIN,
                                   g_param_spec_double ("true-z-min",
                                                        NULL,
                                                        "Minimum true redshift",
                                                        1e-11, G_MAXDOUBLE, 1e-11,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcGalaxySDZProxyGauss:z-lim:
   *
   * Galaxy sample redshift distribution limits.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_Z_LIM,
                                   g_param_spec_boxed ("z-lim",
                                                       NULL,
                                                       "Galaxy sample photometric redshift distribution limits",
                                                       NCM_TYPE_DTUPLE2,
                                                       G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcGalaxySDZProxyGauss:sigma:
   *
   * Proxy redshift standard deviation.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_SIGMA,
                                   g_param_spec_double ("sigma",
                                                        NULL,
                                                        "Proxy redshift standard deviation",
                                                        1.0e-4, G_MAXDOUBLE, 5.0e-2,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  ncm_model_class_check_params_info (model_class);

  sd_position_class->gen            = &_nc_galaxy_sd_z_proxy_gauss_gen;
  sd_position_class->integ          = &_nc_galaxy_sd_z_proxy_gauss_integ;
  sd_position_class->get_true_z_lim = &_nc_galaxy_sd_z_proxy_gauss_get_true_z_lim;
}

static gboolean
_nc_galaxy_sd_z_proxy_gauss_gen (NcGalaxySDZProxy *gsdzp, NcmRNG *rng, const gdouble z, gdouble *gen_zp)
{
  NcGalaxySDZProxyGauss *gsdzpgauss         = NC_GALAXY_SD_Z_PROXY_GAUSS (gsdzp);
  NcGalaxySDZProxyGaussPrivate * const self = nc_galaxy_sd_z_proxy_gauss_get_instance_private (gsdzpgauss);

  gen_zp[0] = ncm_rng_gaussian_gen (rng, z, self->sigma * (1.0 + z));

  return (gen_zp[0] > self->z_min) && (gen_zp[0] < self->z_max);
}

static gdouble
_nc_galaxy_sd_z_proxy_gauss_integ (NcGalaxySDZProxy *gsdzp, const gdouble z, const gdouble zp)
{
  NcGalaxySDZProxyGauss *gsdzpgauss         = NC_GALAXY_SD_Z_PROXY_GAUSS (gsdzp);
  NcGalaxySDZProxyGaussPrivate * const self = nc_galaxy_sd_z_proxy_gauss_get_instance_private (gsdzpgauss);

  return exp (-pow (zp - z, 2.0) / (2.0 * pow (self->sigma * (1.0 + z), 2.0))) / (self->sigma * (1.0 + z) * sqrt (2.0 * M_PI));
  /* return gsl_ran_gaussian_pdf (zp - z, self->sigma * (1.0 + z)); */
}

static void
_nc_galaxy_sd_z_proxy_gauss_get_true_z_lim (NcGalaxySDZProxy *gsdzp, const gdouble zp, gdouble *z_min, gdouble *z_max)
{
  NcGalaxySDZProxyGauss *gsdzpgauss         = NC_GALAXY_SD_Z_PROXY_GAUSS (gsdzp);
  NcGalaxySDZProxyGaussPrivate * const self = nc_galaxy_sd_z_proxy_gauss_get_instance_private (gsdzpgauss);
  const gdouble sigma                       = self->sigma * (1.0 + zp);
  const gdouble max_z                       = zp + 5.0 * sigma;
  const gdouble sigma_max                   = self->sigma * (1.0 + max_z);

  g_assert_nonnull (z_min);
  g_assert_nonnull (z_max);

  *z_min = MAX (zp - 8.0 * sigma_max, self->true_z_min);
  *z_max = zp + 8.0 * sigma_max;
}

/**
 * nc_galaxy_sd_z_proxy_gauss_new:
 * @z_min: the minimum photometric redshift
 * @z_max: the maximum photometric redshift
 * @sigma: the standard deviation
 *
 * Creates a new #NcGalaxySDZProxyGauss
 *
 * Returns: (transfer full): a new NcGalaxySDZProxyGauss.
 */
NcGalaxySDZProxyGauss *
nc_galaxy_sd_z_proxy_gauss_new (const gdouble z_min, const gdouble z_max, const gdouble sigma)
{
  NcmDTuple2 z_lim                  = NCM_DTUPLE2_STATIC_INIT (z_min, z_max);
  NcGalaxySDZProxyGauss *gsdzpgauss = g_object_new (NC_TYPE_GALAXY_SD_Z_PROXY_GAUSS,
                                                    "z-lim", &z_lim,
                                                    "sigma", sigma,
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
 * nc_galaxy_sd_z_proxy_gauss_set_true_z_min:
 * @gsdzpgauss: a #NcGalaxySDZProxyGauss
 * @true_z_min: the minimum true redshift
 *
 * Sets the minimum true redshift.
 *
 */
void
nc_galaxy_sd_z_proxy_gauss_set_true_z_min (NcGalaxySDZProxyGauss *gsdzpgauss, const gdouble true_z_min)
{
  NcGalaxySDZProxyGaussPrivate * const self = nc_galaxy_sd_z_proxy_gauss_get_instance_private (gsdzpgauss);

  self->true_z_min = true_z_min;
}

/**
 * nc_galaxy_sd_z_proxy_gauss_get_true_z_min:
 * @gsdzpgauss: a #NcGalaxySDZProxyGauss
 *
 * Gets the minimum true redshift.
 *
 * Returns: the minimum true redshift.
 */
gdouble
nc_galaxy_sd_z_proxy_gauss_get_true_z_min (NcGalaxySDZProxyGauss *gsdzpgauss)
{
  NcGalaxySDZProxyGaussPrivate * const self = nc_galaxy_sd_z_proxy_gauss_get_instance_private (gsdzpgauss);

  return self->true_z_min;
}

/**
 * nc_galaxy_sd_z_proxy_gauss_set_z_lim:
 * @gsdzpgauss: a #NcGalaxySDZProxyGauss
 * @z_min: the minimum photometric redshift
 * @z_max: the maximum photometric redshift
 *
 * Sets the photometric redshift limits.
 *
 */
void
nc_galaxy_sd_z_proxy_gauss_set_z_lim (NcGalaxySDZProxyGauss *gsdzpgauss, const gdouble z_min, const gdouble z_max)
{
  NcGalaxySDZProxyGaussPrivate * const self = nc_galaxy_sd_z_proxy_gauss_get_instance_private (gsdzpgauss);

  g_assert_cmpfloat (z_min, <, z_max);

  self->z_min = z_min;
  self->z_max = z_max;
}

/**
 * nc_galaxy_sd_z_proxy_gauss_get_z_lim:
 * @gsdzpgauss: a #NcGalaxySDZProxyGauss
 * @z_min: (out): the minimum photometric redshift
 * @z_max: (out): the maximum photometric redshift
 *
 * Gets the photometric redshift limits.
 *
 */
void
nc_galaxy_sd_z_proxy_gauss_get_z_lim (NcGalaxySDZProxyGauss *gsdzpgauss, gdouble *z_min, gdouble *z_max)
{
  NcGalaxySDZProxyGaussPrivate * const self = nc_galaxy_sd_z_proxy_gauss_get_instance_private (gsdzpgauss);

  g_assert_nonnull (z_min);
  g_assert_nonnull (z_max);

  z_min[0] = self->z_min;
  z_max[0] = self->z_max;
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
  NcGalaxySDZProxyGaussPrivate * const self = nc_galaxy_sd_z_proxy_gauss_get_instance_private (gsdzpgauss);

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
  NcGalaxySDZProxyGaussPrivate * const self = nc_galaxy_sd_z_proxy_gauss_get_instance_private (gsdzpgauss);

  return self->sigma;
}


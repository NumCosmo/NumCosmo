/***************************************************************************
 *            nc_halo_density_profile_dk14.c
 *
 *  Tue July 16 15:20:15 2019
 *  Copyright  2014
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * nc_halo_density_profile_dk14.c
 * Copyright (C) 2014 Mariana Penna Lima <pennalima@gmail.com>
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
 * SECTION:nc_halo_density_profile_dk14
 * @title: NcHaloDensityProfileDK14
 * @short_description: Density profile of Diemer \& Kravtsov type.
 *
 * This object implements the #NcHaloDensityProfile class for a Diemer \& Kravtsov (DK14) density profile.
 *
 * The DK14 profile is defined as ... Implementation incomplete! FIXME
 *
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "lss/nc_halo_density_profile_dk14.h"
#include "math/ncm_cfg.h"
#include "math/ncm_util.h"
#include "math/ncm_integrate.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_sf_expint.h>
#endif /* NUMCOSMO_GIR_SCAN */

G_DEFINE_TYPE (NcHaloDensityProfileDK14, nc_halo_density_profile_dk14, NC_TYPE_HALO_DENSITY_PROFILE)

#define VECTOR (NCM_MODEL (dpdk))
#define RT     (ncm_model_orig_param_get (VECTOR, NC_HALO_DENSITY_PROFILE_DK14_RT))
#define BETA   (ncm_model_orig_param_get (VECTOR, NC_HALO_DENSITY_PROFILE_DK14_BETA))

enum
{
  PROP_0,
  PROP_DELTA,
  PROP_SIZE,
};

static void
nc_halo_density_profile_dk14_init (NcHaloDensityProfileDK14 *dpdk)
{
  dpdk->Delta   = 200.0;
  dpdk->r_Delta = 0.0;
}

static void
_nc_halo_density_profile_dk14_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcHaloDensityProfileDK14 *dpdk = NC_HALO_DENSITY_PROFILE_DK14 (object);

  g_return_if_fail (NC_IS_HALO_DENSITY_PROFILE_DK14 (object));

  switch (prop_id)
  {
    case PROP_DELTA:
      dpdk->Delta = g_value_get_double (value);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_halo_density_profile_dk14_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcHaloDensityProfileDK14 *dpdk = NC_HALO_DENSITY_PROFILE_DK14 (object);

  g_return_if_fail (NC_IS_HALO_DENSITY_PROFILE_DK14 (object));

  switch (prop_id)
  {
    case PROP_DELTA:
      g_value_set_double (value, dpdk->Delta);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_halo_density_profile_dk14_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_halo_density_profile_dk14_parent_class)->finalize (object);
}

static void
nc_halo_density_profile_dk14_class_init (NcHaloDensityProfileDK14Class *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);
  /*NcHaloDensityProfileClass *parent_class = NC_HALO_DENSITY_PROFILE_CLASS (klass);*/
  NcmModelClass *model_class = NCM_MODEL_CLASS (klass);

  model_class->set_property = &_nc_halo_density_profile_dk14_set_property;
  model_class->get_property = &_nc_halo_density_profile_dk14_get_property;
  object_class->finalize    = &nc_halo_density_profile_dk14_finalize;

  ncm_model_class_set_name_nick (model_class, "DK14 Density Profile", "DK14");
  ncm_model_class_add_params (model_class, NC_HALO_DENSITY_PROFILE_DK14_SPARAM_LEN, 0, PROP_SIZE);

  /**
   * NcHaloDensityProfileDK14:Delta:
   *
   * Constant that indicates the overdensity with respect to the critical density.
   * FIXME Set correct values (limits)
   */
  g_object_class_install_property (object_class,
                                   PROP_DELTA,
                                   g_param_spec_double ("Delta",
                                                        NULL,
                                                        "Overdensity constant",
                                                        200, 1500, 200,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcHaloDensityProfileDK14:rt:
   *
   * Truncation radius.
   * FIXME Set correct values (limits)
   */
  ncm_model_class_set_sparam (model_class, NC_HALO_DENSITY_PROFILE_DK14_RT, "r_{t}", "rt",
                              0.5,  10.0, 1.0e-1,
                              NC_HALO_DENSITY_PROFILE_DK14_DEFAULT_PARAMS_ABSTOL, NC_HALO_DENSITY_PROFILE_DK14_DEFAULT_RT,
                              NCM_PARAM_TYPE_FIXED);

  /**
   * NcHaloDensityProfileDK14:beta:
   *
   * Sharpness of the steepening.
   * FIXME Set correct values (limits)
   */
  ncm_model_class_set_sparam (model_class, NC_HALO_DENSITY_PROFILE_DK14_BETA, "\\beta", "beta",
                              1.0,  10.0, 4.0,
                              NC_HALO_DENSITY_PROFILE_DK14_DEFAULT_PARAMS_ABSTOL, NC_HALO_DENSITY_PROFILE_DK14_DEFAULT_BETA,
                              NCM_PARAM_TYPE_FIXED);
}

/**
 * nc_halo_density_profile_dk14_new:
 *
 * This function returns a #NcHaloDensityProfile with a #NcHaloDensityProfileDK14 implementation.
 *
 * Returns: A new #NcHaloDensityProfile.
 */
NcHaloDensityProfile *
nc_halo_density_profile_dk14_new ()
{
  return g_object_new (NC_TYPE_HALO_DENSITY_PROFILE_DK14, NULL);
}


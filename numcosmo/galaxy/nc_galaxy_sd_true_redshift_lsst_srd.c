/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 2; tab-width: 2 -*-  */

/***************************************************************************
 *            nc_galaxy_sd_true_redshift_lsst_srd.c
 *
 *  Wed Jul 31 21:45:13 2024
 *  Copyright  2024  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * nc_galaxy_sd_true_redshift_lsst_srd.c
 * Copyright (C) 2024 Caio Lima de Oliveira <caiolimadeoliveira@pm.me>
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
 * SECTION:nc_galaxy_sd_true_redshift_lsst_srd
 * @title: NcGalaxySDTrueRedshiftLSSTSRD
 * @short_description: Class describing galaxy sample redshift distributions as in LSST-SRD.
 * @stability: Unstable
 *
 *
 * Class defining a galaxy sample redshift distribution as described in the LSST
 * Science Roadmap Document,
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math.h"
#include "gsl/gsl_sf.h"
#include "galaxy/nc_galaxy_sd_true_redshift_lsst_srd.h"
#include "galaxy/nc_galaxy_sd_true_redshift.h"
#include "math/ncm_dtuple.h"
#include "math/ncm_model.h"
#include "math/ncm_rng.h"
#include "math/ncm_vector.h"

typedef struct _NcGalaxySDTrueRedshiftLSSTSRDPrivate
{
  gdouble z_min;
  gdouble z_max;
  gdouble z_norm;
  gdouble y0;
} NcGalaxySDTrueRedshiftLSSTSRDPrivate;

struct _NcGalaxySDTrueRedshiftLSSTSRD
{
  NcGalaxySDTrueRedshift parent_instance;
};

enum
{
  PROP_0,
  PROP_LEN,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcGalaxySDTrueRedshiftLSSTSRD, nc_galaxy_sd_true_redshift_lsst_srd, NC_TYPE_GALAXY_SD_TRUE_REDSHIFT);

static void
nc_galaxy_sd_true_redshift_lsst_srd_init (NcGalaxySDTrueRedshiftLSSTSRD *gsdtrlsst)
{
  NcGalaxySDTrueRedshiftLSSTSRDPrivate * const self = nc_galaxy_sd_true_redshift_lsst_srd_get_instance_private (gsdtrlsst);

  self->z_min  = 0.0;
  self->z_max  = 0.0;
  self->z_norm = 0.0;
  self->y0     = 0.0;
}

static void
_nc_galaxy_sd_true_redshift_lsst_srd_dispose (GObject *object)
{
  /* Chain up: end */
  G_OBJECT_CLASS (nc_galaxy_sd_true_redshift_lsst_srd_parent_class)->dispose (object);
}

static void
_nc_galaxy_sd_true_redshift_lsst_srd_finalize (GObject *object)
{
  /* Chain up: end */
  G_OBJECT_CLASS (nc_galaxy_sd_true_redshift_lsst_srd_parent_class)->finalize (object);
}

static gdouble _nc_galaxy_sd_true_redshift_lsst_srd_gen (NcGalaxySDTrueRedshift *gsdtrlsst, NcmRNG *rng);
static gdouble _nc_galaxy_sd_true_redshift_lsst_srd_integ (NcGalaxySDTrueRedshift *gsdtrlsst, gdouble z);
static gboolean _nc_galaxy_sd_true_redshift_lsst_srd_set_lim (NcGalaxySDTrueRedshift *gsdtrlsst, const gdouble z_min, const gdouble z_max);
static gboolean _nc_galaxy_sd_true_redshift_lsst_srd_get_lim (NcGalaxySDTrueRedshift *gsdtrlsst, gdouble *z_min, gdouble *z_max);

static void
nc_galaxy_sd_true_redshift_lsst_srd_class_init (NcGalaxySDTrueRedshiftLSSTSRDClass *klass)
{
  NcGalaxySDTrueRedshiftClass *sd_redshift_class = NC_GALAXY_SD_TRUE_REDSHIFT_CLASS (klass);
  GObjectClass *object_class                     = G_OBJECT_CLASS (klass);
  NcmModelClass *model_class                     = NCM_MODEL_CLASS (klass);

  object_class->dispose  = &_nc_galaxy_sd_true_redshift_lsst_srd_dispose;
  object_class->finalize = &_nc_galaxy_sd_true_redshift_lsst_srd_finalize;

  ncm_model_class_set_name_nick (model_class, "LSST SRD Galaxy Distribution", "LSST_SRD");
  ncm_model_class_add_params (model_class, NC_GALAXY_SD_TRUE_REDSHIFT_LSST_SRD_SPARAM_LEN, 0, PROP_LEN);

  /**
   * NcGalaxySDTrueRedshiftLSSTSRD:alpha:
   *
   * The redshift exponential slope.
   *
   */
  ncm_model_class_set_sparam (model_class, NC_GALAXY_SD_TRUE_REDSHIFT_LSST_SRD_ALPHA, "\\alpha", "alpha",
                              1.0e-8, 1.0, 1.0e-2,
                              NC_GALAXY_SD_TRUE_REDSHIFT_LSST_SRD_DEFAULT_PARAMS_ABSTOL, NC_GALAXY_SD_TRUE_REDSHIFT_LSST_SRD_DEFAULT_ALPHA,
                              NCM_PARAM_TYPE_FIXED);

  /**
   * NcGalaxySDTrueRedshiftLSSTSRD:beta:
   *
   * The redshift power law slope.
   *
   */
  ncm_model_class_set_sparam (model_class, NC_GALAXY_SD_TRUE_REDSHIFT_LSST_SRD_BETA, "\\beta", "beta",
                              1.0e-8, 5.0, 1.0e-1,
                              NC_GALAXY_SD_TRUE_REDSHIFT_LSST_SRD_DEFAULT_PARAMS_ABSTOL, NC_GALAXY_SD_TRUE_REDSHIFT_LSST_SRD_DEFAULT_BETA,
                              NCM_PARAM_TYPE_FIXED);

  /**
   * NcGalaxySDTrueRedshiftLSSTSRD:z0:
   *
   * The redshift pivot.
   *
   */
  ncm_model_class_set_sparam (model_class, NC_GALAXY_SD_TRUE_REDSHIFT_LSST_SRD_Z0, "z_0", "z0",
                              1.0e-8, 1.0, 1.0e-2,
                              NC_GALAXY_SD_TRUE_REDSHIFT_LSST_SRD_DEFAULT_PARAMS_ABSTOL, NC_GALAXY_SD_TRUE_REDSHIFT_LSST_SRD_DEFAULT_Z0,
                              NCM_PARAM_TYPE_FIXED);

  ncm_model_class_check_params_info (model_class);

  sd_redshift_class->gen     = &_nc_galaxy_sd_true_redshift_lsst_srd_gen;
  sd_redshift_class->integ   = &_nc_galaxy_sd_true_redshift_lsst_srd_integ;
  sd_redshift_class->set_lim = &_nc_galaxy_sd_true_redshift_lsst_srd_set_lim;
  sd_redshift_class->get_lim = &_nc_galaxy_sd_true_redshift_lsst_srd_get_lim;
}

#define VECTOR (NCM_MODEL (gsdtr))
#define ALPHA  (ncm_model_orig_param_get (VECTOR, NC_GALAXY_SD_TRUE_REDSHIFT_LSST_SRD_ALPHA))
#define BETA   (ncm_model_orig_param_get (VECTOR, NC_GALAXY_SD_TRUE_REDSHIFT_LSST_SRD_BETA))
#define Z0     (ncm_model_orig_param_get (VECTOR, NC_GALAXY_SD_TRUE_REDSHIFT_LSST_SRD_Z0))

static gdouble
_nc_galaxy_sd_true_redshift_lsst_srd_gen (NcGalaxySDTrueRedshift *gsdtr, NcmRNG *rng)
{
  NcGalaxySDTrueRedshiftLSSTSRD *gsdtrlsst          = NC_GALAXY_SD_TRUE_REDSHIFT_LSST_SRD (gsdtr);
  NcGalaxySDTrueRedshiftLSSTSRDPrivate * const self = nc_galaxy_sd_true_redshift_lsst_srd_get_instance_private (gsdtrlsst);
  const gdouble alpha                               = ALPHA;
  const gdouble beta                                = BETA;
  const gdouble z0                                  = Z0;
  const gdouble y0                                  = pow (z0, alpha);
  const gdouble gamma_a                             = (1.0 + beta) / alpha;
  gdouble z;

  do {
    const gdouble gen_y = ncm_rng_gamma_gen (rng, gamma_a, y0);

    z = pow (gen_y, 1.0 / alpha);
  } while (z < self->z_min || z > self->z_max);

  return z;
}

static gdouble
_nc_galaxy_sd_true_redshift_lsst_srd_integ (NcGalaxySDTrueRedshift *gsdtr, gdouble z)
{
  NcGalaxySDTrueRedshiftLSSTSRD *gsdtrlsst          = NC_GALAXY_SD_TRUE_REDSHIFT_LSST_SRD (gsdtr);
  NcGalaxySDTrueRedshiftLSSTSRDPrivate * const self = nc_galaxy_sd_true_redshift_lsst_srd_get_instance_private (gsdtrlsst);
  const gdouble alpha                               = ALPHA;
  const gdouble beta                                = BETA;
  const gdouble z0                                  = Z0;
  const gdouble y                                   = pow (z0, alpha);

  if ((z >= self->z_min) && (z < self->z_max))
    return pow (z, beta) * exp (-(y / self->y0)) * self->z_norm;

  return 0.0;
}

static gboolean
_nc_galaxy_sd_true_redshift_lsst_srd_set_lim (NcGalaxySDTrueRedshift *gsdtr, const gdouble z_min, const gdouble z_max)
{
  NcGalaxySDTrueRedshiftLSSTSRD *gsdtrlsst          = NC_GALAXY_SD_TRUE_REDSHIFT_LSST_SRD (gsdtr);
  NcGalaxySDTrueRedshiftLSSTSRDPrivate * const self = nc_galaxy_sd_true_redshift_lsst_srd_get_instance_private (gsdtrlsst);
  const gdouble alpha                               = ALPHA;
  const gdouble beta                                = BETA;

  g_assert_cmpfloat (z_min, <, z_max);

  self->z_min  = z_min;
  self->z_max  = z_max;
  self->y0     = pow (Z0, alpha);
  self->z_norm = alpha / pow (Z0, 1 + beta) / (gsl_sf_gamma_inc ((1 + beta) / alpha, pow (z_min / Z0, alpha)) - gsl_sf_gamma_inc ((1 + beta) / alpha, pow (z_max / Z0, alpha)));

  return TRUE;
}

static gboolean
_nc_galaxy_sd_true_redshift_lsst_srd_get_lim (NcGalaxySDTrueRedshift *gsdtr, gdouble *z_min, gdouble *z_max)
{
  NcGalaxySDTrueRedshiftLSSTSRD *gsdtrlsst          = NC_GALAXY_SD_TRUE_REDSHIFT_LSST_SRD (gsdtr);
  NcGalaxySDTrueRedshiftLSSTSRDPrivate * const self = nc_galaxy_sd_true_redshift_lsst_srd_get_instance_private (gsdtrlsst);

  g_assert_nonnull (z_min);
  g_assert_nonnull (z_max);

  *z_min = self->z_min;
  *z_max = self->z_max;

  return TRUE;
}

/**
 * nc_galaxy_sd_true_redshift_lsst_srd_new:
 * @z_min: the minimum redshift
 * @z_max: the maximum redshift
 *
 * Creates a new #NcGalaxySDPositionLSSTSRD, the parameter values correspond to the
 * LSST SRD year 1.
 *
 * Returns: (transfer full): a new #NcGalaxySDPositionLSSTSRD
 */
NcGalaxySDTrueRedshiftLSSTSRD *
nc_galaxy_sd_true_redshift_lsst_srd_new (const gdouble z_min, const gdouble z_max)
{
  NcmDTuple2 lim                           = NCM_DTUPLE2_STATIC_INIT (z_min, z_max);
  NcGalaxySDTrueRedshiftLSSTSRD *gsdtrlsst = g_object_new (NC_TYPE_GALAXY_SD_TRUE_REDSHIFT_LSST_SRD,
                                                           "lim", &lim,
                                                           NULL);

  return gsdtrlsst;
}

/**
 * nc_galaxy_sd_true_redshift_lsst_srd_new_y10:
 * @z_min: the minimum redshift
 * @z_max: the maximum redshift
 *
 * Creates a new #NcGalaxySDPositionLSSTSRD, the parameter values correspond to the
 * LSST SRD year 10.
 *
 * Returns: (transfer full): a new #NcGalaxySDPositionLSSTSRD
 */
NcGalaxySDTrueRedshiftLSSTSRD *
nc_galaxy_sd_true_redshift_lsst_srd_new_y10 (const gdouble z_min, const gdouble z_max)
{
  NcmDTuple2 lim                           = NCM_DTUPLE2_STATIC_INIT (z_min, z_max);
  NcGalaxySDTrueRedshiftLSSTSRD *gsdtrlsst = g_object_new (NC_TYPE_GALAXY_SD_TRUE_REDSHIFT_LSST_SRD,
                                                           "lim", &lim,
                                                           "alpha", NC_GALAXY_SD_TRUE_REDSHIFT_LSST_SRD_Y10_ALPHA,
                                                           "beta", NC_GALAXY_SD_TRUE_REDSHIFT_LSST_SRD_Y10_BETA,
                                                           "z0", NC_GALAXY_SD_TRUE_REDSHIFT_LSST_SRD_Y10_Z0,
                                                           NULL);

  return gsdtrlsst;
}

/**
 * nc_galaxy_sd_true_redshift_lsst_srd_ref:
 * @gsdtrlsst: a #NcGalaxySDTrueRedshiftLSSTSRD
 *
 * Increases the reference count of @gsdtrlsst by one.
 *
 * Returns: (transfer full): @gsdtrlsst.
 */
NcGalaxySDTrueRedshiftLSSTSRD *
nc_galaxy_sd_true_redshift_lsst_srd_ref (NcGalaxySDTrueRedshiftLSSTSRD *gsdtrlsst)
{
  return g_object_ref (gsdtrlsst);
}

/**
 * nc_galaxy_sd_true_redshift_lsst_srd_free:
 * @gsdtrlsst: a #NcGalaxySDTrueRedshiftLSSTSRD
 *
 * Decreases the reference count of @gsdtrlsst by one.
 *
 */
void
nc_galaxy_sd_true_redshift_lsst_srd_free (NcGalaxySDTrueRedshiftLSSTSRD *gsdtrlsst)
{
  g_object_unref (gsdtrlsst);
}

/**
 * nc_galaxy_sd_true_redshift_lsst_srd_clear:
 * @gsdtrlsst: a #NcGalaxySDTrueRedshiftLSSTSRD
 *
 * Decreases the reference count of @gsdtrlsst by one, and sets the
 * pointer @gsdtrlsst to NULL.
 *
 */
void
nc_galaxy_sd_true_redshift_lsst_srd_clear (NcGalaxySDTrueRedshiftLSSTSRD **gsdtrlsst)
{
  g_clear_object (gsdtrlsst);
}

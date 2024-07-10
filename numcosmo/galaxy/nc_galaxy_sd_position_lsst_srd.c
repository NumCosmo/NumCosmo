/***************************************************************************
 *            nc_galaxy_sd_position_lsst_srd.c
 *
 *  Tue June 22 14:59:27 2023
 *  Copyright  2023  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * nc_galaxy_sd_position_lsst_srd.c
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
 * SECTION:nc_galaxy_sd_position_lsst_srd
 * @title: NcGalaxySDPositionLSSTSRD
 * @short_description: Class describing galaxy sample position distributions as in LSST-SRD
 * @stability: Unstable
 *
 * Class defining a galaxy sample position distribution as described in the LSST
 * Science Roadmap Document, it includes the redshift distribution and a angular radius
 * distribution.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

/* #include "nc_enum_types.h" */
#include "galaxy/nc_galaxy_sd_position_lsst_srd.h"
#include "galaxy/nc_galaxy_sd_position.h"
#include "math/ncm_stats_dist1d.h"
#include "math/ncm_stats_dist1d_spline.h"
#include "math/ncm_dtuple.h"
#include "math/ncm_spline.h"
#include "math/ncm_spline_func.h"
#include "math/ncm_spline_cubic.h"
#include "math/ncm_spline_cubic_notaknot.h"
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf.h>

typedef struct _NcGalaxySDPositionLSSTSRDPrivate
{
  gdouble z_lb;
  gdouble z_ub;
  gdouble theta_norm;
  gdouble z_norm;
  gdouble theta_lb;
  gdouble theta_ub;
  gdouble theta_lb2;
  gdouble theta_ub2;
  gdouble y0;
} NcGalaxySDPositionLSSTSRDPrivate;

struct _NcGalaxySDPositionLSSTSRD
{
  NcGalaxySDPosition parent_instance;
};

enum
{
  PROP_0,
  PROP_LEN,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcGalaxySDPositionLSSTSRD, nc_galaxy_sd_position_lsst_srd, NC_TYPE_GALAXY_SD_POSITION);

static void
nc_galaxy_sd_position_lsst_srd_init (NcGalaxySDPositionLSSTSRD *gsdplsst)
{
  NcGalaxySDPositionLSSTSRDPrivate * const self = nc_galaxy_sd_position_lsst_srd_get_instance_private (gsdplsst);

  self->z_lb       = 0.0;
  self->z_ub       = 0.0;
  self->theta_norm = 0.0;
  self->z_norm     = 0.0;
  self->theta_lb   = 0.0;
  self->theta_ub   = 0.0;
  self->theta_lb2  = 0.0;
  self->theta_ub2  = 0.0;
  self->y0         = 0.0;
}

static void
_nc_galaxy_sd_position_lsst_srd_dispose (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_galaxy_sd_position_lsst_srd_parent_class)->dispose (object);
}

static void
_nc_galaxy_sd_position_lsst_srd_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_galaxy_sd_position_lsst_srd_parent_class)->finalize (object);
}

static gdouble _nc_galaxy_sd_position_lsst_srd_gen_theta (NcGalaxySDPosition *gsdp, NcmRNG *rng);
static gdouble _nc_galaxy_sd_position_lsst_srd_gen_z (NcGalaxySDPosition *gsdp, NcmRNG *rng);
static gdouble _nc_galaxy_sd_position_lsst_srd_integ (NcGalaxySDPosition *gsdp, const gdouble theta, const gdouble z);
static void _nc_galaxy_sd_position_lsst_srd_set_z_lim (NcGalaxySDPosition *gsdp, gdouble z_min, gdouble z_max);
static void _nc_galaxy_sd_position_lsst_srd_get_z_lim (NcGalaxySDPosition *gsdp, gdouble *z_min, gdouble *z_max);
static void _nc_galaxy_sd_position_lsst_srd_set_theta_lim (NcGalaxySDPosition *gsdp, gdouble theta_min, gdouble theta_max);
static void _nc_galaxy_sd_position_lsst_srd_get_theta_lim (NcGalaxySDPosition *gsdp, gdouble *theta_min, gdouble *theta_max);

static void
nc_galaxy_sd_position_lsst_srd_class_init (NcGalaxySDPositionLSSTSRDClass *klass)
{
  NcGalaxySDPositionClass *sd_position_class = NC_GALAXY_SD_POSITION_CLASS (klass);
  GObjectClass *object_class                 = G_OBJECT_CLASS (klass);
  NcmModelClass *model_class                 = NCM_MODEL_CLASS (klass);

  object_class->dispose  = &_nc_galaxy_sd_position_lsst_srd_dispose;
  object_class->finalize = &_nc_galaxy_sd_position_lsst_srd_finalize;

  ncm_model_class_set_name_nick (model_class, "LSST SRD Galaxy Distribution", "LSST_SRD");
  ncm_model_class_add_params (model_class, NC_GALAXY_SD_POSITION_LSST_SRD_SPARAM_LEN, 0, PROP_LEN);

  /**
   * NcGalaxySDPositionLSSTSRD:alpha:
   *
   * The redshift exponential slope.
   *
   */
  ncm_model_class_set_sparam (model_class, NC_GALAXY_SD_POSITION_LSST_SRD_ALPHA, "\\alpha", "alpha",
                              1e-8,  1.0, 1.0e-2,
                              NC_GALAXY_SD_POSITION_LSST_SRD_DEFAULT_PARAMS_ABSTOL, NC_GALAXY_SD_POSITION_LSST_SRD_DEFAULT_ALPHA,
                              NCM_PARAM_TYPE_FIXED);

  /**
   * NcGalaxySDPositionLSSTSRD:beta:
   *
   * The redshift power law slope.
   *
   */
  ncm_model_class_set_sparam (model_class, NC_GALAXY_SD_POSITION_LSST_SRD_BETA, "\\beta", "beta",
                              1e-8,  5.0, 1.0e-1,
                              NC_GALAXY_SD_POSITION_LSST_SRD_DEFAULT_PARAMS_ABSTOL, NC_GALAXY_SD_POSITION_LSST_SRD_DEFAULT_BETA,
                              NCM_PARAM_TYPE_FIXED);

  /**
   * NcGalaxySDPositionLSSTSRD:z0:
   *
   * The redshift pivot.
   *
   */
  ncm_model_class_set_sparam (model_class, NC_GALAXY_SD_POSITION_LSST_SRD_Z0, "z_0", "z0",
                              1e-8,  10.0, 1.0e-2,
                              NC_GALAXY_SD_POSITION_LSST_SRD_DEFAULT_PARAMS_ABSTOL, NC_GALAXY_SD_POSITION_LSST_SRD_DEFAULT_Z0,
                              NCM_PARAM_TYPE_FIXED);

  ncm_model_class_check_params_info (model_class);

  sd_position_class->gen_theta     = &_nc_galaxy_sd_position_lsst_srd_gen_theta;
  sd_position_class->gen_z         = &_nc_galaxy_sd_position_lsst_srd_gen_z;
  sd_position_class->integ         = &_nc_galaxy_sd_position_lsst_srd_integ;
  sd_position_class->set_z_lim     = &_nc_galaxy_sd_position_lsst_srd_set_z_lim;
  sd_position_class->get_z_lim     = &_nc_galaxy_sd_position_lsst_srd_get_z_lim;
  sd_position_class->set_theta_lim = &_nc_galaxy_sd_position_lsst_srd_set_theta_lim;
  sd_position_class->get_theta_lim = &_nc_galaxy_sd_position_lsst_srd_get_theta_lim;
}

#define VECTOR (NCM_MODEL (gsdp))
#define ALPHA  (ncm_model_orig_param_get (VECTOR, NC_GALAXY_SD_POSITION_LSST_SRD_ALPHA))
#define BETA   (ncm_model_orig_param_get (VECTOR, NC_GALAXY_SD_POSITION_LSST_SRD_BETA))
#define Z0     (ncm_model_orig_param_get (VECTOR, NC_GALAXY_SD_POSITION_LSST_SRD_Z0))

static gdouble
_nc_galaxy_sd_position_lsst_srd_gen_theta (NcGalaxySDPosition *gsdp, NcmRNG *rng)
{
  NcGalaxySDPositionLSSTSRD *gsdplsst           = NC_GALAXY_SD_POSITION_LSST_SRD (gsdp);
  NcGalaxySDPositionLSSTSRDPrivate * const self = nc_galaxy_sd_position_lsst_srd_get_instance_private (gsdplsst);
  const gdouble cumul_gen                       = ncm_rng_uniform_gen (rng, 0.0, 1.0);

  return sqrt (cumul_gen * 2.0 / self->theta_norm + self->theta_lb2);
}

static gdouble
_nc_galaxy_sd_position_lsst_srd_gen_z (NcGalaxySDPosition *gsdp, NcmRNG *rng)
{
  NcGalaxySDPositionLSSTSRD *gsdplsst           = NC_GALAXY_SD_POSITION_LSST_SRD (gsdp);
  NcGalaxySDPositionLSSTSRDPrivate * const self = nc_galaxy_sd_position_lsst_srd_get_instance_private (gsdplsst);
  const gdouble alpha                           = ALPHA;
  const gdouble beta                            = BETA;
  const gdouble z0                              = Z0;
  const gdouble y0                              = pow (z0, alpha);
  const gdouble gamma_a                         = (1.0 + beta) / alpha;

  gdouble z;

  do {
    const gdouble gen_y = ncm_rng_gamma_gen (rng, gamma_a, y0);

    z = pow (gen_y, 1.0 / alpha);
  } while (z < self->z_lb || z > self->z_ub);

  return z;
}

static gdouble
_nc_galaxy_sd_position_lsst_srd_integ (NcGalaxySDPosition *gsdp, const gdouble theta, const gdouble z)
{
  NcGalaxySDPositionLSSTSRD *gsdplsst           = NC_GALAXY_SD_POSITION_LSST_SRD (gsdp);
  NcGalaxySDPositionLSSTSRDPrivate * const self = nc_galaxy_sd_position_lsst_srd_get_instance_private (gsdplsst);
  const gdouble alpha                           = ALPHA;
  const gdouble beta                            = BETA;
  const gdouble z0                              = Z0;
  const gdouble y                               = pow (z, alpha);

  return pow (z, beta) * exp (-(y / self->y0)) * theta * self->theta_norm * self->z_norm;
  /* return pow (z, beta) * exp (-(y / self->y0)) * theta * self->theta_norm; */
}

static void
_nc_galaxy_sd_position_lsst_srd_set_z_lim (NcGalaxySDPosition *gsdp, gdouble z_min, gdouble z_max)
{
  NcGalaxySDPositionLSSTSRD *gsdplsst           = NC_GALAXY_SD_POSITION_LSST_SRD (gsdp);
  NcGalaxySDPositionLSSTSRDPrivate * const self = nc_galaxy_sd_position_lsst_srd_get_instance_private (gsdplsst);
  const gdouble alpha                           = ALPHA;
  const gdouble beta                            = BETA;

  g_assert_cmpfloat (z_min, <, z_max);

  self->z_lb   = z_min;
  self->z_ub   = z_max;
  self->y0     = pow (Z0, alpha);
  self->z_norm =  alpha / pow (Z0, 1 + beta) / (gsl_sf_gamma_inc ((1 + beta) / alpha, pow (z_min / Z0, alpha)) - gsl_sf_gamma_inc ((1 + beta) / alpha, pow (z_max / Z0, alpha)));
}

static void
_nc_galaxy_sd_position_lsst_srd_get_z_lim (NcGalaxySDPosition *gsdp, gdouble *z_min, gdouble *z_max)
{
  NcGalaxySDPositionLSSTSRD *gsdplsst           = NC_GALAXY_SD_POSITION_LSST_SRD (gsdp);
  NcGalaxySDPositionLSSTSRDPrivate * const self = nc_galaxy_sd_position_lsst_srd_get_instance_private (gsdplsst);

  g_assert_nonnull (z_min);
  g_assert_nonnull (z_max);

  *z_min = self->z_lb;
  *z_max = self->z_ub;
}

static void
_nc_galaxy_sd_position_lsst_srd_set_theta_lim (NcGalaxySDPosition *gsdp, gdouble theta_min, gdouble theta_max)
{
  NcGalaxySDPositionLSSTSRD *gsdplsst           = NC_GALAXY_SD_POSITION_LSST_SRD (gsdp);
  NcGalaxySDPositionLSSTSRDPrivate * const self = nc_galaxy_sd_position_lsst_srd_get_instance_private (gsdplsst);

  g_assert_cmpfloat (theta_min, <, theta_max);

  self->theta_lb = theta_min;
  self->theta_ub = theta_max;

  self->theta_lb2  = self->theta_lb * self->theta_lb;
  self->theta_ub2  = self->theta_ub * self->theta_ub;
  self->theta_norm = 2.0 / (self->theta_ub2 - self->theta_lb2);
}

static void
_nc_galaxy_sd_position_lsst_srd_get_theta_lim (NcGalaxySDPosition *gsdp, gdouble *theta_min, gdouble *theta_max)
{
  NcGalaxySDPositionLSSTSRD *gsdplsst           = NC_GALAXY_SD_POSITION_LSST_SRD (gsdp);
  NcGalaxySDPositionLSSTSRDPrivate * const self = nc_galaxy_sd_position_lsst_srd_get_instance_private (gsdplsst);

  g_assert_nonnull (theta_min);
  g_assert_nonnull (theta_max);

  *theta_min = self->theta_lb;
  *theta_max = self->theta_ub;
}

/**
 * nc_galaxy_sd_position_lsst_srd_new:
 * @z_min: the minimum redshift
 * @z_max: the maximum redshift
 * @theta_min: the minimum angular radius
 * @theta_max: the maximum angular radius
 *
 * Creates a new #NcGalaxySDPositionLSSTSRD, the parameter values correspond to the
 * LSST SRD year 1.
 *
 * Returns: (transfer full): a new NcGalaxySDPositionLSSTSRD.
 */
NcGalaxySDPositionLSSTSRD *
nc_galaxy_sd_position_lsst_srd_new (const gdouble z_min, const gdouble z_max, const gdouble theta_min, const gdouble theta_max)
{
  NcmDTuple2 z_lim                    = NCM_DTUPLE2_STATIC_INIT (z_min, z_max);
  NcmDTuple2 theta_lim                = NCM_DTUPLE2_STATIC_INIT (theta_min, theta_max);
  NcGalaxySDPositionLSSTSRD *gsdplsst = g_object_new (NC_TYPE_GALAXY_SD_POSITION_LSST_SRD,
                                                      "z-lim", &z_lim,
                                                      "theta-lim", &theta_lim,
                                                      NULL);

  return gsdplsst;
}

/**
 * nc_galaxy_sd_position_lsst_srd_new_y10:
 * @z_min: the minimum redshift
 * @z_max: the maximum redshift
 * @theta_min: the minimum angular radius
 * @theta_max: the maximum angular radius
 *
 * Creates a new #NcGalaxySDPositionLSSTSRD, the parameter values correspond to the
 * LSST SRD year 10.
 *
 * Returns: (transfer full): a new NcGalaxySDPositionLSSTSRD.
 */
NcGalaxySDPositionLSSTSRD *
nc_galaxy_sd_position_lsst_srd_new_y10 (const gdouble z_min, const gdouble z_max, const gdouble theta_min, const gdouble theta_max)
{
  NcmDTuple2 z_lim                    = NCM_DTUPLE2_STATIC_INIT (z_min, z_max);
  NcmDTuple2 theta_lim                = NCM_DTUPLE2_STATIC_INIT (theta_min, theta_max);
  NcGalaxySDPositionLSSTSRD *gsdplsst = g_object_new (NC_TYPE_GALAXY_SD_POSITION_LSST_SRD,
                                                      "alpha", NC_GALAXY_SD_POSITION_LSST_SRD_Y10_ALPHA,
                                                      "beta", NC_GALAXY_SD_POSITION_LSST_SRD_Y10_BETA,
                                                      "z0", NC_GALAXY_SD_POSITION_LSST_SRD_Y10_Z0,
                                                      "z-lim", &z_lim,
                                                      "theta-lim", &theta_lim,
                                                      NULL);

  return gsdplsst;
}

/**
 * nc_galaxy_sd_position_lsst_srd_ref:
 * @gsdplsst: a #NcGalaxySDPositionLSSTSRD
 *
 * Increase the reference of @gsdplsst by one.
 *
 * Returns: (transfer full): @gsdplsst.
 */
NcGalaxySDPositionLSSTSRD *
nc_galaxy_sd_position_lsst_srd_ref (NcGalaxySDPositionLSSTSRD *gsdplsst)
{
  return g_object_ref (gsdplsst);
}

/**
 * nc_galaxy_sd_position_lsst_srd_free:
 * @gsdplsst: a #NcGalaxySDPositionLSSTSRD
 *
 * Decrease the reference count of @gsdplsst by one.
 *
 */
void
nc_galaxy_sd_position_lsst_srd_free (NcGalaxySDPositionLSSTSRD *gsdplsst)
{
  g_object_unref (gsdplsst);
}

/**
 * nc_galaxy_sd_position_lsst_srd_clear:
 * @gsdplsst: a #NcGalaxySDPositionLSSTSRD
 *
 * Decrease the reference count of @gsdplsst by one, and sets the pointer *@gsdplsst to
 * NULL.
 *
 */
void
nc_galaxy_sd_position_lsst_srd_clear (NcGalaxySDPositionLSSTSRD **gsdplsst)
{
  g_clear_object (gsdplsst);
}


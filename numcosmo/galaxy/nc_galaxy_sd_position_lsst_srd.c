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
#include "math/ncm_spline.h"
#include "math/ncm_spline_func.h"
#include "math/ncm_spline_cubic.h"
#include "math/ncm_spline_cubic_notaknot.h"
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>

typedef struct _NcGalaxySDPositionLSSTSRDPrivate
{
  NcmVector *z_lim;
  NcmVector *r_lim;
  gdouble z_lb;
  gdouble z_ub;
  gdouble r_norm;
  gdouble r_lb;
  gdouble r_ub;
  gdouble r_lb2;
  gdouble r_ub2;
} NcGalaxySDPositionLSSTSRDPrivate;

struct _NcGalaxySDPositionLSSTSRD
{
  NcGalaxySDPosition parent_instance;
};

enum
{
  PROP_0,
  PROP_Z_LIM,
  PROP_R_LIM,
  PROP_Z_DIST,
  PROP_LEN,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcGalaxySDPositionLSSTSRD, nc_galaxy_sd_position_lsst_srd, NC_TYPE_GALAXY_SD_POSITION);

static void
nc_galaxy_sd_position_lsst_srd_init (NcGalaxySDPositionLSSTSRD *gsdplsst)
{
  NcGalaxySDPositionLSSTSRDPrivate * const self = nc_galaxy_sd_position_lsst_srd_get_instance_private (gsdplsst);

  self->z_lim  = NULL;
  self->r_lim  = NULL;
  self->z_lb   = 0.0;
  self->z_ub   = 0.0;
  self->r_norm = 0.0;
  self->r_lb   = 0.0;
  self->r_ub   = 0.0;
  self->r_lb2  = 0.0;
  self->r_ub2  = 0.0;
}

static void
_nc_galaxy_sd_position_lsst_srd_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcGalaxySDPositionLSSTSRD *gsdplsst = NC_GALAXY_SD_POSITION_LSST_SRD (object);

  g_return_if_fail (NC_IS_GALAXY_SD_POSITION_LSST_SRD (object));

  switch (prop_id)
  {
    case PROP_Z_LIM:
      nc_galaxy_sd_position_lsst_srd_set_z_lim (gsdplsst, g_value_get_object (value));
      break;
    case PROP_R_LIM:
      nc_galaxy_sd_position_lsst_srd_set_r_lim (gsdplsst, g_value_get_object (value));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_galaxy_sd_position_lsst_srd_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcGalaxySDPositionLSSTSRD *gsdplsst = NC_GALAXY_SD_POSITION_LSST_SRD (object);

  g_return_if_fail (NC_IS_GALAXY_SD_POSITION_LSST_SRD (object));

  switch (prop_id)
  {
    case PROP_Z_LIM:
      g_value_set_object (value, nc_galaxy_sd_position_lsst_srd_peek_z_lim (gsdplsst));
      break;
    case PROP_R_LIM:
      g_value_set_object (value, nc_galaxy_sd_position_lsst_srd_peek_r_lim (gsdplsst));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_galaxy_sd_position_lsst_srd_dispose (GObject *object)
{
  NcGalaxySDPositionLSSTSRD *gsdplsst           = NC_GALAXY_SD_POSITION_LSST_SRD (object);
  NcGalaxySDPositionLSSTSRDPrivate * const self = nc_galaxy_sd_position_lsst_srd_get_instance_private (gsdplsst);

  ncm_vector_clear (&self->z_lim);
  ncm_vector_clear (&self->r_lim);

  G_OBJECT_CLASS (nc_galaxy_sd_position_lsst_srd_parent_class)->dispose (object);
}

static void
_nc_galaxy_sd_position_lsst_srd_finalize (GObject *object)
{
  G_OBJECT_CLASS (nc_galaxy_sd_position_lsst_srd_parent_class)->finalize (object);
}

static gdouble _nc_galaxy_sd_position_lsst_srd_gen_r (NcGalaxySDPosition *gsdp, NcmRNG *rng);
static gdouble _nc_galaxy_sd_position_lsst_srd_gen_z (NcGalaxySDPosition *gsdp, NcmRNG *rng);
static gdouble _nc_galaxy_sd_position_lsst_srd_integ (NcGalaxySDPosition *gsdp, const gdouble r, const gdouble z);

static void
nc_galaxy_sd_position_lsst_srd_class_init (NcGalaxySDPositionLSSTSRDClass *klass)
{
  NcGalaxySDPositionClass *sd_position_class = NC_GALAXY_SD_POSITION_CLASS (klass);
  GObjectClass *object_class                 = G_OBJECT_CLASS (klass);
  NcmModelClass *model_class                 = NCM_MODEL_CLASS (klass);

  model_class->set_property = &_nc_galaxy_sd_position_lsst_srd_set_property;
  model_class->get_property = &_nc_galaxy_sd_position_lsst_srd_get_property;
  object_class->dispose     = &_nc_galaxy_sd_position_lsst_srd_dispose;
  object_class->finalize    = &_nc_galaxy_sd_position_lsst_srd_finalize;

  ncm_model_class_set_name_nick (model_class, "LSST SRD Galaxy Distribution", "LSST_SRD");
  ncm_model_class_add_params (model_class, NC_GALAXY_SD_POSITION_LSST_SRD_SPARAM_LEN, 0, PROP_LEN);

  /**
   * NcGalaxySDPositionLSSTSRD:z-lim:
   *
   * Galaxy sample redshift distribution limits.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_Z_LIM,
                                   g_param_spec_object ("z-lim",
                                                        NULL,
                                                        "Galaxy sample redshift distribution limits",
                                                        NCM_TYPE_VECTOR,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcGalaxySDPositionLSSTSRD:r-lim:
   *
   * Galaxy sample radius distribution limits.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_R_LIM,
                                   g_param_spec_object ("r-lim",
                                                        NULL,
                                                        "Galaxy sample radius distribution limits",
                                                        NCM_TYPE_VECTOR,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

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

  sd_position_class->gen_r = &_nc_galaxy_sd_position_lsst_srd_gen_r;
  sd_position_class->gen_z = &_nc_galaxy_sd_position_lsst_srd_gen_z;
  sd_position_class->integ = &_nc_galaxy_sd_position_lsst_srd_integ;
}

#define VECTOR (NCM_MODEL (gsdp))
#define ALPHA  (ncm_model_orig_param_get (VECTOR, NC_GALAXY_SD_POSITION_LSST_SRD_ALPHA))
#define BETA   (ncm_model_orig_param_get (VECTOR, NC_GALAXY_SD_POSITION_LSST_SRD_BETA))
#define Z0     (ncm_model_orig_param_get (VECTOR, NC_GALAXY_SD_POSITION_LSST_SRD_Z0))

static gdouble
_nc_galaxy_sd_position_lsst_srd_gen_r (NcGalaxySDPosition *gsdp, NcmRNG *rng)
{
  NcGalaxySDPositionLSSTSRD *gsdplsst           = NC_GALAXY_SD_POSITION_LSST_SRD (gsdp);
  NcGalaxySDPositionLSSTSRDPrivate * const self = nc_galaxy_sd_position_lsst_srd_get_instance_private (gsdplsst);
  const gdouble cumul_gen                       = ncm_rng_uniform_gen (rng, 0.0, 1.0);

  return sqrt (cumul_gen * 2.0 / self->r_norm + self->r_lb2);
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
_nc_galaxy_sd_position_lsst_srd_integ (NcGalaxySDPosition *gsdp, const gdouble r, const gdouble z)
{
  NcGalaxySDPositionLSSTSRD *gsdplsst           = NC_GALAXY_SD_POSITION_LSST_SRD (gsdp);
  NcGalaxySDPositionLSSTSRDPrivate * const self = nc_galaxy_sd_position_lsst_srd_get_instance_private (gsdplsst);
  const gdouble alpha                           = ALPHA;
  const gdouble beta                            = BETA;
  const gdouble z0                              = Z0;
  const gdouble y0                              = pow (z0, alpha);
  const gdouble gamma_a                         = (1.0 + beta) / alpha;
  const gdouble y                               = pow (z, alpha);

  return gsl_ran_gamma_pdf (y, gamma_a, y0) * alpha * y / z * r * self->r_norm;
}

/**
 * nc_galaxy_sd_position_lsst_srd_new:
 *
 * Creates a new #NcGalaxySDPositionLSSTSRD, the parameter values correspond to the
 * LSST SRD year 1.
 *
 * Returns: (transfer full): a new NcGalaxySDPositionLSSTSRD.
 */
NcGalaxySDPositionLSSTSRD *
nc_galaxy_sd_position_lsst_srd_new ()
{
  NcGalaxySDPositionLSSTSRD *gsdplsst = g_object_new (NC_TYPE_GALAXY_SD_POSITION_LSST_SRD,
                                                      NULL);

  return gsdplsst;
}

/**
 * nc_galaxy_sd_position_lsst_srd_new_y10:
 *
 * Creates a new #NcGalaxySDPositionLSSTSRD, the parameter values correspond to the
 * LSST SRD year 10.
 *
 * Returns: (transfer full): a new NcGalaxySDPositionLSSTSRD.
 */
NcGalaxySDPositionLSSTSRD *
nc_galaxy_sd_position_lsst_srd_new_y10 ()
{
  NcGalaxySDPositionLSSTSRD *gsdplsst = g_object_new (NC_TYPE_GALAXY_SD_POSITION_LSST_SRD,
                                                      "alpha", NC_GALAXY_SD_POSITION_LSST_SRD_Y10_ALPHA,
                                                      "beta", NC_GALAXY_SD_POSITION_LSST_SRD_Y10_BETA,
                                                      "z0", NC_GALAXY_SD_POSITION_LSST_SRD_Y10_Z0,
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

/**
 * nc_galaxy_sd_position_lsst_srd_set_z_lim:
 * @gsdplsst: a #NcGalaxySDPositionLSSTSRD
 * @lim: a #NcmVector
 *
 * Sets the redshift limits @lim.
 */
void
nc_galaxy_sd_position_lsst_srd_set_z_lim (NcGalaxySDPositionLSSTSRD *gsdplsst, NcmVector *lim)
{
  NcGalaxySDPositionLSSTSRDPrivate * const self = nc_galaxy_sd_position_lsst_srd_get_instance_private (gsdplsst);

  g_assert_cmpuint (ncm_vector_len (lim), ==, 2);

  self->z_lb = ncm_vector_get (lim, 0);
  self->z_ub = ncm_vector_get (lim, 1);

  g_assert_cmpfloat (self->z_lb, <, self->z_ub);

  if (self->z_lim != lim)
  {
    ncm_vector_clear (&self->z_lim);
    self->z_lim = ncm_vector_ref (lim);
  }
}

/**
 * nc_galaxy_sd_position_lsst_srd_peek_z_lim:
 * @gsdplsst: a #NcGalaxySDPositionLSSTSRD
 *
 * Gets the redshift limits. The returned vector should not be modified.
 *
 * Returns: (transfer none): the redshift limits.
 */
NcmVector *
nc_galaxy_sd_position_lsst_srd_peek_z_lim (NcGalaxySDPositionLSSTSRD *gsdplsst)
{
  NcGalaxySDPositionLSSTSRDPrivate * const self = nc_galaxy_sd_position_lsst_srd_get_instance_private (gsdplsst);

  return self->z_lim;
}

/**
 * nc_galaxy_sd_position_lsst_srd_set_r_lim:
 * @gsdplsst: a #NcGalaxySDPositionLSSTSRD
 * @lim: a #NcmVector
 *
 * Sets the radius limits @lim.
 */
void
nc_galaxy_sd_position_lsst_srd_set_r_lim (NcGalaxySDPositionLSSTSRD *gsdplsst, NcmVector *lim)
{
  NcGalaxySDPositionLSSTSRDPrivate * const self = nc_galaxy_sd_position_lsst_srd_get_instance_private (gsdplsst);

  g_assert_cmpuint (ncm_vector_len (lim), ==, 2);

  self->r_lb = ncm_vector_get (lim, 0);
  self->r_ub = ncm_vector_get (lim, 1);

  g_assert_cmpfloat (self->r_lb, <, self->r_ub);

  self->r_lb2  = self->r_lb * self->r_lb;
  self->r_ub2  = self->r_ub * self->r_ub;
  self->r_norm = 2.0 / (self->r_ub2 - self->r_lb2);

  if (self->r_lim != lim)
  {
    ncm_vector_clear (&self->r_lim);
    self->r_lim = ncm_vector_ref (lim);
  }
}

/**
 * nc_galaxy_sd_position_lsst_srd_peek_r_lim:
 * @gsdplsst: a #NcGalaxySDPositionLSSTSRD
 *
 * Gets the radius limits. The returned vector should not be modified.
 *
 * Returns: (transfer none): the radius limits.
 */
NcmVector *
nc_galaxy_sd_position_lsst_srd_peek_r_lim (NcGalaxySDPositionLSSTSRD *gsdplsst)
{
  NcGalaxySDPositionLSSTSRDPrivate * const self = nc_galaxy_sd_position_lsst_srd_get_instance_private (gsdplsst);

  return self->r_lim;
}


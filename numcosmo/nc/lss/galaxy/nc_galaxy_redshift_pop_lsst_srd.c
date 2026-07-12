/***************************************************************************
 *            nc_galaxy_redshift_pop_lsst_srd.c
 *
 *  Wed Jul 31 21:45:13 2024
 *  Copyright  2024  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 *  Copyright  2024  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_galaxy_redshift_pop_lsst_srd.c
 * Copyright (C) 2024 Caio Lima de Oliveira <caiolimadeoliveira@pm.me>
 * Copyright (C) 2024 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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
 * NcGalaxyRedshiftPopLSSTSRD:
 *
 * Class describing galaxy sample redshift distributions as in LSST-SRD.
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
#include "nc/lss/galaxy/nc_galaxy_redshift_pop_lsst_srd.h"
#include "nc/lss/galaxy/nc_galaxy_redshift_pop.h"
#include "ncm/core/ncm_dtuple.h"
#include "ncm/model/ncm_model.h"
#include "ncm/core/ncm_rng.h"
#include "ncm/algebra/ncm_vector.h"

typedef struct _NcGalaxyRedshiftPopLSSTSRDPrivate
{
  gdouble z_min;
  gdouble z_max;
  gdouble z_norm;
  gdouble ln_z_norm;
  gdouble y0;
  gdouble alpha;
  gdouble beta;
  gdouble z0;
  gdouble gamma_a;
} NcGalaxyRedshiftPopLSSTSRDPrivate;

struct _NcGalaxyRedshiftPopLSSTSRD
{
  NcGalaxyRedshiftPop parent_instance;
};

enum
{
  PROP_0,
  PROP_LEN,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcGalaxyRedshiftPopLSSTSRD, nc_galaxy_redshift_pop_lsst_srd, NC_TYPE_GALAXY_REDSHIFT_POP);

static void
nc_galaxy_redshift_pop_lsst_srd_init (NcGalaxyRedshiftPopLSSTSRD *gsdrplsst)
{
  NcGalaxyRedshiftPopLSSTSRDPrivate * const self = nc_galaxy_redshift_pop_lsst_srd_get_instance_private (gsdrplsst);

  self->z_min   = 0.0;
  self->z_max   = 0.0;
  self->z_norm  = 0.0;
  self->y0      = 0.0;
  self->alpha   = 0.0;
  self->beta    = 0.0;
  self->z0      = 0.0;
  self->gamma_a = 0.0;
}

static void
_nc_galaxy_redshift_pop_lsst_srd_dispose (GObject *object)
{
  /* Chain up: end */
  G_OBJECT_CLASS (nc_galaxy_redshift_pop_lsst_srd_parent_class)->dispose (object);
}

static void
_nc_galaxy_redshift_pop_lsst_srd_finalize (GObject *object)
{
  /* Chain up: end */
  G_OBJECT_CLASS (nc_galaxy_redshift_pop_lsst_srd_parent_class)->finalize (object);
}

static gdouble _nc_galaxy_redshift_pop_lsst_srd_gen (NcGalaxyRedshiftPop *gsdrplsst, NcmRNG *rng);
static gdouble _nc_galaxy_redshift_pop_lsst_srd_eval (NcGalaxyRedshiftPop *gsdrplsst, gdouble z);
static gdouble _nc_galaxy_redshift_pop_lsst_srd_ln_eval (NcGalaxyRedshiftPop *gsdrplsst, gdouble z);
static void _nc_galaxy_redshift_pop_lsst_srd_set_lim (NcGalaxyRedshiftPop *gsdrplsst, const gdouble z_min, const gdouble z_max);
static void _nc_galaxy_redshift_pop_lsst_srd_get_lim (NcGalaxyRedshiftPop *gsdrplsst, gdouble *z_min, gdouble *z_max);

static void
nc_galaxy_redshift_pop_lsst_srd_class_init (NcGalaxyRedshiftPopLSSTSRDClass *klass)
{
  NcGalaxyRedshiftPopClass *sd_redshift_class = NC_GALAXY_REDSHIFT_POP_CLASS (klass);
  GObjectClass *object_class                     = G_OBJECT_CLASS (klass);
  NcmModelClass *model_class                     = NCM_MODEL_CLASS (klass);

  object_class->dispose  = &_nc_galaxy_redshift_pop_lsst_srd_dispose;
  object_class->finalize = &_nc_galaxy_redshift_pop_lsst_srd_finalize;

  ncm_model_class_set_name_nick (model_class, "LSST SRD Galaxy Distribution", "LSST_SRD");
  ncm_model_class_add_params (model_class, NC_GALAXY_REDSHIFT_POP_LSST_SRD_SPARAM_LEN, 0, PROP_LEN);

  /**
   * NcGalaxyRedshiftPopLSSTSRD:alpha:
   *
   * The redshift exponential slope.
   *
   */
  ncm_model_class_set_sparam (model_class, NC_GALAXY_REDSHIFT_POP_LSST_SRD_ALPHA, "\\alpha", "alpha",
                              1.0e-8, 1.0, 1.0e-2,
                              NC_GALAXY_REDSHIFT_POP_LSST_SRD_DEFAULT_PARAMS_ABSTOL, NC_GALAXY_REDSHIFT_POP_LSST_SRD_DEFAULT_ALPHA,
                              NCM_PARAM_TYPE_FIXED);

  /**
   * NcGalaxyRedshiftPopLSSTSRD:beta:
   *
   * The redshift power law slope.
   *
   */
  ncm_model_class_set_sparam (model_class, NC_GALAXY_REDSHIFT_POP_LSST_SRD_BETA, "\\beta", "beta",
                              1.0e-8, 5.0, 1.0e-1,
                              NC_GALAXY_REDSHIFT_POP_LSST_SRD_DEFAULT_PARAMS_ABSTOL, NC_GALAXY_REDSHIFT_POP_LSST_SRD_DEFAULT_BETA,
                              NCM_PARAM_TYPE_FIXED);

  /**
   * NcGalaxyRedshiftPopLSSTSRD:z0:
   *
   * The redshift pivot.
   *
   */
  ncm_model_class_set_sparam (model_class, NC_GALAXY_REDSHIFT_POP_LSST_SRD_Z0, "z_0", "z0",
                              1.0e-8, 1.0, 1.0e-2,
                              NC_GALAXY_REDSHIFT_POP_LSST_SRD_DEFAULT_PARAMS_ABSTOL, NC_GALAXY_REDSHIFT_POP_LSST_SRD_DEFAULT_Z0,
                              NCM_PARAM_TYPE_FIXED);

  ncm_model_class_check_params_info (model_class);

  sd_redshift_class->gen      = &_nc_galaxy_redshift_pop_lsst_srd_gen;
  sd_redshift_class->eval    = &_nc_galaxy_redshift_pop_lsst_srd_eval;
  sd_redshift_class->ln_eval = &_nc_galaxy_redshift_pop_lsst_srd_ln_eval;
  sd_redshift_class->set_lim  = &_nc_galaxy_redshift_pop_lsst_srd_set_lim;
  sd_redshift_class->get_lim  = &_nc_galaxy_redshift_pop_lsst_srd_get_lim;
}

#define VECTOR (NCM_MODEL (gsdrp))
#define ALPHA  (ncm_model_orig_param_get (VECTOR, NC_GALAXY_REDSHIFT_POP_LSST_SRD_ALPHA))
#define BETA   (ncm_model_orig_param_get (VECTOR, NC_GALAXY_REDSHIFT_POP_LSST_SRD_BETA))
#define Z0     (ncm_model_orig_param_get (VECTOR, NC_GALAXY_REDSHIFT_POP_LSST_SRD_Z0))

static void
_nc_galaxy_redshift_pop_lsst_srd_update (NcGalaxyRedshiftPop *gsdrp)
{
  NcGalaxyRedshiftPopLSSTSRD *gsdrplsst = NC_GALAXY_REDSHIFT_POP_LSST_SRD (gsdrp);
  NcmModel *model                          = NCM_MODEL (gsdrp);

  if (ncm_model_state_is_update (model))
    return;

  {
    NcGalaxyRedshiftPopLSSTSRDPrivate * const self = nc_galaxy_redshift_pop_lsst_srd_get_instance_private (gsdrplsst);

    const gdouble alpha = ALPHA;
    const gdouble beta  = BETA;
    const gdouble z0    = Z0;
    const gdouble y_low = pow (self->z_min, alpha);
    const gdouble y_up  = pow (self->z_max, alpha);

    self->alpha   = alpha;
    self->beta    = beta;
    self->z0      = z0;
    self->gamma_a = (1.0 + beta) / alpha;
    self->y0      = pow (z0, alpha);
    self->z_norm  = alpha / (pow (z0, 1.0 + self->beta) * (gsl_sf_gamma_inc (self->gamma_a, y_low / self->y0) -
                                                           gsl_sf_gamma_inc (self->gamma_a, y_up / self->y0))
    );
    self->ln_z_norm = log (self->z_norm);

    ncm_model_state_set_update (model);
  }
}

static gdouble
_nc_galaxy_redshift_pop_lsst_srd_gen (NcGalaxyRedshiftPop *gsdrp, NcmRNG *rng)
{
  NcGalaxyRedshiftPopLSSTSRD *gsdrplsst          = NC_GALAXY_REDSHIFT_POP_LSST_SRD (gsdrp);
  NcGalaxyRedshiftPopLSSTSRDPrivate * const self = nc_galaxy_redshift_pop_lsst_srd_get_instance_private (gsdrplsst);
  gdouble z;

  _nc_galaxy_redshift_pop_lsst_srd_update (gsdrp);

  do {
    const gdouble gen_y = ncm_rng_gamma_gen (rng, self->gamma_a, self->y0);

    z = pow (gen_y, 1.0 / self->alpha);
  } while (z < self->z_min || z > self->z_max);

  return z;
}

static gdouble
_nc_galaxy_redshift_pop_lsst_srd_ln_eval (NcGalaxyRedshiftPop *gsdrp, gdouble z)
{
  NcGalaxyRedshiftPopLSSTSRD *gsdrplsst          = NC_GALAXY_REDSHIFT_POP_LSST_SRD (gsdrp);
  NcGalaxyRedshiftPopLSSTSRDPrivate * const self = nc_galaxy_redshift_pop_lsst_srd_get_instance_private (gsdrplsst);

  _nc_galaxy_redshift_pop_lsst_srd_update (gsdrp);

  {
    const gdouble y = pow (z, self->alpha);

    return self->beta * log (z) - (y / self->y0) + self->ln_z_norm;
  }
}

static gdouble
_nc_galaxy_redshift_pop_lsst_srd_eval (NcGalaxyRedshiftPop *gsdrp, gdouble z)
{
  NcGalaxyRedshiftPopLSSTSRD *gsdrplsst          = NC_GALAXY_REDSHIFT_POP_LSST_SRD (gsdrp);
  NcGalaxyRedshiftPopLSSTSRDPrivate * const self = nc_galaxy_redshift_pop_lsst_srd_get_instance_private (gsdrplsst);

  _nc_galaxy_redshift_pop_lsst_srd_update (gsdrp);

  {
    const gdouble y = pow (z, self->alpha);

    return pow (z, self->beta) * exp (-(y / self->y0)) * self->z_norm;
  }
}

static void
_nc_galaxy_redshift_pop_lsst_srd_set_lim (NcGalaxyRedshiftPop *gsdrp, const gdouble z_min, const gdouble z_max)
{
  NcGalaxyRedshiftPopLSSTSRD *gsdrplsst          = NC_GALAXY_REDSHIFT_POP_LSST_SRD (gsdrp);
  NcGalaxyRedshiftPopLSSTSRDPrivate * const self = nc_galaxy_redshift_pop_lsst_srd_get_instance_private (gsdrplsst);

  g_assert_cmpfloat (z_min, <, z_max);

  self->z_min = z_min;
  self->z_max = z_max;

  ncm_model_state_mark_outdated (NCM_MODEL (gsdrp));
}

static void
_nc_galaxy_redshift_pop_lsst_srd_get_lim (NcGalaxyRedshiftPop *gsdrp, gdouble *z_min, gdouble *z_max)
{
  NcGalaxyRedshiftPopLSSTSRD *gsdrplsst          = NC_GALAXY_REDSHIFT_POP_LSST_SRD (gsdrp);
  NcGalaxyRedshiftPopLSSTSRDPrivate * const self = nc_galaxy_redshift_pop_lsst_srd_get_instance_private (gsdrplsst);

  g_assert_nonnull (z_min);
  g_assert_nonnull (z_max);

  *z_min = self->z_min;
  *z_max = self->z_max;
}

/**
 * nc_galaxy_redshift_pop_lsst_srd_new:
 *
 * Creates a new #NcGalaxyRedshiftPopLSSTSRD. The default parameter values correspond to the
 * LSST SRD year 1 source parametrization.
 *
 * Returns: (transfer full): a new #NcGalaxyRedshiftPopLSSTSRD
 */
NcGalaxyRedshiftPopLSSTSRD *
nc_galaxy_redshift_pop_lsst_srd_new (void)
{
  NcmDTuple2 lim = NCM_DTUPLE2_STATIC_INIT (NC_GALAXY_REDSHIFT_POP_LSST_SRD_DEFAULT_Z_LOW,
                                            NC_GALAXY_REDSHIFT_POP_LSST_SRD_DEFAULT_Z_HIGH);
  NcGalaxyRedshiftPopLSSTSRD *gsdrplsst = g_object_new (NC_TYPE_GALAXY_REDSHIFT_POP_LSST_SRD,
                                                           "lim", &lim,
                                                           NULL);

  return gsdrplsst;
}

/**
 * nc_galaxy_redshift_pop_lsst_srd_new_y1_source:
 *
 * Creates a new #NcGalaxyRedshiftPopLSSTSRD. The parameter values correspond to the
 * LSST SRD year 1 source parametrization.
 *
 * Returns: (transfer full): a new #NcGalaxyRedshiftPopLSSTSRD
 */
NcGalaxyRedshiftPopLSSTSRD *
nc_galaxy_redshift_pop_lsst_srd_new_y1_source (void)
{
  NcmDTuple2 lim = NCM_DTUPLE2_STATIC_INIT (NC_GALAXY_REDSHIFT_POP_LSST_SRD_DEFAULT_Z_LOW,
                                            NC_GALAXY_REDSHIFT_POP_LSST_SRD_DEFAULT_Z_HIGH);
  NcGalaxyRedshiftPopLSSTSRD *gsdrplsst = g_object_new (NC_TYPE_GALAXY_REDSHIFT_POP_LSST_SRD,
                                                           "lim", &lim,
                                                           "alpha", NC_GALAXY_REDSHIFT_POP_LSST_SRD_Y1_SOURCE_ALPHA,
                                                           "beta", NC_GALAXY_REDSHIFT_POP_LSST_SRD_Y1_SOURCE_BETA,
                                                           "z0", NC_GALAXY_REDSHIFT_POP_LSST_SRD_Y1_SOURCE_Z0,
                                                           NULL);

  return gsdrplsst;
}

/**
 * nc_galaxy_redshift_pop_lsst_srd_new_y1_lens:
 *
 * Creates a new #NcGalaxyRedshiftPopLSSTSRD. The parameter values correspond to the
 * LSST SRD year 1 lens parametrization.
 *
 * Returns: (transfer full): a new #NcGalaxyRedshiftPopLSSTSRD
 */
NcGalaxyRedshiftPopLSSTSRD *
nc_galaxy_redshift_pop_lsst_srd_new_y1_lens (void)
{
  NcmDTuple2 lim = NCM_DTUPLE2_STATIC_INIT (NC_GALAXY_REDSHIFT_POP_LSST_SRD_DEFAULT_Z_LOW,
                                            NC_GALAXY_REDSHIFT_POP_LSST_SRD_DEFAULT_Z_HIGH);
  NcGalaxyRedshiftPopLSSTSRD *gsdrplsst = g_object_new (NC_TYPE_GALAXY_REDSHIFT_POP_LSST_SRD,
                                                           "lim", &lim,
                                                           "alpha", NC_GALAXY_REDSHIFT_POP_LSST_SRD_Y1_LENS_ALPHA,
                                                           "beta", NC_GALAXY_REDSHIFT_POP_LSST_SRD_Y1_LENS_BETA,
                                                           "z0", NC_GALAXY_REDSHIFT_POP_LSST_SRD_Y1_LENS_Z0,
                                                           NULL);

  return gsdrplsst;
}

/**
 * nc_galaxy_redshift_pop_lsst_srd_new_y10_source:
 *
 * Creates a new #NcGalaxyRedshiftPopLSSTSRD. The parameter values correspond to the
 * LSST SRD year 10 source parametrization.
 *
 * Returns: (transfer full): a new #NcGalaxyRedshiftPopLSSTSRD
 */
NcGalaxyRedshiftPopLSSTSRD *
nc_galaxy_redshift_pop_lsst_srd_new_y10_source (void)
{
  NcmDTuple2 lim = NCM_DTUPLE2_STATIC_INIT (NC_GALAXY_REDSHIFT_POP_LSST_SRD_DEFAULT_Z_LOW,
                                            NC_GALAXY_REDSHIFT_POP_LSST_SRD_DEFAULT_Z_HIGH);
  NcGalaxyRedshiftPopLSSTSRD *gsdrplsst = g_object_new (NC_TYPE_GALAXY_REDSHIFT_POP_LSST_SRD,
                                                           "lim", &lim,
                                                           "alpha", NC_GALAXY_REDSHIFT_POP_LSST_SRD_Y10_SOURCE_ALPHA,
                                                           "beta", NC_GALAXY_REDSHIFT_POP_LSST_SRD_Y10_SOURCE_BETA,
                                                           "z0", NC_GALAXY_REDSHIFT_POP_LSST_SRD_Y10_SOURCE_Z0,
                                                           NULL);

  return gsdrplsst;
}

/**
 * nc_galaxy_redshift_pop_lsst_srd_new_y10_lens:
 *
 * Creates a new #NcGalaxyRedshiftPopLSSTSRD. The parameter values correspond to the
 * LSST SRD year 10 lens parametrization.
 *
 * Returns: (transfer full): a new #NcGalaxyRedshiftPopLSSTSRD
 */
NcGalaxyRedshiftPopLSSTSRD *
nc_galaxy_redshift_pop_lsst_srd_new_y10_lens (void)
{
  NcmDTuple2 lim = NCM_DTUPLE2_STATIC_INIT (NC_GALAXY_REDSHIFT_POP_LSST_SRD_DEFAULT_Z_LOW,
                                            NC_GALAXY_REDSHIFT_POP_LSST_SRD_DEFAULT_Z_HIGH);
  NcGalaxyRedshiftPopLSSTSRD *gsdrplsst = g_object_new (NC_TYPE_GALAXY_REDSHIFT_POP_LSST_SRD,
                                                           "lim", &lim,
                                                           "alpha", NC_GALAXY_REDSHIFT_POP_LSST_SRD_Y10_LENS_ALPHA,
                                                           "beta", NC_GALAXY_REDSHIFT_POP_LSST_SRD_Y10_LENS_BETA,
                                                           "z0", NC_GALAXY_REDSHIFT_POP_LSST_SRD_Y10_LENS_Z0,
                                                           NULL);

  return gsdrplsst;
}

/**
 * nc_galaxy_redshift_pop_lsst_srd_new_from_type:
 * @type: a #NcGalaxyRedshiftPopLSSTSRDType
 *
 * Creates a new #NcGalaxyRedshiftPopLSSTSRD using a predefined type.
 * The type determines which LSST SRD parametrization to use (Year 1 or Year 10,
 * source or lens).
 *
 * Returns: (transfer full): a new #NcGalaxyRedshiftPopLSSTSRD
 */
NcGalaxyRedshiftPopLSSTSRD *
nc_galaxy_redshift_pop_lsst_srd_new_from_type (NcGalaxyRedshiftPopLSSTSRDType type)
{
  switch (type)
  {
    case NC_GALAXY_REDSHIFT_POP_LSST_SRD_Y1_SOURCE:
      return nc_galaxy_redshift_pop_lsst_srd_new_y1_source ();

    case NC_GALAXY_REDSHIFT_POP_LSST_SRD_Y1_LENS:
      return nc_galaxy_redshift_pop_lsst_srd_new_y1_lens ();

    case NC_GALAXY_REDSHIFT_POP_LSST_SRD_Y10_SOURCE:
      return nc_galaxy_redshift_pop_lsst_srd_new_y10_source ();

    case NC_GALAXY_REDSHIFT_POP_LSST_SRD_Y10_LENS:
      return nc_galaxy_redshift_pop_lsst_srd_new_y10_lens ();

    default:
      g_error ("nc_galaxy_redshift_pop_lsst_srd_new_from_type: invalid type %d", type);

      return NULL;
  }
}

/**
 * nc_galaxy_redshift_pop_lsst_srd_ref:
 * @gsdrplsst: a #NcGalaxyRedshiftPopLSSTSRD
 *
 * Increases the reference count of @gsdrplsst by one.
 *
 * Returns: (transfer full): @gsdrplsst.
 */
NcGalaxyRedshiftPopLSSTSRD *
nc_galaxy_redshift_pop_lsst_srd_ref (NcGalaxyRedshiftPopLSSTSRD *gsdrplsst)
{
  return g_object_ref (gsdrplsst);
}

/**
 * nc_galaxy_redshift_pop_lsst_srd_free:
 * @gsdrplsst: a #NcGalaxyRedshiftPopLSSTSRD
 *
 * Decreases the reference count of @gsdrplsst by one.
 *
 */
void
nc_galaxy_redshift_pop_lsst_srd_free (NcGalaxyRedshiftPopLSSTSRD *gsdrplsst)
{
  g_object_unref (gsdrplsst);
}

/**
 * nc_galaxy_redshift_pop_lsst_srd_clear:
 * @gsdrplsst: a #NcGalaxyRedshiftPopLSSTSRD
 *
 * Decreases the reference count of @gsdrplsst by one, and sets the
 * pointer @gsdrplsst to NULL.
 *
 */
void
nc_galaxy_redshift_pop_lsst_srd_clear (NcGalaxyRedshiftPopLSSTSRD **gsdrplsst)
{
  g_clear_object (gsdrplsst);
}


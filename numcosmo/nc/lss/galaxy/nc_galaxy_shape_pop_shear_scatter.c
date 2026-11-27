/***************************************************************************
 *            nc_galaxy_shape_pop_shear_scatter.c
 *
 *  Thu Nov 26 2026
 *  Copyright  2026  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * nc_galaxy_shape_pop_shear_scatter.c
 * Copyright (C) 2026 Caio Lima de Oliveira <caiolimadeoliveira@pm.me>
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
 * NcGalaxyShapePopShearScatter:
 *
 * Abstract scatter of a linear shear estimator about the true shear.
 *
 * Unlike the other #NcGalaxyShapePop models, which describe the intrinsic
 * ellipticity $\chi_I$ on the unit disc, a shear-scatter population describes
 * the one-dimensional distribution of the residual
 * $\delta = \hat g_t - \mu_t$ of a galaxy's tangential shear estimate
 * $\hat g_t$ about its expected value $\mu_t$, as in the "P(z) method" of
 * Applegate et al. (2014) (Weighing the Giants III). It lives on the real
 * line, in reduced-shear units, and is consumed only by
 * #NcGalaxyShapeFactorWtG; every other shape factor refuses it, and the
 * intrinsic-ellipticity methods of #NcGalaxyShapePop (eval_p(), gen(),
 * e_rms(), the radial moments and exponents) abort when called on it.
 *
 * The parameters of the concrete models are in reduced-shear units. Each
 * evaluation also receives the per-galaxy measurement-noise variance
 * @var_noise_g, already converted to shear units, which is added in
 * quadrature to the Gaussian core of the kernel when
 * #NcGalaxyShapePopShearScatter:include-noise is %TRUE (the default), and the
 * factor @eps_to_g $= 1/(sR)$ that converts an ellipticity into the shear
 * estimate, for models that take a width from an ellipticity-unit catalog
 * column.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc/lss/galaxy/nc_galaxy_shape_pop_shear_scatter.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <math.h>
#endif /* NUMCOSMO_GIR_SCAN */

typedef struct _NcGalaxyShapePopShearScatterPrivate
{
  gboolean include_noise;
} NcGalaxyShapePopShearScatterPrivate;

enum
{
  PROP_0,
  PROP_INCLUDE_NOISE,
  PROP_LEN,
};

G_DEFINE_ABSTRACT_TYPE_WITH_PRIVATE (NcGalaxyShapePopShearScatter, nc_galaxy_shape_pop_shear_scatter, NC_TYPE_GALAXY_SHAPE_POP)

static void
nc_galaxy_shape_pop_shear_scatter_init (NcGalaxyShapePopShearScatter *gspss)
{
  NcGalaxyShapePopShearScatterPrivate * const self = nc_galaxy_shape_pop_shear_scatter_get_instance_private (gspss);

  self->include_noise = TRUE;
}

static void
_nc_galaxy_shape_pop_shear_scatter_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcGalaxyShapePopShearScatter *gspss = NC_GALAXY_SHAPE_POP_SHEAR_SCATTER (object);

  g_return_if_fail (NC_IS_GALAXY_SHAPE_POP_SHEAR_SCATTER (gspss));

  switch (prop_id)
  {
    case PROP_INCLUDE_NOISE:
      nc_galaxy_shape_pop_shear_scatter_set_include_noise (gspss, g_value_get_boolean (value));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_galaxy_shape_pop_shear_scatter_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcGalaxyShapePopShearScatter *gspss = NC_GALAXY_SHAPE_POP_SHEAR_SCATTER (object);

  g_return_if_fail (NC_IS_GALAXY_SHAPE_POP_SHEAR_SCATTER (gspss));

  switch (prop_id)
  {
    case PROP_INCLUDE_NOISE:
      g_value_set_boolean (value, nc_galaxy_shape_pop_shear_scatter_get_include_noise (gspss));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_galaxy_shape_pop_shear_scatter_finalize (GObject *object)
{
  G_OBJECT_CLASS (nc_galaxy_shape_pop_shear_scatter_parent_class)->finalize (object);
}

/* LCOV_EXCL_START */
#define _SHEAR_SCATTER_REFUSE(name)                                                       \
        g_error ("%s: %s is a shear-scatter population (a distribution of the shear "     \
                 "estimate, not of the intrinsic ellipticity); it can only be used with " \
                 "NcGalaxyShapeFactorWtG.", name, G_OBJECT_TYPE_NAME (gsp))

static gdouble
_nc_galaxy_shape_pop_shear_scatter_eval_p (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data, const gdouble r)
{
  _SHEAR_SCATTER_REFUSE ("nc_galaxy_shape_pop_eval_p");

  return 0.0;
}

static void
_nc_galaxy_shape_pop_shear_scatter_eval_p_array (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data, const GArray *r, GArray **p)
{
  _SHEAR_SCATTER_REFUSE ("nc_galaxy_shape_pop_eval_p_array");
}

static void
_nc_galaxy_shape_pop_shear_scatter_gen (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data, NcmRNG *rng, gdouble *e_int_1, gdouble *e_int_2)
{
  _SHEAR_SCATTER_REFUSE ("nc_galaxy_shape_pop_gen");
}

static gdouble
_nc_galaxy_shape_pop_shear_scatter_e_rms (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data)
{
  _SHEAR_SCATTER_REFUSE ("nc_galaxy_shape_pop_e_rms");

  return 0.0;
}

static gdouble
_nc_galaxy_shape_pop_shear_scatter_exponent (NcGalaxyShapePop *gsp)
{
  _SHEAR_SCATTER_REFUSE ("nc_galaxy_shape_pop_exponent_at_origin/one");

  return 0.0;
}

static gdouble
_nc_galaxy_shape_pop_shear_scatter_moment_2k (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data, const guint k)
{
  _SHEAR_SCATTER_REFUSE ("nc_galaxy_shape_pop_moment_2k");

  return 0.0;
}

static gdouble
_nc_galaxy_shape_pop_shear_scatter_eval_p_reduced (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data, const gdouble r)
{
  _SHEAR_SCATTER_REFUSE ("nc_galaxy_shape_pop_eval_p_reduced");

  return 0.0;
}

static gdouble
_nc_galaxy_shape_pop_shear_scatter_eval_ln_p_default (NcGalaxyShapePopShearScatter *gspss, NcGalaxyShapePopData *data, const gdouble delta_g, const gdouble var_noise_g, const gdouble eps_to_g)
{
  g_error ("_nc_galaxy_shape_pop_shear_scatter_eval_ln_p: not implemented by `%s'.", G_OBJECT_TYPE_NAME (gspss));

  return 0.0;
}

static gdouble
_nc_galaxy_shape_pop_shear_scatter_gen_delta_default (NcGalaxyShapePopShearScatter *gspss, NcGalaxyShapePopData *data, NcmRNG *rng, const gdouble var_noise_g, const gdouble eps_to_g)
{
  g_error ("_nc_galaxy_shape_pop_shear_scatter_gen_delta: not implemented by `%s'.", G_OBJECT_TYPE_NAME (gspss));

  return 0.0;
}

#undef _SHEAR_SCATTER_REFUSE
/* LCOV_EXCL_STOP */

static void
nc_galaxy_shape_pop_shear_scatter_class_init (NcGalaxyShapePopShearScatterClass *klass)
{
  NcGalaxyShapePopClass *gsp_class = NC_GALAXY_SHAPE_POP_CLASS (klass);
  GObjectClass *object_class       = G_OBJECT_CLASS (klass);
  NcmModelClass *model_class       = NCM_MODEL_CLASS (klass);

  model_class->set_property = &_nc_galaxy_shape_pop_shear_scatter_set_property;
  model_class->get_property = &_nc_galaxy_shape_pop_shear_scatter_get_property;
  object_class->finalize    = &_nc_galaxy_shape_pop_shear_scatter_finalize;

  ncm_model_class_set_name_nick (model_class, "Shear-estimator scatter", "ShearScatter");
  ncm_model_class_add_params (model_class, 0, 0, PROP_LEN);

  /**
   * NcGalaxyShapePopShearScatter:include-noise:
   *
   * Whether the per-galaxy measurement noise (the catalog's std_noise column,
   * converted to shear units) is added in quadrature to the Gaussian core of
   * the kernel. Set it to %FALSE when the kernel parameters already describe
   * the whole scatter of the shear estimates, as in Applegate et al. (2014),
   * where the Voigt width is fitted to the observed distribution.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_INCLUDE_NOISE,
                                   g_param_spec_boolean ("include-noise",
                                                         "Include noise",
                                                         "Whether the measurement noise is added to the kernel core",
                                                         TRUE,
                                                         G_PARAM_READWRITE | G_PARAM_STATIC_STRINGS));

  ncm_model_class_check_params_info (model_class);

  gsp_class->eval_p             = &_nc_galaxy_shape_pop_shear_scatter_eval_p;
  gsp_class->eval_p_array       = &_nc_galaxy_shape_pop_shear_scatter_eval_p_array;
  gsp_class->gen                = &_nc_galaxy_shape_pop_shear_scatter_gen;
  gsp_class->e_rms              = &_nc_galaxy_shape_pop_shear_scatter_e_rms;
  gsp_class->exponent_at_origin = &_nc_galaxy_shape_pop_shear_scatter_exponent;
  gsp_class->exponent_at_one    = &_nc_galaxy_shape_pop_shear_scatter_exponent;
  gsp_class->moment_2k          = &_nc_galaxy_shape_pop_shear_scatter_moment_2k;
  gsp_class->eval_p_reduced     = &_nc_galaxy_shape_pop_shear_scatter_eval_p_reduced;

  klass->eval_ln_p = &_nc_galaxy_shape_pop_shear_scatter_eval_ln_p_default;
  klass->gen_delta = &_nc_galaxy_shape_pop_shear_scatter_gen_delta_default;
}

/**
 * nc_galaxy_shape_pop_shear_scatter_ref:
 * @gspss: a #NcGalaxyShapePopShearScatter
 *
 * Increases the reference count of @gspss by one.
 *
 * Returns: (transfer full): @gspss.
 */
NcGalaxyShapePopShearScatter *
nc_galaxy_shape_pop_shear_scatter_ref (NcGalaxyShapePopShearScatter *gspss)
{
  return g_object_ref (gspss);
}

/**
 * nc_galaxy_shape_pop_shear_scatter_free:
 * @gspss: a #NcGalaxyShapePopShearScatter
 *
 * Decreases the reference count of @gspss by one.
 *
 */
void
nc_galaxy_shape_pop_shear_scatter_free (NcGalaxyShapePopShearScatter *gspss)
{
  g_object_unref (gspss);
}

/**
 * nc_galaxy_shape_pop_shear_scatter_clear:
 * @gspss: a #NcGalaxyShapePopShearScatter
 *
 * Decreases the reference count of *@gspss by one, and sets the pointer
 * *@gspss to NULL.
 *
 */
void
nc_galaxy_shape_pop_shear_scatter_clear (NcGalaxyShapePopShearScatter **gspss)
{
  g_clear_object (gspss);
}

/**
 * nc_galaxy_shape_pop_shear_scatter_set_include_noise:
 * @gspss: a #NcGalaxyShapePopShearScatter
 * @include_noise: whether the measurement noise is added to the kernel core
 *
 * Sets #NcGalaxyShapePopShearScatter:include-noise.
 *
 */
void
nc_galaxy_shape_pop_shear_scatter_set_include_noise (NcGalaxyShapePopShearScatter *gspss, gboolean include_noise)
{
  NcGalaxyShapePopShearScatterPrivate * const self = nc_galaxy_shape_pop_shear_scatter_get_instance_private (gspss);

  if (self->include_noise != include_noise)
  {
    self->include_noise = include_noise;

    /* Moves the pkey, which the shape factor hashes, so cached values built
     * with the previous setting are rebuilt. */
    ncm_model_state_mark_outdated (NCM_MODEL (gspss));
  }
}

/**
 * nc_galaxy_shape_pop_shear_scatter_get_include_noise:
 * @gspss: a #NcGalaxyShapePopShearScatter
 *
 * Returns: whether the measurement noise is added to the kernel core.
 */
gboolean
nc_galaxy_shape_pop_shear_scatter_get_include_noise (NcGalaxyShapePopShearScatter *gspss)
{
  NcGalaxyShapePopShearScatterPrivate * const self = nc_galaxy_shape_pop_shear_scatter_get_instance_private (gspss);

  return self->include_noise;
}

/**
 * nc_galaxy_shape_pop_shear_scatter_eval_ln_p:
 * @gspss: a #NcGalaxyShapePopShearScatter
 * @data: a #NcGalaxyShapePopData, prepared with nc_galaxy_shape_pop_prepare()
 * @delta_g: the residual $\hat g_t - \mu_t$ of the shear estimate
 * @var_noise_g: the measurement-noise variance in shear units
 * @eps_to_g: the factor $1/(sR)$ converting an ellipticity into a shear estimate
 *
 * Evaluates the logarithm of the normalized one-dimensional density of
 * @delta_g. @var_noise_g is ignored when
 * #NcGalaxyShapePopShearScatter:include-noise is %FALSE.
 *
 * Returns: $\ln P(\delta_g)$.
 */
gdouble
nc_galaxy_shape_pop_shear_scatter_eval_ln_p (NcGalaxyShapePopShearScatter *gspss, NcGalaxyShapePopData *data, const gdouble delta_g, const gdouble var_noise_g, const gdouble eps_to_g)
{
  NcGalaxyShapePopShearScatterPrivate * const self = nc_galaxy_shape_pop_shear_scatter_get_instance_private (gspss);

  return NC_GALAXY_SHAPE_POP_SHEAR_SCATTER_GET_CLASS (gspss)->eval_ln_p (gspss, data, delta_g, self->include_noise ? var_noise_g : 0.0, eps_to_g);
}

/**
 * nc_galaxy_shape_pop_shear_scatter_eval_p:
 * @gspss: a #NcGalaxyShapePopShearScatter
 * @data: a #NcGalaxyShapePopData, prepared with nc_galaxy_shape_pop_prepare()
 * @delta_g: the residual $\hat g_t - \mu_t$ of the shear estimate
 * @var_noise_g: the measurement-noise variance in shear units
 * @eps_to_g: the factor $1/(sR)$ converting an ellipticity into a shear estimate
 *
 * Returns: $P(\delta_g)$, the exponential of nc_galaxy_shape_pop_shear_scatter_eval_ln_p().
 */
gdouble
nc_galaxy_shape_pop_shear_scatter_eval_p (NcGalaxyShapePopShearScatter *gspss, NcGalaxyShapePopData *data, const gdouble delta_g, const gdouble var_noise_g, const gdouble eps_to_g)
{
  return exp (nc_galaxy_shape_pop_shear_scatter_eval_ln_p (gspss, data, delta_g, var_noise_g, eps_to_g));
}

/**
 * nc_galaxy_shape_pop_shear_scatter_gen_delta:
 * @gspss: a #NcGalaxyShapePopShearScatter
 * @data: a #NcGalaxyShapePopData, prepared with nc_galaxy_shape_pop_prepare()
 * @rng: a #NcmRNG
 * @var_noise_g: the measurement-noise variance in shear units
 * @eps_to_g: the factor $1/(sR)$ converting an ellipticity into a shear estimate
 *
 * Draws a residual $\delta_g$ from the density of
 * nc_galaxy_shape_pop_shear_scatter_eval_ln_p().
 *
 * Returns: the drawn residual.
 */
gdouble
nc_galaxy_shape_pop_shear_scatter_gen_delta (NcGalaxyShapePopShearScatter *gspss, NcGalaxyShapePopData *data, NcmRNG *rng, const gdouble var_noise_g, const gdouble eps_to_g)
{
  NcGalaxyShapePopShearScatterPrivate * const self = nc_galaxy_shape_pop_shear_scatter_get_instance_private (gspss);

  return NC_GALAXY_SHAPE_POP_SHEAR_SCATTER_GET_CLASS (gspss)->gen_delta (gspss, data, rng, self->include_noise ? var_noise_g : 0.0, eps_to_g);
}


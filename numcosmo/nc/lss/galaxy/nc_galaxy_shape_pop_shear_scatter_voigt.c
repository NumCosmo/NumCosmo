/***************************************************************************
 *            nc_galaxy_shape_pop_shear_scatter_voigt.c
 *
 *  Thu Nov 26 2026
 *  Copyright  2026  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * nc_galaxy_shape_pop_shear_scatter_voigt.c
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
 * NcGalaxyShapePopShearScatterVoigt:
 *
 * Voigt scatter of the tangential shear estimate.
 *
 * The residual $\delta = \hat g_t - \mu_t$ follows a Voigt profile, the
 * convolution of a Gaussian of width $\sigma_t$ with a Lorentzian of half
 * width at half maximum $\Gamma$,
 * $$V(\delta; \sigma_t, \Gamma) = \frac{\mathrm{Re}\,w(z)}{\sigma_t\sqrt{2\pi}},
 *   \qquad z = \frac{\delta + i\Gamma}{\sigma_t\sqrt{2}},$$
 * where $w$ is the Faddeeva function. This is the model Applegate et al.
 * (2014) (Weighing the Giants III, their Section 5) adopt for the scatter of
 * the shear estimates, whose wings a Gaussian misses. The Gaussian core is
 * $\sigma_t^2 = \sigma_g^2 + \sigma_{\mathrm{noise},g}^2$, with the
 * measurement noise dropped when #NcGalaxyShapePopShearScatter:include-noise
 * is %FALSE (the Applegate et al. setting, where $\sigma_g$ and $\Gamma$ are
 * fitted to the whole distribution of the shear estimates). Both parameters
 * are in reduced-shear units, as in #NcReducedShearClusterMass.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc/lss/galaxy/nc_galaxy_shape_pop_shear_scatter_voigt.h"
#include "external/misc/Faddeeva.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <math.h>
#include <complex.h>
#endif /* NUMCOSMO_GIR_SCAN */

struct _NcGalaxyShapePopShearScatterVoigt
{
  NcGalaxyShapePopShearScatter parent_instance;
};

enum
{
  PROP_0,
  PROP_LEN,
};

G_DEFINE_TYPE (NcGalaxyShapePopShearScatterVoigt, nc_galaxy_shape_pop_shear_scatter_voigt, NC_TYPE_GALAXY_SHAPE_POP_SHEAR_SCATTER)

static void
nc_galaxy_shape_pop_shear_scatter_voigt_init (NcGalaxyShapePopShearScatterVoigt *gspssv)
{
}

static void
_nc_galaxy_shape_pop_shear_scatter_voigt_finalize (GObject *object)
{
  G_OBJECT_CLASS (nc_galaxy_shape_pop_shear_scatter_voigt_parent_class)->finalize (object);
}

#define VECTOR (NCM_MODEL (gspss))
#define SIGMA  (ncm_model_orig_param_get (VECTOR, NC_GALAXY_SHAPE_POP_SHEAR_SCATTER_VOIGT_SIGMA))
#define GAMMA  (ncm_model_orig_param_get (VECTOR, NC_GALAXY_SHAPE_POP_SHEAR_SCATTER_VOIGT_GAMMA))

static void
_nc_galaxy_shape_pop_shear_scatter_voigt_ldata_noop (NcGalaxyShapePopData *data, NcGalaxyWLObs *obs, const guint i)
{
}

static void
_nc_galaxy_shape_pop_shear_scatter_voigt_ldata_no_columns (NcGalaxyShapePopData *data, GList **columns)
{
}

static void
_nc_galaxy_shape_pop_shear_scatter_voigt_data_init (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data)
{
  data->ldata                  = NULL;
  data->ldata_destroy          = &g_free;
  data->ldata_read_row         = &_nc_galaxy_shape_pop_shear_scatter_voigt_ldata_noop;
  data->ldata_write_row        = &_nc_galaxy_shape_pop_shear_scatter_voigt_ldata_noop;
  data->ldata_required_columns = &_nc_galaxy_shape_pop_shear_scatter_voigt_ldata_no_columns;
}

static void
_nc_galaxy_shape_pop_shear_scatter_voigt_prepare (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data)
{
  /* Nothing to resolve: the parameters are read at evaluation. */
}

static gdouble
_nc_galaxy_shape_pop_shear_scatter_voigt_eval_ln_p (NcGalaxyShapePopShearScatter *gspss, NcGalaxyShapePopData *data, const gdouble delta_g, const gdouble var_noise_g, const gdouble eps_to_g)
{
  const gdouble sigma    = SIGMA;
  const gdouble sigma_t  = sqrt (sigma * sigma + var_noise_g);
  const complex double z = (delta_g + I * GAMMA) / (M_SQRT2 * sigma_t);

  /* Re w (z) > 0 for Im z > 0, and decays as the Lorentzian wing
   * Im z / (sqrt (pi) |z|^2) far from the core, so the logarithm stays finite. */
  return log (creal (Faddeeva_w (z, 0.0))) - log (M_SQRT2 * M_SQRTPI * sigma_t);
}

static gdouble
_nc_galaxy_shape_pop_shear_scatter_voigt_gen_delta (NcGalaxyShapePopShearScatter *gspss, NcGalaxyShapePopData *data, NcmRNG *rng, const gdouble var_noise_g, const gdouble eps_to_g)
{
  const gdouble sigma   = SIGMA;
  const gdouble sigma_t = sqrt (sigma * sigma + var_noise_g);
  const gdouble gauss   = ncm_rng_gaussian_gen (rng, 0.0, sigma_t);
  const gdouble u       = ncm_rng_uniform01_pos_gen (rng);

  /* The Voigt variable is the sum of a Gaussian and a Cauchy variable. */
  return gauss + GAMMA * tan (M_PI * (u - 0.5));
}

static void
nc_galaxy_shape_pop_shear_scatter_voigt_class_init (NcGalaxyShapePopShearScatterVoigtClass *klass)
{
  NcGalaxyShapePopShearScatterClass *gspss_class = NC_GALAXY_SHAPE_POP_SHEAR_SCATTER_CLASS (klass);
  NcGalaxyShapePopClass *gsp_class               = NC_GALAXY_SHAPE_POP_CLASS (klass);
  GObjectClass *object_class                     = G_OBJECT_CLASS (klass);
  NcmModelClass *model_class                     = NCM_MODEL_CLASS (klass);

  object_class->finalize = &_nc_galaxy_shape_pop_shear_scatter_voigt_finalize;

  ncm_model_class_set_name_nick (model_class, "Voigt shear-estimator scatter", "ShearScatterVoigt");
  ncm_model_class_add_params (model_class, NC_GALAXY_SHAPE_POP_SHEAR_SCATTER_VOIGT_SPARAM_LEN, 0, PROP_LEN);

  /**
   * NcGalaxyShapePopShearScatterVoigt:sigma:
   *
   * The width $\sigma_g$ of the Gaussian core of the Voigt scatter of the
   * tangential shear estimate, in reduced-shear units.
   *
   */
  ncm_model_class_set_sparam (model_class,
                              NC_GALAXY_SHAPE_POP_SHEAR_SCATTER_VOIGT_SIGMA,
                              "\\sigma_g",
                              "sigma", 1.0e-3, 1.0, 1.0e-2,
                              NC_GALAXY_SHAPE_POP_SHEAR_SCATTER_VOIGT_DEFAULT_PARAMS_ABSTOL,
                              NC_GALAXY_SHAPE_POP_SHEAR_SCATTER_VOIGT_DEFAULT_SIGMA,
                              NCM_PARAM_TYPE_FIXED);

  /**
   * NcGalaxyShapePopShearScatterVoigt:gamma:
   *
   * The half width at half maximum $\Gamma$ of the Lorentzian wings of the
   * Voigt scatter of the tangential shear estimate, in reduced-shear units.
   * It must be positive; a vanishing $\Gamma$ is the Gaussian of
   * #NcGalaxyShapePopShearScatterGauss.
   *
   */
  ncm_model_class_set_sparam (model_class,
                              NC_GALAXY_SHAPE_POP_SHEAR_SCATTER_VOIGT_GAMMA,
                              "\\Gamma",
                              "gamma", 1.0e-4, 0.5, 1.0e-3,
                              NC_GALAXY_SHAPE_POP_SHEAR_SCATTER_VOIGT_DEFAULT_PARAMS_ABSTOL,
                              NC_GALAXY_SHAPE_POP_SHEAR_SCATTER_VOIGT_DEFAULT_GAMMA,
                              NCM_PARAM_TYPE_FIXED);

  ncm_model_class_check_params_info (model_class);

  gsp_class->data_init = &_nc_galaxy_shape_pop_shear_scatter_voigt_data_init;
  gsp_class->prepare   = &_nc_galaxy_shape_pop_shear_scatter_voigt_prepare;

  gspss_class->eval_ln_p = &_nc_galaxy_shape_pop_shear_scatter_voigt_eval_ln_p;
  gspss_class->gen_delta = &_nc_galaxy_shape_pop_shear_scatter_voigt_gen_delta;
}

/**
 * nc_galaxy_shape_pop_shear_scatter_voigt_new:
 *
 * Creates a new #NcGalaxyShapePopShearScatterVoigt.
 *
 * Returns: (transfer full): a new #NcGalaxyShapePopShearScatterVoigt.
 */
NcGalaxyShapePopShearScatterVoigt *
nc_galaxy_shape_pop_shear_scatter_voigt_new (void)
{
  return g_object_new (NC_TYPE_GALAXY_SHAPE_POP_SHEAR_SCATTER_VOIGT, NULL);
}

/**
 * nc_galaxy_shape_pop_shear_scatter_voigt_ref:
 * @gspssv: a #NcGalaxyShapePopShearScatterVoigt
 *
 * Increases the reference count of @gspssv by one.
 *
 * Returns: (transfer full): @gspssv.
 */
NcGalaxyShapePopShearScatterVoigt *
nc_galaxy_shape_pop_shear_scatter_voigt_ref (NcGalaxyShapePopShearScatterVoigt *gspssv)
{
  return g_object_ref (gspssv);
}

/**
 * nc_galaxy_shape_pop_shear_scatter_voigt_free:
 * @gspssv: a #NcGalaxyShapePopShearScatterVoigt
 *
 * Decreases the reference count of @gspssv by one.
 *
 */
void
nc_galaxy_shape_pop_shear_scatter_voigt_free (NcGalaxyShapePopShearScatterVoigt *gspssv)
{
  g_object_unref (gspssv);
}

/**
 * nc_galaxy_shape_pop_shear_scatter_voigt_clear:
 * @gspssv: a #NcGalaxyShapePopShearScatterVoigt
 *
 * Decreases the reference count of *@gspssv by one, and sets the pointer
 * *@gspssv to NULL.
 *
 */
void
nc_galaxy_shape_pop_shear_scatter_voigt_clear (NcGalaxyShapePopShearScatterVoigt **gspssv)
{
  g_clear_object (gspssv);
}


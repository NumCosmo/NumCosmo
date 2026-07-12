/***************************************************************************
 *            nc_galaxy_shape_pop_gauss.c
 *
 *  Thu Jun 19 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_galaxy_shape_pop_gauss.c
 * Copyright (C) 2026 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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
 * NcGalaxyShapePopGauss:
 *
 * Truncated-Gaussian intrinsic ellipticity distribution.
 *
 * The intrinsic ellipticity follows an isotropic Gaussian of width $\sigma$
 * truncated to the unit disk. The induced density of $x = |\chi_I|^2$ is
 * $P(x) \propto \exp(-x / 2\sigma^2)$ on $[0,1]$.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc/lss/galaxy/nc_galaxy_shape_pop_gauss.h"
#include "nc/lss/galaxy/nc_galaxy_shape_pop_gauss_private.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <math.h>
#include <complex.h>
#endif /* NUMCOSMO_GIR_SCAN */

struct _NcGalaxyShapePopGauss
{
  NcGalaxyShapePop parent_instance;
};

typedef struct _NcGalaxyShapePopGaussLData
{
  gdouble norm;        /* residual normalization */
  gdouble inv_2sigma2; /* 1 / (2 sigma^2) */
  gdouble sigma;       /* truncated-Gaussian width (for sampling) */
} NcGalaxyShapePopGaussLData;

enum
{
  PROP_0,
  PROP_LEN,
};

G_DEFINE_TYPE (NcGalaxyShapePopGauss, nc_galaxy_shape_pop_gauss, NC_TYPE_GALAXY_SHAPE_POP);

static void
nc_galaxy_shape_pop_gauss_init (NcGalaxyShapePopGauss *gspg)
{
}

static void
_nc_galaxy_shape_pop_gauss_finalize (GObject *object)
{
  G_OBJECT_CLASS (nc_galaxy_shape_pop_gauss_parent_class)->finalize (object);
}

/* _data_init, _eval_p and _gen are shared with sibling models (see
 * nc_galaxy_shape_pop_gauss_private.h) so are not static; forward-declared
 * here only to keep definition order readable. */
void _nc_galaxy_shape_pop_gauss_data_init (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data);
gdouble _nc_galaxy_shape_pop_gauss_eval_p (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data, const gdouble x);
void _nc_galaxy_shape_pop_gauss_gen (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data, NcmRNG *rng, gdouble *e_int_1, gdouble *e_int_2);
static void _nc_galaxy_shape_pop_gauss_prepare (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data);
static gdouble _nc_galaxy_shape_pop_gauss_e_rms (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data);

static void
nc_galaxy_shape_pop_gauss_class_init (NcGalaxyShapePopGaussClass *klass)
{
  NcGalaxyShapePopClass *gsp_class = NC_GALAXY_SHAPE_POP_CLASS (klass);
  GObjectClass *object_class       = G_OBJECT_CLASS (klass);
  NcmModelClass *model_class       = NCM_MODEL_CLASS (klass);

  object_class->finalize = &_nc_galaxy_shape_pop_gauss_finalize;

  ncm_model_class_set_name_nick (model_class, "Truncated Gaussian intrinsic ellipticity distribution", "GaussIntrinsic");
  ncm_model_class_add_params (model_class, NC_GALAXY_SHAPE_POP_GAUSS_SPARAM_LEN, 0, PROP_LEN);

  /**
   * NcGalaxyShapePopGauss:sigma:
   *
   * The width $\sigma$ of the truncated Gaussian intrinsic ellipticity.
   *
   */
  ncm_model_class_set_sparam (model_class,
                              NC_GALAXY_SHAPE_POP_GAUSS_SIGMA,
                              "\\sigma",
                              "sigma", 1.0e-2, 1.0, 1.0e-1,
                              NC_GALAXY_SHAPE_POP_GAUSS_DEFAULT_PARAMS_ABSTOL,
                              NC_GALAXY_SHAPE_POP_GAUSS_DEFAULT_SIGMA,
                              NCM_PARAM_TYPE_FIXED);

  ncm_model_class_check_params_info (model_class);

  gsp_class->data_init = &_nc_galaxy_shape_pop_gauss_data_init;
  gsp_class->prepare   = &_nc_galaxy_shape_pop_gauss_prepare;
  gsp_class->eval_p    = &_nc_galaxy_shape_pop_gauss_eval_p;
  gsp_class->gen       = &_nc_galaxy_shape_pop_gauss_gen;
  gsp_class->e_rms     = &_nc_galaxy_shape_pop_gauss_e_rms;
}

#define VECTOR (NCM_MODEL (gsp))
#define SIGMA  (ncm_model_orig_param_get (VECTOR, NC_GALAXY_SHAPE_POP_GAUSS_SIGMA))

static void
_nc_galaxy_shape_pop_gauss_ldata_noop (NcGalaxyShapePopData *data, NcGalaxyWLObs *obs, const guint i)
{
}

static void
_nc_galaxy_shape_pop_gauss_ldata_required_columns (NcGalaxyShapePopData *data, GList **columns)
{
}

void
_nc_galaxy_shape_pop_gauss_data_init (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data)
{
  NcGalaxyShapePopGaussLData *ldata = g_new0 (NcGalaxyShapePopGaussLData, 1);

  data->ldata                  = ldata;
  data->ldata_destroy          = &g_free;
  data->ldata_read_row         = &_nc_galaxy_shape_pop_gauss_ldata_noop;
  data->ldata_write_row        = &_nc_galaxy_shape_pop_gauss_ldata_noop;
  data->ldata_required_columns = &_nc_galaxy_shape_pop_gauss_ldata_required_columns;
  data->ldata_get_sigma        = &_nc_galaxy_shape_pop_gauss_ldata_get_sigma;
}

gdouble
_nc_galaxy_shape_pop_gauss_ldata_get_sigma (NcGalaxyShapePopData *data)
{
  NcGalaxyShapePopGaussLData *ldata = (NcGalaxyShapePopGaussLData *) data->ldata;

  return ldata->sigma;
}

void
_nc_galaxy_shape_pop_gauss_ldata_set_sigma (NcGalaxyShapePopData *data, const gdouble sigma)
{
  NcGalaxyShapePopGaussLData *ldata = (NcGalaxyShapePopGaussLData *) data->ldata;
  const gdouble inv_2sigma2         = 0.5 / (sigma * sigma);

  /* P(x) ∝ exp(-x/2σ²) on [0,1]. */
  ldata->sigma       = sigma;
  ldata->inv_2sigma2 = inv_2sigma2;
  ldata->norm        = inv_2sigma2 / (-expm1 (-inv_2sigma2));
}

static void
_nc_galaxy_shape_pop_gauss_prepare (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data)
{
  _nc_galaxy_shape_pop_gauss_ldata_set_sigma (data, SIGMA);
}

gdouble
_nc_galaxy_shape_pop_gauss_eval_p (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data, const gdouble x)
{
  NcGalaxyShapePopGaussLData *ldata = (NcGalaxyShapePopGaussLData *) data->ldata;

  return ldata->norm * exp (-ldata->inv_2sigma2 * x);
}

static complex double
_gauss_cut_gen (NcmRNG *rng, const gdouble sigma)
{
  gdouble x;
  gdouble y;

  do {
    x = ncm_rng_gaussian_gen (rng, 0.0, sigma);
    y = ncm_rng_gaussian_gen (rng, 0.0, sigma);
  } while (hypot (x, y) > 1.0);

  return x + I * y;
}

void
_nc_galaxy_shape_pop_gauss_gen (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data, NcmRNG *rng, gdouble *e_int_1, gdouble *e_int_2)
{
  NcGalaxyShapePopGaussLData *ldata = (NcGalaxyShapePopGaussLData *) data->ldata;
  const complex double e            = _gauss_cut_gen (rng, ldata->sigma);

  *e_int_1 = creal (e);
  *e_int_2 = cimag (e);
}

static gdouble
_nc_galaxy_shape_pop_gauss_e_rms (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data)
{
  const gdouble sigma  = SIGMA;
  const gdouble lambda = 0.5 / (sigma * sigma);
  const gdouble exp_ml = exp (-lambda);
  /* <x> for the exponential P(x) ∝ exp(-lambda x) truncated to [0,1]. */
  const gdouble mean_x = (1.0 - exp_ml * (1.0 + lambda)) / (lambda * (1.0 - exp_ml));

  return sqrt (0.5 * mean_x);
}

/**
 * nc_galaxy_shape_pop_gauss_new:
 *
 * Creates a new #NcGalaxyShapePopGauss.
 *
 * Returns: (transfer full): a new #NcGalaxyShapePopGauss.
 */
NcGalaxyShapePopGauss *
nc_galaxy_shape_pop_gauss_new (void)
{
  NcGalaxyShapePopGauss *gspg = g_object_new (NC_TYPE_GALAXY_SHAPE_POP_GAUSS,
                                              NULL);

  return gspg;
}

/**
 * nc_galaxy_shape_pop_gauss_ref:
 * @gspg: a #NcGalaxyShapePopGauss
 *
 * Increases the reference count of @gspg by one.
 *
 * Returns: (transfer full): @gspg.
 */
NcGalaxyShapePopGauss *
nc_galaxy_shape_pop_gauss_ref (NcGalaxyShapePopGauss *gspg)
{
  return g_object_ref (gspg);
}

/**
 * nc_galaxy_shape_pop_gauss_free:
 * @gspg: a #NcGalaxyShapePopGauss
 *
 * Decreases the reference count of @gspg by one.
 *
 */
void
nc_galaxy_shape_pop_gauss_free (NcGalaxyShapePopGauss *gspg)
{
  g_object_unref (gspg);
}

/**
 * nc_galaxy_shape_pop_gauss_clear:
 * @gspg: a #NcGalaxyShapePopGauss
 *
 * Decreases the reference count of *@gspg by one, and sets the pointer *@gspg to
 * NULL.
 *
 */
void
nc_galaxy_shape_pop_gauss_clear (NcGalaxyShapePopGauss **gspg)
{
  g_clear_object (gspg);
}


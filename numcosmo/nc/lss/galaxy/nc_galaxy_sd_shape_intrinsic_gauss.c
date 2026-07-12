/***************************************************************************
 *            nc_galaxy_sd_shape_intrinsic_gauss.c
 *
 *  Thu Jun 19 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_galaxy_sd_shape_intrinsic_gauss.c
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
 * NcGalaxySDShapeIntrinsicGauss:
 *
 * Truncated-Gaussian intrinsic ellipticity distribution.
 *
 * The intrinsic ellipticity follows an isotropic Gaussian of width $\sigma$
 * truncated to the unit disk. The induced density of $x = |\chi_I|^2$ is
 * $P(x) \propto \exp(-x / 2\sigma^2)$ on $[0,1]$, so the Gauss-Jacobi weight
 * exponents vanish (plain Gauss-Legendre) and the residual carries the
 * normalized exponential.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc/lss/galaxy/nc_galaxy_sd_shape_intrinsic_gauss.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <math.h>
#include <complex.h>
#endif /* NUMCOSMO_GIR_SCAN */

struct _NcGalaxySDShapeIntrinsicGauss
{
  NcGalaxySDShapeIntrinsic parent_instance;
};

typedef struct _NcGalaxySDShapeIntrinsicGaussLData
{
  gdouble norm;        /* residual normalization */
  gdouble inv_2sigma2; /* 1 / (2 sigma^2) */
  gdouble sigma;       /* truncated-Gaussian width (for sampling) */
} NcGalaxySDShapeIntrinsicGaussLData;

enum
{
  PROP_0,
  PROP_LEN,
};

G_DEFINE_TYPE (NcGalaxySDShapeIntrinsicGauss, nc_galaxy_sd_shape_intrinsic_gauss, NC_TYPE_GALAXY_SD_SHAPE_INTRINSIC);

static void
nc_galaxy_sd_shape_intrinsic_gauss_init (NcGalaxySDShapeIntrinsicGauss *gsig)
{
}

static void
_nc_galaxy_sd_shape_intrinsic_gauss_finalize (GObject *object)
{
  G_OBJECT_CLASS (nc_galaxy_sd_shape_intrinsic_gauss_parent_class)->finalize (object);
}

static void _nc_galaxy_sd_shape_intrinsic_gauss_data_init (NcGalaxySDShapeIntrinsic *gsi, NcGalaxySDShapeIntrinsicData *data);
static void _nc_galaxy_sd_shape_intrinsic_gauss_prepare (NcGalaxySDShapeIntrinsic *gsi, NcGalaxySDShapeIntrinsicData *data);
static gdouble _nc_galaxy_sd_shape_intrinsic_gauss_eval_residual (NcGalaxySDShapeIntrinsic *gsi, NcGalaxySDShapeIntrinsicData *data, const gdouble x);
static void _nc_galaxy_sd_shape_intrinsic_gauss_gen (NcGalaxySDShapeIntrinsic *gsi, NcGalaxySDShapeIntrinsicData *data, NcmRNG *rng, gdouble *e_int_1, gdouble *e_int_2);
static gdouble _nc_galaxy_sd_shape_intrinsic_gauss_e_rms (NcGalaxySDShapeIntrinsic *gsi, NcGalaxySDShapeIntrinsicData *data);

static void
nc_galaxy_sd_shape_intrinsic_gauss_class_init (NcGalaxySDShapeIntrinsicGaussClass *klass)
{
  NcGalaxySDShapeIntrinsicClass *gsi_class = NC_GALAXY_SD_SHAPE_INTRINSIC_CLASS (klass);
  GObjectClass *object_class             = G_OBJECT_CLASS (klass);
  NcmModelClass *model_class             = NCM_MODEL_CLASS (klass);

  object_class->finalize = &_nc_galaxy_sd_shape_intrinsic_gauss_finalize;

  ncm_model_class_set_name_nick (model_class, "Truncated Gaussian intrinsic ellipticity distribution", "GaussIntrinsic");
  ncm_model_class_add_params (model_class, NC_GALAXY_SD_SHAPE_INTRINSIC_GAUSS_SPARAM_LEN, 0, PROP_LEN);

  /**
   * NcGalaxySDShapeIntrinsicGauss:sigma:
   *
   * The width $\sigma$ of the truncated Gaussian intrinsic ellipticity.
   *
   */
  ncm_model_class_set_sparam (model_class,
                              NC_GALAXY_SD_SHAPE_INTRINSIC_GAUSS_SIGMA,
                              "\\sigma",
                              "sigma", 1.0e-2, 1.0, 1.0e-1,
                              NC_GALAXY_SD_SHAPE_INTRINSIC_GAUSS_DEFAULT_PARAMS_ABSTOL,
                              NC_GALAXY_SD_SHAPE_INTRINSIC_GAUSS_DEFAULT_SIGMA,
                              NCM_PARAM_TYPE_FIXED);

  ncm_model_class_check_params_info (model_class);

  gsi_class->data_init     = &_nc_galaxy_sd_shape_intrinsic_gauss_data_init;
  gsi_class->prepare       = &_nc_galaxy_sd_shape_intrinsic_gauss_prepare;
  gsi_class->eval_residual = &_nc_galaxy_sd_shape_intrinsic_gauss_eval_residual;
  gsi_class->gen           = &_nc_galaxy_sd_shape_intrinsic_gauss_gen;
  gsi_class->e_rms         = &_nc_galaxy_sd_shape_intrinsic_gauss_e_rms;
}

#define VECTOR (NCM_MODEL (gsi))
#define SIGMA  (ncm_model_orig_param_get (VECTOR, NC_GALAXY_SD_SHAPE_INTRINSIC_GAUSS_SIGMA))

static void
_nc_galaxy_sd_shape_intrinsic_gauss_ldata_noop (NcGalaxySDShapeIntrinsicData *data, NcGalaxyWLObs *obs, const guint i)
{
}

static void
_nc_galaxy_sd_shape_intrinsic_gauss_ldata_required_columns (NcGalaxySDShapeIntrinsicData *data, GList *columns)
{
}

static void
_nc_galaxy_sd_shape_intrinsic_gauss_data_init (NcGalaxySDShapeIntrinsic *gsi, NcGalaxySDShapeIntrinsicData *data)
{
  NcGalaxySDShapeIntrinsicGaussLData *ldata = g_new0 (NcGalaxySDShapeIntrinsicGaussLData, 1);

  data->ldata                  = ldata;
  data->ldata_destroy          = &g_free;
  data->ldata_read_row         = &_nc_galaxy_sd_shape_intrinsic_gauss_ldata_noop;
  data->ldata_write_row        = &_nc_galaxy_sd_shape_intrinsic_gauss_ldata_noop;
  data->ldata_required_columns = &_nc_galaxy_sd_shape_intrinsic_gauss_ldata_required_columns;
}

static void
_nc_galaxy_sd_shape_intrinsic_gauss_prepare (NcGalaxySDShapeIntrinsic *gsi, NcGalaxySDShapeIntrinsicData *data)
{
  NcGalaxySDShapeIntrinsicGaussLData *ldata = (NcGalaxySDShapeIntrinsicGaussLData *) data->ldata;
  const gdouble sigma                     = SIGMA;
  const gdouble inv_2sigma2               = 0.5 / (sigma * sigma);

  /* P(x) ∝ exp(-x/2σ²) on [0,1]: plain Gauss-Legendre weight (a=b=0) and the
   * residual carries the normalized exponential. */
  data->jacobi_a     = 0.0;
  data->jacobi_b     = 0.0;
  ldata->sigma       = sigma;
  ldata->inv_2sigma2 = inv_2sigma2;
  ldata->norm        = inv_2sigma2 / (-expm1 (-inv_2sigma2));
}

static gdouble
_nc_galaxy_sd_shape_intrinsic_gauss_eval_residual (NcGalaxySDShapeIntrinsic *gsi, NcGalaxySDShapeIntrinsicData *data, const gdouble x)
{
  NcGalaxySDShapeIntrinsicGaussLData *ldata = (NcGalaxySDShapeIntrinsicGaussLData *) data->ldata;

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

static void
_nc_galaxy_sd_shape_intrinsic_gauss_gen (NcGalaxySDShapeIntrinsic *gsi, NcGalaxySDShapeIntrinsicData *data, NcmRNG *rng, gdouble *e_int_1, gdouble *e_int_2)
{
  NcGalaxySDShapeIntrinsicGaussLData *ldata = (NcGalaxySDShapeIntrinsicGaussLData *) data->ldata;
  const complex double e                   = _gauss_cut_gen (rng, ldata->sigma);

  *e_int_1 = creal (e);
  *e_int_2 = cimag (e);
}

static gdouble
_nc_galaxy_sd_shape_intrinsic_gauss_e_rms (NcGalaxySDShapeIntrinsic *gsi, NcGalaxySDShapeIntrinsicData *data)
{
  const gdouble sigma   = SIGMA;
  const gdouble lambda  = 0.5 / (sigma * sigma);
  const gdouble exp_ml  = exp (-lambda);
  /* <x> for the exponential P(x) ∝ exp(-lambda x) truncated to [0,1]. */
  const gdouble mean_x  = (1.0 - exp_ml * (1.0 + lambda)) / (lambda * (1.0 - exp_ml));

  return sqrt (0.5 * mean_x);
}

/**
 * nc_galaxy_sd_shape_intrinsic_gauss_new:
 *
 * Creates a new #NcGalaxySDShapeIntrinsicGauss.
 *
 * Returns: (transfer full): a new #NcGalaxySDShapeIntrinsicGauss.
 */
NcGalaxySDShapeIntrinsicGauss *
nc_galaxy_sd_shape_intrinsic_gauss_new (void)
{
  NcGalaxySDShapeIntrinsicGauss *gsig = g_object_new (NC_TYPE_GALAXY_SD_SHAPE_INTRINSIC_GAUSS,
                                                    NULL);

  return gsig;
}

/**
 * nc_galaxy_sd_shape_intrinsic_gauss_ref:
 * @gsig: a #NcGalaxySDShapeIntrinsicGauss
 *
 * Increases the reference count of @gsig by one.
 *
 * Returns: (transfer full): @gsig.
 */
NcGalaxySDShapeIntrinsicGauss *
nc_galaxy_sd_shape_intrinsic_gauss_ref (NcGalaxySDShapeIntrinsicGauss *gsig)
{
  return g_object_ref (gsig);
}

/**
 * nc_galaxy_sd_shape_intrinsic_gauss_free:
 * @gsig: a #NcGalaxySDShapeIntrinsicGauss
 *
 * Decreases the reference count of @gsig by one.
 *
 */
void
nc_galaxy_sd_shape_intrinsic_gauss_free (NcGalaxySDShapeIntrinsicGauss *gsig)
{
  g_object_unref (gsig);
}

/**
 * nc_galaxy_sd_shape_intrinsic_gauss_clear:
 * @gsig: a #NcGalaxySDShapeIntrinsicGauss
 *
 * Decreases the reference count of *@gsig by one, and sets the pointer *@gsig to
 * NULL.
 *
 */
void
nc_galaxy_sd_shape_intrinsic_gauss_clear (NcGalaxySDShapeIntrinsicGauss **gsig)
{
  g_clear_object (gsig);
}

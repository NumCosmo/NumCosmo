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
 * truncated to the unit disk: the 2D area density is
 * $P_\mathrm{2D}(\chi_I) \propto \exp(-r^2/2\sigma^2)$, $r=|\chi_I|$. The
 * public r-marginal contract (see #NcGalaxyShapePop's own docs)
 * is $P_\mathrm{pop}(r) = 2\pi r\,P_\mathrm{2D}(r) \propto r\exp(-r^2/2\sigma^2)$
 * on $[0,1)$ -- multiplying the natural area-density form by $r$, always
 * safe (never a pole), unlike populations whose natural form is the
 * r-marginal itself (see #NcGalaxyShapePopBeta).
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

/* _data_init, _eval_p and _gen are shared with related models (see
 * nc_galaxy_shape_pop_gauss_private.h) so are not static; forward-declared
 * here only to keep definition order readable. */
void _nc_galaxy_shape_pop_gauss_data_init (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data);
gdouble _nc_galaxy_shape_pop_gauss_eval_p (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data, const gdouble r);
void _nc_galaxy_shape_pop_gauss_gen (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data, NcmRNG *rng, gdouble *e_int_1, gdouble *e_int_2);
static void _nc_galaxy_shape_pop_gauss_prepare (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data);
static gdouble _nc_galaxy_shape_pop_gauss_e_rms (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data);
static gdouble _nc_galaxy_shape_pop_gauss_exponent_at_origin (NcGalaxyShapePop *gsp);
gdouble _nc_galaxy_shape_pop_gauss_moment_2k (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data, const guint k);

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

  gsp_class->data_init          = &_nc_galaxy_shape_pop_gauss_data_init;
  gsp_class->prepare            = &_nc_galaxy_shape_pop_gauss_prepare;
  gsp_class->eval_p             = &_nc_galaxy_shape_pop_gauss_eval_p;
  gsp_class->gen                = &_nc_galaxy_shape_pop_gauss_gen;
  gsp_class->e_rms              = &_nc_galaxy_shape_pop_gauss_e_rms;
  gsp_class->exponent_at_origin = &_nc_galaxy_shape_pop_gauss_exponent_at_origin;
  gsp_class->moment_2k          = &_nc_galaxy_shape_pop_gauss_moment_2k;
  gsp_class->exponent_at_one    = &_nc_galaxy_shape_pop_gauss_exponent_at_one;
  gsp_class->eval_p_reduced     = &_nc_galaxy_shape_pop_gauss_eval_p_reduced;
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

  /* norm/inv_2sigma2 are the x-space (x=r^2) normalization,
   * P(x) = norm*exp(-x/2 sigma^2) on [0,1], used by eval_p() below via the *2r conversion. */
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
_nc_galaxy_shape_pop_gauss_eval_p (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data, const gdouble r)
{
  NcGalaxyShapePopGaussLData *ldata = (NcGalaxyShapePopGaussLData *) data->ldata;

  /* P_pop(r) = 2*pi*r*P_2D(r) = 2*r*norm*exp(-r^2/2sigma^2) -- norm is
   * unchanged from the area-density form (see _ldata_set_sigma()'s own
   * docs), the substitution u=r^2 recovers exactly the integral norm's own
   * definition already satisfies, so no new normalization constant is
   * needed here. */
  return 2.0 * r * ldata->norm * exp (-ldata->inv_2sigma2 * r * r);
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
  /* <x> for the exponential P(x) proportional to exp(-lambda x), truncated to [0,1]. */
  const gdouble mean_x = (1.0 - exp_ml * (1.0 + lambda)) / (lambda * (1.0 - exp_ml));

  return sqrt (0.5 * mean_x);
}

static gdouble
_nc_galaxy_shape_pop_gauss_exponent_at_origin (NcGalaxyShapePop *gsp)
{
  return 1.0;
}

/* The truncated Gaussian is finite and nonzero at r = 1 */
gdouble
_nc_galaxy_shape_pop_gauss_exponent_at_one (NcGalaxyShapePop *gsp)
{
  return 0.0;
}

/* P_pop (r) / r = 2 norm exp (-r^2 / 2 sigma^2), see _nc_galaxy_shape_pop_gauss_eval_p() */
gdouble
_nc_galaxy_shape_pop_gauss_eval_p_reduced (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data, const gdouble r)
{
  NcGalaxyShapePopGaussLData *ldata = (NcGalaxyShapePopGaussLData *) data->ldata;

  return 2.0 * ldata->norm * exp (-ldata->inv_2sigma2 * r * r);
}

/* M_2k = <r^2k> = (2sigma^2)^k k! P(k+1,a)/(1-e^-a), a = inv_2sigma2, in
 * its integer-k form: P(k+1,a) = 1 - e^-a sum_{j=0}^k a^j/j!, and
 * (2sigma^2)^k k! = 1/term_k with term_j = a^j/j!, the value the loop below
 * leaves in @term. For integer k this agrees with gsl_sf_gamma_inc_P() to
 * machine precision. Reads the resolved @data->ldata, as eval_p() does. */
gdouble
_nc_galaxy_shape_pop_gauss_moment_2k (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data, const guint k)
{
  if (k == 0)
  {
    return 1.0;
  }
  else
  {
    NcGalaxyShapePopGaussLData *ldata = (NcGalaxyShapePopGaussLData *) data->ldata;
    const gdouble a                   = ldata->inv_2sigma2;
    const gdouble exp_ma              = exp (-a);
    gdouble term                      = 1.0; /* a^j/j!, j=0 */
    gdouble sum                       = 1.0;
    guint j;

    for (j = 1; j <= k; j++)
    {
      term *= a / j;
      sum  += term;
    }

    return (1.0 - exp_ma * sum) / (term * (1.0 - exp_ma));
  }
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


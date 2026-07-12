/***************************************************************************
 *            nc_galaxy_sd_shape_intrinsic_beta.c
 *
 *  Thu Jun 19 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_galaxy_sd_shape_intrinsic_beta.c
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
 * NcGalaxySDShapeIntrinsicBeta:
 *
 * Beta intrinsic ellipticity distribution (BetaGlobal).
 *
 * The squared intrinsic ellipticity modulus $x = |\chi_I|^2$ follows a Beta
 * distribution $x \sim \mathrm{Beta}(\alpha,\beta)$ with $\alpha = \mu\nu$ and
 * $\beta = (1-\mu)\nu$, where $\mu$ controls the typical ellipticity and $\nu$
 * the concentration. This is the reference exact intrinsic model: it has exact
 * support on $[0,1]$ and maps naturally to Gauss-Jacobi quadrature, with weight
 * exponents $a = \beta-1$ (for $(1-x)$) and $b = \alpha-1$ (for $x$) and a
 * constant residual $1/B(\alpha,\beta)$.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc/lss/galaxy/nc_galaxy_sd_shape_intrinsic_beta.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <math.h>
#include <complex.h>
#include <gsl/gsl_sf_gamma.h>
#endif /* NUMCOSMO_GIR_SCAN */

struct _NcGalaxySDShapeIntrinsicBeta
{
  NcGalaxySDShapeIntrinsic parent_instance;
};

typedef struct _NcGalaxySDShapeIntrinsicBetaLData
{
  gdouble norm;  /* 1 / B(alpha, beta) */
  gdouble alpha; /* mu nu */
  gdouble beta;  /* (1 - mu) nu */
} NcGalaxySDShapeIntrinsicBetaLData;

enum
{
  PROP_0,
  PROP_LEN,
};

G_DEFINE_TYPE (NcGalaxySDShapeIntrinsicBeta, nc_galaxy_sd_shape_intrinsic_beta, NC_TYPE_GALAXY_SD_SHAPE_INTRINSIC);

static void
nc_galaxy_sd_shape_intrinsic_beta_init (NcGalaxySDShapeIntrinsicBeta *gsib)
{
}

static void
_nc_galaxy_sd_shape_intrinsic_beta_finalize (GObject *object)
{
  G_OBJECT_CLASS (nc_galaxy_sd_shape_intrinsic_beta_parent_class)->finalize (object);
}

static void _nc_galaxy_sd_shape_intrinsic_beta_data_init (NcGalaxySDShapeIntrinsic *gsi, NcGalaxySDShapeIntrinsicData *data);
static void _nc_galaxy_sd_shape_intrinsic_beta_prepare (NcGalaxySDShapeIntrinsic *gsi, NcGalaxySDShapeIntrinsicData *data);
static gdouble _nc_galaxy_sd_shape_intrinsic_beta_eval_residual (NcGalaxySDShapeIntrinsic *gsi, NcGalaxySDShapeIntrinsicData *data, const gdouble x);
static void _nc_galaxy_sd_shape_intrinsic_beta_gen (NcGalaxySDShapeIntrinsic *gsi, NcGalaxySDShapeIntrinsicData *data, NcmRNG *rng, gdouble *e_int_1, gdouble *e_int_2);
static gdouble _nc_galaxy_sd_shape_intrinsic_beta_e_rms (NcGalaxySDShapeIntrinsic *gsi, NcGalaxySDShapeIntrinsicData *data);

static void
nc_galaxy_sd_shape_intrinsic_beta_class_init (NcGalaxySDShapeIntrinsicBetaClass *klass)
{
  NcGalaxySDShapeIntrinsicClass *gsi_class = NC_GALAXY_SD_SHAPE_INTRINSIC_CLASS (klass);
  GObjectClass *object_class             = G_OBJECT_CLASS (klass);
  NcmModelClass *model_class             = NCM_MODEL_CLASS (klass);

  object_class->finalize = &_nc_galaxy_sd_shape_intrinsic_beta_finalize;

  ncm_model_class_set_name_nick (model_class, "Beta intrinsic ellipticity distribution", "BetaIntrinsic");
  ncm_model_class_add_params (model_class, NC_GALAXY_SD_SHAPE_INTRINSIC_BETA_SPARAM_LEN, 0, PROP_LEN);

  /**
   * NcGalaxySDShapeIntrinsicBeta:mu:
   *
   * The mean $\mu = \langle x \rangle$ of the Beta distribution of $x = |\chi_I|^2$.
   *
   */
  ncm_model_class_set_sparam (model_class,
                              NC_GALAXY_SD_SHAPE_INTRINSIC_BETA_MU,
                              "\\mu",
                              "mu", 1.0e-3, 1.0 - 1.0e-3, 1.0e-2,
                              NC_GALAXY_SD_SHAPE_INTRINSIC_BETA_DEFAULT_PARAMS_ABSTOL,
                              NC_GALAXY_SD_SHAPE_INTRINSIC_BETA_DEFAULT_MU,
                              NCM_PARAM_TYPE_FIXED);

  /**
   * NcGalaxySDShapeIntrinsicBeta:nu:
   *
   * The concentration $\nu = \alpha + \beta$ of the Beta distribution.
   *
   */
  ncm_model_class_set_sparam (model_class,
                              NC_GALAXY_SD_SHAPE_INTRINSIC_BETA_NU,
                              "\\nu",
                              "nu", 1.0e-2, 1.0e3, 1.0e-1,
                              NC_GALAXY_SD_SHAPE_INTRINSIC_BETA_DEFAULT_PARAMS_ABSTOL,
                              NC_GALAXY_SD_SHAPE_INTRINSIC_BETA_DEFAULT_NU,
                              NCM_PARAM_TYPE_FIXED);

  ncm_model_class_check_params_info (model_class);

  gsi_class->data_init     = &_nc_galaxy_sd_shape_intrinsic_beta_data_init;
  gsi_class->prepare       = &_nc_galaxy_sd_shape_intrinsic_beta_prepare;
  gsi_class->eval_residual = &_nc_galaxy_sd_shape_intrinsic_beta_eval_residual;
  gsi_class->gen           = &_nc_galaxy_sd_shape_intrinsic_beta_gen;
  gsi_class->e_rms         = &_nc_galaxy_sd_shape_intrinsic_beta_e_rms;
}

#define VECTOR (NCM_MODEL (gsi))
#define MU     (ncm_model_orig_param_get (VECTOR, NC_GALAXY_SD_SHAPE_INTRINSIC_BETA_MU))
#define NU     (ncm_model_orig_param_get (VECTOR, NC_GALAXY_SD_SHAPE_INTRINSIC_BETA_NU))

static void
_nc_galaxy_sd_shape_intrinsic_beta_ldata_noop (NcGalaxySDShapeIntrinsicData *data, NcGalaxyWLObs *obs, const guint i)
{
}

static void
_nc_galaxy_sd_shape_intrinsic_beta_ldata_required_columns (NcGalaxySDShapeIntrinsicData *data, GList *columns)
{
}

static void
_nc_galaxy_sd_shape_intrinsic_beta_data_init (NcGalaxySDShapeIntrinsic *gsi, NcGalaxySDShapeIntrinsicData *data)
{
  NcGalaxySDShapeIntrinsicBetaLData *ldata = g_new0 (NcGalaxySDShapeIntrinsicBetaLData, 1);

  data->ldata                  = ldata;
  data->ldata_destroy          = &g_free;
  data->ldata_read_row         = &_nc_galaxy_sd_shape_intrinsic_beta_ldata_noop;
  data->ldata_write_row        = &_nc_galaxy_sd_shape_intrinsic_beta_ldata_noop;
  data->ldata_required_columns = &_nc_galaxy_sd_shape_intrinsic_beta_ldata_required_columns;
}

static void
_nc_galaxy_sd_shape_intrinsic_beta_prepare (NcGalaxySDShapeIntrinsic *gsi, NcGalaxySDShapeIntrinsicData *data)
{
  NcGalaxySDShapeIntrinsicBetaLData *ldata = (NcGalaxySDShapeIntrinsicBetaLData *) data->ldata;
  const gdouble mu                       = MU;
  const gdouble nu                       = NU;
  const gdouble alpha                    = mu * nu;
  const gdouble beta                     = (1.0 - mu) * nu;

  /* P(x) = x^{alpha-1} (1-x)^{beta-1} / B(alpha,beta) = (1-x)^a x^b R, with the
   * Gauss-Jacobi weight (1-x)^a x^b and constant residual R = 1/B(alpha,beta). */
  data->jacobi_a = beta - 1.0;
  data->jacobi_b = alpha - 1.0;
  ldata->norm    = exp (-gsl_sf_lnbeta (alpha, beta));
  ldata->alpha   = alpha;
  ldata->beta    = beta;
}

static gdouble
_nc_galaxy_sd_shape_intrinsic_beta_eval_residual (NcGalaxySDShapeIntrinsic *gsi, NcGalaxySDShapeIntrinsicData *data, const gdouble x)
{
  NcGalaxySDShapeIntrinsicBetaLData *ldata = (NcGalaxySDShapeIntrinsicBetaLData *) data->ldata;

  return ldata->norm;
}

static void
_nc_galaxy_sd_shape_intrinsic_beta_gen (NcGalaxySDShapeIntrinsic *gsi, NcGalaxySDShapeIntrinsicData *data, NcmRNG *rng, gdouble *e_int_1, gdouble *e_int_2)
{
  NcGalaxySDShapeIntrinsicBetaLData *ldata = (NcGalaxySDShapeIntrinsicBetaLData *) data->ldata;
  const gdouble x                        = ncm_rng_beta_gen (rng, ldata->alpha, ldata->beta);
  const gdouble r                        = sqrt (x);
  const gdouble theta                    = ncm_rng_uniform_gen (rng, 0.0, 2.0 * M_PI);

  *e_int_1 = r * cos (theta);
  *e_int_2 = r * sin (theta);
}

static gdouble
_nc_galaxy_sd_shape_intrinsic_beta_e_rms (NcGalaxySDShapeIntrinsic *gsi, NcGalaxySDShapeIntrinsicData *data)
{
  /* <|chi|^2> = <x> = mu, so e_rms = sqrt(mu / 2). */
  return sqrt (0.5 * MU);
}

/**
 * nc_galaxy_sd_shape_intrinsic_beta_new:
 *
 * Creates a new #NcGalaxySDShapeIntrinsicBeta.
 *
 * Returns: (transfer full): a new #NcGalaxySDShapeIntrinsicBeta.
 */
NcGalaxySDShapeIntrinsicBeta *
nc_galaxy_sd_shape_intrinsic_beta_new (void)
{
  NcGalaxySDShapeIntrinsicBeta *gsib = g_object_new (NC_TYPE_GALAXY_SD_SHAPE_INTRINSIC_BETA,
                                                   NULL);

  return gsib;
}

/**
 * nc_galaxy_sd_shape_intrinsic_beta_ref:
 * @gsib: a #NcGalaxySDShapeIntrinsicBeta
 *
 * Increases the reference count of @gsib by one.
 *
 * Returns: (transfer full): @gsib.
 */
NcGalaxySDShapeIntrinsicBeta *
nc_galaxy_sd_shape_intrinsic_beta_ref (NcGalaxySDShapeIntrinsicBeta *gsib)
{
  return g_object_ref (gsib);
}

/**
 * nc_galaxy_sd_shape_intrinsic_beta_free:
 * @gsib: a #NcGalaxySDShapeIntrinsicBeta
 *
 * Decreases the reference count of @gsib by one.
 *
 */
void
nc_galaxy_sd_shape_intrinsic_beta_free (NcGalaxySDShapeIntrinsicBeta *gsib)
{
  g_object_unref (gsib);
}

/**
 * nc_galaxy_sd_shape_intrinsic_beta_clear:
 * @gsib: a #NcGalaxySDShapeIntrinsicBeta
 *
 * Decreases the reference count of *@gsib by one, and sets the pointer *@gsib to
 * NULL.
 *
 */
void
nc_galaxy_sd_shape_intrinsic_beta_clear (NcGalaxySDShapeIntrinsicBeta **gsib)
{
  g_clear_object (gsib);
}
